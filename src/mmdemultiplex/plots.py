#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains utility functions for the demultiplexer package."""
import subprocess
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# import pypipegraph as ppg
import numpy as np
import seaborn as sns

# from mplots import MPPlotJob
from pypipegraph import MultiFileGeneratingJob, Job
from .demultiplex import Demultiplexer
from typing import Union, List, Optional
from pathlib import Path


def make_adapter_heatmap(
    matrix_tsv: Union[str, Path],
    out_file: Union[str, Path],
    value_column: str = "both",
    log10: bool = True,
    figsize: tuple = (10, 8),
    dpi: int = 150,
    cmap: str = "viridis",
    annot: bool = False,
    annot_fmt: str = "g",
    vmax: Optional[float] = None,
):
    """
    Erzeugt eine Heatmap der Barcode-Kombinationen aus der
    adapter_check_matrix_main.tsv.

    Parameter
    ---------
    matrix_tsv : str | Path
        Pfad zu der TSV-Datei, die von _adapter_check_matrix_main.tsv erzeugt wurde.
        Erwartete Spalten: start_barcode, end_barcode, <value_column>, ...
    out_file : str | Path
        Pfad zur Output-Grafik (z.B. 'heatmap.png' oder 'heatmap.pdf').
    value_column : str, default 'both'
        Welche Metrik aus der Matrix-Datei geplottet werden soll
        (z.B. 'both', 'any', 'expected', 'swapped', ...).
    log10 : bool, default True
        Wenn True, werden die Werte mit log10(value+1) transformiert,
        was oft sinnvoll ist, wenn die Counts stark streuen.
    figsize : tuple, default (10, 8)
        Größe der Figure in Zoll.
    dpi : int, default 150
        Auflösung der gespeicherten Grafik.
    cmap : str, default 'viridis'
        Farbskala für die Heatmap.
    annot : bool, default False
        Wenn True, werden die Werte in die Zellen geschrieben.
        (Vorsicht bei großen Matrizen.)
    annot_fmt : str, default 'g'
        Format-String für die Annotation.
    vmax : float | None, default None
        Optionaler Maximalwert für die Farbskalierung.
    """
    matrix_tsv = Path(matrix_tsv)
    out_file = Path(out_file)

    if not matrix_tsv.exists():
        raise FileNotFoundError(f"Matrix file not found: {matrix_tsv}")

    df = pd.read_csv(matrix_tsv, sep="\t")

    required_cols = {"start_barcode", "end_barcode", value_column}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in matrix file: {', '.join(sorted(missing))}")

    # Pivot: Zeilen = start_barcode, Spalten = end_barcode
    mat = df.pivot_table(
        index="start_barcode",
        columns="end_barcode",
        values=value_column,
        aggfunc="sum",
        fill_value=0,
    )

    if log10:
        # log10(value + 1), um mit 0 klarzukommen
        mat_plot = np.log10(mat + 1)
    else:
        mat_plot = mat

    plt.figure(figsize=figsize)
    sns.heatmap(
        mat_plot,
        cmap=cmap,
        annot=annot,
        fmt=annot_fmt,
        cbar=True,
        vmax=vmax,
    )
    plt.title(f"Barcode heatmap ({value_column}{' (log10)' if log10 else ''})")
    plt.xlabel("end_barcode")
    plt.ylabel("start_barcode")
    plt.tight_layout()

    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=dpi)
    plt.close()


def count_raw_input_reads(gz_filename1):
    if gz_filename1.endswith(".gz"):
        p1 = subprocess.Popen(["gunzip", "-c", gz_filename1], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
        x = int(p2.communicate()[0][:-1])
        x = x / 4
        return x
    else:
        t = subprocess.check_output(["wc", "-l", gz_filename1])
        if t == b"":
            return 0
        else:
            x = int(t.split()[0])
            x = x / 4
            return x


def plot_read_counts_callable(demultiplexer: List[Demultiplexer], log: bool = True):

    def __count(filenames, demultiplexer=demultiplexer, log=log):
        outfile, outfile_df = filenames
        samples_to_plot = []
        color = []
        for dm in demultiplexer:
            samples_to_plot.append(dm.input_sample)
            color.append("g")
            dm_samples = dm.get_samples()
            for dm_sample_name in dm_samples:
                samples_to_plot.append(dm_samples[dm_sample_name])
                color.append("c")
        tmp = {"Sample": [], "Count": []}
        for sample in samples_to_plot:
            read_count = count_raw_input_reads(
                str(sample.get_aligner_input_filenames()[0])
            )
            tmp["Sample"].append(sample.name)
            tmp["Count"].append(read_count)
        df = pd.DataFrame(tmp)
        df["Color"] = color
        df.to_csv(outfile_df, sep="\t", index=False)

        fig = plt.figure(figsize=(12, 12))
        x = np.arange(len(samples_to_plot))
        plt.bar(x, df["Count"].values, width=0.8, color=df["Color"].values)
        labels = [str(x.strip()) for x in df["Sample"].values]
        plt.xticks(ticks=x, labels=labels, rotation=75, ha="right")
        if log:
            plt.yscale("log")
        plt.title("Demultiplexed read counts")
        plt.tight_layout()
        fig.savefig(outfile)

    return __count


def plot_read_counts(demultiplexer_or_list, outfile, dependencies=[], log=True):
    outfile_df = outfile + ".tsv"
    demultiplexer = demultiplexer_or_list
    deps = []
    if isinstance(demultiplexer_or_list, Demultiplexer):
        demultiplexer = [demultiplexer]
    for dm in demultiplexer:
        deps.append(dm.do_demultiplex())
        deps.append(dm.input_sample.prepare_input())
        for sample in dm.get_samples().values():
            deps.append(sample.prepare_input())
    return MultiFileGeneratingJob(
        [outfile, outfile_df],
        plot_read_counts_callable(demultiplexer, log),
    ).depends_on(dependencies + deps)


def plot_read_counts_percentage(
    demultiplexer_or_list, outfile, dependencies=[], log=True
):
    demultiplexer = demultiplexer_or_list
    deps = []
    if isinstance(demultiplexer_or_list, Demultiplexer):
        demultiplexer = [demultiplexer]
    for dm in demultiplexer:
        deps.append(dm.do_demultiplex())
        deps.append(dm.input_sample.prepare_input())
        for sample in dm.get_samples().values():
            deps.append(sample.prepare_input())

    def __count():
        samples_to_plot = {}
        color = []
        tmp = {"Sample": [], "Percentage Reads": [], "Color": []}
        for dm in demultiplexer:
            samples_to_plot[dm.input_sample] = [dm.input_sample]
            read_count_100 = count_raw_input_reads(
                str(dm.input_sample.get_aligner_input_filenames()[0])
            )
            tmp["Color"].append("g")
            tmp["Sample"].append(dm.input_sample.name)
            tmp["Percentage Reads"].append(100.0)
            dm_samples = dm.get_samples()
            for dm_sample_name in dm_samples:
                sample = dm_samples[dm_sample_name]
                read_count = count_raw_input_reads(
                    str(sample.get_aligner_input_filenames()[0])
                )
                perc = 100 * (float(read_count) / read_count_100)
                tmp["Color"].append("c")
                tmp["Sample"].append(sample.name)
                tmp["Percentage Reads"].append(perc)
        df = pd.DataFrame(tmp)
        return df

    def __plot(df):
        fig = plt.figure(figsize=(12, 12))
        x = np.arange(df.shape[0])
        plt.barh(x, df["Percentage Reads"].values, height=0.8, color=df["Color"].values)
        for i, v in enumerate(df["Percentage Reads"].values):
            plt.gca().text(
                v + 2,
                i - 0.3,
                f"{v:.2f}%",
                color=df["Color"].values[i],
            )
        labels = [str(x.strip()) for x in df["Sample"].values]
        plt.yticks(ticks=x, labels=labels, ha="right")
        plt.xlim([0, 115])
        plt.title("Demultiplexed read count percentages")
        plt.tight_layout()
        return fig

    return MPPlotJob(outfile, __count, __plot).depends_on(dependencies + deps)
