from pypipegraph import (
    Job,
    # MultiFileGeneratingJob,
    FileGeneratingJob,
    # ParameterInvariant,
)
from typing import List, Optional  # Callable, Dict, Tuple
from pathlib import Path
from .util import dump_matching_reads
from .plots import make_adapter_heatmap


def plot_adapter_heatmap(outfile: Path, matrix_file: Path, dependencies: List[Job],
        value_column: str = "both",
        log10: bool = True,
        figsize: tuple = (10, 8),
        dpi: int = 150,
        cmap: str = "viridis",
        annot: bool = False,
        annot_fmt: str = "g",
        vmax: Optional[float] = None,

    ):
    def __dump(outfile, matrix_file=matrix_file):
        make_adapter_heatmap(
            matrix_tsv=matrix_file,
            out_file=outfile,
            value_column=value_column,
            log10=log10,
            figsize=figsize,
            dpi=dpi,
            cmap=cmap,
            annot=annot,
            annot_fmt=annot_fmt,
            vmax=vmax,
        )
    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def dump_matching_reads_job(
    incoming_reads_file: Path,
    outfile: Path,
    r1: Path,
    r2: Path = None,
    max_reads: Optional[int] = None,
    dependencies: List[Job] = [],
):
    def __dump(
        outfile,
        incoming_reads_file=incoming_reads_file,
        r1=r1,
        r2=r2,
        max_reads=max_reads,
    ):
        dump_matching_reads(incoming_reads_file, Path(outfile), r1, r2, max_reads)

    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)
