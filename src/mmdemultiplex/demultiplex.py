#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
demultiplex.py: This contains a Demultiplexer class that is used for
demultiplexing fastqs. It offers the option to find barcodes within the sequences,
thereby removing trailing sequences. In addition, it can trim fixed length
sequences after the detected barcode. Can be applied to single end paired end
sequencing. Demultiplexing can Be defined on Forward barcodes, reverse barcodes
and a combination of both (in case of paired end reads.
"""

from io import FileIO
import pandas as pd
import mbf_align
import pypipegraph as ppg
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from mbf_align import Sample  # TODO: replace with something else
from .strategies import PE_Decide_On_Start_Trim_Start_End, DemultiplexStrategy
from .util import Fragment, Read, get_fastq_iterator
from pypipegraph import Job


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


from pathlib import Path
import tempfile
import shutil


class TemporaryToPermanent:
    def __init__(self, permanent_file: Path):
        self.permanent_file = permanent_file
        self.tmp_directory = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp_directory.name)
        self.temp_file = self.tmp_path / self.permanent_file

    def open(self, mode: str = "r", buffering: int = -1, encoding=None, errors=None, newline=None):
        return self.temp_file.open(mode, buffering, encoding, errors, newline)

    def close(self):
        self.temp_file.close()
        shutil.move(self.temp_file, self.permanent_file)
        self.tmp_path.cleanup()


class Demultiplexer:
    def __init__(
        self,
        sample: Sample,
        barcode_df_callback: Callable,
        strategy: DemultiplexStrategy = PE_Decide_On_Start_Trim_Start_End,
        output_folder: Path = None,
        filter_func: Callable = None,
        maximal_error_rate: int = 0,
    ):
        self.barcode_df = barcode_df_callback()
        self.name = f"DM_{sample.name}"
        self.library_name = sample.name
        if output_folder is None:
            self.output_folder = Path("cache") / self.name
        else:
            self.output_folder = output_folder
        self.output_folder.mkdir(parents=True, exist_ok=True)
        # if isinstance(barcode_df_or_file, pd.DataFrame):
        #     self.barcodes = barcode_df_or_file
        # elif isinstance(barcode_df_or_file, str):
        #     self.barcodes = pd.read_csv(barcode_df_or_file, sep = '\t')
        # elif isinstance(barcode_df_or_file, Path):
        #     self.barcodes = pd.read_csv(str(barcode_df_or_file.resolve()), sep = '\t')
        # else:
        #     raise ValueError("No flanking barcodes supplied. please enter a barcode dataframe or a file path.")
        self.input_sample = sample
        self.maximal_error_rate = maximal_error_rate
        self.input_files = self.input_sample.get_aligner_input_filenames()
        self.strategy = strategy
        self.is_paired = self.input_sample.is_paired
        self.__initialize_decision_callbacks()

    def get_dependencies(self) -> List[Job]:
        return [
            ppg.ParameterInvariant(self.name + "_demultiplex_params", self.parameters),
            self.input_sample.prepare_input(),
        ]

    def parameters(self) -> List[str]:
        return [self.name, self.input_sample.name, self.strategy.__name__]

    def parameter_string(self) -> str:
        ret = f"{self.name}\nLibrary: {self.input_sample.name}\nDemultiplexStrategy: {self.strategy}\nBarcodeDastaFrame:\n"
        ret += self.barcode_df.to_string()
        return ret

    def __initialize_decision_callbacks(self):
        self.decision_callbacks = {}
        for key, df_1_row in self.barcode_df.iterrows():
            print(type(df_1_row))
            print(df_1_row)
            parameters = df_1_row.to_dict()
            print(parameters)
            self.decision_callbacks[key] = self.strategy(**parameters)

    def get_fastq_iterator(self):
        return get_fastq_iterator(self.is_paired)

    def _decide_on_barcode(self, fragment: Fragment):
        for key in self.decision_callbacks:
            accepted = self.decision_callbacks[key].match_and_trim(fragment)
            if accepted:
                return key, accepted
        return "discard", fragment

    def __write_fragment(self, fragment: Fragment, file_handles: List[FileIO]) -> None:
        for i, read in enumerate(fragment):
            file_handles[i].write(f"@{read.Name}\n{read.Sequence}\n+\n{read.Quality}\n".encode())

    def do_demultiplex(self, dependencies: List[Job] = []):
        deps = dependencies
        files_to_create = {}
        filenames = [sentinel]
        sample_names = [
            f"{self.library_name}_{key}" for key in list(self.decision_callbacks.keys())
        ] + [f"{library_name}_discarded"]
        for sample_name in sample_name:
            files_to_create[sample_name] = [
                self.output_folder / sample_name / f"{sample_name}_R1_.fastq"
            ]
            files_to_create[sample_name][0].parent.mkdir(parents=True, exist_ok=True)

        if self.is_paired:
            for sample_name in sample_names:
                files_to_create[sample_name].append(
                    self.output_folder / sample_name / f"{sample_name}_R2_.fastq"
                )

        def dump():
            # open a bunch of temporary files to write to
            temporary_files = {}
            sentinel = self.result_dir / "done.txt"
            with sentinel.open("w") as done:
                done.write(self.parameter_string())
                for sample_name in files_to_create:
                    temporary_files[sample_name] = [
                        TemporaryToPermanent(f).open("w") for f in files_to_create[sample_name]
                    ]
                # iterate over the input files and decide on each fragment, then write to temporary file
                read_iterator = self.get_fastq_iterator()
                for files_tuple in self.input_files:
                    for fragment in read_iterator(files_tuple):
                        key, accepted = self._decide_on_barcode(fragment)
                        sample_name = f"{self.library_name}_{key}"
                        self.__write_fragment(accepted, temporary_files[sample_name])
                # close open file handle
                for key in temporary_files:
                    for f in temporary_files[key]:
                        f.close()
                done.write("demultiplexing done")

        return ppg.MultiFileGeneratingJob(str(self.result_dir / "done.txt"), dump).depends_on(deps)
