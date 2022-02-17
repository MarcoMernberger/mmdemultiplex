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

import mbf_align
import pypipegraph as ppg
from pathlib import Path
from typing import Callable, List, Dict, Tuple
from mbf_align import Sample
from .strategies import PE_Decide_On_Start_Trim_Start_End, DemultiplexStrategy
from .util import Fragment, get_fastq_iterator, TemporaryToPermanent
from pypipegraph import Job
from .samples import FASTQsFromJobSelect

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class Demultiplexer:
    def __init__(
        self,
        sample: Sample,
        barcode_df_callback: Callable,
        strategy: DemultiplexStrategy = PE_Decide_On_Start_Trim_Start_End,
        output_folder: Path = None,
        maximal_error_rate: int = 0,
        prefix: str = "",
    ):
        self.barcode_df = barcode_df_callback()
        self.name = f"{prefix}{sample.name}"
        self.library_name = sample.name
        if output_folder is None:
            self.output_folder = Path("cache") / self.name
        else:
            self.output_folder = output_folder / self.name
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.input_sample = sample
        self.maximal_error_rate = maximal_error_rate
        self.input_files = self.input_sample.get_aligner_input_filenames()
        self.is_paired = self.input_sample.is_paired
        if isinstance(self.input_files, Tuple):
            self.input_files = [self.input_files]
        elif isinstance(self.input_files, List) and isinstance(self.input_files[0], str):
            self.input_files = [(x,) for x in self.input_files]
        self.strategy = strategy
        self.__initialize_decision_callbacks()

    def get_dependencies(self) -> List[Job]:
        deps = [
            ppg.ParameterInvariant(self.name + "_demultiplex_params", self.parameters),
        ]
        if hasattr(self.input_sample, "prepare_input"):
            deps.append(self.input_sample.prepare_input())
        if hasattr(self.input_sample, "dependencies"):
            deps.extend(self.input_sample.dependencies)
        return deps

    @property
    def parameters(self) -> List[str]:
        return [self.name, self.input_sample.name, self.strategy.__name__]

    def parameter_string(self) -> str:
        ret = f"{self.name}\nLibrary: {self.input_sample.name}\nDemultiplexStrategy: {self.strategy}\nBarcodeDastaFrame:\n"
        ret += self.barcode_df.to_string()
        return ret

    def __initialize_decision_callbacks(self) -> None:
        self.decision_callbacks = {}
        for key, df_1_row in self.barcode_df.iterrows():
            parameters = df_1_row.to_dict()
            self.decision_callbacks[key] = self.strategy(**parameters)

    def get_fastq_iterator(self) -> Callable:
        return get_fastq_iterator(self.is_paired)

    def _decide_on_barcode(self, fragment: Fragment):
        for key in self.decision_callbacks:
            accepted = self.decision_callbacks[key].match_and_trim(fragment)
            if accepted:
                return key, accepted
        return "discarded", fragment

    def __write_fragment(
        self, fragment: Fragment, file_handles: List[TemporaryToPermanent]
    ) -> None:
        for i, read in enumerate(fragment):
            file_handles[i].write(f"@{read.Name}\n{read.Sequence}\n+\n{read.Quality}\n")

    def do_demultiplex(self):
        deps = self.get_dependencies()
        files_to_create = {}
        sample_names = [f"{self.name}_{key}" for key in list(self.decision_callbacks.keys())] + [
            f"{self.name}_discarded"
        ]
        for sample_name in sample_names:
            files_to_create[sample_name] = [
                self.output_folder / sample_name / f"{sample_name}_R1_.fastq"
            ]
            files_to_create[sample_name][0].parent.mkdir(parents=True, exist_ok=True)
        if self.is_paired:
            for sample_name in sample_names:
                files_to_create[sample_name].append(
                    self.output_folder / sample_name / f"{sample_name}_R2_.fastq"
                )
        sentinel = self.output_folder / "done.txt"
        filenames = [sentinel] + [
            filename for files in files_to_create.values() for filename in files
        ]

        def dump():
            # open a bunch of temporary files to write to
            with sentinel.open("w") as done:
                temporary_files = {}
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
                        sample_name = f"{self.name}_{key}"
                        self.__write_fragment(accepted, temporary_files[sample_name])
                # close open file handle
                for key in temporary_files:
                    for f in temporary_files[key]:
                        f.close()
                done.write("\ndemultiplexing done")

        return ppg.MultiFileGeneratingJob(filenames, dump, empty_ok=True).depends_on(deps)

    def _make_samples(self) -> Dict[str, Sample]:
        raw_samples = {}
        pairing = "single"
        if self.is_paired:
            pairing = "paired"
        for key in self.decision_callbacks.keys():
            sample_name = f"{self.name}_{key}"
            raw_samples[sample_name] = Sample(
                sample_name,
                input_strategy=FASTQsFromJobSelect(sample_name, self.do_demultiplex()),
                reverse_reads=False,
                fastq_processor=mbf_align.fastq2.Straight(),
                pairing=pairing,
                vid=None,
            )
        return raw_samples

    def get_samples(self):
        if not hasattr(self, "raw_samples"):
            self.raw_samples = self._make_samples()
        return self.raw_samples
