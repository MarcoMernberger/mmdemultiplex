#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
filterfastq.py: This module has a class to filter fastqs that is used for
demultiplexing fastqs. It offers the option to find barcodes within the sequences,
thereby removing trailing sequences. In addition, it can trim fixed length
sequences after the detected barcode. Can be applied to single end paired end
sequencing. Demultiplexing can Be defined on Forward barcodes, reverse barcodes
and a combination of both (in case of paired end reads.
"""
import mbf

# import pypipegraph as ppg
import pandas as pd
from collections import Counter
from pathlib import Path
from typing import Callable, List, Dict, Tuple, Optional, Union
from mbf.align import Sample
from .strategies import DemultiplexStrategy, PE_Decide_On_Start_Trim_Start_End
from .util import Fragment, get_fastq_iterator, TemporaryToPermanent, len_callback
from pypipegraph import Job, MultiFileGeneratingJob, ParameterInvariant
from .samples import FASTQsFromJobSelect


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def decision_callback_init(input_dataframe: Path) -> Callable:
    """
    This function generates a bunch of decision_callbacks and returns a dictionary
    of those to be used in the FasqDemultiplexer below.

    Parameters
    ----------
    input_dataframe : Path
        The path to a dataframe.

    Returns
    -------
    Dict[str, Callable]
        A dictionary of callables that decide on acceptance of a fragment.
    """

    def __callback_init():
        def decide(sequenced_sgRNAs_sensors):
            def __callback(fragment: Fragment):
                """
                This function decides on acceptance of a fragment based on the sequenced sensors and sgRNAs.
                """
                accepted = (
                    fragment.Read1.Sequence,
                    fragment.Read2.Sequence,
                ) in sequenced_sgRNAs_sensors
                return accepted, fragment

            return __callback

        decision_callbacks = {}
        df = pd.read_csv(input_dataframe, sep="\t")
        for key, df_sub in df.groupby("ID"):
            sequenced_sgRNAs_sensors = dict(
                [
                    ((sg, sensor), key)
                    for sg, sensor in zip(
                        df_sub["sequenced sgRNA"], df_sub["sequenced Sensor"]
                    )
                ]
            )
            decision_callbacks[key] = decide(sequenced_sgRNAs_sensors)
        return decision_callbacks

    return __callback_init


class FastqDemultiplexer:

    def __init__(
        self,
        decision_callback_init: Callable,
        fastq_input: Union[Path, Tuple[Path], Sample],
        output_folder: Optional[Path] = None,
        prefix: str = "",
        filter_callable: Callable = None,
        dependencies: List[Job] = [],
        name: Optional[str] = None,
    ):
        self.name = name
        self.decision_callback_init = decision_callback_init
        self.fastq_input = fastq_input
        self.prefix = prefix
        self.__init_files()
        if output_folder is None:
            self.output_folder = Path("cache") / self.name
        else:
            self.output_folder = output_folder / self.name
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.__initialize_decision_callbacks()
        # self.__check_pairing()
        self.filter_func = filter_callable
        if self.filter_func is None:
            self._decide_for_fragment = self._decide_on_barcode
        else:
            self._decide_for_fragment = self._decide_on_barcode_and_filter
        self.parameters = [fastq_input, output_folder, prefix]
        self.dependencies = dependencies

    def __init_files(self):
        if isinstance(self.fastq_input, Sample):
            self.input_sample = self.fastq_input
            self.input_files = self.input_sample.get_aligner_input_filenames()
            self.is_paired = self.input_sample.is_paired
            name = f"{self.prefix}{self.input_sample.name}"
            self.library_name = self.input_sample.name
            if isinstance(self.input_files, Tuple):
                self.input_files = [self.input_files]
            elif isinstance(self.input_files, List) and isinstance(
                self.input_files[0], str
            ):
                self.input_files = [(x,) for x in self.input_files]
        elif isinstance(self.fastq_input, Path):
            self.input_sample = None
            self.input_files = [self.fastq_input]
            self.is_paired = False
            self.library_name = self.fastq_input.stem
            name = f"{self.prefix}{self.library_name}"
        elif isinstance(self.fastq_input, Tuple):
            self.input_sample = None
            self.input_files = [self.fastq_input]
            self.is_paired = len(self.fastq_input) == 2
            self.library_name = self.fastq_input[0].stem
            name = f"{self.prefix}{self.library_name}"
        else:
            raise NotImplementedError(
                f"Input type {type(self.fastq_input)} not supported. "
                "Please provide a Sample, Path or Tuple of Paths."
            )
        if self.name is None:
            self.name = name

    def get_dependencies(self) -> List[Job]:
        deps = [
            ParameterInvariant(self.name + "_demultiplex_params", self.parameters),
        ] + self.dependencies

        if hasattr(self.input_sample, "prepare_input"):
            deps.append(self.input_sample.prepare_input())
        if hasattr(self.input_sample, "dependencies"):
            deps.extend(self.input_sample.dependencies)
        return deps

    def __initialize_decision_callbacks(self) -> None:
        self.decision_callbacks = self.decision_callback_init()

    def get_fastq_iterator(self) -> Callable:
        return get_fastq_iterator(self.is_paired)

    def _decide_on_barcode(self, fragment: Fragment):
        for key in self.decision_callbacks:
            accepted, fragment = self.decision_callbacks[key](fragment)
            if accepted:
                return key, fragment
        return "discarded", fragment

    def _decide_on_barcode_and_filter(self, fragment: Fragment):
        if not self.filter_func(fragment):
            return "discarded", fragment
        return self._decide_on_barcode(fragment)

    def __write_fragment(
        self, fragment: Fragment, file_handles: List[TemporaryToPermanent]
    ) -> None:
        for i, read in enumerate(fragment):
            file_handles[i].write(f"@{read.Name}\n{read.Sequence}\n+\n{read.Quality}\n")

    def _do_demultiplex_callable(
        self, files_to_create: Dict[str, Path], sentinel: Path
    ):

        def __dump(filenames, self=self, files_to_create=files_to_create. sentinel=sentinel):
            # open a bunch of temporary files to write to
            with sentinel.open("w") as done:
                temporary_files = {}
                # done.write(self.parameter_string())
                for sample_name in files_to_create:
                    temporary_files[sample_name] = [
                        TemporaryToPermanent(f).open("w")
                        for f in files_to_create[sample_name]
                    ]
                # iterate over the input files and decide on each fragment, then write to temporary file
                read_iterator = self.get_fastq_iterator()
                for files_tuple in self.input_files:
                    for fragment in read_iterator(files_tuple):
                        key, accepted = self._decide_for_fragment(fragment)
                        sample_name = f"{self.name}_{key}"
                        self.__write_fragment(accepted, temporary_files[sample_name])
                # close open file handle
                for key in temporary_files:
                    for f in temporary_files[key]:
                        f.close()
                done.write("\ndemultiplexing done")

        return __dump

    def do_demultiplex(self):
        deps = self.get_dependencies()
        files_to_create = {}
        sample_names = [
            f"{self.name}_{key}" for key in list(self.decision_callbacks.keys())
        ] + [f"{self.name}_discarded"]
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

        return MultiFileGeneratingJob(
            filenames, self._do_demultiplex_callable(files_to_create, sentinel), empty_ok=True
        ).depends_on(deps)

    def get_files_to_create(self):
        files_to_create = {}
        sample_names = [
            f"{self.name}_{key}" for key in list(self.decision_callbacks.keys())
        ] + [f"{self.name}_discarded"]
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
        return files_to_create

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
                fastq_processor=mbf.align.fastq2.Straight(),
                pairing=pairing,
                vid=None,
            )
        return raw_samples

    def get_samples(self):
        if not hasattr(self, "raw_samples"):
            self.raw_samples = self._make_samples()
        return self.raw_samples
