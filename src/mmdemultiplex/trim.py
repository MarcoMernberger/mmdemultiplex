#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""trim.py: Contains a class for trimming reads."""

import mbf
import pypipegraph as ppg
from collections import Counter
from pathlib import Path
from typing import Callable, List, Dict, Tuple, Optional
from mbf.align import Sample
from .strategies import (
    DemultiplexStrategy,
    PE_Decide_On_Start_Trim_Start_End,
    SE_Trim_On_Start_Trim_After_X_BP,
)
from .util import Fragment, get_fastq_iterator, TemporaryToPermanent, len_callback
from pypipegraph import Job
from .samples import FASTQsFromJobSelect

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class Trimmer:
    def __init__(
        self,
        sample: Sample,
        strategy: DemultiplexStrategy = SE_Trim_On_Start_Trim_After_X_BP,
        output_folder: Optional[Path] = None,
        prefix: str = "",
        filter_callable: Callable = None,
        parameters: Dict[str, str] = {},
    ):
        self.name = f"{prefix}{sample.name}"
        self.library_name = sample.name
        if output_folder is None:
            self.output_folder = Path("cache") / self.name
        else:
            self.output_folder = output_folder / self.name
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.input_sample = sample
        self.input_files = self.input_sample.get_aligner_input_filenames()
        self.is_paired = self.input_sample.is_paired
        if isinstance(self.input_files, Tuple):
            self.input_files = [self.input_files]
        elif isinstance(self.input_files, List) and isinstance(
            self.input_files[0], str
        ):
            self.input_files = [(x,) for x in self.input_files]
        self.trim_strategy = strategy
        self.__check_pairing()
        self.filter_func = filter_callable
        if self.filter_func is None:
            self._trim_fragment = self._trim_only
        else:
            self._trim_fragment = self._trim_and_filter
        self.strategy_params = parameters
        self.__initialize_strategy()

    def __check_pairing(self):
        if self.is_paired and str(self.trim_strategy.__class__).startswith("SE"):
            raise AttributeError("Strategy is not compatible with paired end reads")
        elif not self.is_paired and str(self.trim_strategy.__class__).startswith("PE"):
            raise AttributeError("Strategy is not compatible with single end reads")

    def __initialize_strategy(self) -> None:
        self.strategy = self.trim_strategy(**self.strategy_params)

    def get_dependencies(self) -> List[Job]:
        deps = [
            ppg.ParameterInvariant(self.name + "_trim_params", self.parameters),
        ]
        if hasattr(self.input_sample, "prepare_input"):
            deps.append(self.input_sample.prepare_input())
        if hasattr(self.input_sample, "dependencies"):
            deps.extend(self.input_sample.dependencies)
        return deps

    @property
    def parameters(self) -> List[str]:
        return [self.name, self.input_sample.name, self.trim_strategy.__name__]

    def parameter_string(self) -> str:
        ret = f"{self.name}\nLibrary: {self.input_sample.name}\nStrategy: {self.trim_strategy}\n"
        return ret

    def get_fastq_iterator(self) -> Callable:
        return get_fastq_iterator(self.is_paired)

    def _trim_only(self, fragment: Fragment) -> Tuple[str, Fragment]:
        trimmed = self.strategy.match_and_trim(fragment)
        if trimmed:
            return "trimmed", trimmed
        return "discarded", fragment

    def _trim_and_filter(self, fragment: Fragment) -> Tuple[str, Fragment]:
        if not self.filter_func(fragment):
            return "discarded", fragment
        return self._trim_only(fragment)

    def __write_fragment(
        self, fragment: Fragment, file_handles: List[TemporaryToPermanent]
    ) -> None:
        for i, read in enumerate(fragment):
            file_handles[i].write(f"@{read.Name}\n{read.Sequence}\n+\n{read.Quality}\n")

    def trim(self):
        deps = self.get_dependencies()
        files_to_create = self.get_files_to_create()
        sentinel = self.output_folder / "done.txt"
        filenames = [sentinel] + [
            filename for files in files_to_create.values() for filename in files
        ]
        for sample_name in files_to_create:
            files_to_create[sample_name][0].parent.mkdir(parents=True, exist_ok=True)

        def dump():
            # open a bunch of temporary files to write to
            with sentinel.open("w") as done:
                temporary_files = {}
                done.write(self.parameter_string())
                for sample_name in files_to_create:
                    temporary_files[sample_name] = [
                        TemporaryToPermanent(f).open("w")
                        for f in files_to_create[sample_name]
                    ]
                # iterate over the input files and decide on each fragment, then write to temporary file
                read_iterator = self.get_fastq_iterator()
                for files_tuple in self.input_files:
                    for fragment in read_iterator(files_tuple):
                        key, trimmed_fragment = self._trim_fragment(fragment)
                        sample_name = f"{self.name}_{key}"
                        self.__write_fragment(
                            trimmed_fragment, temporary_files[sample_name]
                        )
                # close open file handle
                for key in temporary_files:
                    for f in temporary_files[key]:
                        f.close()
                done.write("\ntrimming done")

        return ppg.MultiFileGeneratingJob(filenames, dump, empty_ok=True).depends_on(
            deps
        )

    def get_files_to_create(self):
        files_to_create = {}
        sample_names = [f"{self.name}_trimmed"] + [f"{self.name}_discarded"]
        for sample_name in sample_names:
            files_to_create[sample_name] = [
                self.output_folder / sample_name / f"{sample_name}_R1_.fastq"
            ]
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
        sample_name = f"{self.name}_trimmed"
        raw_samples[sample_name] = Sample(
            sample_name,
            input_strategy=FASTQsFromJobSelect(sample_name, self.trim()),
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
