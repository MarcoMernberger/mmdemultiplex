#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""samples.py: Contains ...."""

from builtins import NotImplementedError
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from mbf_align import fastq2
from mbf_align.strategies import _FASTQsBase
from pypipegraph import Job
import pandas as pd
import pypipegraph as ppg

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class DemultiplexInputSample:
    def __init__(
        self,
        sample_name,
        input_strategy,
        reverse_reads,
        fastq_processor=fastq2.Straight(),
        pairing="paired",
    ):
        self.name = sample_name
        self.input_strategy = input_strategy
        self.reverse_reads = reverse_reads
        self.fastq_processor = fastq_processor
        self.pairing = pairing
        self.dependencies = []
        self.is_paired = self.pairing in ["paired"]
        self._set_input_files()

    def _set_input_files(self):
        input_pairs = self.input_strategy()
        any_r2 = any([len(x) > 1 for x in input_pairs])
        if self.pairing == "single":
            if any_r2:
                raise ValueError(
                    f"{self.name}: paired end lane defined as single end - you need to change the pairing parameter. Available is ('paired', 'single', 'paired_as_single'.)"
                )
            input_filenames = [(str(f[0])) for f in input_pairs]
        elif self.pairing == "paired_as_single":
            input_filenames = [(str(f[0])) for f in input_pairs]
        elif self.pairing == "paired":
            if not any_r2:
                raise ValueError(
                    f"Paired end lane, but no R2 reads found. Found files: {input_pairs}"
                )
            input_filenames = [(str(f[0]), str(f[1])) for f in input_pairs]
        else:
            raise ValueError("unknown pairing")  # pragma: no cover
        if self.pairing == "paired":
            flat_input_filenames = [f for fl in input_pairs for f in fl]
        else:
            flat_input_filenames = input_filenames
        if hasattr(self.input_strategy, "dependencies"):
            deps = self.input_strategy.dependencies
        else:
            deps = [ppg.FileChecksumInvariant(f) for f in flat_input_filenames]
        self.input_filenames = input_filenames
        self.dependencies.extend(deps)

    def get_aligner_input_filenames(self):
        return self.input_filenames


class FASTQsFromJobSelect(_FASTQsBase):
    def __init__(self, sample_name: str, job: Job):
        self.dependencies = job
        self.job = job
        self.sample_name = sample_name

    def __call__(self):
        return self._parse_filenames()

    def _parse_filenames(self):
        correct_files = []
        for filename in self.job.filenames:
            if self.sample_name in filename:
                correct_files.append(Path(filename))
        return super()._parse_filenames(correct_files)

    def __str__(self):
        return f"FASTQsFromJobSelect({self.sample_name} from {self.job})"  # pragma: no cover
