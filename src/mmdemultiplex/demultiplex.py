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
import mbf
import pypipegraph as ppg
from collections import Counter
from pathlib import Path
from typing import Callable, List, Dict, Tuple, Optional
from mbf.align import Sample
from .strategies import DemultiplexStrategy, PE_Decide_On_Start_Trim_Start_End
from .util import Fragment, get_fastq_iterator, TemporaryToPermanent, len_callback
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
        output_folder: Optional[Path] = None,
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
        self.input_files = self.input_sample.get_aligner_input_filenames()
        self.is_paired = self.input_sample.is_paired
        if isinstance(self.input_files, Tuple):
            self.input_files = [self.input_files]
        elif isinstance(self.input_files, List) and isinstance(
            self.input_files[0], str
        ):
            self.input_files = [(x,) for x in self.input_files]
        self.strategy = strategy
        self.__initialize_decision_callbacks()
        self.__check_pairing()

    def __check_pairing(self):
        if self.is_paired and str(self.strategy.__class__).startswith("SE"):
            raise AttributeError("Strategy is not compatible with paired end reads")
        elif not self.is_paired and str(self.strategy.__class__).startswith("PE"):
            raise AttributeError("Strategy is not compatible with single end reads")

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
                        key, accepted = self._decide_on_barcode(fragment)
                        print(fragment.Read1.Name, key, accepted)
                        sample_name = f"{self.name}_{key}"
                        self.__write_fragment(accepted, temporary_files[sample_name])
                # close open file handle
                for key in temporary_files:
                    for f in temporary_files[key]:
                        f.close()
                done.write("\ndemultiplexing done")

        return ppg.MultiFileGeneratingJob(filenames, dump, empty_ok=True).depends_on(
            deps
        )

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

    def count_adapters(self, k: int = 10, most_common: Optional[int] = None):
        """
        count_adapters counts the starting kmers of the reads in the input files.
        This is to check the actual adapters/barcodes used for demultiplexing.
        """
        deps = self.get_dependencies()
        output_file = (
            self.output_folder / f"{self.input_sample.name}_barcode_counts.txt"
        )

        def count():
            read_iterator = self.get_fastq_iterator()
            counter = Counter()
            for files_tuple in self.input_files:
                for fragment in read_iterator(files_tuple):
                    for read in fragment:
                        counter[read.Sequence[:k]] += 1
            with output_file.open("w") as outp:
                for count in counter.most_common(most_common):
                    outp.write(f"{count[0]}\t{count[1]}\n")

        return ppg.FileGeneratingJob(output_file, count, empty_ok=True).depends_on(deps)

    def divide_reads(
        self,
        input_files: List[Tuple[Path, Path]],
        decision_callback=len_callback,
        new_output_folder=Path("results/demultiplexed/divided"),
        dependencies=[],
    ):
        """
        divide_reads counts the starting kmers of the reads in the input files.
        This is to check the actual adapters/barcodes used for demultiplexing.
        """
        deps = self.get_dependencies()
        print(new_output_folder)
        new_output_folder.mkdir(parents=True, exist_ok=True)
        print(new_output_folder.exists())
        outfiles = (
            new_output_folder / input_files[0][0].name,
            new_output_folder / input_files[0][1].name,
        )

        def divide():
            read_iterator = self.get_fastq_iterator()
            temporary_files = [TemporaryToPermanent(f).open("w") for f in outfiles]
            for files_tuple in input_files:
                for fragment in read_iterator(files_tuple):
                    fragment = decision_callback(fragment)
                    self.__write_fragment(fragment, temporary_files)
            # close open file handle
            for f in temporary_files:
                f.close()

        return (
            ppg.MultiFileGeneratingJob(outfiles, divide, empty_ok=True)
            .depends_on(dependencies)
            .depends_on(deps)
        )
