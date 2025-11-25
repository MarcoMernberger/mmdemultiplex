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
import pandas as pd
# import pypipegraph as ppg
from collections import Counter
from pathlib import Path
from typing import Callable, List, Dict, Tuple, Optional, Any, Iterable
from mbf.align import Sample
from mbf.align.fastq import Straight
from .strategies import DemultiplexStrategy, PE_Decide_On_Start_Trim_Start_End
from .util import Fragment, get_fastq_iterator, TemporaryToPermanent, len_callback
from pypipegraph import (
    Job,
    MultiFileGeneratingJob,
    FileGeneratingJob,
    ParameterInvariant,
)
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
        filter_callable: Callable = None,
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
        self.filter_func = filter_callable
        if self.filter_func is None:
            self._decide_for_fragment = self._decide_on_barcode
        else:
            self._decide_for_fragment = self._decide_on_barcode_and_filter

    def __check_pairing(self):
        if self.is_paired and str(self.strategy.__class__).startswith("SE"):
            raise AttributeError("Strategy is not compatible with paired end reads")
        elif not self.is_paired and str(self.strategy.__class__).startswith("PE"):
            raise AttributeError("Strategy is not compatible with single end reads")

    def get_dependencies(self) -> List[Job]:
        deps = [
            ParameterInvariant(self.name + "_demultiplex_params", self.parameters),
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
            accepted = self.decision_callbacks[key].match_and_trim(fragment.copy())
            if accepted:
                return key, accepted
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

        def dump(
            filenames, self=self, sentinel=sentinel, files_to_create=files_to_create
        ):
            # open a bunch of temporary files to write to
            counter = Counter()
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
                        key, accepted = self._decide_for_fragment(fragment)
                        sample_name = f"{self.name}_{key}"
                        counter[sample_name] += 1
                        self.__write_fragment(accepted, temporary_files[sample_name])
                # close open file handle
                done.write("\nKey\tCount\n")
                for key in temporary_files:
                    done.write(f"{key}\t{counter[key]}\n")
                    for f in temporary_files[key]:
                        f.close()
                done.write("\ndemultiplexing done")

        return dump

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
            filenames,
            self._do_demultiplex_callable(files_to_create, sentinel),
            empty_ok=True,
        ).depends_on(deps)

    def _check_on_barcode(self, fragment: "Fragment") -> Dict[Any, dict]:
        """
        For each barcode combination key, call locate_count and get occurences:
   
             returns:
            { key: {
                "b1_in_r1": bool,
                "b1_in_r2": bool,
                "b2_in_r1": bool,
                "b2_in_r2": bool,
                "b1_fwd_r1": bool, "b1_rev_r1": bool,
                "b1_fwd_r2": bool, "b1_rev_r2": bool,
                "b2_fwd_r1": bool, "b2_rev_r1": bool,
                "b2_fwd_r2": bool, "b2_rev_r2": bool,
                "b1_present": bool,
                "b2_present": bool,
            }, ... }
        """
        result: Dict[Any, dict] = {}

        for key, cb in self.decision_callbacks.items():
            print(cb)
            (
                b1_fwd_r1, b2_fwd_r1,
                b1_rev_r1, b2_rev_r1,
                b1_fwd_r2, b2_fwd_r2,
                b1_rev_r2, b2_rev_r2,
            ) = cb.locate_count(fragment.copy())

            b1_in_r1 = b1_fwd_r1 or b1_rev_r1
            b1_in_r2 = b1_fwd_r2 or b1_rev_r2
            b2_in_r1 = b2_fwd_r1 or b2_rev_r1
            b2_in_r2 = b2_fwd_r2 or b2_rev_r2

            b1_present = b1_in_r1 or b1_in_r2
            b2_present = b2_in_r1 or b2_in_r2

            result[key] = {
                "b1_in_r1": b1_in_r1,
                "b1_in_r2": b1_in_r2,
                "b2_in_r1": b2_in_r1,
                "b2_in_r2": b2_in_r2,
                "b1_fwd_r1": b1_fwd_r1,
                "b1_rev_r1": b1_rev_r1,
                "b1_fwd_r2": b1_fwd_r2,
                "b1_rev_r2": b1_rev_r2,
                "b2_fwd_r1": b2_fwd_r1,
                "b2_rev_r1": b2_rev_r1,
                "b2_fwd_r2": b2_fwd_r2,
                "b2_rev_r2": b2_rev_r2,
                "b1_present": b1_present,
                "b2_present": b2_present,
            }

        return result

    def _check_adapters_callable(
        self,
        files_to_create: Iterable[Path],
        sentinel: Path,
    ) -> Callable:
        
        def dump(
            outfiles: Iterable[Path],
            self=self,
            sentinel=sentinel,
            files_to_create=files_to_create,
        ):
            # order as in check_adapters()
            outfile, outfile_matrix, outfile_pairs = outfiles

            # Counter per barcode combination
            adapter_counts: Dict[Any, Counter] = {
                key: Counter() for key in self.decision_callbacks
            }

            # NEW: Co-Occurence-Counts 
            barcode_pair_counts: Counter = Counter()

            global_total_fragments = 0
            global_unassigned_fragments = 0
            global_ambiguous_pairs = 0

            with sentinel.open("w") as done:
                done.write(self.parameter_string())

                read_iterator = self.get_fastq_iterator()

                for files_tuple in self.input_files:
                    for fragment in read_iterator(files_tuple):
                        global_total_fragments += 1

                        # get adapter info for each barcode combination
                        match_info = self._check_on_barcode(fragment)

                        # NEW: barcode-sequenz-sets per fragment         # <<< NEW
                        present_b1 = set()                               # <<< NEW
                        present_b2 = set()                               # <<< NEW
                        for key, info in match_info.items():             # <<< NEW
                            cb = self.decision_callbacks[key]            # <<< NEW
                            b1_seq = str(getattr(cb, "start_barcode", ""))  # <<< NEW
                            b2_seq = str(getattr(cb, "end_barcode", ""))    # <<< NEW
                            if info["b1_present"] and b1_seq:            # <<< NEW
                                present_b1.add(b1_seq)                   # <<< NEW
                            if info["b2_present"] and b2_seq:            # <<< NEW
                                present_b2.add(b2_seq)                   # <<< NEW
                        if present_b1 and present_b2:                    # <<< NEW
                            for b1_seq in present_b1:                    # <<< NEW
                                for b2_seq in present_b2:                # <<< NEW
                                    barcode_pair_counts[(b1_seq, b2_seq)] += 1  # <<< NEW

                        # keys with at least one barcode
                        keys_with_any = [
                            k for k, v in match_info.items()
                            if v["b1_present"] or v["b2_present"]
                        ]
                        # keys with both barcodes
                        keys_with_both = [
                            k for k, v in match_info.items()
                            if v["b1_present"] and v["b2_present"]
                        ]

                        if not keys_with_any:
                            # completely unassignable globally, plus none for all keys
                            global_unassigned_fragments += 1
                            for key in adapter_counts:
                                adapter_counts[key]["none"] += 1
                            continue

                        if len(keys_with_both) > 1:
                            # aAmbiguous, different barcode keys possible
                            global_ambiguous_pairs += 1
                            continue

                        if len(keys_with_both) == 1:
                            # unique pair (B1 & B2)
                            key = keys_with_both[0]
                            info = match_info[key]
                            stats = adapter_counts[key]

                            stats["any"] += 1
                            stats["both"] += 1

                            # positional infos
                            if info["b1_in_r1"]:
                                stats["b1_in_r1"] += 1
                            if info["b1_in_r2"]:
                                stats["b1_in_r2"] += 1
                            if info["b2_in_r1"]:
                                stats["b2_in_r1"] += 1
                            if info["b2_in_r2"]:
                                stats["b2_in_r2"] += 1

                            # strand info
                            if info["b1_fwd_r1"] or info["b1_fwd_r2"]:
                                stats["b1_forward"] += 1
                            if info["b1_rev_r1"] or info["b1_rev_r2"]:
                                stats["b1_reverse"] += 1
                            if info["b2_fwd_r1"] or info["b2_fwd_r2"]:
                                stats["b2_forward"] += 1
                            if info["b2_rev_r1"] or info["b2_rev_r2"]:
                                stats["b2_reverse"] += 1

                            expected = info["b1_in_r1"] and info["b2_in_r2"]
                            swapped = (info["b1_in_r2"] and info["b2_in_r1"]) and not expected

                            if expected:
                                stats["expected"] += 1
                            elif swapped:
                                stats["swapped"] += 1
                            else:
                                stats["other_orientation"] += 1

                        else:
                            # no keys_with_both, but there are keys_with_any
                            # So we have single barcodes present (b1_only / b2_only)
                            # fragments here are counted mutliple times!
                            for key in keys_with_any:
                                v = match_info[key]
                                stats = adapter_counts[key]
                                stats["any"] += 1

                                if v["b1_present"] and not v["b2_present"]:
                                    stats["b1_only"] += 1
                                elif v["b2_present"] and not v["b1_present"]:
                                    stats["b2_only"] += 1
                                else:
                                    # bot be true (should not happen)
                                    stats["other_orientation"] += 1

                # output metrics per barcode combination key
                detail_metrics = [
                    "any",
                    "both",
                    "b1_only",
                    "b2_only",
                    "none",
                    "b1_in_r1",
                    "b1_in_r2",
                    "b2_in_r1",
                    "b2_in_r2",
                    "b1_forward",
                    "b1_reverse",
                    "b2_forward",
                    "b2_reverse",
                    "expected",
                    "swapped",
                    "other_orientation",
                ]

                rows = []
                for key, stats in adapter_counts.items():
                    cb = self.decision_callbacks[key]
                    start_bc = getattr(cb, "start_barcode", "")
                    end_bc = getattr(cb, "end_barcode", "")

                    row = {
                        "sample_name": str(key),
                        "start_barcode": str(start_bc),
                        "end_barcode": str(end_bc),
                    }
                    for m in detail_metrics:
                        row[m] = int(stats.get(m, 0))
                    rows.append(row)

                df_detail = pd.DataFrame(rows)

                # force column order
                df_detail = df_detail[
                    ["sample_name", "start_barcode", "end_barcode"] + detail_metrics
                ]

                # detailed tsv
                df_detail.to_csv(outfile, sep="\t", index=False)

                # matrix file to indicate the number of combinations
                matrix_metrics = [
                    "any",
                    "both",
                    "b1_only",
                    "b2_only",
                    "none",
                    "expected",
                    "swapped",
                    "other_orientation",
                ]
                df_matrix = df_detail[
                    ["sample_name", "start_barcode", "end_barcode"] + matrix_metrics
                ]
                df_matrix.to_csv(outfile_matrix, sep="\t", index=False)

                # NEW: barcode-pair-co-occurence                          # <<< NEW
                pair_rows = []                                            # <<< NEW
                for (b1_seq, b2_seq), cnt in barcode_pair_counts.items(): # <<< NEW
                    pair_rows.append(                                     # <<< NEW
                        {                                                 # <<< NEW
                            "start_barcode": b1_seq,                      # <<< NEW
                            "end_barcode": b2_seq,                        # <<< NEW
                            "both": int(cnt),                             # <<< NEW
                        }                                                 # <<< NEW
                    )                                                     # <<< NEW
                df_pairs = pd.DataFrame(pair_rows)                        # <<< NEW
                if df_pairs.empty:                                        # <<< NEW
                    df_pairs = pd.DataFrame(                              # <<< NEW
                        columns=["start_barcode", "end_barcode", "both"]  # <<< NEW
                    )                                                     # <<< NEW
                df_pairs.to_csv(outfile_pairs, sep="\t", index=False)     # <<< NEW

                # global infos to sentinel
                done.write(
                    f"\nTotal fragments: {global_total_fragments}\n"
                    f"Unassigned fragments (no barcode for any key): {global_unassigned_fragments}\n"
                    f"Ambiguous fragments (multiple barcode pairs with both indices): {global_ambiguous_pairs}\n"
                )
                done.write("Check complete.")

        return dump

    def check_adapters(self):
        deps = self.get_dependencies()
        outfile = self.output_folder / f"{self.name}_adapter_check.tsv"
        outfile_matrix = self.output_folder / f"{self.name}_adapter_check_matrix_main.tsv"
        outfile_pairs = self.output_folder / f"{self.name}_adapter_check_barcode_pairs.tsv"

        files_to_create = [outfile, outfile_matrix, outfile_pairs]
        sentinel = self.output_folder / "adapter_check.txt"

        return MultiFileGeneratingJob(
            files_to_create,
            self._check_adapters_callable(files_to_create, sentinel),
            empty_ok=True,
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
                fastq_processor=Straight(),
                pairing=pairing,
                vid=None,
            )
        return raw_samples

    def get_samples(self):
        if not hasattr(self, "raw_samples"):
            self.raw_samples = self._make_samples()
        return self.raw_samples

    def _count_adapters_callable(self, k: int = 10, most_common: Optional[int] = None):
        """
        count_adapters counts the starting kmers of the reads in the input files.
        This is to check the actual adapters/barcodes used for demultiplexing.
        """

        def __count(output_file, self=self):
            read_iterator = self.get_fastq_iterator()
            counter = Counter()
            for files_tuple in self.input_files:
                for fragment in read_iterator(files_tuple):
                    for read in fragment:
                        counter[read.Sequence[:k]] += 1
            with output_file.open("w") as outp:
                for count in counter.most_common(most_common):
                    outp.write(f"{count[0]}\t{count[1]}\n")

        return __count

    def count_adapters(self, k: int = 10, most_common: Optional[int] = None):
        """
        count_adapters counts the starting kmers of the reads in the input files.
        This is to check the actual adapters/barcodes used for demultiplexing.
        """
        deps = self.get_dependencies()
        output_file = (
            self.output_folder / f"{self.input_sample.name}_barcode_counts.txt"
        )
        return FileGeneratingJob(
            output_file, self._count_adapters_callable(k, most_common), empty_ok=True
        ).depends_on(deps)

    def _divide_reads_callable(
        self,
        decision_callback=len_callback,
    ):
        """
        divide_reads counts the starting kmers of the reads in the input files.
        This is to check the actual adapters/barcodes used for demultiplexing.
        """

        def __divide(outfiles, self=self, decision_callback=decision_callback):
            outfiles[0].parent.mkdir(parents=True, exist_ok=True)
            read_iterator = self.get_fastq_iterator()
            temporary_files = [TemporaryToPermanent(f).open("w") for f in outfiles]
            for files_tuple in input_files:
                for fragment in read_iterator(files_tuple):
                    fragment = decision_callback(fragment)
                    self.__write_fragment(fragment, temporary_files)
            # close open file handle
            for f in temporary_files:
                f.close()

        return __divide

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
        outfiles = (
            new_output_folder / input_files[0][0].name,
            new_output_folder / input_files[0][1].name,
        )
        return (
            MultiFileGeneratingJob(
                outfiles, self._divide_reads_callable(decision_callback), empty_ok=True
            )
            .depends_on(dependencies)
            .depends_on(deps)
        )
