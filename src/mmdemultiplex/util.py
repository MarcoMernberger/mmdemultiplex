#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains utility functions for the demultiplexer package."""

from pathlib import Path
from typing import Optional, Callable, Tuple
from pandas import DataFrame
from dataclasses import dataclass, replace
import tempfile
import shutil
import collections
import gzip, bz2
import re

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


@dataclass
class Read:
    """Data class for sequencing reads"""

    Name: str
    Sequence: str
    Quality: str

    def reverse():
        return Read(
            self.Name, reverse_complement(self.Sequence[::-1]), self.Quality[::-1]
        )

    def __str__(self):
        return f"{self.Name}\n{self.Sequence}\n+\n{self.Quality}\n"


class Fragment:
    """Data class for single-end and paired-end Reads/Fragments."""

    def __init__(self, *reads: Read):
        self.Read1 = reads[0]
        if len(reads) == 2:
            self.Read2 = reads[1]

    @property
    def is_paired(self):
        return hasattr(self, "Read2")

    @property
    def reads(self):
        if self.is_paired:
            return [self.Read1, self.Read2]
        else:
            return [
                self.Read1,
            ]

    def __iter__(self):
        for read in self.reads:
            yield read

    def copy(self):
        return Fragment(replace(self.Read1), replace(self.Read2))

    def __str__(self):
        return f"{self.Read1}\n{self.Read2}\n"


class TemporaryToPermanent:
    def __init__(self, permanent_file: Path):
        self.permanent_file = permanent_file

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback) -> None:
        if exception_type is None:
            self.close()
        else:
            self.file_handle.close()
            self.tmp_directory.cleanup()

    def open(self, *args, **kwargs):
        self.tmp_directory = tempfile.TemporaryDirectory(dir=self.permanent_file.parent)
        self.tmp_path = Path(self.tmp_directory.name)
        self.temp_file = self.tmp_path / self.permanent_file.relative_to(
            self.permanent_file.root
        )
        self.temp_file.parent.mkdir(exist_ok=True, parents=True)
        self.file_handle = self.temp_file.open(*args, **kwargs)
        return self

    def close(self):
        self.file_handle.close()
        shutil.move(self.temp_file, self.permanent_file)
        delattr(self, "file_handle")
        self.tmp_directory.cleanup()
        delattr(self, "tmp_path")
        # delattr(self, "file_handle")

    def write(self, *args, **kwargs):
        self.file_handle.write(*args, **kwargs)

    @property
    def closed(self) -> bool:
        if hasattr(self, "file_handle"):
            return self.file_handle.closed
        return True


def reverse_complement(sequence: str) -> str:
    """
    reverse_complement retuzrns the reverse complement of given sequence.

    Parameters
    ----------
    sequence : str
        Input sequence.

    Returns
    -------
    str
        Reverse complement of input sequence.
    """
    return sequence[::-1].translate(rev_comp_table)


def get_df_callable_for_demultiplexer(
    df_in: DataFrame,
    fw_col_name: str,
    rv_col_name: str,
    sample_col_name: str,
    trim_start_col_name: Optional[str] = None,
    trim_end_col_name: Optional[str] = None,
) -> DataFrame:
    def call():
        df = df_in.copy()
        df[fw_col_name].fillna("", inplace=True)
        df[rv_col_name].fillna("", inplace=True)
        df[fw_col_name] = df[fw_col_name].str.strip()
        df[fw_col_name] = df[fw_col_name].str.upper()
        df[rv_col_name] = df[rv_col_name].str.strip()
        df[rv_col_name] = df[rv_col_name].str.upper()
        whitespace = re.compile(r"\s+")
        assert len((df[fw_col_name] + df[rv_col_name]).unique()) == len(
            df
        )  # check if the barcodes are unique
        df["key"] = df[sample_col_name].str.replace(whitespace, "_", regex=True)
        df = df.set_index("key")

        if trim_start_col_name is None and trim_end_col_name is None:
            df = df.rename(
                columns={
                    fw_col_name: "start_barcode",
                    rv_col_name: "end_barcode",
                }
            )
            return df[["start_barcode", "end_barcode"]]

        else:
            df = df.rename(
                columns={
                    fw_col_name: "start_barcode",
                    rv_col_name: "end_barcode",
                    trim_start_col_name: "trim_after_start",
                    trim_end_col_name: "trim_before_end",
                }
            )
            return df[
                ["start_barcode", "end_barcode", "trim_after_start", "trim_before_end"]
            ]

    return call


def _open_auto(filename: str):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rb")
    if filename.endswith(".bz2"):
        return bz2.open(filename, "rb")
    return open(filename, "rb", buffering=4 * 1024 * 1024)  # großer Buffer


def iterate_fastq(filename: str, reverse_reads: bool) -> Read:
    # op = mbf.align._common.BlockedFileAdaptor(filename)
    op = _open_auto(filename)
    while True:
        try:
            name = op.readline()[1:-1].decode()
            seq = op.readline()[:-1].decode()
            op.readline()
            qual = op.readline()[:-1].decode()
            if reverse_reads:
                seq = seq[::-1].translate(rev_comp_table)
                qual = qual[::-1]
            yield Read(name, seq, qual)
        except StopIteration:
            break


def get_fastq_iterator(paired) -> Callable:
    fastq_iterator = iterate_fastq

    def _iterreads_paired_end(tuple_of_files: Tuple[Path, Path]) -> Fragment:
        for reads in zip(
            fastq_iterator(str(tuple_of_files[0]), reverse_reads=False),
            fastq_iterator(str(tuple_of_files[1]), reverse_reads=False),
        ):
            yield Fragment(*reads)

    def _iterreads_single_end(filetuple) -> Fragment:
        for read in fastq_iterator(str(filetuple[0]), reverse_reads=False):
            yield Fragment(read)

    if paired:
        return _iterreads_paired_end
    else:
        return _iterreads_single_end


def len_callback(fragment):
    if len(fragment.Read1.Sequence) < len(fragment.Read2.Sequence):
        fragment = Fragment(fragment.Read2, fragment.Read1)
    return fragment


def create_fragment(sequences, qualities=None, names=None):
    "make sequences into fragments with default names and qualities"
    if qualities is None:
        qualities = ["F" * len(sequence) for sequence in sequences]
    if names is None:
        names = ["read" + str(i) for i in range(len(sequences))]
    reads = [
        Read(name, sequence, quality)
        for name, sequence, quality in zip(names, sequences, qualities)
    ]
    fragment = Fragment(*reads)
    return fragment
