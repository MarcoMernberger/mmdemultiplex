#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains utility functions for the demultiplexer package."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from dataclasses import dataclass
import pandas as pd
import tempfile
import shutil
import pypipegraph as ppg
import collections
import mbf_align

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


rev_comp_table = maketrans(b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR")


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


@dataclass
class Read:
    Name: str
    Sequence: str
    Quality: str


class Fragment:
    def __init__(self, *reads: Read):
        self.reads = reads
        self.Read1 = self.reads[0]
        if len(reads) == 2:
            self.Read2 = self.reads[1]

    def __iter__(self):
        for read in self.reads:
            yield read


class TemporaryToPermanent:
    def __init__(self, permanent_file: Path):
        self.permanent_file = permanent_file
        self.tmp_directory = tempfile.TemporaryDirectory(dir=permanent_file.parent)
        self.tmp_path = Path(self.tmp_directory.name)
        self.temp_file = self.tmp_path / self.permanent_file.relative_to(self.permanent_file.root)
        self.temp_file.parent.mkdir(exist_ok=True, parents=True)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is None:
            self.close()
        else:
            self.file_handle.close()
            self.tmp_directory.cleanup()

    def open(self, *args, **kwargs):
        self.file_handle = self.temp_file.open(*args, **kwargs)
        return self

    def close(self):
        self.file_handle.close()
        shutil.move(self.temp_file, self.permanent_file)
        self.tmp_directory.cleanup()

    def write(self, *args, **kwargs):
        self.file_handle.write(*args, **kwargs)

    @property
    def closed(self):
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


def iterate_fastq(fn, reverse_reads):
    op = mbf_align._common.BlockedFileAdaptor(fn)
    while True:
        try:
            name = op.readline()[1:-1].decode()
            seq = op.readline()[:-1].decode()
            qual = op.readline()[:-1].decode()
            if reverse_reads:
                seq = seq[::-1].translate(rev_comp_table)
                qual = qual[::-1]
            yield Read(name, seq, qual)
        except StopIteration:
            break


def get_fastq_iterator(paired):
    fastq_iterator = iterate_fastq

    def _iterreads_paired_end(tuple_of_files: Tuple[Path]):
        for reads in zip(
            fastq_iterator(str(tuple_of_files[0]), reverse_reads=False),
            fastq_iterator(str(tuple_of_files[1]), reverse_reads=False),
        ):
            yield Fragment(*reads)

    def _iterreads_single_end(filetuple):
        for read in fastq_iterator(str(filetuple[0]), reverse_reads=False):
            yield Fragment(read)

    if paired:
        return _iterreads_paired_end
    else:
        return _iterreads_single_end
