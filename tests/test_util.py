# -*- coding: utf-8 -*-

import pytest
import pandas as pd
from mmdemultiplex import Demultiplexer, DemultiplexStrategy, PE_Decide_On_Start_Trim_Start_End
from mmdemultiplex.util import Read, Fragment, reverse_complement, TemporaryToPermanent
from pathlib import Path
from unittest.mock import patch
from pypipegraph import Job
from conftest import MockBlockedFileAdapter, iterate_fastq, R1, R2

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_read():
    read = Read("rname", "rsequence", "rquality")
    assert hasattr(read, Name)
    assert hasattr(read, Sequence)
    assert hasattr(read, Quality)
    assert read.Name == "rname"
    assert read.Sequence == "rsequence"
    assert read.Quality == "rquality"


def test_fragment():
    read1 = Read("r1", "seq1", "qual1")
    read2 = Read("r2", "seq2", "qual2")
    read3 = Read("r3", "seq3", "qual3")
    fragment_pe = Fragment(read1, read2)
    fragment_se = Fragment(read3)
    for fragment in [fragment_se, fragment_pe]:
        assert hasattr(fragment, "Read1")
        assert isinstance(fragment.Read1, Read)
        assert hasattr(fragment, "reads")
        assert hasattr(fragment.reads, "__iter__")
    assert hasattr(fragment_pe, "Read2")
    assert fragment_pe.Read1.Name == "r1"
    assert fragment_pe.Read2.Name == "r2"
    assert not hasattr(fragment_se, "Read2")
    assert fragment_se.Read1.Name == "r3"


def test_temporarytopermanent_init(tmp_path):
    permanent_file = tmp_path / "permanent.txt"
    temp_to_perm = TemporaryToPermanent(permanent_file)
    assert temp_to_perm.permanent_file == permanent_file
    # assert nottemp_to_perm.tmp_path.exists()
    # assert str(temp_to_perm.tmp_path) == temp_to_perm.tmp_directory.name
    # assert temp_to_perm.temp_file == temp_to_perm.tmp_path / tmp_path.name / "permanent.txt"
    # assert not temp_to_perm.temp_file.parent.exists()


def test_temporarytopermanent_open(tmp_path):
    permanent_file = tmp_path / "permanent.txt"
    temp_to_perm = TemporaryToPermanent(permanent_file)
    assert not hasattr(temp_to_perm, "file_handle")
    fh = temp_to_perm.open("w")
    assert hasattr(temp_to_perm, "file_handle")
    assert hasattr(temp_to_perm, "tmp_path")
    assert temp_to_perm.tmp_path.exists()
    assert hasattr(temp_to_perm, "tmp_file")
    assert temp_to_perm.temp_file.exists()
    assert not permanent_file.exists()


def test_temporarytopermanent_close():
    permanent_file = tmp_path / "permanent.txt"
    temp_to_perm = TemporaryToPermanent(permanent_file)
    assert not permanent_file.exists()
    temporary_path = ""
    temporary_file = ""
    with temp_to_perm.open("w"):
        temporary_path = str(temp_to_perm.tmp_path)
        temporary_file = str(temp_to_perm.temp_file)
        assert temp_to_perm.tmp_path.exists()
        assert temp_to_perm.temp_file.exists()
    assert not Path(temporary_path).exists()
    assert not Path(temporary_file).exists()
    assert not hasattr(temp_to_perm, "file_handle")
    assert not hasattr(temp_to_perm, "tmp_path")
    assert not hasattr(temp_to_perm, "tmp_file")
    assert permanent_file.exists()


def test_temporarytopermanent_with():
    pass


def test_temporarytopermanent_write():
    permanent_file = tmp_path / "permanent.txt"
    with TemporaryToPermanent(permanent_file).open("w") as fh:
        fh.write("something")
    assert permanent_file.exists()
    with permanent_file.open("") as inp:
        readline = inp.readline()
    assert readline == "something"


def test_temporarytopermanent_interrupt():
    permanent_file = tmp_path / "permanent.txt"
    temp_to_perm = TemporaryToPermanent(permanent_file)
    assert not permanent_file.exists()
    with temp_to_perm.open("w") as fh:
        fh.write("something")
    fh.close()
    assert permanent_file.exists()
    with permanent_file.open("") as inp:
        readline = inp.readline()
    assert readline == "something"


def test_temporarytopermanent_interrupt_no_context_manager():
    pass


def test_reverse_complement(MockBlockedFileAdapter):
    assert reverse_complement("ATCG") == "CGAT"


def test_iterate_fastq():
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        iterator = iterate_fastq("_R1_", reverse_reads=False)
        for name, seq, qual in iterator:
            assert name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
            assert (
                seq
                == "NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
            )
            assert (
                qual
                == "#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFF:FFFFFF"
            )
            break
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter("_R1_")):
        iterator = iterate_fastq("This is not a filehandl", reverse_reads=True)
        for name, seq, qual in iterator:
            assert name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
            assert seq == reverse_complement(
                "NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
            )
            break


def test_get_fastq_iterator(se_sample, pe_sample):
    iterator = get_fastq_iterator(pe_sample.is_paired)
    assert iterator.__name__ == "_iterreads_paired_end"
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        fragments = [fragment for fragment in iterator(pe_sample.filenames)]
        assert fragments[0].Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
        assert fragments[0].Read2.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 2:N:0:TAGCTT"

    iterator = get_fastq_iterator(se_sample.is_paired)
    assert iterator.__name__ == "_iterreads_single_end"
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        fragments = [fragment for fragment in iterator(se_sample.filenames)]
        assert fragments[0].Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
        assert not hasattr(fragments[0], "Read2")
