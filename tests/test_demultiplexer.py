# -*- coding: utf-8 -*-

import pytest
import pandas as pd
import pypipegraph as ppg
from mmdemultiplex import (
    Demultiplexer,
    DemultiplexStrategy,
    PE_Decide_On_Start_Trim_Start_End,
    Fragment,
)
from mmdemultiplex.util import get_fastq_iterator
from pathlib import Path
from unittest.mock import patch
from pypipegraph import Job
from conftest import MockBlockedFileAdapter

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def barcode_df_callback():
    df = pd.DataFrame({"key": ["Sample1"], "start_barcode": ["ATCG"], "end_barcode": ["TTCG"]})
    df = df.set_index("key")
    return df


class MockDecisionCallback:
    def __init__(self):
        self.name = "MockStrategy"

    def match_and_trim(self, fragment):
        if fragment.Read1.Name == "A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT":
            return fragment
        else:
            return False


def test_init(tmp_path, pe_sample, se_sample):
    barcode_df = barcode_df_callback()
    demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
    assert demultiplexer.name == f"DM_{pe_sample.name}"
    assert hasattr(demultiplexer, "barcode_df")
    assert isinstance(demultiplexer.barcode_df, pd.DataFrame)
    assert demultiplexer.library_name == pe_sample.name
    assert demultiplexer.output_folder == tmp_path / demultiplexer.name
    assert demultiplexer.output_folder.exists()
    assert demultiplexer.input_sample == pe_sample
    assert demultiplexer.maximal_error_rate == 0
    assert demultiplexer.input_files == pe_sample.get_aligner_input_filenames()
    assert issubclass(demultiplexer.strategy, DemultiplexStrategy)
    assert demultiplexer.is_paired
    assert hasattr(demultiplexer, "decision_callbacks")
    assert len(demultiplexer.decision_callbacks) == barcode_df.shape[0]
    assert isinstance(
        demultiplexer.decision_callbacks["Sample1"], PE_Decide_On_Start_Trim_Start_End
    )
    # with SE
    demultiplexer = Demultiplexer(se_sample, barcode_df_callback, output_folder=tmp_path)
    assert demultiplexer.name == f"DM_{se_sample.name}"
    assert not demultiplexer.is_paired
    # with default output folder
    with patch("pathlib.Path.mkdir", return_value=None):
        demultiplexer = Demultiplexer(pe_sample, barcode_df_callback)
        assert demultiplexer.output_folder == (Path("cache") / demultiplexer.name)


@pytest.mark.usefixtures("new_pipegraph")
def test_get_dependencies(tmp_path, pe_sample):
    demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
    for job in demultiplexer.get_dependencies():
        assert isinstance(job, Job) or isinstance(job, str)


def test_parameters(tmp_path, pe_sample):
    demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
    parameters = demultiplexer.parameters
    assert parameters[0] == demultiplexer.name
    assert parameters[1] == pe_sample.name
    assert parameters[2] == "PE_Decide_On_Start_Trim_Start_End"
    for item in parameters:
        assert hasattr(item, "__hash__")


def test_parameter_string(tmp_path, pe_sample):
    demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
    assert isinstance(demultiplexer.parameter_string(), str)


def test_fastq_iterator(tmp_path, pe_sample, se_sample):
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
        iterator = demultiplexer.get_fastq_iterator()
        assert callable(iterator)
        for files_tuple in demultiplexer.input_files:
            for fragment in iterator(files_tuple):
                assert isinstance(fragment, Fragment)
                assert hasattr(fragment, "Read1")
                assert hasattr(fragment, "Read2")
                assert hasattr(fragment, "reads")
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        demultiplexer = Demultiplexer(se_sample, barcode_df_callback, output_folder=tmp_path)
        iterator = demultiplexer.get_fastq_iterator()
        assert callable(iterator)
        for files_tuple in demultiplexer.input_files:
            for fragment in iterator(files_tuple):
                assert isinstance(fragment, Fragment)
                assert hasattr(fragment, "Read1")
                assert not hasattr(fragment, "Read2")
                assert hasattr(fragment, "reads")


@pytest.mark.usefixtures("new_pipegraph")
def test_do_demultiplex_pe(tmp_path, pe_sample):
    first_read_sample_name = f"{pe_sample.name}_first_read"
    discarded_sample_name = f"{pe_sample.name}_discarded"
    files_created = {
        "first_read": (
            tmp_path / first_read_sample_name / f"{first_read_sample_name}_R1_.fastq",
            tmp_path / first_read_sample_name / f"{first_read_sample_name}_R2_.fastq",
        ),
        "discarded": (
            tmp_path / discarded_sample_name / f"{discarded_sample_name}_R1_.fastq",
            tmp_path / discarded_sample_name / f"{discarded_sample_name}_R2_.fastq",
        ),
    }
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
        demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
        job = demultiplexer.do_demultiplex()
        print(job.dependencies)
        ppg.run_pipegraph()
        sentinel = demultiplexer.output_folder / "done.txt"
        filepaths = [sentinel] + [
            filename for filetuple in files_created.values() for filename in filetuple
        ]
        for filepath in filepaths:
            assert filepath.exists()
    fastq_iterator = get_fastq_iterator(pe_sample.is_paired)
    for fragment in fastq_iterator(files_created["first_read"]):
        assert fragment.Read1.Name == b"A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
        assert fragment.Read2.Name == b"A01284:56:HNNKWDRXY:1:2101:1524:1000 2:N:0:TAGCTT"
        assert (
            fragment.Read1.Sequence
            == b"NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
        )
        assert (
            fragment.Read2.Sequence
            == b"NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC"
        )
        assert (
            fragment.Read1.Quality
            == b"#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFF:FFFFFF"
        )
        assert (
            fragment.Read2.Quality
            == b"#FF,FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
        )
        break
    for fragment in fastq_iterator(files_created["discarded"]):
        assert fragment.Read1.Name == b"A01284:56:HNNKWDRXY:1:2101:2248:1000 1:N:0:TAGCTT"
        assert fragment.Read2.Name == b"A01284:56:HNNKWDRXY:1:2101:2248:1000 2:N:0:TAGCTT"
        assert (
            fragment.Read1.Sequence
            == b"NTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCA"
        )
        assert (
            fragment.Read2.Sequence
            == b"NAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGC"
        )
        assert (
            fragment.Read1.Quality
            == b"#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFF,FFFFFFFFFF:FFFF"
        )
        assert (
            fragment.Read2.Quality
            == b"#FF:FFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF"
        )
        break


@pytest.mark.usefixtures("new_pipegraph")
def test_do_demultiplex_se(tmp_path, se_sample):
    first_read_sample_name = f"{se_sample.name}_first_read"
    discarded_sample_name = f"{se_sample.name}_discarded"
    files_created = {
        "first_read": (tmp_path / first_read_sample_name / f"{first_read_sample_name}_R1_.fastq",),
        "discarded": (tmp_path / discarded_sample_name / f"{discarded_sample_name}_R1_.fastq",),
    }
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        demultiplexer = Demultiplexer(se_sample, barcode_df_callback, output_folder=tmp_path)
        demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
        demultiplexer.do_demultiplex()
        ppg.run_pipegraph()
        sentinel = demultiplexer.output_folder / "done.txt"
        filepaths = [sentinel] + [
            filename for filetuple in files_created.values() for filename in filetuple
        ]
        for filepath in filepaths:
            assert filepath.exists()
    fastq_iterator = get_fastq_iterator(se_sample.is_paired)
    for fragment in fastq_iterator(files_created["first_read"]):
        assert fragment.Read1.Name == b"A01284:56:HNNKWDRXY:1:2101:1524:1000 1:N:0:TAGCTT"
        assert not hasattr(fragment, "Read2")
        break
    for fragment in fastq_iterator(files_created["discarded"]):
        assert fragment.Read1.Name == b"A01284:56:HNNKWDRXY:1:2101:2248:1000 1:N:0:TAGCTT"
        break


def test_decide_on_barcode(tmp_path, pe_sample):
    with patch("mbf_align._common.BlockedFileAdaptor", MockBlockedFileAdapter):
        fragments = []
        demultiplexer = Demultiplexer(pe_sample, barcode_df_callback, output_folder=tmp_path)
        demultiplexer.decision_callbacks = {"first_read": MockDecisionCallback()}
        iterator = demultiplexer.get_fastq_iterator()
        for files_tuple in demultiplexer.input_files:
            for fragment in iterator(files_tuple):
                fragments.append(fragment)
        first_result = demultiplexer._decide_on_barcode(fragments[0])
        assert first_result[0] == "first_read"
        assert first_result[1] == fragments[0]
        second_result = demultiplexer._decide_on_barcode(fragments[1])
        assert second_result[0] == "discarded"
        assert second_result[1] == fragments[1]


def t__write_fragment(self, fragment, file_handles):
    for i, read in enumerate(fragment):
        file_handles[i].write(f"@{read.Name}\n{read.Sequence}\n+\n{read.Quality}\n".encode())

