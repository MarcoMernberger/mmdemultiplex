# -*- coding: utf-8 -*-

import pytest
import mmdemultiplex
from mmdemultiplex import (
    DemultiplexStrategy,
    PE_Decide_On_Start_Trim_Start_End,
)
from mmdemultiplex.util import Read, Fragment, reverse_complement
from conftest import barcode_df_full_callback, barcode_df_callback


__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_trim_read_front():
    read = Read("name", "ATCTC", "FFGGG")
    index = 2
    new_read = DemultiplexStrategy.trim_read_front(read, index)
    assert new_read.Name.split("_")[-1] == "AT"
    assert new_read.Sequence == "CTC"
    assert new_read.Quality == "GGG"
    read = Read("name", "ATCTC", "FFGGG")
    index = 200
    new_read = DemultiplexStrategy.trim_read_front(read, index)
    assert new_read.Name.split("_")[-1] == "ATCTC"
    assert new_read.Sequence == ""
    assert new_read.Quality == ""


def test_trim_read_back():
    read = Read("name", "ATCTC", "FFFGG")
    index = 3
    new_read = DemultiplexStrategy.trim_read_back(read, index)
    assert new_read.Name.split("_")[-1] == "TC"
    assert new_read.Sequence == "ATC"
    assert new_read.Quality == "FFF"
    read = Read("name", "ATCTC", "FFFGG")
    index = 200
    new_read = DemultiplexStrategy.trim_read_back(read, index)
    assert new_read.Name.split("_")[-1] == ""
    assert new_read.Sequence == "ATCTC"
    assert new_read.Quality == "FFFGG"


def test_PE_Decide_On_Start_Trim_Start_End_init():
    minimal = PE_Decide_On_Start_Trim_Start_End("ATG", "TGA")
    assert minimal.start_barcode == "ATG"
    assert minimal.end_barcode == "TGA"
    assert minimal.trim_before_end == 0
    assert minimal.trim_after_start == 0
    assert minimal.minimal_overlap_start == len(minimal.start_barcode)
    assert minimal.minimal_overlap_end == len(minimal.end_barcode)
    assert minimal.maximal_errors_start == 0
    assert minimal.maximal_errors_end == 0
    params = {
        "start_barcode": "ATTT",
        "end_barcode": "TAAA",
        "trim_after_start": 10,
        "trim_before_end": 12,
        "maximal_errors_start": 2,
        "maximal_errors_end": 0,
        "minimal_overlap_start": 5,
        "minimal_overlap_end": 6,
    }
    full = PE_Decide_On_Start_Trim_Start_End(**params)
    assert full.start_barcode == "ATTT"
    assert full.end_barcode == "TAAA"
    assert full.trim_before_end == 12
    assert full.trim_after_start == 10
    assert full.maximal_errors_start == 2
    assert full.maximal_errors_end == 0
    assert full.minimal_overlap_start == 5
    assert full.minimal_overlap_end == 6
    params = {"key": 1, "not_a_param": 1}
    with pytest.raises(TypeError):
        PE_Decide_On_Start_Trim_Start_End(**params)


def test_PE_Decide_On_Start_Trim_Start_End_init_from_df():
    df = barcode_df_callback()
    for key, row in df.iterrows():
        minimal = PE_Decide_On_Start_Trim_Start_End(**row.to_dict())
        assert minimal.start_barcode == row["start_barcode"]
    df = barcode_df_full_callback()
    for key, row in df.iterrows():
        full = PE_Decide_On_Start_Trim_Start_End(**row.to_dict())
        assert full.start_barcode == row["start_barcode"]

    for key, row in df.iterrows():
        full = PE_Decide_On_Start_Trim_Start_End(**row.to_dict())
        assert full.start_barcode == row["start_barcode"]


def test_get_parameters_PE_Decide_On_Start_Trim_Start_End():
    params = {
        "start_barcode": "ATTT",
        "end_barcode": "TAAA",
        "trim_after_start": 10,
        "trim_before_end": 12,
        "maximal_errors_start": 2,
        "maximal_errors_end": 0,
        "minimal_overlap_start": 5,
        "minimal_overlap_end": 6,
    }
    strat = PE_Decide_On_Start_Trim_Start_End(**params)
    parameters = strat.get_parameters()
    for value in params.values():
        assert value in parameters


def test_init_adapter_PE_Decide_On_Start_Trim_Start_End():
    params = {
        "start_barcode": "ATTT",
        "end_barcode": "TAAA",
        "trim_after_start": 10,
        "trim_before_end": 12,
        "maximal_errors_start": 2,
        "maximal_errors_end": 0,
        "minimal_overlap_start": 5,
        "minimal_overlap_end": 6,
    }
    strat = PE_Decide_On_Start_Trim_Start_End(**params)
    assert hasattr(strat, "adapter_start_forward")
    assert hasattr(strat, "adapter_start_reverse")
    assert hasattr(strat, "adapter_end_forward")
    assert hasattr(strat, "adapter_end_reverse")
    assert strat.adapter_start_forward.adapter_sequence == strat.start_barcode
    assert strat.adapter_start_forward.maximal_number_of_errors == strat.maximal_errors_start
    assert strat.adapter_start_forward.where == mmdemultiplex.adapters.WHERE_START
    assert strat.adapter_start_forward.minimal_overlap == strat.minimal_overlap_start
    assert (
        strat.adapter_start_reverse.adapter_sequence
        == reverse_complement(strat.start_barcode)[::-1]
    )
    assert strat.adapter_start_reverse.maximal_number_of_errors == strat.maximal_errors_start
    assert strat.adapter_start_reverse.where == mmdemultiplex.adapters.WHERE_START
    assert strat.adapter_start_reverse.minimal_overlap == strat.minimal_overlap_start
    assert strat.adapter_end_forward.adapter_sequence == strat.end_barcode
    assert strat.adapter_end_forward.maximal_number_of_errors == strat.maximal_errors_end
    assert strat.adapter_end_forward.where == mmdemultiplex.adapters.WHERE_START
    assert strat.adapter_end_forward.minimal_overlap == strat.minimal_overlap_end
    assert strat.adapter_end_reverse.adapter_sequence == reverse_complement(strat.end_barcode)[::-1]
    assert strat.adapter_end_reverse.maximal_number_of_errors == strat.maximal_errors_end
    assert strat.adapter_end_reverse.where == mmdemultiplex.adapters.WHERE_START
    assert strat.adapter_end_reverse.minimal_overlap == strat.minimal_overlap_end


def test_match_and_trim_PE_Decide_On_Start_Trim_Start_End(paired_fragments):
    df = barcode_df_full_callback()
    amplicon = "ATCGATCG"
    amplicon_rev = "CGATCGAT"
    row = df.loc["Sample2"]
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    accepted = matcher.match_and_trim(paired_fragments["Sample2"])
    assert accepted.Read1.Sequence == amplicon
    assert accepted.Read2.Sequence == amplicon_rev
    accepted = matcher.match_and_trim(paired_fragments["Sample1"])
    assert not accepted
    accepted = matcher.match_and_trim(paired_fragments["Sample3"])
    assert not accepted
    accepted = matcher.match_and_trim(paired_fragments["Sample2-1mismatch"])
    assert accepted.Read1.Sequence == amplicon
    assert accepted.Read2.Sequence == amplicon_rev
    accepted = matcher.match_and_trim(paired_fragments["Sample2-empty"])
    assert not accepted
    accepted = matcher.match_and_trim(paired_fragments["Sample-discard"])
    assert not accepted
    accepted = matcher.match_and_trim(paired_fragments["Sample2-no-adapter"])
    assert not accepted
    accepted = matcher.match_and_trim(paired_fragments["Sample2-reversed"])
    assert accepted.Read1.Sequence == amplicon
    assert accepted.Read2.Sequence == amplicon_rev


def test_match_and_trim_PE_Decide_On_Start_Trim_Start_End_fw_rev_same(paired_fragments):
    df = barcode_df_full_callback()
    amplicon = "ATCGATCG"
    amplicon_rev = "CGATCGAT"
    row = df.loc["Sample4"]
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    accepted = matcher.match_and_trim(paired_fragments["Sample4-fw-rev-same"])
    assert accepted
    assert accepted.Read1.Sequence == amplicon
    assert accepted.Read2.Sequence == amplicon_rev


def test_match_and_trim_PE_Decide_On_Start_truncated_and_error(paired_fragments):
    df = barcode_df_full_callback()
    amplicon = "ATCGATCG"
    amplicon_rev = "CGATCGAT"
    row = df.loc["Sample3"]
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    accepted = matcher.match_and_trim(paired_fragments["Sample3-partial-and-error"])
    assert accepted
    assert accepted.Read1.Sequence == amplicon
    assert accepted.Read2.Sequence == amplicon_rev


def test_match_and_trim_PE_Decide_On_Start_accept_all_no_barcode(paired_fragments):
    df = barcode_df_full_callback()
    row = df.loc["Sample4"]
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    accepted = matcher.match_and_trim(paired_fragments["Sample3-partial-and-error"])
    assert not accepted
    row["start_barcode"] = ""
    row["end_barcode"] = ""
    r1 = paired_fragments["Sample3-partial-and-error"].Read1
    r2 = paired_fragments["Sample3-partial-and-error"].Read2
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    accepted = matcher.match_and_trim(paired_fragments["Sample3-partial-and-error"])
    assert accepted
    assert accepted.Read1.Sequence == r1.Sequence
    assert accepted.Read2.Sequence == r2.Sequence


def test_example_PE_Decide_On_Start_Trim_Start_End():
    row = {
        "start_barcode": "CGTAAACTCACTG",
        "end_barcode": "TCATGTAGCTCTG",
        "trim_after_start": 13,
        "trim_before_end": 13
    }
    matcher = PE_Decide_On_Start_Trim_Start_End(**row)
    r1 = Read(
        'M03491:11:000000000-KFFWP:1:1101:16391:1909 1:N:0:1',
        'CGTAAACTCACTGGTGCTGTGACTGCTTGTAGATGGCCA',
        'AAAAA1FFFFFFGGGGGGEG1FGHG0AFEHHFHF31BFH',
    )
    r2 = Read(
        'M03491:11:000000000-KFFWP:1:1101:16391:1909 2:N:0:1',
        'TCATGTAGCTCTGCCTGTCTTTCAACTCTGTCTCCTTCC',
        '>AAAAFFFFFFFG111G1133ADD33DGH3B3AF11A10',
        )
    fragment = Fragment(r1, r2)
    accepted = matcher.match_and_trim(fragment)
    assert accepted
