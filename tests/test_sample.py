# -*- coding: utf-8 -*-

import pytest
import mbf.align
import os
from mmdemultiplex.samples import DemultiplexInputSample, FASTQsFromJobSelect
from typing import Tuple, List
from pypipegraph import MultiFileGeneratingJob, Job

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


@pytest.mark.usefixtures("new_pipegraph")
def test_demultiplex_input_sample(tmp_path):
    r1 = tmp_path / "bla_R1_.fastq"
    r2 = tmp_path / "bla_R2_.fastq"
    r1.touch()
    r2.touch()
    sample = DemultiplexInputSample(
        "mysample",
        mbf.align.strategies.FASTQsFromFolder(tmp_path),
        reverse_reads=False,
    )
    assert sample.name == "mysample"
    assert isinstance(sample.input_strategy, mbf.align.strategies.FASTQsFromFolder)
    assert not sample.reverse_reads
    assert isinstance(sample.fastq_processor, mbf.align.fastq2.Straight)
    assert sample.pairing == "paired"
    for job in sample.dependencies:
        assert isinstance(job, Job)
    assert sample.is_paired
    sample = DemultiplexInputSample(
        "mysample",
        mbf.align.strategies.FASTQsFromFolder(tmp_path),
        reverse_reads=True,
        pairing="paired_as_single",
    )
    assert sample.pairing == "paired_as_single"
    for job in sample.dependencies:
        assert isinstance(job, Job)
    assert not sample.is_paired
    os.remove(r2)
    strat = mbf.align.strategies.FASTQsFromFolder(tmp_path)
    strat.dependencies = [MultiFileGeneratingJob(["somefiles"], lambda: None)]
    sample = DemultiplexInputSample(
        "mysample",
        strat,
        reverse_reads=True,
        pairing="single",
    )
    assert sample.reverse_reads
    assert sample.pairing == "single"
    for job in sample.dependencies:
        assert isinstance(job, Job)
    assert not sample.is_paired


@pytest.mark.usefixtures("new_pipegraph")
def test_demultiplex_input_files_paired(tmp_path):
    r1 = tmp_path / "bla_R1_.fastq"
    r2 = tmp_path / "bla_R2_.fastq"
    r1.touch()
    r2.touch()
    sample = DemultiplexInputSample(
        "mysample",
        mbf.align.strategies.FASTQsFromFolder(tmp_path),
        reverse_reads=False,
    )
    input_files = sample.get_aligner_input_filenames()
    assert isinstance(input_files, List)
    assert isinstance(input_files[0], Tuple)
    assert input_files[0] == (str(r1), str(r2))
    assert isinstance(input_files[0][0], str)
    assert isinstance(input_files[0][1], str)
    assert r1.exists()
    assert r2.exists()


@pytest.mark.usefixtures("new_pipegraph")
def test_demultiplex_input_files_single(tmp_path):
    r1 = tmp_path / "bla_R1_.fastq"
    r1.touch()
    sample = DemultiplexInputSample(
        "mysample",
        mbf.align.strategies.FASTQsFromFolder(tmp_path),
        reverse_reads=False,
        pairing="single",
    )
    input_files = sample.get_aligner_input_filenames()
    assert isinstance(input_files, List)
    assert isinstance(input_files[0], str)
    assert input_files[0] == str(r1)
    assert r1.exists()


@pytest.mark.usefixtures("new_pipegraph")
def test_demultiplex_input_files_raise_exceptions(tmp_path):
    r1 = tmp_path / "bla_R1_.fastq"
    r2 = tmp_path / "bla_R2_.fastq"
    r1.touch()
    r2.touch()
    with pytest.raises(ValueError):
        DemultiplexInputSample(
            "mysample",
            mbf.align.strategies.FASTQsFromFolder(tmp_path),
            reverse_reads=False,
            pairing="single",
        )
    os.remove(r2)
    with pytest.raises(ValueError):
        DemultiplexInputSample(
            "mysample",
            mbf.align.strategies.FASTQsFromFolder(tmp_path),
            reverse_reads=False,
            pairing="paired",
        )
    with pytest.raises(ValueError):
        DemultiplexInputSample(
            "mysample",
            mbf.align.strategies.FASTQsFromFolder(tmp_path),
            reverse_reads=False,
            pairing="wrong",
        )


@pytest.mark.usefixtures("new_pipegraph")
def test_FASTQsFromJobSelect_init():
    filenames = [
        "sample1_R1_.fastq",
        "sample1_R2_.fastq",
        "sample2_R1_.fastq",
        "sample2_R2_.fastq",
    ]
    job = MultiFileGeneratingJob(filenames, lambda: None)
    sample_name = "sample1"
    select = FASTQsFromJobSelect(sample_name, job)
    assert select.dependencies == job
    assert select.job == job
    assert select.sample_name == sample_name


@pytest.mark.usefixtures("new_pipegraph")
def test_FASTQsFromJobSelect_parse(tmp_path):
    filenames = [
        "sample1_R1_.fastq",
        "sample1_R2_.fastq",
        "sample2_R1_.fastq",
        "sample2_R2_.fastq",
    ]
    job = MultiFileGeneratingJob(filenames, lambda: None)
    sample_name = "sample1"
    select = FASTQsFromJobSelect(sample_name, job)
    filenames = [filepath.name for tuple in select() for filepath in tuple]
    assert "sample1_R1_.fastq" in filenames
    assert "sample1_R2_.fastq" in filenames
    assert "sample2_R1_.fastq" not in filenames
    assert "sample2_R2_.fastq" not in filenames
