# -*- coding: utf-8 -*-
from importlib.metadata import PackageNotFoundError, version

# Change here if project/distribution name differs from package name
dist_name = __name__
try:
    __version__ = version("mmdemultiplex")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.5.0"


from .demultiplex import Demultiplexer
from .trim import Trimmer
from .strategies import (
    DemultiplexStrategy,
    PE_Decide_On_Start_End_Trim_Start_End,
    SE_Decide_On_Start_Trim_Start_End,
    SE_Trim_On_Start_Trim_After_X_BP,
    PE_Trim_On_Start_Trim_After_X_BP,
    PE_Decide_On_Start_End_Trim_Start_End_Force_Barcode_at_Front,
)
from .util import (
    Fragment,
    Read,
    get_df_callable_for_demultiplexer,
    reverse_complement,
    get_fastq_iterator,
    dump_matching_reads,
)
from .filterfastq import FastqDemultiplexer, decision_callback_init
from .samples import DemultiplexInputSample
from .plots import *
from .filter import SequenceFilter
from .jobs import dump_matching_reads_job, plot_adapter_heatmap
from .adapters import Adapter

__all__ = [
    "Demultiplexer",
    "Trimmer",
    "DemultiplexStrategy",
    "PE_Decide_On_Start_End_Trim_Start_End",
    "SE_Decide_On_Start_Trim_Start_End",
    "SE_Trim_On_Start_Trim_After_X_BP",
    "PE_Trim_On_Start_Trim_After_X_BP",
    "PE_Decide_On_Start_End_Trim_Start_End_Force_Barcode_at_Front",
    "Fragment",
    "Read",
    "get_df_callable_for_demultiplexer",
    "reverse_complement",
    "get_fastq_iterator",
    "dump_matching_reads",
    "FastqDemultiplexer",
    "decision_callback_init",
    "DemultiplexInputSample",
    "plot_adapter_heatmap",
    "dump_matching_reads_job",
    "Adapter",
    "SequenceFilter",
]
