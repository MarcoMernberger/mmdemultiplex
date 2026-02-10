# -*- coding: utf-8 -*-
try:
    # Py3.8+
    from importlib.metadata import (
        version as _get_version,
        PackageNotFoundError as _PackageNotFoundError,
    )
except Exception:
    # Fallback for older python
    from importlib_metadata import version as _get_version, PackageNotFoundError as _PackageNotFoundError  # type: ignore

# Change here if project/distribution name differs from package name
dist_name = __name__
try:
    __version__ = _get_version(dist_name)
except _PackageNotFoundError:
    __version__ = "unknown"
finally:
    del _get_version, _PackageNotFoundError


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
