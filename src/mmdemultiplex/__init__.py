# -*- coding: utf-8 -*-
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mmdemultiplex")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.5.0"


from .demultiplex import Demultiplexer
from .trim import Trimmer
from .strategies import (
    PE_Decide_On_Start_Trim_Start_End,
    DemultiplexStrategy,
    PE_Decide_On_Start_End_Trim_Start_End,
    SE_Decide_On_Start_Trim_Start_End,
    SE_Trim_On_Start_Trim_After_X_BP,
)
from .util import (
    Fragment,
    Read,
    get_df_callable_for_demultiplexer,
    reverse_complement,
    get_fastq_iterator,
)

from .samples import DemultiplexInputSample
from .plots import *
from .filter import SequenceFilter
