# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from .demultiplex import Demultiplexer
from .trim import Trimmer
from .strategies import (
    PE_Decide_On_Start_Trim_Start_End,
    DemultiplexStrategy,
    PE_Decide_On_Start_End_Trim_Start_End,
    SE_Decide_On_Start_Trim_Start_End,
    SE_Trim_On_Start_Trim_After_X_BP,
    PE_Trim_On_Start_Trim_After_X_BP,
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
