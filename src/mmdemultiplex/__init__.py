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
from .strategies import PE_Decide_On_Start_Trim_Start_End, DemultiplexStrategy
from .util import Fragment, Read
from .samples import DemultiplexInputSample
from .plots import *