# -*- coding: utf-8 -*-

import pytest
import pandas as pd
from mmdemultiplex import Demultiplexer, DemultiplexStrategy, PE_Decide_On_Start_Trim_Start_End
from pathlib import Path
from unittest.mock import patch
from pypipegraph import Job

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def ttttest_fragment():
    assert hasattr(fragment, "Read1")
    assert hasattr(fragment, "Read2")
    assert hasattr(fragment.Read1, "Name")
    assert hasattr(fragment.Read2, "Name")
    assert hasattr(fragment.Read1, "Sequence")
    assert hasattr(fragment.Read2, "Sequence")
    assert hasattr(fragment.Read1, "Quality")
    assert hasattr(fragment.Read2, "Quality")
