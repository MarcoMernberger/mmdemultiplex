#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""strategies.py: Contains different demultiplexing strategies."""

from abc import abstractmethod, ABC
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union, Literal
from .util import Fragment, Read, reverse_complement
from .adapters import Adapter
import pandas as pd
import pypipegraph as ppg

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class DemultiplexStrategy(ABC):
    @classmethod
    def trim_read_front(self, read: Read, index: int) -> Read:
        read.Name = read.Name + "_" + read.Sequence[:index]
        read.Sequence = read.Sequence[index:]
        read.Quality = read.Quality[index:]
        return read

    @classmethod
    def trim_read_back(self, read: Read, index: int) -> Read:
        read.Name = read.Name + "_" + read.Sequence[index:]
        read.Sequence = read.Sequence[:index]
        read.Quality = read.Quality[:index]
        return read

    @abstractmethod
    def match_and_trim(self, fragment: Fragment) -> Union[Fragment, Literal[False]]:
        pass  # pragma: no cover

    @abstractmethod
    def get_parameters(self) -> List:
        pass  # pragma: no cover


class PE_Decide_On_Start_Trim_Start_End(DemultiplexStrategy):
    """The start barcode alone is relevant"""

    def __init__(
        self,
        start_barcode: str,
        end_barcode: str,
        trim_after_start: int = 0,
        trim_before_end: int = 0,
        maximal_errors_start: int = 0,
        maximal_errors_end: int = 0,
        minimal_overlap_start: int = None,
        minimal_overlap_end: int = None,
    ):
        self.start_barcode = start_barcode
        self.end_barcode = end_barcode
        self.trim_before_end = trim_before_end
        self.trim_after_start = trim_after_start
        self.maximal_errors_start = maximal_errors_start
        self.maximal_errors_end = maximal_errors_end
        self.minimal_overlap_start = minimal_overlap_start
        self.minimal_overlap_end = minimal_overlap_end
        if minimal_overlap_start is None:
            self.minimal_overlap_start = len(self.start_barcode)
        self.minimal_overlap_end = minimal_overlap_end
        if minimal_overlap_end is None:
            self.minimal_overlap_end = len(self.end_barcode)
        self._init_adapter()

    def get_parameters(self) -> List:
        return [
            self.start_barcode,
            self.end_barcode,
            self.trim_before_end,
            self.trim_after_start,
            self.maximal_errors_start,
            self.maximal_errors_end,
            self.minimal_overlap_start,
            self.minimal_overlap_end,
        ]

    def _init_adapter(self) -> None:
        # this is the start of the first read
        self.adapter_start_forward = Adapter(
            self.start_barcode,
            maximal_number_of_errors=self.maximal_errors_start,
            index_adapter_end=True,
            minimal_overlap=self.minimal_overlap_start,
            find_right_most_occurence=False,
        )
        # this is the end of the first read
        self.adapter_end_reverse = Adapter(
            reverse_complement(self.end_barcode),
            maximal_number_of_errors=self.maximal_errors_end,
            index_adapter_end=False,
            minimal_overlap=self.minimal_overlap_end,
            find_right_most_occurence=True,
        )
        # this is the start of the first read
        self.adapter_end_forward = Adapter(
            self.end_barcode,
            maximal_number_of_errors=self.maximal_errors_end,
            index_adapter_end=True,
            minimal_overlap=self.minimal_overlap_end,
            find_right_most_occurence=False,
        )
        # this is the end of the second read
        self.adapter_start_reverse = Adapter(
            reverse_complement(self.start_barcode),
            maximal_number_of_errors=self.maximal_errors_start,
            index_adapter_end=False,
            minimal_overlap=self.minimal_overlap_start,
            find_right_most_occurence=True,
        )

    def match_and_trim(self, fragment: Fragment) -> Union[Fragment, Literal[False]]:
        "Returns None, if the fragment does not match and a trimmed fragment otherwise"
        start_in_r1 = self.adapter_start_forward.locate(fragment.Read1.Sequence)
        start_in_r2 = self.adapter_start_forward.locate(fragment.Read2.Sequence)
        if start_in_r1 is None and start_in_r2 is None:
            # start adapter nowhere to be found, discard
            return False
        elif start_in_r1 is not None and start_in_r2 is None:
            # start adapter in read 1 and not in read 2, nothing to be done
            pass  # pragma: no cover
        elif start_in_r1 is None and start_in_r2 is not None:
            # start adapter in read 2, switch reads
            start_in_r1, start_in_r2 = start_in_r2, start_in_r1
            fragment = Fragment(fragment.Read2, fragment.Read1)  # check in test
        else:  # start_in_r1 and start_in_r2
            if start_in_r2 < start_in_r1:
                # start adapter in read 2 is more likely
                start_in_r1, start_in_r2 = start_in_r2, start_in_r1
                fragment = Fragment(fragment.Read2, fragment.Read1)  # check in test
        # now read 1 ist the one that contains start adapter forward
        fragment.Read1 = self.trim_read_front(
            fragment.Read1, start_in_r1 + self.trim_after_start
        )  # trim start adapter in first read
        end_in_r2 = self.adapter_end_forward.locate(fragment.Read2.Sequence)
        if end_in_r2 is None:
            # end adapter not found, discard
            return False
        fragment.Read2 = self.trim_read_front(fragment.Read2, end_in_r2 + self.trim_before_end)
        # accept fragment, both adapters are found
        # trim reverse adapters at the end of reads
        end_reverse_in_r1 = self.adapter_end_reverse.locate(fragment.Read1.Sequence)
        if end_reverse_in_r1 is not None:
            # that means there is something to trim and it will not leave an empty string (end_reverse_in_r1 == 0)
            fragment.Read1 = self.trim_read_back(
                fragment.Read1, end_reverse_in_r1 - self.trim_before_end
            )
        start_reverse_in_r2 = self.adapter_start_reverse.locate(fragment.Read2.Sequence)
        if start_reverse_in_r2 is not None:
            fragment.Read2 = self.trim_read_back(
                fragment.Read2, start_reverse_in_r2 - self.trim_after_start
            )
        if len(fragment.Read2.Sequence) == 0 or len(fragment.Read1.Sequence) == 0:
            return False
        return fragment
