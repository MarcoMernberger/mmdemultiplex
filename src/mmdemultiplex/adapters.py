#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""adapters.py: Contains adapter classes and an interface to cutadapt."""

from audioop import reverse
from typing import Union, Literal
import cutadapt
import cutadapt.align
import collections

from mmdemultiplex.util import reverse_complement

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


WHERE_START = (
    cutadapt.align.START_WITHIN_SEQ1  # you need this, if you don't want partial matches with the adapter start at the end
    | cutadapt.align.START_WITHIN_SEQ2
    | cutadapt.align.STOP_WITHIN_SEQ2
)
WHERE_END = (
    cutadapt.align.STOP_WITHIN_SEQ1  # you need this, if you want partial matches with the adapter end at the start
    | cutadapt.align.START_WITHIN_SEQ2
    | cutadapt.align.STOP_WITHIN_SEQ2
)
WHERE_BOTH = (
    cutadapt.align.START_WITHIN_SEQ1  # allow partial at front and back
    | cutadapt.align.STOP_WITHIN_SEQ1
    | cutadapt.align.START_WITHIN_SEQ2
    | cutadapt.align.STOP_WITHIN_SEQ2
)
WHERE_NONE = cutadapt.align.START_WITHIN_SEQ2 | cutadapt.align.STOP_WITHIN_SEQ2


class Adapter:
    def __init__(
        self,
        adapter_sequence: str,
        maximal_number_of_errors: int = 0,
        index_adapter_end: bool = True,
        minimal_overlap: int = None,
        find_right_most_occurence: bool = False,
    ):
        """
        This is a wrapper class for cutadapt.Aligner.

        This is wrapper that provides a cutadapt_locate function that returns
        for a specific barcode/adapter sequence the exact position to trim in the
        sequence to be trimmed, depending on the supposed position of the
        adapter sequence (e.g. at the beginning or at the end of the read).
        This is controlled bvia the index_adapter_end variable. This specifies
        if the sequence that should remain is located before or after the
        adapter_sequence. Cutadapt parameter are set accordingly.

        Parameters
        ----------
        adapter_sequence : str
            The adapter sequence to trim.
        maximal_number_of_errors : int, optional
            the maximum number of errors, by default 0.
        index_adapter_end : bool, optional
            Is the sequence of interest before or after the adapter_sequence, by default True.
        minimal_overlap : int, optional
            Minimal overlap for partial adapter occurences, by default adapter length.
        """
        self.adapter_sequence = adapter_sequence
        self.adapter_sequence_length = len(adapter_sequence)
        self.find_right_most = find_right_most_occurence
        self.orient_read = lambda x: x
        self.factor = 1
        self.index_adapter_end = index_adapter_end
        self.where = WHERE_START
        if self.find_right_most:
            self.adapter_sequence = self.adapter_sequence[::-1]
            self.orient_read = lambda x: x[::-1]
            self.factor = -1
            self.index_adapter_end = not self.index_adapter_end
        self.maximal_number_of_errors = maximal_number_of_errors
        self.minimal_overlap = minimal_overlap
        if self.minimal_overlap is None:
            self.minimal_overlap = self.adapter_sequence_length
        self.accept_overlap = self.minimal_overlap < self.adapter_sequence_length
        if self.adapter_sequence_length == 0:
            self.locate = self.locate_null
        else:
            if self.maximal_number_of_errors == 0 and (not self.accept_overlap):
                self.locate = self.exact_locate
                self.error_rate = 0
            else:
                self.locate = self.cutadapt_locate
                self.error_rate = self.maximal_number_of_errors / float(self.minimal_overlap)
            if self.index_adapter_end:
                self.correct_for_adapter_location = (
                    lambda pos: self.adapter_sequence_length + pos
                )  # if the adapter is at the front, we need to trim from the adapter end
                self.get_position_from_adaptermatch = (
                    lambda match: match.rstop
                )  # if the adapter is at the end, we need to trim from the adapter start
            else:
                self.correct_for_adapter_location = (
                    lambda pos: pos
                )  # if its the end adapter, we need to trim at the beginning
                self.get_position_from_adaptermatch = (
                    lambda match: match.rstart
                )  # if its the end adapter, we need to trim at the beginning
                self.adapter = cutadapt.align.Aligner(
                    self.adapter_sequence,
                    self.error_rate,
                    self.where,
                    min_overlap=self.minimal_overlap,
                )
                self.accept_error = self.maximal_number_of_errors > 0

    def cutadapt_match(self, sequence: str) -> Union[int, Literal[False]]:
        """
        cutadapt_match uses the predefined cutadapt.align.Aligner instance
        to locate a partial match of adapter_sequence and returns an
        AdapterMatch object.

        Parameters
        ----------
        sequence : str
            The sequence to search.

        Returns
        -------
        AdapterMatch
            A named tuple containing the start/stop in adapter, start, stop in
            sequuence, the number of mismatches and the error count.
        """
        alignment = self.adapter.locate(sequence)
        if alignment is None:
            return False
        else:
            return self.get_position_from_adaptermatch(AdapterMatch(*alignment))

    def exact_locate(self, sequence: str) -> Union[int, Literal[False]]:
        sequence = self.orient_read(sequence)
        pos = sequence.find(self.adapter_sequence)  # find the exact sequence is faster
        if pos >= 0:
            ret = self.correct_for_adapter_location(pos) * self.factor
        else:
            ret = False
        return ret

    def cutadapt_locate(self, sequence: str) -> Union[int, Literal[False]]:
        """
        locate returns the first occurence of the adapter_sequence in
        sequence.

        Since cutadapt returns the best match, this first tries to find an exact
        match using str.find.
        Otherwise it resorts to cutadat for partial/mismatching occurences,
        Depending on the setting of self.index_adapter_end, this returns the end
        position (if True) or start position (if False) of the adapter occurence
        (0 based). If no adapter_sequence is found, return False.

        Parameters
        ----------
        sequence : str
            The sequence to search.

        Returns
        -------
        Union[int, Literal[False]]
            The position to trim the provided sequence to get rid of the adapter.
            If no adapter is located, return False.
        """
        sequence = self.orient_read(sequence)
        pos = sequence.find(self.adapter_sequence)
        if pos >= 0:
            ret = self.correct_for_adapter_location(pos) * self.factor
        else:
            pos = self.cutadapt_match(sequence)
            if pos:
                ret = pos * self.factor
            else:
                ret = False
        return ret

    def locate_null(self, sequence: str) -> Union[int, Literal[False]]:
        return 0
