#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""adapters.py: Contains adapter classes and an interface to cutadapt."""

from typing import Union, Optional, Tuple
import cutadapt
import cutadapt.align
import collections


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


# To allow skipping of a prefix of the reference at no cost, set the
# |  START_IN_REFERENCE flag.
# |  To allow skipping of a prefix of the query at no cost, set the
# |  START_IN_QUERY flag.
# |  If both are set, a prefix of the reference or the query is skipped,
# |  never both.
# |  Similarly, set STOP_IN_REFERENCE and STOP_IN_QUERY to
# |  allow skipping of suffixes of the reference or of the query. Again, it
# |  is never the case that both suffixes are skipped.
# |  If all flags are set, this results in standard semiglobal alignment.


WHERE_START = 11  # WHERE_START should be start_in_reference & start_in_query & stop_in_query = 0b1011 = 11, that means we can get adapters within the reference, but the start must be present
WHERE_END = 14  # WHERE_END should be stop_in_reference & start_in_query & stop_in_query = 0b1110 = 14, that means we can get adapters within the reference, but the end must be present
WHERE_SEMIGLOBAL = 15
# WHERE_START = (
#     cutadapt.align.START_WITHIN_SEQ1  # you need this, if you don't want partial matches with the adapter start at the end
#     | cutadapt.align.START_WITHIN_SEQ2
#     | cutadapt.align.STOP_WITHIN_SEQ2
# )
# WHERE_END = (
#     cutadapt.align.STOP_WITHIN_SEQ1  # you need this, if you want partial matches with the adapter end at the start
#     | cutadapt.align.START_WITHIN_SEQ2
#     | cutadapt.align.STOP_WITHIN_SEQ2
# )
# WHERE_BOTH = (
#     cutadapt.align.START_WITHIN_SEQ1  # allow partial at front and back
#     | cutadapt.align.STOP_WITHIN_SEQ1
#     | cutadapt.align.START_WITHIN_SEQ2
#     | cutadapt.align.STOP_WITHIN_SEQ2
# )
# WHERE_NONE = cutadapt.align.START_WITHIN_SEQ2 | cutadapt.align.STOP_WITHIN_SEQ2


class Adapter:
    def __init__(
        self,
        adapter_sequence: str,
        maximal_number_of_errors: int = 0,
        index_adapter_end: bool = True,
        minimal_overlap: Optional[int] = None,
        find_right_most_occurence: bool = False,
        find_best_match: bool = True,
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
        self.flags = WHERE_START  # this might be wrong
        self.find_best_match = find_best_match
        self.cutadapt_caller = self.cutadapt_call_first
        if self.find_best_match:
            self.cutadapt_caller = self.cutadapt_call_best
        if self.find_right_most:  # turn everything around
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
                self.error_rate = 0.0
            else:
                self.locate = self.cutadapt_locate
                self.error_rate = self.maximal_number_of_errors / float(
                    self.minimal_overlap
                )
                self.adapter = cutadapt.align.Aligner(
                    reference=self.adapter_sequence,
                    max_error_rate=self.error_rate,
                    flags=self.flags,
                    min_overlap=self.minimal_overlap,
                )
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

    def cutadapt_match(self, sequence: str) -> Optional[int]:
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
        alignment = self.cutadapt_caller(sequence)
        if alignment is None:
            return None
        else:
            return self.get_position_from_adaptermatch(AdapterMatch(*alignment))

    def cutadapt_call_best(self, sequence: str) -> Tuple[int, int, int, int, int, int]:
        """
        Calls cutadapt multiple times to return the best match of an adapter.

        Parameters
        ----------
        sequence : str
            the sequence to search

        Returns
        -------
        Tuple[int, int, int, int, int, int]
            the alignment
        """
        alignment = self.cutadapt_call_first(sequence)
        return self.get_best_from_all(alignment, sequence)

    def cutadapt_call_first(self, sequence: str) -> Tuple[int, int, int, int, int, int]:
        """
        Calls cutadapt once to return the first match of an adapter.

        Parameters
        ----------
        sequence : str
            the sequence to search

        Returns
        -------
        Tuple[int, int, int, int, int, int]
            the alignment
        """
        return self.adapter.locate(sequence)

    def get_best_from_all(
        self, alignment: Tuple[int, int, int, int, int, int], sequence: str
    ) -> Tuple[int, int, int, int, int, int]:
        """
        get_best_from_all calls cudaopt.align.locate multiple times to get all
        possible alignments, then returns the one with the best score.

        Parameters
        ----------
        alignment : Tuple[int, int, int, int, int, int]
            the first found alignment
        sequence : str
            the sequence to search

        Returns
        -------
        Tuple[int, int, int, int, int, int]
            the best scored alignment
        """
        # if we need the rightmost, a new alignment must have the same or better score
        # if we choose the leftmoost, a new alignment must have a better score
        best_alignment = alignment
        current = self.adapter.locate(sequence)
        index = 0
        while current is not None:
            if (
                current[-1] < best_alignment[-1]
            ):  # this one is better, or if finding_right_most, this one is at least equivalent
                best_alignment = (
                    current[:2] + (current[2] + index, current[3] + index) + current[4:]
                )
            sequence = sequence[current[3] :]
            index += current[3]
            current = self.adapter.locate(sequence)
        return best_alignment

    def exact_locate(self, sequence: str) -> Union[int, None]:
        sequence = self.orient_read(sequence)
        pos = sequence.find(self.adapter_sequence)  # find the exact sequence is faster
        if pos >= 0:
            ret = self.correct_for_adapter_location(pos) * self.factor
        else:
            ret = None
        return ret

    def cutadapt_locate(self, sequence: str) -> Union[int, None]:
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
        ret = self.exact_locate(sequence)
        if ret is None:
            optpos = self.cutadapt_match(sequence)
            if optpos is not None:
                ret = optpos * self.factor
            else:
                ret = optpos
        return ret

    def locate_null(self, sequence: str) -> Union[int, None]:
        return 0

    def __str__(self):
        return f"Adapter: {self.adapter_sequence}, errors: {self.maximal_number_of_errors}, overlap: {self.minimal_overlap}, index_adapter_end: {self.index_adapter_end}, find_right_most_occurence: {self.find_right_most}"
