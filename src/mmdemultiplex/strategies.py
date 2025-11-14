#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""strategies.py: Contains different demultiplexing strategies."""

from abc import abstractmethod, ABC
from typing import List, Union, Literal, Optional
from .util import Fragment, Read, reverse_complement
from .adapters import Adapter

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class DemultiplexStrategy(ABC):
    @classmethod
    def trim_read_front(self, read: Read, index: int) -> Read:
        name = read.Name + "_" + read.Sequence[:index]
        sequence = read.Sequence[index:]
        quality = read.Quality[index:]
        return Read(name, sequence, quality)

    @classmethod
    def trim_read_back(self, read: Read, index: int) -> Read:
        name = read.Name + "_" + read.Sequence[index:]
        sequence = read.Sequence[:index]
        quality = read.Quality[:index]
        return Read(name, sequence, quality)

    @abstractmethod
    def match_and_trim(self, fragment: Fragment) -> Union[Fragment, Literal[False]]:
        raise NotImplementedError()  # pragma: no cover

    @abstractmethod
    def get_parameters(self) -> List:
        raise NotImplementedError()  # pragma: no cover

    def trim_fragment_to_length(self, fragment: Fragment) -> Fragment:
        if self.trim_length_r1 > 0:
            fragment.Read1 = self.trim_read_to_length(
                fragment.Read1, self.trim_length_r1
            )  # trim first read to specified length
        if self.trim_length_r2 > 0:
            fragment.Read2 = self.trim_read_to_length(
                fragment.Read2, self.trim_length_r2
            )  # trim second read to specified length
        return fragment

    def trim_read_to_length(self, read: Read, trim_length: int) -> Read:
        if (trim_length > 0) and (len(read.Sequence) >= trim_length):
            read = self.trim_read_back(
                read, trim_length
            )  # trim first read to specified length
        return read


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
        minimal_overlap_start: Optional[int] = None,
        minimal_overlap_end: Optional[int] = None,
        **kwargs
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
        # accept fragment, since only the start adapter matters
        fragment.Read1 = self.trim_read_front(
            fragment.Read1, start_in_r1 + self.trim_after_start
        )  # trim start adapter in first read
        end_in_r2 = self.adapter_end_forward.locate(fragment.Read2.Sequence)
        if end_in_r2 is not None:
            # end adapter found, trim it
            fragment.Read2 = self.trim_read_front(
                fragment.Read2, end_in_r2 + self.trim_before_end
            )
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


class SE_Decide_On_Start_Trim_Start_End(DemultiplexStrategy):
    """Only the start barcode is relevant"""

    def __init__(
        self,
        start_barcode: str,
        end_barcode: str,
        trim_after_start: int = 0,
        trim_before_end: int = 0,
        maximal_errors_start: int = 0,
        maximal_errors_end: int = 0,
        minimal_overlap_start: Optional[int] = None,
        minimal_overlap_end: Optional[int] = None,
        **kwargs
    ):
        """
        Demultiplex and Trim strategy based on start and end adapters.

        This strategy assumes that only the first adapter is relevant for demultiplexing.
        The end adapter is only used for trimming. If the start adapter is not found, the
        fragment is discarded. If the start adapter is found, the fragment is trimmed at the
        start adapter. If the end adapter is found, the fragment is trimmed at the end adapter.
        If the end adapter is not found, the fragment is trimmed at the reverse end adapter.
        If the end adapter is not found in any orienttation, the fragment is not trimmed at the end.
        If the fragment is empty after trimming, it is discarded.

        Parameters
        ----------
        start_barcode : str
            Barcode at the start of read
        end_barcode : str
            Barcode at the end of the read, might be there or not
        trim_after_start : int, optional
            an offset position to trim after the start adapter, by default 0
        trim_before_end : int, optional
            an offset position to trim before the end adapter, by default 0
        maximal_errors_start : int, optional
            max number of errors in start adapter allowed for a match, by default 0
        maximal_errors_end : int, optional
            max number of errors in start adapter allowed for a match, by default 0
        minimal_overlap_start : Optional[int], optional
            the amount of overlap for a truncated start adapter, by default None
        minimal_overlap_end : Optional[int], optional
            the amount of overlap for a truncated end adapter, by default None
        """
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
        # this is the adapter at the start of the read
        self.adapter_start_forward = Adapter(
            self.start_barcode,
            maximal_number_of_errors=self.maximal_errors_start,
            index_adapter_end=True,
            minimal_overlap=self.minimal_overlap_start,
            find_right_most_occurence=False,
        )
        # this is the reverse of the adapter at the end of the read
        self.adapter_end_reverse = Adapter(
            reverse_complement(self.end_barcode),
            maximal_number_of_errors=self.maximal_errors_end,
            index_adapter_end=False,
            minimal_overlap=self.minimal_overlap_end,
            find_right_most_occurence=True,
        )
        # this is the adapter at the end of the read
        self.adapter_end_forward = Adapter(
            self.end_barcode,
            maximal_number_of_errors=self.maximal_errors_end,
            index_adapter_end=True,
            minimal_overlap=self.minimal_overlap_end,
            find_right_most_occurence=False,
        )
        # this is the reverse adapter at the start of the read
        self.adapter_start_reverse = Adapter(
            reverse_complement(self.start_barcode),
            maximal_number_of_errors=self.maximal_errors_start,
            index_adapter_end=False,
            minimal_overlap=self.minimal_overlap_start,
            find_right_most_occurence=True,
        )

    def match_and_trim(self, fragment: Fragment) -> Union[Fragment, Literal[False]]:
        """
        Trims a fragment based on start and end adapters.

        Parameters
        ----------
        fragment : Fragment
            The Fragment to be trimmed/demultiplexed.

        Returns
        -------
        Union[Fragment, Literal[False]]
            A fragment that has the start adapter removed and, if present also
            the end adapter. The sequence will be oriented in forward direction
            from the perspective of the start adapter.
            If the fragment is empty after trimming, False is returned.
            If the start adapter is not found, False is returned.
        """
        start_in_read = self.adapter_start_forward.locate(fragment.Read1.Sequence)
        start_reverse_in_read = self.adapter_start_reverse.locate(
            fragment.Read1.Sequence
        )
        if start_in_read is None and start_reverse_in_read is None:
            # start adapter nowhere to be found, discard
            return False
        elif start_in_read is not None and start_reverse_in_read is None:
            # start adapter in read, nothing to be done
            pass  # pragma: no cover
        elif start_in_read is None and start_reverse_in_read is not None:
            # reverse start adapter in read, turn the thing around
            start_in_read = start_reverse_in_read  # we want to trim the adapter
            fragment = Fragment(
                fragment.Read1.reverse()
            )  # we want to turn around the sequence so all are oriented right
        else:  # start_in_read and start_reverse_in_read
            if start_in_read > start_reverse_in_read:
                # start reverse adapter in read is more likely
                start_in_read = start_reverse_in_read  # we want to trim the adapter
                fragment = Fragment(
                    fragment.Read1.reverse()
                )  # we want to turn around the sequence so all are oriented right
        # now read contains start adapter forward
        fragment.Read1 = self.trim_read_front(
            fragment.Read1, start_in_read + self.trim_after_start
        )  # trim start adapter in first read
        end_in_read = self.adapter_end_forward.locate(fragment.Read1.Sequence)
        if end_in_read is not None:
            # end read found, trim it
            fragment.Read1 = self.trim_read_back(
                fragment.Read1, end_in_read + self.trim_before_end
            )
        # accept fragment, both adapters are found
        # trim reverse adapters at the end of reads
        else:
            end_reverse_in_read = self.adapter_end_reverse.locate(
                fragment.Read1.Sequence
            )
            if end_reverse_in_read is not None:
                fragment.Read1 = self.trim_read_back(
                    fragment.Read1, end_reverse_in_read + self.trim_before_end
                )
            else:
                # no end adapter found, nothing to be done
                pass
        if len(fragment.Read1.Sequence) == 0:
            # we killed it, discard
            return False
        return fragment


class SE_Trim_On_Start_Trim_After_X_BP(DemultiplexStrategy):
    """Only the start barcode is relevant"""

    def __init__(
        self,
        start_barcode: str,
        trim_after_start: int = 20,
        maximal_errors_start: int = 0,
        min_length: int = 0,
        minimal_overlap_start: Optional[int] = None,
        **kwargs
    ):
        """
        Demultiplex and Trim strategy based on start adapter.

        This strategy assumes that there is a single end sequuencing and only
        the first adapter is relevant for trimming.
        The read is trimmed X bp after the start adapter. If the start adapter is not found, the
        The end adapter is only used for trimming. If the start adapter is not found, the
        fragment is discarded.

        Parameters
        ----------
        start_barcode : str
            Barcode at the start of read
        trim_after_start : int, optional
            an offset position to trim after the start adapter, by default 0
        maximal_errors_start : int, optional
            max number of errors in start adapter allowed for a match, by default 0
        min_length : Optional[int], optional
            the minimal remaining read length, by default 0
        minimal_overlap_start : Optional[int], optional
            the amount of overlap for a truncated start adapter, by default None
        """
        self.start_barcode = start_barcode
        self.trim_after_start = trim_after_start
        self.maximal_errors_start = maximal_errors_start
        self.minimal_overlap_start = minimal_overlap_start
        self.min_length = min_length
        if minimal_overlap_start is None:
            self.minimal_overlap_start = len(self.start_barcode)
        self._init_adapter()

    def get_parameters(self) -> List:
        return [
            self.start_barcode,
            self.trim_after_start,
            self.maximal_errors_start,
            self.minimal_overlap_start,
        ]

    def _init_adapter(self) -> None:
        # this is the adapter at the start of the read
        self.adapter_start_forward = Adapter(
            self.start_barcode,
            maximal_number_of_errors=self.maximal_errors_start,
            index_adapter_end=True,
            minimal_overlap=self.minimal_overlap_start,
            find_right_most_occurence=False,
        )

    def match_and_trim(self, fragment: Fragment) -> Union[Fragment, Literal[False]]:
        """
        Trims a fragment based on start and end adapters.

        Parameters
        ----------
        fragment : Fragment
            The Fragment to be trimmed/demultiplexed.

        Returns
        -------
        Union[Fragment, Literal[False]]
            A fragment that has the start adapter removed and, if present also
            the end adapter. The sequence will be oriented in forward direction
            from the perspective of the start adapter.
            If the fragment is empty after trimming, False is returned.
            If the start adapter is not found, False is returned.
        """
        start_in_read = self.adapter_start_forward.locate(fragment.Read1.Sequence)
        if start_in_read is None:
            # start adapter nowhere to be found, discard
            return False
        else:  # start_in_read and start_reverse_in_read
            print(fragment.Read1.Sequence)
            print(self.trim_after_start)
            print(start_in_read)
            fragment.Read1 = self.trim_read_back(
                fragment.Read1, start_in_read + self.trim_after_start
            )  # trim after the adapter
            print(fragment.Read1.Sequence)
            fragment.Read1 = self.trim_read_front(
                fragment.Read1, start_in_read
            )  # trim start adapter in first read
            print(fragment.Read1.Sequence)
            print(len(fragment.Read1.Sequence), self.min_length)
            if len(fragment.Read1.Sequence) < self.min_length:
                # we killed it, discard
                return False
            return fragment


class PE_Decide_On_Start_End_Trim_Start_End(DemultiplexStrategy):
    """The start barcode and end barcode are is relevant"""

    def __init__(
        self,
        start_barcode: str,
        end_barcode: str,
        trim_after_start: int = 0,
        trim_before_end: int = 0,
        trim_length_r1: int = 0,
        trim_length_r2: int = 0,
        maximal_errors_start: int = 0,
        maximal_errors_end: int = 0,
        minimal_overlap_start: Optional[int] = None,
        minimal_overlap_end: Optional[int] = None,
        **kwargs
    ):
        self.start_barcode = start_barcode
        self.end_barcode = end_barcode
        self.trim_after_start = trim_after_start
        self.trim_before_end = trim_before_end
        self.trim_length_r1 = trim_length_r1
        self.trim_length_r2 = trim_length_r2
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
            self.trim_length_r2,
            self.trim_length_r1,
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
        if "M03491:45:000000000-GT86M:1:1101:14198:1505" in fragment.Read1.Name:
            print(
                "Entered",
                self.start_barcode,
                self.end_barcode,
                start_in_r1,
                start_in_r2,
            )
            print(fragment.Read1.Sequence)
            print(fragment.Read2.Sequence)
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
        #        raise ValueError()
        end_in_r2 = self.adapter_end_forward.locate(fragment.Read2.Sequence)
        if end_in_r2 is None:
            # end adapter not found, discard
            return False
        fragment.Read2 = self.trim_read_front(
            fragment.Read2, end_in_r2 + self.trim_before_end
        )
        # accept fragment, both adapters are found
        # trim reverse adapters at the end of reads
        if "M03491:45:000000000-GT86M:1:1101:14198:1505" in fragment.Read1.Name:
            print("Here")

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
        if "M03491:45:000000000-GT86M:1:1101:14198:1505" in fragment.Read1.Name:
            print("not discarded")

        return self.trim_fragment_to_length(fragment)
