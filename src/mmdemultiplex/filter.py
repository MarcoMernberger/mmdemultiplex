from typing import List
from .util import Fragment, reverse_complement
from .adapters import Adapter


class SequenceFilter:

    def __init__(
        self,
        adapter_sequences: List[str],
        present: bool = True,
        maximal_number_of_errors: int = 0,
    ):
        """
        A class for creating a sequence filter for reads processed by the demultiplexer.

        This class accepts a number of sequences that should be present or
         absent in a Read for it to be accepted. If we want to exclude these
         positive reads or exclude negative reads, we can set the present
         parameter accordingly.

        Parameters
        ----------
        adapter_sequences : List[str]
            sub sequences in reads to filter by.
        present . bool, optional
            if we should keep positive reads or discard them, by default True
        maximal_number_of_errors : int, optional
            max errors for cutadapt error matching, by default 0
        """
        self.sequences = adapter_sequences
        self.adapters = []
        self.present = present
        for sequence in self.sequences:
            adapter = Adapter(
                adapter_sequence=sequence,
                maximal_number_of_errors=maximal_number_of_errors,
                index_adapter_end=True,
                minimal_overlap=None,
                find_right_most_occurence=False,
            )
            self.adapters.append(adapter)

    def match(self, fragment: Fragment):
        """
        match is a function that returns True if a sequence was found in a read.

        Parameters
        ----------
        fragment : Fragment
            the fragment to be tested.
        """
        found = False
        for adapter in self.adapters:
            for read in fragment:
                found = (
                    found
                    or (adapter.locate(read.Sequence) is not None)
                    or (adapter.locate(reverse_complement(read.Sequence)) is not None)
                )
            if found:
                break
        return found

    def get_filter_callable(self):

        def accept_positive(fragment: Fragment):
            return self.match(fragment)

        def accept_negative(fragment: Fragment):
            return not self.match(fragment)

        return accept_positive if self.present else accept_negative
