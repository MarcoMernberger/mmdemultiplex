from mmdemultiplex.util import create_fragment
from mmdemultiplex.filter import SequenceFilter
from mmdemultiplex.adapters import Adapter, WHERE_START, WHERE_END, AdapterMatch


def test_SequenceFilter_init():
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], maximal_number_of_errors=2)
    assert sfilter.sequences == ["AGGACATACCAG"]
    assert isinstance(sfilter.adapters, list)
    assert isinstance(sfilter.adapters[0], Adapter)


def test_SequenceFilter_match():
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], maximal_number_of_errors=0)
    example = "GGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)


def test_SequenceFilter_match_error():
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], maximal_number_of_errors=2)
    example1err = "GGGGAGGAGGATGGGGAGTAGGACATACGAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example1err])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example2err = "GGGGAGGAGGATGGGGAGTAGGACATAGGAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example2err])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example3err = "GGGGAGGAGGATGGGGAGTAGGACATTGGAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example3err])
    matcher = sfilter.get_filter_callable()
    assert not matcher(fragment)


def test_SequenceFilter_match_indel():
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], maximal_number_of_errors=2)
    example = "GGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example2del = "GGGGAGGAGGATGGGGAGTAGGACATAAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example2del])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example2in = "GGGGAGGAGGATGGGGAGTAGGACATACCCCAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example2in])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example2in = "GGGGAGGAGGATGGGGAGTAGGACATACCCCCAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example2in])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)


def test_not_accepting_reads_that_mis():
    example = "TCGCTCTGGATTGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGTGGAGGTTGGGGAGTCCAGGGAAGCTGTCCCTCACTGTTGATTTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTTGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACAT"
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], maximal_number_of_errors=2)
    fragment = create_fragment([example])
    matcher = sfilter.get_filter_callable()
    assert not matcher(fragment)


def test_accepting_reads_that_mis():
    sequence = "AGGACATACCAG"
    sfilter = SequenceFilter([sequence], present=False, maximal_number_of_errors=0)
    example = "TCGCTCTGGATTGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGTGGAGGTTGGGGAGTCCAGGGAAGCTGTCCCTCACTGTTGATTTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTTGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACAT"
    fragment = create_fragment([example])
    matcher = sfilter.get_filter_callable()
    assert matcher(fragment)
    example = "GGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTT"
    fragment = create_fragment([example])
    matcher = sfilter.get_filter_callable()
    assert not matcher(fragment)
