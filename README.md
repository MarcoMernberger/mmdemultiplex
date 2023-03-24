# mmdemultiplex

Demultiplex fastqs

## Description

A flexible barcode demultiplexer library for sequencing projects written in python.

Tested using python 3.8 and linux. Requirements are specified in setup.cfg. Does
not require special hardware.

Use 'pip install .' to install. Installation is expected to take less than a
minute.

## Usage

```python
import demultiplex
import pandas as pd
import pypipegraph as pp


ppg.new_pipegraph()

mbf_align_sample = mbf_align.raw.Sample(
	"combined_fastq",
	"incoming/combined.fastq.gz",
	reverse_reads = False,
	pairing="paired",
	vid=None)

dm = demultiplex.Demultiplexer(
	mbf_align_sample,
	pd.DataFrame(
        {
            "key": ["Sample1", "Sample2", "Sample3"],
            "start_barcode": ["CTGGCA", "GGTCGA", "ACTGTG"],
            "end_barcode": ["CAAAAG", "GCCACC", "AAGTGC"],
        }
    ),
	output_dir = "output"
	)

ppg.run_pipegraph()
```

Expected output: 
One `.fastq` per sample in `./output`.


## Demo 
Please refer to tests/conftest.py for demonstration data,
and to the test cases for the expected output.

The tests should complet in a couple of minutes.

To use mmdemultiplex on your data, write a python script as above.
There is no CLI interface.

