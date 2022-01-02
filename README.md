# GeneScanner

# Overview 

GeneScanner reads a raw GeneScan datasheet, removes dirty peaks and calculates the percentage area of each peak among all peaks in the sample.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bionitio-team/bionitio/master/LICENSE).

# Installing

Clone this repository: 
```
git clone https://github.com/CherWeiYuan/GeneScanner.git
```

Move into the repository directory:
```
cd GeneScanner
```

Python 3 is required for this software.

Install inside a virtual environment:
```
pip install virtualenv
virtualenv gsenv -p `which python3`
source gsenv/bin/activate
pip install .
```

# Testing

## Unit tests

Test package is provided in GeneScanner/genescanner/test. You should run the test with the following commands:
```
cd GeneScanner/genescanner
python3 genescanner_test.py
```

GitHub template used here is provided by https://github.com/bionitio-team/bionitio

## Help message
```
$ genescanner -h
Initializing parameters
usage: genescanner [-h] [--outdir OUTDIR] [--prefix PREFIX] [--version] [--peak_gap PEAK_GAP]
                   [--cluster_size CLUSTER_SIZE] [--filter FILTER]
                   input

Reads the output of GeneScan in csv format, remove peaks with small or relatively small areas, and calculates, for
each sample, the percentage of the total area that each peak covers.

positional arguments:
  input                 Input GeneScan datasheet in CSV format

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR       Name of output directory
  --prefix PREFIX       Prefix name of output files
  --version             show program's version number and exit
  --peak_gap PEAK_GAP   DEFAULT = 1.0. A pair of peaks within peak_gap of each other will be processed to give one
                        peak
  --cluster_size CLUSTER_SIZE
                        DEFAULT = 3. The maximum number of peaks within peak_gap of each other that will be processed
                        together. Only one peak with largest area will remain.
  --filter FILTER       DEFAULT = 1.0. Float. Remove all peaks with percentage area lower than filter. Percentage area
                        refers to the area of the peak over the area of all peaks of the same sample.

Example usage: python3 genescanner.py --outdir /mnt/c/genescan/out --prefix test_2 /mnt/c/genescan/in/input.csv
```

