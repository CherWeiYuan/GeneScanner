# GeneScanner

# Overview 

GeneScanner reads a raw GeneScan datasheet, removes dirty peaks, calculates the percentage area of each peak among all peaks in the sample and output a cleaned csv datasheet.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bionitio-team/bionitio/master/LICENSE).

# Installing

Python 3 is required for this software.

Clone this repository: 
```
git clone https://github.com/CherWeiYuan/GeneScanner.git
```

Move into the repository directory:
```
cd GeneScanner
```

Install inside a virtual environment:

Install virtualenv 
```
pip install virtualenv
```

Create and start virtual environment named gsenv
```
virtualenv gsenv -p `which python3`
source gsenv/bin/activate
```

Install GeneScanner inside environment
```
pip install .
```

Deactivate environment once installation is done
```
deactivate
```

# Testing

## Unit tests

Test package is provided in GeneScanner/genescanner/test. You should run the test with the following commands:
```
source gsenv/bin/activate
cd GeneScanner/genescanner
python3 genescanner_test.py
deactivate
```

# Execute

## Quick start
```
cd GeneScanner/genescanner
source gsenv/bin/activate
genescanner \
  --outdir mnt/c/mySamples/output \
  --prefix mySamples \
  /mnt/c/mySamples/mySamples.csv
deactivate
```

## Help message
```
$ genescanner -h
usage: genescanner [-h] [--outdir OUTDIR] [--prefix PREFIX] [--version] [--peak_gap PEAK_GAP]
                   [--cluster_size CLUSTER_SIZE] [--filter FILTER]
                   input

Reads the output of GeneScan in csv format, remove peaks with small area, and calculates, for each sample, the
percentage of the total area that each peak covers.

positional arguments:
  input                 Input GeneScan datasheet in CSV format

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR       Name of output directory
  --prefix PREFIX       Prefix name of output files
  --version             show program's version number and exit
  --peak_gap PEAK_GAP   DEFAULT = 1.7. A pair of peaks within peak_gap of each other will be processed to give one
                        peak
  --cluster_size CLUSTER_SIZE
                        DEFAULT = 3. The maximum number of peaks within peak_gap of each other that will be processed
                        together. Only one peak with largest area will remain.
  --filter FILTER       DEFAULT = 0.0. Float. Remove all peaks with percentage area lower than filter. Percentage area
                        refers to the area of the peak over the area of all peaks of the same sample.

Example usage: genescanner --outdir mnt/c/mySamples/output --prefix mySamples /mnt/c/mySamples/mySamples.csv
```

GitHub template used here is provided by https://github.com/bionitio-team/bionitio
