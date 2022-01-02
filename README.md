# GeneScanner

# Overview 

GeneScanner reads a raw GeneScan datasheet, removes dirty peaks and calculates the percentage area of each peak among all peaks in the sample.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bionitio-team/bionitio/master/LICENSE).

# Installing

## Installing directly from source code

# Test
cd GeneScanner/genescanner
python3 genescanner_test.py

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
