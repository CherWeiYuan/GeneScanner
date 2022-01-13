'''
Unit tests for genescanner
Usage: python3 genescanner_test.py
'''

from genescanner import *
from os import system
from os import remove
from os import path
from glob import glob
import sys

# Set-up temp folder
if path.exists("./test/temp"):
    files=glob('./test/temp/*')
    for filename in files:
        remove(filename)    
if not path.exists("./test/temp"):
    mkdir("./test/temp")

def test_wholeProgram(expected, test_type):
    "Test output of entire program"
    command = system("python3 genescanner.py \
                        --outdir ./test/temp \
                        --prefix WholeProgramTest \
                        --peak_gap 1.7 \
                        --filter 0 \
                        --cluster_size 3 \
                        ./test/test_input/input_basic_test.csv")
    try:
        result = pd.read_csv("./test/temp/WholeProgramTest.csv")
    except PermissionError:
        sys.exit()
    
    expected = pd.read_csv("./test/test_output/output_basic_test.csv")
    pd.testing.assert_frame_equal(result, expected, 
                                  check_dtype=False,
                                  check_index_type = False,)
    print(f"Passed: Whole Level Test 1")

# Tests
test1 = test_wholeProgram("./test/test_output/output_basic_test.csv", "Whole Level Test 1")

