"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) CHER_WEI_YUAN, 29DEC2021 
License     : MIT_LICENSE 
Maintainer  : E0031403@U.NUS.EDU 
Portability : POSIX
This program reads the output of GeneScan in csv format, remove
peaks with small or relatively small areas, and calculates, for 
each sample, the percentage of the total area that each peak covers.

The hard filtering procedure here has low threshold (in ambiguous
situations, many peaks will not be removed) and the decision to keep
or remove peaks is left to the user.
"""

PROGRAM_NAME = "genescanner"
PROGRAM_VERSION = "1.0.1"

EXIT_COLUMN_HEADER_ERROR = 1
DIRECTORY_MISSING_ERROR = 2
PERMISSION_ERROR = 3

import sys
from os import mkdir
from os import path
from logging import basicConfig
from logging import info
from logging import INFO
from logging import error
from argparse import ArgumentParser

import pandas as pd
from itertools import chain
from seaborn import FacetGrid
from matplotlib.pyplot import bar

def parse_args():
    """
    Parse command line arguments.
    """
    description = "Reads the output of GeneScan in csv format, remove peaks with small area, and calculates, for each sample, the percentage of the total area that each peak covers."
    epilog = "Example usage: genescanner \
        --outdir mnt/c/mySamples/output \
        --prefix mySamples \
        /mnt/c/mySamples/mySamples.csv"
    parser = ArgumentParser(description=description,
                            epilog=epilog)
    parser.add_argument('input',
                        type=str,
                        help="Input GeneScan datasheet from a single run in CSV format. A single run can contain many samples but they share the same capillary injection")
    parser.add_argument('--exon_df',
                        type=str,
                        default = None,
                        help='Datasheet containing exon sizes per sample in CSV format')
    parser.add_argument('--outdir',
                        type=str,
                        default='stdout',
                        help='Name of output directory')
    parser.add_argument('--prefix',
                        type=str,
                        default='out',
                        help='Prefix name of output files')
    parser.add_argument('--version',
                        action='version',
                        version=str(PROGRAM_VERSION))
    parser.add_argument('--peak_gap',
                        type=float,
                        default = 1.7,
                        help='DEFAULT = 1.7. A pair of peaks within peak_gap of each other will be processed to give one peak')
    parser.add_argument('--cluster_size',
                        type=int,
                        default = 3,
                        help='DEFAULT = 3. The maximum number of peaks within peak_gap of each other that will be processed together. Only one peak with largest area will remain.')
    parser.add_argument('--filter',
                        metavar='FILTER',
                        type=float,
                        default = 0.0,
                        help='DEFAULT = 0.0. Float. Remove all peaks with percentage area lower than filter. Percentage area refers to the area of the peak over the area of all peaks of the same sample.') 
    parser.add_argument('--shift_range',
                        type=list,
                        default=[-500, 500, 0.5],
                        help='Range of Error landscape')
    return parser.parse_args()

def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)

def loadDf(input_file):
    """
    Parameters
    ----------
    input_file : TYPE String
        DESCRIPTION. Directory and name of GeneScan output sheet in csv.
                     (e.g. ./input/data.csv)

    Returns
    -------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    """
    # Load dataframe
    df = pd.read_csv(input_file)
    
    # Sort dataframe by sample name and size
    df = df.sort_values(by=['Sample File Name', 'Size']).\
        reset_index(drop = True, inplace = False)
        
    # Clean sample names by removing leading and trailing white space
    df["Sample File Name"] = df["Sample File Name"].str.rstrip().str.lstrip()
    
    # Clean df by making sure dtype is correct
    df = df.astype({'Sample File Name': str,
                    'Size': float, 
                    'Height': float, 
                    'Area': float})
    
    # Check if df column names are correct
    expected = ['Sample File Name',
                'Size','Height','Area']
    for i in expected:
        if i not in list(df.columns):
            exit_with_error(f"Unexpected column header detected. Rename columns to {expected} and retry", 
                            EXIT_COLUMN_HEADER_ERROR)    
    return df
       
def getSampleNames(df):
    """
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.

    Returns
    -------
    sample_names : TYPE List
        DESCRIPTION. List of unique sample names from df.
    """
    sample_names = list(set(df.loc[:,"Sample File Name"]))
    return sample_names

def findPeakCluster(index, build_list, df, peak_gap):
    """
    Recursively finds members of a peak cluster starting from the peak with
    the smallest size.
    
    Parameters
    ----------
    index : TYPE Integer
        DESCRIPTION. The index of df that corresponds to the
                     rows (i.e. peaks) that are clustered (within peak_gap
                     of each other) and awaiting to be processed to give 
                     fewer peaks.
    build_list : TYPE List
        DESCRIPTION. List of index of peaks in peak clusters
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    peak_gap : TYPE Integer
        DESCRIPTION. User-supplied. A pair of peaks within peak_gap of
                     each other will be processed to give one peak.

    Returns
    -------
    TYPE List
    A list of index corresponding to peaks in a peak cluster. 
    """    
    # Return build_list if we reach the end of dataframe
    if index == max(df.index):
        return list(set(build_list))
    # Stop recursion when next peak is not within peak_gap of current peak
    elif df.loc[index + 1, "Size"] - df.loc[index, "Size"] > peak_gap:
        return list(set(build_list))
    # Recursion to next peak
    else:
        build_list += [index, index + 1]
        return findPeakCluster(index + 1, build_list, df, peak_gap)

def findMountainRanges(df, sample_names, peak_gap):
    """
    Finds a list of peak_cluster (each peak cluster is a list of index
    corresponding to peaks within peak_gap of each other) called
    mountain_ranges (i.e. the collection of all peak clusters in df).
    
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    sample_names : TYPE List
        DESCRIPTION. List of unique sample names from df.

    Returns
    -------
    mountain_ranges : TYPE List
        DESCRIPTION. A master list containing list of indexes of 
                     continuous peaks within peak_gap of each other.
    """
    mountain_ranges = []
    # Loop through samples
    for i in sample_names:
        dfSample = df.loc[df.iloc[:,1] == i, :]
        # Loop through rows in subset df
        # and run findPeakCluster
        for index, row in dfSample.iterrows():
            # Skip rows checked by findPeakCluster
            if index in chain(*mountain_ranges):
                pass
            else:
                mountain_ranges += [findPeakCluster(index, 
                                                      [], 
                                                      dfSample, 
                                                      peak_gap)]
    # Remove duplicate empty nested lists
    mountain_ranges = [x for x in mountain_ranges if x != []]
    return mountain_ranges
            
def cleanMountainRanges(df, mountain_ranges, cluster_size):
    """
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    clean_mountain_ranges : TYPE List
        DESCRIPTION. A master list containing list of indexes of 
                     continuous peaks within peak_gap of each other

    Returns
    -------
    remove: TYPE List
        DESCRIPTION. List of index to remove from df
    """
    remove = []
    for mountains in mountain_ranges:
        if len(mountains) == 2:
            a = df.loc[mountains[0], "Area"]
            b = df.loc[mountains[1], "Area"]
            if a > b:
                remove += [[mountains[1]]]
            elif b > a:
                remove += [[mountains[0]]]
            # Keep both peaks if they are equal
            elif a == b:
                pass
        elif len(mountains) > 2 and len(mountains) <= cluster_size:
            # Keep list of all areas of peaks in cluster
            lst = []
            for i in range(len(mountains)):
                lst += [df.loc[mountains[i], "Area"]]
            # Remove index corresponding to largest area
            # and add remaining index to remove
            mountains.pop(lst.index(max(lst)))
            remove += [ mountains ]
            
    remove = list(chain(*remove))
    return remove

def RemoveArtefacts(df, remove):
    """
    Create df for output to user
    
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    remove : TYPE List
        DESCRIPTION. Flattened list of index to remove from df

    Returns
    -------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet with dirty peaks
                     removed.
    """
    return df.drop(labels = remove, axis = 0, inplace = False)

def labelArtefacts(df, remove):
    """
    Create df for plotting
    
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    remove : TYPE List
        DESCRIPTION. Flattened list of index to remove from df

    Returns
    -------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet with dirty peaks
                     to be removed marked as "Removed".
    """
    df["Status"] = "Kept"
    if remove == []:
        pass
    else:
        for i in remove:
            df.loc[i, "Status"] = "Removed"
    return df
        
def AddPercentage(processed_df, sample_names):
    """
    Parameters
    ----------
    processed_df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet from RemoveArtefacts.
    sample_names : TYPE List
        DESCRIPTION. List of unique sample names from df.

    Returns
    -------
    df_out : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned and processed GeneScan datasheet.
    """
    data_list = []
    # Loop through samples
    for i in sample_names:
        dfSample = processed_df.loc[processed_df.iloc[:,1] == i, :]
        
        # Calculate percentage per peak (row)
        # and keep only useful rows
        total_area = sum(dfSample.loc[:,"Area"])
        for index, row in dfSample.iterrows():
            name = row["Sample File Name"]
            size = row["Size"]
            height = row["Height"]
            area = row["Area"]
            percentage = int(round((row["Area"]/ total_area) * 100, 0))
            data_list.append([name, size, height, area, percentage])

    # Make output dataframe
    df_out = pd.DataFrame(data_list, columns = \
                          ["Sample File Name", 
                           "Size", 
                           "Height", 
                           "Area", 
                           "Percentage"])
    return df_out

def filterAreaPercent(df_out, filter_threshold):
    return df_out.query(f"Percentage >= {filter_threshold}").\
        sort_values(by=['Sample File Name', 'Size'], ascending = True)

def plot(df_before, df_after, prefix, outdir):
    # Sort dataframe by sample name and size
    df_before = df_before.sort_values(by=['Sample File Name', 'Size'])
    df_after = df_after.sort_values(by=['Sample File Name', 'Size'])
        
    # Add new column to differentiate before and after
    df_before["Processed"] = "Before"
    df_after["Processed"] = "After"
    
    # Combine df 
    df_combine = pd.concat([df_before, df_after]).\
        reset_index(drop = True, inplace = False)

    # Change "Size" to categorical variable
    df_combine.loc[:,"Size"] = df_combine.loc[:,"Size"].astype('category')
    
    # Grid plot
    grid = FacetGrid(df_combine, 
                     hue='Status',
                     hue_order=["Kept", "Removed"],
                     row='Sample File Name',
                     col='Processed',
                     sharey=False, 
                     sharex=False,
                     legend_out=True)
    grid.map(bar, 
             'Size', 
             'Area', 
             width = 10)
    
    grid.set_titles('{row_name} ({col_name})')
    grid.figure.subplots_adjust(wspace=1.2)
    grid.set_axis_labels("Size (bp)", "Area of peak")
    grid.add_legend()
    grid.savefig(f"{outdir}/{prefix}.png")

def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        basicConfig(filename=log_filename,
                            level=INFO,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z",
                            force=True)
        info('program started')
        info('command line: %s', ' '.join(sys.argv))

def main():
    # Get user-defined parameters
    args = parse_args()
    input_file = args.input
    prefix = args.prefix
    outdir = args.outdir
    peak_gap = args.peak_gap
    filter_threshold = args.filter
    cluster_size = args.cluster_size
    exon_df = args.exon_dif
    shift_range = args.shift_range
    shift_start = shift_range[0]
    shift_end  = shift_range[1]
    shift_step = shift_range[2]
    
    # Ensure output file is not opened
    try:
        temp = pd.read_csv(f"{outdir}/{prefix}.csv")
        temp.to_csv(f"{outdir}/{prefix}.csv", index = False)
    except PermissionError:
        exit_with_error(f"Please close {outdir}/{prefix}.csv and try again", 
                        PERMISSION_ERROR)
    except FileNotFoundError:
        pass
                        
    # Check if output directory is available
    try:
        if not path.exists(f"{outdir}"):
            mkdir(f"{outdir}")
    except:
        exit_with_error("Output directory does not exist", 
                        DIRECTORY_MISSING_ERROR)
        
    # Initiate logging
    log_filename = f"{prefix}_log.txt"
    init_logging(f"{outdir}/{log_filename}")
    
    # Load and clean input
    df = loadDf(input_file)
    
    # Get all sample names
    sample_names = getSampleNames(df)
    
    # Find all peak clusters
    # Peak clusters collectively form mountain ranges
    mountain_ranges = findMountainRanges(df, sample_names, peak_gap)
    
    # For every peak cluster of cluster_size, pick peak with largest area
    remove = cleanMountainRanges(df, mountain_ranges, cluster_size)
    processed_df = RemoveArtefacts(df, remove)
    
    # Calculate percentage area of each peak across total area of its sample
    processed_df = AddPercentage(processed_df, sample_names)
    
    # Remove percentages below filter threshold
    processed_df = filterAreaPercent(processed_df, filter_threshold)
    
    # Save output df to csv
    processed_df.to_csv(f"{outdir}/{prefix}_cleanPeaks.csv", 
                        index = False)

    # Plot 
    plot(labelArtefacts(df, remove), 
         labelArtefacts(processed_df, []), 
         prefix, 
         outdir)
    
    # Get dict of shift against Error
    Error_dict = drawErrorLandscape(processed_df, 
                                  exon_df, 
                                  sample_names,
                                  shift_start, shift_end, shift_step)
    
    # Calculate shift
    shift = findShift(Error_dict)
    
    # Assign exon combinations to processed_df
    exon_processed_df = processDF(processed_df, 
                                  exon_df, 
                                  sample_names, 
                                  shift, 
                                  "AssignExonCombinations")
    
    # Save processed_df to csv
    exon_processed_df.to_csv(f"{outdir}/{prefix}_AssignedExons.csv", 
                             index = False)
    
    # Inform completion
    info(f"Output saved to {outdir}")
    info("Complete")
    print(f"Complete. Output saved to {outdir}")

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()

##### Development for isoform assignment to peaks
# Note all sizes are in base pair (bp)
# Log shift and command inputs
# Plot error-shift line graph as visual quality control
# Position of df is important: pos 1 is "Sample File Name"
# Need to read exon_df

from itertools import combinations
from bisect import bisect_left
from numpy import arange

def findClosestError(Error_list, shift):
    Error_list = sorted(Error_list)
    pos = bisect_left(Error_list, shift)
    if pos == 0:
        return Error_list[0]
    if pos == len(Error_list):
        return Error_list[-1]
    before = Error_list[pos - 1]
    after = Error_list[pos]
    if (shift - before) > (after - shift):
        return after
    else:
        return before

def findAllExonCombinations(exon_sizes, size):
    # Finds all exon combinations 
    # and calculates Error for each combination
    # Returns dictionary where keys are Error values
    # and values are list of exon_combinations: 
    #      [ [exon set 1], [exon set 2]... ]
    perm_list = []
    Error_list = []
    
    # Find all exon combinations
    for i in range(len(exon_sizes)):
        perm_list += [ list(x) for x in combinations(exon_sizes, i+1) ]
    
    # Find Error, i.e. the difference between peak size and sum of exon size
    Error_list = [ sum(x)-size for x in perm_list ]
    
    # Make dictionary of Error: Exon combination
    errorExonMap = {}
    for i in range(len(Error_list)):
        try:
            errorExonMap[Error_list[i]] += [ perm_list[i] ]
        except KeyError:
            errorExonMap[Error_list[i]] = [ perm_list[i] ]
    
    return errorExonMap

def SelectExonCombinations(errorExonMap, size):
    # errorExonMap = { ErrorA = [ [exon set 1], [exon set 2]... ],
    #                  ErrorB = [ [exon set 3], [exon set 4]... ] }
    
    # Find exon combination closest to peak_size closest to size,
    # whether Error is positive or negative 
    closest_Error = findClosestError(list(errorExonMap.keys()), size)
    return { closest_Error: errorExonMap[closest_Error]}

def processDF(processed_df, 
              exon_df, 
              sample_names,
              shift,
              function):
    """
    exon df structure
    Sample File Name    Exon    Exon Size
    """

    processed_df["Exon Combination"] = None
    processed_df["Error"] = None
    
    # Loop through processed_df sample by sample
    for i in sample_names:
        dfSample = processed_df.loc[processed_df.iloc[:, 0] == i, :]
        
        # Get all possible exon sizes for each sample
        dfExonSample = exon_df.loc[exon_df.iloc[:, 0] == i, :]
        exon_sizes = list(dfExonSample["Exon Size"])
        
        # Assign Error, exon combinations to each peak in sample 
        for index, row in dfSample.iterrows():
            # Get exon combination with error closest to zero
            peak_size = dfSample.loc[index, "Size"]
            errorExonMap = findAllExonCombinations(exon_sizes, peak_size + shift) 
            out = SelectExonCombinations(errorExonMap, 0) 
            # Assign exon combinations
            for key, value in out.items():
                try:
                    processed_df.at[index, "Exon Combination"] += [value]
                except TypeError:
                    processed_df.at[index, "Exon Combination"] = [value]
                    processed_df.at[index, "Error"] = round(key, 2)  
                    
    if function == "findLowestError":
        return sum(processed_df["Error"])
    
    elif function == "AssignExonCombinations":
        return processed_df

def drawErrorLandscape(processed_df, 
                       exon_df, 
                       sample_names,
                       shift_start, shift_end, shift_step):
    # Create new df
    Error_dict = {}
    
    # Initialize parameters
    count = 0
    total_iterations = (shift_end - shift_start)/shift_step
    
    # Add Shift, Error to Error_dict
    for shift in arange(shift_start, shift_end, shift_step):
        Error = processDF(processed_df, 
                          exon_df, 
                          sample_names,
                          shift,
                          "findLowestError")
        Error_dict[Error] = shift
        count += 1
        print(f"Iteration {count}/ {total_iterations} completed")
        
    return Error_dict

def findShift(Error_dict):
    Error_list = list(Error_dict.keys())
    Error_near_zero = findClosestError(Error_list, 0)
    shift = Error_dict[Error_near_zero]
    return shift

def translateSizeToExon(processed_df, exon_df):
    for index, row in processed_df.iterrows():
        

# TEST
import pandas as pd
exon_df = pd.read_csv("test_exon_input.csv")
df = loadDF("input_basic_test.csv")
processDF(processed_df, 
              exon_df, 
              sample_names,
              shift,
              "findLowestError")
exon_processed_df1 = processDF(processed_df, 
                                  exon_df, 
                                  sample_names, 
                                  -273, 
                                  "AssignExonCombinations")
    
exon_processed_df2 = processDF(processed_df, 
                                  exon_df, 
                                  sample_names, 
                                  0, 
                                  "AssignExonCombinations")    
