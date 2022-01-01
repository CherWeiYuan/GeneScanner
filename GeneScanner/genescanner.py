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
PROGRAM_VERSION = 1.0

EXIT_COLUMN_HEADER_ERROR = 1

import sys
from os import mkdir
from os import path
from logging import basicConfig
from logging import info
from logging import DEBUG
from logging import error
from argparse import ArgumentParser
from argparse import BooleanOptionalAction

import pandas as pd
from itertools import chain
from seaborn import FacetGrid
from matplotlib.pyplot import bar

def parse_args():
    """
    Parse command line arguments.
    """
    description = "Reads the output of GeneScan in csv format, remove \
        peaks with small or relatively small areas, and calculates, for \
        each sample, the percentage of the total area that each peak covers."
    epilog = "Usage: python3 genescanner.py --outdir /mnt/c/genescan/out \
        --prefix test_2 \
        /mnt/c/genescan/in/input.csv"
    parser = ArgumentParser(description=description,
                            epilog=epilog)
    parser.add_argument('input',
                        type=str,
                        help='Input GeneScan datasheet in CSV format')
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
                        type=int,
                        default = 1,
                        help='DEFAULT = 1. A pair of peaks within peak_gap \
                            of each other will be processed to give one peak')
    parser.add_argument('--threshold',
                        type=float,
                        default = 0.2,
                        help='DEFAULT = 0.2. Values from 0.0 to 1.0. The \
                            fraction is the area of a peak within a peak \
                            cluster divided by the area of the largest peak \
                            in the same cluster. If the fraction of a peak \
                            is lower than the threshold, the peak will not \
                            appear in the final output.') 
    parser.add_argument('--filter',
                        metavar='FILTER',
                        type=float,
                        default = 1.0,
                        help='DEFAULT = 1.0. Float. Remove all peaks with \
                            percentage area lower than filter. Percentage \
                            area refers to the area of the peak over the \
                            area of all peaks of the same sample.') 
    parser.add_argument('--plot', 
                        default=False, 
                        action=BooleanOptionalAction,
                        help='DEFAULT = False. If True, make barplots of \
                            GeneScan data before and after removing peaks.')
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
    
    # Check if df column names are correct
    expected = ['Dye/Sample Peak','Sample File Name',
                'Marker','Allele','Size','Height','Area',
                'Data Point']
    if list(df.columns) != expected:
        exit_with_error(f"Unexpected column header detected. \
                        Rename columns to {expected} and retry", 
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
    if index == max(df.index):
        return list(set(build_list))
    elif df.loc[index + 1, "Size"] - df.loc[index, "Size"] > peak_gap:
        return list(set(build_list))
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

def resolvePeakCluster(df, mountain_ranges, threshold):
    """

    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    mountain_ranges : TYPE List
        DESCRIPTION. A master list containing list of indexes of 
                     continuous peaks within peak_gap of each other.
    threshold : TYPE Float
        DESCRIPTION. 0.0 to 1.0. The fraction is the area of a peak
                     within a peak cluster divided by the area of the
                     largest peak in the same cluster. If the fraction
                     of a peak is lower than the threshold, the peak
                     will not appear in the final output.

    Returns
    -------
    clean_mountain_ranges : TYPE List
        DESCRIPTION. A master list containing list of indexes of 
                     continuous peaks within peak_gap of each other
                     after the removal of peaks lower than threshold

    """
    # Remove dirty peaks within peak cluster if possible
    clean_mountain_ranges = []
    # Loop through each peak cluster
    for peak_cluster in mountain_ranges:
        area_list = []
        clean_peak_cluster = []
        # Examine clusters if number of peaks exceeds three
        if len(peak_cluster) >= 3:
            for i in peak_cluster:
                area_list += [ df.iloc[i, 6] ]
            max_area = max(area_list)
            for k in range(len(area_list)):
                # Keep index in peak_cluster if its peak
                # is above user-defined threshold
                if (area_list[k]/ max_area) >= threshold:
                    clean_peak_cluster += [peak_cluster[k]]
            clean_mountain_ranges += [ clean_peak_cluster ]
        else:
            clean_mountain_ranges += [ peak_cluster ]
            
    return clean_mountain_ranges
            
def cleanMountainRanges(df, clean_mountain_ranges):
    """
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    clean_mountain_ranges : TYPE List
        DESCRIPTION. A master list containing list of indexes of 
                     continuous peaks within peak_gap of each other
                     after the removal of peaks lower than threshold

    Returns
    -------
    remove: TYPE List
        DESCRIPTION. List of index to remove from df
    """
    remove = []
    for mountains in clean_mountain_ranges:
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
    remove = list(chain(*remove))
    return remove

def RemoveArtefacts(df, remove):
    return df.drop(labels = remove, axis = 0, inplace = False)

def labelArtefacts(df, remove):
    df["Status"] = "Kept"
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
            percentage = (row["Area"]/ total_area) * 100
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
    return df_out.query(f"Percentage >= {filter_threshold}")

def plot(df, file_prefix, subplot_prefix, plot_bool, outdir):
    # Sort dataframe by sample name and size
    df = df.sort_values(by=['Sample File Name', 'Size']).\
        reset_index(drop = True, inplace = False)
    if plot_bool == True:
        # Change "Size" to categorical variable
        plot_df = df
        plot_df.loc[:,"Size"] = plot_df.loc[:,"Size"].astype('category')
        if "Status" in plot_df.columns:
            # Make barplot for data before processing
            colors = {"Removed": 'r', "Kept": 'g'}
            grid = FacetGrid(plot_df, 
                             col='Sample File Name', 
                             col_wrap=5, 
                             sharey=False, 
                             sharex=False)
            grid.map(bar, 
                     'Size', 
                     'Area', 
                     width = 10, 
                     color = [colors[i] for i in plot_df['Status']])
            grid.set_titles('{col_name}')

            grid.savefig(f"{outdir}/{subplot_prefix}_{file_prefix}.png")
        else:
            # Make barplot for data after processing
            grid = FacetGrid(plot_df,
                             col='Sample File Name', 
                             col_wrap=5, 
                             sharey=False, 
                             sharex=False)
            grid.map(bar, 
                     'Size', 
                     'Area', 
                     width = 10)
            grid.set_titles('{col_name}')

            grid.savefig(f"{outdir}/{subplot_prefix}_{file_prefix}.png")
    else:
        pass

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
                            level=DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        info('program started')
        info('command line: %s', ' '.join(sys.argv))

def main():
    print("Initializing parameters")
    # Get user-defined parameters
    args = parse_args()
    input_file = args.input
    prefix = args.prefix
    outdir = args.outdir
    peak_gap = args.peak_gap
    threshold = args.threshold
    filter_threshold = args.filter
    plot_bool = args.plot
    info("User parameters")
    info("---------------")
    info(f"Input file: {input_file}")
    info(f"Peak gap: {peak_gap}")
    info(f"Threshold: {threshold}")

    # Check if output directory is available
    if not path.exists(f"{outdir}"):
        mkdir(f"{outdir}")
        
    # Initiate logging
    log_filename = f"{prefix}_log.txt"
    init_logging(f"{outdir}/{log_filename}")
    
    print("Loading and cleaning data")
    # Load and clean input
    df = loadDf(input_file)
    
    # Get all sample names
    sample_names = getSampleNames(df)
    
    print("Processing calculations")
    # Find all peak clusters
    # Peak clusters collectively form mountain ranges
    mountain_ranges = findMountainRanges(df, sample_names, peak_gap)
    
    # Remove small peaks in peak cluster
    clean_mountain_ranges = resolvePeakCluster(df, mountain_ranges, threshold)
    
    # For every peak pair, pick peak with larger area
    remove = cleanMountainRanges(df, clean_mountain_ranges)
    processed_df = RemoveArtefacts(df, remove)
    
    # Calculate percentage area of each peak across total area of its sample
    processed_df = AddPercentage(processed_df, sample_names)
    
    # Remove percentages below filter threshold
    processed_df = filterAreaPercent(processed_df, filter_threshold)
    
    # Save output df to csv
    processed_df.to_csv(f"{outdir}/{prefix}.csv")

    # Plot if plot_bool is True
    if plot_bool:
        print("Plotting")
    
    plot(labelArtefacts(df, remove), 
         prefix, 
         "before", 
         plot_bool, 
         outdir)
    plot(processed_df, 
         prefix, 
         "after", 
         plot_bool, 
         outdir)
    
    # Inform completion
    info(f"Output saved to {outdir}")
    info("Complete")
    print(f"Complete. Output saved to {outdir}")

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()

