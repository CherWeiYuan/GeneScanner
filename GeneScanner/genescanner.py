'''
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
'''

import pandas as pd

input_file = "./test/test_input/input_basic_test.csv"
peak_gap = 1 # <= peak_gap will be processed
threshold = 0.2

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

def findPeakCluster(index, build_list, df, peak_gap, peak_cluster):
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
    build_list : TYPE
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
    peak_gap : TYPE Integer
        DESCRIPTION. User-supplied. A pair of peaks within peak_gap of
                     each other will be processed to give one peak.
    peak_cluster : TYPE List
        DESCRIPTION. A list of index of df that corresponds to the
                     rows (i.e. peaks) that are clustered (within peak_gap
                     of each other) and awaiting to be processed to give 
                     fewer peaks.

    Returns
    -------
    TYPE List
    A list of index corresponding to peaks in a peak cluster. 
    """    
    if index in peak_cluster:
        return []
    if df.loc[index + 1, "Size"] - df.loc[index, "Size"] > peak_gap:
        return list(set(build_list))
    else:
        build_list += [index, index + 1]
        return findPeakCluster(index + 1, build_list, df, peak_gap)

def findMountainRanges(df, sample_names):
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
            if index in mountain_ranges:
                pass
            else:
                mountain_ranges += [ [findPeakCluster(index, [], df, peak_gap)] ]
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
                area_list += df.loc[i, "Area"]
            max_area = max(area_list)
            for k in range(len(area_list)):
                # Keep index in peak_cluster if its peak
                # is above user-defined threshold
                if (area_list[k]/ max_area) >= threshold:
                    clean_peak_cluster += [peak_cluster[k]]
        clean_mountain_ranges += [ clean_peak_cluster ]
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
                remove += [mountains[1]]
            elif b > a:
                remove += [mountains[0]]
            # Keep both peaks if they are equal
            elif a == b:
                pass
    return remove

def RemoveArtefacts(df, remove):
    return df.drop(labels = remove, axis = 0, inplace = False)
        
def AddPercentage(df, sample_names):
    """
    Parameters
    ----------
    df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe of cleaned GeneScan datasheet.
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
        dfSample = df.loc[df.iloc[:,1] == i, :]
        
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

