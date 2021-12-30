'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) CHER_WEI_YUAN, 29DEC2021 
License     : MIT_LICENSE 
Maintainer  : E0031403@U.NUS.EDU 
Portability : POSIX
This program reads the output of GeneScan in csv format and 
calculates, for each sample, the percentage of the total area 
that each peak covers.

Features:
    For peaks within PEAK_GAP size of each other, the peak with 
    a smaller area will be ignored if mode is select_larger_area,
    or the two peaks will be combined if mode is combine_areas.
    A filter can be applied to remove peaks that cover below 
    a certain percentage of total area.
'''

import pandas as pd

input_file = "./test/test_input/input_basic_test.csv"
peak_gap = 1 # <= peak_gap will be processed
method = "pick_larger_peak" # other methods: "keep_separated_peaks" / "combine_peaks"/ "hybrid"

# Load dataframe
df = pd.read_csv(input_file)

# Sort dataframe by sample name and size
df = df.sort_values(by=['Sample File Name', 'Size']).\
    reset_index(drop = True, inplace = False)

# Clean sample names by removing leading and trailing white space
df["Sample File Name"] = df["Sample File Name"].str.rstrip().str.lstrip()
        
# Get sample names
sample_names = list(set(df.loc[:,"Sample File Name"]))

## Find mountain ranges: cluster of peaks within peak_gap along x-axis (size)
## Store index of peaks
def advance(index, build_list, df, peak_gap, mountain_ranges):       
    if index in mountain_ranges:
        return []
    if df.loc[index + 1, "Size"] - df.loc[index, "Size"] > peak_gap:
        return list(set(build_list))
    else:
        build_list += [index, index + 1]
        return advance(index + 1, build_list, df, peak_gap)

mountain_ranges = []
for i in sample_names:
    dfSample = df.loc[df.iloc[:,1] == i, :]
    for index, row in dfSample.iterrows():
        if index in mountain_ranges:
            pass
        else:
            mountain_ranges += [ [advance(index, [], df, peak_gap)] ]
        
# Keep largest peak and remove smaller peak(s)
for mountains in mountain_range:
    ## continue


# Loop through samples
data_list = []
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

