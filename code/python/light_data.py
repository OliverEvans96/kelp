# light_data.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Wed 01 Jun 2016 10:11:00 AM EDT
# Last Edited: Thu 16 Jun 2016 10:46:58 AM CEST

# Parse data from HOBO loggers

## IMPORTANT NOTE ##
# The data contains some zero values in the 1-kelp and 2-kelp datasets,
# presumably because the light intensity fell below the minimum intensity 
# threshold of the HOBO loggers used 
# (http://www.onsetcomp.com/products/data-loggers/ua-002-08)
# In the 1-kelp and 2-kelp, the lowest non-zero value is 10.8.
# Shane recommended using ~10% of the lowest value, so I will replace
# all zero-values with a value of 1.

##############
## Imports  ## 
##############

from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
from datetime import datetime
import pickle
from keyboard import keyboard
from scipy import io
import os

###############
## Functions ##
###############

# Function to format dates from file
def date_parser(date_str,time_str):
    return datetime.strptime("{} {}".format(date_str,time_str),"%m.%d.%y %I:%M:%S%p")

# Convert array of Python datetime object to Matlab datenum format
# Matlab starts from the date "0/0/0"
# 86400 seconds per day
def datenum(date_list):
    datenum_list = zeros_like(date_list,dtype=float)
    ii = 0
    for date in date_list:
        dt = date.to_datetime() - datetime(1,1,1)
        datenum_list[ii] = dt.days + 367 + dt.seconds/86400
        ii += 1
    return datenum_list

# Given a 24-hour time string in the format 'HH:MM',
# return a datetime object of that time on May 18th, 2016
def str_to_datetime(time):
    date = '5.18.16 '
    date_time = date+time
    return datetime.strptime(date_time,'%m.%d.%y %H:%M')

# Extract datasets from specified intervals for one string
# time should be a 1d array of datetime objects of size n_steps
# intensities should be a 1d array of floats of size n_steps
# interval_list should be a list, the length of which corresponds
# to the number of datasets contained in the corresponding string
# returns a Pandas DataFrame fo size n_steps x 2,
# The columns are time and intensity
def extract_datasets(time,intensities,interval_list):
    # Initialize list of datasets
    dataset_list = []

    # Number of datasets in this string
    n_datasets = len(interval_list)//2

    # Loop through datasets in the string
    for dataset_num in range(n_datasets):
        
        # Extract start and end times
        start_time = str_to_datetime(interval_list[2*dataset_num])
        end_time = str_to_datetime(interval_list[2*dataset_num+1])

        # Boolean array of which timesteps to consider
        time_range = logical_and(start_time<time,time<end_time)

        # Extract time and data in this interval
        dataset_time = time[time_range]
        dataset_intensities = intensities[time_range]

        # Save to array
        dataset_arr = array([dataset_time,dataset_intensities],dtype=object).T

        # Convert to DataFrame
        dataset_df = pd.DataFrame(dataset_arr,columns=['date_time','intensity'])
        dataset_list.append(dataset_df)

    # Return filled list
    return dataset_list

##########
## Main ##
##########

# Directory where data files are located
data_dir = "../data/Data HOBO/txt/"

# Dictionaries in which to save variables for Matlab
matlab_dict_strings = {}
matlab_dict_datasets = {}

# Directory names for each string
str_names = ["Streng_1","Streng_2"]

# Number of strings
n_strings = len(str_names)

# Depths at which light was measured
depth_list = []
depth_labels = []

# List of two dictionaries of data arrays: One for each string
str_dicts = [{} for ii in range(n_strings)];

# Directory for each string
str_dir = [data_dir+name+"/" for name in str_names]

# Lists of files for each string
str_files = [[f for f in os.listdir(data_dir+name) if not f.startswith('.')] for name in str_names]

# Sort file lists
for file_list in str_files:
    file_list.sort()

# Initialize list of datasets
all_datasets = []

# Time intervals to extract each dataset from
# Odd entries are the beginnings of intervals
# Even entries are the ends of intervals
# (Number of entries must be even)
# 24-hour clock on May 18th, 2016
dataset_intervals = [['13:00','16:30'],['14:15','15:30','16:10','16:45']]

# Number of datasets
n_datasets = sum(len(interval)//2 for interval in dataset_intervals)

# Dataset counter
running_dataset_count = 0

# Dataset names for display
# Streng_2 actually contains 2 datasets: 1 kelp and 2 kelp ropes
dataset_names = ["Control","1 Kelp","2 Kelp"]

# Names to be used for saving files
dataset_filenames = ["control","1_kelp","2_kelp"]

# List of dictionaries of data arrays: One for each dataset
dataset_dicts = [{} for ii in range(n_datasets)];

# Loop through strings
for str_num,file_list in enumerate(str_files):
    print("String #{}".format(str_num))

    # Loop through files for this string
    for filename in file_list:
        print("Filename: {}".format(filename))

        # Depth of this file (remove .txt from filename)
        depth = filename[:-4]

        # Save depth
        if(str_num == 0):
            depth_labels.append(depth)
            depth_list.append(float(depth[:-1]))

        # Load data file to a Pandas DataFrame (df_strings)
        # tab delimited
        # skip first two lines
        # don't use headers on second line (too messy)
        # explicitly specify column headers
        # Date information is stored in columns 1 and 2
        # Use above defined date parsing function
        # ignore bad lines causing mysterious error
        # comma is used as thousands separator in this data
        # Specify data type for each column
        df_strings = pd.read_csv(
            str_dir[str_num]+filename,
            delim_whitespace=True,
            skiprows=2,
            header=None,
            names=['index','date','time','intensity','eof'],
            parse_dates=[[1,2]],
            date_parser=date_parser,
            error_bad_lines=False,
            thousands=',',
            dtype={'index':int,'date':str,'time':str,'intensity':float,'eof':str})

        # During date parsing, date and time columns are deleted
        # and date_time is inserted as the first column

        # Strip second column (index)
        # Strip last column (Only entry in this column is the word 
        # 'Logged' on last line + header string on 2nd line)
        df_strings = df_strings.iloc[:,[0,2]]

        # Isolate data from May 18th, 2016
        start_time = datetime(2016,5,18)
        end_time = datetime(2016,5,19)
        df_strings = df_strings.loc[logical_and(start_time<df_strings['date_time'],df_strings['date_time']<end_time)]

        # Convert Python datetime to Matlab datenum
        matlab_df_strings = df_strings.copy()
        matlab_df_strings['date_time'] = datenum(df_strings['date_time'])

        # Save string DataFrame to dictionary as NumPy array for Matlab and Python
        matlab_name_strings = str_names[str_num] + "_" + depth
        matlab_dict_strings[matlab_name_strings] = array(matlab_df_strings)
        str_dicts[str_num][depth] = array(df_strings)

        # Extract datasets from strings
        str_dataset_list = extract_datasets(
            time=df_strings['date_time'],
            intensities=df_strings['intensity'],
            interval_list=dataset_intervals[str_num])

        # Append datasets to list
        # str_dataset_num is only counting datasets in this string
        for str_dataset_num,df_datasets in enumerate(str_dataset_list):
            # Current total dataset number
            total_dataset_num = running_dataset_count + str_dataset_num

            # Convert zero-values to ones. See important note at top of this file
            df_datasets.loc[df_datasets.iloc[:,1]==0,[False,True]] = 1

            # Convert Python datetime to Matlab datenum
            matlab_df_datasets = df_datasets.copy()
            matlab_df_datasets['date_time'] = datenum(df_datasets['date_time'])

            # Save dataset DataFrame to dictionary as NumPy array for Matlab and Python
            matlab_name_datasets = ( dataset_filenames[total_dataset_num] 
                + "_" + depth )
            matlab_dict_datasets[matlab_name_datasets] = array(matlab_df_datasets)
            dataset_dicts[total_dataset_num][depth] = df_datasets

    # Count number of datasets found in all strings counted thus far
    running_dataset_count += (str_dataset_num + 1)

# Number of depths (number of files)
n_depths = len(depth_list)

# Set up plot canvas
fig = figure(1,figsize=[16,6])
logfig = figure(2,figsize=[16,6])

#Set font
font = {'family':'serif','size':10}
mpl.rc('font',**font)

# Arrange plots on a grid of the following size
plot_grid = [1,n_datasets]

# Loop through datasets
print()
print("Before")
for dataset_num,dataset_dict in enumerate(dataset_dicts):
    print("Dataset #{}".format(dataset_num))

    # Set up subplot
    figure(1)
    ax = subplot(*plot_grid,dataset_num+1)
    title(dataset_names[dataset_num])

    figure(2)
    logax = subplot(*plot_grid,dataset_num+1)
    title(dataset_names[dataset_num])

    # Loop through depths
    for depth_num,depth in enumerate(depth_labels):
        print("Depth: {}".format(depth))
        # Plot data
        figure(1)
        dataset_dict[depth].plot(x='date_time',y='intensity',ax=ax,label=depth)
        xlabel('Time')
        ylabel('Intensity')

        figure(2)
        dataset_dict[depth].plot(x='date_time',y='intensity',ax=logax,label=depth,logy=True)
        xlabel('Time')
        ylabel('Intensity (log)')

print("After")

# Number of timesteps for each array
n_steps_list_strings = [str_dicts[ii]['1m'].shape[0] 
    for ii in range(len(str_dicts))]
n_steps_list_datasets = [dataset_dicts[ii]['1m'].shape[0] 
    for ii in range(len(dataset_dicts))]

# For each string, convert dictionary of 2d arrays to a 3d array
str_list = []
for str_num,n_steps in enumerate(n_steps_list_strings):
    arr = zeros([n_steps,n_depths,2],dtype=object)
    for depth_num,depth_label in enumerate(depth_labels):
        arr[:,depth_num,:] = str_dicts[str_num][depth_label]
    str_list.append(arr)

# For each dataset, convert dictionary of 2d arrays to a 3d array
dataset_list = []
for dataset_num,n_steps in enumerate(n_steps_list_datasets):
    arr = zeros([n_steps,n_depths,2],dtype=object)
    for depth_num,depth_label in enumerate(depth_labels):
        arr[:,depth_num,:] = dataset_dicts[dataset_num][depth_label]
    dataset_list.append(arr)

# Save .mat files (Matlab)
mat_filename = '../data/Data HOBO/light_attenuation_data_strings.mat'
io.savemat(mat_filename,matlab_dict_strings)
mat_filename = '../data/Data HOBO/light_attenuation_data_datasets.mat'
io.savemat(mat_filename,matlab_dict_datasets)

# Save .pickle file (Python)
pickle_filename = '../data/Data HOBO/light_attenuation_data_strings.pickle'
with open(pickle_filename,'wb') as pickle_file:
    pickle.dump(str_list,pickle_file)
pickle_filename = '../data/Data HOBO/light_attenuation_data_datasets.pickle'
with open(pickle_filename,'wb') as pickle_file:
    pickle.dump(dataset_list,pickle_file)

# Save plots
figure(1)
tight_layout()
savefig('../plots/light_data.eps')
savefig('../plots/light_data.png')

figure(2)
tight_layout()
savefig('../plots/light_data_log.eps')
savefig('../plots/light_data_log.png')

print("light_data.py done!")

