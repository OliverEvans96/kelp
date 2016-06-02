# light_data.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Wed 01 Jun 2016 10:11:00 AM EDT
# Last Edited: Thu 02 Jun 2016 11:12:25 AM EDT

from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
from datetime import datetime
from scipy import io
import os

# Function to format dates from file
def date_parser(date_str,time_str):
    return datetime.strptime("{} {}".format(date_str,time_str),"%m.%d.%y %I:%M:%S%p")

# Convert array of Python datetime object to Matlab datenum format
# Matlab starts from "0/0/0"
# Python starts from year 1
# 86400 seconds per day
def datenum(date_list):
    datenum_list = zeros_like(date_list,dtype=float)
    ii = 0
    for date in date_list:
        dt = date.to_datetime() - datetime(1,1,1)
        datenum_list[ii] = dt.days + 367 + dt.seconds/86400
        ii += 1
    return datenum_list

# Directory where data files are located
data_dir = "../data/Data HOBO/txt/"

# List of two dictionaries of data arrays: One for each string
str_data = [{},{}];

# Directory names for each string
str_names = ["Streng_1","Streng_2"]

# Number of strings
n_strings = len(str_names)

# Directory for each string
str_dir = [data_dir+name+"/" for name in str_names]

# Lists of files for each string
str_files = [[f for f in os.listdir(data_dir+name) if not f.startswith('.')] for name in str_names]

# Sort file lists
for file_list in str_files:
    file_list.sort()

# Set up plot canvas
ion()
fig = figure(1,figsize=[16,6])
logfig = figure(2,figsize=[16,6])
# Arrange plots on a grid of the following size
plot_grid = [1,n_strings]
# Interactively display plots


#Set font
font = {'family':'serif','size':10}
mpl.rc('font',**font)

# Loop through both strings
for str_num,file_list in enumerate(str_files):
    print("String #{}".format(str_num))

    # Set up subplot
    figure(1)
    ax = subplot(*plot_grid,str_num+1)
    title(str_names[str_num])

    figure(2)
    logax = subplot(*plot_grid,str_num+1)
    title(str_names[str_num])

    # Loop through files for this string
    for filename in file_list:
        print("Filename: {}".format(filename))
        # Depth of this file (remove .txt from filename)
        depth = filename[:-4]

        # Load data file to a Pandas DataFrame (df)
        # tab delimited
        # skip first two lines
        # don't use headers on second line (too messy)
        # explicitly specify column headers
        # Date information is stored in columns 1 and 2
        # Use above defined date parsing function
        # ignore bad lines causing mysterious error
        # comma is used as thousands separator in this data
        # Specify data type for each column
        df = pd.read_csv(
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

        # Strip last column (Only entry in this column is 'Logged' on last line + header string on 2nd line)
        df = df.iloc[:,:-1]

        # Isolate data from May 18th, 2016
        start_time = datetime(2016,5,18)
        end_time = datetime(2016,5,19)
        df = df.loc[logical_and(start_time<df['date_time'],df['date_time']<end_time)]

        # Plot data
        figure(1)
        df.plot(x='date_time',y='intensity',ax=ax,label=depth)
        xlabel('Time')
        ylabel('Intensity')

        figure(2)
        df.plot(x='date_time',y='intensity',ax=logax,label=depth,logy=True)
        xlabel('Time')
        ylabel('Intensity (log)')

        # Convert Python datetime to Matlab datenum
        df['date_time'] = datenum(df['date_time'])

        # Save DataFrame to dictionary as NumPy array
        key_name = str_names[str_num] + "_" + depth
        str_data[str_num][key_name] = array(df)

# Combine dictionaries from both strings
overall_dict = {**str_data[0],**str_data[1]}

# Save .mat files
mat_filename = '../data/Data HOBO/mat/light_attenuation_data.mat'
io.savemat(mat_filename,overall_dict)

# Final plot formatting
tight_layout()

# Save plots
figure(1)
savefig('../plots/light_data.eps')
savefig('../plots/light_data.png')

figure(2)
savefig('../plots/light_data_log.eps')
savefig('../plots/light_data_log.png')


