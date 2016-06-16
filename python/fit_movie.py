# fit_movie.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Tue 14 Jun 2016 09:42:56 AM CEST
# Last Edited: Thu 16 Jun 2016 02:26:42 PM CEST

# Create one figure with three subplots on a 1x3 grid, one for each
# dataset, plotting intensity as a function of depth for a particular
# time. Make a movie for all times. Control should play the whole time,
# and other two should turn on and off at the appropriate times.

from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
from datetime import datetime,timedelta
from keyboard import keyboard
import pickle
import time
from os import system
from sys import argv

# Convert float to nice scientific notation string for latex parsing
def sci_not(xx,digits=2):
    fmt_str = ('{:.'+'{}'.format(digits)+'e}')
    base = fmt_str.format(xx).split('e')[0]
    expnt = str(int(floor(log10(xx))))
    return r'{}\cdot 10^{}'.format(base,expnt)

# Name of this run of fit_movie.py
# Useful for trying different parameters 
# and saving cases separately
try:
    run_name = argv[1]
except IndexError:
    run_name = 'default'

run_dir = run_name + '/'

# Make sure run_dir exists
system('mkdir -p ../results/attenuation/{}movies/img'.format(run_dir))

# Colors to use
dataset_colors = ['b','g','r']
depth_colors = ['b','g','r','brown','m','c','y','k']

# Set font
font = {'family':'serif','size':10}
mpl.rc('font',**font)

# Interactive plotting
#ion()

# Create figure
fig = figure(figsize=[12,8])

# Load experimental data from light_data.py
with open('../data/Data HOBO/light_attenuation_data_datasets.pickle','rb') as data_file:
    dataset_array = pickle.load(data_file)

# Load fit parameters from attenuation.py
with open('../results/attenuation/{}fit_parameters.pickle'.format(run_dir),'rb') as param_file:
    parameters = pickle.load(param_file)

# Load fit parameters from attenuation.py
with open('../results/attenuation/{}run_info.pickle'.format(run_dir),'rb') as info_file:
    run_info = pickle.load(info_file)

# Number of datasets
n_datasets = len(dataset_array)

# Names of datasets
dataset_labels = ['Control','1 Kelp','2 Kelp']
dataset_filenames = ['control','1_kelp','2_kelp']

# Depth Labels
depth_labels = ['{}m'.format(x) for x in range(1,9)]
n_depths = len(depth_labels)

# Excluded depths
exclude_depths = run_info['exclude_depths']

# Continuous plotting variable for depth
zz = linspace(0,9,801)

# Discrete depth variable
zd = arange(1,n_depths+1)

# Same, but with excluded depths for some datasets
zd_ex = []

# Which depths to use for each dataset
depths_to_use = array([[True for depth in range(n_depths)] for ii in range(n_datasets)])

# Generate step numbers relative to first measurement taken
offsets = []

# Number of timesteps for each dataset
n_steps_list = [dataset.shape[0] for dataset in dataset_array]

# Find first and last timesteps overall
first_ts = min([dataset[:,0,0].min() for dataset in dataset_array])
last_ts = max([dataset[:,0,0].max() for dataset in dataset_array])

# Number of steps spanning whole data range
total_n_steps = int((last_ts-first_ts).value/1e10)

# Array of datetimes spanning whole collection period
dt_array = array([first_ts+timedelta(seconds=10)*ii
    for ii in range(total_n_steps)])

# Only use data collected before 18:00
for dataset in dataset_array:
    # Generate step number offsets from control
    offsets.append(int((dataset[0,0,0]-first_ts).value/1e10))

# Time bars to progress through overview plots
time_bars = []

###################
## Set up figure ##
###################

# Set figure margins
subplots_adjust(
    bottom = 0.15,
    left = 0.10,
    top = 0.85,
    right = 0.90,
    hspace = 0.25,
    wspace = 0.10)

# Title figure with time
fig_title = figtext(0.5,0.95,'',
    horizontalalignment='center',
    fontsize='16')

# Axis limits for intensity plots
log_limits = (1e0,1e6)

# Arrays for data & fit Line2D objects
data_lines = []
fit_lines = []

# Markers for time bar
time_markers = []

# Parameter text
param_text = []

# Loop through datasets
for dataset_num,dataset in enumerate(dataset_array):

    ## One timestep ##
    subplot(2,3,dataset_num+1)

    # Find first and last timesteps overall
    first_dataset_ts = dataset[0,0,0]
    last_dataset_ts = dataset[-1,0,0]

    # Axis Title
    title(dataset_labels[dataset_num],fontweight='bold',y=1.1)

    # X-axis label
    xlabel('Light Intensity (Lux)')

    # Set axis limits
    xlim(*log_limits)
    ylim(-9,0)

    # Axis labels
    if(dataset_num == 0):
        yticks(-zd,zd)
        ylabel('Depth (m)')
    else:
        yticks(visible=False)

    # Display parameters
    param_text.append(text(0,1.01,'',
        transform=gca().transAxes,
        fontsize='12',
        ha='left',va='bottom'))

    # Exclude certain depths
    zd_ex.append(delete(zd,exclude_depths[dataset_num]))

    # Which depths to use for this dataset
    depths_to_use[dataset_num,exclude_depths[dataset_num]] = False

    # Create plots
    data_lines.append(semilogx(
        zeros_like(zd_ex[dataset_num]),-zd_ex[dataset_num],'o',
        label='data',color=dataset_colors[dataset_num])[0])
    fit_lines.append(semilogx(zeros_like(zz),-zz,'--',
        label='model',color=dataset_colors[dataset_num])[0])
    
    # Create legend
    legend(bbox_to_anchor=(1,1),loc='lower right',borderpad=0)

    ## Overview plots ##
    # Plot data
    subplot(2,3,dataset_num+4)
    for depth_num,depth_label in enumerate(depth_labels):
        semilogy(dataset[:,depth_num,0],dataset[:,depth_num,1],
            color=depth_colors[depth_num],label=depth_label)

    # Set axis range
    xlim(first_dataset_ts,last_dataset_ts)
    ylim(log_limits)

    # X-axis labels
    xlabel('Time')

    # Y-axis labels only on left side
    if(dataset_num == 0):
        ylabel('Light Intensity (Lux)')
    else:
        yticks(visible=False)

    # Legend only on right side
    if(dataset_num == 2):
        legend(bbox_to_anchor=(1,1),loc='upper left')

    # Rotate x axis labels, horizontally align to right side of text
    xticks(rotation=45,ha='right')

    # Find extreme values of intensity
    I_min = dataset[:,:,1].min()
    I_max = dataset[:,:,1].max()

    # Set up time bars & save line object to list
    time_bar = semilogy([first_dataset_ts,first_dataset_ts],log_limits,'k--')
    time_bars.append(time_bar[0])

    # Markers
    time_markers.append(
        [semilogy(dataset[0,depth_num,0],dataset[0,depth_num,1],'o',
            mfc=depth_colors[depth_num],mec='w',mew=1)[0]
        for depth_num in range(n_depths)])

# Turn off all plot lines to start
for dataset_num in range(n_datasets):
    data_lines[dataset_num].set_visible(False)
    fit_lines[dataset_num].set_visible(False)
    time_bars[dataset_num].set_visible(False)
    for depth_num in range(n_depths):
        time_markers[dataset_num][depth_num].set_visible(False)
    param_text[dataset_num].set_visible(False)

# Whether each dataset has become active yet
active_flag = [False,False,False]

# Whether this is the first timestep a dataset is active
first_flag = [False,False,False]

# Whether this is the last timestep a dataset is active
last_flag = [False,False,False]

###############
## Plot Loop ##
###############

# Loop through timesteps (of control dataset)
for step_num,step_dt in enumerate(dt_array):

    #step_num+=700

    # Don't keep running after all datasets have ended
    if(step_num>=total_n_steps):
        break

    # Print step number
    print("Timestep {:04d}".format(step_num))
    #print("active_flag: {}".format(active_flag))
    #print("first_flag: {}".format(first_flag))
    #print("last_flag: {}".format(last_flag))
    #print()

    # Update figure title
    fig_title.set_text('Light Intensity Data - '
        + str(step_dt).split()[1])

    # Loop through datasets
    for dataset_num,dataset in enumerate(dataset_array):

        # Relative step number
        rsn = step_num - offsets[dataset_num]
        
        # Check whether this dataset is active, adjust flags
        if(0 <= rsn and rsn < n_steps_list[dataset_num]):
            if active_flag[dataset_num]:
                first_flag[dataset_num] = False
            else:
                first_flag[dataset_num] = True
                active_flag[dataset_num] = True
        else:
            if active_flag[dataset_num]:
                active_flag[dataset_num] = False
                last_flag[dataset_num] = True

        # Check for last step
        if(rsn == n_steps_list[dataset_num]):
            last_flag[dataset_num] = True
        else:
            last_flag[dataset_num] = False


        # Plot if active
        if active_flag[dataset_num]:

            # Activate
            if first_flag[dataset_num]:
                time_bars[dataset_num].set_visible(True)
                data_lines[dataset_num].set_visible(True)
                fit_lines[dataset_num].set_visible(True)
                for depth_num in range(n_depths):
                    time_markers[dataset_num][depth_num].set_visible(True)
                param_text[dataset_num].set_visible(True)

            # Extract parameters
            kk = parameters[dataset_num][rsn,0]
            I0 = parameters[dataset_num][rsn,1]

            # Extract coefficient of determination
            r_squared = parameters[dataset_num][rsn,2]

            # Update data (I(z))
            data_lines[dataset_num].set_xdata(dataset[rsn,depths_to_use[dataset_num],1])

            # Update exponential model
            fit_lines[dataset_num].set_xdata(I0*exp(-kk*zz))

            # Update overview plot time bars
            tt = dataset[rsn,0,0]
            time_bars[dataset_num].set_xdata([tt,tt])

            # Update time markers
            for depth_num in range(n_depths):
                time_markers[dataset_num][depth_num].set_xdata(tt)
                time_markers[dataset_num][depth_num].set_ydata(dataset[rsn,depth_num,1])

            # Display parameters
            param_text[dataset_num].set_text(('$k={:.2f}$'
                + '\n' + r'$I_0={}$   $R^2={:.2f}$')
                .format(kk,sci_not(I0),r_squared))

        else:
            # Deactivate
            if last_flag[dataset_num]:
                time_bars[dataset_num].set_visible(False)
                data_lines[dataset_num].set_visible(False)
                fit_lines[dataset_num].set_visible(False)
                for depth_num in range(n_depths):
                    time_markers[dataset_num][depth_num].set_visible(False)
                param_text[dataset_num].set_visible(False)

    # Save images
    draw()

    savefig('../results/attenuation/{}movies/img/{:04d}.png'.format(run_dir,step_num))

