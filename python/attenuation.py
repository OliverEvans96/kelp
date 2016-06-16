# attenuation.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Thu 02 Jun 2016 11:54:11 AM EDT
# Last Edited: Thu 16 Jun 2016 02:34:53 PM CEST

# Use data from light_data.py to calculate attenuation in water with and without kelp
# For each time step, use least squares to perform exponential regression to model
# intensity as a function of depth. The two parameters from the fit are k (kk),
# the exponential decay factor and I0, surface intensity. For each time step, plot
# a point in k-I0 space to show the results of the fit. Color the point with the 
# timestep. Also, plot the distributions of k and I0. Also, plot the r_squareds 
# as a function of time

from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import pickle
from datetime import datetime
from keyboard import keyboard
from sys import argv
from os import system

# Linear least squares
# N data points
# y = ax + b
# return [a,b]
def least_squares_fit(x,y):
    # Number of points
    N = len(x)

    # Remove non-numerical values
    good_vals = logical_and(isfinite(x),isfinite(y))
    x = x[good_vals]
    y = y[good_vals]

    # Set up and solve system
    A = array([[sum(x**2),sum(x)],[sum(x),N]])
    b = array([sum(x*y),sum(y)])

    # Note that linalg.lstsq is matrix equation solver from NumPy 
    return linalg.lstsq(A,b)[0]

# Name of this run of attenuation.py
# Useful for trying different parameters 
# and saving cases separately
try:
    run_name = argv[1]
except IndexError:
    run_name = 'default'

run_dir = run_name + '/'

# Make sure run_dir exists
system('mkdir -p ../results/attenuation/{}plots'.format(run_dir))

print("Run name: '{}'".format(run_name))

# Load data from light_data.py
with open('../data/Data HOBO/light_attenuation_data_datasets.pickle','rb') as pickle_file:
    dataset_array = pickle.load(pickle_file)

#Set font
font = {'family':'serif','size':10}
mpl.rc('font',**font)

# Depth keys
depth_range = arange(1,9)
depth_list = ['{}m'.format(zz) for zz in depth_range]
n_depths = len(depth_range)

# Number of datasets
n_datasets = len(dataset_array)

# Depth values to exclude for each dataset
# e.g. 0 => exclude 1m, 6 => exclude 7m
exclude_depths=[[],[],[6]]

# Names of datasets
dataset_labels = ['Control','1 Kelp','2 Kelp']
dataset_filenames = ['control','1_kelp','2_kelp']

# Only use data collected before 18:00
for arr in dataset_array:
    end_time = datetime(2016,5,18,18,0)
    arr = arr[(arr[:,0,0]<end_time),:,:]

# Number of timesteps for each dataset
n_steps_list = [dataset.shape[0] for dataset in dataset_array]

# Number of parameters
n_params = 3

# Initialize array to save calculated parameter values for each string
parameters = [zeros([n_steps,3]) for n_steps in n_steps_list]
good_parameters = [zeros([n_steps,3]) for n_steps in n_steps_list]

# Quantities whose distributions to plot
n_quantities = 4

# Axes limits for each quantity
dist_limits = [[0,1.0],[0,150000],[0,1],[0,50000]]

# Axes limits for calculated parameters for each dataset
param_limits = [[[0.1,0.4],[0,50000]],
                [[0.0,0.7],[0,50000]],
                [[-1.,1.0],[0,10]]]

## Number of bins to use in joint plot marginals
#n_param_bins = 10
#
## Bins to use for each parameter for each dataset
#param_bins = [[arange(*param_limits[ii][jj],n_param_bins) 
#    for jj in range(2)] for ii in range(n_datasets)]

# Colors to use for datasets
dataset_colors = ['b','g','r']

# Loop through datasets
for dataset_num,data in enumerate(dataset_array):

    # Loop through timesteps
    for step_num in range(n_steps_list[dataset_num]):
        # Light intensity
        II = data[step_num,:,1].astype(float)
        # Depth (starts at 1m, increments by 1m)
        zz = depth_range

        # Remove excluded depths
        II = delete(II,exclude_depths[dataset_num])
        zz = delete(zz,exclude_depths[dataset_num])

        # Solve for exponential fit using least squares
        # I = I_0 * exp(-kz)
        # => ln(I) = ln(I_0) - kz
        X = least_squares_fit(zz,log(II))
        kk = -X[0]
        I0 = exp(X[1])

        # Model values for ln(I)
        yy = log(I0) - kk*zz
        
        # Sum of squares error
        SSE = sum((log(II) - yy)**2)
        # Sum of squares total
        SST = sum((log(II) - mean(log(II)))**2)

        # Coefficient of determination
        # "Proportion of observed y variation that can be 
        # explained by the simple linear regression model"
        # - Devore: Prob & Stats for Engineering..., p. 485
        r_squared = 1 - SSE/SST

        # Save values
        parameters[dataset_num][step_num,:] = [kk,I0,r_squared]

    # Remove timesteps where fit fails (>=1 parameter == inf)
    #non_numeric_rows = where(logical_or(isinf(parameters[dataset_num]),isnan(parameters[dataset_num])))[0]
    #good_parameters[dataset_num] = delete(parameters[dataset_num],non_numeric_rows,axis=0)

    # Extract calculated variables
    kk = parameters[dataset_num][:,0]
    I0 = parameters[dataset_num][:,1]
    r_squared = parameters[dataset_num][:,2]

    # Intensity for all timesteps and all depths
    II_all = data[:,:,1].flatten()

    # Which quantities to plot
    plot_quantities = [kk,I0,r_squared,II_all]
    fig_titles=['k','$I_0$','$R^2$','Intensity']
    fig_filenames=['k','I0','r_squared','intensity']

    # Report statistical data
    print("{} Stats:".format(dataset_labels[dataset_num]))
    for ii,qq in enumerate(plot_quantities):
        print("{}:".format(fig_titles[ii]))
        print("mean={:.2g}".format(mean(qq)))
        print("std={:.2g}".format(std(qq)))
        print()

    # Plot k & I0 results
    figure(1,figsize=[8,6])
    plot(kk,I0,'o',label=dataset_labels[dataset_num],alpha=0.6)

    # Time plot of parameters
    time=data[:,0,0]
    # Plot distributions 
    for fig_num in range(n_quantities-1):
        # Time plots
        figure(fig_num*2+2,figsize=[7,3])
        plot(time,parameters[dataset_num][:,fig_num],color=dataset_colors[dataset_num],alpha=0.8,
            label=dataset_labels[dataset_num])
        title('Time plot: {}'.format(fig_titles[fig_num]))
        xlabel('time')
        ylabel(fig_titles[fig_num])
        legend()
        savefig('../results/attenuation/{}plots/time_{}.png'.format(run_dir,fig_filenames[fig_num]))
        savefig('../results/attenuation/{}plots/time_{}.eps'.format(run_dir,fig_filenames[fig_num]))

        # Distribution plots
        figure(fig_num*2+3,figsize=[7,3])
        qq = plot_quantities[fig_num]
        partial_data = qq[logical_and(dist_limits[fig_num][0]<qq,qq<dist_limits[fig_num][1])]
        sns.distplot(partial_data,
            label=dataset_labels[dataset_num],
            kde_kws={"shade": False})
        title('Dist plot: {}'.format(fig_titles[fig_num]))
        xlabel('value')
        ylabel('occurrence')
        gca().set_xlim(*dist_limits[fig_num])
        legend()
        savefig('../results/attenuation/{}plots/dist_{}.png'.format(run_dir,fig_filenames[fig_num]))
        savefig('../results/attenuation/{}plots/dist_{}.eps'.format(run_dir,fig_filenames[fig_num]))

    # Jointplots
    (sns.jointplot(kk,I0,space=0,kind='kde',stat_func=None,
        color=dataset_colors[dataset_num],
        joint_kws={
            'title':dataset_labels[dataset_num],
            'shade':False})
        .set_axis_labels('k','$I_0$'))
    title(dataset_labels[dataset_num])

    tight_layout()
    savefig('../results/attenuation/{}plots/joint_{}.png'.format(run_dir,dataset_filenames[dataset_num]))
    savefig('../results/attenuation/{}plots/joint_{}.eps'.format(run_dir,dataset_filenames[dataset_num]))

# Save calculated parameter data from fitting
with open('../results/attenuation/{}fit_parameters.pickle'.format(run_dir),'wb') as param_file:
    pickle.dump(parameters,param_file)

# Save other information about this run
run_info = {'exclude_depths':exclude_depths}
with open('../results/attenuation/{}run_info.pickle'.format(run_dir),'wb') as info_file:
    pickle.dump(run_info,info_file)

# Title, etc. for k-I0 plot
figure(1,figsize=[8,6])
title("k vs. $I_0$")
xlabel("k")
ylabel("$I_0$")
legend()
tight_layout()
savefig('../results/attenuation/{}plots/k_I0.png'.format(run_dir))
savefig('../results/attenuation/{}plots/k_I0.eps'.format(run_dir))

# Show plots
show(block=False)

print("attenuation.py done!")
