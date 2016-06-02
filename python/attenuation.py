# attenuation.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Thu 02 Jun 2016 11:54:11 AM EDT
# Last Edited: Thu 02 Jun 2016 03:37:29 PM EDT

# Use data from light_data.py to calculate attenuation in water with and without kelp
# For each time step, use least squares to perform exponential regression to model
# intensity as a function of depth. The two parameters from the fit are k (kk),
# the exponential decay factor and I0, surface intensity. For each time step, plot
# a point in k-I0 space to show the results of the fit. Color the point with the 
# timestep. Also, plot the distributions of k and I0. Also, plot the residuals 
# as a function of time

from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import pickle
from keyboard import keyboard

# Load data from light_data.py
with open('../data/Data HOBO/light_attenuation_data.pickle','rb') as pickle_file:
    str_array = pickle.load(pickle_file)

# Depth keys
depth_range = arange(1,9)
depth_list = ['{}m'.format(zz) for zz in depth_range]
n_depths = len(depth_range)

# Names of strings
str_names = ['Streng_1','Streng_2']

# Number of timesteps
n_steps = str_array[0].shape[0]

# Initialize array to save calculated parameter values for each string
parameters = [zeros([n_steps,3]),zeros([n_steps,3])]

# Figure counter
num_figs = 4

# Loop through strings
for str_num,data in enumerate(str_array):
    # Loop through timesteps
    for step_num in range(n_steps):
        
        # Initialize matrix system of equations
        A = zeros([2,2])
        b = zeros(2)

        # Light intensity
        II = data[step_num,:,2]
        # Depth (starts at 1m, increments by 1m)
        zz = depth_range + 1

        # Set values for Least squares system
        A[0,0] = sum(zz**2)
        A[0,1] = sum(zz)
        A[1,0] = sum(zz)
        A[1,1] = n_depths
        b[0] = sum(-zz*II)
        b[1] = sum(-II)

        # Solve system
        [kk,lnI0] = linalg.solve(A,b)
        I0 = exp(lnI0)

        # Calculate residual (sum of squared distances between data & model)
        res = sum((II + kk*zz - lnI0)**2)

        # Save values
        parameters[str_num][step_num,:] = [kk,I0,res]
    
    # Remove timesteps where fit fails (>=1 parameter == inf)
    inf_rows = where(parameters[str_num]==inf)[0]
    parameters[str_num] = delete(parameters[str_num],inf_rows,axis=0)

    # Extract calculated variables for this timestep
    kk = parameters[str_num][:,0]
    I0 = parameters[str_num][:,1]
    res = parameters[str_num][:,2]
    
    # Intensity for all timesteps and all depths
    II_all = data[:,:,2].flatten()

    # Which quantities to plot
    plot_quantities = [kk,I0,res,II_all]
    fig_titles=['k','$I_0$','residuals','Intensity']

    # Plot distributions 
    for fig_num in range(num_figs):
        figure(fig_num)
        sns.distplot(plot_quantities[fig_num],label=str_names[str_num])

# Figure titles and legends
for fig_num in range(num_figs):
    figure(fig_num)
    title(fig_titles[fig_num])
    legend()

# Show plots
show(block=False)

