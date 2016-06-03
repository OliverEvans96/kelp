# attenuation.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Thu 02 Jun 2016 11:54:11 AM EDT
# Last Edited: Thu 02 Jun 2016 08:55:33 PM EDT

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
from datetime import datetime
from keyboard import keyboard

# Linear least squares
# N data points
# y = ax + b
# return [a,b]
def least_squares_fit(x,y):
    N = len(x)
    A = array([[sum(x**2),sum(x)],[sum(x),N]])
    b = array([sum(x*y),sum(y)])
    # Note that linalg.lstsq is matrix equation solver from NumPy 
    return linalg.lstsq(A,b)[0]

# Load data from light_data.py
with open('../data/Data HOBO/light_attenuation_data.pickle','rb') as pickle_file:
    str_array = pickle.load(pickle_file)

# Depth keys
depth_range = arange(1,9)
depth_list = ['{}m'.format(zz) for zz in depth_range]
n_depths = len(depth_range)

# Names of strings
str_names = ['Streng_1','Streng_2']
# Streng_1 has no kelp
# Streng_2 has kelp
str_labels = ['Control','Kelp']

# Only use data collected before 18:00
for arr in str_array:
    end_time = datetime(2016,5,18,18,0)
    arr = arr[(arr[:,0,0]<end_time),:,:]

# Number of timesteps
n_steps = str_array[0].shape[0]

# Initialize array to save calculated parameter values for each string
parameters = [zeros([n_steps,3]),zeros([n_steps,3])]
good_parameters = [zeros([n_steps,3]),zeros([n_steps,3])]

# Quantities whose distributions to plot
n_quantities = 4

# Axes limits for each quantity
dist_limits = [[-.75,.75],[0,150000],[0,1e10],[0,50000]]

# Color maps to use
colors = ['b','g']
cmaps=["Blues","Greens"]

#Set font
font = {'family':'serif','size':10}
mpl.rc('font',**font)

# Loop through strings
for str_num,data in enumerate(str_array):
    # Loop through timesteps
    for step_num in range(n_steps):
        
        # Initialize matrix system of equations
        A = zeros([2,2])
        b = zeros(2)

        # Light intensity
        II = data[step_num,:,2].astype(float)
        # Depth (starts at 1m, increments by 1m)
        zz = depth_range + 1

        # Solve for exponential fit using least squares
        X = least_squares_fit(zz,log(II))
        kk = -X[0]
        I0 = exp(X[1])

        # Calculate residual (sum of squared distances between data & model)
        res = sum((II + kk*zz - log(I0))**2)

        # Save values
        parameters[str_num][step_num,:] = [kk,I0,res]
    
    # Remove timesteps where fit fails (>=1 parameter == inf)
    inf_rows = where(logical_or(isinf(parameters[str_num]),isnan(parameters[str_num])))[0]
    good_parameters[str_num] = delete(parameters[str_num],inf_rows,axis=0)

    # Extract calculated variables
    kk = good_parameters[str_num][:,0]
    I0 = good_parameters[str_num][:,1]
    res = good_parameters[str_num][:,2]

    # Intensity for all timesteps and all depths
    II_all = data[:,:,2].flatten()

    # Which quantities to plot
    plot_quantities = [kk,I0,res,II_all]
    fig_titles=['k','$I_0$','residuals','Intensity']

    # Report statistical data
    print("{} Stats:".format(str_labels[str_num]))
    for ii,qq in enumerate(plot_quantities):
        print("{}:".format(fig_titles[ii]))
        print("mean={:.2g}".format(mean(qq)))
        print("std={:.2g}".format(std(qq)))
        print()

    # Normalized timestep indices
    ind = data[:,0,1]/data[-1,0,1]

    # Plot k & I0 results
    figure(1,figsize=[8,6])
    plot(kk,I0,'o',label=str_labels[str_num],alpha=0.6)

    # Time plot of parameters
    figure(2,figsize=[16,16])
    time=data[:,0,0]
    # Plot distributions 
    for fig_num in range(n_quantities-1):
        # Time plots
        subplot(3,2,fig_num*2+1)
        plot(time,parameters[str_num][:,0],color=colors[str_num],alpha=0.8,
            label=str_labels[str_num])
        title('Time plot: {}'.format(fig_titles[fig_num]))
        xlabel('time')
        ylabel(fig_titles[fig_num])
        legend()

        # Distribution plots
        subplot(3,2,fig_num*2+2)
        qq = plot_quantities[fig_num]
        partial_data = qq[logical_and(dist_limits[fig_num][0]<qq,qq<dist_limits[fig_num][1])]
        sns.distplot(partial_data,
            label=str_labels[str_num],
            kde_kws={"shade": False})
        title('Dist plot: {}'.format(fig_titles[fig_num]))
        xlabel('value')
        ylabel('occurrence')
        gca().set_xlim(*dist_limits[fig_num])
        legend()

    tight_layout()
    savefig('../plots/attenuation/parameters.png')
    savefig('../plots/attenuation/parameters.eps')

    # Jointplots
    jplot = sns.JointGrid(kk,I0,space=0)
    jplot = jplot.plot_joint(sns.kdeplot,cmap=cmaps[str_num],fill=True)
    xlabel('k')
    ylabel('$I_0$')
    title(str_labels[str_num])
    jplot = jplot.plot_marginals(sns.kdeplot,color=colors[str_num])
    gca().set_xlim(-0.75,0.75)
    gca().set_ylim(0,50000)
    tight_layout()
    savefig('../plots/attenuation/joint_{}.png'.format(str_labels[str_num]))
    savefig('../plots/attenuation/joint_{}.eps'.format(str_labels[str_num]))

# Title, etc. for k-I0 plot
figure(1,figsize=[8,6])
title("k vs. $I_0$")
xlabel("k")
ylabel("$I_0$")
legend()
tight_layout()
savefig('../plots/attenuation/k_I0.png')
savefig('../plots/attenuation/k_I0.eps')

# Show plots
show(block=False)

