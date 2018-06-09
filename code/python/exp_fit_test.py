# exp_fit_test.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Thu 02 Jun 2016 03:57:45 PM EDT
# Last Edited: Thu 02 Jun 2016 05:13:33 PM EDT


from numpy import *
from matplotlib.pyplot import *
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import pickle
from keyboard import keyboard

# Linear least squares
# N data points
# y = ax + b
# return [a,b]
def least_squares_fit(x,y):
    N = len(x)
    A = array([[sum(x**2),sum(x)],[sum(x),N]])
    b = array([sum(x*y),sum(y)])
    print("x: {}".format(x))
    print("y: {}".format(y))
    print("A: {}".format(A))
    print("b: {}".format(b))
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

# Number of timesteps
n_steps = str_array[0].shape[0]

# Initialize array to save calculated parameter values for each string
parameters = [zeros([n_steps,3]),zeros([n_steps,3])]

# Quantities whose distributions to plot
n_quantities = 4

# Total number of plot figures
n_figs = 0

str_num = 0
data = str_array[0]
step_num=2000

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

# Extract calculated variables
kk = parameters[str_num][step_num,0]
I0 = parameters[str_num][step_num,1]
res = parameters[str_num][step_num,2]

# Intensity for all timesteps and all depths
II_all = data[step_num,:,2].astype(float)

# Continuous plotting variable
z = linspace(1,8,1001)

# Plot
figure(1)
scatter(depth_range,II_all)
plot(z,I0*exp(-kk*z),'r')
title('regular: {}'.format(data[step_num,0,0]))

# Log plot
figure(2)
scatter(depth_range,log(II_all))
plot(z,log(I0)-kk*z)
title('log: {}'.format(data[step_num,0,0]))

# Show plots
show(block=False)

