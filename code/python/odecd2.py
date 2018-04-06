# odecd2.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Wed 03 Aug 2016 11:19:50 AM EDT
# Last Edited: Wed 03 Aug 2016 02:42:34 PM EDT

# Solve dy/dx = -cy
# using iterative method w/ CD2
# and BD2 for final point

from numpy import *
from matplotlib.pyplot import *

# Create grid
xmin = 0
xmax = 3
dx = 5e-1
nx = int((xmax-xmin)/dx)
x = linspace(xmin,xmax,nx)
y = zeros_like(x)

# Initial guess
y[:] = 1

# Boundary condition
y0 = 1
y[0] = y0

# Decay parameter
c = 1

# Grouping
r_tilde = 1/(2*dx*c)

# Iteration parameters
maxiter = 5
abs_tol = 1e-3

# Plot actual solution
plot(x,exp(-c*x),label='Exact')

"""
## ITERATION METHOD ##
# Loop through iterations
for k in range(maxiter):
    # Reset error
    err = 0

    # Loop through points
    for i in range(1,nx):
        # CD2 for 0<=i<nx-1
        if i < nx-1:
            ynew = -(y[i+1] - y[i-1])/(2*dx*c)
        # BD2 for i=nx-1
        else:
            ynew = 2*dx*(4*y[i-1] - y[i-2]) / (3*(1+c))

        # Accumulate error & update
        err += abs(y[i]-ynew)
        print("{:.2e} -> {:.2e}: err[{}] = {:.2e}".format(y[i],ynew,i,err))
        y[i] = ynew

    # Plot
    plot(x,y,label='k={}'.format(k))

    # Report error
    print("k={}: err={:.2e}".format(k,err/dx))

    # Stopping criteria
    if err/nx < abs_tol:
        print("Break!")
        break
"""

## MATRIX METHOD ##
A = zeros([nx,nx])
b = zeros([nx,1])

# BC
A[0,0] = 1
b[0] = y0

# PDE (CD2)
for i in range(1,nx-1):
    A[i,i] = 1
    A[i,i-1] = -r_tilde
    A[i,i+1] = r_tilde

# PDE (BD2)
A[nx-1,nx-3] = -r_tilde
A[nx-1,nx-2] = 4*r_tilde
A[nx-1,nx-1] = 1+3*r_tilde

# Solve
y_mat = linalg.lstsq(A,b)[0]
print(y_mat)

# Plot
plot(x,y_mat,label='Matrix')

# Show plot
legend()
show()
