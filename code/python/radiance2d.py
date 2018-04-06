# radiance2d.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Fri 29 Jul 2016 02:04:11 PM EDT
# Last Edited: Wed 03 Aug 2016 09:19:27 AM EDT

# Solve radiance PDE in 2D.
# x = rcos(theta)
# y = rsin(Theta)

# Imports
from numpy import *
from matplotlib.pyplot import *
from keyboard import keyboard

# Find bin of x in [a,b] w/ n intervals
def find_bin(a,b,n,x):
    return floor((x-a)/(b-a)*n)

# Binary search

# Set up grid
dx = 1e-1
dy = 1e-1
dtheta = pi/50

xmin = 0
xmax = 10
ymin = 0
ymax = 10

thetamin = dtheta
thetamax = 2*pi

x = arange(xmin,xmax,dx)
y = arange(ymin,ymax,dx)
theta = arange(thetamin,thetamax,dx)

nx = len(x)
ny = len(y)
ntheta = len(theta)

L = zeros([nx,ny,ntheta])

# Simulation parameters
a = 0.5 # Absorptance (1/m)
b = 0.3 # Scatterance (1/m)
c = a + b # Attenuance (1/m)
theta_s = 3*pi/2 # Angle of sun (rad)
k_theta_s = int(digitize(theta_s,theta))
maxiter=10 # Max # of iterations

# Set up VSF
vsf_arr = loadtxt('../data/vsf/nuc_vsf.txt')
vsf_th = arange(0,pi,dtheta)
vsf = interp(vsf_th,*vsf_arr.T)

"""
# Plot vsf
plot(vsf_th,vsf,'ob-')
plot(*vsf_arr.T,'or-')
semilogy()
show()
"""

# Initial conditions
# Light only from direction of sun
print("digit(theta_s={:.2f}) = {:2d}".format(theta_s,int(digitize(theta_s,theta))))
L[:,-1,int(digitize(theta_s,theta))] = 10

# Solve PDE
for m in range(maxiter):
    resid = 0
    # Start from layer below top
    for j in range(ny-2,-1,-1):
        print("############")
        print("## j = {:2d} ##".format(j))
        print("############")
        print()
        for i in range(nx):
            for k in range(ntheta-2):
                dLdx = (L[(i+1)%(nx),j,k] - L[(i-1)%nx,j,k]) / (2*dx)
                dLdy = (L[i,j+1,k] - L[i,j-1,k]) / (2*dy)

                # Gain from scattering
                L_star = 0
                """
                for l in range(ntheta):
                    if l != k:
                        L_star += vsf[abs(k-l)%len(vsf)]*L[i,j,l]
                print("L_star =",L_star)
                """

                #L_new = (dLdx*cos(theta[k]) + dLdy*sin(theta[k]) - L_star) / c
                # PDE
                # CD2 for j > 1
                if j > 1:
                    L_new = (L_star - dLdx*cos(theta[k]) - dLdy*sin(theta[k])) / c
                # FD2 for j = 1
                else:
                    L_new = ( (-4*L[i,j+1,k] + L[i,j+2,k]) / (2*dy) 
                            + (L_star - dLdx*cos(theta[k])) / sin(theta[k]) ) \
                            / (-3/(2*dy) + c/sin(theta[k]))

                #L_new = (L[i,j-1,k] + dy/sin(theta[k])*(L_star - dLdx)) / (1 - c*dy/sin(theta[k]))

                # Calculate residual & update

                if isfinite(resid):
                    resid += abs(L_new-L[i,j,k])/(dx*dy*dtheta)
                else:
                    print("NF:",resid)

                # Troubleshooting
                if(k == k_theta_s) and (j == ny+3):
                    print("({},{},{}): {} -> {}".format(i,j,k,L[i,j,k],L_new))
                    print("dLdx =",dLdx)
                    print("dLdy =",dLdy)
                    #print("dLdy =",dLdy)
                    print()
                if(k == k_theta_s and j>1):
                    print("({},{}): {:9.2e} -> {:9.2e}"
                            .format(i,j,L[i,j,k],L_new))
                    print("|-------------------------------| theta = {:.2f}"
                            .format(theta[k]))
                    print("|           {:9.2e}           | dLdx = {:9.2e}"
                            .format(L[i,j+1,k],dLdx))
                    print("| {:9.2e} {:9.2e} {:9.2e} | dLdy = {:9.2e}"
                            .format(L[(i-1)%nx,j,k],L[i,j,k],L[(i+1)%nx,j,k],dLdy))
                    print("|           {:9.2e}           | cos(theta) = {:.2f}"
                            .format(L[i,j-1,k],cos(theta[k])))
                    print("|-------------------------------| sin(theta) = {:.2f}"
                            .format(sin(theta[k])))
                    print()
                L[i,j,k] = L_new
            #print("| {:.2e} ".format(L[i,j,k_theta_s]),end='')
#        print("|")
#        print("|-------------------------------------------------------------")

    print("################")
    print("m = {}: resid = {}".format(m,resid))
    print("################")
    break
