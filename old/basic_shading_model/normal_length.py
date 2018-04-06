#!/usr/local/miniconda3/bin/python

#shading.py
#Oliver Evans
#4-11-16
#Kelp Research

from numpy import *
from matplotlib.pyplot import *
from scipy.special import erf

##############
## 2D model ##
##############

#Assume plant length is normally distributed with mean mu, std sigma
mu=6
sigma=3

#Box limits
zmax=10
xmax=10

#Percent light absorbed
alpha=0.1
#Initial irradiance (just below surface)
I0=5
#Attentuation coefficient
Kd=.1
#Number of kelp plants
n=10
#Spacing between kelp plants
dz=zmax/n

#Expected Number of shading plants
def N(k,x):
    return k/2*(1-erf((x-mu)/(sqrt(2)*sigma)))

#Irradiance as a function of depth and horizontal position
def I(z,x):
    y=zeros_like(z)
    for k in range(1,n+1):
        vals=logical_and((k-1)*dz<z,z<=k*dz)
        print(vals)
        zk=z[vals]
        xk=x[vals]
        y[vals]=(1-alpha)**N(k,x)*I0*exp(-Kd*zk)
    return y

#Depth array
Z=arange(0,zmax,1e-2)
X=arange(0,xmax,1e-2)

z,x=meshgrid(Z,X)


#Plot
contour(z,I(z,x))

show()
