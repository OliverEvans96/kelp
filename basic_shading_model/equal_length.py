#!/usr/local/miniconda3/bin/python

#shading.py
#Oliver Evans
#4-11-16
#Kelp Research

from numpy import *
from matplotlib.pyplot import *

##############
## 1D model ##
##############

#Percent light absorbed
alpha=0.1
#Initial irradiance (just below surface)
I0=5
#Attentuation coefficient
Kd=.1
#Maximum z
zmax=10
#Number of kelp plants
n=1000
#Spacing between kelp plants
dz=zmax/n

#Irradiance as a function of depth
def I(z):
    y=zeros_like(z)
    for k in range(1,n+1):
        vals=logical_and((k-1)*dz<z,z<=k*dz)
        zk=z[vals]
        y[vals]=(1-alpha)**k*I0*exp(-Kd*zk)
    return y

#Depth array
z=arange(0,zmax,1e-2)

#Plot
plot(z,I(z))

show()
