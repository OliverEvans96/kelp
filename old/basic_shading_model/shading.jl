#!/usr/local/bin/julia

#shading.jl
#Oliver Evans
#4-11-16
#Kelp Research

using PyPlot

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
n=40
#Spacing between kelp plants
dz=zmax/n

#Irradiance as a function of depth
function I(z)
	y=zeros(z)
	for k=1:n
		vals=(k-1)*dz.<z.<=k*dz
		zk=z[vals]
		y[vals]=(1-alpha)^k*I0*exp(-Kd*zk)
	end
	return y
end

#Depth array
z=0:1e-2:zmax

#Plot
plot(z,I(z))

show()
