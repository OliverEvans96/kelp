#!/usr/local/miniconda3/bin/python

#normal_length.jl
#Oliver Evans
#4-11-16
#Kelp Research

println("Loading...")
using PyPlot
println("Loaded!")

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
function N(k,x)
    return (k-1)/2*(1-erf((x-mu)/(sqrt(2)*sigma)))
end

#Irradiance as a function of depth and horizontal position
function I(z,x)
    result=zeros(length(z),length(x))
    for k=1:n
        vals=(k-1)*dz.<=z.<=k*dz
        zk=z[vals]
        result[vals,:]=(1-alpha).^N(k,x)'*I0.*exp(-Kd*zk)
	end
    return result
end

#Depth array
z=0:1e-2:zmax
x=0:1e-2:xmax

#Contour Plot
contourf(x,-z,I(z,x),10)
colorbar()
title("Irradiance over depth & length")
xlabel("\$x\$")
ylabel("\$z\$",rotation=true)

#Plot kelp plants
for k=1:n
	plot([0,xmax],[-k*dz,-k*dz],"--k")
end

savefig("normal_length.eps")
#show()
