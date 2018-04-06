include("test_definitions.jl")
using KelpTest

xmin = 0
ymin = 0
zmin = 0
xmax = 10
ymax = 10
zmax = 10
lims = [xmin ymin zmin; xmax ymax zmax]

nx = 10
ny = 10
nz = 10
nθ = 10
nϕ = 10
nums = [nx ny nz nθ nϕ]

i = 10
j = 10
k = 6
l = 8
m = 9
indices = [i, j, k, l, m]

s, ds, ã, gₙ, rad_scatter = test_traverse(lims..., nums..., indices...)
