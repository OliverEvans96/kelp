include("test_definitions.jl")
using KelpTest

xmin = 0
ymin = 0
zmin = 0
xmax = 1
ymax = 1
zmax = 1
lims = [xmin ymin zmin; xmax ymax zmax]

nx = 1
ny = 1
nz = 10
nθ = 1
nϕ = 2
nums = [nx ny nz nθ nϕ]

i = 1
j = 1
k = 6
l = 1
m = 1
indices = [i, j, k, l, m]

s, ds, ã, gₙ, rad_scatter = test_traverse(lims..., nums..., indices...)
