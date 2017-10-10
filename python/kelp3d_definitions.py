import numpy as np
import ipyvolume as ipv
import matplotlib.pyplot as plt

from kelp3d_objs import *
from kelp_widgets import *

grid = Grid()
rope = Rope(grid)
frond = Frond()
params = Params()
iops = OpticalProperties(grid)
kelp = Kelp(grid, rope, frond, params)
bc = BoundaryCondition(grid)
light = Light(kelp, iops, bc, params)
