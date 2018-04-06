import numpy as np
import ipyvolume as ipv
import matplotlib.pyplot as plt
import IPython
import warnings

from kelp3d_objs import *
from kelp_widgets import *

warnings.filterwarnings('ignore', category=RuntimeWarning)

def init_objects():
    global grid, rope, frond, params, iops, kelp, bc, light
    grid = Grid()
    rope = Rope(grid)
    frond = Frond()
    params = Params()
    iops = OpticalProperties(grid)
    kelp = Kelp(grid, rope, frond, params)
    bc = BoundaryCondition(grid)
    light = Light(kelp, iops, bc, params)

def init_widgets():
    global fw, rw, iw, bcw, gw, pw, vpw, radw
    fw = FrondWidget(frond)
    rw = RopeWidget(rope)
    iw = IOPWidget(iops)
    bcw = BCWidget(bc)
    gw = GridWidget(grid)
    pw = ParamWidget(params)
    vpw = VolumePlotWidget(kelp, light)
    radw = RadianceWidget(light)

def stats(arr):
    print("min: {}".format(np.nanmin(arr)))
    print("max: {}".format(np.nanmax(arr)))
    print("mean: {}".format(np.nanmean(arr)))
    print("std: {}".format(np.nanstd(arr)))

init_objects()
init_widgets()

if __name__ == '__main__':
    kelp.gen_kelp()
    IPython.embed()
    light.calculate_light_field()
