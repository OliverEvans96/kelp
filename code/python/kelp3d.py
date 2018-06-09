# Imports
import numpy as np
import subprocess
from hdf_kelp import *


# Initialization
def write_all(grid, rope, frond, params):
    gridfile = '../hdf5/kelp3d/grid.hdf5'
    ropefile = '../hdf5/kelp3d/rope.hdf5'
    frondfile = '../hdf5/kelp3d/frond.hdf5'
    paramfile = '../hdf5/kelp3d/param.hdf5'

    hdf_write_grid(gridfile, grid)
    hdf_write_rope(ropefile, rope, grid)
    hdf_write_frond(frondfile, frond)
    hdf_write_params(paramfile, params)

def calculate_kelp(grid, rope, frond, params):
    result = subprocess.run(
        'cd .. && ./bin/pykelp3d',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout = result.stdout.decode()
    stderr = result.stderr.decode()
    if len(stdout) > 0:
        print("OUT:")
        print(stdout)
    if len(stderr) > 0:
        print("ERR:")
        print(stderr)

def read_all(kelp, grid):
    kelpfile = '../hdf5/kelp3d/kelp.hdf5'
    radfile = '../hdf5/kelp3d/rad.hdf5'
    irradfile = '../hdf5/kelp3d/irrad.hdf5'

    hdf_read_kelp(kelpfile, kelp, grid)
    hdf_read_rad(radfile, light, grid)
    hdf_read_irrad(irradfile, light, grid)
