import traitlets as tr
import numpy as np

def num_from_spacing(xmin, xmax, dx):
    return int(np.floor((xmax - xmin) / dx))
def spacing_from_num(xmin, xmax, nx):
    return (xmax - xmin) / nx

class Grid(tr.HasTraits):
    xmin = tr.Float()
    xmax = tr.Float()
    ymin = tr.Float()
    ymax = tr.Float()
    zmin = tr.Float()
    zmax = tr.Float()

    dx = tr.Float()
    dy = tr.Float()
    dz = tr.Float()

    nx = tr.Int()
    ny = tr.Int()
    nz = tr.Int()

    def __init__(self):
        self.xmin, self.xmax, self.dx = -1, 1, 2.5e-2
        self.ymin, self.ymax, self.dy = -1, 1, 2.5e-2
        self.zmin, self.zmax, self.dz = 0, 10, 1
        self.ntheta = 20
        self.nphi = 10

        self.num_from_spacing()
        self.assign_linspace()


    def num_from_spacing(self):
        self.nx = num_from_spacing(self.xmin, self.xmax, self.dx)
        self.ny = num_from_spacing(self.ymin, self.ymax, self.dy)
        self.nz = num_from_spacing(self.zmin, self.zmax, self.dz)

    def assign_linspace(self):
        self.x = np.arange(self.xmin, self.xmax, self.dx)
        self.y = np.arange(self.ymin, self.ymax, self.dy)
        self.z = np.arange(self.zmin, self.zmax, self.dz)


class Rope(tr.HasTraits):
    frond_lengths = tr.Any()
    frond_stds = tr.Any()
    water_speeds = tr.Any()
    water_angles = tr.Any()

    grid = tr.Any()

    def __init__(self, grid):
        self.grid = grid

        a = 2e-1
        b = 1e-1
        c = 5e-1
        d = np.pi / 4

        z = grid.z

        #self.frond_lengths = np.exp(-a*z) * np.sin(z) ** 2
        self.frond_lengths = .5 * z**2 * np.exp(1-z)
        self.frond_stds = b * np.ones_like(z)
        self.water_speeds = c * np.ones_like(z)
        #self.water_angles = 2*np.pi / grid.zmax * z
        self.water_angles = d * np.ones_like(z)

class Frond(tr.HasTraits):
    fs = tr.Float()
    fr = tr.Float()

    def __init__(self):
        self.fs = .5
        self.fr = 2

class Params(tr.HasTraits):
    quadrature_degree = tr.Int()

    def __init__(self):
        self.quadrature_degree = 5

class Kelp(tr.HasTraits):
    p_kelp = tr.Any()

    def __init__(self):
        pass
