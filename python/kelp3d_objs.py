import traitlets as tr
import numpy as np
from numpy.polynomial.legendre import leggauss

class SpaceDim(tr.HasTraits):
    minval = tr.Float()
    maxval = tr.Float()
    num = tr.Int()
    vals = tr.Any()

    def __init__(self, minval, maxval, num):
        super().__init__()
        self.minval = minval
        self.maxval = maxval
        self.num = num

        self.init_logic()
        self.update()

    # Deal only with num for simplicity

    def spacing_from_num(self):
        self.spacing = (self.maxval - self.minval) / self.num

    def assign_linspace(self):
        self.vals = np.arange(self.minval, self.maxval, self.spacing)

    def update(self, *args):
        self.spacing_from_num()
        self.assign_linspace()

    def init_logic(self):
        self.observe(self.update, names=['minval', 'maxval', 'num'])

class AngleDim(tr.HasTraits):
    num = tr.Int()
    vals = tr.Any()

    def __init__(self, minval, maxval, num):
        super().__init__()
        self.num = num
        self.assign_legendre()

    def assign_legendre(self):
        self.vals = leggauss(self.num)

    def init_logic(self):
        self.observe(self.assign_linspace, names=['num'])


class Grid(tr.HasTraits):
    x = tr.Any()
    y = tr.Any()
    z = tr.Any()
    theta = tr.Any()
    phi = tr.Any()

    def __init__(self):
        super().__init__()
        self.x = SpaceDim(-1, 1, 10)
        self.y = SpaceDim(-1, 1, 10)
        self.z = SpaceDim(-1, 1, 10)
        self.theta = AngleDim(0, 2*np.pi, 10)
        self.phi = AngleDim(0, np.pi, 10)

class Rope(tr.HasTraits):
    frond_lengths = tr.Any()
    frond_stds = tr.Any()
    water_speeds = tr.Any()
    water_angles = tr.Any()

    grid = tr.Any()

    def __init__(self, grid):
        super().__init__()
        self.grid = grid

        a = 2e-1
        b = 1e-1
        c = 5e-1
        d = np.pi / 4

        z = grid.z.vals

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
        super().__init__()
        self.fs = .5
        self.fr = 2

class Params(tr.HasTraits):
    quadrature_degree = tr.Int()

    def __init__(self):
        super().__init__()
        self.quadrature_degree = 5

class Kelp(tr.HasTraits):
    p_kelp = tr.Any()

    def __init__(self):
        super().__init__()
        pass

class Light(tr.HasTraits):
    radiance = tr.Any()
    irradiance = tr.Any()

    def __init__(self):
        super().__init__()
        pass

class OpticalProperties(tr.HasTraits):
    a_kelp = tr.Float()
    b_kelp = tr.Float()
    a_water = tr.Float()
    b_water = tr.Float()
    grid = tr.Any()
    vsf = tr.Any()

    def __init__(self, grid):
        super().__init__()
        vsf = None

class BoundaryCondition(tr.HasTraits):

    def __init__(self):
        super().__init__()
        bc = None
        pass
