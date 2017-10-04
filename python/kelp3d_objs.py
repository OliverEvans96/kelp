import traitlets as tr
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simps
import ipyvolume as ipv

from fortran_wrappers.pykelp3d_wrap import f90wrap_py_gen_kelp as gen_kelp_f90

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
    minval = tr.Float()
    maxval = tr.Float()
    vals = tr.Any()

    def __init__(self, minval, maxval, num):
        super().__init__()
        self.minval = minval
        self.maxval = maxval
        self.num = num
        self.update()
        self.init_logic()

    def assign_legendre(self):
        self.vals = self.affine_transform(
            vals=leggauss(self.num)[0],
            oldmin=-1,
            oldmax=1,
            newmin=self.minval,
            newmax=self.maxval
        )

    def affine_transform(self, vals, oldmin, oldmax, newmin, newmax):
        return newmin + (newmin-newmax)/(oldmin-oldmax) * (vals-oldmin)

    def update(self, *args):
        self.assign_legendre()

    def init_logic(self):
        self.observe(self.update, names='num')

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
        self.z = SpaceDim(0, 1, 10)
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

    def __init__(self, grid, rope, frond, params):
        super().__init__()
        self.grid = grid
        self.rope = rope
        self.frond = frond
        self.params = params

        self.gen_kelp()

    def volume_plot(self):
        "Transform so that the new y is the old -z for IPyVolume."
        return ipv.quickvolshow(np.swapaxes(self.p_kelp[:,:,::-1], 1, 2)) 

    def gen_kelp(self):
        # Grid
        xmin = self.grid.x.minval
        xmax = self.grid.x.maxval
        nx = self.grid.x.num
        ymin = self.grid.y.minval
        ymax = self.grid.y.maxval
        ny = self.grid.y.num
        zmin = self.grid.z.minval
        zmax = self.grid.z.maxval
        nz = self.grid.z.num

        # Rope
        frond_lengths = self.rope.frond_lengths
        frond_stds = self.rope.frond_stds
        water_speeds = self.rope.water_speeds
        water_angles = self.rope.water_angles

        # Frond
        fs = self.frond.fs
        fr = self.frond.fr

        # Params
        quadrature_degree = self.params.quadrature_degree

        # Kelp
        self.p_kelp = np.asfortranarray(np.zeros([nx, ny, nz]))

        # Call fortran
        gen_kelp_f90(
            xmin, xmax, nx,
            ymin, ymax, ny,
            zmax, nz,
            frond_lengths,
            frond_stds,
            water_speeds,
            water_angles,
            fs, fr,
            self.p_kelp
        )

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
        self.grid = grid
        self.init_vals()

    def init_vals(self):
        self.a_kelp = 1
        self.b_kelp = 1
        self.a_water = 1
        self.b_water = 1

        z = self.grid.z.vals
        self.set_vsf(np.exp(-z))

    def set_vsf(self, vsf):
        self.vsf = self.normalize(vsf)

    def normalize(self, vsf):
        norm = simps(x=self.grid.phi.vals, y=vsf)
        try:
            return vsf / norm
        except ZeroDivisionError:
            return vsf

class BoundaryCondition(tr.HasTraits):

    def __init__(self):
        super().__init__()
        bc = None
        pass
