import traitlets as tr
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simps
import ipyvolume as ipv

from fortran_wrappers.pykelp3d_wrap import f90wrap_py_gen_kelp as gen_kelp_f90
from fortran_wrappers.pyrte3d_wrap import f90wrap_py_calculate_light_field as calculate_light_field_f90

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

    #def angular_integral(self):
    #   pass

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
        self.x = SpaceDim(-1, 1, 6)
        self.y = SpaceDim(-1, 1, 6)
        self.z = SpaceDim(0, 1, 6)
        self.theta = AngleDim(0, 2*np.pi, 6)
        self.phi = AngleDim(0, np.pi, 6)

class Rope(tr.HasTraits):
    frond_lengths = tr.Any()
    frond_stds = tr.Any()
    water_speeds = tr.Any()
    water_angles = tr.Any()

    grid = tr.Any()

    def __init__(self, grid):
        super().__init__()
        self.grid = grid

        a = 0 * 2e-1
        b = 0 * 1e-1
        c = 0 * 5e-1
        d = 0 * np.pi / 4

        z = grid.z.vals

        self.frond_lengths = a * np.exp(-a*z) * np.sin(z) ** 2
        #self.frond_lengths = .5 * z**2 * np.exp(1-z)
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
    maxiter_inner = tr.Int()
    maxiter_outer = tr.Int()
    tol_abs = tr.Float()
    tol_rel = tr.Float()

    def __init__(self):
        super().__init__()
        self.init_vals()

    def init_vals(self):
        self.quadrature_degree = 5
        self.maxiter_inner = 10
        self.maxiter_outer = 10
        self.tol_abs = 1e-3
        self.tol_rel = 1e-3

class Kelp(tr.HasTraits):
    p_kelp = tr.Any()

    def __init__(self, grid, rope, frond, params):
        super().__init__()
        self.grid = grid
        self.rope = rope
        self.frond = frond
        self.params = params

    def volume_plot(self):
        "Transform so that the new y is the old -z for IPyVolume."
        return ipv.quickvolshow(np.swapaxes(self.p_kelp[:,:,::-1], 1, 2)) 

    def gen_kelp(self, *args):
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

    def __init__(self, kelp, iops, bc, params):
        super().__init__()
        self.kelp = kelp
        self.grid = kelp.grid
        self.iops = iops
        self.bc = bc
        self.params = params

        self.init_vals()

    def init_vals(self):
        nx = self.grid.x.num
        ny = self.grid.y.num
        nz = self.grid.z.num
        ntheta = self.grid.theta.num
        nphi = self.grid.phi.num

        self.irradiance = np.asfortranarray(np.zeros([nx, ny, nz]))
        self.radiance = np.asfortranarray(np.zeros([nx, ny, nz, ntheta, nphi]))

    def calculate_light_field(self, *args):
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
        ntheta = self.grid.theta.num
        nphi = self.grid.phi.num

        # IOPs
        num_vsf = self.iops.num_vsf
        vsf_angles = self.iops.vsf_angles
        vsf_vals = self.iops.vsf_vals
        p_kelp = self.kelp.p_kelp
        a_water = self.iops.a_water
        a_kelp = self.iops.a_kelp
        b_water = self.iops.b_water
        b_kelp = self.iops.b_kelp

        # Boundary Condition
        theta_s = self.bc.theta_s
        phi_s = self.bc.phi_s
        max_rad = self.bc.max_rad
        decay = self.bc.decay

        # Params
        tol_abs = self.params.tol_abs
        tol_rel = self.params.tol_rel
        maxiter_inner = self.params.maxiter_inner
        maxiter_outer = self.params.maxiter_outer

        # Kelp
        p_kelp = self.kelp.p_kelp

        # Call fortran
        calculate_light_field_f90(
            xmin, xmax, nx,
            ymin, ymax, ny,
            zmax, nz,
            ntheta, nphi,
            a_water, a_kelp, b_water, b_kelp,
            num_vsf, vsf_angles, vsf_vals,
            theta_s, phi_s, max_rad, decay,
            tol_abs, tol_rel, maxiter_inner, maxiter_outer,
            p_kelp, self.radiance, self.irradiance
        )

    def volume_plot(self):
        "Transform so that the new y is the old -z for IPyVolume."
        return ipv.quickvolshow(np.swapaxes(self.irradiance[:,:,::-1], 1, 2)) 


class OpticalProperties(tr.HasTraits):
    a_kelp = tr.Float()
    b_kelp = tr.Float()
    a_water = tr.Float()
    b_water = tr.Float()
    grid = tr.Any()
    num_vsf = tr.Int()
    vsf_vals = tr.Any()
    vsf_angles = tr.Any()

    def __init__(self, grid):
        super().__init__()
        self.grid = grid
        self.init_vals()

    def init_vals(self):
        self.a_kelp = 2
        self.b_kelp = 1
        self.a_water = 1
        self.b_water = 1

        self.num_vsf = self.grid.phi.num
        self.vsf_angles = self.grid.phi.vals
        self.set_vsf(np.exp(-self.vsf_angles))

    def init_logic(self):
        tr.link(
            (self, 'vsf_angles'),
            (self.grid.theta, 'vals')
        )

    def set_vsf(self, vsf_vals):
        self.vsf_vals = self.normalize(self.vsf_angles, vsf_vals)

    def normalize(self, vsf_angles, vsf_vals):
        norm = np.abs(simps(x=vsf_angles, y=vsf_vals))
        try:
            return vsf_vals / norm
        except ZeroDivisionError:
            return vsf_vals

class BoundaryCondition(tr.HasTraits):
    "2D Angular Gaussian of light centered around sun position."

    # Direction of light ray from sun
    theta_s = tr.Float()
    phi_s = tr.Float()

    # Maximum intensity
    max_rad = tr.Float()

    # Exponential decay in other directions
    decay = tr.Float()

    grid = tr.Any()

    def __init__(self, grid):
        self.grid = grid
        super().__init__()
        self.init_values()

    def calc_rad(self, theta, phi):
        delta = angle_diff(theta, phi, self.theta_s, self.phi_s)
        return self.max_rad * np.exp(-decay * delta)

    def angle_diff(self, theta1, phi1, theta2, phi2):
        "Angle between vectors pointing in two directions"
        return np.arccos(
            np.sin(theta1)*np.sin(theta2)
            * (np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2))
            + np.cos(theta1)*np.cos(theta2)
        )

    def create_surf_rad_grid(self):
        theta_grid, phi_grid = np.meshgrid(
            self.grid.theta.vals,
            # Only take first half of phi vals (downwelling light)
            self.grid.phi.vals[:self.grid.phi.num/2],
            indexing='ij'
        )
        return calc_rad(theta_grid, phi_grid)

    def init_values(self):
        self.theta_s = 0
        self.phi_s = 0
        self.max_rad = 1
        self.decay = 1


