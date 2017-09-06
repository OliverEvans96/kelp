import h5py

# HDF Read/Write Functions
def hdf_write_trait(hdf, obj, dset_name, dims, dtype):
    try:
        dset = hdf.create_dataset(dset_name, dims, dtype)
    except RuntimeError:
        dset = hdf.require_dataset(dset_name, dims, dtype)

    dset[...] = getattr(obj, dset_name)

def hdf_write_grid(filename, grid):
    grid.bounds = (grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax)
    grid.nums = (grid.nx, grid.ny, grid.nz, grid.ntheta, grid.nphi)
    #grid.spacings = (grid.dx, grid.dy, grid.dz, grid.dtheta, grid.dphi)

    with h5py.File(filename) as file:
        hdf_write_trait(file, grid, 'bounds', (6,), 'd')
        hdf_write_trait(file, grid, 'nums', (5,), 'i')
        #hdf_write_trait(file, grid, 'spacings', (5,), 'd')

def hdf_write_rope(filename, rope, grid):
    nz = grid.nz

    with h5py.File(filename) as file:
        hdf_write_trait(file, rope, 'frond_lengths', (nz,), 'd')
        hdf_write_trait(file, rope, 'frond_stds', (nz,), 'd')
        hdf_write_trait(file, rope, 'water_speeds', (nz,), 'd')
        hdf_write_trait(file, rope, 'water_angles', (nz,), 'd')

def hdf_write_frond(filename, frond):
    frond.frond_arr = (frond.fs, frond.fr)
    with h5py.File(filename) as file:
        hdf_write_trait(file, frond, 'frond_arr', (2,), 'd')

def hdf_write_params(filename, params):
    params.param_arr = (params.quadrature_degree,)
    with h5py.File(filename) as file:
        hdf_write_trait(file, params, 'param_arr', (1,), 'd')

def hdf_read_kelp(filename, kelp, grid):
    nums = (grid.nx, grid.ny, grid.nz)

    with h5py.File(filename) as file:
        # Transpose indices since data is written in fortran order (column-major)
        dset = file.require_dataset('p_kelp', nums[::-1], 'd')
        kelp.p_kelp = dset[...].T


