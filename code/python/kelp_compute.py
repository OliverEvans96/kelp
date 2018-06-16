# stdlib
import sqlite3
from datetime import datetime
import subprocess
import multiprocessing
import itertools as it
import os

# 3rd party
import netCDF4 as nc
import tempfile
import ipyparallel as ipp

###############
# Misc. utils #
###############

def get_random_unused_filename(dir='.', prefix='', suffix=''):
    with tempfile.NamedTemporaryFile(dir=dir, prefix=prefix, suffix=suffix) as fh:
        filename = fh.name

    return filename

def get_git_commit_hash():
    """Get hash of current git branch."""
    return subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()
    return filename


#######################
# Kelp-specific funcs #
#######################

def create_table(table_name, prefix='.'):
    """
    Create sqlite db file, create table, and return connection.
    Assume table_name is unique (e.g. includes date.)
    """

    create_table_template = '''
    CREATE TABLE {table_name} (

        /* Inputs */
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        absorptance_kelp REAL,
        a_water REAL,
        b REAL,
        ns REAL,
        nz REAL,
        na REAL,
        num_dens REAL,
        kelp_dist CHAR(32),
        fs REAL,
        fr REAL,
        ft REAL,
        max_length REAL,
        length_std REAL,
        zmax REAL,
        rope_spacing REAL,
        I0 REAL,
        phi_s REAL,
        theta_s REAL,
        decay REAL,
        num_cores INTEGER,
        num_scatters INTEGER,
        fd_flag INTEGER,
        lis_opts CHAR(256),
        date CHAR(64),
        git_commit CHAR(40),
        compute_time REAL,
        lis_iter INTEGER,
        lis_time REAL,
        lis_resid REAL,
        /* Array data (NETCDF) */
        data_path CHAR(256)
        )
    '''

    base_dir = os.path.join(prefix, table_name)
    data_dir = os.path.join(base_dir, 'data')
    db_path = os.path.join(base_dir, table_name+'.db')
    # Will fail if `table_name` is not unique in `prefix`
    os.mkdir(base_dir)
    os.mkdir(data_dir)

    conn = sqlite3.connect(db_path)
    conn.execute(create_table_template.format(table_name=table_name))
    conn.close()

    return base_dir, db_path, table_name

def get_table_names(conn):
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    return cursor.fetchall()

def insert_run(db_path, table_name=None, **params):
    insert_template = '''
    INSERT INTO {table_name} VALUES (
        NULL, /* id (autoincrement) */
        {absorptance_kelp},
        {a_water},
        {b},
        {ns},
        {nz},
        {na},
        {num_dens},
        '{kelp_dist}',
        {fs},
        {fr},
        {ft},
        {max_length},
        {length_std},
        {zmax},
        {rope_spacing},
        {I0},
        {phi_s},
        {theta_s},
        {decay},
        {num_cores},
        {num_scatters},
        {fd_flag},
        '{lis_opts}',
        '{date}',
        '{git_commit}',
        {compute_time},
        {lis_iter},
        {lis_time},
        {lis_resid},
        '{data_path}'
        )
    '''

    # If table_name is not provided,
    # just use the first one we find
    if not table_name:
        table_name = get_table_names(conn)[0]

    insert_cmd = insert_template.format(table_name=table_name, **params)

    # Execute SQL command
    conn = sqlite3.connect(db_path)
    cursor = conn.execute(insert_cmd)
    # Save changes
    conn.commit()
    print("Executed SQL:")
    print(insert_cmd)
    conn.close()

def create_nc(db_path, **results):
    """Create netCDF file to store results"""
    # Use directory containing DB
    base_dir = os.path.dirname(db_path)

    # Generate random name in correct directory
    data_path = get_random_unused_filename(
        dir=os.path.join(base_dir, 'data'),
        suffix='.nc'
    )

    # Create file
    rootgrp = nc.Dataset(data_path, 'w', format='NETCDF4')

    # Get dimension sizes
    nx = results['nx']
    ny = results['ny']
    nz = results['nz']
    ntheta = results['ntheta']
    nphi = results['nphi']
    nomega = results['nomega']

    # Create Dimensions
    x = rootgrp.createDimension('x', nx)
    y = rootgrp.createDimension('y', ny)
    z = rootgrp.createDimension('z', nz)
    theta = rootgrp.createDimension('theta', ntheta)
    phi = rootgrp.createDimension('phi', nphi)
    omega = rootgrp.createDimension('omega', nomega)

    # Create dimension variables
    # xvals = rootgrp.createVariable('x','f8', ('x',))
    # yvals = rootgrp.createVariable('y','f8', ('y',))
    # zvals = rootgrp.createVariable('z','f8', ('z',))
    # thetavals = rootgrp.createVariable('theta','f8', ('theta',))
    # phivals = rootgrp.createVariable('phi','f8', ('phi',))

    # Assign dimension variables
    # xvals[:] = results.pop('x')
    # yvals[:] = results.pop('y')
    # zvals[:] = results.pop('z')
    # thetavals[:] = results.pop('theta')
    # phivals[:] = results.pop('phi')

    # Create results variables
    rad = rootgrp.createVariable('rad', 'f8', ('x', 'y', 'z', 'omega'))
    irrad = rootgrp.createVariable('irrad', 'f8', ('x', 'y', 'z'))
    avg_irrad = rootgrp.createVariable('avg_irrad', 'f8', ('z',))
    perc_irrad = rootgrp.createVariable('perc_irrad', 'f8', ('z',))
    p_kelp = rootgrp.createVariable('p_kelp', 'f8', ('x', 'y', 'z'))

    # Assign results variables
    rad[:] = results.pop('rad')
    irrad[:] = results.pop('irrad')
    avg_irrad[:] = results.pop('avg_irrad')
    perc_irrad[:] = results.pop('perc_irrad')
    p_kelp[:] = results.pop('p_kelp')

    # Assume all others are scalar metadata
    for var_name, val in results.items():
        print("{}: {} ({})".format(var_name, val, type(val)))
        var = rootgrp.createVariable(var_name, type(val))
        var[...] = val

    # Close file
    rootgrp.close()

    return data_path

def grid_centers(xmin, xmax, nx):
    """Evenly spaced grid centers."""
    import numpy as np
    dx = (xmax - xmin) / nx
    x = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
    return x

def get_kelp_dist(kelp_dist, max_length, zmin, zmax, nz):
    import numpy as np
    # TODO: Scale by zmax
    # Kelp distribution profiles
    # Normalize by max val

    z = grid_centers(zmin, zmax, nz)

    maxval = 12*np.exp(-2) + 0.5
    frond_length_funcs = {
        'top-heavy': (3.0 * z**2 * np.exp(-z) + 0.5) / maxval,
        'bottom-heavy': (3.0 * (zmax-z)**2 * np.exp((z-zmax)) + 0.5) / maxval,
        'uniform': 0*z + 1.0,
        'none': 0*z
    }
    frond_lengths = max_length * frond_length_funcs[kelp_dist]
    frond_stds = 0.0 * np.ones_like(z)
    water_speeds = 0.5 * np.ones_like(z)
    water_angles = 2*np.pi / zmax * (z-zmin)

    return frond_lengths, frond_stds, water_speeds, water_angles



###############

# TODO: This is not right anymore, but it doesn't matter much
def nokelp_calculate(a_water, b, ns, na, const, num_threads=1):
    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    # Extract constants
    (rope_spacing, zmin, zmax, nz, I0, phi_s, theta_s, decay, xmin, xmax, ymin, ymax, absortpance_kelp,
         num_scatters, gmres_flag, gmres_flag, lis_options) = const
    a_kelp = a_water

    num_vsf = na
    vsf_angles = np.linspace(0,np.pi, na)
    vsf_vals = 0*vsf_angles + 1/(4*np.pi)
    ns = int(ns)

    nomega = int(na*(na-2)+2)
    p_kelp = np.asfortranarray(np.zeros([ns,ns,nz]))
    radiance = np.asfortranarray(np.zeros([ns, ns, nz, nomega]))
    irradiance = np.asfortranarray(np.zeros([ns, ns, nz]))

    # Start timer
    tic = time.time()

    # Calculate light field
    f90.calculate_light_field(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        na, na,
        a_water, a_kelp, b,
        vsf_angles, vsf_vals,
        theta_s, phi_s, I0, decay,
        p_kelp, radiance, irradiance,
        num_scatters, gmres_flag, lic_options,
    )

    # End timer
    toc = time.time()
    date = datetime.now().ctime()

    return {
        'duration': toc - tic,
        'date': date,
        'radiance': radiance,
        'irradiance': irradiance,
        'p_kelp': p_kelp
    }



def kelp_calculate_full(base_dir, db_path, table_name, absorptance_kelp, a_water, b, ns, nz, na, num_dens, kelp_dist, fs, fr, ft, max_length, length_std, zmax, rope_spacing, I0, phi_s, theta_s, decay, num_cores, num_scatters, fd_flag, lis_opts):
    # TODO: num_cores doesn't do anything yet.

    # TODO: Remove imports?
    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    zmin = 0
    dz = (zmax-zmin)/nz

    ns = int(ns)
    nz = int(nz)
    na = int(na)

    nx = ny = ns

    num_vsf = na
    vsf_angles = np.linspace(0,np.pi, na)
    vsf_vals = 0*vsf_angles + 1/(4*np.pi)

    ntheta = na
    nphi = na

    # nphi = int(na/2)
    # if nphi % 2 != 0:
    #     nphi += 1

    nomega = int(ntheta*(nphi-2)+2)
    p_kelp = np.asfortranarray(np.zeros([nx, ny, nz]))
    rad = np.asfortranarray(np.zeros([nx, ny, nz, nomega]))
    irrad = np.asfortranarray(np.zeros([nx, ny, nz]))

    # Start timer
    tic = time.time()

    # Rope spacing determines horizontal bounds
    xmin = ymin = -rope_spacing/2
    xmax = ymax = rope_spacing/2

    # Kelp shape
    frond_lengths, frond_stds, water_speeds, water_angles = get_kelp_dist(kelp_dist, max_length, zmin, zmax, nz)

    # Number of fronds in each depth layer from density
    # Keep constant over depth for now
    num_fronds = num_dens * dz * np.ones(nz)

    # Absorptance = % of light absorbed for whole frond (units 1).
    # a_kelp = absorption coefficient = %/m (units 1/m).
    a_kelp = absorptance_kelp / ft

    print("xmin = {}".format(xmin))
    print("xmax = {}".format(xmax))
    print("nx = {}".format(nx))
    print("ymin = {}".format(ymin))
    print("ymax = {}".format(ymax))
    print("ny = {}".format(ny))
    print("zmin = {}".format(zmin))
    print("zmax = {}".format(zmax))
    print("nz = {}".format(nz))
    print("frond_lengths.shape = {}".format(frond_lengths.shape))
    print("frond_stds.shape = {}".format(frond_stds.shape))
    print("num_fronds.shape = {}".format(num_fronds.shape))
    print("water_speeds.shape = {}".format(water_speeds.shape))
    print("water_angles.shape = {}".format(water_angles.shape))
    print("fs = {}".format(fs))
    print("fr = {}".format(fr))
    print("ft = {}".format(ft))
    print("p_kelp.shape = {}".format(p_kelp.shape))

    # Generate kelp
    f90.gen_kelp(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        frond_lengths,
        frond_stds,
        num_fronds,
        water_speeds,
        water_angles,
        fs, fr, ft,
        p_kelp
    )

    # Create arrays to hand mutable values to fortran
    lis_iter = np.array([0], dtype=int)
    lis_time = np.array([0], dtype=float)
    lis_resid = np.array([0], dtype=float)

    # Calculate light field
    f90.calculate_light_field(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        ntheta, nphi,
        a_water, a_kelp, b,
        vsf_angles, vsf_vals,
        theta_s, phi_s, I0, decay,
        p_kelp, rad, irrad,
        num_scatters, fd_flag, lis_opts,
        lis_iter, lis_time, lis_resid
    )

    avg_irrad = np.mean(irrad, axis=(0,1))
    perc_irrad = np.sum(p_kelp*irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))

    # End timer
    toc = time.time()
    date = datetime.now().ctime()
    git_commit = get_git_commit_hash()
    compute_time = toc - tic

    # Extract values from arrays
    lis_iter = lis_iter[0]
    lis_time = lis_time[0]
    lis_resid = lis_resid[0]

    # TODO: Separate DB stuff into another function/wrapper?
    params = {
        'absorptance_kelp': absorptance_kelp,
        'a_water': a_water,
        'b': b,
        # TODO: Should I choose either ns or nx, not both?
        'ns': ns,
        'na': na,
        'nx': nx,
        'ny': ny,
        'nz': nz,
        'ntheta': ntheta,
        'nphi': nphi,
        'nomega': nomega,
        'num_dens': num_dens,
        'kelp_dist': kelp_dist,
        'fs': fs,
        'fr': fr,
        'ft': ft,
        'max_length': max_length,
        'length_std': length_std,
        'zmax': zmax,
        'rope_spacing': rope_spacing,
        'I0': I0,
        'phi_s': phi_s,
        'theta_s': theta_s,
        'decay': decay,
        'num_cores': num_cores,
        'num_scatters': num_scatters,
        'fd_flag': fd_flag,
        'lis_opts': lis_opts,
        'date': date,
        'git_commit': git_commit,
        'compute_time': compute_time,
        'lis_iter': lis_iter,
        'lis_time': lis_time,
        'lis_resid': lis_resid,
    }

    # SQL has no bool type
    for var_name, val in params.items():
        if type(val) == bool:
            params[var_name] = int(val)

    # Create data file containing results
    nc_dict = {
        'p_kelp': p_kelp,
        'rad': rad,
        'irrad': irrad,
        'avg_irrad': avg_irrad,
        'perc_irrad': perc_irrad,
        **params,
    }
    data_path = create_nc(db_path, **nc_dict)


    # Save to DB
    db_dict = {
        'data_path': data_path,
        **params
    }
    insert_run(db_path, table_name, **db_dict)

def kelp_calculate(base_dir, db_path, table_name, a_water, b, ns, nz, na, kelp_dist, num_scatters, fd_flag, lis_opts='', num_cores=None):
    """kelp_calculate_full, but with some sensible defaults"""

    # TODO: Remove imports?
    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    absorptance_kelp = 0.8

    # Broch 2013
    # 150 individuals/meter
    num_dens = 150

    fs = 0.5
    # Handa figure 5
    fr = 5.0
    # From Solveig Foldal's Master's Thesis
    ft = 4e-4

    # From Solveig's Master's Thesis
    max_length = 6.0
    length_std = 0.2 * max_length

    zmax = 10 # Max. vertical
    rope_spacing = 10 # Horizontal

    # Fairly sunny day
    I0 = 50.0
    # Light from directly above
    phi_s = 0.0
    theta_s = 0.0
    decay = 1.0

    if not num_cores:
        num_cores = multiprocessing.cpu_count()

    return kelp_calculate_full(
        base_dir, db_path, table_name,
        absorptance_kelp, a_water, b,
        ns, nz, na, num_dens, kelp_dist,
        fs, fr, ft, max_length, length_std,
        zmax, rope_spacing,
        I0, phi_s, theta_s, decay,
        num_cores, num_scatters,
        fd_flag, lis_opts
    )

###############################
# Grid study

def grid_study_compute(study_name, a_water, b, kelp_dist, ns_list, nz_list, na_list, lis_opts, study_dir='.'):
    """
    Do grid study with Cartesian product of given resolutions.

    Store results in sqlite database
    `study_name` will be database name
    Database file located at {study_dir}/{study_name}/{study_name.db}
    Other data files located at {study_dir}/{study_name}
    """
    # Create database
    base_dir, db_path, table_name = create_table(study_name, study_dir)

    # TODO: Take `lv` as arg?
    lv = ipp.Client().load_balanced_view()

    # TODO: Should I choose either ns or nx, not both?
    # Constant values
    # One scatter before FD
    num_scatters = 1
    # TODO: THIS DEFEATS THE POINT OF THE GRID STUDY
    fd_flag = False

    # Futures from all computations
    # TODO: Are they necessary? (maybe for tracking progress)
    fut_dict = {}

    # Loop through grid
    # TODO: Is db_path necessary
    for ns, nz, na in it.product(ns_list, nz_list, na_list):
        fut_dict[(ns, nz, na)] = lv.apply(
            kelp_calculate,
            base_dir, db_path, table_name,
            a_water, b,
            ns, nz, na,
            kelp_dist, num_scatters,
            fd_flag, lis_opts)

    return fut_dict, db_path, table_name
