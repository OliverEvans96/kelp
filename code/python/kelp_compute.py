# stdlib
import sqlite3
from datetime import datetime
import subprocess
import multiprocessing
import os

# 3rd party
import netCDF4 as nc
import tempfile

###############
# Misc. utils #
###############

def get_random_unused_filename(dir='.', prefix='', suffix=''):
    with tempfile.NamedTemporaryFile(dir=dir, prefix=prefix, suffix=suffix) as fh:
        filename = fh.name

def get_git_commit_hash():
    """Get hash of current git branch."""
    return subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()
    return filename


#######################
# Kelp-specific funcs #
#######################

def create_table(run_name, prefix='.'):
    """
    Create sqlite db file, create table, and return connection.
    Assume run_name is unique (e.g. includes date.)
    """

    create_table_template = '''
    CREATE TABLE {db_name} (

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
        lis_time INTEGER,
        lis_resid REAL,
        /* Array data (NETCDF) */
        data_path CHAR(256)
        )
    '''

    base_dir = os.path.join(prefix, run_name)
    data_dir = os.path.join(base_dir, 'data')
    db_file = os.path.join(base_dir, run_name+'.db')
    # Will fail if `run_name` is not unique in `prefix`
    os.mkdir(base_dir)
    os.mkdir(data_dir)

    conn = sqlite3.connect(db_file)
    conn.execute(create_table_template.format(db_name=run_name))

    return conn


def get_first_db_name(conn):
    cursor = conn.execute("PRAGMA database_list;")
    # TODO: Might need to add another index here
    return cursor.fetch()

def insert_run(conn, **params):
    insert_template = '''
    INSERT INTO {db_name} VALUES (
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
    # If db_name is not provided,
    # just use the first one we find
    if 'db_name' not in params.keys():
        params['db_name'] = get_first_db_name(conn)

    insert_cmd = insert_template.format(**params)
    return conn.execute(insert_cmd)

def create_nc(data_path, **results):
    # Create file
    rootgrp = nc.Dataset(data_path, 'w', format='NETCDF4')

    # Get dimension sizes
    nx = len(results['x'])
    ny = len(results['y'])
    nz = len(results['z'])
    ntheta = len(results['theta'])
    nphi = len(results['phi'])

    # Create Dimensions
    x = rootgrp.createDimension('x', nx)
    y = rootgrp.createDimension('y', ny)
    z = rootgrp.createDimension('z', nz)
    theta = rootgrp.createDimension('theta', ntheta)
    phi = rootgrp.createDimension('phi', nphi)

    # Create dimension variables
    xvals = rootgrp.createVariable('x','f8', ('x',))
    yvals = rootgrp.createVariable('y','f8', ('y',))
    zvals = rootgrp.createVariable('z','f8', ('z',))
    thetavals = rootgrp.createVariable('theta','f8', ('theta',))
    phivals = rootgrp.createVariable('phi','f8', ('phi',))

    # Assign dimension variables
    xvals[:] = results.pop('x')
    yvals[:] = results.pop('y')
    zvals[:] = results.pop('z')
    thetavals[:] = results.pop('theta')
    phivals[:] = results.pop('phi')

    # Create results variables
    rad = rootgrp.createVariable('rad', 'f8', ('x', 'y', 'z', 'theta', 'phi'))
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
        var = rootgrp.createVariable(var_name, type(val))
        var[...] = val

    # Close file
    rootgrp.close()

def get_kelp_dist(kelp_dist, max_length):
    # TODO: Scale by zmax
    # Kelp distribution profiles
    # Normalize by max val
    maxval = 12*np.exp(-2) + 0.5
    frond_length_funcs = {
        'top-heavy': (3.0 * z**2 * np.exp(-z) + 0.5) / maxval,
        'bottom-heavy': (3.0 * (zmax-z)**2 * np.exp((z-zmax)) + 0.5) / maxval,
        'uniform': 0*z + 1.0,
        'none': 0*z
    }
    frond_lengths = max_length * frond_length_funcs[kelp_profile]
    frond_stds = 0.0 * np.ones(nz)
    water_speeds = 0.5 * np.ones(nz)
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



def kelp_calculate_full(conn, absorptance_kelp, a_water, b, ns, nz, na, num_dens, kelp_dist, fs, fr, ft, max_length, length_std, zmax, rope_spacing, I0, phi_s, theta_s, decay, num_cores, num_scatters, fd_flag, lis_opts):
    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    dz = (zmax-zmin)/nz

    num_vsf = na
    vsf_angles = np.linspace(0,np.pi, na)
    vsf_vals = 0*vsf_angles + 1/(4*np.pi)
    ns = int(ns)

    ntheta = na
    nphi = int(na/2)
    if nphi % 2 != 0:
        nphi += 1

    nomega = int(ntheta*(nphi-2)+2)
    p_kelp = np.asfortranarray(np.zeros([ns,ns,nz]))
    radiance = np.asfortranarray(np.zeros([ns, ns, nz, nomega]))
    irradiance = np.asfortranarray(np.zeros([ns, ns, nz]))

    # z grid centers
    z = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz)

    # Start timer
    tic = time.time()

    # Kelp shape
    frond_lengths, frond_stds, water_speeds, water_angles = get_kelp_dist(kelp_dist, max_length)

    # Number of fronds in each depth layer from density
    # Keep constant over depth for now
    num_fronds = num_dens * dz * np.ones(nz)

    # Absorptance = % of light absorbed for whole frond (units 1).
    # a_kelp = absorption coefficient = %/m (units 1/m).
    a_kelp = absorptance_kelp / ft

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

    # Calculate light field
    f90.calculate_light_field(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        ntheta, nphi,
        a_water, a_kelp, b,
        vsf_angles, vsf_vals,
        theta_s, phi_s, I0, decay,
        p_kelp, radiance, irradiance,
        num_scatters, fd_flag, lis_opts,
        lis_iter, lis_time, lis_resid
    )

    # End timer
    toc = time.time()
    date = datetime.now().ctime()
    git_commit = get_git_commit_hash()
    compute_time = toc - tic


    params = {
        'absorptance_kelp': absorptance_kelp,
        'a_water': a_water,
        'b': b,
        'ns': ns,
        'nz': nz,
        'na': na,
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
        'data_path': data_path
    }

    # Save to DB
    insert_run(conn, **params)

def kelp_calculate(conn, absorptance_kelp, a_water, b, ns, nz, na, kelp_dist, num_scatters, fd_flag, lis_opts='', num_cores=None):
    """kelp_calculate_full, but with some sensible defaults"""

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

    if num_cores = None:
        num_cores = multiprocessing.cpu_count()

    return kelp_calculate_full(
        conn, absorptance_kelp, a_water, b,
        ns, nz, na, num_dens, kelp_dist,
        fs, fr, ft, max_length, length_std,
        zmax, rope_spacing,
        I0, phi_s, theta_s, decay,
        num_cores, num_scatters,
        fd_flag, lis_opts
    )

###############################
# Grid study

def grid_study_compute(a_water, absorptance_kelp, kelp_profile, ns_list, nz_list, na_list):
    ns_max = max(ns_list)
    nz_max = max(nz_list)
    na_max = max(na_list)

    print("ns_max = {}".format(ns_max))
    print("nz_max = {}".format(nz_max))
    print("na_max = {}".format(na_max))

    # best solution (no kelp, no scattering)
    best_results = lv.apply(kelp_param.kelp_calculate,
        a_water,
        b,
        ns_max,
        na_max,
        nz_max,
        kelp_profile,
        absorptance_kelp,
        gmres_flag=True,
        num_scatters=4,
        const=const
    ).result()

    perceived_irrad_dict = {}
    duration_dict = {}
    abs_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])
    rel_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])

    p_kelp = best_results['p_kelp']
    best_irrad = best_results['irradiance']
    best_perceived_irrad = np.sum(p_kelp*best_irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    perceived_irrad_dict[(ns_max,nz_max,na_max)] = best_perceived_irrad

    # Vary nz
    ns = ns_max
    na = na_max
    for i, nz in enumerate(nz_list[:-1]):
        perceived_irrad, abs_err, rel_err, duration = compute_err(ns, nz, na, best_perceived_irrad, a_water, b, kelp_profile, absorptance_kelp, const)
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        duration_dict[(ns, nz, na)] = duration
        abs_err_arr[i, -1, -1] = abs_err
        rel_err_arr[i, -1, -1] = rel_err

    # Vary ns
    nz = nz_max
    na = na_max
    for i, ns in enumerate(ns_list[:-1]):
        perceived_irrad, abs_err, rel_err, duration = compute_err(ns, nz, na, best_perceived_irrad, a_water, b, kelp_profile, absorptance_kelp, const)
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        duration_dict[(ns, nz, na)] = duration
        abs_err_arr[-1, i, -1] = abs_err
        rel_err_arr[-1, i, -1] = rel_err

    # Vary na
    ns = ns_max
    nz = nz_max
    for i, na in enumerate(na_list[:-1]):
        perceived_irrad, abs_err, rel_err, duration = compute_err(ns, nz, na, best_perceived_irrad, a_water, b, kelp_profile, absorptance_kelp, const)
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        duration_dict[(ns, nz, na)] = duration
        abs_err_arr[-1, -1, i] = abs_err
        rel_err_arr[-1, -1, i] = rel_err

    return perceived_irrad_dict, abs_err_arr, rel_err_arr, duration_dict
