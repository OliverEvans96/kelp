# stdlib
import sqlite3
from datetime import datetime
import subprocess
import multiprocessing
import itertools as it
import functools as ft
import threading
import time
import os
import re

# 3rd party
import netCDF4 as nc
import tempfile
import ipyparallel as ipp
import numpy as np

###############
# Misc. utils #
###############

def get_random_unused_filename(dir='.', prefix='', suffix=''):
    with tempfile.NamedTemporaryFile(dir=dir, prefix=prefix, suffix=suffix, delete=False) as fh:
        filename = fh.name

    return filename

def get_git_commit_hash():
    """Get hash of current git branch."""
    return subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()
    return filename

def print_call(run_func, run_args, run_kwargs):
    func_name = run_func.__name__
    args_str = ', '.join([arg.__repr__() for arg in run_args])
    kwargs_str = ', '.join(['{}={}'.format(k,v.__repr__()) for k,v in run_kwargs.items()])
    sig_str = ', '.join([args_str, kwargs_str])
    print("Calling {}({})".format(func_name, sig_str))


#######################
# Kelp-specific funcs #
#######################

def create_table(conn, table_name, prefix='.'):
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

    conn.execute(create_table_template.format(table_name=table_name))
    conn.commit()

def create_dirs(study_dir):
    os.mkdir(study_dir)
    os.mkdir(os.path.join(study_dir, 'data'))

def get_table_names(conn):
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    return cursor.fetchall()

def insert_run(conn, table_name=None, **params):
    insert_template = '''
    INSERT INTO {table_name} VALUES (
        NULL, /* id (autoincrement) */
        :absorptance_kelp,
        :a_water,
        :b,
        :ns,
        :nz,
        :na,
        :num_dens,
        :kelp_dist,
        :fs,
        :fr,
        :ft,
        :max_length,
        :length_std,
        :zmax,
        :rope_spacing,
        :I0,
        :phi_s,
        :theta_s,
        :decay,
        :num_cores,
        :num_scatters,
        :fd_flag,
        :lis_opts,
        :date,
        :git_commit,
        :compute_time,
        :lis_iter,
        :lis_time,
        :lis_resid,
        :data_path
        );'''

    # If table_name is not provided,
    # just use the first one we find
    if not table_name:
        table_name = get_table_names(conn)[0]

    insert_command = insert_template.format(table_name=table_name)
    try:
        conn.execute(insert_command, params)
    except sqlite3.OperationalError as e:
        print('FAILURE WITH:')
        print("'{}'".format(insert_command))
        raise e

def combine_dbs(study_dir, table_name):
    data_dir = os.path.join(study_dir, 'data')
    dbs = [
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir) if re.match('.*\.db$', f)
    ]

    combined_db = os.path.join(study_dir, '{}.db'.format(table_name))
    combined_conn = sqlite3.connect(combined_db)
    create_table(combined_conn, table_name)

    for db in dbs:
        # Read from individual tables
        conn = sqlite3.connect(db)
        print("Combining {} (tables: {})".format(db, get_table_names(conn)))
        cursor = conn.execute('SELECT * FROM {}'.format(table_name))
        print("read.")
        columns = [tup[0] for tup in cursor.description]
        # Write to combined table
        for row_tuple in cursor:
            row_dict = dict(zip(columns, row_tuple))
            insert_run(combined_conn, table_name, **row_dict)
        conn.close()

    combined_conn.commit()
    combined_conn.close()

def create_nc(data_path, **results):
    """Create netCDF file to store results"""

    # Create file
    rootgrp = nc.Dataset(data_path, 'w', format='NETCDF4_CLASSIC')

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
        var_type = type(val)
        print("{}: {} ({})".format(var_name, val, var_type))
        # String handling for interoperability with fortran:
        # https://stackoverflow.com/a/37091930
        if var_type == str:
            # Strings treated as character arrays
            dim_name = '{}_dim'.format(var_name)
            char_dim = rootgrp.createDimension(dim_name, len(val))

            type_str = 'S{}'.format(len(val))
            char_val = nc.stringtochar(np.array([val], type_str))
            print("CREATE STR VAR: ('{}', '{}', '{}')".format(var_name, type_str, dim_name))
            var = rootgrp.createVariable(var_name, 'S1', (dim_name,))
            var[...] = char_val
        else:
            # Scalar int or float
            if var_type == float:
                type_str = 'f4'
            elif var_type == int:
                type_str = 'i4'
            var = rootgrp.createVariable(var_name, type_str)
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

def study_decorator(study_func):
    """Create directories before execution and
    merge dbs afterwards.

    Should be applied to functions like grid_study, etc.

    Function sholud return func_list, args_list, kwargs_list.
    Each of these will be applied via executor.apply().

    executor should be something like `ipp.Client().load_balanced_view()`
    with an `apply` method.

    Store results in sqlite database
    `study_name` will be database name
    Database file located at {base_dir}/{study_name}/{study_name.db}
    Other data files located at {base_dir}/{study_name}/data/*.nc
    """
    @ft.wraps(study_func)
    def wrapper(study_name, *study_args, base_dir=os.curdir, executor=None, **study_kwargs):

        study_dir = os.path.join(base_dir, study_name)
        if not executor:
            executor = ipp.Client().load_balanced_view()

        create_dirs(study_dir)
        study_calls = study_func(*study_args, **study_kwargs)

        run_futures = []

        # Execute function calls from study function
        for run_func, run_args, run_kwargs in it.zip_longest(*study_calls):
            # Make args, kwargs optional
            if not run_args:
                run_args = ()
            if not run_kwargs:
                run_kwargs = {}

            run_kwargs = {
                'study_dir': study_dir,
                'study_name': study_name,
                **run_kwargs
            }

            #print_call(run_func, run_args, run_kwargs)
            run_futures.append(executor.apply(run_func, *run_args, **run_kwargs))

        # Once all functions have run, combine the results
        def wait_and_combine():
            for future in run_futures:
                future.wait()
                print("{} futures done.".format(sum([f.done() for f in run_futures])))
            combine_dbs(study_dir, study_name)
        combine_thread = threading.Thread(target=wait_and_combine)
        combine_thread.start()

        # This will be returned immediately.
        # It can be probed with combined_thread.isAlive()
        return combine_thread, run_futures

    return wrapper

def run_decorator(run_func):
    """Run function and save results to database.
    The function is expected to return a dictionary
    of all values to be stored in the database.

    The .nc file will be written to data/`filename`.nc,
    and the temporary .db file will be written to `filename`.db,
    where `filename` is randomly generated.
    All .db files should later be merged into one db per directory.
    """
    @ft.wraps(run_func)
    def wrapper(*args, study_dir=None, study_name=None, **kwargs):
        scalar_params, results = run_func(*args, **kwargs)

        if not study_dir:
            raise ValueError("kwarg `study_dir` required for functions wrapped by `run_decorator`")
        if not study_name:
            raise ValueError("kwarg `study_name` required for functions wrapped by `run_decorator`")

        # Generate random name in correct directory
        data_path = get_random_unused_filename(
            dir=os.path.join(study_dir, 'data'),
            suffix='.nc'
        )
        # Change suffix from .nc to .db for db path.
        db_path = re.sub('\.nc$', '.db', data_path)

        # SQL has no bool type
        for var_name, val in scalar_params.items():
            if type(val) == bool:
                scalar_params[var_name] = int(val)

        # Create data file containing results
        nc_dict = {
            **scalar_params,
            **results
        }
        data_path = create_nc(data_path, **nc_dict)

        # Save to DB
        db_dict = {
            'data_path': data_path,
            **scalar_params
        }

        conn = sqlite3.connect(db_path)
        create_table(conn, study_name)
        insert_run(conn, study_name, **db_dict)
        conn.commit()
        conn.close()

    return wrapper

###############

def kelp_calculate_full(absorptance_kelp, a_water, b, ns, nz, na, num_dens, kelp_dist, fs, fr, ft, max_length, length_std, zmax, rope_spacing, I0, phi_s, theta_s, decay, num_cores, num_scatters, fd_flag, lis_opts):
    # TODO: num_cores doesn't do anything yet.

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
    avg_irrad = np.asfortranarray(np.zeros([nz], dtype=np.float32))
    perc_irrad = np.asfortranarray(np.zeros([nz], dtype=np.float32))

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

    # TODO: Remove?
    print("a_kelp = {}".format(a_kelp))
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
        p_kelp, rad, irrad, avg_irrad, perc_irrad,
        num_scatters, fd_flag, lis_opts,
        lis_iter, lis_time, lis_resid
    )

    # End timer
    toc = time.time()
    date = datetime.now().ctime()
    git_commit = get_git_commit_hash()
    compute_time = toc - tic

    # Extract values from arrays
    lis_iter = int(lis_iter)
    lis_time = float(lis_time)
    lis_resid = float(lis_resid)

    scalar_params = {
        'absorptance_kelp': absorptance_kelp,
        'a_water': a_water,
        'b': b,
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

    results = {
        'p_kelp': p_kelp,
        'rad': rad,
        'irrad': irrad,
        'avg_irrad': avg_irrad,
        'perc_irrad': perc_irrad,
    }

    return scalar_params, results

@run_decorator
def kelp_calculate(a_water, b, ns, nz, na, kelp_dist, num_scatters, fd_flag, lis_opts='', num_cores=None):
    """kelp_calculate_full, but with some sensible defaults, saving results to .db/.nc due to wrapper"""

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
    rope_spacing = 15 # Horizontal

    # Fairly sunny day
    I0 = 50.0
    # Light from directly above
    phi_s = 0.0
    theta_s = 0.0
    decay = 1.0

    if not num_cores:
        num_cores = multiprocessing.cpu_count()

    return kelp_calculate_full(
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

@study_decorator
def grid_study_compute(a_water, b, kelp_dist, ns_list, nz_list, na_list, lis_opts):
    """
    Do grid study with Cartesian product of given resolutions.
    """

    # One scatter before FD
    num_scatters = 1

    fd_flag = True

    # Actual calling will be performed by decorator.
    # Functions to be called
    func_list = []
    # Arguments to be passed
    args_list = [] # tuples/lists
    kwargs_list = [] # dictionaries

    # Loop through grid
    for ns, nz, na in it.product(ns_list, nz_list, na_list):
        func_list.append(kelp_calculate)
        args_list.append((
            a_water, b,
            ns, nz, na,
            kelp_dist, num_scatters,
            fd_flag, lis_opts
        ))

    return func_list, args_list, kwargs_list
