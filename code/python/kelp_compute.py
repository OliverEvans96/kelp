# stdlib
import sqlite3
from datetime import datetime
import subprocess
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
        date CHAR(64),
        compute_time REAL,
        git_commit CHAR(40),
        num_cores INTEGER,
        num_scatters INTEGER,
        fd_flag INTEGER,
        lis_opts CHAR(256),
        lis_iter INTEGER,
        lis_time INTEGER,
        lis_resid REAL,

        /* Array data (NETCDF)
        Contains:
          - radiance (ns, ns, nz, na, na)
          - irradiance (ns, ns, nz)
          - p_kelp (ns, ns, nz)
          - avg_irrad (nz)
          - perc_irrad (nz)
        */
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
        '{date}',
        {compute_time},
        '{git_commit}',
        {num_cores},
        {num_scatters},
        {fd_flag},
        '{lis_opts}',
        {lis_iter},
        {lis_time},
        {lis_resid},
        '{data_path}'
        )
    '''
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

###############

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
        'irradiance': irradiance
    }

def nokelp_visualize(duration, date, radiance, irradiance):
    print("Computation finished {}".format(date))
    print("Took {:.2f} seconds".format(duration))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig = ipv.figure()
        ipv.volshow(irradiance.T, controls=False)

    # Set initial angle
    fig.anglex = 0.85 * np.pi
    fig.angley = -0.30 * np.pi
    fig.anglez = -0.67 * np.pi
    ipv.show()

def kelp_calculate(a_water, b, ns, na, nz, kelp_profile, absorptance_kelp, gmres_flag, num_scatters, const):
    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    # Extract constants
    (rope_spacing, zmin, zmax, I0, phi_s, theta_s, decay, xmin, xmax, ymin, ymax,
         lic_options) = const

    dz = (zmax-zmin)/nz
    print("dz = {}".format(dz))

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


    # Kelp distribution profiles
    frond_length_funcs = {
        'top-heavy': 3.0 * z**2 * np.exp(-z) + 0.5,
        'bottom-heavy': 3.0 * (zmax-z)**2 * np.exp((z-zmax)) + 0.5,
        'uniform': 0*z + 1.0,
        'none': 0*z
    }

    # Top-heavy
    frond_lengths = frond_length_funcs[kelp_profile]
    #frond_lengths = np.ones(nz)
    frond_stds = 0.0 * np.ones(nz)
    num_fronds = 10 * np.ones(nz)
    water_speeds = 0.5 * np.ones(nz)
    water_angles = 2*np.pi / zmax * (z-zmin)

    fs = 0.5
    fr = 0.5
    ft = 1e-2

    print("theoretical max_kelp = {}".format(num_fronds.max()*ft/dz))

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

    print("Max kelp = {}".format(p_kelp.max()))

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
