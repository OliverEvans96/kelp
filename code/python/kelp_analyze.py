# stdlib
import os

# 3rd party
import sqlite3
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# local
import uneven_diff

def table_to_df(conn, table_name):
    select_cmd = 'SELECT * FROM {table_name}'.format(table_name=table_name)
    cursor = conn.execute(select_cmd)
    columns = [col[0] for col in cursor.description]
    data = [x for x in cursor]

    df = pd.DataFrame(data, columns=columns)
    return df

def query_results(conn, table_name, **kwargs):
    for key, val in kwargs.items():
        # Sanitize string values
        if isinstance(val, str):
            kwargs[key] = '"{}"'.format(val)
        
        # Sanitize bool values (convert to int)
        if isinstance(val, bool):
            kwargs[key] = '{}'.format(int(val))

    # Form SQL WHERE clause
    where_condition = ' AND '.join([
        '{}={}'.format(key, value)
        for key, value in kwargs.items()
    ])
    if where_condition:
        where_clause = 'WHERE ' + where_condition
    else:
        where_clause = ''
        
    # Form entire query
    query = ' '.join([
        'SELECT data_path FROM {table_name}',
        '{where_clause}'
    ]).format(
        table_name=table_name,
        where_clause=where_clause
    )
    
    #print("query: '{}'".format(query))

    # Execute query
    cursor = conn.execute(query)
    results_list = cursor.fetchall()

    # Extract matching datasets
    datasets = [nc.Dataset(results[-1]) for results in results_list]
    return datasets

def nokelp_visualize(compute_time, date, radiance, irradiance):
    print("Computation finished {}".format(date))
    print("Took {:.2f} seconds".format(compute_time))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig = ipv.figure()
        ipv.volshow(irradiance.T, controls=False)

    # Set initial angle
    fig.anglex = 0.85 * np.pi
    fig.angley = -0.30 * np.pi
    fig.anglez = -0.67 * np.pi
    ipv.show()


def kelp_visualize(compute_time, date, radiance, irradiance, p_kelp):
    print("Computation finished {}".format(date))
    print("Took {:.2f} seconds".format(compute_time))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kelpfig = ipv.figure()
        ipv.volshow(p_kelp.T, controls=False)
        # set initial angle
        kelpfig.anglex = 0.85 * np.pi
        kelpfig.angley = -0.30 * np.pi
        kelpfig.anglez = -0.67 * np.pi

        irradfig = ipv.figure()
        ipv.volshow(irradiance.T, controls=False)
        # set initial angle
        irradfig.anglex = 0.85 * np.pi
        irradfig.angley = -0.30 * np.pi
        irradfig.anglez = -0.67 * np.pi


    rad_text = "Rad Stats<br>---------<br>min: {:.3e}<br>max: {:.3e}<br>mean: {:.3e}<br>std: {:.3e}".format(
        radiance.min(),
        radiance.max(),
        radiance.mean(),
        radiance.std()
    )

    irrad_text = "Irrad Stats<br>-----------<br>min: {:.3e}<br>max: {:.3e}<br>mean: {:.3e}<br>std: {:.3e}".format(
        irradiance.min(),
        irradiance.max(),
        irradiance.mean(),
        irradiance.std()
    )

    display(ipw.HBox([
        kelpfig,
        irradfig,
        ipw.Box(layout=ipw.Layout(width='10px')),
        ipw.HTML("<tt>{}<br><br>{}</tt>".format(rad_text, irrad_text))
    ]))


def block_mean(large_arr, small_shape):
    """Calculate an array of block means of `large_arr` which has the shape `small_shape`"""

    if all(n_large % n_small == 0 for n_small, n_large in zip(small_shape, large_arr.shape)):
        # Try to abstract over number of dimensions
        avg = np.zeros(small_shape)
        block_sizes = [n_large // n_small for n_small, n_large in zip(small_shape, large_arr.shape)]
        #print("block_sizes = {}".format(block_sizes))

        # Loop over all combinations of indices
        for inds in it.product(*(range(n) for n in small_shape)):
            #print("inds = {}".format(inds))
            startinds = [ind * block_size for ind, block_size in zip(inds, block_sizes)]
            stopinds = [(ind+1) * block_size for ind, block_size in zip(inds, block_sizes)]
            slices = tuple(slice(startind, stopind) for startind, stopind in zip(startinds, stopinds))
            #print("startinds = {}".format(startinds))
            #print("stopinds = {}".format(stopinds))
            avg[inds] = large_arr[slices].mean()

        return avg

    else:
        raise IndexError("`small_shape` must divide `large_arr.shape` elementwise.")


################
## Grid Study ##
################

def calculate_perceived_irrad(p_kelp, irrad):
    perceived_irrad = np.sum(p_kelp*irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    # TODO: Not true.
    # If p_kelp is 0, then so is perceied_irrad
    perceived_irrad[np.isnan(perceived_irrad)] = 0
    
    return perceived_irrad

def compute_err(conn, table_name, best_perceived_irrad, **run_dict):
    #print("ns={}, nz={}, na={}".format(ns,nz,na))

    compute_results = query_results(
        conn,
        table_name,
        **run_dict
    )[0]
    # This one (not best)
    nz = int(compute_results['nz'][:])
    zmin = 0
    zmax = compute_results['zmax'][:]
    p_kelp = compute_results['p_kelp'][:]
    irrad = compute_results['irrad'][:]
    perceived_irrad = calculate_perceived_irrad(p_kelp, irrad)
    dz = (zmax - zmin) / nz
    z = dz * (np.arange(nz) + 0.5)
    
    # Best
    nz_best = len(best_perceived_irrad)
    dz_best = (zmax - zmin) / nz_best
    z_best = dz_best * (np.arange(nz_best) + 0.5)
    
    # Piecewise-constant error
    #abs_err, rel_err = uneven_diff.discrete_err(zmin, zmax, best_perceived_irrad, perceived_irrad)
    
    # Linear interpolation error
    abs_err, rel_err = uneven_diff.err_linear(
        zmin, zmax,
        z_best, z,
        best_perceived_irrad, 
        perceived_irrad
    )
    
    
    return perceived_irrad, abs_err, rel_err, compute_results['compute_time']

def cori_plot_two_avg_irrads(study_name, grid1, grid2):
    """
    Plot both perceived_irrads. rel_err is w.r.t. tup1.
    """
    base_dir = os.path.join(os.environ['SCRATCH'], 'kelp-results')
    study_dir = os.path.join(base_dir, study_name)
    db_path = os.path.join(study_dir, '{}.db'.format(study_name))
    
    conn = sqlite3.connect(db_path)
    
    results1 = query_results(
        conn, 
        study_name, 
        **grid1
    )[0]
    results2 = query_results(
        conn, 
        study_name, 
        **grid2
    )[0]
    
    conn.close()
    
    zmin = 0
    zmax = results1['zmax'][:]
    
    p_kelp1 = results1['p_kelp'][:]
    irrad1 = results1['irrad'][:]
    perceived_irrad1 = calculate_perceived_irrad(p_kelp1, irrad1)
    
    p_kelp2 = results2['p_kelp'][:]
    irrad2 = results2['irrad'][:]
    perceived_irrad2 = calculate_perceived_irrad(p_kelp2, irrad2)
    
    nz1 = grid1['nz']
    dz1 = (zmax - zmin) / nz1
    z1 = dz1 * (np.arange(nz1) + 0.5)
    
    nz2 = grid2['nz']
    dz2 = (zmax - zmin) / nz2
    z2 = dz2 * (np.arange(nz2) + 0.5)
    
    # # Piecewise-constant interpolation
    # abs_err, rel_err = uneven_diff.discrete_err(
    #     zmin, zmax,
    #     perceived_irrad1,
    #     perceived_irrad2,
    #     verbose=True
    # )
    
    # Linear interpolation
    abs_err, rel_err = uneven_diff.err_linear(
        zmin, zmax,
        z1, z2,
        perceived_irrad1,
        perceived_irrad2,
    )
    
    plt.plot(z1, perceived_irrad1, 'o-', label='perc_irrad1')
    plt.plot(z2, perceived_irrad2, 'o-', label='perc_irrad2')
    plt.xlabel('z')
    plt.ylabel('perceived irradiance')
    plt.legend()
    
    print("abs_err = {:.3e}".format(abs_err))
    print("rel_err = {:.3e}".format(rel_err))
    
def cori_plot_two_avg_irrads_onespace(study_name, grid1, grid2):
    grid1['nz'] = grid1['ns']
    grid2['nz'] = grid2['ns']
    cori_plot_two_avg_irrads(study_name, grid1, grid2)
    
def plot_slice_1d(n_list, best_ind, irrad_dict, pos, label, zmin, zmax):
    ns_max, nz_max, na_max = best_ind
    n_max = best_ind[pos]
    best_perceived_irrad = irrad_dict[(ns_max,nz_max,na_max)]

    dz = (zmax-zmin)/nz_max
    z_best = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz_max)

    print(label)

    plt.figure(figsize=[8,6])

    # Log Average Irradiance
    ind = best_ind[:]
    plt.subplot(2,2,1)
    for n in n_list[:-1]:
        ind[pos] = n
        nz = ind[1]
        dz = (zmax-zmin)/nz
        z = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz)
        irrad = irrad_dict[tuple(ind)]
        plt.plot(z, irrad, '-o', label='{}={}'.format(label,n))
    plt.plot(z_best, best_perceived_irrad, '-s', lw=3, label='{}={}'.format(label,n_max))
    plt.yscale('log')
    plt.xlim(zmin, zmax)
    #plt.ylim(1e-2, I0)
    #plt.legend()
    plt.xlabel('z (depth)')
    plt.ylabel(r'Average Irradiance ($\mathrm{W/m}^2$)')

    # Linear Average Irradiance
    ind = best_ind[:]
    plt.subplot(2,2,2)
    for n in n_list[:-1]:
        ind[pos] = n
        nz = ind[1]
        dz = (zmax-zmin)/nz
        z = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz)
        irrad = irrad_dict[tuple(ind)]
        plt.plot(z, irrad, '-o', label='{}={}'.format(label,n))
    plt.plot(z_best, best_perceived_irrad, '-s', lw=3, label='{}={}'.format(label,n_max))
    plt.yscale('linear')
    plt.xlim(zmin, zmax)
    #plt.ylim(0, I0)
    plt.legend(loc='upper left', bbox_to_anchor=(1.1,1))
    plt.xlabel('z (depth)')
    plt.ylabel(r'Average Irradiance ($\mathrm{W/m}^2$)')

    # Absolute Error vs depth
    ind = best_ind[:]
    plt.subplot(2,2,3)
    for n in n_list[:-1]:
        ind[pos] = n
        nz = ind[1]
        dz = (zmax-zmin)/nz
        z = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz)
        irrad = irrad_dict[tuple(ind)]
        this_interp = interp1d(z, irrad, fill_value='extrapolate')
        best_interp = interp1d(z_best, best_perceived_irrad)
        diff_interp = lambda z: this_interp(z) - best_interp(z)
        z_both = np.sort(np.concatenate([z, z_best]))
        abs_err = np.abs(best_interp(z_both)-this_interp(z_both))
        plt.plot(z_both, abs_err, '-o', label='{}={}'.format(label,n))
    plt.yscale('log')
    plt.xlim(zmin, zmax)
    #plt.legend()
    plt.xlabel('z (depth)')
    plt.ylabel('Absolute Error')

    # Relative Error vs depth
    ind = best_ind[:]
    plt.subplot(2,2,4)
    for n in n_list[:-1]:
        ind[pos] = n
        nz = ind[1]
        dz = (zmax-zmin)/nz
        z = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz)
        irrad = irrad_dict[tuple(ind)]
        this_interp = interp1d(z, irrad, fill_value='extrapolate')
        best_interp = interp1d(z_best, best_perceived_irrad)
        z_both = np.sort(np.concatenate([z, z_best]))
        rel_err = np.abs((best_interp(z_both)-this_interp(z_both))/best_interp(z_both))
        plt.plot(z_both, rel_err, '-o', label='{}={}'.format(label,n))
    plt.yscale('log')
    plt.xlim(zmin, zmax)
    #plt.legend()
    plt.xlabel('z (depth)')
    plt.ylabel('Relative Error')

    plt.tight_layout()
    plt.show()

def grid_study_analyze_full(db_path, table_name):
    """
    Analyze results from grid_study_compute.
    Compare all to best, calculate errors.
    """
    conn = sqlite3.connect(db_path)

    cursor = conn.execute('''
    SELECT ns, nz, na
    FROM {table_name}
    '''.format(table_name=table_name))

    # Get all unique values of ns, nz, na
    ns_list, nz_list, na_list = (
        sorted(map(int, set(z))) for z in zip(*cursor.fetchall())
    )

    ns_max = max(ns_list)
    nz_max = max(nz_list)
    na_max = max(na_list)

    # Assume there's only one run
    # that matches largest grid
    best_results = query_results(
        conn,
        table_name,
        ns=ns_max,
        nz=nz_max,
        na=na_max,
        fd_flag=True
    )[0]

    perceived_irrad_dict = {}
    compute_time_dict = {}
    abs_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])
    rel_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])

    p_kelp = best_results['p_kelp'][:]
    best_irrad = best_results['irrad'][:]
    best_perceived_irrad = np.sum(p_kelp*best_irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    perceived_irrad_dict[(ns_max,nz_max,na_max)] = best_perceived_irrad

    for i, ns in enumerate(ns_list):
        for j, nz in enumerate(nz_list):
            for k, na in enumerate(na_list):
                perceived_irrad, abs_err, rel_err, compute_time = compute_err(
                    conn, 
                    table_name, 
                    best_perceived_irrad,
                    ns=ns,
                    nz=nz,
                    na=na
                )
                perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
                compute_time_dict[(ns, nz, na)] = compute_time
                abs_err_arr[i, j, k] = abs_err
                rel_err_arr[i, j, k] = rel_err
    return perceived_irrad_dict, abs_err_arr, rel_err_arr, compute_time_dict

def grid_study_analyze_onespace(db_path, table_name):
    """
    Analyze results from grid_study_compute.
    Compare all to best, calculate errors.
    """
    conn = sqlite3.connect(db_path)

    cursor = conn.execute('''
    SELECT ns, na
    FROM {table_name}
    '''.format(table_name=table_name))

    # Get all unique values of ns, nz, na
    ns_list, na_list = (
        sorted(map(int, set(z))) for z in zip(*cursor.fetchall())
    )

    ns_max = max(ns_list)
    na_max = max(na_list)

    # Assume there's only one run
    # that matches largest grid + FD
    best_results = query_results(
        conn,
        table_name,
        ns=ns_max,
        na=na_max,
        fd_flag=True
    )[0]

    perceived_irrad_dict = {}
    compute_time_dict = {}
    abs_err_arr = np.zeros([len(ns_list),len(na_list)])
    rel_err_arr = np.zeros([len(ns_list),len(na_list)])

    p_kelp = best_results['p_kelp'][:]
    best_irrad = best_results['irrad'][:]
    best_perceived_irrad = np.sum(p_kelp*best_irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    perceived_irrad_dict[(ns_max,na_max)] = best_perceived_irrad

    for i, ns in enumerate(ns_list):
        for k, na in enumerate(na_list):
            perceived_irrad, abs_err, rel_err, compute_time = compute_err(
                conn, 
                table_name, 
                best_perceived_irrad,
                ns=ns,
                na=na,
                fd_flag=True
            )
            perceived_irrad_dict[(ns, na)] = perceived_irrad
            compute_time_dict[(ns, na)] = compute_time
            abs_err_arr[i, k] = abs_err
            rel_err_arr[i, k] = rel_err
            
    # Set errors for finest grid to NaN 
    # (since there's nothing to compare it to)
    abs_err_arr[-1, -1] = np.nan
    rel_err_arr[-1, -1] = np.nan
    
    return perceived_irrad_dict, abs_err_arr, rel_err_arr, compute_time_dict

def grid_study_analyze_edges(db_path, table_name):
    """
    Analyze results from grid_study_compute.
    Compare all to best, calculate errors.
    """
    conn = sqlite3.connect(db_path)

    cursor = conn.execute('''
    SELECT ns, nz, na
    FROM {table_name}
    '''.format(table_name=table_name))

    # Get all unique values of ns, nz, na
    ns_list, nz_list, na_list = (
        sorted(map(int, set(z))) for z in zip(*cursor.fetchall())
    )

    ns_max = max(ns_list)
    nz_max = max(nz_list)
    na_max = max(na_list)

    # Assume there's only one run
    # that matches largest grid
    best_results = query_results(
        conn,
        table_name,
        ns=ns_max,
        nz=nz_max,
        na=na_max
    )[0]

    perceived_irrad_dict = {}
    compute_time_dict = {}
    abs_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])
    rel_err_arr = np.zeros([len(ns_list),len(nz_list),len(na_list)])

    p_kelp = best_results['p_kelp'][:]
    best_irrad = best_results['irrad'][:]
    best_perceived_irrad = np.sum(p_kelp*best_irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    perceived_irrad_dict[(ns_max,nz_max,na_max)] = best_perceived_irrad

    # Vary nz
    ns = ns_max
    na = na_max
    for i, nz in enumerate(nz_list[:-1]):
        perceived_irrad, abs_err, rel_err, compute_time = compute_err(
            conn, 
            table_name, 
            best_perceived_irrad,
            ns=ns,
            nz=nz,
            na=na,
            fd_flag=True
        )
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        compute_time_dict[(ns, nz, na)] = compute_time
        abs_err_arr[i, -1, -1] = abs_err
        rel_err_arr[i, -1, -1] = rel_err

    # Vary ns
    nz = nz_max
    na = na_max
    for i, ns in enumerate(ns_list[:-1]):
        perceived_irrad, abs_err, rel_err, compute_time = compute_err(
            conn, 
            table_name, 
            best_perceived_irrad,
            ns=ns,
            nz=nz,
            na=na
        )
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        compute_time_dict[(ns, nz, na)] = compute_time
        abs_err_arr[-1, i, -1] = abs_err
        rel_err_arr[-1, i, -1] = rel_err

    # Vary na
    ns = ns_max
    nz = nz_max
    for i, na in enumerate(na_list[:-1]):
        perceived_irrad, abs_err, rel_err, compute_time = compute_err(
            conn, 
            table_name, 
            best_perceived_irrad,
            ns=ns,
            nz=nz,
            na=na
        )
        perceived_irrad_dict[(ns, nz, na)] = perceived_irrad
        compute_time_dict[(ns, nz, na)] = compute_time
        abs_err_arr[-1, -1, i] = abs_err
        rel_err_arr[-1, -1, i] = rel_err

    return perceived_irrad_dict, abs_err_arr, rel_err_arr, compute_time_dict

# Plot Convergence Curves
def grid_study_plot(irrad_dict, abs_err_arr, rel_err_arr, compute_time_dict, zmin, zmax):
    # TODO: plot compute_time_dict

    # Get all unique values of ns, nz, na from dict keys
    ns_list, nz_list, na_list = (
        sorted(map(int, set(z))) for z in zip(*irrad_dict.keys())
    )

    ns_max = max(ns_list)
    nz_max = max(nz_list)
    na_max = max(na_list)

    best_ind = [ns_max, nz_max, na_max]
    for pos, (label, n_list) in enumerate(zip(['ns','nz','na'],(ns_list,nz_list,na_list))):
        plot_slice_1d(n_list, best_ind, irrad_dict, pos, label, zmin, zmax)

    print("together")
    # Relative error vs resolution
    fig, ax1 = plt.subplots()
    ax1.plot(ns_list[:-1], rel_err_arr[:-1,-1,-1], 'o-', label='ns')
    ax1.plot(nz_list[:-1], rel_err_arr[-1,:-1,-1], 'o-', label='nz')
    ax1.plot(na_list[:-1], rel_err_arr[-1,-1,:-1], 'o-', label='na')
    ax1.set_xlabel('Number of grid points')
    ax1.set_ylabel('Relative Error')
    ax1.set_yscale('log')
    ax1.legend(loc='upper left')
    plt.tight_layout()

    plt.show()
