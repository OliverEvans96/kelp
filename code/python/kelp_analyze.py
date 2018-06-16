# 3rd party
import netCDF4 as nc
import pandas as pd

def table_to_df(conn, table_name):
    select_cmd = 'SELECT * FROM {table_name}'.format(table_name=table_name)
    cursor = conn.execute(select_cmd)
    columns = [col[0] for col in cursor.description]
    data = [x for x in cursor]

    df = pd.DataFrame(data, columns=columns)
    return df

def query_results(conn, table_name, **kwargs):
    # Sanitize string values
    for key, val in kwargs.items():
        if isinstance(val, str):
            kwargs[key] = '"{}"'.format(val)

    # Form SQL WHERE clause
    where_clause = ' AND '.join([
        '{}={}'.format(key, value)
        for key, value in kwargs.items()
    ])

    # Form entire query
    query = '''SELECT data_path FROM {table_name}
    WHERE {where_clause}
    '''.format(
        table_name=table_name,
        where_clause=where_clause
    )

    # Execute query
    cursor = conn.execute(query)
    results_list = [x for x in cursor]

    # Extract matching datasets
    datasets = [nc.Dataset(results[-1]) for results in results_list]
    return datasets

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


def kelp_visualize(duration, date, radiance, irradiance, p_kelp):
    print("Computation finished {}".format(date))
    print("Took {:.2f} seconds".format(duration))
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


#######################
## Grid Study

def merge_diff_grid(xmin, xmax, y1, y2):
    """Cast discrete quantities with different spacings to a common grid.
    Piecewise constant interpolation.
    """
    n1 = len(y1)
    n2 = len(y2)
    d1 = (xmax-xmin)/n1
    d2 = (xmax-xmin)/n2
    # Grid centers
    x1 = np.linspace(xmin+d1/2, xmax-d1/2, n1)
    x2 = np.linspace(xmin+d2/2, xmax-d2/2, n2)

    # Edges
    e1 = np.linspace(xmin, xmax, n1+1)
    e2 = np.linspace(xmin, xmax, n2+1)

    i_sort = np.argsort(np.concatenate([x1,x2]))
    i1 = i_sort<n1
    i2 = i_sort>=n1

    x3 = np.zeros(n1+n2)
    x3[i1] = x1
    x3[i2] = x2

    z1 = np.zeros_like(x3)
    z2 = np.zeros_like(x3)

    #z1[i1] = y1
    #z2[i2] = y2

    z1[x3<(x1[0]+x1[1])/2]
    for i in range(n1):
        wh = np.logical_and(
            e1[i]<=x3,
            x3<e1[i+1]
        )
        z1[wh] = y1[i]

    for i in range(n2):
        wh = np.logical_and(
            e2[i]<=x3,
            x3<e2[i+1]
        )
        z2[wh] = y2[i]

    return x3, z1, z2

def abs_err_uneven(xmin, xmax, y1, y2):
    "Integral of absolute difference between discrete quantities with different spacings."
    x, z1, z2 = merge_diff_grid(xmin, xmax, y1, y2)

    return np.sum(np.abs(
        (z1-z2) * np.diff(np.concatenate([x,[xmax]]))
    ))

def rel_err_uneven(xmin, xmax, y1, y2):
    "Integral of absolute difference between discrete quantities with different spacings."
    x, z1, z2 = merge_diff_grid(xmin, xmax, y1, y2)

    return np.sum(np.abs(
        (z1-z2)/z1 * np.diff(np.concatenate([x,[xmax]]))
    ))


def compute_block_err(ns, nz, na, best_irrad, a_water, b, kelp_profile, absorptance_kelp, const):
    print("ns={}, nz={}, na={}".format(ns,nz,na))
    compute_results = lv.apply(kelp_param.kelp_calculate,
        a_water,
        b,
        ns,
        na,
        nz,
        kelp_profile,
        absorptance_kelp,
        gmres_flag=True,
        num_scatters=4,
        const=const
    ).result()

    irrad = compute_results['irradiance'].mean(axis=(0,1))
    block_best_irrad = kelp_param.block_mean(best_irrad, irrad.shape)
    abs_err = np.mean(np.abs(irrad-block_best_irrad))
    rel_err = np.mean(np.abs((irrad-block_best_irrad)/block_best_irrad))
    return irrad, abs_err, rel_err

def compute_err(ns, nz, na, best_perceived_irrad, a_water, b, kelp_profile, absorptance_kelp, const):
    print("ns={}, nz={}, na={}".format(ns,nz,na))

    compute_results = lv.apply(kelp_param.kelp_calculate,
        a_water,
        b,
        ns,
        na,
        nz,
        kelp_profile,
        absorptance_kelp,
        gmres_flag=True,
        num_scatters=4,
        const=const
    ).result()

    p_kelp = compute_results['p_kelp']
    irrad = compute_results['irradiance']
    # Perceived irradiance for each depth layer
    perceived_irrad = np.sum(p_kelp*irrad, axis=(0,1)) / np.sum(p_kelp, axis=(0,1))
    # If p_kelp is 0, then so is perceied_irrad
    perceived_irrad[np.isnan(perceived_irrad)] = 0
    abs_err = abs_err_uneven(zmin, zmax, best_perceived_irrad, perceived_irrad)
    rel_err = rel_err_uneven(zmin, zmax, best_perceived_irrad, perceived_irrad)
    return perceived_irrad, abs_err, rel_err, compute_results['duration']

def slice_1d(n_list, best_ind, irrad_dict, pos, label):
    ns_max, nz_max, na_max = best_ind
    n_max = best_ind[pos]
    best_perceived_irrad = irrad_dict[(ns_max,nz_max,na_max)]

    dz = (zmax-zmin)/nz_max
    z_best = np.linspace(zmin+0.5*dz, zmax-0.5*dz, nz_max)

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
    plt.ylabel(r'Average Irradiance ($\mbox{W/m}^2)')

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
    plt.ylabel(r'Average Irradiance ($\mbox{W/m}^2)')

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

# TODO: args
def grid_study_analyze():
    """
    Analyze results from grid_study_compute.
    Compare all to best, calculate errors.
    """
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



# Plot Convergence Curves
def grid_study_plot(ns_list, nz_list, na_list, irrad_dict, abs_err_arr, rel_err_arr):
    ns_max = max(ns_list)
    nz_max = max(nz_list)
    na_max = max(na_list)

    nz = ns_max
    ns = ns_max

    best_ind = [ns_max, nz_max, na_max]
    for pos, (label,n_list) in enumerate(zip(['ns','nz','na'],(ns_list,nz_list,na_list))):
        slice_1d(n_list, best_ind, irrad_dict, pos, label)

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
