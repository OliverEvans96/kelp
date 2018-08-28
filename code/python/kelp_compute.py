"""kelp_compute.py"""

# stdlib
from datetime import datetime
import itertools as it
import multiprocessing
import subprocess
import threading
import json
import inspect
import sqlite3
import time
import os
import re

# 3rd party
import netCDF4 as nc
import tempfile
import ipyparallel as ipp
import numpy as np
import sympy as sp
import functools
# NOTE: Using boltons.funcutils.wraps in place
# of functools.wraps in order to correctly preserve
# signature of wrapper for the sake of f2py
# see https://hynek.me/articles/decorators/
import boltons.funcutils as fu

# local
import run_utils as ru
import mms

## Kelp-specific funcs ##

def sym_to_num(fun, *args):
    """
    Convert sympy function to numpy function,
    with the output shape broadcasted to the shape
    of the sum of all arguments.

    This is required in case one or more arguments
    are not used explicity in the formula.
    """

    f = fu.wraps(fun)(sp.lambdify(
        args,
        fun(*args),
        modules=("numpy",)
    ))

    @fu.wraps(fun)
    def wrapper(*inner_args):
        """
        Reshape output to always match broadcasted
        sum of inputs, even if they are not all
        explicitly used in the function.
        """
        array_args = map(np.array, inner_args)
        shape = np.shape(sum(array_args))
        ans = f(*inner_args)
        return np.broadcast_to(ans, shape)

    return wrapper

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
        'huge': 0*z + 20,
        'none': 0*z
    }
    frond_lengths = max_length * frond_length_funcs[kelp_dist]
    frond_stds = 0.0 * np.ones_like(z)
    water_speeds = 0.5 * np.ones_like(z)
    if kelp_dist == 'huge':
        water_speeds *= 20

    water_angles = 2*np.pi / zmax * (z-zmin)

    return frond_lengths, frond_stds, water_speeds, water_angles


## Run Functions ##

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
    vsf_angles = np.linspace(0, np.pi, na)
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
        p_kelp, rad, irrad,
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

@ru.run_decorator
def kelp_calculate(a_water, b, ns, nz, na, kelp_dist, num_scatters, fd_flag, lis_opts='', num_cores=None):
    """kelp_calculate_full, but with some sensible defaults, saving results to .db/.nc due to wrapper"""

    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time

    absorptance_kelp = 0.7

    # Broch 2013
    # 150 individuals/meter
    num_dens = 120

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
        absorptance_kelp, a_water, b,
        ns, nz, na, num_dens, kelp_dist,
        fs, fr, ft, max_length, length_std,
        zmax, rope_spacing,
        I0, phi_s, theta_s, decay,
        num_cores, num_scatters,
        fd_flag, lis_opts
    )


def solve_rte_with_callbacks_full(ns, nz, ntheta, nphi, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict, num_scatters, fd_flag, lis_opts):
    # TODO: num_cores doesn't do anything yet.

    from kelp3d_objs import f90
    import numpy as np
    from datetime import datetime
    import time
    import sympy as sp

    # Get symbolic expressions for callbacks as strings
    space = sp.var('x, y, z')
    x, y, z = space
    angle = sp.var('theta, phi')
    theta, phi = angle
    delta = sp.var('Delta')

    # NOTE: sol_expr is not actually used in solution procedure,
    # it's just required so that it can be stored for future reference.
    # If the true solution is not known, just pass an empty string.
    sol_expr_str = str(sol_expr)
    abs_expr_str = str(abs_expr)
    source_expr_str = str(source_expr)
    bc_expr_str = str(bc_expr)
    vsf_expr_str = str(vsf_expr)

    # Convert parameter values dictionary to string
    param_dict_str = json.dumps(param_dict)

    # Convert expressions to sympy functions
    source_sym = mms.symify(source_expr, *space, *angle, **param_dict)
    abs_sym = mms.symify(abs_expr, *space, **param_dict)
    bc_sym = mms.symify(bc_expr, *angle, **param_dict)
    vsf_sym = mms.symify(vsf_expr, delta, **param_dict)

    # Convert sympy functions to numpy functions
    abs_func_N = sym_to_num(abs_sym, *space)
    source_func_N = sym_to_num(source_sym, *space, *angle)
    bc_func_N = sym_to_num(bc_sym, *angle)
    vsf_func_N = sym_to_num(vsf_sym, delta)

    # Calculate source expansion
    source_expansion_N = mms.gen_series_N(source_expr, num_scatters, **param_dict)

    # Assign grid variables
    zmin = 0
    dz = (zmax-zmin)/nz

    ns = int(ns)
    nz = int(nz)
    ntheta = int(ntheta)
    nphi = int(nphi)

    #na = int(na)
    #ntheta = na
    #nphi = na

    nx = ny = ns

    nomega = int(ntheta*(nphi-2)+2)

    # Initialize solution arrays
    rad = np.asfortranarray(np.zeros([nx, ny, nz, nomega]))
    irrad = np.asfortranarray(np.zeros([nx, ny, nz]))

    # Start timer
    tic = time.time()

    # Rope spacing determines horizontal bounds
    xmin = ymin = -rope_spacing/2
    xmax = ymax = rope_spacing/2

    # Create arrays to hand mutable values to fortran
    lis_iter = np.array([0], dtype=int)
    lis_time = np.array([0], dtype=float)
    lis_resid = np.array([0], dtype=float)

    # Calculate light field
    f90.solve_rte_with_callbacks(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        ntheta, nphi,
        b, abs_func_N, source_func_N, source_expansion_N, bc_func_N, vsf_func_N,
        rad, irrad,
        num_scatters, fd_flag, lis_opts,
        lis_iter, lis_time, lis_resid
    )

    # End timer
    toc = time.time()
    date = datetime.now().ctime()
    git_commit = ru.get_git_commit_hash()
    compute_time = toc - tic

    # Extract values from arrays
    lis_iter = int(lis_iter)
    lis_time = float(lis_time)
    lis_resid = float(lis_resid)

    scalar_params = {
        'b': b,
        'ns': ns,
        'ntheta': nphi,
        'nphi': nphi,
        'nx': nx,
        'ny': ny,
        'nz': nz,
        'ntheta': ntheta,
        'nphi': nphi,
        'nomega': nomega,
        'zmax': zmax,
        'rope_spacing': rope_spacing,
        'num_scatters': num_scatters,
        'fd_flag': fd_flag,
        'lis_opts': lis_opts,
        'date': date,
        'git_commit': git_commit,
        'compute_time': compute_time,
        'lis_iter': lis_iter,
        'lis_time': lis_time,
        'lis_resid': lis_resid,
        'sol_expr': sol_expr_str,
        'abs_expr': abs_expr_str,
        'source_expr': source_expr_str,
        'bc_expr': bc_expr_str,
        'vsf_expr': vsf_expr_str,
        'param_dict': param_dict_str,
    }

    results = {
        'rad': rad,
        'irrad': irrad,
    }

    return scalar_params, results

@ru.run_decorator
def solve_rte_with_callbacks(ns, nz, ntheta, nphi, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict, num_scatters, fd_flag):
    lis_opts = '-i gmres -restart 100'

    return solve_rte_with_callbacks_full(ns, nz, ntheta, nphi, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict, num_scatters, fd_flag, lis_opts)


## Study Functions ##

@ru.study_decorator
def grid_study_compute(a_water, b, kelp_dist, ns_list, nz_list, na_list, lis_opts):
    """
    Do grid study with Cartesian product of given resolutions.
    """

    # One scatter before FD
    num_scatters = 0

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

@ru.study_decorator
def grid_study_compute_onespace(a_water, b, kelp_dist, ns_list, na_list, lis_opts):
    """
    Do grid study with Cartesian product of given resolutions.
    """

    # One scatter before FD
    num_scatters = 0

    fd_flag = True

    # Actual calling will be performed by decorator.
    # Functions to be called
    func_list = []
    # Arguments to be passed
    args_list = [] # tuples/lists
    kwargs_list = [] # dictionaries

    # Loop through grid
    # using nz = ns
    for ns, na in it.product(ns_list, na_list):
        # FD
        func_list.append(kelp_calculate)
        args_list.append((
            a_water, b,
            ns, ns, na,
            kelp_dist
        ))
        kwargs_list.append({
            'num_scatters': num_scatters,
            'fd_flag': True,
            'lis_opts': lis_opts
        })

        # No scattering
        func_list.append(kelp_calculate)
        args_list.append((
            a_water, b,
            ns, ns, na,
            kelp_dist
        ))
        kwargs_list.append({
            'num_scatters': 0,
            'fd_flag': False,
            'lis_opts': lis_opts
        })

    return func_list, args_list, kwargs_list

@ru.study_decorator
def asymptotics_study_compute(a_water_list, b_list, kelp_dist, fd_ns, fd_nz, fd_na, as_ns, as_nz, as_na, num_scatters_list, lis_opts):
    """
    For a grid of IOPs (`a_water` and `b` values), compute FD solution
    and compare to asymptotics solution for a range of `num_scatters`.
    """

    # One scatter before FD
    pre_fd_num_scatters = 0

    # Actual calling will be performed by decorator.
    # Functions to be called
    func_list = []
    # Arguments to be passed
    args_list = [] # tuples/lists
    kwargs_list = [] # dictionaries

    # Loop through IOP grid
    for a_water in a_water_list:
        for b in b_list:
            # FD solution
            func_list.append(kelp_calculate)
            args_list.append((
                a_water, b,
                fd_ns, fd_nz, fd_na,
                kelp_dist
            ))
            kwargs_list.append({
                'num_scatters': pre_fd_num_scatters,
                'fd_flag': True,
                'lis_opts': lis_opts
            })

            # Asymptotics solutions
            for num_scatters in num_scatters_list:
                func_list.append(kelp_calculate)
                args_list.append((
                    a_water, b,
                    as_ns, as_nz, as_na,
                    kelp_dist
                ))
                kwargs_list.append({
                    'num_scatters': int(num_scatters),
                    'fd_flag': False,
                    'lis_opts': lis_opts
                })

    return func_list, args_list, kwargs_list

@ru.study_decorator
def fd_verify_compute(ns_list, nz_list, ntheta_list, nphi_list, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict):
    """
    Given a list of resolutions in each dimension,
    loop over each list while holding all others at
    the highest resolution in order to test the convergence
    order of the FD algorithm.
    """

    # Sort and collect all dimensions
    dim_names = ('ns', 'nz', 'ntheta', 'nphi')
    dim_resolutions = map(
        sorted,
        (ns_list, nz_list, ntheta_list, nphi_list)
    )
    max_res_list = map(
        max,
        (ns_list, nz_list, ntheta_list, nphi_list)
    )
    dim_dict = dict(zip(dim_names, dim_resolutions))

    # One scatter before FD
    num_scatters = 0
    fd_flag = True

    # Arguments which do not change between runs
    const_args = (
        rope_spacing, zmax, b,
        sol_expr, abs_expr, source_expr, bc_expr, vsf_expr,
        param_dict, num_scatters, fd_flag
    )

    # Actual calling will be performed by decorator.
    # Functions to be called
    func_list = []
    # Arguments to be passed
    args_list = [] # iterables
    kwargs_list = [] # dictionaries

    # Run the largest grid once
    ns, nz, ntheta, nphi = max_res_list
    func_list.append(solve_rte_with_callbacks)
    args_list.append((ns, nz, ntheta, nphi, *const_args))
    kwargs_list.append({})

    # Loop over dimensions
    for dim_num, dim_name in enumerate(dim_names):
        # List of resolutions in the current dimension
        current_dim = dim_dict[dim_name]
        # Set all resolutions to their maximum values
        current_res_list = max_res_list

        # Loop over all smaller resolutions in the current dimension
        for res in current_dim[:-1]:
            current_res_list[dim_num] = res
            ns, nz, ntheta, nphi = current_res_list

            # Run this grid size
            func_list.append(solve_rte_with_callbacks)
            args_list.append((ns, nz, ntheta, nphi, *const_args))
            kwargs_list.append({})

    return func_list, args_list, kwargs_list
