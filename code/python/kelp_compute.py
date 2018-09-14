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
from IPython.display import display
import matplotlib
import numpy as np
import sympy as sp
import functools
# NOTE: Using boltons.funcutils.wraps in place
# of functools.wraps in order to correctly preserve
# signature of wrapper for the sake of f2py
# see https://hynek.me/articles/decorators/
import boltons.funcutils as fu

matplotlib.use('Agg')

# local
import run_utils as ru
import mms

## Kelp-specific funcs ##

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
    # If the true solution is not known, just pass 0.
    # However, the expression is evaluated and its results are stored
    # in the .nc file along with the numerical results.
    sol_expr_str = str(sol_expr)
    abs_expr_str = str(abs_expr)
    source_expr_str = str(source_expr)
    bc_expr_str = str(bc_expr)
    vsf_expr_str = str(vsf_expr)

    # Convert parameter values dictionary to string
    param_dict_str = json.dumps(param_dict)

    # Convert expressions to sympy functions
    sol_sym = mms.symify(sol_expr, *space, *angle, **param_dict)
    source_sym = mms.symify(source_expr, *space, *angle, **param_dict)
    abs_sym = mms.symify(abs_expr, *space, **param_dict)
    bc_sym = mms.symify(bc_expr, *angle, **param_dict)
    vsf_sym = mms.symify(vsf_expr, delta, **param_dict)

    # Convert sympy functions to numpy functions
    sol_func_N = mms.expr_to_num(sol_expr, *space, *angle)
    source_func_N = mms.expr_to_num(source_expr, *space, *angle)
    abs_func_N = mms.expr_to_num(abs_expr, *space)
    bc_func_N = mms.expr_to_num(bc_expr, *angle)
    vsf_func_N = mms.expr_to_num(vsf_expr, delta)

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
    true_rad = np.asfortranarray(np.zeros([nx, ny, nz, nomega]))
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

    # Calculate true light field
    grid = mms.gen_grid(ns, nz, ntheta, nphi, rope_spacing, zmax)
    true_rad = sol_func_N(*grid)

    # Calculate approximate light field
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
        'true_rad': true_rad,
    }

    return scalar_params, results

@ru.run_decorator
def solve_rte_with_callbacks(ns, nz, ntheta, nphi, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict, num_scatters, fd_flag):
    lis_opts = '-i gmres -restart 100 -maxiter 5000'

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
def verify_single_space_compute(ns_list, ntheta, nphi, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, num_scatters, fd_flag, param_dict):
    """
    Maintain constant ntheta, nphi while looping
    over spatial resolutions together
    (body diagonal in resolution-space)
    """

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

    # Loop over all smaller resolutions in the current dimension
    for ns in ns_list:
        nz = ns

        # Run this grid size
        print("Running grid ({:2d},{:2d},{:2d},{:2d})".format(ns, nz, ntheta, nphi))
        func_list.append(solve_rte_with_callbacks)
        args_list.append((ns, nz, ntheta, nphi, *const_args))
        kwargs_list.append({})

    return func_list, args_list, kwargs_list

@ru.study_decorator
def verify_compute(ns_list, nz_list, ntheta_list, nphi_list, rope_spacing, zmax, b, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, num_scatters, fd_flag, param_dict):
    """
    Given a list of resolutions in each dimension,
    loop over each list while holding all others at
    the highest resolution in order to test the convergence
    order of the FD algorithm.
    """

    # Sort and collect all dimensions
    dim_names = ('ns', 'nz', 'ntheta', 'nphi')
    dim_resolutions = list(map(
        sorted,
        (ns_list, nz_list, ntheta_list, nphi_list)
    ))
    max_res_list = list(map(
        max,
        (ns_list, nz_list, ntheta_list, nphi_list)
    ))
    dim_dict = dict(zip(dim_names, dim_resolutions))

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
    print("Running grid ({:2d},{:2d},{:2d},{:2d})".format(ns, nz, ntheta, nphi))
    func_list.append(solve_rte_with_callbacks)
    args_list.append((ns, nz, ntheta, nphi, *const_args))
    kwargs_list.append({})

    # Loop over dimensions
    for dim_num, dim_name in enumerate(dim_names):
        # List of resolutions in the current dimension
        current_dim = dim_dict[dim_name]
        # Set all resolutions to their maximum values
        # Use [:] to just copy values and not modify list.
        current_res_list = max_res_list[:]

        # Resolutions for current dimension except largest
        smaller_res_current = [
            res for res in current_dim
            if res != max_res_list[dim_num]
        ]

        # Loop over all smaller resolutions in the current dimension
        for res in smaller_res_current:
            current_res_list[dim_num] = res
            ns, nz, ntheta, nphi = current_res_list

            # Run this grid size
            print("Running grid ({:2d},{:2d},{:2d},{:2d})".format(ns, nz, ntheta, nphi))
            func_list.append(solve_rte_with_callbacks)
            args_list.append((ns, nz, ntheta, nphi, *const_args))
            kwargs_list.append({})

    return func_list, args_list, kwargs_list

@ru.study_decorator
def verify_asym_compute(b_list, num_scatters_list, ns, nz, ntheta, nphi, rope_spacing, zmax, sol_expr, abs_expr, source_expr, bc_expr, vsf_expr, param_dict):
    """
    Maintain constant grid,
    loop over:
    - b
    - num_scatters

    without fd
    """

    fd_flag = False

    const_kwargs = {
        'ns': ns,
        'nz': nz,
        'ntheta': ntheta,
        'nphi': nphi,
        'rope_spacing': rope_spacing,
        'zmax': zmax,
        'sol_expr': sol_expr,
        'abs_expr': abs_expr,
        'source_expr': source_expr,
        'bc_expr': bc_expr,
        'vsf_expr': vsf_expr,
        'fd_flag': fd_flag
    }

    # Actual calling will be performed by decorator.
    # Functions to be called
    func_list = []
    # Arguments to be passed
    args_list = [] # iterables
    kwargs_list = [] # dictionaries

    for b in b_list:
        param_dict['b'] = b
        for num_scatters in num_scatters_list:
            print("Running asym.: ({:.2f}, {:2d})".format(b, num_scatters))
            run_kwargs = {
                'b': b,
                'num_scatters': num_scatters,
                'param_dict': param_dict.copy(),
                **const_kwargs
            }

            func_list.append(solve_rte_with_callbacks)
            args_list.append([])
            kwargs_list.append(run_kwargs)

    return func_list, args_list, kwargs_list

@ru.study_decorator
def verify_asym_noscat_1d_compute(nz_list, rope_spacing, zmax, abs_expr, bc_expr, param_dict):
    ns = 1
    ntheta = 1
    nphi = 2
    num_scatters = 0

    b = 0
    param_dict_copy = param_dict.copy()
    param_dict_copy['b'] = b

    # Having trouble with normal integration
    # for some reason, so computing antiderivative
    # and evaluating manually
    abs_antideriv = sp.integrate(abs_expr, sp.Symbol('z'))
    sol_expr = bc_expr * sp.exp(
        abs_antideriv.subs('z', 0) - abs_antideriv
    )

    source_expr = 0
    vsf_expr = 0

    fd_flag = False

    const_kwargs = {
        'ns': ns,
        'ntheta': ntheta,
        'nphi': nphi,
        'rope_spacing': rope_spacing,
        'zmax': zmax,
        'sol_expr': sol_expr,
        'abs_expr': abs_expr,
        'source_expr': source_expr,
        'bc_expr': bc_expr,
        'vsf_expr': vsf_expr,
        'fd_flag': fd_flag,
        'b': b,
        'num_scatters': num_scatters,
        'param_dict': param_dict_copy
    }

    func_list = []
    args_list = []
    kwargs_list = []

    for nz in nz_list:
        run_kwargs = {'nz': nz, **const_kwargs}

        func_list.append(solve_rte_with_callbacks)
        args_list.append([])
        kwargs_list.append(run_kwargs)

    return func_list, args_list, kwargs_list

@ru.study_decorator
def verify_asym_noscat_const_abs_and_source_compute(ns, nz, ntheta, nphi, a, sigma, rope_spacing, zmax, param_dict):

    num_scatters = 0
    b = 0
    param_dict_copy = param_dict.copy()
    param_dict_copy['b'] = b

    abs_expr = a
    source_expr = sigma
    bc_expr = 0

    # Path length from origin
    s = sp.Piecewise(
        (
            sp.Symbol('z') / sp.cos(sp.Symbol('phi')),
            sp.Symbol('phi') < sp.pi/2
        ),
        (
            (sp.Symbol('z') - zmax) / sp.cos(sp.Symbol('phi')),
            True
        )
    )

    # Analytical solution
    sol_expr = sigma/a * (1 - sp.exp(-a*s))
    vsf_expr = 0

    fd_flag = False

    const_kwargs = {
        'ntheta': ntheta,
        'nphi': nphi,
        'rope_spacing': rope_spacing,
        'zmax': zmax,
        'sol_expr': sol_expr,
        'abs_expr': abs_expr,
        'source_expr': source_expr,
        'bc_expr': bc_expr,
        'vsf_expr': vsf_expr,
        'fd_flag': fd_flag,
        'b': b,
        'num_scatters': num_scatters,
        'param_dict': param_dict_copy
    }

    func_list = []
    args_list = []
    kwargs_list = []

    run_kwargs = {
        'ns': ns,
        'nz': nz,
        'ntheta': ntheta,
        'nphi': nphi,
        **const_kwargs
    }

    func_list.append(solve_rte_with_callbacks)
    args_list.append([])
    kwargs_list.append(run_kwargs)

    return func_list, args_list, kwargs_list

@ru.study_decorator
def verify_ss_asym_noscat_compute(ns_list, ntheta, nphi, abs_expr, source_expr, bc_expr, rope_spacing, zmax, param_dict):

    num_scatters = 0
    b = 0
    param_dict_copy = param_dict.copy()
    param_dict_copy['b'] = b

    vsf_expr = 0
    fd_flag = False

    abs_sym = mms.symify(abs_expr, *mms.space)
    source_sym = mms.symify(source_expr, *mms.space, *mms.angle)
    bc_sym = mms.symify(bc_expr, *mms.angle)
    vsf_sym = mms.symify(0, sp.Symbol('Delta'))

    s = sp.Symbol('s')
    s_p = sp.Symbol('s_p')
    s_pp = sp.Symbol('s_{pp}')

    x0, y0 = sp.var('x_0, y_0')
    z_hat = sp.Matrix([0,0,1])
    vec_x0 = sp.Matrix([
        x0,
        y0,
        sp.Piecewise(
            (0, mms.dot(mms.vec_om, z_hat) > 0),
            (zmax, True)
        )
    ])

    vec_l0_p = sp.simplify(mms.vec_l0(vec_x0, mms.vec_om, s_p, zmax))
    vec_l0_pp = sp.simplify(mms.vec_l0(vec_x0, mms.vec_om, s_pp, zmax))

    print("Evaluating abs & source")
    a_tilde = abs_sym(
        *vec_l0_pp
    )
    sigma_tilde = source_sym(
        *vec_l0_p,
        *mms.angle
    )

    print("a_tilde:")
    display(a_tilde)

    print("sigma_tilde:")
    display(sigma_tilde)

    print("inner")

    inner = sp.integrate(
        a_tilde,
        (s_pp, s_p, s)
    )
    display(inner)

    print("inner_prod")
    inner_prod = sp.simplify(
        sigma_tilde * sp.exp(-inner),
    )
    display(inner_prod)

    print("outer")
    outer = sp.integrate(
        inner_prod,
        (s_p, 0, s)
    )
    display(outer)

    u0_source_expr = outer

    print("Double integral")
    # Integrate light from distributed source

    print("BC")
    # Integrate light from boundary condition
    u0_bc_expr = (
        bc_sym(*mms.angle) * sp.exp(
            -sp.integrate(
                a_tilde,
                (s_pp, 0, s)
            )
        )
    )

    print("Combine")
    # Superpose source and bc solutions
    u0_s_expr = u0_source_expr + u0_bc_expr

    print("Lambdify u0")
    display(u0_s_expr)
    # Convert to funcion of x, y, z

    # Manually extract upwelling and downwelling
    # pieces because sympy can't do piecewise lambdify
    u0_s_down, u0_s_up = mms.split_piecewise(u0_s_expr)
    u0_down_func = sp.lambdify(
        ('s', 'x_0', 'y_0'),
        u0_s_down,
        modules=("sympy",)
    )
    u0_up_func = sp.lambdify(
        ('s', 'x_0', 'y_0'),
        u0_s_up,
        modules=("sympy",)
    )
    print("Evaluate L(x,y,z)")
    ph = sp.Symbol('phi')
    th = sp.Symbol('theta')
    s_tilde_sym = sp.Piecewise(
        (sp.Symbol('z')/sp.cos(ph), sp.cos(ph) > 0),
        ((sp.Symbol('z')-zmax)/sp.cos(ph), True)
    )
    s_down, s_up = mms.split_piecewise(s_tilde_sym)
    x0_sym = sp.simplify(sp.Symbol('x') - s_tilde_sym * sp.sin(ph)*sp.cos(th))
    y0_sym = sp.simplify(sp.Symbol('y') - s_tilde_sym * sp.sin(ph)*sp.sin(th))

    print("x0:")
    display(x0_sym)
    print("y0:")
    display(y0_sym)

    x0_down, x0_up = mms.split_piecewise(x0_sym)
    y0_down, y0_up = mms.split_piecewise(y0_sym)

    print("evaluate sol.")
    sol_down_expr = u0_down_func(s_down, x0_down, y0_down)
    sol_up_expr = u0_up_func(s_up, x0_up, y0_up)
    print("sol_down_expr:")
    display(sol_down_expr)
    print("sol_up_expr:")
    display(sol_up_expr)

    print("combined")
    sol_expr = sp.Piecewise(
        (sol_down_expr, sp.cos(sp.Symbol('phi')) > 0),
        (sol_up_expr, True)
    )

    display(sol_expr)

    print("symify")
    sol_down_sym = mms.symify(sol_down_expr, *mms.space, *mms.angle)
    sol_up_sym = mms.symify(sol_up_expr, *mms.space, *mms.angle)

    print("Checking down")
    down_diff = mms.check_sol(sol_down_sym, b, abs_sym, vsf_sym, source_sym)
    print("diff:")
    display(down_diff)
    num_down_diff = sp.lambdify(
        (*mms.space, *mms.angle),
        down_diff,
        modules=("numpy",)
    )
    max_down_diff = np.max(num_down_diff(*mms.gen_grid(10, 10, 10, 10, 1, 1)))
    print("max_down_diff = {:.2e}".format(max_down_diff))

    print("Checking up")
    up_diff = mms.check_sol(sol_up_sym, b, abs_sym, vsf_sym, source_sym)
    print("diff:")
    display(up_diff)
    num_up_diff = sp.lambdify(
        (*mms.space, *mms.angle),
        up_diff,
        modules=("numpy",)
    )
    max_up_diff = np.max(num_up_diff(*mms.gen_grid(10, 10, 10, 10, 1, 1)))
    print("max_up_diff = {:.2e}".format(max_up_diff))

    print("Done with sym. calc.")

    const_kwargs = {
        'ntheta': ntheta,
        'nphi': nphi,
        'rope_spacing': rope_spacing,
        'zmax': zmax,
        'sol_expr': sol_expr,
        'abs_expr': abs_expr,
        'source_expr': source_expr,
        'bc_expr': bc_expr,
        'vsf_expr': vsf_expr,
        'fd_flag': fd_flag,
        'b': b,
        'num_scatters': num_scatters,
        'param_dict': param_dict_copy
    }

    func_list = []
    args_list = []
    kwargs_list = []

    for ns in ns_list:
        nz = ns
        run_kwargs = {
            'ns': ns,
            'nz': nz,
            'ntheta': ntheta,
            'nphi': nphi,
            **const_kwargs
        }

        func_list.append(solve_rte_with_callbacks)
        args_list.append([])
        kwargs_list.append(run_kwargs)

    return func_list, args_list, kwargs_list
