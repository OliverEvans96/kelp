# Grid Study
# Oliver Evans
# November 07, 2017

import numpy as np

def single_grid_case(ns, na, ratio, current, kelp):
    mat_time = None
    solve_time = None
    niter = None
    rad = None

    return mat_time, solve_time, niter, rad

def run_grid_study(ns_vals, na_vals, coef_ratios, current_profiles, kelp_shapes):
    for ns in ns_vals:
        for na in na_vals:
            for r, ratio in enumerate(coef_ratios):
                for c, current in enumerate(current_profiles):
                    for k, kelp in enumerate(kelp_shapes):
                        mat_time, solve_time, niter, rad = single_grid_case(ns, na, ratio, current, kelp)

if __name__ == '__main__':
    ns_vals = [8,12,16,20]
    na_vals = [8,12,16,20]
    coef_ratios = [1e-2,0,1e2]
    kelp_shapes = [
        lambda z, zmax, rmax: rmax*(1-z/zmax),
        lambda z, zmax, rmax: rmax*z/zmax,
        lambda z, zmax, rmax: rmax
    ]
    current_profiles = [
        {
            'speed': lambda z: 0,
            'angle': lambda z: 0
        },
        {
            'speed': lambda z: 1,
            'angle': lambda z: 0,
        },
        {
            'speed': lambda z: 1,
            'angle': lambda z, zmax: 2*np.pi*zmax
        }
    ]

    run_grid_study(ns_vals, na_vals, coef_ratios, current_profiles, kelp_shapes):
