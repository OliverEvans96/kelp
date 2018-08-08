import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def maxind1d(bool_arr):
    """Return largest True index"""
    return max(np.where(bool_arr)[0])

## Piecewise-constant interpolation

def merge_diff_grids(xmin, xmax, y1, y2):
    n1 = len(y1)
    n2 = len(y2)
    d1 = (xmax-xmin)/n1
    d2 = (xmax-xmin)/n2

    # Edges
    e1 = np.linspace(xmin, xmax, n1+1)
    e2 = np.linspace(xmin, xmax, n2+1)

    # Merge e1 and e2
    i_sort = np.argsort(np.concatenate([e1, e2]))
    i1 = i_sort<n1+1
    i2 = i_sort>=n1+1
    e3 = np.zeros(n1+n2+2)
    e3[i1] = e1
    e3[i2] = e2
    # Remove duplicate extreme endpoints (xmin, xmax)
    e3 = e3[1:-1]

    # Heights of new edges
    z1 = np.zeros(n1+n2-1)
    z2 = np.zeros(n1+n2-1)

    # Loop over edges, excluding xmax
    for j in range(n1+n2-1):
        k1 = maxind1d(e1<=e3[j])
        k2 = maxind1d(e2<=e3[j])
        z1[j] = y1[k1]
        z2[j] = y2[k2]

    return e3, z1, z2

def plot_discrete_diff(xmin, xmax, y1, y2):
    plt.figure(figsize=[10,8])

    n1 = len(y1)
    n2 = len(y2)

    e1 = np.linspace(xmin, xmax, n1+1)
    e2 = np.linspace(xmin, xmax, n2+1)

    x1 = e1[:-1] + np.diff(e1)/2
    x2 = e2[:-1] + np.diff(e2)/2

    e3, z1, z2 = merge_diff_grids(xmin, xmax, y1, y2)

    # Difference
    plt.bar(
        e3[:-1], np.abs(z2-z1), np.diff(e3),
        bottom=np.min([z2,z1], axis=0),
        align='edge', ec=(0,0,1),
        fc=(0,0,1,0.5), lw=2, label='diff'
    )


    # Original
    plt.bar(
        e1[:-1], y1, np.diff(e1),
        align='edge',
        fc=(0,0,0,0), ec='C1', lw=2, label='y1'
    )
    plt.bar(
        e2[:-1], y2, np.diff(e2),
        align='edge',
        fc=(0,0,0,0), ec='C2', lw=2, label='y2'
    )

    # Points
    plt.plot(x1, y1, 'C1o')
    plt.plot(x2, y2, 'C2o')

    # Combined (1st attempt)
    # plt.plot(x3, z1, 'C3o-', label='z1')
    # plt.plot(x3, z2, 'C4o-', label='z2')

    plt.legend()

def get_bin_centers(xmin, xmax, n1, n2):
    """n1, n2 are number of bins (n1+1, n2+1 edges)"""
    # Edges
    e1 = np.linspace(xmin, xmax, n1+1)
    e2 = np.linspace(xmin, xmax, n2+1)

    # Bin centers
    x1 = e1[:-1] + np.diff(e1)/2
    x2 = e2[:-1] + np.diff(e2)/2

    return x1, x2

def discrete_err(xmin, xmax, y1, y2, verbose=False, net=False):
    """
    Positive area between two piecewise constant curves
    of different uniform partitions of [xmin, xmax].
    Relative error is divided by norm of y1.

    If net=True, calculate net area between curves.
    Otherwise, use absolute area between curves.
    """
    e3, z1, z2 = merge_diff_grids(xmin, xmax, y1, y2)
    diff = z1-z2 if net else np.abs(z1-z2)
    interval_widths = np.diff(e3)
    if(verbose):
        print("diff, widths:\n", '\n'.join(
            map(
                lambda tup: '({:.2e}, {:.2e})'.format(*tup),
                zip(diff, interval_widths)
            )
        ), sep='')
    tot_abs_err = abs(np.sum(diff * interval_widths))
    tot_rel_err = abs_err / np.sum(z1*np.diff(e3))

    avg_abs_err, avg_rel_err = map(
        lambda tot_err: tot_err / np.sum(interval_widths),
        (tot_abs_err, tot_rel_err)
    )

    return abs_err, rel_err

## Linear interpolation

def merge_linear(x1, x2, y1, y2):
    """
    Linearly interpolate y1 and y2 onto
    the union of x1 and x2.
    """

    # Use interp1d instead of griddata to enable 'extrapolate'
    x3 = np.sort([*x1, *x2])
    f1 = np.vectorize(interp1d(x1, y1, kind='linear', fill_value='extrapolate'))
    f2 = np.vectorize(interp1d(x2, y2, kind='linear', fill_value='extrapolate'))
    z1 = f1(x3)
    z2 = f2(x3)
    return x3, z1, z2

def err_linear(xmin, xmax, x1, x2, y1, y2, net=False):
    """
    If net=True, calculate net area between curves.
    Otherwise, use absolute area between curves.

    absolute error: average difference between curves
    relative error: average difference between curves
                    divided by average height of first curve

    """
    x3, z1, z2 = merge_linear(x1, x2, y1, y2)

    diff = z1-z2 if net else np.abs(z1-z2)
    avg_abs_err = abs(np.trapz(
        x=x3,
        y=diff
    ) / (xmax - xmin))

    # Global version:
    avg_rel_err = abs(np.trapz(
        x=x3,
        y=diff
    ) / np.trapz(
        x=x1,
        y=y1
    ))

    # # Pointwise version:
    # avg_rel_err = abs(np.trapz(
    #     x=x3,
    #     y=(diff / z1)
    # ))


    return avg_abs_err, avg_rel_err
