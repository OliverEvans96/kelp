# 3rd-party
import numpy as np
import sympy as sp

# local
import mms

def lin_fit(x, y, x0, x1):
    x_arr = np.array(x)
    y_arr = np.array(y)
    which_inds = np.logical_and(
        x_arr>=x0,
        x_arr<=x1
    )
    x_fit = x_arr[which_inds]
    y_fit = y_arr[which_inds]

    def resid(args):
        m, b = args
        res = np.sum((m*x_fit + b - y_fit) ** 2)
        return res

    m0 = 1
    b0 = 0
    res = minimize(resid, (m0, b0))
    m, b = res.x

    return m, b

def plot_lin_fit(x, y, x0, x1, xlabel='x', ylabel='y'):
    xmin = np.min(x)
    ymin = np.min(y)
    xmax = np.max(x)
    ymax = np.max(y)

    plt.plot(x, y, 'o-')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.vlines((x0, x1), ymin, ymax, colors='k', linestyles='dashed')

    m, b = lin_fit(x, y, x0, x1)
    label = 'm={:.2f}, b={:.2f}'.format(m, b)
    plt.plot([xmin, xmax], [m*xmin + b, m*xmax + b], '--')
    plt.title(label)
    plt.show()

def max_derivs(expr, rope_spacing, zmax, do_space=True, do_angle=True, **param_vals):
    dims = []
    names = []

    grid = mms.gen_grid(10, 10, 10, 10, rope_spacing, zmax)
    inds = np.zeros_like(grid[0], dtype=bool)

    if do_space:
        dims += mms.space
        names += ['x', 'y', 'z']

        if do_angle:
            inds[...] = True
        else:
            inds[:,:,:,0] = True
            grid = grid[:3]
    if do_angle:
        dims += mms.angle
        names += ['theta', 'phi']

        if not do_space:
            inds[0,0,0,:]
            grid = grid[:-2]

    max_val_dict = {}
    for name, dim in zip(names, dims):
        deriv_expr = sp.diff(expr, dim)
        deriv_sym = mms.symify(deriv_expr, *dims, **param_vals)
        deriv_N = mms.sym_to_num(deriv_sym, *dims)

        max_val = np.max(np.abs(deriv_N(*[g[inds] for g in grid])))
        max_val_dict[name] = max_val

    return max_val_dict
