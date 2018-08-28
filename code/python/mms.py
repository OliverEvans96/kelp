# Functions for the Method of Manufactured Solutions

import functools as ft
from boltons import funcutils as fu

import numpy as np
import sympy as sp
from IPython.display import display

import discrete_plot

from fortran_wrappers.light_utils_wrap import light_utils_wrap as light_utils

space = sp.var('x, y, z')
x, y, z = space
vec_x = sp.Matrix(space)

angle = sp.var(r'theta, phi')
th, ph = angle
vec_om = sp.Matrix([sp.sin(ph)*sp.cos(th), sp.sin(ph)*sp.sin(th), sp.cos(ph)])
delta = sp.var(r'Delta')

angle_prime = sp.var(r'theta_p, phi_p')
thp, php = angle_prime
vec_omp = sp.Matrix([sp.sin(php)*sp.cos(thp), sp.sin(php)*sp.sin(thp), sp.cos(php)])

def grad(expr, space=space):
    return sp.Matrix([sp.diff(expr, d) for d in space])

def sphere_integral(expr, angle):
    theta, phi = angle
    return sp.integrate(expr*sp.sin(phi), (theta, 0, 2*sp.pi), (phi, 0, sp.pi))

def dot(a, b):
    return sum(a.T * b)

def p_hat(l, m, na):
    nomega = na*(na-2)+2
    if m == 0:
        p = 0
    elif m == na:
        p = nomega
    else:
        p = (m-1)*na + l + 1

    return p

def vec_l(vec_x, vec_omega, s, zmax):
    z_hat = sp.Matrix([0,0,1])
    z0 = sp.Piecewise(
        (0, dot(vec_omega, z_hat) > 0),
        (zmax, True)
    )
    s_tilde = (dot(vec_x, z_hat) - z0) / dot(vec_omega, z_hat)
    vec_x0 = vec_x - s_tilde * vec_omega
    vec_l = (s*vec_x + (s_tilde-s)*vec_x0) / s_tilde
    return vec_l

def vec_l0(vec_x0, vec_omega, s, zmax):
    z_hat = sp.Matrix([0,0,1])
    z0 = sp.Piecewise(
        (0, dot(vec_omega, z_hat) > 0),
        (zmax, True)
    )
    s_tilde = (zmax - z0) / dot(vec_omega, z_hat)
    vec_x1 = vec_x0 + s_tilde * vec_omega
    vec_l = (s*vec_x1 + (s_tilde-s)*vec_x0) / s_tilde
    return vec_l

def gen_grid(ns, nz, na, rope_spacing, zmax):
    ds = rope_spacing/ns
    dz = zmax/nz

    x = y = -rope_spacing/2 + ds * (np.arange(ns) + 1/2)
    z = dz * (np.arange(nz) + 1/2)

    ntheta = nphi = na
    nomega = ntheta*(nphi-2) + 2

    dtheta = 2*np.pi/ntheta
    dphi = np.pi/(nphi-1)

    theta = dtheta * np.arange(ntheta)
    phi = dphi * np.arange(nphi)

    l = np.arange(ntheta)
    m = np.arange(nphi)
    p = np.arange(nomega)

    L, M = np.meshgrid(l, m, indexing='ij')

    phat = (M-1)*ntheta + L + 1
    phat[:,0] = 0
    phat[:,-1] = nomega-1

    theta_p = np.zeros(nomega)
    phi_p = np.zeros(nomega)
    theta_p[phat] = theta[L]
    theta_p[0] = theta_p[-1] = 0
    phi_p[phat] = phi[M]

    # A bit redundant, but seems to work
    X, Y, Z, Theta = np.meshgrid(x, y, z, theta_p, indexing='ij')
    _, _, _, Phi = np.meshgrid(x, y, z, phi_p, indexing='ij')

    return X, Y, Z, Theta, Phi

def display_eq(lhs_name, rhs_expr):
    """
    Use IPython.display.display to display a symbolic equation
    """
    display(sp.Eq(sp.Symbol(lhs_name), rhs_expr))

def subs_dict(expr, param_dict):
    """
    Substitute a dictionary of values into a symbolic expression.
    """
    # Multiply by a placeholder symbol to the expression in case it is originally
    # just a numeric constant, in which case .subs fails.
    placeholder = sp.Symbol('PLACEHOLDER')
    new_expr = expr * placeholder
    for key, val in param_dict.items():
        new_expr = new_expr.subs(key, val)

    # Set the placeholder symobol to one.
    new_expr = new_expr.subs(placeholder, 1)

    return new_expr

def symify(expr, *args, **subs):
    """
    Create symbolic sympy function from `expr`,
    with `args` as free variables,
    and `subs` values as fixed constants.
    """
    return sp.lambdify(
        args,
        subs_dict(expr, subs),
        modules=("sympy",)
    )

# ---

## Calculation functions

def calculate_bc(L, params=()):
    z = space[-1]
    zmin = 0
    return L(*space, *angle, *params).subs(z, zmin)

def calculate_source(L, b, a, beta, params=()):
    L_om = L(*space, *angle, *params)
    L_omp = L(*space, *angle_prime, *params)

    deriv = dot(vec_om, grad(L_om))
    atten = (a(*space)+b)*L_om

    scat_integrand = beta(dot(vec_om, vec_omp)) * L_omp
    scat = b * sphere_integral(scat_integrand, angle=angle_prime)

    source = deriv + atten - scat

    return source

def check_sol(L, b, a, beta, sigma):
    L_om = L(*space, *angle)
    L_omp = L(*space, *angle_prime)

    deriv = dot(vec_om, grad(L_om))
    atten = (a(*space)+b)*L_om

    scat_integrand = beta(dot(vec_om, vec_omp)) * L_omp
    scat = b * sphere_integral(scat_integrand, angle=angle_prime)

    source = sigma(*space, *angle)

    return sp.simplify(deriv + atten - scat - source)

# ---

## Series functions

def series_expr(expr, n):
    b = sp.Symbol('b')
    return sp.Poly(expr.series(b, n=n+1).removeO(), b).all_coeffs()[::-1]

def gen_series_N(expr, N, **param_vals):
    terms = [
        subs_dict(t, param_vals)
        for t in series_expr(expr, N)
    ]

    def lambda_N(x, y, z, theta, phi, n):
        args = (x, y, z, theta, phi)
        array_args = map(np.array, args)
        shape = np.shape(sum(array_args))
        ans = sp.lambdify(
            (*space, *angle),
            terms[n],
            modules=("numpy",)
        )(x, y, z, theta, phi)
        return np.broadcast_to(ans, shape)

    return lambda_N

def gen_series_sym(expr, N, **param_vals):
    terms = [
        subs_dict(t, param_vals)
        for t in series_expr(expr, N)
    ]

    def lambda_sym(x, y, z, theta, phi, n):
        return sp.lambdify(
            (*space, *angle),
            terms[n],
            modules=("sympy",)
        )(x, y, z, theta, phi)

    return lambda_sym