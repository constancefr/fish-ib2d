"""
Fluid-free calibration of the triangulated wing (Wing_Geom.py): measures
`stiffness_scale`, the single multiplier on every spring, so that the IB2d
wing bends like the reference.

Why it is needed: tributary-area springs (mesh_utils.tributary_springs) put
each triangle's area into the springs exactly once, but a random spring
network with stiffness_scale = 1 is only ~1/3 as stiff as a sheet of modulus E
(each spring resists stretching along its own direction only). See NOTES.md,
section 3.

Static tests, small loads (linear), stiffness_scale = 1, so that
deflection is exactly proportional to 1 / stiffness_scale:

  1. strip: a straight, slender strip (0.14 x 0.007 m, L/t = 20, so shear is
     < 0.3 %) meshed the same way, clamped at one end, uniform load, compared
     with Euler-Bernoulli -> the "lattice" scale (what the spring network
     needs to behave like modulus E).
  2. wing: uniform transverse load q on the flexible link, centreline
     deflection compared with
       a) the reference's own model (Wing_3order.jl: Ritz with a cubic
          Legendre curvature basis, unshearable)  <- calibration target
       b) Euler-Bernoulli with the true thickness profile (for context)
     scale fitted by least squares over 0 < s/L3 <= FIT_UP_TO. The last 10 %
     is left out on purpose: there the trailing edge is < 2 mm thick and the
     reference's cubic curvature cannot follow it (model, not material).
  3. IB2d cross-check: IB2d's own spring + target force routines evaluated
     at the solved shape must balance the load (same units and forces as the
     real run).
  4. clamp: tip deflection with target_stiffness_factor x10.
  5. density sweep of 1 and 2 -> the default fill_density_multiplier is the
     lowest one whose scale is within 2 % of the next density's.

The static problem is solved with the linear stiffness matrix of the springs
and targets (k n n^T per spring), assembled sparse: thousands of points rule
out the finite-difference Newton solve the old centreline model used.

Run:  python verify_structure.py           (writes verify_structure.png)
"""

import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Wing_Geom as WG
import wing_case
import mesh_utils as mu
from ib_units import ib_ds

MULTIPLIERS = [1, 4, 8, 16, 32]
Q = 0.5             # N/m per unit depth: small enough to stay linear
FIT_UP_TO = 0.9     # fit stiffness_scale over 0 < s/L3 <= this
STRIP_L, STRIP_T = 0.14, 0.007


# --------------------------------------------------------------------------
# Linear static solve
# --------------------------------------------------------------------------

def stiffness_matrix(x, y, connections, edge_k, target_idx, k_target):
    '''Sparse linear stiffness (IB2d units) of the springs (k n n^T each) and
    the target tethers, at the rest configuration (x, y). DOF order x0, y0, x1, ...'''
    ij = np.array([(i, j) for i, j, _ in connections])
    k = np.array([edge_k[(i, j)] for i, j, _ in connections])
    d = np.column_stack([x[ij[:, 1]] - x[ij[:, 0]], y[ij[:, 1]] - y[ij[:, 0]]])
    n = d / np.linalg.norm(d, axis=1)[:, None]
    rows, cols, vals = [], [], []
    for a in range(2):
        for b in range(2):
            kab = k * n[:, a] * n[:, b]
            for p, q, s in ((0, 0, 1), (1, 1, 1), (0, 1, -1), (1, 0, -1)):
                rows.append(2 * ij[:, p] + a); cols.append(2 * ij[:, q] + b); vals.append(s * kab)
    t = np.asarray(target_idx)
    for a in range(2):
        rows.append(2 * t + a); cols.append(2 * t + a); vals.append(np.full(len(t), k_target))
    N = 2 * len(x)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))


def solve(x, y, connections, edge_k, target_idx, k_target, fy):
    '''Displacement (ux, uy) under the transverse force densities fy (IB2d units).'''
    K = stiffness_matrix(x, y, connections, edge_k, target_idx, k_target)
    f = np.zeros(2 * len(x)); f[1::2] = fy
    u = spla.spsolve(K.tocsc(), f)
    return u[0::2], u[1::2]


def surface_load(xs, ids, q, x_start, ds):
    '''Uniform load q [N/m] along x on the points `ids` (one surface, x-sorted
    or not) from x_start on: each point gets q times its tributary x-length,
    as an IB2d force density (F / ds).'''
    order = np.argsort(xs[ids])
    xo = np.concatenate([[x_start], xs[ids][order]])
    mid = 0.5 * (xo[1:] + xo[:-1])
    edges = np.concatenate([[x_start], mid[1:], [xo[-1]]])
    trib = np.diff(edges)
    fy = np.zeros(len(xs))
    fy[np.asarray(ids)[order]] = q * trib / ds
    return fy


def centreline_deflection(xr, yr, triangles, uy, xq, yq):
    '''uy interpolated (linear, on the rest triangulation) at the points (xq, yq).'''
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(xr, yr, triangles)
    return np.asarray(mtri.LinearTriInterpolator(tri, uy)(xq, yq))


def fit_scale(u1, u_ref):
    '''stiffness_scale s minimising sum (u1/s - u_ref)^2 (u1 = deflection at scale 1).'''
    return float(np.sum(u1 * u1) / np.sum(u1 * u_ref))


# --------------------------------------------------------------------------
# Reference solutions
# --------------------------------------------------------------------------

def reference_solutions(E, rigid_len, L3, q, n=40001):
    '''Deflection of the flexible link under a uniform load q [N/m], clamped
    at s = 0: Euler-Bernoulli with the true EI(s) = E w(s)^3 / 12, and the
    Ritz solution in the reference's cubic Legendre curvature basis (what
    Wing_3order.jl solves). The thickness fit crosses zero ~0.09 mm before
    the trailing edge (it is -4e-5 m at x' = 0.2), so EB needs a small
    thickness floor or that zero-stiffness sliver dominates the tip: with a
    1 um floor the tip is ~5.7 mm (the old NOTES' "4.9 mm, log-divergent"),
    with 10 um it converges to ~0.83 mm. The Ritz solution is insensitive.'''
    s = np.linspace(0.0, L3, n)
    w = np.maximum(WG.thickness(rigid_len + s), 1e-5)
    EI = E * w**3 / 12.0

    def cumint(f):
        return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(s))])
    y_eb = cumint(cumint(q * (L3 - s)**2 / 2.0 / EI))

    P = [np.polynomial.legendre.Legendre.basis(k)(2 * s / L3 - 1) / L3 for k in range(4)]
    Y = [cumint(cumint(p)) for p in P]
    K = np.array([[np.trapezoid(EI * a * b, s) for b in P] for a in P])
    F = np.array([q * np.trapezoid(yk, s) for yk in Y])
    y_gvs = sum(ck * yk for ck, yk in zip(np.linalg.solve(K, F), Y))
    return s, y_eb, y_gvs


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

def strip_test(ds, mult, E):
    '''Clamped slender strip, uniform load on its top edge -> fitted scale vs EB.'''
    L, t = STRIP_L, STRIP_T
    rect = np.array([[0.0, 0.0], [L, 0.0], [L, t], [0.0, t]])
    h = ds / np.sqrt(max(mult, 1.0))
    loop = mu.resample_loop(rect, h)
    fill = mu.fill_polygon_with_halton(rect, ds, density_multiplier=mult, min_boundary_gap=0.5 * h)
    pts, _, tris = mu.triangulate_with_holes(loop, [], loop, fill)
    x, y = pts[:, 0], pts[:, 1]
    conn, ek = mu.tributary_springs(x, y, tris, E, ds, 1.0)
    clamp = np.where(x < 1e-12)[0]
    kT = 15.0 * mu.max_incident_stiffness(clamp.tolist(), ek)
    top = np.where((np.abs(y - t) < 1e-12) & (x > 1e-12))[0]
    ux, uy = solve(x, y, conn, ek, clamp, kT, surface_load(x, top, Q, 0.0, ds))
    s = np.linspace(0.05, FIT_UP_TO, 18) * L
    u1 = centreline_deflection(x, y, tris, uy, s, np.full(len(s), t / 2))
    eb = Q * s**2 * (6 * L**2 - 4 * L * s + s**2) / (24 * E * t**3 / 12)
    return fit_scale(u1, eb), len(x)


def wing_test(mult, target_factor=15.0):
    '''Wing at stiffness_scale 1 under the uniform load on its flexible upper
    skin. Returns (w, s, u1): the structure, stations s along the flexible
    link and the centreline deflection there (m).'''
    w = WG.build_wing(fill_density_multiplier=mult, stiffness_scale=1.0,
                      target_stiffness_factor=target_factor, init_q=[0, 0, 0, 0])
    xp, yp = w['xr'] - w['x_LE'], w['yr'] - w['y0']
    rl = w['rigid_len']
    upper = np.where((np.arange(len(xp)) < w['n_outline']) & (yp >= 0) & (xp > rl))[0]
    fy = surface_load(xp, upper, Q, rl, w['ds'])
    ux, uy = solve(xp, yp, w['connections'], w['edge_k'], w['target_idx'], w['k_Target'], fy)
    L3 = WG.CHORD - rl
    s = np.linspace(0.0, L3, 57)
    u1 = centreline_deflection(xp, yp, w['triangles'], uy, rl + s[:-1], np.zeros(len(s) - 1))
    u1 = np.append(u1, uy[w['te_index']])              # the tip is the TE point itself
    return w, s, u1, (ux, uy, fy)


def ib2d_residual(w, ux, uy, fy, eps):
    '''IB2d's own spring + target forces at the solved shape, plus the load,
    with displacement and load both scaled by eps: max |residual| / max |load|.
    A linear mismatch (units, sign, a missing spring) would give the same
    ratio at any eps; the geometric nonlinearity the linear solve leaves out
    shrinks like eps. The wing is far stiffer in stretching than in bending,
    so at full load that second-order stretch already balances forces of the
    order of the load itself.'''
    ux, uy, fy = eps * ux, eps * uy, eps * fy
    sys.path.insert(0, os.path.join(HERE, '..', '..', 'external', 'IB2d', 'pyIB2d', 'IBM_Blackbox'))
    from please_Find_Lagrangian_Forces_On_Eulerian_grid import (
        give_Me_Spring_Lagrangian_Force_Densities as spring_f,
        give_Me_Target_Lagrangian_Force_Densities as target_f)
    xr, yr = w['xr'], w['yr']
    springs = np.array([(i, j, w['edge_k'][(i, j)], np.hypot(xr[j] - xr[i], yr[j] - yr[i]), 1.0)
                        for i, j, _ in w['connections']])
    t = w['target_idx']
    targets = np.column_stack([t, xr[t], yr[t], np.full(len(t), w['k_Target'])])
    x, y = xr + ux, yr + uy
    Lbig = 1e3                                  # no periodic wrap
    fx1, fy1 = spring_f(w['ds'], len(x), x, y, springs, Lbig, Lbig)
    fx2, fy2 = target_f(w['ds'], x, y, targets, Lbig, Lbig)
    r = np.concatenate([fx1 + fx2, fy1 + fy2 + fy])
    return np.abs(r).max() / np.abs(fy).max()


def main():
    p = wing_case.load_case_params()
    g = p['geom']
    E, rl = g['E_flex'], g['rigid_length']
    L3 = WG.CHORD - rl
    ds = ib_ds(wing_case.read_input2d('Lx'), int(wing_case.read_input2d('Nx')))
    s_ref, y_eb, y_gvs = reference_solutions(E, rl, L3, Q)

    print('Reference, uniform load q = %.2f N/m on the flexible link (deflection, mm):' % Q)
    print('   s/L3    beam (true thickness)   reference (cubic curvature)')
    for f in (0.5, 0.75, 0.9, 0.95, 1.0):
        print('   %.2f    %8.4f                %8.4f'
              % (f, 1e3 * np.interp(f * L3, s_ref, y_eb), 1e3 * np.interp(f * L3, s_ref, y_gvs)))

    print('\nDensity sweep (stiffness_scale fitted over s/L3 <= %.1f; tip = mesh tip / reference tip at that scale):'
          % FIT_UP_TO)
    print('   mult   strip pts  strip scale   wing pts  wing scale  vs beam   tip/ref  tip/beam   time[s]')
    mult0 = g.get('fill_density_multiplier', 8.0)
    rows = {}
    for mult in sorted(set(MULTIPLIERS) | {mult0}):
        t0 = time.perf_counter()
        a_strip, n_strip = strip_test(ds, mult, E)
        w, s, u1, _ = wing_test(mult)
        fit = (s > 0) & (s <= FIT_UP_TO * L3 + 1e-12)
        ref = np.interp(s, s_ref, y_gvs); eb = np.interp(s, s_ref, y_eb)
        a_ref, a_eb = fit_scale(u1[fit], ref[fit]), fit_scale(u1[fit], eb[fit])
        rows[mult] = (w, s, u1, a_ref)
        print('   %4g   %8d   %9.3f    %8d   %9.3f  %8.3f  %8.3f  %8.3f   %7.1f'
              % (mult, n_strip, a_strip, len(w['xr']), a_ref, a_eb,
                 u1[-1] / a_ref / ref[-1], u1[-1] / a_eb / eb[-1], time.perf_counter() - t0))

    w, s, u1, a = rows[mult0]
    print('\nAt the case density x%g:' % mult0)
    _, _, u1c, _ = wing_test(mult0, target_factor=150.0)
    print('   clamp: tip deflection with target_stiffness_factor 150 / 15 = %.4f' % (u1c[-1] / u1[-1]))
    w1, _, _, (ux, uy, fy) = wing_test(mult0)
    print('   IB2d force routines at the solved shape, max|residual| / max|load| at load x1, x1e-2, x1e-4:'
          ' %s  (falls ~eps: linear part exact)'
          % ', '.join('%.1e' % ib2d_residual(w1, ux, uy, fy, e) for e in (1.0, 1e-2, 1e-4)))
    ref = np.interp(s, s_ref, y_gvs)
    print('   profile at stiffness_scale = %.3f (mm):  s/L3   mesh   reference   beam' % a)
    for k in range(0, len(s), 7):
        print('        %.2f  %7.4f  %7.4f  %7.4f' % (s[k] / L3, 1e3 * u1[k] / a, 1e3 * ref[k],
                                                  1e3 * np.interp(s[k], s_ref, y_eb)))
    print('        %.2f  %7.4f  %7.4f  %7.4f   (tip)' % (1.0, 1e3 * u1[-1] / a, 1e3 * ref[-1], 1e3 * y_eb[-1]))
    print('   -> set geom.stiffness_scale = %.3f in case_params.json (currently %s)'
          % (a, g.get('stiffness_scale')))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    ax[0].plot(s_ref / L3, 1e3 * y_gvs, 'k-', label='reference (cubic curvature)')
    ax[0].plot(s_ref / L3, 1e3 * y_eb, 'k--', label='beam, true thickness')
    for mult, (_, ss, uu, aa) in rows.items():
        ax[0].plot(ss / L3, 1e3 * uu / aa, '.-', ms=3, lw=0.8, label='mesh x%g (scale %.2f)' % (mult, aa))
    ax[0].axvline(FIT_UP_TO, color='0.6', lw=0.8)
    ax[0].set_xlabel('s / L3'); ax[0].set_ylabel('deflection [mm]'); ax[0].legend(fontsize=7)
    ax[0].set_title('uniform load %.1f N/m, fitted scale' % Q)
    for mult, (_, ss, uu, aa) in rows.items():
        ref = np.interp(ss, s_ref, y_gvs)
        ax[1].plot(ss[1:] / L3, (uu[1:] / aa) / ref[1:] - 1, '.-', ms=3, lw=0.8, label='x%g' % mult)
    ax[1].axvline(FIT_UP_TO, color='0.6', lw=0.8); ax[1].axhline(0, color='k', lw=0.5)
    ax[1].set_xlabel('s / L3'); ax[1].set_ylabel('mesh / reference - 1'); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, 'verify_structure.png'), dpi=140)


if __name__ == '__main__':
    main()
