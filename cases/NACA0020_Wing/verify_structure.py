"""
Fluid-free check of the IB2d wing structure (the calibration NOTES.md asks for).

Uses IB2d's own force routines (IBM_Blackbox) and the wing's user-force hook,
so the EI -> k_theta mapping in Wing_Geom.bend_k is tested against the code that
will actually run:

  1. IB2d's stock beam force is NOT -grad(U) and carries a net force (why the
     wing uses the user-force bending instead); the user-force bending is
     checked against a finite-difference gradient and has zero net force
  2. static cantilever under a small uniform transverse load q on the flexible
     link (Newton solve on IB2d's forces), compared against
       a) Euler-Bernoulli with the same EI(x)              (continuum truth)
       b) Rayleigh-Ritz with the GVS cubic-curvature basis (what Wing_3order.jl solves)
  3. same load applied on the upper *skin* instead of the centerline (tests the ties)

Run:  python verify_structure.py            (needs numpy only)
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'external', 'IB2d', 'pyIB2d', 'IBM_Blackbox'))

import Wing_Geom as WG
import wing_case
from please_Find_Lagrangian_Forces_On_Eulerian_grid import (
    give_Me_Spring_Lagrangian_Force_Densities as spring_f,
    give_Me_Beam_Lagrangian_Force_Densities as ib2d_beam_f,
    give_Me_Target_Lagrangian_Force_Densities as target_f)
from give_Me_General_User_Defined_Force_Densities import give_Me_General_User_Defined_Force_Densities as bend_hook


def make_force(w):
    Lbig = 1e3    # periodic-wrap tests never trigger
    ids = w['targets'][:, 0].astype(int)
    tgt = np.column_stack([ids, w['xr'][ids], w['yr'][ids], w['targets'][:, 1]])   # clamp at the REST position
    Nb = len(w['x'])
    gen = w['beams'].T                     # columns L M R k_theta theta0 c_theta -> rows, as the driver loads them

    def force(x, y):
        fx, fy = spring_f(w['ds'], Nb, x, y, w['springs'], Lbig, Lbig)
        bx, by = bend_hook(w['ds'], Nb, x, y, x, y, 1.0, 0.0, gen)
        tx, ty = target_f(w['ds'], x, y, tgt, Lbig, Lbig)
        return fx + bx + tx, fy + by + ty
    return force


def check_beams(w):
    """(a) stock IB2d beam: net force on a kinked triple; (b) wing bending hook: gradient + net force."""
    x3 = np.array([0.0, 1.0, 2.0]); y3 = np.array([0.0, 0.1, -0.05])
    fx, fy = ib2d_beam_f(1.0, 3, x3, y3, np.array([[0, 1, 2, 1.0, 0.0]]), 1e3, 1e3)
    print('stock IB2d beam, kinked triple:   net force = (%.4f, %.4f)   [should be 0]' % (fx.sum(), fy.sum()))

    gen = np.array([[0], [1], [2], [3.0], [0.0], [0.0]])
    fx, fy = bend_hook(1.0, 3, x3, y3, x3, y3, 1.0, 0.0, gen)
    def U(x, y):
        a = (x[1] - x[0], y[1] - y[0]); b = (x[2] - x[1], y[2] - y[1])
        return 0.5 * 3.0 * np.arctan2(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1])**2
    num = np.zeros(6); v = np.concatenate([x3, y3]); h = 1e-7
    for i in range(6):
        vp = v.copy(); vp[i] += h; vm = v.copy(); vm[i] -= h
        num[i] = -(U(vp[:3], vp[3:]) - U(vm[:3], vm[3:])) / (2 * h)
    err = np.abs(np.concatenate([fx, fy]) - num).max()
    print('wing bending hook, same triple:   net force = (%.2e, %.2e),  max |f - (-grad U)| = %.2e' % (fx.sum(), fy.sum(), err))


def solve_static(w, ext_fx, ext_fy, nsteps=8):
    '''Newton solve of  f_struct(X) + f_ext = 0 on the free DOFs, load ramped in nsteps.'''
    force = make_force(w)
    x, y = w['xr'].copy(), w['yr'].copy()     # start from the straight rest state
    Nb = len(x)
    fixed = set(w['targets'][:, 0].astype(int))
    free = np.array([i for i in range(Nb) if i not in fixed])
    idx = np.concatenate([free, free + Nb])

    def resid(v, lam):
        xx, yy = v[:Nb], v[Nb:]
        fx, fy = force(xx, yy)
        return np.concatenate([fx + lam * ext_fx, fy + lam * ext_fy])[idx]

    v = np.concatenate([x, y])
    for step in range(1, nsteps + 1):
        lam = step / nsteps
        for it in range(40):
            r = resid(v, lam)
            if np.linalg.norm(r, np.inf) < 1e-7 * (1 + np.abs(ext_fy).max()):
                break
            J = np.zeros((len(idx), len(idx)))
            h = 1e-8
            for c, ci in enumerate(idx):
                vp = v.copy(); vp[ci] += h
                J[:, c] = (resid(vp, lam) - r) / h
            dv = np.linalg.lstsq(J, -r, rcond=None)[0]
            v = v.copy(); v[idx] += dv
    return v[:Nb], v[Nb:]


def reference_solutions(E, span, rigid_len, L3, floor, q, init=None):
    '''EB (variable EI) and GVS-style cubic-curvature Ritz tip deflection under uniform load q [N/m] on the link.'''
    s = np.linspace(0.0, L3, 4001)
    xp = rigid_len + s
    wth = np.maximum(WG.thickness(xp), floor)
    EI = E * span * wth**3 / 12.0

    def cumint(f):
        return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(s))])
    M = q * (L3 - s)**2 / 2.0
    y_eb = cumint(cumint(M / EI))

    P = [np.polynomial.legendre.Legendre.basis(k)(2 * s / L3 - 1) / L3 for k in range(4)]
    Y = [cumint(cumint(p)) for p in P]
    K = np.array([[np.trapezoid(EI * a * b, s) for b in P] for a in P])
    F = np.array([q * np.trapezoid(y, s) for y in Y])
    c = np.linalg.solve(K, F)
    y_gvs = sum(ci * yi for ci, yi in zip(c, Y))
    return s, y_eb, y_gvs


def main():
    p = wing_case.load_case_params()
    g = p['geom']
    E, span, rigid_len, floor = g['E_flex'], g['span'], g['rigid_length'], g['EI_floor_thickness']
    L3 = WG.CHORD - rigid_len

    w = WG.build_wing(bent=False)
    print('IB2d units: ds=%.5f  Nb=%d' % (w['ds'], len(w['x'])))
    check_beams(w)

    q = 0.5   # N/m per unit depth; small enough to stay linear (tip deflection << L3)
    iM, jr, nc, s_c = w['iM'], w['jr'], w['nc'], w['s_c']
    s_ref, y_eb, y_gvs = reference_solutions(E, span, rigid_len, L3, floor, q)

    def report(tag, y_chain):
        # chain tip deflection vs references
        print('  %-28s tip y = %.5f m    (EB %.5f, GVS-cubic %.5f)   chain/EB = %.3f   chain/GVS = %.3f'
              % (tag, y_chain, y_eb[-1], y_gvs[-1], y_chain / y_eb[-1], y_chain / y_gvs[-1]))

    # -- centerline loading (tests beams + axial springs + clamp) --
    print('Uniform load q = %.2f N/m on the flexible link, EI floor thickness = %.1f mm' % (q, 1e3 * floor))
    ext_fx = np.zeros(len(w['x'])); ext_fy = np.zeros(len(w['x']))
    for j in range(jr, nc + 1):
        wgt = 0.5 if j in (jr, nc) else 1.0
        ext_fy[iM[j]] += q * s_c * wgt / w['ds']          # F_phys / ds -> IB2d force density
    x, y = solve_static(w, ext_fx, ext_fy)
    report('load on centerline', y[iM[nc]] - w['y0'])

    # -- skin loading (tests the tie triangulation as well) --
    ext_fx[:] = 0; ext_fy[:] = 0
    xs = w['xs']
    flex = [k for k in range(1, w['n']) if xs[k] > rigid_len]
    dxk = np.gradient(xs)
    for k in flex:
        ext_fy[w['iU'][k]] += q * dxk[k] / w['ds']
    x2, y2 = solve_static(w, ext_fx, ext_fy)
    report('load on upper skin', y2[iM[nc]] - w['y0'])

    # -- profile comparison along the centerline --
    print('  centerline deflection profile [mm]:  s/L3   chain   EB    GVS-cubic')
    for j in range(jr, nc + 1, 2):
        sj = s_c * (j - jr)
        print('      %.2f   %7.3f  %7.3f  %7.3f' % (sj / L3, 1e3 * (y[iM[j]] - w['y0']),
              1e3 * np.interp(sj, s_ref, y_eb), 1e3 * np.interp(sj, s_ref, y_gvs)))

    # -- what the un-floored (true) EI would do at the tip --
    _, y_eb0, y_gvs0 = reference_solutions(E, span, rigid_len, L3, 1e-6, q)
    print('  (no EI floor:  EB tip %.5f m  [log-divergent], GVS-cubic tip %.5f m)' % (y_eb0[-1], y_gvs0[-1]))


if __name__ == '__main__':
    main()
