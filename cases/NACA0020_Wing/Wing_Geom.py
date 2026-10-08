"""
Builds the Lagrangian structure of the NACA0020 wing from the WaterLily_FSI
example (gvs_wings/Wing_3order.jl) and writes the IB2d input files
(.vertex, .spring, .user_force, .target, .geo_connect).

Structural model (see NOTES.md for the full GVS -> IB2d mapping):

  * Reference model: link 2 (0 - 0.06 m from the leading edge) is rigid and,
    because heave/pitch are frozen in the coupled run, fixed in space; link 3
    (0.06 - 0.20 m) is an inextensible, unshearable Euler-Bernoulli beam
    (curvature-only strain, cubic Legendre basis) of variable thickness
    w(x), E = 1e4 Pa (its Kelvin-Voigt damping eta = 2e3 Pa s is NOT
    reproduced, see NOTES.md).
  * IB2d: bending lives on a coarse *centerline* chain of vertices with a
    discrete Euler-Bernoulli energy 0.5 (EI/s) theta^2 per vertex, EI(x) =
    E * span * w(x)^3 / 12. IB2d's stock `.beam` is NOT used because its force
    is not the gradient of its energy (net force != 0, see
    give_Me_General_User_Defined_Force_Densities.py), so the bending (and the
    bending goes through `user_force_model` / wing.user_force. Inextensibility comes from axial springs and the clamp
    from stiff target points.
  * The fluid must see the real airfoil, so the outline is a *skin* of
    vertices (upper + lower surface, spacing ds) that is rigidly triangulated
    onto the centerline by two stiff "tie" springs per skin vertex. The skin
    carries no bending stiffness of its own (only a very weak outline spring
    that keeps the vertex spacing even and draws the outline in VisIt).
  * Everything is per unit depth (span = 1): the coupled WaterLily model scales
    structure (h = 0.1 m) and fluid load (w = 0.1) by the same span, so the
    response is span-independent.

IB2d unit convention (important): the driver spreads `f * ds` with the
*global* ds = Lx / (2 Nx), so a physical nodal force F [N per unit depth]
needs f = F / ds. Every stiffness below is therefore computed in physical
units and divided by ds via to_ib().
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wing_case

CHORD = 0.2   # m; the thickness fit below is only valid for this chord


# --------------------------------------------------------------------------
# Airfoil thickness + material -> IB2d constants
# --------------------------------------------------------------------------

def thickness(xp):
    '''Full airfoil thickness [m] at distance xp [m] from the leading edge.
    Verbatim from Wing_3order.jl (thickness_from_physical_x); zero at the tip.'''
    five_x = 5.0 * np.asarray(xp, dtype=float)
    w = 0.4 * (0.0322 * five_x**3 - 0.0333 * np.cbrt(five_x) - 0.1269 * five_x**2
               + 0.3797 * np.sqrt(five_x) - 0.2518 * five_x)
    return np.maximum(w, 0.0)


def to_ib(K_phys, ds):
    '''Physical stiffness [N/m per unit depth] -> IB2d spring/target constant.'''
    return K_phys / ds


def bend_k(EI, s, ds):
    """
    k_theta stored in wing.user_force for bending stiffness EI [N m per unit
    depth] on a chain of spacing s. Physical discrete energy per vertex is
    0.5 (EI/s) theta^2; the driver spreads f*ds, so k_theta = EI / (s * ds).
    """
    return EI / (s * ds)


def surface_stations(ds_skin):
    '''x' [m] of skin stations spaced ~ds_skin in *arc length* along the
    surface (the leading edge is nearly vertical, so uniform x' would leave
    huge gaps there). Both surfaces share the same stations.'''
    u = np.linspace(0.0, 1.0, 200001)
    xp = CHORD * u**2
    hw = 0.5 * thickness(xp)
    S = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xp), np.diff(hw)))])
    n = int(round(S[-1] / ds_skin))
    return np.interp(np.linspace(0.0, S[-1], n + 1), S, xp)


def centerline_shape(init_q, L3, s_nodes):
    '''
    Bent centerline of the flexible link for a curvature field
    kappa(X) = (1/L3) * sum_k q_k P_k(2X-1), X = s/L3 (same Legendre basis and
    1/L scaling as the GVS model), sampled at arc lengths s_nodes.
    Returns (x, y, theta) relative to the clamp, tangent = +x at s = 0.
    '''
    fine = np.linspace(0.0, L3, 20001)
    kappa = np.polynomial.legendre.legval(2.0 * fine / L3 - 1.0, init_q) / L3
    theta = np.concatenate([[0.0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * np.diff(fine))])
    xf = np.concatenate([[0.0], np.cumsum(0.5 * (np.cos(theta[1:]) + np.cos(theta[:-1])) * np.diff(fine))])
    yf = np.concatenate([[0.0], np.cumsum(0.5 * (np.sin(theta[1:]) + np.sin(theta[:-1])) * np.diff(fine))])
    return np.interp(s_nodes, fine, xf), np.interp(s_nodes, fine, yf), np.interp(s_nodes, fine, theta)


# --------------------------------------------------------------------------
# Structure builder (importable by verify_structure.py)
# --------------------------------------------------------------------------

def build_wing(params=None, Nx=None, Lx=None, bent=True):
    '''
    Returns a dict with the vertex coordinates (bent initial state if
    bent=True, straight rest state otherwise) plus the spring / beam / target
    tables in IB2d units, and index bookkeeping used for verification.
    '''
    g = (params or wing_case.load_case_params()).get('geom', {})
    Nx = Nx if Nx is not None else int(wing_case.read_input2d('Nx'))
    Lx = Lx if Lx is not None else wing_case.read_input2d('Lx')
    ds = Lx / (2.0 * Nx)          # IB2d's global Lagrangian spacing (hard-coded in the driver)

    E = g.get('E_flex', 1.0e4)
    span = g.get('span', 1.0)
    rigid_len = g.get('rigid_length', 0.06)
    s_c = g.get('centerline_spacing', 0.01)
    w_floor = g.get('EI_floor_thickness', 0.008)
    init_q = list(g.get('init_q', [0, 0, 0, 0]))
    K_tie, K_skin, K_target = g.get('K_tie', 1e4), g.get('K_skin', 50.0), g.get('K_target', 1e4)
    px, py = g.get('pitch_axis_world', [0.4, 0.6])
    x_LE = px - g.get('pitch_axis_from_LE', 0.03)
    y0 = py

    nc = int(round(CHORD / s_c))
    jr = int(round(rigid_len / s_c))
    assert abs(nc * s_c - CHORD) < 1e-9 and abs(jr * s_c - rigid_len) < 1e-9, \
        'centerline_spacing must divide both the chord and rigid_length'
    L3 = CHORD - rigid_len

    # ---- rest (straight) coordinates in the wing frame (x' from LE, y up) ----
    xs = surface_stations(ds)                 # k = 0 (LE) ... n (TE)
    n = len(xs) - 1
    hw = 0.5 * thickness(xs)

    # vertex order: centerline M_0..M_nc  (M_0 = LE point, M_nc = TE point),
    # then upper surface k=1..n-1, then lower surface k=1..n-1
    xM = s_c * np.arange(nc + 1)
    iM = np.arange(nc + 1)
    iU = {k: nc + k for k in range(1, n)}
    iL = {k: nc + (n - 1) + k for k in range(1, n)}
    Nb = nc + 1 + 2 * (n - 1)

    xr = np.zeros(Nb); yr = np.zeros(Nb)
    xr[iM] = xM
    for k in range(1, n):
        xr[iU[k]], yr[iU[k]] = xs[k], +hw[k]
        xr[iL[k]], yr[iL[k]] = xs[k], -hw[k]

    # which centerline segment each skin vertex hangs off
    seg = {k: int(np.clip(np.floor(xs[k] / s_c + 1e-12), jr, nc - 1)) for k in range(1, n)}
    fixed_skin = [k for k in range(1, n) if xs[k] <= rigid_len + 1e-12]
    flex_skin = [k for k in range(1, n) if xs[k] > rigid_len + 1e-12]

    # ---- bent initial coordinates ----
    x, y = xr.copy(), yr.copy()
    if bent and any(abs(q) > 0 for q in init_q):
        jf = np.arange(jr, nc + 1)
        cx, cy, cth = centerline_shape(init_q, L3, s_c * (jf - jr))
        x[iM[jf]] = rigid_len + cx
        y[iM[jf]] = cy
        for k in flex_skin:
            j = seg[k]
            b = np.array([x[iM[j + 1]] - x[iM[j]], y[iM[j + 1]] - y[iM[j]]])
            t_hat = b / np.hypot(*b)
            n_hat = np.array([-t_hat[1], t_hat[0]])
            along = xs[k] - xM[j]
            for idx, off in ((iU[k], +hw[k]), (iL[k], -hw[k])):
                x[idx] = x[iM[j]] + along * t_hat[0] + off * n_hat[0]
                y[idx] = y[iM[j]] + along * t_hat[1] + off * n_hat[1]

    # ---- springs: (i, j, k_ib, rest_length, 1.0) ----
    springs = []
    def add_spring(i, j, K_phys):
        rl = float(np.hypot(xr[j] - xr[i], yr[j] - yr[i]))
        springs.append((i, j, to_ib(K_phys, ds), rl, 1.0))

    # axial (inextensibility) springs along the flexible centerline
    for j in range(jr, nc):
        w_mid = max(float(thickness(0.5 * (xM[j] + xM[j + 1]))), w_floor)
        add_spring(iM[j], iM[j + 1], E * w_mid * span / s_c)
    # ties: skin vertex -> both ends of its centerline segment (rigid triangle)
    for k in flex_skin:
        for idx in (iU[k], iL[k]):
            add_spring(idx, iM[seg[k]], K_tie)
            add_spring(idx, iM[seg[k] + 1], K_tie)
    # weak outline springs (LE -> upper -> TE -> lower -> LE), closed loop
    outline = [iM[0]] + [iU[k] for k in range(1, n)] + [iM[nc]] + [iL[k] for k in range(n - 1, 0, -1)]
    for a, b in zip(outline, outline[1:] + outline[:1]):
        add_spring(a, b, K_skin)
    springs = np.array(springs)

    # ---- bending rows for wing.user_force: (L, M, R, k_theta, theta0) ----
    # one per flexible centerline vertex (the first sits on the clamp, so its
    # left neighbour is a fixed vertex -> fixes the root slope)
    beams = []
    for j in range(jr, nc):
        # compliance-averaged (harmonic-mean) EI over the vertex's dual cell; the
        # clamp vertex only owns the flexible half of its cell. Point-sampling EI
        # instead over-softens the chain by ~8% on this tapering section.
        xa = max(xM[j] - 0.5 * s_c, rigid_len); xb = min(xM[j] + 0.5 * s_c, CHORD)
        w_cell = np.maximum(thickness(np.linspace(xa, xb, 41)), w_floor)
        I = span / (12.0 * np.mean(1.0 / w_cell**3))
        L_, M_, R_ = iM[j - 1], iM[j], iM[j + 1]
        a = np.array([xr[M_] - xr[L_], yr[M_] - yr[L_]]); b = np.array([xr[R_] - xr[M_], yr[R_] - yr[M_]])
        th0 = float(np.arctan2(a[0] * b[1] - a[1] * b[0], a @ b))
        # the clamp vertex's rotation stands for the curvature integral over only HALF a
        # cell (theta_0 = kappa*s/2), so matching 0.5 EI kappa^2 (s/2) needs k = 2 EI/s;
        # this removes the O(s) excess compliance of the chain.
        root = 2.0 if j == jr else 1.0
        beams.append((L_, M_, R_, root * bend_k(E * I, s_c, ds), th0))
    beams = np.array(beams)

    # ---- targets: clamp everything in the rigid leading-edge segment ----
    target_ids = [int(iM[j]) for j in range(0, jr + 1)]
    for k in fixed_skin:
        target_ids += [iU[k], iL[k]]
    targets = np.array([(i, to_ib(K_target, ds)) for i in target_ids])

    return dict(x=x + x_LE, y=y + y0, xr=xr + x_LE, yr=yr + y0, springs=springs, beams=beams,
                targets=targets, outline=outline, ds=ds, s_c=s_c, iM=iM, jr=jr, nc=nc,
                x_LE=x_LE, y0=y0, iU=iU, iL=iL, n=n, xs=xs, seg=seg)


# --------------------------------------------------------------------------
# File writers (same formats as FinRay_Geom.py)
# --------------------------------------------------------------------------

def write_files(w, struct_name='wing'):
    with open(struct_name + '.vertex', 'w') as f:
        f.write('%d\n' % len(w['x']))
        for xi, yi in zip(w['x'], w['y']):
            f.write('%1.16e %1.16e\n' % (xi, yi))

    with open(struct_name + '.spring', 'w') as f:
        f.write('%d\n' % len(w['springs']))
        for i, j, k, rl, a in w['springs']:
            f.write('%d %d %1.16e %1.16e %1.16e\n' % (i, j, k, rl, a))

    with open(struct_name + '.user_force', 'w') as f:
        f.write('%d\n' % len(w['beams']))
        for L_, M_, R_, k, th0 in w['beams']:
            f.write('%d %d %d %1.16e %1.16e\n' % (L_, M_, R_, k, th0))

    with open(struct_name + '.target', 'w') as f:
        f.write('%d\n' % len(w['targets']))
        for i, k in w['targets']:
            f.write('%d %1.16e\n' % (i, k))

    with open(struct_name + '.geo_connect', 'w') as f:
        for i, j, *_ in w['springs']:
            f.write('%d %d\n%d %d\n' % (i, j, j, i))


def plot_preview(w, path='geometry_preview.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 4))
    for a, xx, yy, title in ((ax[0], w['xr'], w['yr'], 'rest (straight)'), (ax[1], w['x'], w['y'], 'initial')):
        for i, j, *_ in w['springs']:
            a.plot([xx[int(i)], xx[int(j)]], [yy[int(i)], yy[int(j)]], 'c-', lw=0.3)
        oi = w['outline'] + w['outline'][:1]
        a.plot(xx[oi], yy[oi], 'b-', lw=1.0)
        a.plot(xx[w['iM']], yy[w['iM']], 'r.-', ms=4, lw=0.8)
        tid = w['targets'][:, 0].astype(int)
        a.plot(xx[tid], yy[tid], 'k.', ms=3)
        a.set_aspect('equal'); a.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=140)


def main():
    w = build_wing()
    write_files(w)
    print('Nb=%d vertices, %d springs, %d bending vertices, %d target points (ds=%.5f)'
          % (len(w['x']), len(w['springs']), len(w['beams']), len(w['targets']), w['ds']))
    if '--plot' in sys.argv:
        plot_preview(w)


if __name__ == '__main__':
    main()
