'''
The SoRoLily reference (WaterLily_FSI, WaterLily_2Way_Coupling_Example:
gvs_tests.jl + gvs_wings/Wing_3order.jl), ported to Python where the
comparison needs it: constants, the wing's kinematics and the fluid panels.

Everything here is transcribed from the Julia source (file:line given), so the
IB2d results can be compared with exactly what the reference computed:

  * kinematics (SoRoLily2D.jl extract_points_2d, exp_SE2): the carriage is
    frozen (gvs_functions.jl dynamics_gvs sets the heave/pitch accelerations
    to 0 and they start at 0), so link 2 (rigid, 0.06 m) starts at
    pivot - 0.03 m (Wing_3order.jl center_offset) and link 3 (0.14 m) is
    integrated from its base station by station with the midpoint rule:
    g <- g exp(h (kappa(X_mid), 1, 0)), kappa(X) = sum_k q_k P_k(2X - 1) / L3
    (Phi_LegendrePolynomial, soft_basis_at: curvature row scaled by 1/L3);
  * the fluid body (gvs_functions.jl geometry_gvs): a closed polygon (a
    degree-1 NURBS, body_simulations.jl simulation_initialisation) through
    every other station of 11 rigid and 21 soft stations, upper and lower
    surfaces at +-w(X)/2 along the local normal (sorolily_helper_functions.jl
    centerline_points); 31 points, 30 panels, starting and ending at the TE;
  * forces (gvs_functions.jl get_forces_gvs): F = -p n l per panel, scaled
    from simulation to physical units by rho U^2 (= 40 Pa per unit of p) and
    by the span w = 0.1 m, plus a viscous term. force_history rows (Julia,
    1-based) are [sum of panel torques about their own reference points,
    sum Fx, generalized forces (heave, pitch, q1..q4)]: total Fy is NOT
    stored and is rebuilt here from the stored (filtered) panel pressures.
'''

import numpy as np

# --- constants (gvs_tests.jl, Wing_3order.jl) --------------------------------
U = 0.2                     # m/s, U_p
CHORD = 0.2                 # m, L_p
NU = 1.0023e-6              # m^2/s, nu_p
RHO = 1000.0                # kg/m^3 (get_forces_gvs: 0.5*1000*U_p^2)
SPAN = 0.1                  # m, w passed to get_forces_gvs and h in Wing_3order.jl
L2, L3 = 0.06, 0.14         # m, link lengths (ld)
PIVOT = np.array([0.4, 0.6])  # m: offset [2*2^5, 3*2^5] cells / (L_c / L_p = 160 cells/m)
CELLS_PER_CHORD = 32        # L_c
DOMAIN = (1.2, 1.2)         # m: sz = (6*2^5, 6*2^5) cells
DT = 0.01                   # s, fixed fluid (and coupling) step
P_SCALE = RHO * U ** 2      # Pa per unit of simulation pressure (U_c = 1, rho_c = 1)
E, ETA, RHO_S, POISSON = 1e4, 2e3, 1064.6, 0.5   # link 3 material
N_RIGID, N_SOFT = 11, 21    # stations: 2*f_n+1, 2*t_n+1 with f_n = 5, t_n = 10
FORCE_CUTOFF_HZ = 1.5       # Butterworth order 6: pressure at cutoff = 3 Hz, forces at 0.5*cutoff
PRESSURE_CUTOFF_HZ = 3.0


def thickness(x):
    '''Full thickness [m] at x [m] from the LE (Wing_3order.jl
    thickness_from_physical_x, unclipped as in the reference: -4e-5 m at x = 0.2).'''
    f = 5.0 * np.asarray(x, dtype=float)
    return 0.4 * (0.0322 * f ** 3 - 0.0333 * np.cbrt(f) - 0.1269 * f ** 2
                  + 0.3797 * np.sqrt(f) - 0.2518 * f)


def curvature(q, X):
    '''kappa [1/m] of link 3 at X = s / L3 for the 4 curvature coordinates q.'''
    return np.polynomial.legendre.legval(2.0 * np.asarray(X) - 1.0, q) / L3


def _exp_se2(w, vx, vy):
    '''exp_SE2 (SoRoLily2D.jl): (R, p) of the twist (w, vx, vy).'''
    c, s = np.cos(w), np.sin(w)
    if abs(w) <= 1e-2:
        V = np.array([[1 - w * w / 6, -w / 2], [w / 2, 1 - w * w / 6]])
    else:
        V = np.array([[s / w, -(1 - c) / w], [(1 - c) / w, s / w]])
    return np.array([[c, -s], [s, c]]), V @ np.array([vx, vy])


def stations(q, n_rigid=N_RIGID, n_soft=N_SOFT):
    '''
    Poses along the wing for curvature coordinates q (4,), world frame [m]:
    returns (rigid, soft), each a dict with 'c' (n, 2) backbone points,
    'n' (n, 2) unit normals (local y axis), 'r' (n,) half-thickness, 'X'.
    Same station placement and midpoint-rule integration as extract_points_2d.
    '''
    out = []
    R = np.eye(2)
    p = PIVOT + np.array([-L2 / 2, 0.0])          # link 2 base: gi offset by -ld[2]/2
    for link, n, L, x0 in (('rigid', n_rigid, L2, 0.0), ('soft', n_soft, L3, L2)):
        X = np.linspace(0.0, 1.0, n)
        c, nrm = np.zeros((n, 2)), np.zeros((n, 2))
        for k in range(n):
            c[k], nrm[k] = p, R[:, 1]
            if k < n - 1:
                xm, h = 0.5 * (X[k] + X[k + 1]), (X[k + 1] - X[k]) * L
                kap = curvature(q, xm) if link == 'soft' else 0.0
                dR, dp = _exp_se2(h * kap, h, 0.0)
                p, R = p + R @ dp, R @ dR
        out.append(dict(c=c, n=nrm, r=0.5 * thickness(x0 + L * X), X=X))
    return out[0], out[1]


def trailing_edge(q):
    '''TE position [m] (the soft backbone's last station).'''
    return stations(q)[1]['c'][-1]


def straight_te():
    return PIVOT + np.array([-L2 / 2 + CHORD, 0.0])


def panels(q):
    '''
    The fluid body of geometry_gvs: polygon points (31, 2) [m] (TE -> upper
    surface -> LE -> lower surface -> TE), panel midpoints, lengths, outward
    normals ([0 1; -1 0] * tangent), and for each panel a label
    ('upper'/'lower') and the chordwise rest position x' [m] of its midpoint
    (for matching panels between codes).
    '''
    rig, sof = stations(q)
    up = lambda s: s['c'] + s['r'][:, None] * s['n']
    lo = lambda s: s['c'] - s['r'][:, None] * s['n']
    # Julia (1-based): [reverse(soft_upper[:,3:2:end]) reverse(rigid_upper[:,1:2:end])
    #                   rigid_lower[:,3:2:end] soft_lower[:,3:2:end-2] soft_upper[:,end]]
    su, ru, rl, sl = up(sof), up(rig), lo(rig), lo(sof)
    pts = np.vstack([su[2::2][::-1], ru[0::2][::-1], rl[2::2], sl[2:-1:2], su[-1:]])
    xs_s = L2 + L3 * sof['X']
    xs_r = L2 * rig['X']
    xr = np.concatenate([xs_s[2::2][::-1], xs_r[0::2][::-1], xs_r[2::2], xs_s[2:-1:2], xs_s[-1:]])
    d = np.diff(pts, axis=0)
    length = np.hypot(d[:, 0], d[:, 1])
    t = d / length[:, None]
    normal = np.column_stack([t[:, 1], -t[:, 0]])
    mid = pts[:-1] + 0.5 * d
    x_mid = 0.5 * (xr[1:] + xr[:-1])
    n_up = len(su[2::2]) + len(ru[0::2]) - 1      # panels from the TE to the LE
    side = np.array(['upper'] * n_up + ['lower'] * (len(mid) - n_up))
    return dict(points=pts, mid=mid, length=length, normal=normal, side=side, x_rest=x_mid)


def pressure_force(q, p_sim):
    '''Pressure force (Fx, Fy) [N, over the 0.1 m span] from the panel
    pressures p_sim (simulation units), as get_forces_gvs computes F_p.'''
    pn = panels(q)
    f = -(np.asarray(p_sim)[:, None] * pn['normal'] * pn['length'][:, None]).sum(axis=0)
    return f * P_SCALE * SPAN


def butterworth(cutoff_hz, fs=1.0 / DT, order=6, output='ba'):
    '''The reference's causal low-pass filter (DSP.jl
    digitalfilter(Lowpass(fc), Butterworth(6)); scipy.signal.butter is the same
    bilinear-transform design): (b, a), or second-order sections with
    output='sos' (needed at high sampling rates, where (b, a) is ill-conditioned).'''
    from scipy.signal import butter
    return butter(order, cutoff_hz, btype='low', fs=fs, output=output)
