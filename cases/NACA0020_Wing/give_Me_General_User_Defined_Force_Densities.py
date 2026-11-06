'''
IB2d `user_force_model = 1` hook for the NACA0020 wing (called by the driver
every step with the half-step Lagrangian positions). Two uses, chosen by the
contents of wing.user_force (written by the geometry script):

1. Current model (Wing_Geom.py writes one number: the mode).
   Mode 1 (geom.spring_force = "hook", input2d springs = 0): the springs of
   wing.spring, computed here vectorised with exactly IB2d's formula
   (give_Me_Spring_Lagrangian_Force_Densities: f = 0.5 (a+1) k (|d| - RL)^a d/|d|
   on the leader, minus that on the follower, same periodic-wrap rule).
   IB2d's version is a Python loop over springs (~17 us per spring per step
   on the cluster's nodes: half the step time at 192^2); the result is the
   same up to the order of floating-point summation (spring_force_check.py).
   Mode 0 (geom.spring_force = "ib2d", springs = 1): IB2d computes the
   springs; this hook adds no force.
   Both modes record every Lagrangian point and the hydrodynamic force on the wing
   every diagnostics.interval seconds into viz_IB2d/wing_history.npz (see
   wing_diagnostics.py):
       pos (n, Nb, 2) float32 [m]   half-step positions at t = step*dt + dt/2
       force (n, 2) [N/m]           force of the fluid on the wing per unit
                                    depth, from the clamp tethers:
                                    F = sum_i k_Target,i * ds * (x_i - x_target,i)
       force_every_step (n_steps+1, 2)   the same at every step (step k at
                                    t = k dt + dt/2), so it can be low-pass
                                    filtered without aliasing the clamp's
                                    ~2 ms ringing
   The wing is held only by the target tethers and the IB structure is
   massless, so F is the whole hydrodynamic load (pressure + viscous) on the
   wing, rigid block included.

2. Legacy bending force (Wing_Geom_centerline.py; wing.user_force rows
   L M R k_theta theta0), kept so the centreline model's runs can be
   reproduced. The rest of this docstring describes it.

Why this exists instead of IB2d's stock `.beam`: IB2d's beam force
(give_Me_Beam_Lagrangian_Force_Densities, python and matlab alike) is not the
gradient of its energy 0.5 k (cross - C)^2 -- the left/right-vertex forces have
the wrong sign, so a triple carries a net force (and a uniform-curvature arc,
which is a pure bending moment, feels a spurious force on every interior
vertex). That makes k_beam impossible to calibrate against a physical EI, so
bending is implemented here as a conservative, momentum- and moment-free
turning-angle energy per vertex triple (L, M, R):

    U_phys = 0.5 * (EI / s) * (theta - theta0)^2,   theta = angle(R-M) - angle(M-L)

which is the discrete Euler-Bernoulli energy 0.5 * EI * kappa^2 * s with
kappa = theta / s. The driver spreads `f * ds` (ds = Lx / 2Nx), so the constant
stored in wing.user_force is k_theta = EI / (s * ds) and the force density is

    f = -k_theta * (theta - theta0) * d(theta)/dX.

The reference's Kelvin-Voigt damping (Eta = 2e3 Pa s) is deliberately absent:
an explicit-force dashpot on massless Lagrangian vertices has a loop gain
c*dt*g/m ~ 1e2-1e3 (per vertex or modal: a modal moment ends up as forces on the
few end vertices, whose local inertia is ~rho*dx^2), so it is either unstable or
bounded to ~1 % of the intended value. It needs implicit fluid-structure
coupling; see NOTES.md.

wing.user_force columns (one row per bending vertex):
    L  M  R  k_theta  theta0
'''

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def _angle(ax, ay, bx, by):
    return np.arctan2(ax * by - ay * bx, ax * bx + ay * by)


_rec = None
_wing = None
_springs = None


def spring_forces(ds, Nb, xLag, yLag, springs, Lx, Ly):
    '''Vectorised give_Me_Spring_Lagrangian_Force_Densities (same arguments).'''
    i, j = springs[:, 0].astype(int), springs[:, 1].astype(int)
    k, rl, a = springs[:, 2], springs[:, 3], springs[:, 4]
    dx, dy = xLag[j] - xLag[i], yLag[j] - yLag[i]
    dx = np.where(np.abs(dx) > Lx / 2, np.sign(dx) * (Lx - np.sign(dx) * dx), dx)
    dy = np.where(np.abs(dy) > Ly / 2, np.sign(dy) * (Ly - np.sign(dy) * dy), dy)
    r = np.sqrt(dx ** 2 + dy ** 2)
    mag = 0.5 * (a + 1.0) * k * np.power(r - rl, a)
    sx, sy = mag * (dx / r), mag * (dy / r)
    fx = np.bincount(i, sx, Nb) - np.bincount(j, sx, Nb)
    fy = np.bincount(i, sy, Nb) - np.bincount(j, sy, Nb)
    return fx, fy


def _diagnostics(ds, Nb, xLag, yLag, dt, current_time):
    global _rec, _wing
    if _rec is None:
        import wing_diagnostics
        m = np.load(os.path.join(_HERE, 'wing_mesh.npz'))
        _wing = dict(tid=m['target_idx'], kT=m['k_Target'] * ds,
                     xT=m['xr'][m['target_idx']], yT=m['yr'][m['target_idx']],
                     xc=float(m['xr'].mean()), yc=float(m['yr'].mean()), L=1.0)   # |offset| < 1 m: inside any sane state
        _rec = wing_diagnostics.Recorder('wing_history', {'pos': (Nb, 2), 'force': (2,)},
                                         dtypes={'pos': np.float32})
        _wing['force_every_step'] = np.full((_rec.n_steps + 1, 2), np.nan)
        _rec.meta = dict(time_offset=0.5 * dt, force_every_step=_wing['force_every_step'])
    step = _rec.step_of(current_time)
    t, kT = _wing['tid'], _wing['kT']
    F = (float(np.sum(kT * (xLag[t] - _wing['xT']))), float(np.sum(kT * (yLag[t] - _wing['yT']))))
    if step < len(_wing['force_every_step']):
        _wing['force_every_step'][step] = F
    # diverged: non-finite, or (after a fluid overflow) finite but absurd
    # positions, which IB2d turns into garbage grid indices (IndexError)
    finite = (np.isfinite(F[0]) and np.isfinite(F[1])
              and np.all(np.abs(xLag - _wing['xc']) < _wing['L']) and np.all(np.abs(yLag - _wing['yc']) < _wing['L']))
    if _rec.wants(step) or (not finite and not _rec.stopped):
        # a NaN anywhere crashes IB2d later in this step (IndexError when it
        # spreads the forces), so record it now and save both histories
        _rec.record(current_time + 0.5 * dt, pos=np.column_stack([xLag, yLag]), force=F)
    if not finite:
        _rec.save()
        fluid = sys.modules.get('please_Compute_External_Forcing')
        if fluid is not None and getattr(fluid, '_rec', None) is not None:
            fluid._rec.save()
    _rec.maybe_save(step)
    return np.zeros(Nb), np.zeros(Nb)


def _hook_springs(ds, Nb, xLag, yLag):
    global _springs
    if _springs is None:
        import wing_case
        with open(os.path.join(_HERE, 'wing.spring')) as f:
            f.readline()
            _springs = (np.loadtxt(f, ndmin=2), wing_case.read_input2d('Lx'), wing_case.read_input2d('Ly'))
    sp, Lx, Ly = _springs
    return spring_forces(ds, Nb, xLag, yLag, sp, Lx, Ly)


def give_Me_General_User_Defined_Force_Densities(ds, Nb, xLag, yLag,
                                                 xLag_P, yLag_P, dt, current_time, general_force):
    if np.ndim(general_force) < 2:               # wing.user_force holds the mode (0 or 1)
        fx, fy = _diagnostics(ds, Nb, xLag, yLag, dt, current_time)
        if int(general_force) == 1:
            fx, fy = _hook_springs(ds, Nb, xLag, yLag)
        return fx, fy
    L = general_force[0, :].astype(int)
    M = general_force[1, :].astype(int)
    R = general_force[2, :].astype(int)
    k_th = general_force[3, :]
    th0 = general_force[4, :]

    ax, ay = xLag[M] - xLag[L], yLag[M] - yLag[L]
    bx, by = xLag[R] - xLag[M], yLag[R] - yLag[M]
    la2, lb2 = ax * ax + ay * ay, bx * bx + by * by

    theta = _angle(ax, ay, bx, by)
    moment = -k_th * (theta - th0)                 # generalized force conjugate to theta

    # d(theta)/dX for the three vertices (gradients sum to zero -> no net force/moment)
    ga = np.stack([-ay, ax]) / la2                  # d angle(a) / d a
    gb = np.stack([-by, bx]) / lb2                  # d angle(b) / d b
    dL, dM, dR = ga, -ga - gb, gb                   # d theta / d X_L, X_M, X_R

    fx = np.zeros(Nb)
    fy = np.zeros(Nb)
    for idx, d in ((L, dL), (M, dM), (R, dR)):
        np.add.at(fx, idx, moment * d[0])
        np.add.at(fy, idx, moment * d[1])
    return fx, fy
