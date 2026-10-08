'''
User-defined fibre force for the NACA0020 wing (IB2d `user_force_model = 1`,
called by the driver every step; rows come from wing.user_force).

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

import numpy as np


def _angle(ax, ay, bx, by):
    return np.arctan2(ax * by - ay * bx, ax * bx + ay * by)


def give_Me_General_User_Defined_Force_Densities(ds, Nb, xLag, yLag,
                                                 xLag_P, yLag_P, dt, current_time, general_force):
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
