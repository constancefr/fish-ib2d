'''
Free-stream forcing for the NACA0020 wing (IB2d `arb_ext_force = 1`).

IB2d is doubly periodic, so the WaterLily setup (uniform inflow of speed U at
angle alpha(t) on the left, free outflow on the right) is reproduced with a
"fringe" (sponge) region: in a strip at the downstream end of the domain the
fluid is relaxed toward the free-stream velocity,

    f = -rho * lambda * chi(x) * (u - U (cos alpha(t), sin alpha(t)))     [N/m^3]

Anything leaving the test section (wake, vortices) is damped in the fringe and
the fluid re-enters at x = 0 (periodic wrap) as clean free stream. The angle of
attack is changed by rotating the free-stream vector, exactly as the WaterLily
example rotates its inflow (the wing itself never pitches).

The fringe alone changes the flow direction only where it acts; the turned
fluid then has to travel round the periodic box to reach the wing (~2 s here),
whereas an inflow BC in a bounded domain changes the far-field potential flow
almost instantly. A uniform body force on the domain-mean velocity (the k = 0
mode), f = -rho * lambda_mean * (<u> - u_target), restores that: it rotates the
whole free stream at once and leaves vorticity to be advected as usual.

Tunables live in case_params.json -> "flow".
'''

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wing_case

_flow = wing_case.load_case_params()['flow']
_rho = wing_case.read_input2d('rho')


def _mask(x, Lx):
    '''Smooth 0..1 profile: ramps up over fringe_ramp at fringe_x[0], ramps down at fringe_x[1].'''
    x0, x1 = _flow['fringe_x']
    r = _flow['fringe_ramp']
    up = np.array([wing_case.smoothstep((xi - x0) / r) for xi in x])
    down = np.array([1.0 - wing_case.smoothstep((xi - (x1 - r)) / r) for xi in x])
    return up * down


def please_Compute_External_Forcing(dt, current_time, x, y, grid_Info, uX, uY, first, inds):
    if first == 1:
        inds = _mask(x, grid_Info[2])[None, :]      # (1, Nx) fringe weights, cached between calls
        first = 0

    alpha = wing_case.aoa_rad(current_time, _flow)
    uX_tar = _flow['U'] * np.cos(alpha)
    uY_tar = _flow['U'] * np.sin(alpha)
    k = _rho * _flow['fringe_rate'] * inds
    km = _rho * _flow.get('mean_rate', 0.0)

    return (-k * (uX - uX_tar) - km * (uX.mean() - uX_tar),
            -k * (uY - uY_tar) - km * (uY.mean() - uY_tar), first, inds)
