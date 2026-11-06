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

Initial state (flow.start_ramp): IB2d starts the fluid at rest
(IBM_Driver.py: U = np.zeros), whereas WaterLily starts it at the free stream
everywhere (Flow(): u0 = ic_function(uBC), i.e. u = U at t = 0, body
included). With start_ramp = tau > 0 the domain-mean velocity is driven
along U * smoothstep(t / tau) exactly (deadbeat: the uniform force
rho * (u_target(t + dt) - <u>) / dt each step, which in IB2d's periodic
solver changes only the k = 0 mode, by exactly f dt / rho), and the fringe
relaxes towards the same ramped target; after tau the normal fringe + mean
controller take over. A one-step impulse (tau = 0) is the literal WaterLily
start, but it excites the clamp's ringing mode (the wing's interior fluid on
the stiff target tethers, period ~2 ms) to +-4000 N/m, undamped over 0.1 s;
tau = 0.02 s (~10 ringing periods, 2 % of the convective time c/U = 1 s)
starts the flow as fast as the clamp allows (IMPLEMENTATION_DECISIONS.md, O2).

Diagnostics: every diagnostics.interval seconds the fluid's domain-mean
velocity, kinetic energy per unit depth KE = rho/2 sum(u^2+v^2) dx^2, its
grid-scale part KE_gs (the same sum over u minus its 3x3 periodic box
average: the measure of cases/FinRay/postprocess.py) and the peak speed go
to viz_IB2d/fluid_history.npz (wing_diagnostics.py).

Tunables live in case_params.json -> "flow".
'''

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wing_case
import wing_diagnostics

_flow = wing_case.load_case_params()['flow']
_rho = wing_case.read_input2d('rho')


def _mask(x, Lx):
    '''Smooth 0..1 profile: ramps up over fringe_ramp at fringe_x[0], ramps down at fringe_x[1].'''
    x0, x1 = _flow['fringe_x']
    r = _flow['fringe_ramp']
    up = np.array([wing_case.smoothstep((xi - x0) / r) for xi in x])
    down = np.array([1.0 - wing_case.smoothstep((xi - (x1 - r)) / r) for xi in x])
    return up * down


_rec = None


def _box3(u):
    '''3x3 box average on a periodic grid.'''
    s = sum(np.roll(np.roll(u, i, axis=0), j, axis=1) for i in (-1, 0, 1) for j in (-1, 0, 1))
    return s / 9.0


def _record(dt, current_time, uX, uY, dx):
    global _rec
    if _rec is None:
        _rec = wing_diagnostics.Recorder('fluid_history', {
            'u_mean': (2,), 'ke': (), 'ke_gs': (), 'u_max': ()})
    step = _rec.step_of(current_time)
    if _rec.wants(step):
        ru, rv = uX - _box3(uX), uY - _box3(uY)
        _rec.record(current_time, u_mean=(uX.mean(), uY.mean()),
                    ke=0.5 * _rho * np.sum(uX ** 2 + uY ** 2) * dx ** 2,
                    ke_gs=0.5 * _rho * np.sum(ru ** 2 + rv ** 2) * dx ** 2,
                    u_max=np.sqrt(uX ** 2 + uY ** 2).max())
    _rec.maybe_save(step)


def please_Compute_External_Forcing(dt, current_time, x, y, grid_Info, uX, uY, first, inds):
    if first == 1:
        inds = _mask(x, grid_Info[2])[None, :]      # (1, Nx) fringe weights, cached between calls
        first = 0
    _record(dt, current_time, uX, uY, grid_Info[4])

    alpha = wing_case.aoa_rad(current_time, _flow)
    tau = _flow.get('start_ramp', 0.0)
    k = _rho * _flow['fringe_rate'] * inds
    km = _rho * _flow.get('mean_rate', 0.0)
    ramp = 1.0
    if tau > 0 and current_time < tau:
        ramp = wing_case.smoothstep((current_time + dt) / tau)   # target at the end of this step
        km = _rho / dt                                           # deadbeat on the mean
    uX_tar = ramp * _flow['U'] * np.cos(alpha)
    uY_tar = ramp * _flow['U'] * np.sin(alpha)

    return (-k * (uX - uX_tar) - km * (uX.mean() - uX_tar),
            -k * (uY - uY_tar) - km * (uY.mean() - uY_tar), first, inds)
