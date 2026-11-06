'''
Per-step diagnostics for the NACA0020 wing, recorded from inside IB2d through
its hooks (IB2d itself is unchanged):

  * give_Me_General_User_Defined_Force_Densities.py (user_force_model = 1)
    gets the half-step Lagrangian positions every step -> every Lagrangian
    point and the hydrodynamic force on the wing (from the clamp tethers);
  * please_Compute_External_Forcing.py (arb_ext_force = 1) gets the fluid
    velocity every step -> domain-mean velocity, kinetic energy, its
    grid-scale part and the peak speed.

Both sample every `diagnostics.interval` seconds of simulated time (0.01 s by
default, the reference's sampling), independently of print_dump, so the .vtk
field output can be sparse. Each writes an .npz into viz_IB2d/ (which
run_case.py archives with the run), rewritten every `diagnostics.save_every`
seconds and at the last step, so a run that diverges or is killed keeps its
history up to then. Sampling stops at the first non-finite value.
'''

import os

import numpy as np

import wing_case

_HERE = os.path.dirname(os.path.abspath(__file__))


class Recorder:
    '''Fixed-interval sampler of named arrays, saved to viz_IB2d/<name>.npz.'''

    def __init__(self, name, shapes, dtypes=None):
        params = wing_case.load_case_params()
        cfg = params.get('diagnostics', {})
        self.dt = wing_case.read_input2d('dt')
        self.t_final = wing_case.read_input2d('Tfinal')
        self.every = max(1, int(round(cfg.get('interval', 0.01) / self.dt)))
        self.save_steps = max(1, int(round(cfg.get('save_every', 1.0) / self.dt)))
        self.n_steps = int(round(self.t_final / self.dt))
        n = self.n_steps // self.every + 1
        dtypes = dtypes or {}
        self.data = {k: np.full((n,) + tuple(s), np.nan, dtype=dtypes.get(k, np.float64))
                     for k, s in shapes.items()}
        self.t = np.full(n, np.nan)
        self.k = 0
        self.stopped = False
        self.path = os.path.join(_HERE, 'viz_IB2d', name + '.npz')
        self.meta = {}

    def step_of(self, current_time):
        return int(round(current_time / self.dt))

    def wants(self, step):
        return not self.stopped and step % self.every == 0 and self.k < len(self.t)

    def record(self, t, **values):
        for key, v in values.items():
            self.data[key][self.k] = v
        self.t[self.k] = t
        self.k += 1
        if not all(np.all(np.isfinite(np.asarray(v, dtype=float))) for v in values.values()):
            self.stopped = True          # NaN never clears: keep the first bad sample, then stop
            self.save()

    def maybe_save(self, step):
        if step % self.save_steps == 0 or step >= self.n_steps - 1:
            self.save()

    def save(self):
        tmp = self.path + '.tmp.npz'
        np.savez(tmp, t=self.t[:self.k], interval=self.every * self.dt,
                 **{k: v[:self.k] for k, v in self.data.items()}, **self.meta)
        os.replace(tmp, self.path)
