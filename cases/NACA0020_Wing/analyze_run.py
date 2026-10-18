'''
Post-processes a finished run (viz_IB2d/lagsPts.*.vtk, written every
print_dump steps): trailing-edge deflection and the net hydrodynamic force on
the wing, and saves history.npz + history.png into the run folder.

The hydrodynamic force is taken from the clamp: the wing is held only by the
stiff target tethers, so in the quasi-static limit the force the fluid exerts on
the wing (per unit depth) is minus the tether force,
    F_fluid = sum_i k_Target,i * (x_i - x_i,0) * ds     [N/m]
(x_i,0 = the target positions, i.e. the first frame; the .target constants are
in IB2d units, physical k = k_IB * ds). Lift/drag are its components
normal/parallel to the *current* free stream. Multiply by the span (0.1 m in
the WaterLily run) to compare with its forces.

The TE displacement is measured from the straight (rest) wing's trailing
edge, so at t = 0 it is the pre-bend (64.4 mm up for init_q = [1, 0, 0, 0]).

Reads from the run folder: input2d, case_params.json, wing.target,
run_info.json (mesh.te_index). run_case.py runs it on every new run.
Run by hand (from the example folder):
    python analyze_run.py RUN_FOLDER
'''

import glob
import json
import os
import sys

import numpy as np
import vtk
from vtk.util import numpy_support as ns

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import wing_case
from Wing_Geom import CHORD
from run_records import read_input2d
from ib_units import ib_ds


def read_lag(f):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(f)
    r.Update()
    return ns.vtk_to_numpy(r.GetOutput().GetPoints().GetData())[:, :2]


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()
    files = sorted(glob.glob(os.path.join(case, 'viz_IB2d', 'lagsPts.*.vtk')))
    if not files:
        sys.exit('no viz_IB2d/lagsPts.*.vtk in ' + case)
    inp = read_input2d(os.path.join(case, 'input2d'))
    params = json.load(open(os.path.join(case, 'case_params.json')))
    flow, g = params['flow'], params['geom']
    te = json.load(open(os.path.join(case, 'run_info.json')))['mesh']['te_index']
    tgt = np.loadtxt(os.path.join(case, 'wing.target'), skiprows=1, ndmin=2)
    tid, k_t = tgt[:, 0].astype(int), tgt[:, 1] * ib_ds(inp['Lx'], inp['Nx'])   # physical [N/m per m]

    px, py = g.get('pitch_axis_world', [0.4, 0.6])
    te_rest = np.array([px - g.get('pitch_axis_from_LE', 0.03) + CHORD, py])
    P0 = read_lag(files[0])

    t = np.arange(len(files)) * inp['dt'] * inp['print_dump']
    tip = np.zeros((len(files), 2)); F = np.zeros((len(files), 2))
    for n, f in enumerate(files):
        P = read_lag(f)
        tip[n] = P[te] - te_rest
        F[n] = (k_t[:, None] * (P[tid] - P0[tid])).sum(axis=0)

    alpha = np.array([wing_case.aoa_rad(tt, flow) for tt in t])
    drag = F[:, 0] * np.cos(alpha) + F[:, 1] * np.sin(alpha)
    lift = -F[:, 0] * np.sin(alpha) + F[:, 1] * np.cos(alpha)
    np.savez(os.path.join(case, 'history.npz'), t=t, tip=tip, F=F, lift=lift, drag=drag, alpha=alpha)

    print('   t[s]   TE dx[mm]   TE dy[mm]   lift[N/m]  drag[N/m]')
    for n in range(0, len(t), max(1, len(t) // 20)):
        print('%7.3f  %9.2f  %9.2f  %10.4f  %9.4f' % (t[n], 1e3 * tip[n, 0], 1e3 * tip[n, 1], lift[n], drag[n]))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.5))
    ax[0].plot(t, 1e3 * tip[:, 1]); ax[0].set_title('TE y-displacement from straight [mm]')
    ax[1].plot(t, lift); ax[1].set_title('lift per unit depth [N/m]')
    ax[2].plot(t, drag); ax[2].set_title('drag per unit depth [N/m]')
    for a in ax:
        a.set_xlabel('t [s]'); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(case, 'history.png'), dpi=140)


if __name__ == '__main__':
    main()
