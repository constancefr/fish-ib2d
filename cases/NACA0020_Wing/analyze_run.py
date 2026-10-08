'''
Post-processes a finished run (viz_IB2d/lagsPts.*.vtk, written every
print_dump steps): tip deflection, tip angle-of-attack-relative motion and the
net hydrodynamic force on the wing, and saves history.npz + history.png.

The hydrodynamic force is taken from the clamp: the wing is held only by the
stiff target tethers, so in the quasi-static limit the force the fluid exerts on
the wing (per unit depth) is minus the tether force,
    F_fluid = sum_i K_target * (x_i - x_i,rest)     [N/m]
Lift/drag are its components normal/parallel to the *current* free stream.
Multiply by the span (0.1 m in the WaterLily run) to compare with its Lift (N).

Run (from the example folder, after main2d.py):
    python analyze_run.py [case_dir]
'''

import glob
import os
import sys

import numpy as np
import vtk
from vtk.util import numpy_support as ns

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Wing_Geom as WG
import wing_case


def read_lag(f):
    r = vtk.vtkUnstructuredGridReader()
    r.SetFileName(f)
    r.Update()
    return ns.vtk_to_numpy(r.GetOutput().GetPoints().GetData())[:, :2]


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()
    files = sorted(glob.glob(os.path.join(case, 'viz_IB2d', 'lagsPts.*.vtk')))
    dt = wing_case.read_input2d('dt')
    dump = wing_case.read_input2d('print_dump')
    flow = wing_case.load_case_params()['flow']
    g = wing_case.load_case_params()['geom']

    w = WG.build_wing(bent=False)
    tid = w['targets'][:, 0].astype(int)
    K_t = g['K_target']
    te = int(w['iM'][w['nc']])
    rest = np.column_stack([w['xr'], w['yr']])

    t = np.arange(len(files)) * dt * dump
    tip = np.zeros((len(files), 2)); F = np.zeros((len(files), 2))
    for n, f in enumerate(files):
        P = read_lag(f)
        tip[n] = P[te] - rest[te]
        F[n] = K_t * (P[tid] - rest[tid]).sum(axis=0)

    alpha = np.array([wing_case.aoa_rad(tt, flow) for tt in t])
    drag = F[:, 0] * np.cos(alpha) + F[:, 1] * np.sin(alpha)
    lift = -F[:, 0] * np.sin(alpha) + F[:, 1] * np.cos(alpha)
    np.savez(os.path.join(case, 'history.npz'), t=t, tip=tip, F=F, lift=lift, drag=drag, alpha=alpha)

    print('   t[s]   tip dx[mm]  tip dy[mm]   lift[N/m]  drag[N/m]')
    for n in range(0, len(t), max(1, len(t) // 20)):
        print('%7.3f  %9.2f  %9.2f  %10.4f  %9.4f' % (t[n], 1e3 * tip[n, 0], 1e3 * tip[n, 1], lift[n], drag[n]))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.5))
    ax[0].plot(t, 1e3 * tip[:, 1]); ax[0].set_title('TE y-displacement [mm]')
    ax[1].plot(t, lift); ax[1].set_title('lift per unit depth [N/m]')
    ax[2].plot(t, drag); ax[2].set_title('drag per unit depth [N/m]')
    for a in ax:
        a.set_xlabel('t [s]'); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(case, 'history.png'), dpi=140)


if __name__ == '__main__':
    main()
