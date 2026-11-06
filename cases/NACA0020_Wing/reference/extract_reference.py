'''
Extracts the SoRoLily reference run (data_gvs/gvs_tests.h5, written by
gvs_tests.jl) into sorolily_reference.npz, in physical units, with everything
the IB2d comparison needs, and plots reference_overview.png.

    python reference/extract_reference.py [path/to/gvs_tests.h5 [OUT.npz]]

(default: ~/finray_fsi/WaterLily_FSI/WaterLily_2Way_Coupling_Example/data_gvs/gvs_tests.h5
-> reference/sorolily_reference.npz, the committed reference, so comparisons
don't need the .h5 or h5py). With OUT.npz (e.g. a SoRoLily re-run with other
settings) the overview plot goes next to it as <OUT>_overview.png.

Contents (time t every 0.01 s, 0 .. 10.01 s; n = 1002 samples):
  t, aoa [rad]                 time, angle of the free stream
  q (n, 4), qd (n, 4)          curvature coordinates of link 3 and their rates
  te (n, 2) [m]                TE position (world); te_disp = te - straight TE
  centreline_s (21,) [m]       arc length along link 3 of the 21 soft stations
  centreline (n, 21, 2) [m]    backbone at those stations (world)
  fx_stored (n,) [N]           force_history row 2: sum Fx, pressure + viscous,
                               filtered (Butterworth 6, 1.5 Hz), 0.1 m span
  fx_p, fy_p (n,) [N]          pressure force rebuilt from the stored panel
                               pressures (filtered at 3 Hz), 0.1 m span, then
                               filtered like force_history (1.5 Hz)
  gen_force (n, 6)             force_history rows 3-8 (generalized forces)
  p_panels (n, 30) [Pa]        panel pressures x rho U^2 (gauge: the reference
                               subtracts p at the domain's corner cell)
  panel_x_rest (30,), panel_side (30,)   chordwise rest position and side of
                               each panel (fixed through the run)
  per-unit-depth forces: divide the [N] values by SPAN = 0.1 m.
'''

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sorolily as sl

DEFAULT_H5 = os.path.expanduser('~/finray_fsi/WaterLily_FSI/WaterLily_2Way_Coupling_Example/data_gvs/gvs_tests.h5')


def main(h5_path=DEFAULT_H5, out_path=None):
    import h5py
    from scipy.signal import lfilter
    with h5py.File(h5_path, 'r') as f:
        state = f['state_history'][:]          # (n, 12): q (6) then qd (6)
        force = f['force_history'][:]          # (n, 8)
        t = f['time_history'][:]
        p = f['pressure_history'][:]           # (n, 30), simulation units
        aoa = f['aoa_history'][:]

    q, qd = state[:, 2:6], state[:, 8:12]
    assert np.allclose(state[:, [0, 1, 6, 7]], 0.0), 'carriage not frozen?'

    n = len(t)
    cl = np.zeros((n, sl.N_SOFT, 2))
    fxy = np.zeros((n, 2))
    for i in range(n):
        cl[i] = sl.stations(q[i])[1]['c']
        fxy[i] = sl.pressure_force(q[i], p[i])
    # force_history is the 1.5 Hz-filtered panel sum of the 3 Hz-filtered
    # pressure; apply the same second filter to the rebuilt pressure force
    b, a = sl.butterworth(sl.FORCE_CUTOFF_HZ)
    fx_p, fy_p = lfilter(b, a, fxy[:, 0]), lfilter(b, a, fxy[:, 1])
    te = cl[:, -1]
    pn = sl.panels(q[0])

    out = out_path or os.path.join(HERE, 'sorolily_reference.npz')
    np.savez_compressed(
        out, t=t, aoa=aoa, q=q, qd=qd, te=te, te_disp=te - sl.straight_te(),
        centreline_s=sl.L3 * np.linspace(0, 1, sl.N_SOFT), centreline=cl,
        fx_stored=force[:, 1], fx_p=fx_p, fy_p=fy_p, gen_force=force[:, 2:],
        p_panels=p * sl.P_SCALE, panel_x_rest=pn['x_rest'], panel_side=pn['side'],
        source=os.path.abspath(h5_path))

    # check the rebuilt pressure force against the stored total Fx (which also
    # holds the viscous term): mean and RMS difference after the start-up
    m = t > 1.0
    print('wrote', out)
    print('Fx stored vs rebuilt from pressure (t > 1 s): mean %.4f vs %.4f N, rms diff %.4f N'
          % (force[m, 1].mean(), fx_p[m].mean(), np.sqrt(np.mean((force[m, 1] - fx_p[m]) ** 2))))
    print('gen_force[heave] vs rebuilt Fy (t > 1 s): mean %.4f vs %.4f N'
          % (force[m, 2].mean(), fy_p[m].mean()))
    print('TE y displacement at t = 0, 0.5, 1, 2 s: %s mm'
          % ', '.join('%.2f' % (1e3 * np.interp(tt, t, te[:, 1] - sl.straight_te()[1])) for tt in (0, 0.5, 1, 2)))
    png = (os.path.splitext(out)[0] + '_overview.png') if out_path else os.path.join(HERE, 'reference_overview.png')
    plot(t, aoa, te - sl.straight_te(), force, fx_p, fy_p, q, png, h5_path)


def plot(t, aoa, ted, force, fx_p, fy_p, q, png, h5_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    ax[0, 0].plot(t, 1e3 * ted[:, 1], 'k'); ax[0, 0].set_ylabel('TE y displacement [mm]')
    ax[0, 0].set_ylim(-10, 70)
    ax[0, 1].plot(t, q); ax[0, 1].set_ylabel('curvature coordinates q'); ax[0, 1].legend(['q1', 'q2', 'q3', 'q4'])
    ax[1, 0].plot(t, force[:, 1] / sl.SPAN, 'k', label='Fx stored (pressure + viscous)')
    ax[1, 0].plot(t, fx_p / sl.SPAN, 'C0--', label='Fx rebuilt (pressure)')
    ax[1, 0].plot(t, fy_p / sl.SPAN, 'C3', label='Fy rebuilt (pressure)')
    ax[1, 0].set_ylabel('force per unit depth [N/m]'); ax[1, 0].legend(fontsize=8); ax[1, 0].set_ylim(-5, 10)
    ax[1, 1].plot(t, np.degrees(aoa), 'k'); ax[1, 1].set_ylabel('AoA [deg]')
    for a in ax[1]:
        a.set_xlabel('t [s]')
    fig.suptitle('SoRoLily (%s)' % os.path.relpath(h5_path, os.path.expanduser('~')))
    fig.tight_layout(); fig.savefig(png, dpi=130)


if __name__ == '__main__':
    main(*sys.argv[1:])
