'''
Cross-sweep summary of the NACA0020 study (IMPLEMENTATION_DECISIONS.md §4):
reads every run_info.json under the data root and writes into
<data root>/SUMMARY/:

  stability_map.png    every full-schedule / diagnostic run in the (Nx, mu)
                       plane: stable, or the time it diverged; cell Re
                       (U dx / nu) contours
  best_vs_reference.png  TE, Fx, Fy of the stable runs at the lowest viscosity
                       per grid, against SoRoLily
  runs.csv             one row per run (sweep, case, Nx, mu/water, dt, cell
                       Re, stable, divergence time, key measures)

    python summarise.py [DATA_ROOT]     (default ~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing)
'''

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / 'reference'))
import sorolily as sl
from postprocess import REF_FILE, reference_metrics

WATER_MU = sl.RHO * sl.NU
SWEEPS = ('stage3_dt', 'stage4_ladder', 'stage4_diag', 'stage4_refine', 'stage4_mesh', 'stage5_sensitivity',
          'stage5_bend_sign')
COLS = ['te.y_mean_0deg', 'te.y_std_0deg', 'te.y_mean_15deg', 'te.y_std_15deg', 'te.t50',
        'te.rms_diff_release', 'force.fx_mean_0deg', 'force.fy_mean_0deg', 'force.fx_mean_15deg',
        'force.fy_mean_15deg', 'shape.rms_dist_mean_0deg', 'shape.rms_dist_mean_15deg',
        'pressure.rms_diff_mean_0deg', 'pressure.rms_diff_mean_15deg', 'clamp_drift_max']


def get(d, path):
    for p in path.split('.'):
        if not isinstance(d, dict) or p not in d:
            return None
        d = d[p]
    return d


def load(root):
    rows = []
    for sw in SWEEPS:
        for p in sorted((root / sw / 'results').glob('*/run_info.json')):
            info = json.loads(p.read_text())
            i, g, a = info['inputs']['input2d'], info['inputs']['geom'], info.get('analysis', {})
            dx = i['Lx'] / i['Nx']
            stable = a.get('stable')
            if stable is None:
                stable = info['outcome']['status'] == 'completed' and not info['outcome'].get('exit_code')
            rows.append(dict(sweep=sw, case=info['case_name'], dir=p.parent, Nx=int(i['Nx']),
                             mu_w=round(i['mu'] / WATER_MU, 3), dt=i['dt'], Tfinal=i['Tfinal'],
                             cell_re=sl.U * dx * i['rho'] / i['mu'], ktf=g['target_stiffness_factor'],
                             init_q0=g['init_q'][0], mult=g['fill_density_multiplier'], stable=bool(stable),
                             t_div=a.get('first_nonfinite_time'), ms_step=info['timing'].get('sim_time_per_step_ms'),
                             **{c: get(a, c) for c in COLS}))
    return rows


def main(root):
    root = Path(root)
    out = root / 'SUMMARY'
    out.mkdir(exist_ok=True)
    rows = load(root)
    keys = [k for k in rows[0] if k != 'dir']
    with open(out / 'runs.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, keys, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # stability map: runs that use the production setup (ktf 3, pre-bend), full or diagnostic length
    fig, ax = plt.subplots(figsize=(8, 5.5))
    sel = [r for r in rows if r['sweep'] != 'stage3_dt' and r['ktf'] == 3.0]
    for r in sel:
        shift = {1.0: 0.95, 8.0: 1.0}.get(r['mult'], 1.0) * (1.05 if r['init_q0'] == 0 else 1.0)
        x, y = r['Nx'] * shift, r['mu_w']
        if r['stable']:
            ax.plot(x, y, 'o', color='C2', ms=9)
        else:
            ax.plot(x, y, 'x', color='C3', ms=9, mew=2)
            ax.annotate(f"{r['t_div']:.2g} s" if r['t_div'] else '', (x, y), textcoords='offset points',
                        xytext=(6, 4), fontsize=7, color='C3')
    nx = np.geomspace(80, 460, 50)
    for cre in (20, 42, 62, 125, 1250):
        ax.plot(nx, sl.U * 1.2 / nx / (cre * sl.NU), 'k:', lw=0.6)
        ax.text(nx[-1], sl.U * 1.2 / nx[-1] / (cre * sl.NU), f' cell Re {cre}', fontsize=7, va='center')
    ax.set_xscale('log'); ax.set_yscale('log')
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([96, 192, 384]); ax.set_xticklabels(['96', '192', '384'])
    ax.set_xlabel('Nx (cells per chord = Nx / 6)'); ax.set_ylabel('viscosity / water')
    ax.set_title('IB2d runs (target factor 3): o stable for the whole run, x diverged at the time shown\n'
                 'marker left of the column: structural mesh x1 (else x4/x8); right: no pre-bend', fontsize=9)
    fig.tight_layout(); fig.savefig(out / 'stability_map.png', dpi=130); plt.close(fig)

    # best stable runs vs reference: the lowest mu that is stable for the full 10 s, per Nx
    ref = np.load(REF_FILE)
    full = [r for r in rows if r['stable'] and r['Tfinal'] >= 10 and r['sweep'] in ('stage4_ladder', 'stage4_refine')]
    best = {}
    for r in sorted(full, key=lambda r: r['mu_w']):
        best.setdefault(r['Nx'], r)
    fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    ax[0].plot(ref['t'], 1e3 * ref['te_disp'][:, 1], 'k', lw=2, label='SoRoLily (water, Re 4.0e4)')
    ax[1].plot(ref['t'], ref['fx_stored'] / sl.SPAN, 'k', lw=2)
    ax[2].plot(ref['t'], ref['fy_p'] / sl.SPAN, 'k', lw=2)
    for k, (Nxb, r) in enumerate(sorted(best.items())):
        d = np.genfromtxt(r['dir'] / 'timeseries.csv', delimiter=',', names=True)
        lab = f"IB2d Nx {Nxb}, {r['mu_w']:g}x water (Re {sl.U * sl.CHORD / (r['mu_w'] * sl.NU):.3g})"
        ax[0].plot(d['t'], 1e3 * d['te_dy'], color=f'C{k}', label=lab)
        ax[1].plot(d['t'], d['fx_filtered'], color=f'C{k}')
        ax[2].plot(d['t'], d['fy_filtered'], color=f'C{k}')
    ax[0].set_ylabel('TE y displacement [mm]'); ax[0].legend(fontsize=8)
    ax[1].set_ylabel('Fx [N/m], 1.5 Hz low-pass'); ax[2].set_ylabel('Fy [N/m], 1.5 Hz low-pass')
    ax[2].set_xlabel('t [s]')
    for a in ax:
        a.axvline(5.0, color='0.7', lw=0.8)
    fig.tight_layout(); fig.savefig(out / 'best_vs_reference.png', dpi=130); plt.close(fig)

    rm = reference_metrics(ref)
    (out / 'reference_metrics.json').write_text(json.dumps(rm, indent=2) + '\n')
    print('wrote', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing'))
