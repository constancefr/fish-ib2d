'''
Summarises one sweep folder (as made by slurm/submit_sweep.sh) into
<sweep>/report/: summary.md (tables), and figures. Reads only the run
folders' run_info.json, timeseries.csv and viz_IB2d/wing_history.npz, so it
can be re-run at any time.

    python report.py <sweep folder> [...]

  * dt sweeps (inputs.sweep starts with "stage3"): stability map per
    (Nx, target_stiffness_factor, mu) against dt, and dt_crit = the largest dt
    of the ladder that is stable together with every smaller dt.
  * every sweep: table of the comparison measures (postprocess.py) next to
    the reference's, and overlays of the TE displacement and the filtered
    forces against SoRoLily (te_overlay.png, force_overlay.png).
'''

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / 'reference'))
import sorolily as sl
from postprocess import REF_FILE, reference_metrics

WATER_MU = sl.RHO * sl.NU


def load_runs(sweep_dir):
    runs = []
    for p in sorted((Path(sweep_dir) / 'results').glob('*/run_info.json')):
        info = json.loads(p.read_text())
        info['_dir'] = p.parent
        runs.append(info)
    return runs


def key(info):
    i, g = info['inputs']['input2d'], info['inputs']['geom']
    return int(i['Nx']), g['target_stiffness_factor'], round(i['mu'] / WATER_MU, 3)


def stable(info):
    a = info.get('analysis', {})
    if 'stable' in a and a['stable'] is not None:
        return bool(a['stable'])
    return info['outcome']['status'] == 'completed' and not info['outcome'].get('exit_code')


def fmt(x, scale=1.0, nd=3):
    if x is None:
        return '-'
    return f'{x * scale:.{nd}g}'


def dt_report(runs, out, lines):
    groups = {}
    for r in runs:
        groups.setdefault(key(r), []).append(r)
    lines += ['## Critical time step (1 s runs: start-up ramp + release of the pre-bend)', '',
              'Each cell: S = stable through 1 s, or the time [s] of the first non-finite value. '
              'ms/step in brackets (wall time / steps, on whichever node ran it).', '']
    dts = sorted({r['inputs']['input2d']['dt'] for r in runs}, reverse=True)
    lines.append('| Nx | k_T factor | mu / water | ' + ' | '.join(f'{d:g}' for d in dts) + ' | dt_crit |')
    lines.append('|' + '---|' * (len(dts) + 4))
    crit = {}
    for k in sorted(groups):
        row = {r['inputs']['input2d']['dt']: r for r in groups[k]}
        cells = []
        for d in dts:
            r = row.get(d)
            if r is None:
                cells.append('')
                continue
            ms = r['timing'].get('sim_time_per_step_ms')
            s = 'S' if stable(r) else fmt(r.get('analysis', {}).get('first_nonfinite_time'), nd=2)
            cells.append(f'{s} ({ms:.0f})' if ms else s)
        dc = None
        for d in sorted(row):                      # ascending: stop at the first failure
            if stable(row[d]):
                dc = d
            else:
                break
        crit[k] = dc
        lines.append(f'| {k[0]} | {k[1]:g} | {k[2]:g} | ' + ' | '.join(cells) + f' | {fmt(dc)} |')
    lines.append('')
    (out / 'dt_crit.json').write_text(json.dumps(
        [dict(Nx=k[0], target_stiffness_factor=k[1], mu_over_water=k[2], dt_crit=v) for k, v in crit.items()],
        indent=2) + '\n')
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.5))
    labels = []
    for y, k in enumerate(sorted(groups)):
        for r in groups[k]:
            d = r['inputs']['input2d']['dt']
            ax.plot(d, y, 'o' if stable(r) else 'x', color='C2' if stable(r) else 'C3', ms=8)
        labels.append(f'Nx {k[0]}, k_T x{k[1]:g}, mu {k[2]:g}x')
    ax.set_xscale('log'); ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('dt [s]'); ax.set_title('1 s runs: o stable, x diverged')
    fig.tight_layout(); fig.savefig(out / 'dt_stability.png', dpi=120); plt.close(fig)
    return crit


MEASURES = [
    ('te.y_initial', 'TE y(0) [mm]', 1e3), ('te.t50', 't50 [s]', 1), ('te.t10', 't10 [s]', 1),
    ('te.y_min_release', 'TE min, release [mm]', 1e3),
    ('te.y_mean_0deg', 'TE mean 0° [mm]', 1e3), ('te.y_std_0deg', 'TE std 0° [mm]', 1e3),
    ('te.y_mean_15deg', 'TE mean 15° [mm]', 1e3), ('te.y_std_15deg', 'TE std 15° [mm]', 1e3),
    ('force.fx_mean_0deg', 'Fx mean 0° [N/m]', 1), ('force.fy_mean_0deg', 'Fy mean 0° [N/m]', 1),
    ('force.fx_mean_15deg', 'Fx mean 15° [N/m]', 1), ('force.fy_mean_15deg', 'Fy mean 15° [N/m]', 1),
    ('force.fy_std_15deg', 'Fy std 15° [N/m]', 1),
]
DIFFS = [
    ('te.rms_diff_release', 'TE rms diff release [mm]', 1e3), ('te.rms_diff_0deg', 'TE rms diff 0° [mm]', 1e3),
    ('te.rms_diff_15deg', 'TE rms diff 15° [mm]', 1e3),
    ('force.rms_diff_fx_15deg', 'Fx rms diff 15° [N/m]', 1), ('force.rms_diff_fy_15deg', 'Fy rms diff 15° [N/m]', 1),
    ('shape.rms_dist_t0p5', 'centreline rms dist t=0.5 [mm]', 1e3),
    ('shape.rms_dist_mean_15deg', 'mean centreline rms dist 15° [mm]', 1e3),
    ('pressure.rms_diff_mean_15deg', 'surface p rms diff 15° [Pa]', 1),
    ('fluid.gs_ke_growth_rate', 'grid-scale KE growth [1/s]', 1), ('fluid.cell_re_max', 'max cell Re', 1),
]


def get(d, path):
    for p in path.split('.'):
        if not isinstance(d, dict) or p not in d:
            return None
        d = d[p]
    return d


def label(r):
    i = r['inputs']['input2d']
    return f"Nx {int(i['Nx'])}, dt {i['dt']:g}, mu {i['mu'] / WATER_MU:g}x" + (
        '' if stable(r) else ' (diverged)')


def physics_report(runs, out, lines):
    ref = np.load(REF_FILE)
    rm = reference_metrics(ref)
    cols = [r for r in runs if r.get('analysis')]
    lines += ['## Comparison with SoRoLily', '',
              '| measure | SoRoLily | ' + ' | '.join(r['case_name'] for r in cols) + ' |',
              '|' + '---|' * (len(cols) + 2)]
    lines.append('| status | | ' + ' | '.join(
        'stable' if stable(r) else f"NaN at {fmt(get(r['analysis'], 'first_nonfinite_time'), nd=3)} s"
        for r in cols) + ' |')
    lines.append('| Re (chord) | %.3g | ' % (sl.U * sl.CHORD / sl.NU) + ' | '.join(
        fmt(get(r['analysis'], 're_chord'), nd=3) for r in cols) + ' |')
    for path, name, s in MEASURES:
        lines.append(f'| {name} | {fmt(get(rm, path), s)} | ' + ' | '.join(
            fmt(get(r['analysis'], path), s) for r in cols) + ' |')
    for path, name, s in DIFFS:
        lines.append(f'| {name} | | ' + ' | '.join(fmt(get(r['analysis'], path), s) for r in cols) + ' |')
    lines.append('| ms / step | | ' + ' | '.join(fmt(r['timing'].get('sim_time_per_step_ms'), nd=3)
                                                  for r in cols) + ' |')
    lines.append('')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax[0].plot(ref['t'], 1e3 * ref['te_disp'][:, 1], 'k', lw=2, label='SoRoLily')
    fig2, ax2 = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax2[0].plot(ref['t'], ref['fx_stored'] / sl.SPAN, 'k', lw=2, label='SoRoLily')
    ax2[1].plot(ref['t'], ref['fy_p'] / sl.SPAN, 'k', lw=2, label='SoRoLily')
    for k, r in enumerate(cols):
        ts = r['_dir'] / 'timeseries.csv'
        if not ts.exists():
            continue
        d = np.genfromtxt(ts, delimiter=',', names=True)
        if d.size < 2:
            continue
        c = f'C{k % 10}'
        ax[0].plot(d['t'], 1e3 * d['te_dy'], color=c, lw=1, label=label(r))
        ax[1].plot(d['t'], 1e3 * (d['te_dy'] - d['te_dy_ref']), color=c, lw=1)
        ax2[0].plot(d['t'], d['fx_filtered'], color=c, lw=1, label=label(r))
        ax2[1].plot(d['t'], d['fy_filtered'], color=c, lw=1)
    ax[0].set_ylabel('TE y displacement [mm]'); ax[0].legend(fontsize=7, ncol=2)
    ax[1].set_ylabel('IB2d - SoRoLily [mm]'); ax[1].set_xlabel('t [s]'); ax[1].axhline(0, color='k', lw=0.5)
    fig.tight_layout(); fig.savefig(out / 'te_overlay.png', dpi=120); plt.close(fig)
    ax2[0].set_ylabel('Fx [N/m] (filtered 1.5 Hz)'); ax2[0].legend(fontsize=7, ncol=2)
    ax2[1].set_ylabel('Fy [N/m] (filtered 1.5 Hz)'); ax2[1].set_xlabel('t [s]')
    for a in ax2:
        a.set_ylim(-6, 10)
    fig2.tight_layout(); fig2.savefig(out / 'force_overlay.png', dpi=120); plt.close(fig2)


def main(dirs):
    for sd in dirs:
        sd = Path(sd)
        runs = load_runs(sd)
        out = sd / 'report'
        out.mkdir(exist_ok=True)
        sweep = runs[0]['inputs'].get('sweep', '') if runs else ''
        lines = [f'# Sweep {sd.name} ({sweep}): {len(runs)} runs', '']
        if sweep.startswith('stage3'):
            dt_report(runs, out, lines)
        else:
            physics_report(runs, out, lines)
        (out / 'summary.md').write_text('\n'.join(lines) + '\n')
        print('wrote', out / 'summary.md')


if __name__ == '__main__':
    main(sys.argv[1:])
