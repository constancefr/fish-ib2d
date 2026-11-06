'''
Post-processing of NACA0020 runs: compares one IB2d run with the SoRoLily
reference (reference/sorolily_reference.npz) on the four outputs of
IMPLEMENTATION_DECISIONS.md decision 6, adds an "analysis" section to the
run's run_info.json and writes comparison.png, flow_snapshots.png and
timeseries.csv into the run folder.

    python postprocess.py <run folder> [<run folder> ...]
    python postprocess.py --ref OTHER.npz --name NAME <run folder> ...   (another SoRoLily run)

run_case.py calls it after every run. It reads only the run folder
(viz_IB2d/wing_history.npz, fluid_history.npz, P.*.vtk, Omega.*.vtk,
wing_mesh.npz, case_params.json, input2d), so it can be re-run on old runs.

Windows (s): release [0, 2], 0 deg [2, 5], 15 deg [5.5, 10] (the AoA switch is
at 5 s). Every measure is computed the same way for the reference
(reference_metrics()).

analysis.te        trailing-edge displacement from the straight wing, y [m]
  y_initial, t50, t10   pre-bend, and the first times it falls to 50 % / 10 %
  y_min_release         minimum over the release window (overshoot if < 0)
  y_mean_*, y_std_*     mean and std over the 0 deg / 15 deg windows
  rms_diff_*            RMS of (IB2d - reference) over each window
analysis.force     hydrodynamic force per unit depth [N/m], world axes (Fx
                   along the 0 deg free stream). IB2d: from the clamp tethers
                   (wing_history.npz), then the reference's own causal filter
                   (Butterworth 6, 1.5 Hz at 100 Hz), as force_history was.
                   Reference: Fx stored, Fy rebuilt from the panel pressures.
  f{x,y}_mean_*, f{x,y}_std_*, rms_diff_f{x,y}_*
analysis.shape     centreline of the flexible link at the reference's 21
                   stations (same material points: IB2d's rest centreline
                   mapped through the deformed triangulation)
  rms_dist_t0p25/t0p5/t1   RMS distance to the reference [m] at those times
  rms_dist_mean_*          the same between the window-mean shapes
analysis.pressure  surface pressure at the reference's 30 panel midpoints
                   (mapped like the centreline), sampled supp/2 cells outside
                   the deformed outline (outside IB2d's smeared pressure jump,
                   as the reference's BDIM value is), gauge = p at
                   (3.125, 3.125) mm as the reference (its corner cell), every
                   field frame
  rms_diff_mean_*   RMS over panels of the window-mean difference [Pa]
  ref_rms_mean_*    RMS over panels of the reference's window mean [Pa] (scale)
analysis.fluid     (fluid_history.npz) u_mean_end: domain-mean velocity at the
                   end [m/s]; gs_ke_growth_rate: slope of ln(KE_gs) over
                   t >= 1 s [1/s] (> 0: grid-scale noise grows, as in
                   cases/FinRay/postprocess.py); gs_ke_fraction_end: KE_gs/KE
                   over the last second; u_max, cell_re_max = u_max dx rho/mu
analysis.stable            no NaN and IB2d exited normally (IB2d usually stops
                           with an exception on NaN positions instead of
                           writing NaN frames, so outcome.status alone can say
                           "incomplete"); first_nonfinite_time = the first step
                           with a non-finite position or tether force (else
                           the last sample of a crashed run)
analysis.clamp_drift_max   largest displacement of a target point [m]
analysis.force_ringing_rms RMS of the per-step tether force minus its 1.5 Hz
                           low-pass [N/m]: the clamp's ringing (IB artefact)
'''

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE / 'reference', HERE.parents[1] / 'common'):
    sys.path.insert(0, str(p))
import run_records as rr
import sorolily as sl

REF_FILE = HERE / 'reference' / 'sorolily_reference.npz'
WINDOWS = {'release': (0.0, 2.0), '0deg': (2.0, 5.0), '15deg': (5.5, 10.0)}
SHAPE_TIMES = {'t0p25': 0.25, 't0p5': 0.5, 't1': 1.0}
GAUGE_POINT = (0.5 * sl.CHORD / sl.CELLS_PER_CHORD,) * 2


def _r(x, digits=6):
    if x is None:
        return None
    x = float(x)
    return float(f'{x:.{digits}g}') if np.isfinite(x) else None


def _win(t, name):
    a, b = WINDOWS[name]
    return (t >= a - 1e-9) & (t <= b + 1e-9)


def _full_window(t, name):
    '''True if the samples t cover the window.'''
    a, b = WINDOWS[name]
    return len(t) > 0 and t[0] <= a + 0.06 and t[-1] >= b - 0.06      # one field interval (0.05 s) of slack


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------

def read_scalar_vtk(path):
    '''(Ny, Nx) array of an IB2d STRUCTURED_POINTS scalar file (x fastest,
    origin 0, spacing dx) -- via VTK, whatever writer IB2d used.'''
    import vtk
    from vtk.util import numpy_support as ns
    r = vtk.vtkStructuredPointsReader()
    r.SetFileName(str(path))
    r.Update()
    out = r.GetOutput()
    nx, ny, _ = out.GetDimensions()
    return ns.vtk_to_numpy(out.GetPointData().GetScalars()).reshape(ny, nx)


def sample_periodic(field, Lx, Ly, px, py):
    '''Bilinear interpolation of a (Ny, Nx) node field (x_i = i dx) at points.'''
    ny, nx = field.shape
    fx, fy = (np.asarray(px) / Lx * nx) % nx, (np.asarray(py) / Ly * ny) % ny
    i0, j0 = np.floor(fx).astype(int), np.floor(fy).astype(int)
    ax, ay = fx - i0, fy - j0
    i1, j1 = (i0 + 1) % nx, (j0 + 1) % ny
    return ((1 - ax) * (1 - ay) * field[j0, i0] + ax * (1 - ay) * field[j0, i1]
            + (1 - ax) * ay * field[j1, i0] + ax * ay * field[j1, i1])


def outward_normals(outline, pts):
    '''Unit outward normals of the closed outline (LE -> upper -> TE -> lower:
    clockwise) at its vertex nearest to each point.'''
    k = np.argmin(np.hypot(outline[None, :, 0] - pts[:, 0, None], outline[None, :, 1] - pts[:, 1, None]), axis=1)
    tng = np.roll(outline, -1, axis=0)[k] - np.roll(outline, 1, axis=0)[k]
    nrm = np.column_stack([-tng[:, 1], tng[:, 0]])
    return nrm / np.linalg.norm(nrm, axis=1)[:, None]


class MaterialMap:
    '''Material points given by rest positions: barycentric weights in the
    rest triangulation, applied to any deformed configuration.'''

    def __init__(self, xr, yr, triangles, qx, qy):
        import matplotlib.tri as mtri
        tri = mtri.Triangulation(xr, yr, triangles)
        k = np.asarray(tri.get_trifinder()(qx, qy))
        cen = np.column_stack([xr[triangles].mean(1), yr[triangles].mean(1)])
        for m in np.where(k < 0)[0]:             # on the outline: nearest triangle (weights extrapolate slightly)
            k[m] = np.argmin(np.hypot(cen[:, 0] - qx[m], cen[:, 1] - qy[m]))
        v = triangles[k]
        x1, y1, x2, y2, x3, y3 = xr[v[:, 0]], yr[v[:, 0]], xr[v[:, 1]], yr[v[:, 1]], xr[v[:, 2]], yr[v[:, 2]]
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        w1 = ((y2 - y3) * (qx - x3) + (x3 - x2) * (qy - y3)) / det
        w2 = ((y3 - y1) * (qx - x3) + (x1 - x3) * (qy - y3)) / det
        self.v, self.w = v, np.column_stack([w1, w2, 1 - w1 - w2])

    def __call__(self, pos):
        '''pos (..., N, 2) -> (..., n, 2)'''
        return np.einsum('...kjd,kj->...kd', pos[..., self.v, :], self.w)


# --------------------------------------------------------------------------
# Measures (same code for IB2d and the reference)
# --------------------------------------------------------------------------

def te_metrics(t, y, t_ref=None, y_ref=None):
    out = dict(y_initial=_r(y[0]) if len(y) else None)
    y0 = y[0] if len(y) else np.nan
    for frac, key in ((0.5, 't50'), (0.1, 't10')):
        below = np.where(y <= frac * y0)[0]
        out[key] = _r(t[below[0]]) if len(below) else None
    m = _win(t, 'release')
    out['y_min_release'] = _r(y[m].min()) if m.any() and _full_window(t, 'release') else None
    for w in ('0deg', '15deg'):
        m = _win(t, w)
        ok = m.any() and _full_window(t, w)
        out[f'y_mean_{w}'] = _r(y[m].mean()) if ok else None
        out[f'y_std_{w}'] = _r(y[m].std()) if ok else None
    if t_ref is not None:
        for w in WINDOWS:
            m = _win(t, w)
            ok = m.any() and _full_window(t, w)
            out[f'rms_diff_{w}'] = _r(np.sqrt(np.mean((y[m] - np.interp(t[m], t_ref, y_ref)) ** 2))) if ok else None
    return out


def force_metrics(t, fx, fy, t_ref=None, fx_ref=None, fy_ref=None):
    out = {}
    for w in ('0deg', '15deg'):
        m = _win(t, w)
        ok = m.any() and _full_window(t, w)
        for name, f, fr in (('fx', fx, fx_ref), ('fy', fy, fy_ref)):
            out[f'{name}_mean_{w}'] = _r(f[m].mean()) if ok else None
            out[f'{name}_std_{w}'] = _r(f[m].std()) if ok else None
            if t_ref is not None:
                out[f'rms_diff_{name}_{w}'] = (_r(np.sqrt(np.mean((f[m] - np.interp(t[m], t_ref, fr)) ** 2)))
                                               if ok else None)
    return out


def reference_metrics(ref=None):
    '''The same measures for the reference (per unit depth).'''
    ref = ref if ref is not None else np.load(REF_FILE)
    t = ref['t']
    return dict(te=te_metrics(t, ref['te_disp'][:, 1]),
                force=force_metrics(t, ref['fx_stored'] / sl.SPAN, ref['fy_p'] / sl.SPAN))


def window_mean(t, values, name):
    m = _win(t, name)
    return values[m].mean(axis=0) if m.any() and _full_window(t, name) else None


# --------------------------------------------------------------------------
# One run
# --------------------------------------------------------------------------

def analyse(run_dir, ref_file=REF_FILE, tag=''):
    '''Comparison with the SoRoLily run in ref_file; tag (e.g. '_nu30')
    suffixes the files it writes (the flow snapshots only for the default reference).'''
    from scipy.signal import lfilter, sosfilt
    run_dir = Path(run_dir)
    viz = run_dir / 'viz_IB2d'
    inp = rr.read_input2d(run_dir / 'input2d')
    ref = np.load(ref_file)
    mesh = np.load(run_dir / 'wing_mesh.npz')
    wh = np.load(viz / 'wing_history.npz')
    t = wh['t']
    good = np.all(np.isfinite(wh['pos']), axis=(1, 2)) & np.all(np.isfinite(wh['force']), axis=1)
    n_ok = int(np.argmin(good)) if not good.all() else len(t)   # samples before the first non-finite one
    t, pos, F = t[:n_ok], wh['pos'][:n_ok].astype(float), wh['force'][:n_ok]
    xr, yr, te = mesh['xr'], mesh['yr'], int(mesh['te_index'])
    result = dict(t_end=_r(t[-1]) if len(t) else None,
                  re_chord=_r(inp['rho'] * sl.U * sl.CHORD / inp['mu']))
    # stability: the last step with finite tether forces (IB2d usually crashes
    # rather than writing NaN frames, so this is the reliable divergence time)
    if 'force_every_step' in wh.files:
        fe = wh['force_every_step']
        ok = np.all(np.isfinite(fe), axis=1)
        n_steps = int(round(inp['Tfinal'] / inp['dt']))
        bad = np.where(~ok[:n_steps])[0]
        info = json.loads((run_dir / 'run_info.json').read_text())
        crashed = bool(info['outcome'].get('exit_code'))
        first_bad = (bad[0] + 0.5) * inp['dt'] if len(bad) else None
        result['stable'] = (not crashed) and first_bad is None
        if crashed:
            # the hook saves the history when it sees the divergence; if it
            # missed it (runs before 2026-10-18 20:00: overflowed but finite
            # positions), the last good field frame bounds the time from below
            last_good = info['outcome'].get('last_output_time') or 0.0
            first_bad = max(first_bad or 0.0, last_good)
        result['first_nonfinite_time'] = _r(first_bad) if first_bad is not None else None
    if not len(t):
        return dict(result, note='no finite samples')

    tid = mesh['target_idx']
    result['clamp_drift_max'] = _r(np.hypot(pos[:, tid, 0] - xr[tid], pos[:, tid, 1] - yr[tid]).max())

    # trailing edge
    y_te = pos[:, te, 1] - yr[te]
    x_te = pos[:, te, 0] - xr[te]
    result['te'] = te_metrics(t, y_te, ref['t'], ref['te_disp'][:, 1])

    # forces: the reference's causal filter (Butterworth 6, 1.5 Hz), applied to
    # the force at every step (so the clamp's ~2 ms ringing is filtered, not
    # aliased into the 0.01 s samples), then sampled at the sample times
    if 'force_every_step' in wh.files:
        fe = wh['force_every_step']
        n_fe = int(np.argmin(np.all(np.isfinite(fe), axis=1))) if not np.all(np.isfinite(fe)) else len(fe)
        fe = fe[:n_fe]
        t_fe = (np.arange(n_fe) + 0.5) * inp['dt']
        sos = sl.butterworth(sl.FORCE_CUTOFF_HZ, fs=1.0 / inp['dt'], output='sos')
        fx_f = np.interp(t, t_fe, sosfilt(sos, fe[:, 0]))
        fy_f = np.interp(t, t_fe, sosfilt(sos, fe[:, 1]))
        result['force_ringing_rms'] = _r(np.sqrt(np.mean(np.sum((fe - np.column_stack(
            [sosfilt(sos, fe[:, 0]), sosfilt(sos, fe[:, 1])])) ** 2, 1))))
    else:
        b, a = sl.butterworth(sl.FORCE_CUTOFF_HZ, fs=1.0 / float(wh['interval']))
        fx_f, fy_f = lfilter(b, a, F[:, 0]), lfilter(b, a, F[:, 1])
    fx_ref, fy_ref = ref['fx_stored'] / sl.SPAN, ref['fy_p'] / sl.SPAN
    result['force'] = force_metrics(t, fx_f, fy_f, ref['t'], fx_ref, fy_ref)
    result['force'].update(fx_mean_15deg_unfiltered=_r(window_mean(t, F[:, 0], '15deg')),
                           fy_mean_15deg_unfiltered=_r(window_mean(t, F[:, 1], '15deg')))

    # centreline at the reference's stations (rest: straight centreline, world)
    s = ref['centreline_s']
    qx = mesh['x_LE'] + mesh['rigid_len'] + s
    qy = np.full(len(s), float(mesh['y_axis']))
    cmap = MaterialMap(xr, yr, mesh['triangles'], qx[:-1], qy[:-1])
    cl = np.concatenate([cmap(pos), pos[:, te:te + 1, :]], axis=1)       # last station = the TE vertex
    shape = {}
    for key, tt in SHAPE_TIMES.items():
        if t[-1] >= tt:
            k, kr = np.argmin(np.abs(t - tt)), np.argmin(np.abs(ref['t'] - tt))
            shape[f'rms_dist_{key}'] = _r(np.sqrt(np.mean(np.sum((cl[k] - ref['centreline'][kr]) ** 2, 1))))
        else:
            shape[f'rms_dist_{key}'] = None
    for w in ('0deg', '15deg'):
        a_ib, a_ref = window_mean(t, cl, w), window_mean(ref['t'], ref['centreline'], w)
        shape[f'rms_dist_mean_{w}'] = (_r(np.sqrt(np.mean(np.sum((a_ib - a_ref) ** 2, 1))))
                                       if a_ib is not None else None)
    result['shape'] = shape

    # surface pressure at the reference's panel midpoints, every field frame
    pn0 = sl.panels(np.zeros(4))
    pmap = MaterialMap(xr, yr, mesh['triangles'], pn0['mid'][:, 0], pn0['mid'][:, 1])
    # IB2d's pressure jumps across the interface, smeared over the delta
    # function's support (supp/2 cells each side): on the surface it mixes the
    # inside (which balances the structure's own stresses: +-800 Pa) with the
    # outside. The reference's BDIM surface value is the outside pressure, so
    # sample just outside the smearing, supp/2 * dx along the outward normal.
    n_out = int(mesh['n_outline'])
    p_offset = 0.5 * inp['supp'] * inp['Lx'] / inp['Nx']
    frames = sorted(viz.glob('P.*.vtk'), key=lambda p: int(p.name.split('.')[1]))
    t_frame = inp['print_dump'] * inp['dt']
    tp, pp = [], []
    for f in frames:
        k = int(f.name.split('.')[1])
        tk = k * t_frame
        if tk > t[-1] + 1e-9:
            break
        field = read_scalar_vtk(f)
        if not np.all(np.isfinite(field)):
            break
        i = np.argmin(np.abs(t - tk))
        pts = pmap(pos[i])
        pts = pts + p_offset * outward_normals(pos[i, :n_out], pts)
        p = sample_periodic(field, inp['Lx'], inp['Ly'], pts[:, 0], pts[:, 1])
        p -= sample_periodic(field, inp['Lx'], inp['Ly'], [GAUGE_POINT[0]], [GAUGE_POINT[1]])[0]
        tp.append(tk); pp.append(p)
    tp, pp = np.array(tp), np.array(pp)
    pres = dict(n_frames=len(tp), sample_offset_m=_r(p_offset))
    for w in ('0deg', '15deg'):
        m_ib = window_mean(tp, pp, w) if len(tp) else None
        m_ref = window_mean(ref['t'], ref['p_panels'], w)
        pres[f'rms_diff_mean_{w}'] = _r(np.sqrt(np.mean((m_ib - m_ref) ** 2))) if m_ib is not None else None
        pres[f'ref_rms_mean_{w}'] = _r(np.sqrt(np.mean(m_ref ** 2)))
    result['pressure'] = pres

    # fluid
    fluid = {}
    fh_path = viz / 'fluid_history.npz'
    if fh_path.exists():
        fh = np.load(fh_path)
        tf, ke, kg, um, umean = fh['t'], fh['ke'], fh['ke_gs'], fh['u_max'], fh['u_mean']
        ok = np.isfinite(ke) & np.isfinite(kg) & np.isfinite(um)
        tf, ke, kg, um, umean = tf[ok], ke[ok], kg[ok], um[ok], umean[ok]
        dxg = inp['Lx'] / inp['Nx']
        if len(tf):
            fluid['u_mean_end'] = _r(np.hypot(*umean[-1]))
            fluid['u_max'] = _r(um.max())
            fluid['cell_re_max'] = _r(um.max() * dxg * inp['rho'] / inp['mu'])
            m = (tf >= 1.0) & (kg > 0)
            fluid['gs_ke_growth_rate'] = _r(np.polyfit(tf[m], np.log(kg[m]), 1)[0]) if m.sum() >= 3 else None
            last = tf >= tf[-1] - 1.0
            fluid['gs_ke_fraction_end'] = _r(kg[last].sum() / ke[last].sum()) if ke[last].sum() > 0 else None
    result['fluid'] = fluid

    # time series + figures
    with open(run_dir / f'timeseries{tag}.csv', 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['t', 'te_dx', 'te_dy', 'te_dy_ref', 'fx', 'fy', 'fx_filtered', 'fy_filtered',
                     'fx_ref', 'fy_ref'])
        for k in range(len(t)):
            wr.writerow(['%.6g' % v for v in (
                t[k], x_te[k], y_te[k], np.interp(t[k], ref['t'], ref['te_disp'][:, 1]), F[k, 0], F[k, 1],
                fx_f[k], fy_f[k], np.interp(t[k], ref['t'], fx_ref), np.interp(t[k], ref['t'], fy_ref))])
    try:
        plot_comparison(run_dir, ref, t, y_te, F, fx_f, fy_f, cl, tp, pp, pn0, fh_path, inp, tag)
        if not tag:
            plot_flow(run_dir, inp, t, pos, mesh)
    except Exception as e:              # a figure must never lose the analysis
        result['plot_error'] = repr(e)
    return result


def plot_comparison(run_dir, ref, t, y_te, F, fx_f, fy_f, cl, tp, pp, pn0, fh_path, inp, tag=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 2, figsize=(14, 12))
    rt = ref['t']
    a = ax[0, 0]
    a.plot(rt, 1e3 * ref['te_disp'][:, 1], 'k', label='SoRoLily')
    a.plot(t, 1e3 * y_te, 'C3', label='IB2d')
    a.set_ylabel('TE y displacement [mm]'); a.set_xlabel('t [s]'); a.legend(); a.set_title('trailing edge')
    a = ax[0, 1]
    a.plot(rt, 1e3 * ref['te_disp'][:, 1], 'k'); a.plot(t, 1e3 * y_te, 'C3')
    a.set_xlim(0, 2.5); a.set_xlabel('t [s]'); a.set_title('release (zoom)')
    a = ax[1, 0]
    a.plot(rt, ref['fx_stored'] / sl.SPAN, 'k', label='Fx SoRoLily')
    a.plot(rt, ref['fy_p'] / sl.SPAN, 'k--', label='Fy SoRoLily')
    a.plot(t, F[:, 0], 'C0', lw=0.4, alpha=0.4); a.plot(t, F[:, 1], 'C3', lw=0.4, alpha=0.4)
    a.plot(t, fx_f, 'C0', label='Fx IB2d (filtered)'); a.plot(t, fy_f, 'C3', label='Fy IB2d (filtered)')
    lim = np.nanpercentile(np.abs(np.concatenate([fx_f, fy_f, ref['fy_p'] / sl.SPAN])), 99.5) * 1.3
    a.set_ylim(-lim, lim); a.set_ylabel('force per unit depth [N/m]'); a.set_xlabel('t [s]')
    a.legend(fontsize=8); a.set_title('hydrodynamic force (faint: IB2d unfiltered)')
    a = ax[1, 1]
    for k, tt in enumerate((0.0, 0.25, 0.5, 1.0)):
        if t[-1] >= tt:
            i, ir = np.argmin(np.abs(t - tt)), np.argmin(np.abs(rt - tt))
            a.plot(ref['centreline'][ir, :, 0], 1e3 * (ref['centreline'][ir, :, 1] - sl.PIVOT[1]), 'k.-', ms=3)
            a.plot(cl[i, :, 0], 1e3 * (cl[i, :, 1] - sl.PIVOT[1]), '.-', color=f'C{k}', ms=3, label=f't = {tt:g} s')
    a.set_xlabel('x [m]'); a.set_ylabel('centreline y [mm]'); a.legend(fontsize=8)
    a.set_title('centreline (black: SoRoLily)')
    a = ax[2, 0]
    x_rest = pn0['x_rest']
    up = pn0['side'] == 'upper'
    for w, c in (('0deg', 'C0'), ('15deg', 'C3')):
        m_ref = window_mean(rt, ref['p_panels'], w)
        m_ib = window_mean(tp, pp, w) if len(tp) else None
        for side, ls in ((up, '-'), (~up, '--')):
            a.plot(x_rest[side], m_ref[side], 'k' + ls, lw=0.8)
            if m_ib is not None:
                a.plot(x_rest[side], m_ib[side], c + ls, label=f'IB2d {w} ' + ('upper' if side is up else 'lower'))
    a.set_xlabel("x' from LE [m]"); a.set_ylabel('window-mean surface pressure [Pa]'); a.legend(fontsize=7)
    a.set_title('surface pressure (black: SoRoLily; solid upper, dashed lower)')
    a = ax[2, 1]
    if fh_path.exists():
        fh = np.load(fh_path)
        a.semilogy(fh['t'], fh['ke'], label='KE'); a.semilogy(fh['t'], fh['ke_gs'], label='KE grid-scale')
        a2 = a.twinx(); a2.plot(fh['t'], np.hypot(*fh['u_mean'].T), 'k', lw=0.8); a2.set_ylabel('|<u>| [m/s]')
        a.legend(fontsize=8)
    a.set_xlabel('t [s]'); a.set_ylabel('fluid KE [J/m]'); a.set_title('fluid')
    fig.suptitle('%s   (Nx %d, dt %g, mu %.3g x water)%s' % (run_dir.name, inp['Nx'], inp['dt'],
                                                         inp['mu'] / (sl.RHO * sl.NU),
                                                         f'   vs SoRoLily run {tag[1:]}' if tag else ''))
    fig.tight_layout(); fig.savefig(run_dir / f'comparison{tag}.png', dpi=110); plt.close(fig)


def plot_flow(run_dir, inp, t, pos, mesh):
    '''Vorticity snapshots (IB2d only: the reference stored no fields).'''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    viz = run_dir / 'viz_IB2d'
    frames = {int(p.name.split('.')[1]): p for p in viz.glob('Omega.*.vtk')}
    if not frames:
        return
    t_frame = inp['print_dump'] * inp['dt']
    want = [tt for tt in (0.5, 2.0, 4.95, 6.0, 8.0, 10.0) if tt <= t[-1] + 1e-9]
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    nout = int(mesh['n_outline'])
    for a, tt in zip(axs.ravel(), want):
        k = min(frames, key=lambda kk: abs(kk * t_frame - tt))
        w = read_scalar_vtk(frames[k])
        if not np.all(np.isfinite(w)):
            continue
        ext = (0, inp['Lx'], 0, inp['Ly'])
        a.imshow(np.clip(w, -10, 10), origin='lower', extent=ext, cmap='RdBu_r', vmin=-10, vmax=10)
        i = np.argmin(np.abs(t - k * t_frame))
        o = np.vstack([pos[i, :nout], pos[i, :1]])
        a.plot(o[:, 0], o[:, 1], 'k-', lw=0.8)
        a.set_xlim(0.2, 1.2); a.set_ylim(0.3, 0.9); a.set_aspect('equal')
        a.set_title('vorticity, t = %.2f s' % (k * t_frame))
    for a in axs.ravel()[len(want):]:
        a.axis('off')
    fig.tight_layout(); fig.savefig(run_dir / 'flow_snapshots.png', dpi=100); plt.close(fig)


def main(run_dirs, ref_file=None, name=None):
    '''Default: the water reference -> run_info['analysis']. With --ref FILE
    --name NAME: that SoRoLily run -> run_info['comparisons'][NAME], files
    suffixed _NAME; 'analysis' is left as it is.'''
    for d in run_dirs:
        d = Path(d)
        if not (d / 'viz_IB2d' / 'wing_history.npz').exists():
            print(f'{d}: no viz_IB2d/wing_history.npz, skipped')
            continue
        info_path = d / 'run_info.json'
        info = json.loads(info_path.read_text())
        if ref_file:
            a = analyse(d, Path(ref_file), f'_{name}')
            a['reference_file'] = str(Path(ref_file).resolve())
            info.setdefault('comparisons', {})[name] = a
        else:
            a = info['analysis'] = analyse(d)
        rr.write_run_info(d, info)
        print(f"{d.name}: t_end {a.get('t_end')}, TE rms diff (release) {a.get('te', {}).get('rms_diff_release')}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    args = sys.argv[1:]
    ref_file = name = None
    if '--ref' in args:
        i = args.index('--ref'); ref_file = args[i + 1]; del args[i:i + 2]
    if '--name' in args:
        i = args.index('--name'); name = args[i + 1]; del args[i:i + 2]
    if bool(ref_file) != bool(name):
        sys.exit('--ref and --name go together')
    main(args, ref_file, name)
