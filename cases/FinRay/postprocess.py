'''
Post-processing of FinRay runs: adds an "analysis" section to each run's
run_info.json, and writes two time series into the run folder.

    python postprocess.py <run folder> [<run folder> ...]

run_case.py calls this after every run. It only reads the run folder
(viz_IB2d/, case_params.json, input2d, run_info.json), so it can be re-run on
old runs without re-simulating; rebuild the index afterwards
(`python run_case.py index [--all] DIR`).

Tip (analysis.tip; tip_trajectory.csv)
    The tip is the outer tip vertex: the point with the largest x at rest
    (frame 0). Its displacement from rest (dx, dy) is measured over the
    measured cycles: run_length.measured_cycles whole cycles after
    run_length.transient_cycles (one cycle = 1/FREQUENCY, both sides push and
    return). If the run is shorter (or went NaN), the whole cycles that fit
    after the transient are used, and n_cycles says how many.
      y_amplitude     mean over cycles of (max dy - min dy) / 2 [m]
      y_amplitude_cv  std / mean of that per-cycle amplitude (periodicity check)
      y_mean          mean dy over the measured cycles [m]
      x_amplitude     the same as y_amplitude, for dx [m]
      phase_lag_deg   lag of dy behind the actuation, at the actuation
                      frequency (first Fourier harmonic over the measured
                      cycles). Actuation reference: d_bottom - d_top, the
                      bottom patch's push minus the top patch's push, which
                      moves the tip in +y when the structure follows it
                      quasi-statically, so ~0 = in phase, > 0 = the tip lags.
    tip_trajectory.csv holds t, x, y, dx, dy and the reference for every frame,
    for comparing whole trajectories between runs.

Fluid (analysis.fluid; fluid_energy.csv)
    Per frame, from uX/uY: kinetic energy per unit depth
    KE = rho/2 * sum(u^2 + v^2) * dx^2, and its grid-scale part KE_gs, the
    same sum over the residual u - (3x3 box average of u) (periodic grid).
      gs_ke_growth_rate   slope of ln(KE_gs) against t [1/s], fitted over the
                          measured cycles (or from STIFFNESS_RAMP_TIME on if
                          no whole cycle fits); > 0 means grid-scale noise grows
      gs_ke_fraction_end  KE_gs / KE over the last cycle of that window
      u_max               peak fluid speed over the run (before any NaN) [m/s]
      cell_re_max         u_max * dx * rho / mu
'''

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import run_records as rr


def _floats_after(text, header, n):
    '''The n numbers that follow the line starting with `header` in a VTK file.'''
    i = text.index("\n", text.index(header)) + 1
    vals = np.array(text[i:].split(maxsplit=n)[:n], dtype=float)
    if len(vals) != n:
        raise ValueError(f"expected {n} values after {header!r}, found {len(vals)}")
    return vals


def read_lag_points(path):
    '''(n, 2) Lagrangian point positions from a lagPtsConnect.*.vtk frame.'''
    text = Path(path).read_text()
    i = text.index("POINTS ")
    n = int(text[i:].split()[1])
    return _floats_after(text, "POINTS ", 3 * n).reshape(n, 3)[:, :2]


def read_grid_scalar(path):
    '''(Ny, Nx) scalar field from an IB2d STRUCTURED_POINTS .vtk frame.'''
    text = Path(path).read_text()
    i = text.index("DIMENSIONS ")
    nx, ny = (int(v) for v in text[i:].split()[1:3])
    return _floats_after(text, "LOOKUP_TABLE", nx * ny).reshape(ny, nx)


def actuation(t, frequency, amplitude):
    '''(d_top, d_bottom): the two patches' pushes in +x, as in
    update_Target_Point_Positions.py (top first, sin^2 bump, then bottom).'''
    phase = (2 * np.pi * frequency * np.asarray(t)) % (2 * np.pi)
    first = phase < np.pi
    d_top = np.where(first, amplitude * np.sin(phase) ** 2, 0.0)
    d_bottom = np.where(first, 0.0, amplitude * np.sin(phase - np.pi) ** 2)
    return d_top, d_bottom


def _box3(u):
    '''3x3 box average on a periodic grid.'''
    s = sum(np.roll(np.roll(u, i, axis=0), j, axis=1) for i in (-1, 0, 1) for j in (-1, 0, 1))
    return s / 9.0


def _frames(viz, prefix):
    '''{frame number: path} of viz/<prefix>.NNNN.vtk.'''
    return {int(p.name.split(".")[1]): p for p in viz.glob(f"{prefix}.*.vtk")}


def _first_harmonic(signal, t, frequency):
    '''Complex Fourier coefficient at `frequency` (whole cycles, uniform samples).'''
    return np.mean(signal * np.exp(-2j * np.pi * frequency * t))


def _r(x, digits=6):
    return None if x is None or not np.isfinite(x) else float(f"{x:.{digits}g}")


def analyse(run_dir):
    run_dir = Path(run_dir)
    viz = run_dir / "viz_IB2d"
    params = json.loads((run_dir / "case_params.json").read_text())
    inp = rr.read_input2d(run_dir / "input2d")
    act = params.get("actuation", {})
    frequency, amplitude = act.get("FREQUENCY", 1.0), act.get("AMPLITUDE", 0.01)
    ramp = act.get("STIFFNESS_RAMP_TIME", 0.2)
    rl = params.get("run_length", {})
    n_transient, n_measured = rl.get("transient_cycles", 1), rl.get("measured_cycles", 3)
    t_frame = inp["print_dump"] * inp["dt"]
    period = 1.0 / frequency

    # --- tip: every frame up to the first NaN
    lag = _frames(viz, "lagPtsConnect")
    if 0 not in lag:
        raise RuntimeError(f"{viz}: no frame 0 (rest positions)")
    rest = read_lag_points(lag[0])
    tip = int(np.argmax(rest[:, 0]))
    ts, xy = [], []
    for k in sorted(lag):
        p = read_lag_points(lag[k])[tip]
        if not np.all(np.isfinite(p)):
            break
        ts.append(k * t_frame)
        xy.append(p)
    t = np.round(np.array(ts), 9)
    xy = np.array(xy)
    dx, dy = xy[:, 0] - rest[tip, 0], xy[:, 1] - rest[tip, 1]
    d_top, d_bottom = actuation(t, frequency, amplitude)
    ref = d_bottom - d_top

    with open(run_dir / "tip_trajectory.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "x", "y", "dx", "dy", "actuation_ref"])
        for row in zip(t, xy[:, 0], xy[:, 1], dx, dy, ref):
            w.writerow([f"{v:.9g}" for v in row])

    # whole measured cycles that fit in the data
    t0 = n_transient * period
    n_fit = int(np.floor((t[-1] - t0) / period + 1e-6)) if t[-1] > t0 else 0
    n_cycles = min(n_measured, n_fit)
    tip_info = dict(index=tip, rest_x=float(rest[tip, 0]), rest_y=float(rest[tip, 1]),
                    n_cycles=n_cycles, window=None, y_amplitude=None, y_amplitude_cv=None,
                    y_mean=None, x_amplitude=None, phase_lag_deg=None)
    if n_cycles:
        t1 = t0 + n_cycles * period
        amp_y, amp_x = [], []
        for c in range(n_cycles):
            m = (t >= t0 + c * period - 1e-9) & (t <= t0 + (c + 1) * period + 1e-9)
            amp_y.append((dy[m].max() - dy[m].min()) / 2)
            amp_x.append((dx[m].max() - dx[m].min()) / 2)
        win = (t >= t0 - 1e-9) & (t < t1 - 1e-9)   # whole cycles, end point excluded
        cy = _first_harmonic(dy[win], t[win], frequency)
        cr = _first_harmonic(ref[win], t[win], frequency)
        lag_deg = np.degrees(np.angle(cr * np.conj(cy)))
        tip_info.update(window=[round(t0, 9), round(t1, 9)],
                        y_amplitude=_r(np.mean(amp_y)),
                        y_amplitude_cv=_r(np.std(amp_y) / np.mean(amp_y)) if np.mean(amp_y) else None,
                        y_mean=_r(dy[win].mean()), x_amplitude=_r(np.mean(amp_x)),
                        phase_lag_deg=_r(lag_deg, 5))

    # --- fluid: kinetic energy and its grid-scale part, up to the first NaN
    ux, uy = _frames(viz, "uX"), _frames(viz, "uY")
    rho, mu = inp["rho"], inp["mu"]
    dxg = inp["Lx"] / inp["Nx"]
    rows = []
    for k in sorted(set(ux) & set(uy)):
        u, v = read_grid_scalar(ux[k]), read_grid_scalar(uy[k])
        if not (np.all(np.isfinite(u)) and np.all(np.isfinite(v))):
            break
        ke = 0.5 * rho * np.sum(u ** 2 + v ** 2) * dxg ** 2
        ru, rv = u - _box3(u), v - _box3(v)
        ke_gs = 0.5 * rho * np.sum(ru ** 2 + rv ** 2) * dxg ** 2
        rows.append((round(k * t_frame, 9), ke, ke_gs, float(np.sqrt(u ** 2 + v ** 2).max())))
    with open(run_dir / "fluid_energy.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "ke", "ke_gs", "u_max"])
        for row in rows:
            w.writerow([f"{v:.9g}" for v in row])

    fluid = dict(gs_ke_growth_rate=None, gs_ke_fraction_end=None, fit_window=None,
                 u_max=None, cell_re_max=None)
    if rows:
        tf, ke, ke_gs, umax = (np.array(c) for c in zip(*rows))
        fluid["u_max"] = _r(umax.max())
        fluid["cell_re_max"] = _r(umax.max() * dxg * rho / mu)
        a, b = (t0, t0 + n_cycles * period) if n_cycles else (ramp, tf[-1])
        m = (tf >= a - 1e-9) & (tf <= b + 1e-9) & (ke_gs > 0)
        if m.sum() >= 3:
            fluid["gs_ke_growth_rate"] = _r(np.polyfit(tf[m], np.log(ke_gs[m]), 1)[0])
            fluid["fit_window"] = [round(float(tf[m][0]), 9), round(float(tf[m][-1]), 9)]
            last = m & (tf >= tf[m][-1] - period + 1e-9)
            if ke[last].sum() > 0:
                fluid["gs_ke_fraction_end"] = _r(ke_gs[last].sum() / ke[last].sum())

    return dict(tip=tip_info, fluid=fluid)


def main(run_dirs):
    failed = 0
    for d in run_dirs:
        d = Path(d)
        try:
            analysis = analyse(d)
        except Exception as e:   # one bad run shouldn't stop a batch
            print(f"{d}: FAILED ({type(e).__name__}: {e})")
            failed += 1
            continue
        info_path = d / "run_info.json"
        if info_path.exists():
            info = json.loads(info_path.read_text())
            info["analysis"] = analysis
            rr.write_run_info(d, info)
        else:
            print(f"{d}: no run_info.json (backfill it first); analysis printed only")
        tip, fl = analysis["tip"], analysis["fluid"]
        print(f"{d.name}: {tip['n_cycles']} cycle(s), tip y amplitude {tip['y_amplitude']} m, "
              f"lag {tip['phase_lag_deg']} deg; grid-scale KE growth {fl['gs_ke_growth_rate']} /s, "
              f"cell Re max {fl['cell_re_max']}")
    return failed


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(1 if main(sys.argv[1:]) else 0)
