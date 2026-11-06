'''
Sweep definitions for the NACA0020 wing's run_case.py (same conventions as
cases/FinRay/sweeps.py): one function per sweep, chosen by name.

    python run_case.py --sweep <name> --list         # the cases, nothing run
    python run_case.py --sweep <name> --case N       # run case N
    slurm/submit_sweep.sh <sweep_folder> <name>      # every case on the cluster

A sweep function takes the base case (case_params.json, a fresh deep copy) and
returns a list of cases. Parameters are named by their path in case_params,
e.g. "geom.target_stiffness_factor", "input2d.dt", "input2d.mu". `grid()`
builds the cartesian product and names each case after the values it varies.
Register a sweep with @sweep. The stages refer to IMPLEMENTATION_DECISIONS.md §4.

Derived values (finalize, every case):
  * input2d.print_dump = output.field_interval / dt (the .vtk field output);
    the hooks' diagnostics sample every diagnostics.interval regardless.
    Both must be whole numbers of steps: pick dt with 1/dt a multiple of 100.
  * input2d.springs = 0 when geom.spring_force = "hook" (springs computed,
    vectorised, in the user-force hook), else 1.
  * Nx = Ny always (dx = dy); vary them together with ("input2d.Nx", "input2d.Ny").
  * geom.stiffness_scale = null means "measured for this mesh" (Wing_Geom.py),
    so sweeps over Nx or the density never carry a stale calibration.
'''

import copy
import itertools
import math

SWEEPS = {}

WATER_MU = 1.0023e-3      # Pa s: rho * nu of the reference (nu = 1.0023e-6, rho = 1000)

ALIASES = {
    "geom.fill_density_multiplier": "mult",
    "geom.stiffness_scale": "scale",
    "geom.target_stiffness_factor": "ktf",
    "input2d.Tfinal": "T",
    "flow.fringe_rate": "fringe",
    "flow.start_ramp": "ramp",
}

# dt ladders per Nx (1/dt a multiple of 100). Nx = 192 brackets the archived
# 10 s run (dt 2e-4 completed); the others are shifted by about 2x per grid
# doubling (the scaling is measured, not assumed).
DT_LADDER = {
    96: [1.25e-3, 1e-3, 6.25e-4, 5e-4, 4e-4, 3.125e-4, 2.5e-4, 2e-4],
    192: [6.25e-4, 5e-4, 4e-4, 3.125e-4, 2.5e-4, 2e-4, 1.5625e-4, 1.25e-4, 1e-4, 7.8125e-5, 6.25e-5],
    384: [3.125e-4, 2.5e-4, 2e-4, 1.5625e-4, 1.25e-4, 1e-4, 7.8125e-5, 6.25e-5],
}


def sweep(fn):
    SWEEPS[fn.__name__] = fn
    return fn


def set_(case, path, value):
    *parents, key = path.split(".")
    d = case
    for k in parents:
        d = d.setdefault(k, {})
    if key not in d and parents and parents[0] in ("geom", "flow", "input2d", "output", "diagnostics"):
        known = ", ".join(sorted(d))
        raise KeyError(f"{path}: not in the base case (known keys: {known}); add it to "
                       f"case_params.json first if it is new")
    d[key] = value


def _label(path):
    return ALIASES.get(path, path.split(".")[-1])


def _fmt(path, value):
    if path == "input2d.mu":
        return f"mu{value / WATER_MU:g}w"            # multiples of water
    if isinstance(value, bool):
        return f"{_label(path)}{int(value)}"
    if isinstance(value, (int, float)):
        return f"{_label(path)}{value:g}"
    return f"{_label(path)}{value}"


def grid(base, axes):
    '''Cartesian product over `axes` ({path or tuple of paths: [values]}).'''
    keys = list(axes)
    cases = []
    for values in itertools.product(*(axes[k] for k in keys)):
        case = copy.deepcopy(base)
        parts = []
        for key, value in zip(keys, values):
            paths = key if isinstance(key, tuple) else (key,)
            for p in paths:
                set_(case, p, value)
            parts.append(_fmt(paths[0], value))
        case["case_name"] = "_".join(parts) or "base"
        cases.append(case)
    return cases


def _whole_steps(seconds, dt, what):
    n = seconds / dt
    if abs(n - round(n)) > 1e-6 or round(n) < 1:
        raise ValueError(f"{what} = {seconds} s is not a whole number of steps at dt = {dt}")
    return int(round(n))


def finalize(case, sweep_name):
    inp = case["input2d"]
    # geom.spring_force = "hook": the vectorised springs in the user-force hook,
    # so IB2d's own (looped) spring force must be off, and vice versa
    inp["springs"] = 0 if case["geom"].get("spring_force", "ib2d") == "hook" else 1
    inp["print_dump"] = _whole_steps(case["output"]["field_interval"], inp["dt"], "output.field_interval")
    _whole_steps(case["diagnostics"]["interval"], inp["dt"], "diagnostics.interval")
    _whole_steps(inp["Tfinal"], inp["dt"], "input2d.Tfinal")
    if not math.isclose(inp["Lx"] / inp["Nx"], inp["Ly"] / inp["Ny"]):
        raise ValueError(f"{case['case_name']}: dx != dy; vary Nx and Ny together")
    case["sweep"] = sweep_name
    return case


def build(name, base):
    if name not in SWEEPS:
        raise SystemExit(f"unknown sweep '{name}'; defined: {', '.join(SWEEPS)}")
    cases = [finalize(c, name) for c in SWEEPS[name](copy.deepcopy(base))]
    names = [c["case_name"] for c in cases]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise ValueError(f"sweep '{name}': duplicate case names {dupes}")
    return cases


# --- Sweeps -------------------------------------------------------------------

@sweep
def base(b):
    '''The base case alone (case_params.json as it is).'''
    return grid(b, {})


@sweep
def smoke(b):
    '''The base case for 0.1 s: the code runs end to end (geometry, hooks, post-processing).'''
    b["input2d"]["Tfinal"] = 0.1
    cases = grid(b, {})
    cases[0]["case_name"] = "smoke"
    return cases


def _dt_cases(b, Nx, factors, mus):
    out = []
    for c in grid(b, {"geom.target_stiffness_factor": factors,
                      "input2d.mu": [m * WATER_MU for m in mus],
                      "input2d.dt": DT_LADDER[Nx]}):
        c["input2d"]["Nx"] = c["input2d"]["Ny"] = Nx
        c["case_name"] = f"Nx{Nx}_" + c["case_name"]
        out.append(c)
    return out


@sweep
def stage3_dt(b):
    '''
    Stage 3, critical dt: 1 s runs (the start-up ramp and the release of the
    pre-bend, the most violent part of the schedule); pass = no NaN.
      Nx 192: target_stiffness_factor 3 and 15 (O4) x mu 30x and 1x water
      Nx 96, 384: both factors at mu 30x (dt_crit per grid for stage 4)
    '''
    b["input2d"]["Tfinal"] = 1.0
    return (_dt_cases(b, 192, [3.0, 15.0], [30, 1])
            + _dt_cases(b, 96, [3.0, 15.0], [30])
            + _dt_cases(b, 384, [3.0, 15.0], [30]))


# Production time steps (stage 3, sweep stage3_dt, target factor 3): the largest
# ladder value <= 0.7 dt_crit. dt_crit (1 s runs, mu 30x; at Nx 192 the same
# with water viscosity) = 1e-3, 3.125e-4, 1.5625e-4 at Nx 96, 192, 384.
DT_PRODUCTION = {96: 6.25e-4, 192: 2e-4, 384: 1e-4}
MU_LEVELS = [1, 3, 10, 30]          # x water (O5 default)


def _ladder_case(b, Nx, mu):
    c = copy.deepcopy(b)
    c["input2d"]["Nx"] = c["input2d"]["Ny"] = Nx
    c["input2d"]["dt"] = DT_PRODUCTION[Nx]
    c["input2d"]["mu"] = mu * WATER_MU
    c["geom"]["target_stiffness_factor"] = 3.0
    c["case_name"] = f"Nx{Nx}_mu{mu:g}w"
    return c


@sweep
def stage4_ladder(b):
    '''Stage 4: Nx x mu ladder, full 10 s schedule, dt = DT_PRODUCTION[Nx].'''
    return [_ladder_case(b, Nx, mu) for Nx in sorted(DT_PRODUCTION) for mu in MU_LEVELS]


@sweep
def stage4_diag(b):
    '''
    Stage 4 diagnosis (2 s runs): why runs below 30x water diverge. Halving dt
    tests a time-step (IB coupling) limit; init_q = 0 (no pre-bend, so no
    violent release) tests whether the release motion triggers the near-wall
    grid-scale instability.
    '''
    out = []
    for Nx, mu, dt, q, tag in ((192, 10, 1e-4, 1.0, "dt_half"), (384, 10, 5e-5, 1.0, "dt_half"),
                               (192, 10, None, 0.0, "nobend"), (384, 10, None, 0.0, "nobend"),
                               (192, 1, None, 0.0, "nobend")):
        c = _ladder_case(b, Nx, mu)
        if dt:
            c["input2d"]["dt"] = dt
        c["geom"]["init_q"] = [q, 0.0, 0.0, 0.0]
        c["input2d"]["Tfinal"] = 2.0
        c["case_name"] += "_" + tag
        out.append(c)
    return out


@sweep
def stage5_sensitivity(b):
    '''
    Stage 5 at the lowest stable setting found in stage 4 (Nx 192, 30x water,
    dt 2e-4), full 10 s. Planned: fringe/mean rates x0.5 and x2 (O3). Added
    after stage 4 showed a sustained ~22 mm TE offset at 0 deg (SoRoLily ~0):
    no pre-bend (spontaneous or a remnant of the release?), 100x water (Re),
    no start ramp (start-up), target factor 15 (clamp), density x8 (mesh).
    '''
    base_ = _ladder_case(b, 192, 30)
    variants = [
        ("fringe_x0.5", {"flow.fringe_rate": 20.0, "flow.mean_rate": 10.0}),
        ("fringe_x2", {"flow.fringe_rate": 80.0, "flow.mean_rate": 40.0}),
        ("nobend", {"geom.init_q": [0.0, 0.0, 0.0, 0.0]}),
        ("mu100w", {"input2d.mu": 100 * WATER_MU}),
        ("noramp", {"flow.start_ramp": 0.0}),
        ("ktf15_dt1e-4", {"geom.target_stiffness_factor": 15.0, "input2d.dt": 1e-4}),
        ("mult8", {"geom.fill_density_multiplier": 8.0}),
    ]
    out = []
    for tag, changes in variants:
        c = copy.deepcopy(base_)
        for path, v in changes.items():
            set_(c, path, v)
        c["case_name"] = "Nx192_mu100w" if tag == "mu100w" else base_["case_name"] + "_" + tag
        out.append(c)
    return out


@sweep
def stage4_refine(b):
    '''
    Stage 4 refinement: the ladder showed a threshold in cell Re = U dx / nu
    between 42 (Nx 192, 30x: stable) and 62 (Nx 384, 10x: diverged), not a dt
    limit (stage4_diag). These bracket it and find the lowest stable mu:
    Nx 192 at 20x (cell Re 62), Nx 384 at 15x (42) and 20x (31).
    '''
    return [_ladder_case(b, 192, 20), _ladder_case(b, 384, 15), _ladder_case(b, 384, 20)]


@sweep
def stage4_mesh(b):
    '''
    Stage 4, structural mesh vs FSI stability (added 2026-10-18): the
    instability below the stable viscosities is a grid-scale (odd-even)
    pressure/vorticity checkerboard fed by the IB forces. At multiplier x4 the
    Lagrangian spacing is h = dx/4, half IB2d's own ds = dx/2; at x1, h = ds.
    Full 10 s runs at x1 (stiffness_scale measured for that mesh).
    '''
    out = []
    for Nx, mu, q in ((192, 30, 1.0), (192, 30, 0.0), (192, 10, 1.0), (192, 3, 1.0), (192, 1, 1.0),
                      (384, 10, 1.0), (384, 3, 1.0), (384, 1, 1.0)):
        c = _ladder_case(b, Nx, mu)
        c["geom"]["fill_density_multiplier"] = 1.0
        c["geom"]["init_q"] = [q, 0.0, 0.0, 0.0]
        c["case_name"] = f"Nx{Nx}_mu{mu:g}w_mult1" + ("" if q else "_nobend")
        out.append(c)
    return out


@sweep
def stage5_bend_sign(b):
    '''
    Added 2026-10-18: the sustained 0 deg TE offset (+22 mm, robust to fringe,
    start-up, Re 400-1330) -- hysteresis of the undamped release (then it flips
    with the pre-bend's sign) or a systematic asymmetry (then it doesn't)?
    Pre-bend q1 = -1 and a half-size pre-bend q1 = +0.5, Nx 192, 30x water.
    '''
    out = []
    for q, tag in ((-1.0, "bend_neg"), (0.5, "bend_half")):
        c = _ladder_case(b, 192, 30)
        c["geom"]["init_q"] = [q, 0.0, 0.0, 0.0]
        c["case_name"] += "_" + tag
        out.append(c)
    return out
