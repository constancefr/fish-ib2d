'''
Sweep definitions for run_case.py: one function per sweep, chosen by name.

    python run_case.py --sweep <name> --list         # the cases, nothing run
    python run_case.py --sweep <name> --case N       # run case N
    slurm/submit_sweep.sh <sweep_folder> <name>      # every case on the cluster

A sweep function takes the base case (case_params.json, a fresh deep copy) and
returns a list of cases. Parameters are named by their path in case_params,
e.g. "geom.ribbon_fill_density_multiplier", "input2d.dt", "actuation.FREQUENCY".
`grid()` builds the cartesian product of some parameters and names each case
after the values it varies; a sweep can then change cases one by one (e.g. a dt
per Nx). Register a sweep with @sweep.

The regime (actuation.FREQUENCY, actuation.AMPLITUDE) and the run length
(run_length.transient_cycles + measured_cycles) are read from the base case,
so editing case_params.json and resubmitting a sweep reruns it consistently in
the new regime. input2d.Tfinal = null (the default) means Tfinal is derived:
(transient_cycles + measured_cycles) / FREQUENCY, rounded up to a whole number
of steps. A sweep that needs a fixed length (e.g. 1 s dt-stability runs) sets
input2d.Tfinal itself.
'''

import copy
import itertools
import math

SWEEPS = {}

# Short names for case names (anything not listed uses its last path component)
ALIASES = {
    "geom.ribbon_fill_density_multiplier": "mult",
    "geom.ribbon_fill_stiffness_scale": "scale",
    "geom.halton_seed": "seed",
    "geom.target_stiffness_factor": "ktf",
    "geom.E_material": "E",
    "actuation.FREQUENCY": "f",
    "actuation.AMPLITUDE": "A",
    "actuation.STIFFNESS_RAMP_TIME": "ramp",
    "input2d.Tfinal": "T",
}


def sweep(fn):
    '''Registers fn under its own name.'''
    SWEEPS[fn.__name__] = fn
    return fn


def get(case, path):
    d = case
    for k in path.split("."):
        d = d[k]
    return d


def set_(case, path, value):
    *parents, key = path.split(".")
    d = case
    for k in parents:
        d = d.setdefault(k, {})
    if key not in d and parents and parents[0] in ("geom", "actuation", "input2d", "run_length"):
        # a typo would otherwise be silently ignored by the scripts that read case_params
        known = ", ".join(sorted(d))
        raise KeyError(f"{path}: not in the base case (known keys: {known}); add it to "
                       f"case_params.json first if it is new")
    d[key] = value


def _label(path):
    return ALIASES.get(path, path.split(".")[-1])


def grid(base, axes):
    '''
    Cartesian product over `axes`, a dict {path: [values]}. A key can be a tuple
    of paths that all take the same value, e.g. ("input2d.Nx", "input2d.Ny").
    Each case is named after the values it varies, in the order given, e.g.
    "Nx64_mu0.003".
    '''
    keys = list(axes)
    cases = []
    for values in itertools.product(*(axes[k] for k in keys)):
        case = copy.deepcopy(base)
        parts = []
        for key, value in zip(keys, values):
            paths = key if isinstance(key, tuple) else (key,)
            for p in paths:
                set_(case, p, value)
            parts.append(f"{_label(paths[0])}{value:g}" if isinstance(value, (int, float))
                         else f"{_label(paths[0])}{value}")
        case["case_name"] = "_".join(parts) or "base"
        cases.append(case)
    return cases


def finalize(case, sweep_name):
    '''Fills in derived values and checks the case (done for every sweep).'''
    inp = case["input2d"]
    if inp.get("Tfinal") is None:
        rl = case["run_length"]
        n_cycles = rl["transient_cycles"] + rl["measured_cycles"]
        t = n_cycles / case["actuation"]["FREQUENCY"]
        inp["Tfinal"] = round(math.ceil(t / inp["dt"] - 1e-9) * inp["dt"], 12)
    if not math.isclose(inp["Lx"] / inp["Nx"], inp["Ly"] / inp["Ny"]):
        raise ValueError(f"{case['case_name']}: dx != dy (Nx={inp['Nx']}, Ny={inp['Ny']}); "
                         f"vary them together: grid(base, {{('input2d.Nx', 'input2d.Ny'): [...]}})")
    case["sweep"] = sweep_name
    return case


def build(name, base):
    '''The finished case list of sweep `name`.'''
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
    '''The base case for 0.05 s: checks that the code runs end to end.'''
    b["input2d"]["Tfinal"] = 0.05
    cases = grid(b, {})
    cases[0]["case_name"] = "smoke"
    return cases


@sweep
def density(b):
    '''Density sweep at the base settings (the former make_sweep()).'''
    return grid(b, {"geom.ribbon_fill_density_multiplier": [0, 1, 8, 32, 64, 128, 256]})
