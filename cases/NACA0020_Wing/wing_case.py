'''
Shared helpers for the NACA0020 wing example: reads case_params.json and
input2d (the single source of truth for the grid / fluid values), and
holds the angle-of-attack schedule used by the free-stream forcing.
'''

import json
import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))


def load_case_params():
    path = os.path.join(_HERE, 'case_params.json')
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def read_input2d(key):
    '''Value of `key = <number>` in input2d (comments after '%' ignored).'''
    with open(os.path.join(_HERE, 'input2d')) as f:
        for line in f:
            m = re.match(r'\s*' + re.escape(key) + r'\s*=\s*([-+0-9.eE]+)', line.split('%')[0] + ' ')
            if m:
                return float(m.group(1))
    raise KeyError(key)


def smoothstep(s):
    '''0 for s<=0, 1 for s>=1, s^2(3-2s) between (same kernel as WaterLily_FSI value_schedule :cubic).'''
    s = min(max(s, 0.0), 1.0)
    return s * s * (3.0 - 2.0 * s)


def aoa_rad(t, flow):
    '''Free-stream angle at time t: aoa_deg[0] -> aoa_deg[1] at aoa_switch_time over aoa_ramp.'''
    a0, a1 = flow['aoa_deg']
    frac = smoothstep((t - flow['aoa_switch_time']) / flow['aoa_ramp'])
    return (a0 + (a1 - a0) * frac) * 3.141592653589793 / 180.0
