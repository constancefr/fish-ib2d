'''
Called by IB2d's driver every timestep as:
    target_info = update_Target_Point_Positions(dt, current_time, target_info)

target_info columns: [Lag Pt ID, x-target, y-target, stiffness]. Only
columns 1 and 2 (x, y) get overwritten here (ID and stiffness are left
untouched).

Implements a two-corner alternating push for tail actuation: the base-bottom
corner pushes out and back (first half of each cycle) while the
base-top corner holds at baseline, then the roles swap for the second
half, repeating in a loop.

For this FinRay setup, "left" is defined as the TOP-LEFT base corner and
"right" as the BOTTOM-LEFT base corner.

Each "corner" is actually a small PATCH of several nearby Lagrangian
points (see _actuation_patch_indices in FinRay_Geom.py), not a single
point -- every point in a patch moves by the SAME rigid displacement,
approximating a real actuator pressing over a small contact area rather
than concentrating the whole push on one infinitesimal point (a single
point being pushed while its immediate neighbors held fully rigid was
showing up as an unrealistic local twist right at that point).
'''

import numpy as np
import json, os


_case_params = {}
if os.path.exists('case_params.json'):
    with open('case_params.json') as f:
        _case_params = json.load(f)
_act = _case_params.get('actuation', {})

# --- Tunable actuation parameters
# NOTE: case_params.json's "actuation" keys (and run_case.py's BASE_CASE)
# are uppercase (FREQUENCY/AMPLITUDE/...) -- these .get() calls must match
# that casing exactly, or they silently miss and fall back to the defaults
# below regardless of what's configured (this previously used lowercase
# keys, so every case_params.json override of these was a silent no-op).
FREQUENCY = _act.get('FREQUENCY', 1.0)  # Hz, flapping frequency
AMPLITUDE = _act.get('AMPLITUDE', 0.01)  # m, how far a corner is pushed perpendicular to the base
STIFFNESS_RAMP_TIME = _act.get('STIFFNESS_RAMP_TIME', 0.2)  # s, smoothly ramps target stiffness from 0 to nominal
# Target stiffness itself is set once, in FinRay_Geom.py (geom key
# 'target_stiffness_factor'); the .target file value is used as-is here.

# Active-corner sequencing over each cycle.
# The row->patch mapping is inferred from baseline y-coordinates (every
# row closer to the top corner than the bottom one is in the "left"
# patch, and vice versa) so this stays correct regardless of .target row
# order OR how many points _actuation_patch_indices put in each patch.
LEFT_ACTIVE_FIRST_HALF = True

# Cached baseline (rest) positions
_baseline = {'x': None, 'y': None}
_corner_rows = {'left': None, 'right': None}  # each an array of row indices (a patch), not a single row
_baseline_k = None


def update_Target_Point_Positions(dt, current_time, target_info):
    global _baseline_k

    # Cache the rest positions on the first call before any actuation
    if _baseline['x'] is None:
        _baseline['x'] = target_info[:, 1].copy()
        _baseline['y'] = target_info[:, 2].copy()
        _baseline_k = target_info[:, 3].copy()

        # In this geometry, both patches lie on the left edge, split by
        # height around the midline -- robust to patch size since every
        # point in the top patch is, by construction (a small radius
        # around the true corner), much closer to the top corner than to
        # the bottom one, and vice versa.
        midline_y = (_baseline['y'].max() + _baseline['y'].min()) / 2.0
        _corner_rows['left'] = np.where(_baseline['y'] > midline_y)[0]
        _corner_rows['right'] = np.where(_baseline['y'] <= midline_y)[0]

    left_rows = _corner_rows['left']
    right_rows = _corner_rows['right']

    T = 1.0 / FREQUENCY
    phase = (2 * np.pi * current_time / T) % (2 * np.pi)  # 0 to 2*pi each cycle

    # Smooth (zero velocity at both ends of each active half-cycle) bump
    # using sin^2, rather than plain sin to avoid velocity discontinuity
    disp_first = AMPLITUDE * np.sin(phase) ** 2 if phase < np.pi else 0.0
    disp_second = AMPLITUDE * np.sin(phase - np.pi) ** 2 if phase >= np.pi else 0.0

    if LEFT_ACTIVE_FIRST_HALF:
        disp_left = disp_first
        disp_right = disp_second
    else:
        disp_left = disp_second
        disp_right = disp_first

    # Perpendicular to base displacement = x-direction. Every point in a
    # patch gets the SAME rigid displacement (moves together as one small
    # block), not an independent one -- that's what makes it behave like a
    # patch instead of several unrelated point pushes.
    target_info[left_rows, 1] = _baseline['x'][left_rows] + disp_left
    target_info[left_rows, 2] = _baseline['y'][left_rows]

    target_info[right_rows, 1] = _baseline['x'][right_rows] + disp_right
    target_info[right_rows, 2] = _baseline['y'][right_rows]

    # Smoothly turn on target stiffness to avoid an impulsive force spike
    # in the first few steps of explicit coupling.
    if STIFFNESS_RAMP_TIME > 0.0:
        ramp = min(max(current_time / STIFFNESS_RAMP_TIME, 0.0), 1.0)
    else:
        ramp = 1.0
    target_info[:, 3] = _baseline_k * ramp

    # if current_time % 0.02 < dt:
    #     print(f"t={current_time:.3f}s: left patch x={target_info[left_rows, 1]}, right patch x={target_info[right_rows, 1]}")
    #     print(f"t={current_time:.3f}s: left patch y={target_info[left_rows, 2]}, right patch y={target_info[right_rows, 2]}")

    return target_info