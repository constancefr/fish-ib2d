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
'''

import numpy as np
import json, os


_case_params = {}
if os.path.exists('case_params.json'):
    with open('case_params.json') as f:
        _case_params = json.load(f)
_act = _case_params.get('actuation', {})

# --- Tunable actuation parameters
FREQUENCY = _act.get('frequency', 1.0)  # Hz, flapping frequency
AMPLITUDE = _act.get('amplitude', 0.01)  # m, how far a corner is pushed perpendicular to the base
STIFFNESS_RAMP_TIME = _act.get('stiffness_ramp_time', 0.2)  # s, smoothly ramps target stiffness from 0 to nominal
TARGET_STIFFNESS_SCALE = _act.get('target_stiffness_scale', 0.2)  # 1.0 uses file value; lower values reduce destabilizing forcing

# Active-corner sequencing over each cycle.
# The row mapping for LEFT/RIGHT is inferred from baseline y-coordinates
# so this stays correct even if .target row order changes.
LEFT_ACTIVE_FIRST_HALF = True

# Cached baseline (rest) positions
_baseline = {'x': None, 'y': None}
_corner_rows = {'left': None, 'right': None}
_baseline_k = None


def update_Target_Point_Positions(dt, current_time, target_info):
    global _baseline_k

    # Cache the rest positions on the first call before any actuation
    if _baseline['x'] is None:
        _baseline['x'] = target_info[:, 1].copy()
        _baseline['y'] = target_info[:, 2].copy()
        _baseline_k = target_info[:, 3].copy()

        # In this geometry, both target points lie on the left edge
        _corner_rows['left'] = int(np.argmax(_baseline['y']))
        _corner_rows['right'] = int(np.argmin(_baseline['y']))

    left_row = _corner_rows['left']
    right_row = _corner_rows['right']

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

    # Perpendicular to base displacement = x-direction
    target_info[left_row, 1] = _baseline['x'][left_row] + disp_left
    target_info[left_row, 2] = _baseline['y'][left_row]

    target_info[right_row, 1] = _baseline['x'][right_row] + disp_right
    target_info[right_row, 2] = _baseline['y'][right_row]

    # Smoothly turn on target stiffness to avoid an impulsive force spike
    # in the first few steps of explicit coupling.
    if STIFFNESS_RAMP_TIME > 0.0:
        ramp = min(max(current_time / STIFFNESS_RAMP_TIME, 0.0), 1.0)
    else:
        ramp = 1.0
    target_info[:, 3] = _baseline_k * TARGET_STIFFNESS_SCALE * ramp

    # if current_time % 0.02 < dt:
    #     print(f"t={current_time:.3f}s: left corner x={target_info[left_row, 1]:.4f}, right corner x={target_info[right_row, 1]:.4f}")
    #     print(f"t={current_time:.3f}s: left corner y={target_info[left_row, 2]:.4f}, right corner y={target_info[right_row, 2]:.4f}")

    return target_info