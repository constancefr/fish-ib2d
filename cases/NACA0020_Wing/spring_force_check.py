'''
Checks the vectorised spring force of the user-force hook
(give_Me_General_User_Defined_Force_Densities.spring_forces) against IB2d's
own give_Me_Spring_Lagrangian_Force_Densities, on the current wing.spring /
wing.vertex with random perturbations (plus one spring wrapped across the
periodic boundary), and times both.

    python Wing_Geom.py && python spring_force_check.py
'''

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'external', 'IB2d', 'pyIB2d', 'IBM_Blackbox'))
from please_Find_Lagrangian_Forces_On_Eulerian_grid import give_Me_Spring_Lagrangian_Force_Densities as ib2d_springs
from give_Me_General_User_Defined_Force_Densities import spring_forces
import wing_case


def main():
    Lx, Ly = wing_case.read_input2d('Lx'), wing_case.read_input2d('Ly')
    ds = Lx / (2 * wing_case.read_input2d('Nx'))
    x, y = np.loadtxt(os.path.join(HERE, 'wing.vertex'), skiprows=1, unpack=True)
    sp = np.loadtxt(os.path.join(HERE, 'wing.spring'), skiprows=1, ndmin=2)
    rng = np.random.default_rng(0)
    worst = 0.0
    for trial in range(3):
        xx = x + 1e-3 * rng.standard_normal(len(x))
        yy = y + 1e-3 * rng.standard_normal(len(y))
        if trial == 2:                       # move one spring end across the periodic boundary
            j = int(sp[0, 1]); xx[j] += Lx
        t0 = time.perf_counter(); a = ib2d_springs(ds, len(x), xx, yy, sp, Lx, Ly); t1 = time.perf_counter()
        b = spring_forces(ds, len(x), xx, yy, sp, Lx, Ly); t2 = time.perf_counter()
        scale = max(np.abs(a[0]).max(), np.abs(a[1]).max())
        err = max(np.abs(a[0] - b[0]).max(), np.abs(a[1] - b[1]).max()) / scale
        worst = max(worst, err)
        print('trial %d: max |hook - IB2d| / max |f| = %.2e   time IB2d %.1f ms, hook %.2f ms (%d springs)'
              % (trial, err, 1e3 * (t1 - t0), 1e3 * (t2 - t1), len(sp)))
    print('PASS' if worst < 1e-12 else 'FAIL', '(worst %.2e)' % worst)
    return worst < 1e-12


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
