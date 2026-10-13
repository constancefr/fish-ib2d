"""
IB2d unit conventions shared by every case.

IB2d is 2D: every force, stiffness and mass is per unit depth (span = 1 m),
in SI units. The driver also treats Lagrangian forces as densities and
spreads `f * ds` onto the grid with a *global* ds = Lx / (2 Nx)
(please_Find_Lagrangian_Forces_On_Eulerian_grid.py), whatever the actual
spacing of the Lagrangian points. A physical nodal force F [N per unit depth]
therefore needs f = F / ds, and every spring/target constant written to an
input file is K_phys / ds.
"""


def ib_ds(Lx, Nx):
    '''IB2d's global force-spreading length ds = Lx / (2 Nx) [m]. Not the
    actual point spacing, which may be finer (e.g. FinRay's interior fill).'''
    return Lx / (2.0 * Nx)


def to_ib(K_phys, ds):
    '''Physical stiffness [N/m per unit depth] -> IB2d spring/target constant.'''
    return K_phys / ds
