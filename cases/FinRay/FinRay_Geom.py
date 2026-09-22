"""
Creates the Lagrangian structure geometry for a finray tail and writes the
associated input files (.vertex, .geo_connect, .spring, .beam, .target).

Structure model: each originally-1-point-thick edge (the outline + the 4
internal rays) is built as a two-rail "ribbon" of the given wall_thickness,
with the rails tied together by rung + diagonal (Warren-truss) springs.
Bending resistance of the tail emerges from this truss geometry itself (one
rail stretches while the other compresses) rather than from an explicit
EI-derived torsional spring on a single centerline, which is a much closer
approximation of how a solid silicone cross-section actually deforms. The
legacy single-point-thick implementation is kept below (commented out) for
reference.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as _MplPath
import json, os, sys

_EXAMPLES_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _EXAMPLES_DIR not in sys.path:
    sys.path.insert(0, _EXAMPLES_DIR)
import _mesh_utils as mu

_case_params = {}
if os.path.exists('case_params.json'):
    with open('case_params.json') as f:
        _case_params = json.load(f)
_geom = _case_params.get('geom', {})


def FinRay_Geom():
    # Eulerian grid params (MAKE SURE MATCHES THOSE IN input2d !!!)
    Nx = 32           # MUST BE EVEN!!!
    Lx = 1.0          # length of grid in x-direction
    dx = Lx / Nx      # spatial resolution

    # Immersed structure geometric / dynamic params
    ds = 0.5 * dx               # Lagrangian point spacing (2x resolution of Eulerian grid!!!)
    struct_name = 'finray'    # name for .vertex, .spring, etc files (must match what's in 'input2d')

    #-------------------------------------------------------------------
    # Material properties - trying to mimick real system from dimensions and Dragon Skin 30
    # silicone properties
    #-------------------------------------------------------------------
    E_material = _geom.get('E_material', 0.74e6)
    wall_thickness = _geom.get('wall_thickness', 0.00396)   # m (3.96mm in real build)
    extrude_depth = _geom.get('extrude_depth', 0.05)        # m

    # Ribbon interior-fill params (Halton-sampled points inside each
    # ribbon's actual 2D area, meshed via Delaunay triangulation). See
    # fill_polygon_with_halton()/triangulate_with_holes() in _mesh_utils.py.
    #   ribbon_fill_spacing: target point spacing; None -> falls back to ds.
    #   ribbon_fill_density_multiplier: scales the density-derived point
    #     count (0 = no fill, exactly today's sparse-truss geometry; 1 =
    #     "same density as ds"; push higher, e.g. 8-32, to densify the mesh).
    #   ribbon_fill_stiffness_scale: scales the interior 'fill' spring
    #     stiffness formula (see build_Tail_Ribbon_Connections) -- flagged
    #     as an approximation to calibrate empirically, see NOTES.md.
    fill_spacing = _geom.get('ribbon_fill_spacing', None)
    fill_density_multiplier = _geom.get('ribbon_fill_density_multiplier', 1.0)
    fill_stiffness_scale = _geom.get('ribbon_fill_stiffness_scale', 1.0)

    # Construct geometry: each edge/ray is now a 2-rail ribbon of width
    # wall_thickness (instead of a single 1-point-thick line). Also returns
    # which rail indices belong to the outline vs. each ray, needed to wire
    # up connections correctly below, and ribbonPolys (each ribbon's
    # boundary/interior points) needed for the interior mesh fill.
    xLag, yLag, outlineRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys = \
        _Build_Tail_Geometry_Ribbon(ds, wall_thickness,
                                     fill_spacing=fill_spacing,
                                     fill_density_multiplier=fill_density_multiplier)

    # Build connections: rail springs ('outline'/'ray'), cross-thickness
    # 'rung' + 'diag' springs that give each ribbon its shear stiffness,
    # 'attach' springs tying each ray's rail-ends to the outline, and 'fill'
    # springs from Delaunay-triangulating each ribbon's interior Halton
    # points against its own boundary.
    connections, fill_edge_k = build_Tail_Ribbon_Connections(
        xLag, yLag, outlineRailIdx, rayRailIdxList, ribbonPolys,
        E_material=E_material, extrude_depth=extrude_depth,
        fill_stiffness_scale=fill_stiffness_scale)

    # Plot geometry to test
    # NOTE: colors must be passed via color= (not as the positional fmt
    # string) since 'fill' needs an RGBA tuple -- plt.plot(x, y, a_tuple)
    # does NOT set the color, it silently treats the tuple as a *second*
    # y-dataset (with auto x=[0,1,2,...]), which draws a bogus line across
    # the whole axes and blows out the plot's autoscale.
    plot_colors = {'outline': 'b', 'ray': 'b', 'rung': 'c', 'diag': 'y',
                    'fill': (0.6, 0.6, 0.6, 0.6)}
    for (i, j, kind) in connections:
        lw = 0.4 if kind == 'fill' else 0.75
        plt.plot([xLag[i], xLag[j]], [yLag[i], yLag[j]],
                  color=plot_colors.get(kind, 'k'), linewidth=lw)
    for rail in outlineRailIdx:
        plt.plot(xLag[rail], yLag[rail], 'r*')
    for railPair in rayRailIdxList:
        for rail in railPair:
            plt.plot(xLag[rail], yLag[rail], 'g*')
    for ribbon in ribbonPolys:
        if ribbon['interior_idx']:
            plt.plot(xLag[ribbon['interior_idx']], yLag[ribbon['interior_idx']], 'k.', markersize=2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('square')
    plt.savefig('geometry_preview.png', dpi=150)
    plt.show()

    # Print .vertex file
    print_Lagrangian_Vertices(xLag, yLag, struct_name)

    # Print .geo_connect file
    print_Geometry_Connections(connections, struct_name)

    # Print .spring file
    #   k_Spring/k_Rung: nominal (ds-based) reference stiffnesses -- still
    #     used below for k_Beam/k_Target, which want a single representative
    #     magnitude, not a per-spring value.
    #   Each individual 'outline'/'ray'/'rung'/'diag' spring's actual
    #     stiffness is now computed per-spring from its OWN rest length
    #     (k = E*A / L_own, the standard truss-element formula -- see the
    #     length-aware stiffness note in NOTES.md), not from a single
    #     constant shared by every spring of that kind. Recovering the E*A
    #     (or G*A) product from k_Spring/k_Rung avoids redefining
    #     spring_Stiffness_From_Material/shear_Stiffness_From_Material:
    #     since k_Spring = E*A/ds, E*A = k_Spring*ds.
    #   Ray-to-outline junctions no longer get a separate spring constant --
    #     each ray is now welded directly into the outline's own rail point
    #     sequence (see weld_to_inner_rail in _Build_Tail_Geometry_Ribbon),
    #     so that joint is just an ordinary 'outline'/'ray' rail spring.
    k_Spring = spring_Stiffness_From_Material(E_material, wall_thickness / 2.0, extrude_depth, ds)
    k_Rung = shear_Stiffness_From_Material(E_material, wall_thickness, extrude_depth, ds)
    EA_rail = k_Spring * ds   # = E * (wall_thickness/2 * extrude_depth)
    GA_web = k_Rung * ds      # = G * (wall_thickness * extrude_depth)
    # 'fill' (interior mesh) springs are resolved through fill_edge_k
    # instead; EA_rail here is just a defensive fallback.
    EA_by_kind = {'outline': EA_rail, 'ray': EA_rail, 'rung': GA_web, 'diag': GA_web,
                  'fill': EA_rail}
    print_Lagrangian_Springs(xLag, yLag, connections, EA_by_kind, struct_name, fill_edge_k=fill_edge_k)

    # Print .beam file
    #   NOTE: the tail's bending resistance comes from the ribbon truss
    #   above (rail springs + rungs/diagonals), NOT from these beams. These
    #   are a light per-rail regularization only, to stop an individual
    #   rail's point-chain from zig-zagging under local compression (a
    #   standard IBM issue for curves discretized purely with springs) -
    #   hence why k_Beam is just a small fraction of k_Spring rather than a
    #   value derived from a real EI.
    beam_reg_fraction = 0.05
    k_Beam = beam_reg_fraction * k_Spring
    beams = build_Tail_Ribbon_Beams(outlineRailIdx, rayRailIdxList)
    print_Lagrangian_Beams(xLag, yLag, beams, k_Beam, struct_name)

    # Print .target file
    #   k_Target is set to very high to hold the points nearly rigid.
    #   Actuation pushes each base corner as a small rigid PATCH of nearby
    #   points (not just the single corner vertex) -- mimicking a real
    #   actuator, which presses over a small contact area rather than an
    #   infinitesimal point. A single-point push otherwise concentrates all
    #   the forcing on one Lagrangian point, which showed up as an
    #   unrealistic local twist right at that point rather than a
    #   distributed push (see NOTES.md).
    k_Target = 1000.0 * k_Spring
    actuation_patch_radius = _geom.get('actuation_patch_radius', None)
    if actuation_patch_radius is None:
        # Tied to wall_thickness, NOT ds: the two rails at a corner are
        # only ~wall_thickness apart, so this reliably captures "both
        # rails of this corner" (the actual physical extent of the
        # material there) regardless of mesh resolution. A ds-based
        # default was tried first and rejected -- for this geometry, ds
        # happens to be close to the distance to the *adjacent* slant
        # edge's nearest point, so a ds-sized radius pulled in a point
        # from a different edge (not more of the base itself), which
        # measurably destabilized the sim (NaN) at high
        # ribbon_fill_density_multiplier -- see NOTES.md. Going much
        # bigger than wall_thickness risks the same failure mode again.
        actuation_patch_radius = 1.0 * wall_thickness
    target_indices = _actuation_patch_indices(xLag, yLag, baseCornerIdx, actuation_patch_radius,
                                               ribbonPolys=ribbonPolys)
    print_Lagrangian_Target_Pts(target_indices, k_Target, struct_name)



def print_Lagrangian_Vertices(xLag, yLag, struct_name):
    '''
    Prints the Lagrangian points to .vertex file.

    xLag: x-coordinates of Lagrangian points
    yLag: y-coordinates of Lagrangian points
    '''
    N = len(xLag)  # total num of Lagrangian points

    with open(struct_name + '.vertex', 'w') as f:
        f.write('%d\n' % N)
        # Loops over all points
        for s in range(N):
            f.write('%1.16e %1.16e\n' % (xLag[s], yLag[s]))

def print_Geometry_Connections(connections, struct_name):
    '''
    Prints the connections to .geo_connect file.

    connections: list of (i, j, kind) generated by build_Tail_Ribbon_Connections()
    '''
    with open(struct_name + '.geo_connect', 'w') as f:
        for (i, j, kind) in connections:
            s1, s2 = i, j
            f.write('%d %d\n' % (s1, s2))
            f.write('%d %d\n' % (s2, s1))

# --- Legacy (single-point-thick) implementation, kept for reference -------
#
# def print_Lagrangian_Springs(xLag, yLag, connections, k_Spring, k_Attach, struct_name, deg_NL=1.0):
#     '''
#     Prints springs to .spring file.
#
#     connections: list of (i, j, kind) generated by build_Tail_Connections()
#     k_Spring: spring constant for the material
#     k_Attach: spring constant for the attachments
#     struct_name: name of the structure
#     deg_NL: degree of nonlinearity
#     '''
#     N = len(connections)
#
#     with open(struct_name + '.spring', 'w') as f:
#         f.write('%d\n' % N)   # Print # of springs
#
#         for (i, j, kind) in connections:
#             ds_Rest = float(np.hypot(xLag[j] - xLag[i], yLag[j] - yLag[i])) # set initial distance as resting length
#             k = k_Attach if kind == 'attach' else k_Spring # fused junctions get stiffer springs
#             f.write('%d %d %1.16e %1.16e %1.16e\n' % (i, j, k, ds_Rest, deg_NL))
#
# ---------------------------------------------------------------------------

def print_Lagrangian_Springs(xLag, yLag, connections, EA_by_kind, struct_name,
                              fill_edge_k=None, deg_NL=1.0):
    '''
    Prints springs to .spring file.

    connections: list of (i, j, kind) generated by build_Tail_Ribbon_Connections()
    EA_by_kind: dict mapping each connection 'kind' (see build_Tail_Ribbon_Connections)
        to that kind's modulus x area product (E*A for 'outline'/'ray' rail
        springs, G*A for 'rung'/'diag' web springs -- NOT yet divided by a
        length). Each spring's actual stiffness is computed here as
        k = EA_by_kind[kind] / ds_Rest, i.e. from THAT SPRING'S OWN rest
        length, not a single nominal-ds-based constant shared by every
        spring of that kind -- see the "length-aware stiffness" note in
        NOTES.md for why this matters (in short: k=EA/L is the standard
        truss-element formula, and treating L as a fixed nominal ds instead
        of each spring's actual length makes any short spring -- e.g. one
        created by weld_to_inner_rail's junction splicing -- silently too
        soft for its length).
    struct_name: name of the structure
    fill_edge_k: dict mapping (min(i,j), max(i,j)) -> stiffness for 'fill'
        (interior mesh) connections, whose stiffness is computed directly
        per-edge in build_Tail_Ribbon_Connections (a different, tributary-
        area-based formula, not the EA/L one here). Falls back to
        EA_by_kind['fill'] / ds_Rest if None.
    deg_NL: degree of nonlinearity
    '''
    N = len(connections)

    with open(struct_name + '.spring', 'w') as f:
        f.write('%d\n' % N)   # Print # of springs

        for (i, j, kind) in connections:
            ds_Rest = float(np.hypot(xLag[j] - xLag[i], yLag[j] - yLag[i])) # set initial distance as resting length
            if kind == 'fill' and fill_edge_k is not None:
                k = fill_edge_k[(min(i, j), max(i, j))]
            else:
                k = EA_by_kind[kind] / ds_Rest
            f.write('%d %d %1.16e %1.16e %1.16e\n' % (i, j, k, ds_Rest, deg_NL))

def _actuation_patch_indices(xLag, yLag, baseCornerIdx, radius, ribbonPolys=None):
    '''
    Expands baseCornerIdx's 2 corners (each already 2 points -- both rails)
    into a "patch" of every Lagrangian point within `radius` of that
    corner's true position -- so update_Target_Point_Positions.py can move
    a whole small rigid cluster of nearby points together instead of just
    the single corner vertex, approximating a real actuator's finite
    contact area rather than a point force.

    baseCornerIdx: (railA_bottom, railB_bottom, railA_top, railB_top), as
        returned by _Build_Tail_Geometry_Ribbon.
    radius: patch radius, same units as xLag/yLag (meters). Should stay
        well under half the base width, or the two patches could merge.
    ribbonPolys: if given, restricts the patch to BOUNDARY (rail) points
        only, excluding every ribbon's Halton-fill interior points -- a
        real actuator grips the material's outer surface, not its bulk
        interior, and forcing many densely-packed fill points into exact
        rigid lockstep (their count scales with ribbon_fill_density_
        multiplier, so this can mean dozens of them at high density) via
        the very stiff target springs was found to destabilize the
        simulation (NaN) at high multiplier/stiffness-scale combinations
        that were otherwise stable -- see NOTES.md. None restricts to
        nothing (every point is eligible), matching the pre-fix search.

    Returns a sorted list of point indices (the union of both patches;
    always includes the original 4 corner points regardless of radius).
    '''
    xLag, yLag = np.asarray(xLag), np.asarray(yLag)
    bottom_ref = np.array([xLag[baseCornerIdx[0]] + xLag[baseCornerIdx[1]],
                            yLag[baseCornerIdx[0]] + yLag[baseCornerIdx[1]]]) / 2.0
    top_ref = np.array([xLag[baseCornerIdx[2]] + xLag[baseCornerIdx[3]],
                         yLag[baseCornerIdx[2]] + yLag[baseCornerIdx[3]]]) / 2.0

    pts = np.column_stack((xLag, yLag))
    near_bottom = np.linalg.norm(pts - bottom_ref, axis=1) <= radius
    near_top = np.linalg.norm(pts - top_ref, axis=1) <= radius
    eligible = near_bottom | near_top
    if ribbonPolys is not None:
        interior_idx = set()
        for ribbon in ribbonPolys:
            interior_idx.update(ribbon['interior_idx'])
        for idx in interior_idx:
            eligible[idx] = False
    patch = set(np.where(eligible)[0].tolist())
    patch.update(baseCornerIdx)  # guarantee the true corners are always included
    return sorted(patch)


def print_Lagrangian_Target_Pts(target_indices, k_Target, struct_name):
    '''
    Prints the target points to .target file.

    target_indices: list of 0-indexed points to actuate
    k_Target: spring constant for the target points
    '''
    N = len(target_indices)

    with open(struct_name + '.target', 'w') as f:
        f.write('%d\n' % N)
        for idx in target_indices:
            f.write('%d %1.16e\n' % (idx, k_Target))

def print_Lagrangian_Beams(xLag, yLag, beams, k_Beam, struct_name):
    '''
    Prints the beams to .beam file.

    beams: list of (L, M, R) 0-indexed point-index triples generated by build_Tail_Ribbon_Beams()
    k_Beam: beam stiffness
    '''
    N = len(beams)

    with open(struct_name + '.beam', 'w') as f:
        f.write('%d\n' % N)

        for (L, M, R) in beams:
            Xp, Yp = xLag[L], yLag[L]
            Xq, Yq = xLag[M], yLag[M]
            Xr, Yr = xLag[R], yLag[R]
            C = (Xr - Xq) * (Yq - Yp) - (Yr - Yq) * (Xq - Xp)  # curvature (signed area of triangle formed by the three points L, M, R)
            f.write('%d %d %d %1.16e %1.16e\n' % (L, M, R, k_Beam, C))



def spring_Stiffness_From_Material(E, thickness, depth, ds):
    """
    Converts a material's Young's modulus into IB2d's discrete spring
    stiffness, treating each spring as an axial rod segment of rest
    length ds and cross-sectional area (thickness x depth):

        k = E * A / ds      (standard truss-element axial stiffness)

    E:  Young's modulus, Pa. NOTE: silicone is hyperelastic (nonlinear
        stress-strain), so a single E is a linear approximation.
    thickness: effective in-plane thickness carried by this spring, m.
        For the ribbon's rail ('outline'/'ray') springs, pass
        wall_thickness/2 - each rail is one "flange" of a 2-flange truss
        idealization of the wall's full cross-section, so it only carries
        half the material.
    depth: out-of-plane extrusion depth, m.
    ds: current Lagrangian point spacing (hence why k_Spring is recomputed
        if Nx/Lx is changed, rather than hardcoded).
    """
    A = thickness * depth
    return E * A / ds


def shear_Stiffness_From_Material(E, thickness, depth, ds, nu=0.499):
    """
    Converts a material's Young's modulus into the discrete spring
    stiffness for a ribbon's cross-thickness ('rung'/'diag') springs,
    which stand in for the "web" shear resistance of the ribbon's
    idealized 2-flange beam cross-section:

        G = E / (2*(1+nu))     (shear modulus; nu ~ 0.499 for near-
                                 incompressible silicone rubber)
        k = G * A / ds,   A = thickness * depth

    thickness: full ribbon wall thickness (unlike the flange/rail springs
        in spring_Stiffness_From_Material, the web spans the whole gap
        between the two rails, not just half of it).
    depth: out-of-plane extrusion depth, m.
    ds: current Lagrangian point spacing.
    """
    G = E / (2.0 * (1.0 + nu))
    A = thickness * depth
    return G * A / ds


# beam_Stiffness_From_Material() is unused by the current ribbon model
# (bending comes from the truss geometry itself - see build_Tail_Ribbon_Connections),
# but is kept as a reference for the legacy single-centerline EI-beam formula.
def beam_Stiffness_From_Material(E, thickness, depth, ds):
    """
    Converts a material's Young's modulus into IB2d's discrete beam
    (torsional spring) stiffness, using the standard scaling for a
    discretized beam approximating a continuum flexural rigidity EI:

        I = thickness^3 * depth / 12   (rectangular cross-section,
                                         thickness in the bending direction)
        k_Beam = E * I / ds

    NOTE: IB2d's torsional spring energy is defined via a signed-area
    form of "curvature" rather than a normalized curvature. Might have to
    verify that the scaling of k_Beam with E and I is correct for IB2d's
    discrete implementation.
    """
    I = (thickness ** 3) * depth / 12.0
    return E * I / ds


# --- Legacy (single-point-thick) implementation, kept for reference -------
#
# def build_Tail_Beams(outlineIdx, rayIdxList):
#     '''
#     Builds the list of (left, middle, right) 0-indexed point-index
#     triples that get bending resistance.
#     Same piece-aware logic as build_Tail_Connections: the outline is
#     a closed loop (wraps around), each ray is an open chain (does not),
#     and rays get no bending connection to each other.
#
#     NOTE: this does not add a beam spanning the
#     ray-to-outline junction itself, so bending moment isn't transmitted
#     across that connection (only the axial attachment spring force is
#     via build_Tail_Connections). The junction currently behaves like a hinge
#     instead of a fused connection.
#     '''
#
#     beams = []
#
#     # Outline: closed loop, one beam centered at every point
#     No = len(outlineIdx)
#     for i in range(No):
#         L = outlineIdx[(i - 1) % No]
#         M = outlineIdx[i]
#         R = outlineIdx[(i + 1) % No]
#         beams.append((L, M, R))
#
#     # Each ray: open chain, only interior points get a beam
#     for ray in rayIdxList:
#         for i in range(1, len(ray) - 1):
#             beams.append((ray[i - 1], ray[i], ray[i + 1]))
#
#     return beams
#
# ---------------------------------------------------------------------------

def build_Tail_Ribbon_Beams(outlineRailIdx, rayRailIdxList):
    '''
    Builds light per-rail regularization beams for the ribboned structure.

    NOTE: these are NOT the source of the tail's bending stiffness (that
    comes from the ribbon truss itself - see build_Tail_Ribbon_Connections);
    they only prevent an individual rail's point-chain from zig-zagging
    under local compression, so FinRay_Geom() gives them a small stiffness
    (beam_reg_fraction * k_Spring) rather than one derived from a real EI.

    Same piece-aware logic as the legacy build_Tail_Beams: the outline's
    two rails are closed loops (wrap around), each ray's two rails are open
    chains (do not), and rays get no bending connection to each other.
    '''
    beams = []

    def add_rail_beams(railIdx, closed):
        n = len(railIdx)
        if closed:
            for i in range(n):
                beams.append((railIdx[(i - 1) % n], railIdx[i], railIdx[(i + 1) % n]))
        else:
            for i in range(1, n - 1):
                beams.append((railIdx[i - 1], railIdx[i], railIdx[i + 1]))

    railA, railB = outlineRailIdx
    add_rail_beams(railA, True)
    add_rail_beams(railB, True)

    for (rA, rB) in rayRailIdxList:
        add_rail_beams(rA, False)
        add_rail_beams(rB, False)

    return beams


def give_Me_Immersed_Boundary_Geometry(ds, Nx):
    '''
    Returns the Lagrangian coordinates of the immersed boundary geometry.

    ds: Lagrangian point spacing
    Nx: Eulerian grid resolution
    '''
    wall_thickness = _geom.get('wall_thickness', 0.00396)
    xLag, yLag, _, _, _, _ = _Build_Tail_Geometry_Ribbon(
        ds, wall_thickness,
        fill_spacing=_geom.get('ribbon_fill_spacing', None),
        fill_density_multiplier=_geom.get('ribbon_fill_density_multiplier', 1.0))
    return xLag, yLag


# --- Legacy (single-point-thick) implementation, kept for reference -------
#
# def build_Tail_Connections(xLag, yLag, outlineIdx, rayIdxList):
#     '''
#     Builds the list of (i, j, kind) triples.
#     kind is one of:
#       'outline': the closed outline loop (last point wraps to first).
#       'ray':     each ray's open internal chain.
#       'attach':  each ray's two endpoints connected to their nearest
#                  outline point.
#
#     kind is returned explicitly so the caller can give 'attach' springs
#     a different (stiffer) constant.
#     '''
#
#     connections = []
#
#     # Outline: closed loop
#     No = len(outlineIdx)
#     for k in range(No):
#         i = outlineIdx[k]
#         j = outlineIdx[(k + 1) % No]     # wraps last -> first, closing the loop
#         connections.append((i, j, 'outline'))
#
#     # Each ray: open chain (no wraparound and no cross-ray links)
#     for ray in rayIdxList:
#         for k in range(len(ray) - 1):
#             connections.append((ray[k], ray[k + 1], 'ray'))
#
#     # Attach each ray's endpoints to their nearest outline point
#     outline_xy = np.column_stack((xLag[outlineIdx], yLag[outlineIdx]))
#     for ray in rayIdxList:
#         for endpoint in (ray[0], ray[-1]):
#             p = np.array([xLag[endpoint], yLag[endpoint]])
#             dists = np.linalg.norm(outline_xy - p, axis=1)
#             nearest_idx = outlineIdx[int(np.argmin(dists))]
#             connections.append((endpoint, nearest_idx, 'attach'))
#
#     return connections
#
# ---------------------------------------------------------------------------

def build_Tail_Ribbon_Connections(xLag, yLag, outlineRailIdx, rayRailIdxList,
                                   ribbonPolys=None, E_material=None,
                                   extrude_depth=None, fill_stiffness_scale=1.0):
    '''
    Builds the list of (i, j, kind) triples for the ribboned structure.
    kind is one of:
      'outline': along each of the outline's two rails (closed loop).
      'ray':     along each ray's two rails (open chain).
      'rung':    cross-thickness spring tying corresponding points on a
                 ribbon's two rails together (its "web").
      'diag':    alternating cross-thickness diagonal (Warren truss) that
                 gives the ribbon shear stiffness - without these, a ribbon
                 built from 'rung' springs alone is a mechanism that can
                 freely shear into a parallelogram.
      'fill':    edges from Delaunay-triangulating a ribbon's interior
                 Halton-sampled points against its own boundary (see
                 ribbonPolys / _mesh_utils.triangulate_with_holes below).
                 Only added if ribbonPolys is given.

    kind is returned explicitly so the caller can give each a different
    stiffness (e.g. 'attach' stiffer/softer than the material itself).

    ribbonPolys: optional list of dicts (see _Build_Tail_Geometry_Ribbon),
        one per ribbon, each with 'boundary_idx'/'exterior'/'holes'/
        'interior_idx'. When given (together with E_material and
        extrude_depth), each ribbon's interior points are triangulated and
        added as 'fill' connections.

    Returns (connections, fill_edge_k): fill_edge_k maps each 'fill' edge's
    (min(i,j), max(i,j)) pair to its stiffness (see the tributary-area
    formula below); it's empty if ribbonPolys is None.
    '''

    connections = []

    def add_rail_chain(railIdx, closed, kind):
        n = len(railIdx)
        stop = n if closed else n - 1
        for k in range(stop):
            connections.append((railIdx[k], railIdx[(k + 1) % n], kind))

    def add_rungs_and_diagonals(railA, railB, closed):
        n = len(railA)
        for k in range(n):
            connections.append((railA[k], railB[k], 'rung'))
        stop = n if closed else n - 1
        for k in range(stop):
            kNext = (k + 1) % n
            if k % 2 == 0:            # alternate diagonal direction each cell
                connections.append((railA[k], railB[kNext], 'diag'))
            else:
                connections.append((railB[k], railA[kNext], 'diag'))

    # Outline ribbon: closed loop
    railA, railB = outlineRailIdx
    add_rail_chain(railA, True, 'outline')
    add_rail_chain(railB, True, 'outline')
    add_rungs_and_diagonals(railA, railB, True)

    # Each ray's ribbon: open chain (no wraparound and no cross-ray links)
    for (rA, rB) in rayRailIdxList:
        add_rail_chain(rA, False, 'ray')
        add_rail_chain(rB, False, 'ray')
        add_rungs_and_diagonals(rA, rB, False)

    # NOTE: rays no longer get a separate 'attach' spring to the outline --
    # each ray rail's endpoints are now welded directly into the outline's
    # own inner-rail point sequence at geometry-construction time (see
    # weld_to_inner_rail in _Build_Tail_Geometry_Ribbon), so they're already
    # covered by the 'outline'/'ray'/'rung'/'diag' connections built above:
    # a welded endpoint is simultaneously a member of the outline's rail
    # chain (gets 'outline' springs to its outline neighbors) and of its
    # ray's rail chain (gets 'ray'/'rung'/'diag' springs there too).

    # Two different rail endpoints can legitimately weld to the very same
    # existing outline vertex (or two adjacent rung/diag pairs can
    # coincidentally match an already-built 'outline' pair) -- drop any
    # resulting self-loop (i == j, which would otherwise become a
    # zero-length, divide-by-zero spring) and any exact duplicate (i, j)
    # pair, keeping the first (lowest-priority-kind) occurrence.
    seen_pairs = set()
    deduped = []
    for (i, j, kind) in connections:
        if i == j:
            continue
        pair = (min(i, j), max(i, j))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        deduped.append((i, j, kind))
    connections = deduped

    # Interior mesh fill: Delaunay-triangulate each ribbon's boundary +
    # interior Halton points, and add any triangulation edge not already
    # covered by the rail/rung/diag/attach logic above as a 'fill' spring.
    # Stiffness is derived per-edge (not per-kind, unlike everything else
    # here) via a tributary-area lattice-spring approximation -- see
    # NOTES.md for the derivation and its caveats.
    fill_edge_k = {}
    if ribbonPolys:
        existing_pairs = {(min(i, j), max(i, j)) for (i, j, _k) in connections}
        for ribbon in ribbonPolys:
            b_idx, i_idx = ribbon['boundary_idx'], ribbon['interior_idx']
            if not i_idx:
                # No interior fill points for this ribbon (e.g.
                # density_multiplier=0, or a ribbon too thin for even one
                # point) -- skip triangulation entirely rather than
                # triangulating the boundary points alone, which would add
                # spurious boundary-to-boundary 'fill' diagonals the sparse
                # truss never had.
                continue
            boundary_pts = np.column_stack((xLag[b_idx], yLag[b_idx]))
            interior_pts = np.column_stack((xLag[i_idx], yLag[i_idx]))
            all_pts, edges_local, tris = mu.triangulate_with_holes(
                ribbon['exterior'], ribbon['holes'], boundary_pts, interior_pts)

            # Per-triangle area split evenly (1/3) across its 3 edges, so
            # each edge's stiffness reflects the material tributary to it.
            edge_area = {}
            for simplex in tris:
                p0, p1, p2 = all_pts[simplex]
                area = 0.5 * abs((p1[0] - p0[0]) * (p2[1] - p0[1])
                                  - (p2[0] - p0[0]) * (p1[1] - p0[1]))
                i, j, k = simplex
                for a, b in ((i, j), (j, k), (k, i)):
                    key = (min(a, b), max(a, b))
                    edge_area[key] = edge_area.get(key, 0.0) + area / 3.0

            local_to_global = list(b_idx) + list(i_idx)
            for (a, b) in edges_local:
                gi, gj = local_to_global[a], local_to_global[b]
                pair = (min(gi, gj), max(gi, gj))
                if pair in existing_pairs:
                    continue  # already covered by outline/ray/rung/diag/attach
                existing_pairs.add(pair)
                L = float(np.hypot(xLag[gi] - xLag[gj], yLag[gi] - yLag[gj]))
                A = edge_area.get((a, b), 0.0)
                # k = fill_stiffness_scale * E * depth * A_tributary / L^2
                # (equates discrete spring energy to continuum strain energy
                # over the edge's tributary volume A_tributary*depth)
                k = fill_stiffness_scale * E_material * extrude_depth * A / max(L ** 2, 1e-30)
                connections.append((gi, gj, 'fill'))
                fill_edge_k[pair] = k

    return connections, fill_edge_k


def _offset_polyline(pts, halfThickness, closed):
    '''
    Offsets an ordered polyline (or closed polygon, if closed=True) to
    either side by halfThickness, using per-vertex mitered normals, to
    build the two rails of a ribbon around it.

    pts: (n,2) array of centerline points, in order
    halfThickness: perpendicular offset distance for each rail
    closed: True if pts forms a closed loop (last point implicitly
        connects back to the first), False for an open chain

    Returns (railA, railB), each an (n,2) array of offset points, in the
    same order/index correspondence as pts.
    '''
    pts = np.asarray(pts, dtype=float)
    n = len(pts)
    normals = np.zeros((n, 2))

    def unit_normal(seg):
        L = np.hypot(seg[0], seg[1])
        if L < 1e-14:
            return None
        return np.array([-seg[1], seg[0]]) / L   # rotate direction 90 deg

    for i in range(n):
        if closed:
            iprev, inext = (i - 1) % n, (i + 1) % n
        else:
            iprev, inext = max(i - 1, 0), min(i + 1, n - 1)

        nPrev = unit_normal(pts[i] - pts[iprev])
        nNext = unit_normal(pts[inext] - pts[i])
        nPrev = nPrev if nPrev is not None else nNext
        nNext = nNext if nNext is not None else nPrev

        bis = nPrev + nNext
        bisNorm = np.linalg.norm(bis)
        if bisNorm < 1e-10:      # ~180 degree turn: fall back to one side's normal
            bis, miterScale = nPrev, 1.0
        else:
            bis = bis / bisNorm
            cosHalf = np.dot(bis, nPrev)
            # clamp so sharp corners (e.g. the fin tip) don't spike out too far
            miterScale = 1.0 / np.clip(cosHalf, 0.25, 1.0)

        normals[i] = bis * miterScale

    railA = pts + halfThickness * normals
    railB = pts - halfThickness * normals
    return railA, railB


# --- Legacy (single-point-thick) implementation, kept for reference -------
#
# def _Build_Tail_Geometry(ds):
#     # !Provide real tail geometry here! (same units as Lx, Ly in input2d)
#     L = 0.136   # Tail length
#     W = 0.051   # Base width
#
#     # Distance of each ray from the base, measured along the tail's length.
#     rayPositions = [0.025, 0.047, 0.069, 0.091]
#
#     # Position of the base's midpoint within the computational domain.
#     x0 = 0.3    # x-location of the base
#     y0 = 0.5    # y-location of the base's midpoint
#
#     # Half-width of the triangle at distance d from the base
#     def halfWidth(d):
#         return (W / 2) * (1 - d / L)
#
#     xLag = []
#     yLag = []
#
#     # Outline: base -> top slanted edge -> bottom slanted edge
#     # (traversed so shared corner points aren't duplicated; the loop is
#     # closed later in build_Tail_Connections)
#     # --- Base (short side), vertical segment at x = x0
#     nBase = max(round(W / ds), 2)
#     yBase = np.linspace(y0 - W / 2, y0 + W / 2, nBase)
#     xBase = x0 * np.ones_like(yBase)
#     xLag.extend(xBase)
#     yLag.extend(yBase)
#
#     # --- Top slanted edge: base-top corner -> tip
#     topStart = np.array([x0, y0 + W / 2])
#     tip = np.array([x0 + L, y0])
#     nTop = max(round(np.linalg.norm(tip - topStart) / ds), 2)
#     xTop = np.linspace(topStart[0], tip[0], nTop)
#     yTop = np.linspace(topStart[1], tip[1], nTop)
#     xLag.extend(xTop[1:])        # drop shared corner
#     yLag.extend(yTop[1:])
#
#     # --- Bottom slanted edge: tip -> base-bottom corner
#     botEnd = np.array([x0, y0 - W / 2])
#     nBot = max(round(np.linalg.norm(botEnd - tip) / ds), 2)
#     xBot = np.linspace(tip[0], botEnd[0], nBot)
#     yBot = np.linspace(tip[1], botEnd[1], nBot)
#     xLag.extend(xBot[1:-1])      # drop both shared corners
#     yLag.extend(yBot[1:-1])
#
#     outlineIdx = list(range(len(xLag)))   # every point so far belongs to the outline
#
#     # The two base corners, by global index, which are needed for actuation
#     # (target points) later. baseBottomIdx is the very first point added
#     # (start of the base segment); baseTopIdx is the last point of the
#     # base segment, i.e. index nBase-1, before the top edge begins.
#     baseBottomIdx = 0
#     baseTopIdx = nBase - 1
#     baseCornerIdx = (baseBottomIdx, baseTopIdx)
#
#     # 4 internal rays, each parallel to the base.
#     # Each ray's indices are tracked separately (rayIdxList) so that
#     # build_Tail_Connections can keep them as independent open chains
#     # instead of accidentally chaining ray -> next ray.
#     rayIdxList = []
#     for d in rayPositions:
#         hw = halfWidth(d)
#         nRay = max(round((2 * hw) / ds), 2)
#         yRay = np.linspace(y0 - hw, y0 + hw, nRay)
#         xRay = (x0 + d) * np.ones_like(yRay)
#
#         startIdx = len(xLag)
#         xLag.extend(xRay)
#         yLag.extend(yRay)
#         rayIdxList.append(list(range(startIdx, len(xLag))))
#
#     xLag = np.array(xLag)
#     yLag = np.array(yLag)
#
#     return xLag, yLag, outlineIdx, rayIdxList, baseCornerIdx
#
# ---------------------------------------------------------------------------

def _n_points_for_length(length, ds):
    '''
    Number of points to discretize an OPEN segment of the given length at
    ~ds spacing. Rounds to the nearest number of SEGMENTS (not points) and
    derives point count from that (points = segments + 1) -- for an open
    chain, n points give n-1 segments, and rounding the POINT count
    directly (as every one of base/top-slant/bottom-slant/ray used to)
    silently drops that -1: negligible for a long edge with many segments,
    but severe for a short one. E.g. a ray needing "2.14 segments" rounded
    the POINT count to 2 (best case per the old formula), i.e. exactly 1
    segment spanning the ray's FULL length -- 214% of ds, not the ~107%
    the segment-based rounding here gives. At least 1 segment (2 points)
    is always used, since a line needs at least 2 points to exist at all.
    '''
    n_segments = max(round(length / ds), 1)
    return n_segments + 1


def _Build_Tail_Geometry_Ribbon(ds, wall_thickness, fill_spacing=None,
                                 fill_density_multiplier=1.0):
    '''
    Builds the ribboned tail geometry: the same outline + 4 internal-ray
    centerlines as the legacy _Build_Tail_Geometry, except each centerline
    is offset into two parallel rails (see _offset_polyline), wall_thickness
    apart, so the structure has a real physical thickness instead of
    encoding it only through the spring/beam stiffness formulas.

    Each ribbon (the outline, and each of the 4 rays) also gets its 2D
    interior filled with Halton-sampled points (see
    _mesh_utils.fill_polygon_with_halton), so the ribbon can be meshed as a
    proper continuum strip rather than just its two rail lines -- see
    ribbonPolys below and build_Tail_Ribbon_Connections.

    fill_spacing: target interior-point spacing; None falls back to ds.
    fill_density_multiplier: scales the fill point count (0 disables fill).

    Returns:
      xLag, yLag: all Lagrangian point coordinates
      outlineRailIdx: (railA_idx, railB_idx) for the closed outline ribbon
      rayRailIdxList: list of (railA_idx, railB_idx) for each open ray ribbon
      baseCornerIdx: the 4 point indices at the base's two corners (both
          rails), used as target/actuation points
      ribbonPolys: list of dicts, one per ribbon (outline + each ray), each
          with 'kind', 'boundary_idx' (rail point indices, exterior loop
          followed by any hole loop), 'exterior', 'holes' (point arrays,
          for _mesh_utils.triangulate_with_holes), and 'interior_idx' (the
          new Halton-fill point indices)
    '''
    fill_spacing = ds if fill_spacing is None else fill_spacing
    # !Provide real tail geometry here! (same units as Lx, Ly in input2d)
    L = 0.136   # Tail length
    W = 0.051   # Base width

    # Distance of each ray from the base, measured along the tail's length.
    rayPositions = [0.025, 0.047, 0.069, 0.091]

    # Position of the base's midpoint within the computational domain.
    x0 = 0.3    # x-location of the base
    y0 = 0.5    # y-location of the base's midpoint

    halfT = wall_thickness / 2.0

    # Half-width of the triangle at distance d from the base
    def halfWidth(d):
        return (W / 2) * (1 - d / L)

    # --- Outline centerline: base -> top slanted edge -> tip -> bottom
    # slanted edge (same construction as the legacy version; it's offset
    # into two rails below instead of being used directly)
    cx, cy = [], []

    nBase = _n_points_for_length(W, ds)
    yBase = np.linspace(y0 - W / 2, y0 + W / 2, nBase)
    xBase = x0 * np.ones_like(yBase)
    cx.extend(xBase); cy.extend(yBase)

    topStart = np.array([x0, y0 + W / 2])
    tip = np.array([x0 + L, y0])
    nTop = _n_points_for_length(np.linalg.norm(tip - topStart), ds)
    xTop = np.linspace(topStart[0], tip[0], nTop)
    yTop = np.linspace(topStart[1], tip[1], nTop)
    cx.extend(xTop[1:]); cy.extend(yTop[1:])          # drop shared corner

    botEnd = np.array([x0, y0 - W / 2])
    nBot = _n_points_for_length(np.linalg.norm(botEnd - tip), ds)
    xBot = np.linspace(tip[0], botEnd[0], nBot)
    yBot = np.linspace(tip[1], botEnd[1], nBot)
    cx.extend(xBot[1:-1]); cy.extend(yBot[1:-1])      # drop both shared corners

    # Same corner bookkeeping as the legacy version, but against the
    # centerline (rail indices are derived from these below)
    baseBottomCenterIdx = 0
    baseTopCenterIdx = nBase - 1
    centerline = np.column_stack((cx, cy))

    railA, railB = _offset_polyline(centerline, halfT, closed=True)

    xLag, yLag = [], []

    def add_points(arr):
        start = len(xLag)
        xLag.extend(arr[:, 0]); yLag.extend(arr[:, 1])
        return list(range(start, len(xLag)))

    railAIdx = add_points(railA)
    railBIdx = add_points(railB)
    outlineRailIdx = (railAIdx, railBIdx)

    baseCornerIdx = (railAIdx[baseBottomCenterIdx], railBIdx[baseBottomCenterIdx],
                      railAIdx[baseTopCenterIdx], railBIdx[baseTopCenterIdx])

    ribbonPolys = []

    # The outline ribbon is a hollow "picture frame": one rail is the outer
    # loop, the other the inner (hole) loop. Which is which depends only on
    # containment, not on which offset direction _offset_polyline happened
    # to call railA/railB. Rays only ever WELD to the INNER rail (they live
    # in the hollow interior) -- but both rails are kept as plain (mutable)
    # lists and spliced IN LOCKSTEP (see weld_to_inner_rail below), so they
    # stay the same length with the same per-position correspondence that
    # add_rungs_and_diagonals relies on (rail A's k-th point assumed to sit
    # directly across the ribbon from rail B's k-th point). An earlier
    # version only spliced the inner rail, which desynced that
    # correspondence from the first splice onward -- add_rungs_and_diagonals
    # then paired far-apart points by matching (now-meaningless) array
    # index instead of physical position, producing long spurious 'rung'/
    # 'diag' springs that cut straight across the tail's open "windows".
    if _MplPath(railA, closed=True).contains_point(railB[0]):
        ext_pts, ext_idx, hole_pts, hole_idx = list(railA), railAIdx, list(railB), railBIdx
    else:
        ext_pts, ext_idx, hole_pts, hole_idx = list(railB), railBIdx, list(railA), railAIdx

    def weld_to_inner_rail(p_xy):
        '''
        Finds where p_xy (a ray rail's endpoint) actually meets the
        outline's inner rail, and welds it there: reuses the nearest
        existing inner-rail vertex if p_xy's projection lands within EPS of
        one (avoids creating a near-duplicate point on top of an existing
        one -- and the near-zero-length spring that would come with it),
        otherwise splices a brand-new point into hole_pts/hole_idx (in
        place, so it becomes a genuine member of the outline's own rail
        chain -- automatically picked up by add_rail_chain/add_rungs_and_
        diagonals below, no separate 'attach' spring needed).

        Whenever a new point IS spliced into the inner rail, a partner
        point is also spliced into the outer rail at the same segment index
        and the same parametric position t (i.e. the equivalent point on
        the corresponding outer-rail segment) -- this partner has no
        special role of its own (nothing welds to it), it exists purely to
        keep both rails the same length / positionally in sync. Since
        wall_thickness is tiny relative to the outline's feature sizes (the
        same justification used elsewhere for the offset/triangulation
        approximations), using the same t on the parallel outer segment is
        an accurate stand-in for that segment's true corresponding point.

        Returns (idx, xy) of the final, shared point -- xy may differ
        slightly from p_xy (it's snapped exactly onto the inner rail), so
        callers should use it (not the original p_xy) for the ray's own
        rail coordinates too, or the ray and outline would no longer
        actually touch.
        '''
        foot, seg_i, t = mu.project_point_to_polyline(p_xy, np.asarray(hole_pts), closed=True)
        n = len(hole_pts)
        EPS = 0.02  # fraction of the segment length
        if t <= EPS:
            return hole_idx[seg_i], hole_pts[seg_i]
        if t >= 1.0 - EPS:
            j = (seg_i + 1) % n
            return hole_idx[j], hole_pts[j]

        new_idx = add_points(foot.reshape(1, 2))[0]
        hole_pts.insert(seg_i + 1, foot)
        hole_idx.insert(seg_i + 1, new_idx)

        n_ext = len(ext_pts)
        a_ext, b_ext = ext_pts[seg_i], ext_pts[(seg_i + 1) % n_ext]
        foot_ext = a_ext + t * (b_ext - a_ext)
        new_ext_idx = add_points(foot_ext.reshape(1, 2))[0]
        ext_pts.insert(seg_i + 1, foot_ext)
        ext_idx.insert(seg_i + 1, new_ext_idx)

        return new_idx, foot

    # 4 internal rays, each parallel to the base, ribboned the same way.
    # Each ray's rail indices are tracked separately (rayRailIdxList) so
    # that build_Tail_Ribbon_Connections can keep them as independent open
    # chains instead of accidentally chaining ray -> next ray.
    #
    # Each rail's two ends are welded directly onto the outline's inner rail
    # (weld_to_inner_rail above) rather than left as free-floating points
    # tied to the nearest existing outline point by a separate 'attach'
    # spring -- that old approach left a small standoff gap (a spring only
    # resists distance changes, not angle, so with no beam spanning the
    # joint it behaved like a hinge rather than a fused connection; see
    # NOTES.md / build_Tail_Beams). Only each rail's INTERIOR points (not
    # its welded endpoints) get fresh point indices.
    rayRailIdxList = []
    for d in rayPositions:
        hw = halfWidth(d)
        nRay = _n_points_for_length(2 * hw, ds)
        yRay = np.linspace(y0 - hw, y0 + hw, nRay)
        xRay = (x0 + d) * np.ones_like(yRay)
        rayCenterline = np.column_stack((xRay, yRay))

        rA, rB = _offset_polyline(rayCenterline, halfT, closed=False)

        botA_idx, rA[0] = weld_to_inner_rail(rA[0])
        topA_idx, rA[-1] = weld_to_inner_rail(rA[-1])
        botB_idx, rB[0] = weld_to_inner_rail(rB[0])
        topB_idx, rB[-1] = weld_to_inner_rail(rB[-1])

        rA_interior_idx = add_points(rA[1:-1]) if nRay > 2 else []
        rB_interior_idx = add_points(rB[1:-1]) if nRay > 2 else []
        rAIdx = [botA_idx] + rA_interior_idx + [topA_idx]
        rBIdx = [botB_idx] + rB_interior_idx + [topB_idx]
        rayRailIdxList.append((rAIdx, rBIdx))

        # Close the open two-rail strip into a simple polygon (up one rail,
        # back down the other) so it can be Halton-filled/triangulated.
        ring_pts = np.vstack([rA, rB[::-1]])
        ring_idx = rAIdx + rBIdx[::-1]
        ray_interior_xy = mu.fill_polygon_with_halton(
            ring_pts, holes=None, spacing=fill_spacing,
            density_multiplier=fill_density_multiplier)
        ray_interior_idx = add_points(ray_interior_xy)
        ribbonPolys.append(dict(kind='ray', boundary_idx=ring_idx,
                                 exterior=ring_pts, holes=[],
                                 interior_idx=ray_interior_idx))

    # The outline's own Halton fill happens LAST, using the final
    # hole_pts/hole_idx (after every ray has welded its endpoints in) --
    # otherwise its interior mesh would be triangulated against a boundary
    # that's missing the very junction points the rays just added.
    outline_interior_xy = mu.fill_polygon_with_halton(
        ext_pts, holes=[hole_pts], spacing=fill_spacing,
        density_multiplier=fill_density_multiplier)
    outline_interior_idx = add_points(outline_interior_xy)
    ribbonPolys.append(dict(kind='outline', boundary_idx=ext_idx + hole_idx,
                             exterior=ext_pts, holes=[hole_pts],
                             interior_idx=outline_interior_idx))

    xLag = np.array(xLag)
    yLag = np.array(yLag)

    return xLag, yLag, outlineRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys


if __name__ == "__main__":
    FinRay_Geom()
