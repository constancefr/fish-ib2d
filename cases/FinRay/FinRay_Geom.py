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
import json, os

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

    # Construct geometry: each edge/ray is now a 2-rail ribbon of width
    # wall_thickness (instead of a single 1-point-thick line). Also returns
    # which rail indices belong to the outline vs. each ray, needed to wire
    # up connections correctly below.
    xLag, yLag, outlineRailIdx, rayRailIdxList, baseCornerIdx = \
        _Build_Tail_Geometry_Ribbon(ds, wall_thickness)

    # Build connections: rail springs ('outline'/'ray'), cross-thickness
    # 'rung' + 'diag' springs that give each ribbon its shear stiffness, and
    # 'attach' springs tying each ray's rail-ends to the outline
    connections = build_Tail_Ribbon_Connections(xLag, yLag, outlineRailIdx, rayRailIdxList)

    # Plot geometry to test
    plot_colors = {'outline': 'b-', 'ray': 'b-', 'rung': 'c-', 'diag': 'y-', 'attach': 'm-'}
    for (i, j, kind) in connections:
        plt.plot([xLag[i], xLag[j]], [yLag[i], yLag[j]], plot_colors.get(kind, 'k-'), linewidth=0.75)
    for rail in outlineRailIdx:
        plt.plot(xLag[rail], yLag[rail], 'r*')
    for railPair in rayRailIdxList:
        for rail in railPair:
            plt.plot(xLag[rail], yLag[rail], 'g*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('square')
    plt.show()

    # Print .vertex file
    print_Lagrangian_Vertices(xLag, yLag, struct_name)

    # Print .geo_connect file
    print_Geometry_Connections(connections, struct_name)

    # Print .spring file
    #   k_Spring: per-rail axial stiffness. Each ribbon is idealized as a
    #     2-flange truss approximating a solid EI beam, so each rail only
    #     carries half the wall's cross-sectional area (A/2), separated by
    #     ~wall_thickness from its partner rail.
    #   k_Rung:   cross-thickness ('rung' + 'diag') springs approximate the
    #     "web" of that idealized beam, i.e. the material's shear
    #     resistance across the wall thickness.
    #   k_Attach: ray-to-outline junction springs - play around to see how
    #     tail flap is affected by making these stiffer/softer than the
    #     material itself.
    k_Spring = spring_Stiffness_From_Material(E_material, wall_thickness / 2.0, extrude_depth, ds)
    k_Rung = shear_Stiffness_From_Material(E_material, wall_thickness, extrude_depth, ds)
    k_Attach = 50.0 * k_Spring
    k_by_kind = {'outline': k_Spring, 'ray': k_Spring, 'rung': k_Rung, 'diag': k_Rung, 'attach': k_Attach}
    print_Lagrangian_Springs(xLag, yLag, connections, k_by_kind, struct_name)

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
    k_Target = 1000.0 * k_Spring
    print_Lagrangian_Target_Pts(list(baseCornerIdx), k_Target, struct_name)



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

def print_Lagrangian_Springs(xLag, yLag, connections, k_by_kind, struct_name, deg_NL=1.0):
    '''
    Prints springs to .spring file.

    connections: list of (i, j, kind) generated by build_Tail_Ribbon_Connections()
    k_by_kind: dict mapping each connection 'kind' (see build_Tail_Ribbon_Connections)
        to its spring constant
    struct_name: name of the structure
    deg_NL: degree of nonlinearity
    '''
    N = len(connections)

    with open(struct_name + '.spring', 'w') as f:
        f.write('%d\n' % N)   # Print # of springs

        for (i, j, kind) in connections:
            ds_Rest = float(np.hypot(xLag[j] - xLag[i], yLag[j] - yLag[i])) # set initial distance as resting length
            k = k_by_kind[kind]
            f.write('%d %d %1.16e %1.16e %1.16e\n' % (i, j, k, ds_Rest, deg_NL))

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
    xLag, yLag, _, _, _ = _Build_Tail_Geometry_Ribbon(ds, wall_thickness)
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

def build_Tail_Ribbon_Connections(xLag, yLag, outlineRailIdx, rayRailIdxList):
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
      'attach':  each ray's two rail-ends connected to their nearest
                 outline rail point.

    kind is returned explicitly so the caller can give each a different
    stiffness (e.g. 'attach' stiffer/softer than the material itself).
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

    # Attach each ray's rail-ends to their nearest outline rail point
    outlineIdx = list(outlineRailIdx[0]) + list(outlineRailIdx[1])
    outline_xy = np.column_stack((xLag[outlineIdx], yLag[outlineIdx]))
    for (rA, rB) in rayRailIdxList:
        for endpoint in (rA[0], rB[0], rA[-1], rB[-1]):
            p = np.array([xLag[endpoint], yLag[endpoint]])
            dists = np.linalg.norm(outline_xy - p, axis=1)
            nearest_idx = outlineIdx[int(np.argmin(dists))]
            connections.append((endpoint, nearest_idx, 'attach'))

    return connections


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

def _Build_Tail_Geometry_Ribbon(ds, wall_thickness):
    '''
    Builds the ribboned tail geometry: the same outline + 4 internal-ray
    centerlines as the legacy _Build_Tail_Geometry, except each centerline
    is offset into two parallel rails (see _offset_polyline), wall_thickness
    apart, so the structure has a real physical thickness instead of
    encoding it only through the spring/beam stiffness formulas.

    Returns:
      xLag, yLag: all Lagrangian point coordinates
      outlineRailIdx: (railA_idx, railB_idx) for the closed outline ribbon
      rayRailIdxList: list of (railA_idx, railB_idx) for each open ray ribbon
      baseCornerIdx: the 4 point indices at the base's two corners (both
          rails), used as target/actuation points
    '''
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

    nBase = max(round(W / ds), 2)
    yBase = np.linspace(y0 - W / 2, y0 + W / 2, nBase)
    xBase = x0 * np.ones_like(yBase)
    cx.extend(xBase); cy.extend(yBase)

    topStart = np.array([x0, y0 + W / 2])
    tip = np.array([x0 + L, y0])
    nTop = max(round(np.linalg.norm(tip - topStart) / ds), 2)
    xTop = np.linspace(topStart[0], tip[0], nTop)
    yTop = np.linspace(topStart[1], tip[1], nTop)
    cx.extend(xTop[1:]); cy.extend(yTop[1:])          # drop shared corner

    botEnd = np.array([x0, y0 - W / 2])
    nBot = max(round(np.linalg.norm(botEnd - tip) / ds), 2)
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

    # 4 internal rays, each parallel to the base, ribboned the same way.
    # Each ray's rail indices are tracked separately (rayRailIdxList) so
    # that build_Tail_Ribbon_Connections can keep them as independent open
    # chains instead of accidentally chaining ray -> next ray.
    rayRailIdxList = []
    for d in rayPositions:
        hw = halfWidth(d)
        nRay = max(round((2 * hw) / ds), 2)
        yRay = np.linspace(y0 - hw, y0 + hw, nRay)
        xRay = (x0 + d) * np.ones_like(yRay)
        rayCenterline = np.column_stack((xRay, yRay))

        rA, rB = _offset_polyline(rayCenterline, halfT, closed=False)
        rAIdx = add_points(rA)
        rBIdx = add_points(rB)
        rayRailIdxList.append((rAIdx, rBIdx))

    xLag = np.array(xLag)
    yLag = np.array(yLag)

    return xLag, yLag, outlineRailIdx, rayRailIdxList, baseCornerIdx


if __name__ == "__main__":
    FinRay_Geom()
