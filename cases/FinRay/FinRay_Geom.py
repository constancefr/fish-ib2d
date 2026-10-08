"""
Creates the Lagrangian structure geometry for a finray tail and writes the
associated input files (.vertex, .geo_connect, .spring, .target).

Structure model: each originally-1-point-thick edge (the frame + the 4
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

_COMMON_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'common')
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)
import mesh_utils as mu

_case_params = {}
if os.path.exists('case_params.json'):
    with open('case_params.json') as f:
        _case_params = json.load(f)
_geom = _case_params.get('geom', {})


def FinRay_Geom():
    #-------------------------------------------------------------------
    # Eulerian grid params (MAKE SURE MATCHES THOSE IN input2d !!!)
    Nx = 32           # MUST BE EVEN!!!
    Lx = 1.0          # length of grid in x-direction
    dx = Lx / Nx      # spatial resolution

    #-------------------------------------------------------------------
    # Immersed structure params
    ds = 0.5 * dx               # maximum Lagrangian point spacing (2x resolution of Eulerian grid!!!)
    struct_name = 'finray'    # name for .vertex, .spring, etc files (must match what's in 'input2d')

    #-------------------------------------------------------------------
    # Material properties - trying to mimick real system from dimensions and Dragon Skin 30 properties
    E_material = _geom.get('E_material', 0.74e6)
    wall_thickness = _geom.get('wall_thickness', 0.00396)   # m (3.96mm in real build)
    extrude_depth = _geom.get('extrude_depth', 0.05)        # m

    #-------------------------------------------------------------------
    # Halton sampling params - see fill_polygon_with_halton() + triangulate_with_holes() 
    # in common/mesh_utils.py
    fill_spacing = _geom.get('ribbon_fill_spacing', ds) # default to ds
        # note: fill_spacing doesn't have to be >=ds, as it isn't the boundary
    fill_density_multiplier = _geom.get('ribbon_fill_density_multiplier', 1.0) 
        # scales number of interior points (0 = no fill, 1 = same density as ds, >1 = denser)
        # TODO: what does it actually scale? fill_spacing? is 1 really the same density as ds??
    stiffness_scale = _geom.get('ribbon_fill_stiffness_scale', 1.0)
        # scales the stiffness of EVERY spring (rails included) -- the single
        # calibration multiplier alpha (1.0 = the formula's literal value,
        # 0.5 = half as stiff, etc.). The key keeps its old name so existing
        # case_params.json / run_case.py sweeps still work. See NOTES.md,
        # "single stiffness formula".

    #-------------------------------------------------------------------
    # Construct geometry
    xLag, yLag, frameRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys = \
        _Build_Tail_Geometry_Ribbon(ds, wall_thickness,
                                     fill_spacing=fill_spacing,
                                     fill_density_multiplier=fill_density_multiplier)
    #-------------------------------------------------------------------
    # Build springs
    connections, edge_k = build_Tail_Ribbon_Connections(
        xLag, yLag, frameRailIdx, rayRailIdxList, ribbonPolys,
        E_material=E_material, extrude_depth=extrude_depth,
        stiffness_scale=stiffness_scale)

    #-------------------------------------------------------------------
    # Set actuation / target stiffness
    actuation_patch_radius = _geom.get('actuation_patch_radius', None)
    if actuation_patch_radius is None:
        actuation_patch_radius = 1.0 * wall_thickness
    target_indices = _actuation_patch_indices(xLag, yLag, baseCornerIdx, actuation_patch_radius,
                                               ribbonPolys=ribbonPolys)

    target_stiffness_factor = _geom.get('target_stiffness_factor', 15.0)
    k_Target = target_stiffness_factor * _max_incident_stiffness(target_indices, edge_k)

    #-------------------------------------------------------------------
    # Plot geometry to test
    plot_colors = {'frame': 'b', 'ray': 'b', 'fill': (0.6, 0.6, 0.6, 0.6)}
    for (i, j, kind) in connections:
        lw = 0.4 if kind == 'fill' else 0.75
        plt.plot([xLag[i], xLag[j]], [yLag[i], yLag[j]],
                  color=plot_colors.get(kind, 'k'), linewidth=lw)
    for rail in frameRailIdx:
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

    # Write to IB2d input files
    print_Lagrangian_Vertices(xLag, yLag, struct_name) # .vertex
    print_Geometry_Connections(connections, struct_name) # .geo_connect
    print_Lagrangian_Springs(xLag, yLag, connections, edge_k, struct_name) # .spring
    print_Lagrangian_Target_Pts(target_indices, k_Target, struct_name) # .target



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

def print_Lagrangian_Springs(xLag, yLag, connections, edge_k, struct_name, deg_NL=1.0):
    '''
    Prints springs to .spring file.

    connections: list of (i, j, kind) generated by build_Tail_Ribbon_Connections()
    edge_k: dict mapping (min(i,j), max(i,j)) -> stiffness for every
        connection, computed per-edge in build_Tail_Ribbon_Connections
        (tributary-area formula, the same for every kind).
    struct_name: name of the structure
    deg_NL: degree of nonlinearity
    '''
    N = len(connections)

    with open(struct_name + '.spring', 'w') as f:
        f.write('%d\n' % N)   # Print # of springs

        for (i, j, kind) in connections:
            ds_Rest = float(np.hypot(xLag[j] - xLag[i], yLag[j] - yLag[i])) # set initial distance as resting length
            k = edge_k[(min(i, j), max(i, j))]
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

def _max_incident_stiffness(point_indices, edge_k):
    '''
    Largest total spring stiffness attached to any one of point_indices,
    i.e. sum of k over every spring touching that point, maxed over points.
    Used as the reference magnitude for k_Target.

    edge_k: dict (min(i,j), max(i,j)) -> k, from build_Tail_Ribbon_Connections.
    '''
    incident = {p: 0.0 for p in point_indices}
    for (i, j), k in edge_k.items():
        if i in incident:
            incident[i] += k
        if j in incident:
            incident[j] += k
    return max(incident.values())


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

def build_Tail_Ribbon_Connections(xLag, yLag, frameRailIdx, rayRailIdxList,
                                   ribbonPolys, E_material, extrude_depth,
                                   stiffness_scale=1.0):
    '''
    Builds every spring of the ribboned structure from ONE triangulated mesh
    per ribbon, and gives every spring the same tributary-area stiffness:

        k = stiffness_scale * E * depth * A_tributary / L**2

    A_tributary is 1/3 of the area of each triangle sharing that edge
    (summed over every ribbon the edge belongs to -- e.g. a frame inner-rail
    segment that is also a ray's welded end). This equates the springs'
    stored energy to the continuum strain energy of the triangles they
    span, so the springs together represent the wall's material exactly
    once, at any sampling density (see NOTES.md, "single stiffness
    formula"). Consequences:
      - a rail (edge) spring touches triangles on one side only, so it
        gets roughly half the stiffness of a similar interior spring;
      - stiffness_scale scales ALL springs (it is the single calibration
        multiplier alpha), not just the interior ones.

    kind (for plotting / bookkeeping only -- it no longer changes k):
      'frame': along each of the frame's two rails (closed loop).
      'ray':   along each ray's two rails (open chain).
      'fill':  every other triangulation edge (cross-wall and interior).

    ribbonPolys: list of dicts (see _Build_Tail_Geometry_Ribbon), one per
        ribbon, each with 'boundary_idx'/'exterior_loop'/'inner_loops'/
        'interior_idx'. Each ribbon is triangulated even with no interior
        points (density multiplier 0), so its cross-wall edges still exist.

    Returns (connections, edge_k): connections is a list of (i, j, kind);
    edge_k maps each connection's (min(i,j), max(i,j)) pair to its k.
    '''

    connections = []

    def add_rail_chain(railIdx, closed, kind):
        n = len(railIdx)
        stop = n if closed else n - 1
        for k in range(stop):
            connections.append((railIdx[k], railIdx[(k + 1) % n], kind))

    for railIdx in frameRailIdx:
        add_rail_chain(railIdx, True, 'frame')
    for (rA, rB) in rayRailIdxList:
        add_rail_chain(rA, False, 'ray')
        add_rail_chain(rB, False, 'ray')

    # Ray rail endpoints are welded into the frame's inner rail (see
    # weld_to_inner_rail in _Build_Tail_Geometry_Ribbon), so the same
    # vertex pair can appear twice, or two endpoints can weld onto one
    # vertex. Drop self-loops (zero-length, divide-by-zero springs) and
    # duplicate pairs, keeping the first occurrence's kind.
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

    # Triangulate every ribbon and accumulate each edge's tributary area
    # (1/3 of every triangle sharing it), keyed by GLOBAL point pair so an
    # edge shared by two ribbons collects area from both.
    edge_area = {}
    for ribbon in ribbonPolys:
        b_idx, i_idx = ribbon['boundary_idx'], ribbon['interior_idx']
        boundary_pts = np.column_stack((xLag[b_idx], yLag[b_idx]))
        interior_pts = np.column_stack((xLag[i_idx], yLag[i_idx])) if i_idx else np.empty((0, 2))
        all_pts, _edges, tris = mu.triangulate_with_holes(
            ribbon['exterior_loop'], ribbon['inner_loops'], boundary_pts, interior_pts)
        local_to_global = list(b_idx) + list(i_idx)
        for simplex in tris:
            p0, p1, p2 = all_pts[simplex]
            area = 0.5 * abs((p1[0] - p0[0]) * (p2[1] - p0[1])
                              - (p2[0] - p0[0]) * (p1[1] - p0[1]))
            g = [local_to_global[v] for v in simplex]
            for a, b in ((g[0], g[1]), (g[1], g[2]), (g[2], g[0])):
                if a == b:
                    continue
                pair = (min(a, b), max(a, b))
                edge_area[pair] = edge_area.get(pair, 0.0) + area / 3.0

    # Every rail segment must be a triangulation edge (triangulate_with_holes
    # forces boundary segments in); one that isn't would get no material.
    missing = [(i, j) for (i, j, _k) in connections if (min(i, j), max(i, j)) not in edge_area]
    if missing:
        raise RuntimeError(f"{len(missing)} rail segment(s) are not edges of the ribbon "
                           f"triangulation, e.g. {missing[:3]} -- they would carry no material")

    for pair in sorted(edge_area):
        if pair not in seen_pairs:
            seen_pairs.add(pair)
            connections.append((pair[0], pair[1], 'fill'))

    edge_k = {}
    for pair, A in edge_area.items():
        L = float(np.hypot(xLag[pair[0]] - xLag[pair[1]], yLag[pair[0]] - yLag[pair[1]]))
        edge_k[pair] = stiffness_scale * E_material * extrude_depth * A / max(L ** 2, 1e-30)

    return connections, edge_k


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

def _n_points_for_length(length, ds):
    '''
    Number of points to discretize an OPEN segment of the given length at
    ~ds spacing. Rounds to the nearest number of SEGMENTS (not points) and
    derives point count from that (points = segments + 1).
    At least 1 segment (2 points) since a line needs at least 2 points to exist.
    '''
    n_segments = max(round(length / ds), 1)
    return n_segments + 1


def _Build_Tail_Geometry_Ribbon(ds, wall_thickness, fill_spacing=None,
                                 fill_density_multiplier=1.0):
    '''
    Builds tail geometry with Halton sampling.
    The structure is currently hardcoded with 4 evenly spaced inner rays.
        TODO: implement variable structure parameters (e.g. # rays, length, base width)

    Input:
        ds: Lagrangian point spacing (m)
        wall_thickness: ribbon wall thickness (m)
        fill_spacing: interior-point spacing; None falls back to ds.
        fill_density_multiplier: scales the fill point count (0 disables fill).

    Returns:
        xLag, yLag: 1D arrays of all Lagrangian point coordinates
        frameRailIdx: two index lists (railA_idx, railB_idx) - one for each edge of the triangular frame
        rayRailIdxList: list of index lists (railA_idx, railB_idx) for each ray ribbon
        baseCornerIdx: 4 point indices at the base's two corners as seeds for actuation (two corners x two rails)
        ribbonPolys: list of dicts, one per ribbon (frame + each ray), each with:
            - 'kind'
            - 'boundary_idx' (rail point indices, i.e. exterior loop followed by any inner loop)
            - 'exterior_loop'
            - 'inner_loops' (point arrays of inner loops i.e. around holes, for mesh_utils.triangulate_with_holes)
            - 'interior_idx' (the new Halton fill point indices)

    NOTE: the tail is placed with its base vertical and the tip pointing right in the domain.
    '''
    fill_spacing = ds if fill_spacing is None else fill_spacing
    # TODO: Provide real tail geometry here! (same units as Lx, Ly in input2d)
    L = 0.136   # tail length
    W = 0.051   # base width

    # Distance of each ray from the base, measured along the tail's length
    rayPositions = [0.025, 0.047, 0.069, 0.091]

    # Position of the base's midpoint within the computatioxnal domain
    x0 = 0.3
    y0 = 0.5

    halfT = wall_thickness / 2.0

    # Half-width of the triangle at distance d from the base
    def halfWidth(d):
        return (W / 2) * (1 - d / L)

    # Outline centerline
    cx, cy = [], []
    # Base
    nBase = _n_points_for_length(W, ds)
    yBase = np.linspace(y0 - W / 2, y0 + W / 2, nBase) # evenly space points vertically
    xBase = x0 * np.ones_like(yBase)
    cx.extend(xBase); cy.extend(yBase)
    # Top slanted edge: base-top corner -> tip
    topStart = np.array([x0, y0 + W / 2])
    tip = np.array([x0 + L, y0])
    nTop = _n_points_for_length(np.linalg.norm(tip - topStart), ds)
    xTop = np.linspace(topStart[0], tip[0], nTop)
    yTop = np.linspace(topStart[1], tip[1], nTop)
    cx.extend(xTop[1:]); cy.extend(yTop[1:])          # drop shared corner
    # Bottom slanted edge: tip -> base-bottom corner
    botEnd = np.array([x0, y0 - W / 2])
    nBot = _n_points_for_length(np.linalg.norm(botEnd - tip), ds)
    xBot = np.linspace(tip[0], botEnd[0], nBot)
    yBot = np.linspace(tip[1], botEnd[1], nBot)
    cx.extend(xBot[1:-1]); cy.extend(yBot[1:-1])      # drop both shared corners

    # Corner bookkeeping against the centerline (rail indices are derived from these below)
    baseBottomCenterIdx = 0
    baseTopCenterIdx = nBase - 1
    centerline = np.column_stack((cx, cy))

    # Build to rails by offsetting the centreline - doesn't determine which is exterior vs interior
    railA, railB = _offset_polyline(centerline, halfT, closed=True)

    xLag, yLag = [], []

    def add_points(arr):
        start = len(xLag)
        xLag.extend(arr[:, 0]); yLag.extend(arr[:, 1])
        return list(range(start, len(xLag)))

    railAIdx = add_points(railA)
    railBIdx = add_points(railB)
    frameRailIdx = (railAIdx, railBIdx)

    baseCornerIdx = (railAIdx[baseBottomCenterIdx], railBIdx[baseBottomCenterIdx],
                      railAIdx[baseTopCenterIdx], railBIdx[baseTopCenterIdx])

    ribbonPolys = []

    # Indentify external vs interal rail
    if mu.closed_path(railA).contains_point(railB[0]): # railB[0] lies inside area of railA
        ext_pts, ext_idx, inner_pts, inner_idx = list(railA), railAIdx, list(railB), railBIdx
    else:
        ext_pts, ext_idx, inner_pts, inner_idx = list(railB), railBIdx, list(railA), railAIdx

    def weld_to_inner_rail(p_xy):
        '''
        Finds where p_xy (a ray rail's endpoint) meets the frame's inner rail and welds it there.
        Reuses the nearest existing inner-rail vertex if p_xy's projection lands within EPS of
        one (avoids creating a near-duplicate point on top of an existing one and a near-zero-length spring).
        Otherwise splices a new point into inner_pts/inner_idx.

        Returns (idx, xy) of the shared point.
        NOTE: xy may differ slightly from p_xy so callers should use xy for the ray's own
        rail coordinates too, or the ray and frame would no longer actually touch.
        '''
        foot, seg_i, t = mu.project_point_to_polyline(p_xy, np.asarray(inner_pts), closed=True)
        n = len(inner_pts)
        EPS = 0.02  # fraction of the segment length
        if t <= EPS:
            return inner_idx[seg_i], inner_pts[seg_i]
        if t >= 1.0 - EPS:
            j = (seg_i + 1) % n
            return inner_idx[j], inner_pts[j]

        new_idx = add_points(foot.reshape(1, 2))[0]
        inner_pts.insert(seg_i + 1, foot)
        inner_idx.insert(seg_i + 1, new_idx)

        n_ext = len(ext_pts)
        a_ext, b_ext = ext_pts[seg_i], ext_pts[(seg_i + 1) % n_ext]
        foot_ext = a_ext + t * (b_ext - a_ext)
        new_ext_idx = add_points(foot_ext.reshape(1, 2))[0]
        ext_pts.insert(seg_i + 1, foot_ext)
        ext_idx.insert(seg_i + 1, new_ext_idx)

        return new_idx, foot

    # 4 internal rays parallel to the base.
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

        # Close the open two-rail strip into a simple polygon
        ring_pts = np.vstack([rA, rB[::-1]])
        ring_idx = rAIdx + rBIdx[::-1]
        
        # Halton sampling of ray's interior
        ray_interior_xy = mu.fill_polygon_with_halton(
            ring_pts, inner_loops=None, spacing=fill_spacing,
            density_multiplier=fill_density_multiplier)
        ray_interior_idx = add_points(ray_interior_xy)
        ribbonPolys.append(dict(kind='ray', boundary_idx=ring_idx,
                                 exterior_loop=ring_pts, inner_loops=[],
                                 interior_idx=ray_interior_idx))

    # Halton sampling of frame's interior
    frame_interior_xy = mu.fill_polygon_with_halton(
        ext_pts, inner_loops=[inner_pts], spacing=fill_spacing,
        density_multiplier=fill_density_multiplier)
    frame_interior_idx = add_points(frame_interior_xy)
    ribbonPolys.append(dict(kind='frame', boundary_idx=ext_idx + inner_idx,
                             exterior_loop=ext_pts, inner_loops=[inner_pts],
                             interior_idx=frame_interior_idx))

    xLag = np.array(xLag)
    yLag = np.array(yLag)

    return xLag, yLag, frameRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys


if __name__ == "__main__":
    FinRay_Geom()
