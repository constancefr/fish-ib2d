"""
Creates the Lagrangian structure geometry for a finray tail and writes the
associated input files (.vertex, .geo_connect, .spring, .target).

Structure model: the tail (triangular frame + 4 internal rays, every wall
wall_thickness wide) is ONE polygon with 5 holes (the open windows between
base, rays and tip). It is meshed once: points along every edge, Halton fill
inside, and one constrained Delaunay triangulation (_Build_Tail_Mesh). Every
triangle edge is a spring with the tributary-area stiffness
(build_Tail_Connections), so bending stiffness comes from the wall's material
itself (one side stretches while the other compresses) rather than from a
separate bending spring. The legacy single-point-thick model is in
FinRay_Geom_zero_thickness.py.

IB2d unit convention (important, see common/ib_units.py): IB2d is 2D, so
every stiffness is per unit depth (span = 1 m), and the driver spreads
`f * ds` with the global ds = Lx / (2 Nx). Spring constants are therefore
computed per unit depth in SI units and divided by ds via to_ib(). The real
fin's depth (geom.extrude_depth) does NOT enter the simulation: it is only
for turning per-depth forces into rough 3D forces (F_3D ~ F_2D * depth).
"""

import numpy as np
import matplotlib.pyplot as plt
import json, os, sys, time

_COMMON_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'common')
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)
import mesh_utils as mu
from ib_units import ib_ds, to_ib

_case_params = {}
if os.path.exists('case_params.json'):
    with open('case_params.json') as f:
        _case_params = json.load(f)
_geom = _case_params.get('geom', {})


def FinRay_Geom():
    #-------------------------------------------------------------------
    # Eulerian grid params, from case_params.json["input2d"] (run_case.py
    # writes the same values into input2d); defaults match input2d.
    _grid = _case_params.get('input2d', {})
    Nx = int(_grid.get('Nx', 32))     # MUST BE EVEN!!!
    Lx = float(_grid.get('Lx', 1.0))  # length of grid in x-direction

    #-------------------------------------------------------------------
    # Immersed structure params
    ds = ib_ds(Lx, Nx)          # = dx/2: boundary point spacing AND IB2d's force-spreading length
    struct_name = 'finray'    # name for .vertex, .spring, etc files (must match what's in 'input2d')

    #-------------------------------------------------------------------
    # Material properties - trying to mimick real system from dimensions and Dragon Skin 30 properties
    E_material = _geom.get('E_material', 0.74e6)
    wall_thickness = _geom.get('wall_thickness', 0.00396)   # m (3.96mm in real build)
    # extrude_depth (geom, real fin depth in m, 0.027) is not read here: the 2D
    # model is per unit depth. Use it only to scale output forces to 3D.

    #-------------------------------------------------------------------
    # Halton sampling params - see fill_polygon_with_halton() + triangulate_with_holes() 
    # in common/mesh_utils.py
    fill_spacing = _geom.get('ribbon_fill_spacing') or ds # default to ds (the json holds null)
        # note: fill_spacing doesn't have to be >=ds, as it isn't the boundary
    fill_density_multiplier = _geom.get('ribbon_fill_density_multiplier', 1.0) 
        # scales number of interior points (0 = no fill, 1 = same density as ds, >1 = denser)
        # TODO: what does it actually scale? fill_spacing? is 1 really the same density as ds??
    stiffness_scale = _geom.get('ribbon_fill_stiffness_scale', 1.0)
        # scales the stiffness of EVERY spring (outline edges included) -- the single
        # calibration multiplier alpha (1.0 = the formula's literal value,
        # 0.5 = half as stiff, etc.). The key keeps its old name so existing
        # case_params.json / run_case.py sweeps still work. See NOTES.md,
        # "single stiffness formula".


    #-------------------------------------------------------------------
    # Build the mesh: one outline, one fill, one triangulation
    t_start = time.perf_counter()
    xLag, yLag, triangles, n_boundary, baseCornerIdx, outer, windows = \
        _Build_Tail_Mesh(ds, wall_thickness, fill_spacing, fill_density_multiplier)

    #-------------------------------------------------------------------
    # Build springs
    connections, edge_k = build_Tail_Connections(
        xLag, yLag, triangles, E_material=E_material, ds=ds, stiffness_scale=stiffness_scale)

    #-------------------------------------------------------------------
    # Set actuation / target stiffness
    actuation_patch_radius = _geom.get('actuation_patch_radius', None)
    if actuation_patch_radius is None:
        actuation_patch_radius = 1.0 * wall_thickness
    target_indices = _actuation_patch_indices(xLag, yLag, baseCornerIdx, actuation_patch_radius)

    target_stiffness_factor = _geom.get('target_stiffness_factor', 15.0)
    k_Target = target_stiffness_factor * _max_incident_stiffness(target_indices, edge_k)
    build_time = time.perf_counter() - t_start

    #-------------------------------------------------------------------
    # Mesh metrics for the run record (run_case.py copies them into run_info.json)
    tail_area = mu._shoelace_area(outer) - sum(mu._shoelace_area(w) for w in windows)
    x, y = xLag[triangles], yLag[triangles]
    mesh_area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
                             - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])).sum()
    k = np.array(list(edge_k.values()))
    mesh_info = dict(
        ds=ds, fill_spacing=fill_spacing, fill_density_multiplier=fill_density_multiplier,
        edge_spacing=fill_spacing / np.sqrt(max(fill_density_multiplier, 1.0)),
        n_edge_points=n_boundary,
        n_fill_points_requested=int(round(fill_density_multiplier * tail_area / fill_spacing ** 2)),
        n_fill_points=len(xLag) - n_boundary,
        n_points=len(xLag), n_triangles=len(triangles), n_springs=len(connections),
        k_median=float(np.median(k)), k_max=float(k.max()),
        n_target_points=len(target_indices), k_Target=k_Target,
        tail_area_m2=tail_area, mesh_area_fraction=mesh_area / tail_area,
        build_time_s=round(build_time, 3),
    )
    with open('mesh_info.json', 'w') as f:
        json.dump(mesh_info, f, indent=2)

    #-------------------------------------------------------------------
    # Plot geometry to test
    plot_colors = {'boundary': 'b', 'fill': (0.6, 0.6, 0.6, 0.6)}
    for (i, j, kind) in connections:
        lw = 0.4 if kind == 'fill' else 0.75
        plt.plot([xLag[i], xLag[j]], [yLag[i], yLag[j]],
                  color=plot_colors[kind], linewidth=lw)
    plt.plot(xLag[n_boundary:], yLag[n_boundary:], 'k.', markersize=2)  # fill points
    plt.plot(xLag[target_indices], yLag[target_indices], 'r*')          # actuation patches
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

    connections: list of (i, j, kind) generated by build_Tail_Connections()
    '''
    with open(struct_name + '.geo_connect', 'w') as f:
        for (i, j, kind) in connections:
            s1, s2 = i, j
            f.write('%d %d\n' % (s1, s2))
            f.write('%d %d\n' % (s2, s1))

def print_Lagrangian_Springs(xLag, yLag, connections, edge_k, struct_name, deg_NL=1.0):
    '''
    Prints springs to .spring file.

    connections: list of (i, j, kind) generated by build_Tail_Connections()
    edge_k: dict mapping (min(i,j), max(i,j)) -> stiffness for every
        connection, computed per-edge in build_Tail_Connections
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


def _actuation_patch_indices(xLag, yLag, baseCornerIdx, radius):
    '''
    The two actuation patches, matching the real build: the tips of two
    linear actuators are attached to the base's OUTER face, one at each end;
    each pushes then pulls (+x and back) in turn and holds its corner in
    place while the other one moves. A patch is every point on that face
    (x = the outer corners' x) within `radius` of its outer corner, corner
    included: the attachment footprint. No window (inner face), slanted-edge
    or fill point is driven. update_Target_Point_Positions.py moves each
    patch as one rigid block (stiff target springs, so it pulls and holds as
    well as pushes).

    baseCornerIdx: (outer_bottom, outer_top), as returned by _Build_Tail_Mesh.
    radius: attachment reach along the base face (m). Should stay well under half
        the base width, or the two patches could merge.

    Returns a sorted list of point indices (the union of both patches).
    '''
    pts = np.column_stack((xLag, yLag))
    on_face = np.isclose(xLag, xLag[baseCornerIdx[0]], rtol=0.0, atol=1e-12)
    near = np.zeros(len(pts), dtype=bool)
    for c in baseCornerIdx:
        near |= np.linalg.norm(pts - pts[c], axis=1) <= radius
    return sorted(np.where(on_face & near)[0].tolist())


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

    edge_k: dict (min(i,j), max(i,j)) -> k, from build_Tail_Connections.
    '''
    incident = {p: 0.0 for p in point_indices}
    for (i, j), k in edge_k.items():
        if i in incident:
            incident[i] += k
        if j in incident:
            incident[j] += k
    return max(incident.values())

def build_Tail_Connections(xLag, yLag, triangles, E_material, ds, stiffness_scale=1.0):
    '''
    One spring per edge of the tail's triangulation, each with the same
    tributary-area stiffness:

        k = to_ib(stiffness_scale * E * A_tributary / L**2, ds)
          = stiffness_scale * E * A_tributary / (L**2 * ds)

    i.e. a physical stiffness per unit depth [N/m per m], divided by IB2d's
    force-spreading length ds (see module docstring / common/ib_units.py).

    A_tributary is 1/3 of the area of each triangle sharing that edge. This
    equates the springs' stored energy to the continuum strain energy of the
    triangles they span, so the springs represent the wall's material exactly
    once, at any sampling density (see NOTES.md, "single stiffness formula").
    An outline edge borders one triangle only, so it gets roughly half the
    stiffness of a similar interior edge. stiffness_scale scales ALL springs
    (the single calibration multiplier alpha).

    kind (plotting only): 'boundary' for outline edges (one triangle),
    'fill' for every other edge.

    Returns (connections, edge_k): connections is a list of (i, j, kind);
    edge_k maps each connection's (min(i,j), max(i,j)) pair to its k.
    '''
    edge_area, edge_ntri = {}, {}
    for tri in triangles:
        x, y = xLag[tri], yLag[tri]
        area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            pair = (int(min(a, b)), int(max(a, b)))
            edge_area[pair] = edge_area.get(pair, 0.0) + area / 3.0
            edge_ntri[pair] = edge_ntri.get(pair, 0) + 1

    connections = [(i, j, 'boundary' if edge_ntri[(i, j)] == 1 else 'fill')
                   for (i, j) in sorted(edge_area)]
    edge_k = {}
    for (i, j), A in edge_area.items():
        L = np.hypot(xLag[i] - xLag[j], yLag[i] - yLag[j])
        edge_k[(i, j)] = to_ib(stiffness_scale * E_material * A / L ** 2, ds)
    return connections, edge_k


def _tail_outline(wall_thickness):
    '''
    The tail as one polygon with holes: the outer outline and the 5 open
    windows (between base, rays and tip), each an (n, 2) vertex array (first
    vertex not repeated). Every wall is wall_thickness wide, centered on the
    triangle's edges and on each ray, with sharp (mitered) corners.
    The structure is currently hardcoded with 4 evenly spaced inner rays.
        TODO: implement variable structure parameters (e.g. # rays, length, base width)

    NOTE: the tail is placed with its base vertical and the tip pointing right in the domain.

    Returns (outer, windows): outer is the outline polygon, windows is a list of
    the 5 open windows (holes) inside it.
    '''
    L = 0.136   # tail length
    W = 0.051   # base width
    rayPositions = [0.025, 0.047, 0.069, 0.091]  # distance of each ray from the base
    x0, y0 = 0.3, 0.5   # position of the base's midpoint in the domain
    halfT = wall_thickness / 2.0

    s = np.hypot(1.0, W / (2 * L))   # 1/cos of the slanted edges' slope

    def half_width(x, off):   # half-width at x between the slanted edges, each moved outward by off
        return W / 2 * (1 - (x - x0) / L) + off * s

    def x_tip(off):           # where half_width(x, off) = 0
        return x0 + L * (1 + 2 * off * s / W)

    def quad(xa, xb, off):    # the region xa <= x <= xb between the moved slanted edges
        pts = [(xa, y0 - half_width(xa, off)), (xb, y0 - half_width(xb, off)),
               (xb, y0 + half_width(xb, off)), (xa, y0 + half_width(xa, off))]
        if half_width(xb, off) < 1e-12:   # xb is the tip: both right corners coincide
            del pts[2]
        return np.array(pts)

    outer = quad(x0 - halfT, x_tip(halfT), halfT)
    faces = ([x0 + halfT]                                            # inner face of the base
             + [x0 + d + side * halfT for d in rayPositions for side in (-1, 1)]
             + [x_tip(-halfT)])                                      # inner tip
    windows = [quad(xa, xb, -halfT) for xa, xb in zip(faces[::2], faces[1::2])]
    return outer, windows


def _Build_Tail_Mesh(ds, wall_thickness, fill_spacing, fill_density_multiplier):
    '''
    Meshes the whole tail at once: points every ~h along the outline and
    every window, Halton fill inside, one constrained triangulation over all
    of it. h = fill_spacing / sqrt(multiplier) is the fill's mean spacing, so
    edge and interior points are equally dense (h = fill_spacing for
    multiplier <= 1).

    Returns:
        xLag, yLag: point coordinates; outline points first (indices
            0..n_boundary-1), then the fill points
        triangles: (nt, 3) array of point indices
        n_boundary: number of outline points
        baseCornerIdx: (outer_bottom, outer_top), the base's two outer corner
            points (seeds for the actuation patches)
        outer, windows: the polygons from _tail_outline
    '''
    outer, windows = _tail_outline(wall_thickness)
    h = fill_spacing / np.sqrt(max(fill_density_multiplier, 1.0))
    loops = [mu.resample_loop(p, h) for p in [outer] + windows]
    fill = mu.fill_polygon_with_halton(outer, fill_spacing, inner_loops=windows,
                                       density_multiplier=fill_density_multiplier)
    pts, _, triangles = mu.triangulate_with_holes(loops[0], loops[1:], np.vstack(loops), fill)

    corners = (outer[0], outer[-1])   # outer outline's base corners (bottom, top)
    baseCornerIdx = tuple(int(np.argmin(np.linalg.norm(pts - c, axis=1))) for c in corners)
    return pts[:, 0], pts[:, 1], triangles, sum(len(l) for l in loops), baseCornerIdx, outer, windows


if __name__ == "__main__":
    FinRay_Geom()
