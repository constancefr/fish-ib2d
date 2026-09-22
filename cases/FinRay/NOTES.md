In FinRay_Geom.py -----------

_actuation_patch_indices / actuation_patch_radius (base actuation as a
patch, not a single point):
    Each base corner used to have its actuation displacement applied to a
    SINGLE Lagrangian point (and, because of the update_Target_Point_
    Positions.py argmax/argmin row lookup, actually only ONE of that
    corner's 2 rail points -- the other sat completely fixed, never
    updated). Concentrating the whole prescribed motion onto one point
    while its immediate ~wall_thickness-away neighbor stayed rigid was
    producing a visible local twist right at the corner, distinct from the
    real Fin Ray shear/twist that alternately driving the two corners is
    supposed to produce (see update_Target_Point_Positions.py's own
    docstring -- the alternating-corner scheme itself is intentional and
    unchanged; this only addresses how each corner's own patch is driven).

    Fixed by _actuation_patch_indices: expands each corner into every
    Lagrangian point within actuation_patch_radius of that corner's true
    position, restricted to BOUNDARY (rail) points only (ribbonPolys'
    interior_idx points are excluded). update_Target_Point_Positions.py
    then classifies every target row into the "top" or "bottom" patch by
    baseline y vs. the midline, and moves every point in a patch by the
    identical rigid displacement -- approximating a real actuator's finite
    contact area. Verified directly: at a given phase, every row in a
    patch gets bit-identical dx (np.allclose), confirming rigid (not
    independent) motion.

    STABILITY WARNING (confirmed, not hypothetical): the natural first
    choice, actuation_patch_radius defaulting to ~ds, is WRONG for this
    geometry and was rejected. At ds-scale radius, the patch pulled in the
    nearest point on the ADJACENT top/bottom slant edge (~13-14mm away)
    rather than more of the base edge itself (whose own next point is
    ~17mm away, i.e. *further* than the adjacent edge's point) -- ds
    happens to be roughly comparable to both distances for this
    geometry's proportions, so nothing guarantees a ds-sized radius stays
    within the same edge. Forcing that adjacent-edge point into rigid
    lockstep with the corner (a location it doesn't physically belong to)
    measurably destabilized the simulation: reproduced going fully NaN
    (case: E_material=0.74e6, dt=1e-3, ribbon_fill_density_multiplier=32,
    ribbon_fill_stiffness_scale=1.0, FREQUENCY=0.1, Tfinal=5.0) with an
    8-point patch (ds-radius, boundary-only), while the exact same case
    stayed clean through the full run with a 4-point patch (both rails of
    the true corner only, no extra points).

    Fixed by defaulting actuation_patch_radius to 1.0*wall_thickness
    instead of ds -- ties the radius to the corner's own actual physical
    scale (the 2 rails are only ~wall_thickness apart) rather than the
    mesh resolution, which for this geometry keeps it safely under the
    ~13mm distance to the nearest unrelated point on an adjacent edge.
    Re-verified stable (same Tfinal=5.0 case) at this default.

    Real, NOT-yet-solved limitation: this is a Euclidean-distance patch,
    not an edge-aware (arc-length-along-the-same-rail-chain) one, so
    future geometry changes (different wall_thickness/rayPositions/L/W) could in
    principle put another edge's point closer than intended again --
    there's no general guarantee a given actuation_patch_radius stays
    within the base edge for arbitrary geometry. If you change
    actuation_patch_radius (or the base geometry) from these defaults,
    re-check stability the same way this was found: run to Tfinal and
    check run_case.py's own divergence warning (or grep the last
    viz_IB2d/lagPtsConnect.*.vtk frame for "nan").

_n_points_for_length (boundary density non-uniformity, FOUND + FIXED):
    The base, top/bottom slant edges, and each of the 4 rays used to each
    independently compute their own point count as
    n_points = max(round(length/ds), 2). This looks reasonable but is
    wrong for an OPEN chain: n points give n-1 SEGMENTS, and rounding the
    POINT count directly silently drops that -1. Negligible when there are
    many segments (the slant edges, 8-9 points, ended up ~110% of ds --
    fine), but severe when there are few (the base and rays, 2-3 points):
    measured actual segment lengths (at this example's default ds) ranged
    up to 214% of ds (ray at d=0.047 -- 2 points forced by rounding down
    from 2.14, giving exactly 1 segment spanning its FULL length), with
    the base itself at 163% of ds. Boundary density was NOT actually
    uniform across the structure, worst specifically at the base and the
    inner rays -- exactly the regions carrying the most actuation-driven
    deformation (see the triangle-inversion note below).

    Fixed by rounding to the nearest number of SEGMENTS first
    (n_segments = max(round(length/ds), 1)), then deriving point count
    from that (n_points = n_segments + 1) -- applied uniformly to the
    base, both slant edges, and every ray via the new
    _n_points_for_length() helper. Result: base segments now 94-124% of
    ds (was 163%), rays now 70-96% of ds (was up to 214%), slant edges
    ~98% of ds (was ~110%, also improved). Verified: no degenerate
    springs at any tested ribbon_fill_density_multiplier, and re-ran the
    triangle-inversion check below (same case: mult=32, scale=1.0,
    FREQUENCY=0.1, Tfinal=5.0) -- inversions are delayed and reduced but
    NOT eliminated (outline: first inversion t=0.44s->0.62s, peak
    5->4/378 triangles; ray0: t=0.9s->1.58s, peak 4->1/41 triangles).
    Consistent with the inversion note's own conclusion: this fix removes
    unnecessary sparseness (one contributing factor), but doesn't address
    the deeper limitation that a pure spring network has no mechanism to
    actually PREVENT inversion, only resist it proportionally to local
    stiffness.

'fill' triangle inversion under large deformation (KNOWN, not a bug):
    During simulation (not the static geometry_preview.png), some 'fill'
    triangles in the outline and ray meshes can invert (flip orientation)
    under large actuation-driven bending -- visible as an interior Halton
    point ending up geometrically on the wrong side of a nearby boundary
    edge. Confirmed by tracking each fill triangle's signed area over a
    full simulation run (case: E_material=0.74e6, dt=1e-3,
    ribbon_fill_density_multiplier=32, ribbon_fill_stiffness_scale=1.0,
    FREQUENCY=0.1, Tfinal=5.0) and comparing sign to its rest
    configuration. Before the _n_points_for_length fix above: outline's
    first inversion at t=0.44s (peak 5/372 triangles concurrently
    inverted), ray0's (the widest ray) at t=0.9s (peak 4/38). Both are
    small, localized, and self-correcting -- 'fill' springs are linear, so
    an inverted triangle isn't permanently stuck, it gets pulled back.

    Root cause: a pure mass-spring lattice has no explicit constraint
    against a triangle folding over -- springs only RESIST it in
    proportion to local stiffness and mesh density, they don't PREVENT it
    (unlike e.g. a co-rotational FEM formulation or an explicit
    area-preservation/inversion-barrier term, neither of which this model
    has). The base is hit first/worst because it's simultaneously (a) the
    most forced part of the structure -- actuation drives it directly via
    the target points -- and (b) was also the most sparsely-resolved
    boundary in the whole structure (see _n_points_for_length above),
    meaning larger, less-braced triangles there.

    _n_points_for_length measurably helps (delays onset, reduces peak
    concurrent inversions -- see that note for numbers) by removing the
    unnecessary sparseness, but does not eliminate inversion entirely,
    since the underlying "nothing prevents it" limitation is still there.
    If it becomes more than cosmetic, further options (not yet
    implemented): push ribbon_fill_density_multiplier/
    ribbon_fill_stiffness_scale higher specifically in the base region, or
    a more invasive fix like an explicit anti-inversion penalty term.

print_Lagrangian_Springs (length-aware stiffness):
    Each 'outline'/'ray'/'rung'/'diag' spring's stiffness is computed
    per-spring as k = EA_by_kind[kind] / ds_Rest (ds_Rest = that spring's
    OWN measured rest length), not a single k_Spring/k_Rung constant
    shared by every spring of that kind. This is the standard truss-element
    formula: for a rod of Young's modulus E, cross-section A, rest length
    L, under strain eps=dL/L, force = E*A*eps = (E*A/L)*dL, so
    k = E*A/L -- a SHORTER piece of the same material is stiffer (less
    material to stretch across), a LONGER piece is softer.

    Why this matters here specifically: FinRay_Geom() used to compute
    k_Spring = E*A/ds ONCE using the nominal ds, and hand that same number
    to every 'outline'/'ray' spring regardless of its own actual length.
    Most rail segments ARE close to ds, so this was a fine approximation
    for them -- but weld_to_inner_rail's junction splicing creates some
    springs as short as ~6% of ds. Under the old constant-k scheme those
    were silently ~17x too soft for their length (k should scale as 1/L,
    but they got the same k as a full-ds spring). Fixed by having
    FinRay_Geom() pass EA_by_kind (E*A or G*A products, NOT yet divided by
    a length -- recovered from the existing k_Spring/k_Rung via
    E*A = k_Spring*ds, so spring_Stiffness_From_Material/
    shear_Stiffness_From_Material didn't need to change) and computing
    k = EA_by_kind[kind] / ds_Rest inside print_Lagrangian_Springs, per
    spring. k_Spring/k_Rung themselves are still used as-is for k_Beam/
    k_Target, which want one representative magnitude, not a per-spring
    value. Verified: springs near 100% of ds keep ~the same k as before;
    a 6%-of-ds spring's k rose ~17x; a spring at ~190% of ds dropped to
    ~half -- all consistent with k=EA/L. Re-ran the fluid dynamics check
    (Tfinal=5.0, ribbon_fill_density_multiplier=1.0,
    ribbon_fill_stiffness_scale=0.1) after this change: still clean, no
    NaN, despite several springs getting substantially stiffer.

beam_Stiffness_From_Material:
    IB2d's torsional spring energy is defined via a signed-area
    form of "curvature" rather than a normalized curvature. Might have to
    verify that the scaling of k_Beam with E and I is correct for IB2d's
    discrete implementation.

    Recommended calibration before trusting comparisons between designs:
    build an isolated single strip (no ribs, no fluid), fix one end with
    stiff target points, apply a known force at the free tip, and check
    the resulting deflection against the analytical cantilever formula
        delta = F * L^3 / (3 * E * I)
    Adjust k_Beam until simulated and analytical deflection agree.

build_Tail_Beams:
    This does not add a beam spanning the
    ray-to-outline junction itself, so bending moment isn't transmitted
    across that connection (only the axial attachment spring force is
    via build_Tail_Connections). The junction currently behaves like a hinge
    instead of a fused connection.

    RESOLVED for the current ribbon model: _Build_Tail_Geometry_Ribbon's
    weld_to_inner_rail welds each ray rail's two ends directly into the
    outline's own inner-rail point sequence (a shared vertex, found via
    _mesh_utils.project_point_to_polyline -- the true point where the ray's
    edge meets the outline, not just the nearest existing outline point),
    instead of adding a separate 'attach' spring to a nearby-but-different
    point. The 'attach' kind and its k_Attach=50*k_Spring stiffness no
    longer exist. This still isn't a full moment-transmitting weld (no beam
    spans across the junction, so relative rotation between the ray and the
    outline right at that shared point is still only resisted by whatever
    rung/diag/fill springs happen to reach it, not by a beam), but the
    standoff gap + isolated one-spring joint is gone.

    BUG FOUND + FIXED during this change: weld_to_inner_rail originally
    only spliced new junction points into the outline's INNER rail
    (hole_pts/hole_idx), leaving the OUTER rail (ext_pts/ext_idx) at its
    original length. add_rungs_and_diagonals (in
    build_Tail_Ribbon_Connections) pairs the two rails by matching array
    position k, assuming len(railA) == len(railB) and that railA[k] sits
    directly across the ribbon from railB[k] -- true before any splicing
    (both rails come from the same centerline via _offset_polyline), but
    silently violated from the first inner-rail-only splice onward: railB
    now had extra points the array-index pairing didn't account for, so
    'rung'/'diag' springs past that point connected increasingly
    mismatched (physically far-apart) pairs -- visible as long spurious
    diagonal lines cutting straight across the tail's open "windows" in
    geometry_preview.png. Fixed by having weld_to_inner_rail splice a
    matching "partner" point into the outer rail too, at the same segment
    index and parametric position t (accurate given wall_thickness is tiny
    relative to the outline's feature sizes -- same justification used for
    the Delaunay centroid-filter approximation elsewhere in this file) --
    the partner point has no role of its own, it exists purely to keep
    both rails the same length with the same per-position correspondence.

build_Tail_Ribbon_Connections (interior 'fill' springs):
    Each ribbon's actual 2D area (not just its two rails) is now filled
    with Halton-sampled interior points (_mesh_utils.fill_polygon_with_halton)
    and Delaunay-triangulated (_mesh_utils.triangulate_with_holes) into
    'fill' springs, controlled by case_params.json's geom.ribbon_fill_*
    keys. Their stiffness is NOT the rail/rung/diag flange-web formula used
    elsewhere in this file -- it's a separate tributary-area lattice-spring
    approximation:

        k_edge = ribbon_fill_stiffness_scale * E_material * extrude_depth
                 * A_edge / L_edge**2

    where A_edge is 1/3 of the area of each triangle sharing that edge
    (equal-split lumping), derived by equating discrete spring energy to
    continuum strain energy over the edge's tributary volume
    (A_edge * extrude_depth).

    Known limitations, not yet resolved:
      - A uniform triangular central-force lattice has a fixed effective
        Poisson's ratio around 1/3 and cannot independently target the
        near-incompressible nu=0.499 assumed by shear_Stiffness_From_Material
        for the rung/diag springs -- the ribbon is now a hybrid of two
        different stiffness derivations at its rail boundary vs. its
        interior fill, not one consistent continuum model.
      - Only genuinely new edges (touching an interior point, or a
        boundary-boundary diagonal not already covered by the existing
        outline/ray/rung/diag logic) get this formula; existing rail
        edges are untouched.

    Recommended calibration before trusting comparisons across
    ribbon_fill_density_multiplier: same approach as the beam calibration
    above -- isolate a single ribbon (no fluid, no other rays), clamp one
    end via target points, apply a known tip force, and compare deflection
    to the analytical cantilever formula. Sweep ribbon_fill_stiffness_scale
    until simulated and analytical deflection agree. See
    verify_ribbon_fill.py for a fluid-free geometry/mesh sanity check
    (degenerate-spring / runaway-count check across several density
    multipliers) -- it does not do this force calibration.

    STABILITY WARNING (confirmed, not hypothetical): at this example's
    defaults (dt=1e-3, E_material=0.74e6, ribbon_fill_density_multiplier=1.0),
    ribbon_fill_stiffness_scale=1.0 (the formula's literal, uncalibrated
    value) makes the simulation go fully NaN partway through the run --
    reproduced by running to Tfinal=5.0 and inspecting the .vtk point data:

        scale=1.0, 0.5:  diverges to NaN around t=3.3-3.6s
        scale=0.1, 0.01: clean (no NaN) through the full Tfinal=5.0

    BASE_CASE in run_case.py now defaults to scale=0.1 because of this.
    The likely mechanism: some 'fill' edges reach k on the order of the
    model's previously-stiffest element (k_Attach), and adding dozens of
    comparably-stiff edges pushes the explicit time-stepping scheme (fixed
    dt, no dt/k coupling) past its stability margin -- consistent with a
    slowly-growing resonance rather than an immediate blow-up (it takes
    several oscillation periods of the FREQUENCY=1.0 Hz actuation to
    appear). This 0.1 default is an empirically-found *starting point*,
    not a substitute for the tip-deflection calibration above, and it will
    need re-checking (rerun a case to Tfinal and grep its last .vtk frame
    for "nan") whenever E_material, dt, or ribbon_fill_density_multiplier
    change materially from these defaults.

    UPDATE after the ray-to-outline weld fix (build_Tail_Beams note above):
    re-ran the same check (Tfinal=5.0, ribbon_fill_density_multiplier=1.0)
    with the welded-junction geometry and scale=1.0 no longer diverges --
    clean through the full run, where it previously went NaN by t~3.3-3.6s.
    Consistent with the likely mechanism above: k_Attach (50*k_Spring, the
    old model's stiffest element by far) is gone entirely now, so there's
    no longer an isolated very-stiff spring for the explicit integrator to
    struggle with at that joint. NOT yet re-tested at higher
    ribbon_fill_density_multiplier (8.0+) -- BASE_CASE's scale=0.1 default
    is left as-is (still valid, just conservative) rather than bumped back
    to 1.0, since multiplier=1.0 stability doesn't guarantee it at higher
    multipliers. Worth re-sweeping scale at higher multipliers now that the
    joint itself is no longer the likely bottleneck.