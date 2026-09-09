In FinRay_Geom.py -----------

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