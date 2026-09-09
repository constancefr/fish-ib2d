'''
VisIt batch movie script for the FinRay IB2d example.

Usage (from a normal shell, not VisIt's GUI):
    visit -cli -nowin -s make_movie.py -- <case_dir> <output.mp4> [fps]
EXAMPLE FOR BASE CASE amp0.010_E7.40e+05:
    /Applications/VisIt.app/Contents/Resources/bin/visit -quiet -cli -nowin -s make_movie.py -- results/amp0.010_E7.40e+05 output_sim.mp4

<case_dir> should contain a viz_IB2d/ folder, as produced by run_case.py
or a plain `python main2d.py` run. Renders one PNG per saved timestep,
then encodes them into an mp4 with ffmpeg.
'''

import glob
import os
import subprocess
import sys

# args after the "--" VisIt passes through to sys.argv
argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
case_dir = os.path.abspath(argv[0])
movie_out = os.path.abspath(argv[1])
fps = int(argv[2]) if len(argv) > 2 else 20

viz_dir = os.path.join(case_dir, "viz_IB2d")
frame_dir = os.path.join(case_dir, "_frames")
print(f"CASE_DIR: {case_dir}")
print(f"VIZ_DIR: {viz_dir}\nFRAME_DIR: {frame_dir}\nMOVIE_OUT: {movie_out}\nFPS: {fps}")
os.makedirs(frame_dir, exist_ok=True)


def db(pattern):
    """
    Build the 'family' database string VisIt uses to treat a numbered
    vtk sequence (lagsPts.0000.vtk, lagsPts.0001.vtk, ...) as one
    time-varying database.
    """
    matches = sorted(glob.glob(os.path.join(viz_dir, pattern)))
    if not matches:
        raise FileNotFoundError(f"No files matching {pattern} in {viz_dir}")

    visit_file = os.path.join(viz_dir, f"_{pattern.replace('*', 'series')}.visit")
    with open(visit_file, "w") as f:
        f.write("\n".join(matches) + "\n")
    return visit_file


def build_plots():
    DeleteAllPlots()

    # --- Raw Lagrangian points as small dots on top of the mesh
    OpenDatabase(db("lagsPts.*.vtk"))
    AddPlot("Mesh", "mesh")
    
    m1 = MeshAttributes()
         # argument 1 would indicate we want to return the current attributes, not the default ones
    m1.legendFlag = 1
    m1.lineWidth = 0
    m1.meshColor = (0, 0, 0, 255)
    m1.meshColorSource = m1.MeshCustom  # Foreground, MeshCustom, MeshRandom
    m1.opaqueColorSource = m1.Background  # Background, OpaqueCustom, OpaqueRandom
    m1.opaqueMode = m1.Auto  # Auto, On, Off
    m1.pointSize = 0.05
    m1.opaqueColor = (255, 255, 255, 255)
    m1.smoothingLevel = m1.NONE  # NONE, Fast, High
    m1.pointSizeVarEnabled = 0
    m1.pointSizeVar = "default"
    m1.pointType = m1.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    m1.showInternal = 0
    m1.showGenerated = 0
    m1.pointSizePixels = 5
    m1.opacity = 1
    SetActivePlots(0)
    SetPlotOptions(m1)

    # --- Structure: connected Lagrangian mesh (springs/beams/target pts)
    OpenDatabase(db("lagPtsConnect.*.vtk"))
    AddPlot("Mesh", "mesh")
    SetActivePlots(1) # just in case, since we now have two mesh plots open

    m2 = MeshAttributes()
    m2.legendFlag = 1
    m2.lineWidth = 9 # for some reason this has to be -1 what you would input in the GUI hmm...
    m2.meshColor = (128, 128, 128, 255)
    m2.meshColorSource = m2.MeshCustom  # Foreground, MeshCustom, MeshRandom
    m2.opaqueColorSource = m2.Background  # Background, OpaqueCustom, OpaqueRandom
    m2.opaqueMode = m2.Auto  # Auto, On, Off
    m2.pointSize = 0.05
    m2.opaqueColor = (255, 255, 255, 255)
    m2.smoothingLevel = m2.NONE  # NONE, Fast, High
    m2.pointSizeVarEnabled = 0
    m2.pointSizeVar = "default"
    m2.pointType = m2.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    m2.showInternal = 0
    m2.showGenerated = 0
    m2.pointSizePixels = 2
    m2.opacity = 1
    SetPlotOptions(m2)


    # --- Velocity magnitude colormap
    OpenDatabase(db("uMag.*.vtk"))
    AddPlot("Pseudocolor", "uMag")

    p = PseudocolorAttributes()
    p.scaling = p.Linear  # Linear, Log, Skew
    p.skewFactor = 1
    p.limitsMode = p.OriginalData  # OriginalData, ActualData
    p.minFlag = 0
    p.min = 0
    p.useBelowMinColor = 0
    p.belowMinColor = (0, 0, 0, 255)
    p.maxFlag = 0
    p.max = 1
    p.useAboveMaxColor = 0
    p.aboveMaxColor = (0, 0, 0, 255)
    p.centering = p.Natural  # Natural, Nodal, Zonal
    p.colorTableName = "Blues"
    p.invertColorTable = 0
    p.opacityType = p.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    p.opacityVariable = ""
    p.opacity = 1
    p.opacityVarMin = 0
    p.opacityVarMax = 1
    p.opacityVarMinFlag = 0
    p.opacityVarMaxFlag = 0
    p.pointSize = 0.05
    p.pointType = p.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    p.pointSizeVarEnabled = 0
    p.pointSizeVar = "default"
    p.pointSizePixels = 2
    p.lineType = p.Line  # Line, Tube, Ribbon
    p.lineWidth = 0
    p.tubeResolution = 10
    p.tubeRadiusSizeType = p.FractionOfBBox  # Absolute, FractionOfBBox
    p.tubeRadiusAbsolute = 0.125
    p.tubeRadiusBBox = 0.005
    p.tubeRadiusVarEnabled = 0
    p.tubeRadiusVar = ""
    p.tubeRadiusVarRatio = 10
    p.tailStyle = p.NONE  # NONE, Spheres, Cones
    p.headStyle = p.NONE  # NONE, Spheres, Cones
    p.endPointRadiusSizeType = p.FractionOfBBox  # Absolute, FractionOfBBox
    p.endPointRadiusAbsolute = 0.125
    p.endPointRadiusBBox = 0.05
    p.endPointResolution = 10
    p.endPointRatio = 5
    p.endPointRadiusVarEnabled = 0
    p.endPointRadiusVar = ""
    p.endPointRadiusVarRatio = 10
    p.renderSurfaces = 1
    p.renderWireframe = 0
    p.renderPoints = 0
    p.smoothingLevel = 0
    p.legendFlag = 1
    p.lightingFlag = 1
    p.wireframeColor = (0, 0, 0, 0)
    p.wireframeColorByVar = 0
    p.pointColor = (0, 0, 0, 0)
    p.pointColorByVar = 0
    SetActivePlots(2)
    SetPlotOptions(p)

    # --- Velocity vectors
    OpenDatabase(db("u.*.vtk"))
    AddPlot("Vector", "u")

    v = VectorAttributes()
    v.glyphLocation = v.UniformInSpace  # AdaptsToMeshResolution, UniformInSpace
    v.useStride = 0
    v.nVectors = 10000
    v.stride = 1
    v.origOnly = 1
    v.limitsMode = v.OriginalData  # OriginalData, CurrentPlot
    v.minFlag = 0
    v.min = 0
    v.maxFlag = 0
    v.max = 1
    v.colorByMagnitude = 1
    v.colorTableName = "difference"
    v.invertColorTable = 0
    v.vectorColor = (0, 0, 0, 255)
    v.useLegend = 1
    v.scale = 0.15
    v.scaleByMagnitude = 1
    v.autoScale = 1
    v.glyphType = v.Arrow  # Arrow, Ellipsoid
    v.headOn = 1
    v.headSize = 0.15
    v.lineStem = v.Line  # Cylinder, Line
    v.lineWidth = 0
    v.stemWidth = 0.08
    v.vectorOrigin = v.Tail  # Head, Middle, Tail
    v.geometryQuality = v.Fast  # Fast, High
    v.animationStep = 0
    SetActivePlots(3)
    SetPlotOptions(v)

    DrawPlots()

    # --- View + annotations ---]
    view = View2DAttributes()
    view.windowCoords = (0.208638, 0.675145, 0.25991, 0.726417)
    view.viewportCoords = (0.1, 0.9, 0.1, 0.9)
    view.fullFrameActivationMode = view.Auto  # On, Off, Auto
    view.fullFrameAutoThreshold = 100
    view.xScale = view.LINEAR  # LINEAR, LOG
    view.yScale = view.LINEAR  # LINEAR, LOG
    view.windowValid = 1
    SetView2D(view)
    
    annot = AnnotationAttributes()
    annot.axes2D.visible = 0
    annot.axes2D.autoSetTicks = 1
    annot.axes2D.autoSetScaling = 1
    annot.axes2D.lineWidth = 0
    annot.axes2D.tickLocation = annot.axes2D.Outside  # Inside, Outside, Both
    annot.axes2D.tickAxes = annot.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
    annot.axes2D.xAxis.title.visible = 1
    annot.axes2D.xAxis.title.font.font = annot.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
    annot.axes2D.xAxis.title.font.scale = 1
    annot.axes2D.xAxis.title.font.useForegroundColor = 1
    annot.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
    annot.axes2D.xAxis.title.font.bold = 1
    annot.axes2D.xAxis.title.font.italic = 1
    annot.axes2D.xAxis.title.userTitle = 0
    annot.axes2D.xAxis.title.userUnits = 0
    annot.axes2D.xAxis.title.title = "X-Axis"
    annot.axes2D.xAxis.title.units = ""
    annot.axes2D.xAxis.label.visible = 1
    annot.axes2D.xAxis.label.font.font = annot.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
    annot.axes2D.xAxis.label.font.scale = 1
    annot.axes2D.xAxis.label.font.useForegroundColor = 1
    annot.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
    annot.axes2D.xAxis.label.font.bold = 1
    annot.axes2D.xAxis.label.font.italic = 1
    annot.axes2D.xAxis.label.scaling = 0
    annot.axes2D.xAxis.tickMarks.visible = 1
    annot.axes2D.xAxis.tickMarks.majorMinimum = 0
    annot.axes2D.xAxis.tickMarks.majorMaximum = 1
    annot.axes2D.xAxis.tickMarks.minorSpacing = 0.02
    annot.axes2D.xAxis.tickMarks.majorSpacing = 0.2
    annot.axes2D.xAxis.grid = 0
    annot.axes2D.yAxis.title.visible = 1
    annot.axes2D.yAxis.title.font.font = annot.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
    annot.axes2D.yAxis.title.font.scale = 1
    annot.axes2D.yAxis.title.font.useForegroundColor = 1
    annot.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
    annot.axes2D.yAxis.title.font.bold = 1
    annot.axes2D.yAxis.title.font.italic = 1
    annot.axes2D.yAxis.title.userTitle = 0
    annot.axes2D.yAxis.title.userUnits = 0
    annot.axes2D.yAxis.title.title = "Y-Axis"
    annot.axes2D.yAxis.title.units = ""
    annot.axes2D.yAxis.label.visible = 1
    annot.axes2D.yAxis.label.font.font = annot.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
    annot.axes2D.yAxis.label.font.scale = 1
    annot.axes2D.yAxis.label.font.useForegroundColor = 1
    annot.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
    annot.axes2D.yAxis.label.font.bold = 1
    annot.axes2D.yAxis.label.font.italic = 1
    annot.axes2D.yAxis.label.scaling = 0
    annot.axes2D.yAxis.tickMarks.visible = 1
    annot.axes2D.yAxis.tickMarks.majorMinimum = 0
    annot.axes2D.yAxis.tickMarks.majorMaximum = 1
    annot.axes2D.yAxis.tickMarks.minorSpacing = 0.02
    annot.axes2D.yAxis.tickMarks.majorSpacing = 0.2
    annot.axes2D.yAxis.grid = 0

    annot.databaseInfoFlag = 0
    annot.timeInfoFlag = 1
    annot.databaseInfoFont.font = annot.databaseInfoFont.Arial  # Arial, Courier, Times
    annot.databaseInfoFont.scale = 1
    annot.databaseInfoFont.useForegroundColor = 1
    annot.databaseInfoFont.color = (0, 0, 0, 255)
    annot.databaseInfoFont.bold = 0
    annot.databaseInfoFont.italic = 0
    annot.databaseInfoExpansionMode = annot.File  # File, Directory, Full, Smart, SmartDirectory
    annot.databaseInfoTimeScale = 1
    annot.databaseInfoTimeOffset = 0
    annot.legendInfoFlag = 0
    annot.backgroundColor = (255, 255, 255, 255)
    annot.foregroundColor = (0, 0, 0, 255)
    annot.gradientBackgroundStyle = annot.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
    annot.gradientColor1 = (0, 0, 255, 255)
    annot.gradientColor2 = (0, 0, 0, 255)
    annot.backgroundMode = annot.Solid  # Solid, Gradient, Image, ImageSphere
    annot.backgroundImage = ""
    annot.imageRepeatX = 1
    annot.imageRepeatY = 1
    annot.axesArray.visible = 1
    annot.axesArray.ticksVisible = 1
    annot.axesArray.autoSetTicks = 1
    annot.axesArray.autoSetScaling = 1
    annot.axesArray.lineWidth = 0
    annot.axesArray.axes.title.visible = 1
    annot.axesArray.axes.title.font.font = annot.axesArray.axes.title.font.Arial  # Arial, Courier, Times
    annot.axesArray.axes.title.font.scale = 1
    annot.axesArray.axes.title.font.useForegroundColor = 1
    annot.axesArray.axes.title.font.color = (0, 0, 0, 255)
    annot.axesArray.axes.title.font.bold = 0
    annot.axesArray.axes.title.font.italic = 0
    annot.axesArray.axes.title.userTitle = 0
    annot.axesArray.axes.title.userUnits = 0
    annot.axesArray.axes.title.title = ""
    annot.axesArray.axes.title.units = ""
    annot.axesArray.axes.label.visible = 1
    annot.axesArray.axes.label.font.font = annot.axesArray.axes.label.font.Arial  # Arial, Courier, Times
    annot.axesArray.axes.label.font.scale = 1
    annot.axesArray.axes.label.font.useForegroundColor = 1
    annot.axesArray.axes.label.font.color = (0, 0, 0, 255)
    annot.axesArray.axes.label.font.bold = 0
    annot.axesArray.axes.label.font.italic = 0
    annot.axesArray.axes.label.scaling = 0
    annot.axesArray.axes.tickMarks.visible = 1
    annot.axesArray.axes.tickMarks.majorMinimum = 0
    annot.axesArray.axes.tickMarks.majorMaximum = 1
    annot.axesArray.axes.tickMarks.minorSpacing = 0.02
    annot.axesArray.axes.tickMarks.majorSpacing = 0.2
    annot.axesArray.axes.grid = 0

    annot.userInfoFlag = 0
    SetAnnotationAttributes(annot)


def save_frames():
    swatts = SaveWindowAttributes()
    swatts.family = 0
    swatts.format = swatts.PNG
    swatts.width = 1280
    swatts.height = 960
    swatts.resConstraint = swatts.NoConstraint

    n_states = TimeSliderGetNStates()
    for state in range(n_states):
        TimeSliderSetState(state)
        swatts.fileName = os.path.join(frame_dir, f"frame_{state:04d}.png")
        SetSaveWindowAttributes(swatts)
        SaveWindow()
    return n_states


def encode(n_states):
    # ffmpeg gives more predictable control than visit_utils' own encoder
    pattern = os.path.join(frame_dir, "frame_%04d.png")
    subprocess.run([
        "ffmpeg", "-y", "-framerate", str(fps),
        "-i", pattern,
        "-pix_fmt", "yuv420p", movie_out,
    ], check=True)
    print(f"Wrote {movie_out} ({n_states} frames @ {fps} fps)")


build_plots()
_n_states = save_frames()
encode(_n_states)
sys.exit()