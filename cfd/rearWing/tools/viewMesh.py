# viewMesh.py -- automatic mesh images for the Ahmed body case
#
# Usage, from the case directory (Allrun.pre runs it after snappyHexMesh):
#     pvbatch viewMesh.py
# On a node without a display:
#     xvfb-run -a pvbatch viewMesh.py
#
# Output, in report/mesh/fig/ (the figures of the mesh report, scripts/report):
#     mesh_symmetry.png   mid-plane mesh: body, refinement boxes, near wake
#     mesh_rear.png       close-up of the trailing edges and near wake (layers)
#     mesh_surface.png    surface mesh on the walls
#     mesh_layers.png     number of layers extruded on each face (snappy)
#
# Cameras are not tuned by hand: they are placed from the body sizes in
# system/meshVariables (scripts/geometryInfo), so they follow any STL.
# The mid-plane section is a "crinkle" slice at y = cutY (scripts/postSetup):
# it shows whole cells, not a cut through them.

import os

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

# ------------------------------------------------------------ data
caseDir = os.getcwd()
# <case>.foam: empty file that the OpenFOAM reader of ParaView opens (as
# paraFoam does); the same file can be opened in the ParaView window
foamFile = os.path.join(caseDir, os.path.basename(caseDir) + '.foam')
outDir = os.path.join(caseDir, 'report', 'mesh', 'fig')
os.makedirs(outDir, exist_ok=True)
if not os.path.exists(foamFile):                # empty marker file
    open(foamFile, 'w').close()


def variables(path):
    """name value; lines of an OpenFOAM dictionary -> dict: numbers as
    floats, words as strings, lists "(a b)" as lists of strings"""
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.split('//')[0].strip()
        if not line.endswith(';') or len(line.split()) < 2:
            continue
        name, value = line[:-1].split(None, 1)
        value = value.strip()
        if value.startswith('(') and value.endswith(')'):
            out[name] = value[1:-1].split()
            continue
        try:
            out[name] = float(value)
        except ValueError:
            if len(value.split()) == 1:
                out[name] = value.strip('"')
    return out


var = variables(os.path.join(caseDir, 'system', 'meshVariables'))
var.update(variables(os.path.join(caseDir, 'system', 'meshParameters')))
x0, x1 = var['stlXmin'], var['stlXmax']         # nose, base
y0 = var.get('yCentre', var['stlYmin'])         # symmetry / centre plane
HALF = var.get('halfModel', 'yes') == 'yes'     # half model: mirror images
WALLS = var.get('wallPatches', ['body', 'legs'])  # the STL regions
L, H = var['bodyLength'], var['bodyHeight']
zTop, zGround = var['bodyTop'], var['ground']
# lowest z of the views: the ground, or half a height H below the STL when
# the ground is far (groundPlane none, isolated wing)
zGround = max(zGround, var['stlZmin'] - 0.5 * H)
# wing elements (refRegion): mid-height and trailing edge, the centre of the
# section views (the Ahmed-body views started from the ground, which for the
# isolated wing leaves the wing out of the image)
zMid = zTop - 0.5 * H
c = var.get('refLength', L)
xTE = var.get('xCofR', x0 + 0.25 * c) + 0.75 * c
caseType = ('Decomposed Case' if os.path.isdir('processor0')
            else 'Reconstructed Case')

# ------------------------------------------------------------ view
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1920, 1080]
view.OrientationAxesVisibility = 1
try:                                            # white background
    view.UseColorPaletteForBackground = 0       # (property names differ
    view.BackgroundColorMode = 'Single Color'   #  between versions)
except AttributeError:
    pass
view.Background = [1.0, 1.0, 1.0]


def side_camera(xa, xb, zc, resolution):
    """parallel view along +y of x in [xa, xb], centred on z = zc"""
    aspect = resolution[0] / resolution[1]
    half = 0.5 * (xb - xa) / aspect             # half of the visible height
    xc = 0.5 * (xa + xb)
    view.CameraParallelProjection = 1
    view.CameraPosition = [xc, y0 - 10 * L, zc]
    view.CameraFocalPoint = [xc, y0, zc]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelScale = half


def screenshot(name, resolution):
    SaveScreenshot(os.path.join(outDir, name), view,
                   ImageResolution=resolution)
    print('  written report/mesh/fig/' + name)


# ================================================== 1-2: mid-plane mesh
volume = OpenFOAMReader(registrationName='volume', FileName=foamFile)
volume.CaseType = caseType
volume.MeshRegions = ['internalMesh']
volume.CellArrays = []                          # only the cells, no fields
volume.Decomposepolyhedra = 0                   # keep polyhedra whole
volume.UpdatePipeline()

section = Slice(registrationName='section', Input=volume)
section.SliceType = 'Plane'
section.SliceType.Origin = [x0, var.get('cutY', y0 + 0.01 * L), 0.0]
section.SliceType.Normal = [0.0, 1.0, 0.0]
section.Crinkleslice = 1
d = Show(section, view)
d.SetRepresentationType('Surface With Edges')
d.ColorArrayName = ['POINTS', '']               # solid colour: only the cells
d.AmbientColor = d.DiffuseColor = [0.85, 0.85, 0.85]
d.EdgeColor = [0.1, 0.1, 0.3]

# 1: from the nose - H/2 to the end of the near box (base + 2H): body,
#    near box, wake box and the transitions between the levels
res = [2700, 900]
side_camera(x0 - 0.5 * H, x1 + 2.0 * H, zMid, res)
screenshot('mesh_symmetry.png', res)

# 2: trailing edges of mainplane and flap and the near wake (from half the
#    chord before the trailing edge to H behind it): the layers are visible
res = [1920, 1080]
side_camera(xTE - 0.5 * c, xTE + H, zMid, res)
screenshot('mesh_rear.png', res)
Hide(section, view)

# ===================================== 3-4: body surface and its layers
walls = OpenFOAMReader(registrationName='walls', FileName=foamFile)
walls.CaseType = caseType
walls.MeshRegions = ['patch/' + w for w in WALLS]
walls.Decomposepolyhedra = 0
walls.SkipZeroTime = 0                          # the field is in 0/: list
ReloadFiles(walls)                              # time 0 too (re-scan)
walls.CellArrays = ['nSurfaceLayers']           # written by snappyHexMesh
walls.UpdatePipeline(0.0)
view.ViewTime = 0.0
full = (Reflect(Input=walls, Plane='Y', Center=y0, CopyInput=1)  # whole body
        if HALF else walls)

# from behind, beside and above: upper surfaces, flap and endplates
view.CameraParallelProjection = 0
view.CameraPosition = [x1 + 1.2 * L, y0 - 1.4 * L, zTop + 0.9 * L]
view.CameraFocalPoint = [0.5 * (x0 + x1), y0, zMid]
view.CameraViewUp = [0.0, 0.0, 1.0]

d = Show(full, view)
d.SetRepresentationType('Surface With Edges')
d.ColorArrayName = ['POINTS', '']   # solid colour
screenshot('mesh_surface.png', [1920, 1080])

# number of layers per face: nLayers where the whole stack was extruded
if 'nSurfaceLayers' in walls.CellData.keys():
    d.SetRepresentationType('Surface')
    ColorBy(d, ('CELLS', 'nSurfaceLayers'))
    lut = GetColorTransferFunction('nSurfaceLayers')
    lut.ApplyPreset('Turbo', True)
    lut.RescaleTransferFunction(0.0, var.get('nLayers', 1.0))
    d.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = 'layer'
    bar.ComponentTitle = ''
    bar.TitleColor = bar.LabelColor = [0, 0, 0]
    screenshot('mesh_layers.png', [1920, 1080])
