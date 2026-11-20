# images.py -- result images for the report (run by scripts/report)
#
# Usage, from the case directory:
#     pvbatch tools/report/images.py <time> <outDir>
#
# The views were designed in the ParaView window and recorded with
# Tools > Start Trace (trace_cp.py, trace_yplus.py, trace_midplane.py,
# trace_cp0.py, trace_streamlines.py). Here the traces are cleaned and made
# general: every position, size and camera distance recorded by hand is
# written in terms of the case quantities (system/meshVariables and
# meshParameters), so the images follow any STL and flow speed.
#   - data: ParaView never reads the volume mesh. Sections, planes,
#     iso-surfaces and streamlines are cut by OpenFOAM, in parallel, from the
#     last time (VTK files in postProcessing/<function>/<time>/):
#       wingSurfaces      walls, midPlane, wakePlane  (during the run,
#                                                      system/functions)
#       planes, isoSurfaces, noseStreamlines, wakeStreamlines
#                         (after the run, system/reportFunctions, started by
#                          scripts/report)
#     so the time to make the images does not depend on the mesh size;
#   - colour scales: physical, as chosen in the traces: Cp from -3 to 1
#     on the wing (1 = stagnation; -3: the suction peak of a high-lift wing,
#     -1 for the Ahmed body), Cp0 from 0 to 1 (1 = free stream), |U| from 0
#     to 1.5 Uinf, y+ in log scale one decade around the design value;
#   - cameras: view directions as recorded; distances, frames and focal
#     points scaled with the body (L = length, H = height).
#
# Output (PNG): cp_side, cp_top, cp_bottom, cp_rear, yplus_body, yplus_1_3,
#               yplus_45_55, wss_side, wss_top, u_midplane, cp_midplane,
#               cp0_wake, streamlines, streamlines_top, streamlines_bottom,
#               streamlines_side, q_iso, cp0_iso,
#               plane_<x|y|z>_<U|k|p|Cp0>_<01..N> (nPlanes per direction,
#               meshParameters; skipped with REPORT_PLANES=0)
#        (TeX): planes.tex, the tables of the planes for the report
#        (CSV): cp_centreline (x, z, Cp on the wing section y = cutY),
#               from the surfaces written by the function "wingSurfaces"
#  A missing file skips its images (with a message), not the whole report.

import glob
import math
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

time, outDir = sys.argv[1], sys.argv[2]
os.makedirs(outDir, exist_ok=True)
RES = [1920, 1080]                       # image size, as in the traces
ASPECT = RES[0] / RES[1]


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


var = variables('system/meshVariables')
var.update(variables('system/meshParameters'))
var.update(variables('constant/caseSettings'))       # QStar
x0, x1 = var['stlXmin'], var['stlXmax']          # nose, base
y0 = var.get('yCentre', var['stlYmin'])          # symmetry / centre plane
HALF = var.get('halfModel', 'yes') == 'yes'      # half model: mirror images
WALLS = ['patch/' + w for w in var.get('wallPatches', ['body', 'legs'])]
L, H, W = var['bodyLength'], var['bodyHeight'], var['bodyHalfWidth']
zTop, ground = var['bodyTop'], var['ground']
# lowest z of the views: the ground, or half a height H below the STL when
# the ground is far (groundPlane none, isolated wing), as in tools/viewMesh.py
ground = max(ground, var['stlZmin'] - 0.5 * H)
zMid = zTop - 0.5 * H                            # mid-height of the wing elements
zBottom = zTop - H                               # underbody
uInf = var['uInf']
# centre of the body (legs included): focal point of the body views
centre = [0.5 * (x0 + x1), y0, 0.5 * (zTop + ground)]

# ----------------------------------------------------------------- helpers
view = CreateView('RenderView')
view.ViewSize = RES
view.OrientationAxesVisibility = 0
try:
    view.UseColorPaletteForBackground = 0
    view.BackgroundColorMode = 'Single Color'
except AttributeError:
    pass
view.Background = [1.0, 1.0, 1.0]



def isTime(name):
    try:
        float(name)
        return True
    except ValueError:
        return False


def vtk(function, name):
    """VTK file postProcessing/<function>/<time>/<name>.vtk (time: the one
    given, or the last one of that function); None if missing"""
    folder = os.path.join('postProcessing', function)
    # ESI: the streamLine function writes in postProcessing/sets/<function>/
    if not os.path.isdir(folder):
        folder = os.path.join('postProcessing', 'sets', function)
    times = [time] if os.path.isdir(os.path.join(folder, time)) else \
        sorted((t for t in os.listdir(folder) if isTime(t)), key=float) \
        if os.path.isdir(folder) else []
    path = os.path.join(folder, times[-1], name + '.vtk') if times else ''
    # ESI: the streamLine function names its files itself (track0.vtp):
    # for 'tracks' the first VTK file of the folder is used
    if times and not os.path.exists(path) and name == 'tracks':
        found = sorted(glob.glob(os.path.join(folder, times[-1], '*.vt[kp]')))
        path = found[0] if found else path
    if not os.path.exists(path):
        print('  missing %s: skipped' % (path or folder))
        return None
    r = OpenDataFile(path)                       # legacy .vtk or XML .vtp
    r.UpdatePipeline()
    return r


def mirror(source):
    """whole body / domain: reflection on the symmetry plane y = y0 (half
    model; a whole body is already whole)"""
    if not HALF:
        return source
    m = Reflect(Input=source, Plane='Y', Center=y0, CopyInput=1)
    m.UpdatePipeline()
    return m


def mirrorLines(source):
    """as mirror, for lines (streamlines): the Reflect filter scrambles the
    points of long polylines, so the copy is made by a transform y -> 2 y0 - y
    and appended to the original"""
    if not HALF:
        return source
    t = Transform(Input=source)
    t.Transform = 'Transform'
    t.Transform.Translate = [0.0, 2 * y0, 0.0]
    t.Transform.Scale = [1.0, -1.0, 1.0]
    m = AppendDatasets(Input=[source, t])
    m.UpdatePipeline()
    return m


def body(source, grey):
    """the walls (surface of the function "wingSurfaces") in a solid grey"""
    d = Show(source, view)
    d.SetRepresentationType('Surface')
    d.ColorArrayName = ['POINTS', '']           # no field: solid colour
    d.AmbientColor = d.DiffuseColor = [grey] * 3
    return source


def hasPreset(name):
    """colour map available in this ParaView ('Fast', the default map of
    ParaView 6, is missing in older versions: 'Cool to Warm' instead)"""
    try:
        from paraview.modules.vtkRemotingViews import \
            vtkSMTransferFunctionPresets
        return vtkSMTransferFunctionPresets.GetInstance().HasPreset(name)
    except ImportError:
        return True


def colour(source, assoc, field, preset, vmin, vmax, black=True, log=False,
           component=None, below=None, above=None):
    """show source coloured by field, with a fixed scale and its legend
    (only this one: the legends of the previous fields are hidden)"""
    for f in ('Cp', 'Cp0', 'U', 'yPlus', 'Cf', 'Umag', 'k', 'p'):
        GetScalarBar(GetColorTransferFunction(f), view).Visibility = 0
    d = Show(source, view)
    d.SetRepresentationType('Surface')
    ColorBy(d, (assoc, field) if component is None
            else (assoc, field, component))
    lut = GetColorTransferFunction(field)
    if preset:
        lut.ApplyPreset(preset if hasPreset(preset) else 'Cool to Warm', True)
    lut.AutomaticRescaleRangeMode = 'Never'
    lut.RescaleTransferFunction(vmin, vmax)
    GetOpacityTransferFunction(field).RescaleTransferFunction(vmin, vmax)
    if log:
        lut.MapControlPointsToLogSpace()
    lut.UseLogScale = 1 if log else 0
    # values out of the scale in a grey (below, above), not in the end colours
    lut.UseBelowRangeColor = 1 if below else 0
    lut.UseAboveRangeColor = 1 if above else 0
    if below:
        lut.BelowRangeColor = below
    if above:
        lut.AboveRangeColor = above
    d.SetScalarBarVisibility(view, True)
    if black:
        bar = GetScalarBar(lut, view)
        bar.TitleColor = bar.LabelColor = [0.0, 0.0, 0.0]
    return d


def parallel(focal, direction, up, half):
    """parallel projection: view along direction, half the image height"""
    view.CameraParallelProjection = 1
    view.CameraFocalPoint = focal
    view.CameraPosition = [f - 10 * L * c for f, c in zip(focal, direction)]
    view.CameraViewUp = up
    view.CameraParallelScale = half


def perspective(focal, direction, distance, up, angle=30.0):
    """perspective: camera at distance from focal, looking along direction"""
    n = math.sqrt(sum(c * c for c in direction))
    view.CameraParallelProjection = 0
    view.CameraFocalPoint = focal
    view.CameraPosition = [f - distance * c / n for f, c in zip(focal, direction)]
    view.CameraViewUp = up
    view.CameraViewAngle = angle


def save(name):
    SaveScreenshot(os.path.join(outDir, name + '.png'), view,
                   ImageResolution=RES)
    print('  written ' + name + '.png')


def clear():
    for s in GetSources().values():
        d = GetDisplayProperties(s, view)
        if d.ColorArrayName[1]:                  # coloured by a field
            d.SetScalarBarVisibility(view, False)
        Hide(s, view)


# the walls (function "wingSurfaces": face values of Cp, y+, tau_w), the whole
# body for the images, the half model for the sections seen from the side
walls = vtk('wingSurfaces', 'walls')
if walls is None:
    sys.exit('no walls surface (function "wingSurfaces"): no images')
wallsFull = mirror(walls)

# ======================================== 1-3: Cp on the wing (trace_cp.py)
CPMIN = -3.0                                     # suction peak of the wing
colour(wallsFull, 'CELLS', 'Cp', 'Turbo', CPMIN, 1.0)
# side (+Y) and top (-Z) views, parallel: the body fills 3/4 of the image
# width from the side and 1/2 from above (as framed in the trace)
parallel(centre, [0, 1, 0], [0, 0, 1], L / (2 * 0.75 * ASPECT))
save('cp_side')
parallel(centre, [0, 0, -1], [0, 1, 0], max(L / (2 * 0.5 * ASPECT), 1.1 * W))
save('cp_top')
# from below: the suction side of the inverted wing
parallel(centre, [0, 0, 1], [0, 1, 0], max(L / (2 * 0.5 * ASPECT), 1.1 * W))
save('cp_bottom')
# from behind, above and beside: slant and base (direction from the trace,
# distance about twice the body length)
perspective(centre, [-1.597, -1.399, -0.747], 2.15 * L, [-0.272, -0.193, 0.943])
save('cp_rear')
clear()

# --------------- Cp along the mid-plane of the body (data for a plot)
# the body surface cut by the mid-plane section: Cp of the faces along the
# centre line of nose, roof, slant, base and underbody
cut = Slice(Input=walls, SliceType='Plane')
cut.SliceType.Origin = [centre[0], var.get('cutY', y0), centre[2]]
cut.SliceType.Normal = [0.0, 1.0, 0.0]
SaveData(os.path.join(outDir, 'cp_centreline.csv'),
         proxy=CellCenters(Input=cut), PointDataArrays=['Cp'])
print('  written cp_centreline.csv')

# ===================================== 4: y+ on the body (trace_yplus.py)
yp = var.get('yPlus', 1.0)
# log scale, one decade below and above the design y+
colour(wallsFull, 'CELLS', 'yPlus', 'Viridis (matplotlib)', yp / 10, yp * 10,
       log=True)
# from below, in front and beside: nose, side, underbody and legs
perspective(centre, [0.459, 0.839, 0.295], 1.6 * L,
            [-0.085, -0.289, 0.954])
save('yplus_body')
# the two ranges of the post-processing guidelines: 1-3 (first cell well
# inside the viscous sublayer) and 4.5-5.5 (edge of the viscous sublayer);
# out of range in dark grey (below) and light grey (above)
DARK, LIGHT = [0.3, 0.3, 0.3], [0.85, 0.85, 0.85]
colour(wallsFull, 'CELLS', 'yPlus', 'Viridis (matplotlib)', 1.0, 3.0,
       below=DARK, above=LIGHT)
save('yplus_1_3')
colour(wallsFull, 'CELLS', 'yPlus', 'Viridis (matplotlib)', 4.5, 5.5,
       below=DARK, above=LIGHT)
save('yplus_45_55')
clear()

# ========================= 4b: wall shear stress (function wallShearStress)
# skin-friction coefficient Cf = |tau_w| / (0.5 Uinf^2) (kinematic tau_w),
# from 0 to twice the value assumed for the near-wall mesh (meshParameters)
cf = Calculator(Input=wallsFull, AttributeType='Cell Data',
                ResultArrayName='Cf',
                Function='mag(wallShearStress)/%g' % (0.5 * uInf ** 2))
colour(cf, 'CELLS', 'Cf', 'Viridis (matplotlib)', 0.0, 2 * var.get('cf', 0.003),
       above=[0.97, 0.97, 0.97])
parallel(centre, [0, 1, 0], [0, 0, 1], L / (2 * 0.75 * ASPECT))
save('wss_side')
parallel(centre, [0, 0, -1], [0, 1, 0], L / (2 * 0.5 * ASPECT))
save('wss_top')
clear()

# ============== 5-6: mid-plane section, |U| and Cp (trace_midplane.py)
# section y = cutY (function "wingSurfaces"), seen from the side with the body
section = vtk('wingSurfaces', 'midPlane')
if section is not None:
    body(walls, 0.39)
    colour(section, 'POINTS', 'U', 'Turbo', 0.0, 1.5 * uInf, black=False,
           component='Magnitude')
    # side view from the ground up: from 0.75 H before the nose to 4 H
    # behind the base (near wake), as framed in the trace
    xa, xb = x0 - 0.75 * H, x1 + 4.0 * H
    half_mid = (xb - xa) / (2 * ASPECT)
    parallel([0.5 * (xa + xb), y0, zMid], [0, 1, 0], [0, 0, 1], half_mid)
    save('u_midplane')
    colour(section, 'POINTS', 'Cp', 'Fast', CPMIN, 1.0, black=False)
    save('cp_midplane')
    clear()

# ================== 7: total pressure in the wake (trace_cp0.py)
# section x = cutX (function "wingSurfaces"), whole width, seen from behind
wake = vtk('wingSurfaces', 'wakePlane')
if wake is not None:
    body(wallsFull, 0.37)
    colour(mirror(wake), 'POINTS', 'Cp0', 'Viridis (matplotlib)', 0.0, 1.0)
    # from behind (-X), centred on the symmetry plane, from the ground up to
    # three body heights
    parallel([x1, y0, zMid + 0.5 * H], [-1, 0, 0], [0, 0, 1], 1.5 * H)   # the wake goes up
    save('cp0_wake')
    clear()

# ================== 8-9: streamlines (trace_streamlines.py)
# traced by OpenFOAM (system/noseStreamlines: seeds two surface cells in
# front of the nose, over 90 % of the frontal section, one line every H/10,
# both ways: back to the inlet and forward round the body into the wake)
zc = 0.5 * (zTop + zBottom)
nose = vtk('noseStreamlines', 'tracks')
if nose is not None:
    lines = mirrorLines(nose)
    body(wallsFull, 0.34)
    colour(lines, 'POINTS', 'U', 'Fast', 0.0, 1.5 * uInf,
           component='Magnitude')
    # 3D view from the front, above and beside (direction from the trace,
    # almost parallel: recorded with a narrow view angle); image centre
    # moved downstream to include the near wake
    parallel([centre[0] + 0.31 * L, y0 - 0.07 * L, centre[2] + 0.13 * L],
             [0.403, 0.770, -0.495], [0.172, 0.467, 0.867], 0.58 * L)
    save('streamlines')
    # from above (-Z): body and wake, centred 0.4 L behind the body centre
    parallel([centre[0] + 0.4 * L, y0, centre[2]], [0, 0, -1], [0, 1, 0],
             0.63 * L)
    save('streamlines_top')
    # from below (+Z, the ground is not drawn): underbody, legs and the flow
    # under the body; same frame as from above, flow again from left to
    # right. In perspective (30 deg, distance for the same image height):
    # the legs are vertical, and in a parallel view along their axis they
    # would collapse into circles of zero area and disappear
    half = 0.63 * L
    perspective([centre[0] + 0.4 * L, y0, centre[2]], [0, 0, 1],
                half / math.tan(math.radians(15.0)), [0, -1, 0])
    save('streamlines_bottom')

    # --------------------------- 10: the wake vortices, side view
    # the recirculation behind the base is closed: the lines from the nose
    # pass round it. The lines of system/wakeStreamlines (seeds H/4 behind
    # the base, over its lower half and the middle half of the width) draw
    # the two recirculating vortices; seen from the side (+Y) together with
    # the lines from the nose
    base = vtk('wakeStreamlines', 'tracks')
    if base is not None:
        colour(mirrorLines(base), 'POINTS', 'U', 'Fast', 0.0, 1.5 * uInf,
               component='Magnitude')
    # side view: from H/2 before the nose to 2.5 H behind the base (near
    # wake), centred on the body
    xa, xb = x0 - 0.5 * H, x1 + 2.5 * H
    parallel([0.5 * (xa + xb), y0, zMid], [0, 1, 0], [0, 0, 1],
             (xb - xa) / (2 * ASPECT))
    save('streamlines_side')
    clear()


# ================== 11-12: vortical structures, iso-surfaces
# cut by OpenFOAM (system/isoSurfaces), mirrored; coloured by |U|/Uinf
def iso(name):
    """iso-surface <name> of system/isoSurfaces (whole body), with
    |U|/Uinf; None if missing"""
    s = vtk('isoSurfaces', name)
    if s is None:
        return None
    return Calculator(Input=mirror(s), ResultArrayName='Umag',
                      Function='mag(U)/%g' % uInf)


# from behind, above and beside (as cp_rear), farther away to include the
# near wake
rear = ([x1, y0, zc], [-1.597, -1.399, -0.747], 2.6 * L,
        [-0.272, -0.193, 0.943])
for name, image in (('Q', 'q_iso'), ('Cp0', 'cp0_iso')):
    surface = iso(name)
    if surface is None:
        continue
    body(wallsFull, 0.6)
    colour(surface, 'POINTS', 'Umag', 'Turbo', 0.0, 1.5)
    perspective(*rear)
    save(image)
    clear()

# ================== 13: planes (guidelines of the course: 16 per direction)
# cut by OpenFOAM (system/planes) across the near box (scripts/postSetup:
# planeXmin..planeZmax), at the centres of nPlanes equal intervals; |U|, k,
# p, Cp0 on each; x and z planes on the whole domain (mirrored), y planes on
# the half model (y >= y0)
N = int(var.get('nPlanes', 16))
if os.environ.get('REPORT_PLANES', '1') != '0' and N > 0:
    SMALL = [960, 540]
    lo = {d: var['plane%smin' % d.upper()] for d in 'xyz'}
    hi = {d: var['plane%smax' % d.upper()] for d in 'xyz'}
    pos = {d: [lo[d] + (i + 0.5) * (hi[d] - lo[d]) / N for i in range(N)]
           for d in 'xyz'}
    plane = {(d, i): vtk('planes', '%s%02d' % (d, i + 1))
             for d in 'xyz' for i in range(N)}
    # k: one scale for all the planes, from 0 to its maximum on the x planes
    kMax = max([plane['x', i].PointData['k'].GetRange()[1]
                for i in range(N) if plane['x', i] is not None] or [0.0]) \
        or uInf ** 2
    q = 0.5 * uInf ** 2            # p from -q/2 to q/2: Cp from -0.5 to 0.5
    fields = [('U', 'Turbo', 0.0, 1.5 * uInf, 'Magnitude', '$|U|$ [m/s]'),
              ('k', 'Viridis (matplotlib)', 0.0, kMax, None,
               '$k$ [m$^2$/s$^2$]'),
              ('p', 'Cool to Warm', -0.5 * q, 0.5 * q, None,
               '$p/\\rho$ [m$^2$/s$^2$]'),
              ('Cp0', 'Viridis (matplotlib)', 0.0, 1.0, None,
               '$C_{p,\\mathrm{tot}}$')]
    # cameras: parallel, the box (mirrored for x and z) fills the image
    xc = 0.5 * (lo['x'] + hi['x'])
    asp = SMALL[0] / SMALL[1]
    # half heights of the x and y views; their bottom edge on the ground
    hx = 1.05 * max(0.5 * (hi['z'] - lo['z']), (hi['y'] - y0) / asp)
    hy = 1.05 * max(0.5 * (hi['z'] - lo['z']), 0.5 * (hi['x'] - lo['x']) / asp)
    cam = {
        'x': lambda c: ([c, y0, lo['z'] + hx], [-1, 0, 0], [0, 0, 1], hx),
        'y': lambda c: ([xc, c, lo['z'] + hy], [0, 1, 0], [0, 0, 1], hy),
        'z': lambda c: ([xc, y0, c], [0, 0, -1], [0, 1, 0],
                        1.05 * max(hi['y'] - y0,
                                   0.5 * (hi['x'] - lo['x']) / asp))}
    label = {'x': lambda c: '$x/L$ = %.2f' % ((c - x0) / L),
             'y': lambda c: '$y/W$ = %.2f' % ((c - y0) / W),
             'z': lambda c: '$z/H$ = %.2f' % ((c - ground) / H)}
    # colour scales: set once for all the planes (applying a colour map is
    # slow, and the scales are the same on every plane)
    for f, preset, vmin, vmax, comp, _ in fields:
        lut = GetColorTransferFunction(f)
        lut.ApplyPreset(preset if hasPreset(preset) else 'Cool to Warm', True)
        lut.AutomaticRescaleRangeMode = 'Never'
        lut.RescaleTransferFunction(vmin, vmax)
        lut.UseLogScale = 0
        lut.UseBelowRangeColor = lut.UseAboveRangeColor = 0
        bar = GetScalarBar(lut, view)
        bar.TitleColor = bar.LabelColor = [0.0, 0.0, 0.0]
        bar.Visibility = 0
    done = []                      # (direction, plane) with images
    for d in 'xyz':
        for i, c in enumerate(pos[d]):
            if plane[d, i] is None:
                continue
            cut = plane[d, i] if d == 'y' else mirror(plane[d, i])
            rep = Show(cut, view)
            rep.SetRepresentationType('Surface')
            parallel(*cam[d](c))
            for f, _, vmin, vmax, comp, _ in fields:
                ColorBy(rep, ('POINTS', f) if comp is None
                        else ('POINTS', f, comp))
                # ColorBy may rescale to the data of this plane: the fixed
                # scale again, and only the legend of this field
                lut = GetColorTransferFunction(f)
                lut.RescaleTransferFunction(vmin, vmax)
                GetOpacityTransferFunction(f).RescaleTransferFunction(vmin,
                                                                      vmax)
                for g in fields:
                    if g[0] != f:
                        GetScalarBar(GetColorTransferFunction(g[0]),
                                     view).Visibility = 0
                rep.SetScalarBarVisibility(view, True)
                name = 'plane_%s_%s_%02d' % (d, f, i + 1)
                SaveScreenshot(os.path.join(outDir, name + '.png'), view,
                               ImageResolution=SMALL)
            Hide(cut, view)
            done.append((d, i))
        print('  written plane_%s_*.png (%d planes)'
              % (d, sum(1 for e in done if e[0] == d)))
    for g in fields:
        GetScalarBar(GetColorTransferFunction(g[0]), view).Visibility = 0
    # LaTeX: one page per direction and field, 4 images per row
    with open(os.path.join(outDir, 'planes.tex'), 'w') as tex:
        for d in 'xyz':
            ids = [i for e, i in done if e == d]
            if not ids:
                continue
            for f, _, _, _, _, title in fields:
                tex.write('\\begin{figure}[p]\\centering\\setlength\\tabcolsep{1pt}\n'
                          '\\begin{tabular}{cccc}\n')
                for r in range(0, len(ids), 4):
                    row = ids[r:r + 4]
                    tex.write(' & '.join(
                        '\\includegraphics[width=0.245\\textwidth]'
                        '{fig/plane_%s_%s_%02d.png}' % (d, f, i + 1)
                        for i in row) + ' \\\\\n')
                    tex.write(' & '.join('\\scriptsize ' + label[d](pos[d][i])
                                         for i in row) + ' \\\\[2pt]\n')
                tex.write('\\end{tabular}\n\\caption{%s on %d planes normal '
                          'to $%s$.}\n\\end{figure}\n\\clearpage\n'
                          % (title, len(ids), d))
    print('  written planes.tex')
