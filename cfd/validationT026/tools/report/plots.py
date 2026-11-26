# plots.py -- charts of the report, made with ParaView (run by scripts/report)
#
# Usage, from the case directory:
#     pvbatch tools/report/plots.py <dataDir> <figDir> <from> <CdMean> <ClMean>
#                                   [CdReference]
#
# In <dataDir> (CSV files prepared by scripts/report):
#     coeffs.csv          iteration, Cd, Cl, ...   (postProcessing/forceCoeffs)
#     residuals.csv       iteration, p, Ux, ...    (postProcessing/solverInfo)
#     cp_section.csv      xc = x/c, minusCp = -Cp on the wing section
#                         y = cutY (from tools/report/images.py)
#     inlet.csv           iteration, CpInlet       (postProcessing/inletPressure)
#     profile1..5.csv     zH, uU, wU, Cp0 on the vertical wake lines
#     centreline.csv      xH, uU, Cp on the wake line (profiles)
#
# Output (PNG) in <figDir>: residuals, coeffs, cp_centreline, inlet,
#     profiles_u, profiles_w, centreline
# (the probe charts of the Ahmed-body case have been removed)
# (a chart whose CSV file is missing is skipped)
#
# Axis ranges come from the data: the history of the coefficients is shown
# around the averaging window (from <from> on), enlarged by its own width on
# each side, because the first iterations are far off.

import csv
import math
import os
import sys

from paraview.simple import *

dataDir, figDir = sys.argv[1], sys.argv[2]
start, cdMean, clMean = (float(x) for x in sys.argv[3:6])
cdRef = float(sys.argv[6]) if len(sys.argv) > 6 else None
colours = [[0.12, 0.47, 0.71], [0.84, 0.15, 0.16], [0.17, 0.63, 0.17],
           [1.00, 0.50, 0.05], [0.58, 0.40, 0.74], [0.55, 0.34, 0.29],
           [0.89, 0.47, 0.76]]


def columns(name):
    """column names and values of a CSV file of dataDir"""
    with open(os.path.join(dataDir, name)) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    data = {h: [float(r[i]) for r in rows[1:] if r[i] not in ('', 'N/A')]
            for i, h in enumerate(head)}
    return head, data


def chart(xTitle, yTitle, log=False):
    view = CreateView('XYChartView')
    view.ViewSize = [1600, 800]
    view.LegendLocation = 'TopRight'
    view.BottomAxisTitle = xTitle
    view.LeftAxisTitle = yTitle
    view.LeftAxisLogScale = 1 if log else 0
    for a in ('LeftAxis', 'BottomAxis', 'RightAxis'):
        setattr(view, a + 'TitleFontSize', 22)
        setattr(view, a + 'LabelFontSize', 18)
    view.LegendFontSize = 20
    return view


def series(view, csvName, x, ys, labels=None, corner=None, points=False,
           colour=None, dashed=()):
    """plot the columns ys against x from a CSV file of dataDir; colour: one
    colour for all ys (otherwise the next ones of the list); dashed: the ys
    drawn with a dashed line"""
    reader = CSVReader(FileName=[os.path.join(dataDir, csvName)])
    d = Show(reader, view, 'XYChartRepresentation')
    d.UseIndexForXAxis = 0
    d.XArrayName = x
    d.SeriesVisibility = ys
    labels = labels or ys
    d.SeriesLabel = [v for y, l in zip(ys, labels) for v in (y, l)]
    d.SeriesColor = [str(v) for i, y in enumerate(ys)
                     for v in [y] + (colour or colours[
                         (len(view.Representations) - 1 + i) % len(colours)])]
    d.SeriesLineThickness = [v for y in ys for v in (y, '3')]
    d.SeriesLineStyle = [v for y in ys for v in (y, '2' if y in dashed else '1')]
    if points:
        d.SeriesLineStyle = [v for y in ys for v in (y, '0')]
        d.SeriesMarkerStyle = [v for y in ys for v in (y, '4')]
    if corner is not None:
        d.SeriesPlotCorner = [v for y in ys for v in (y, str(corner))]
    return d


def save(view, name):
    SaveScreenshot(os.path.join(figDir, name + '.png'), view,
                   ImageResolution=view.ViewSize)
    print('  written ' + name + '.png')


# ------------------------------------------------------------ residuals
# A logarithmic axis cannot show zero, negative or non-numeric values (a
# solver that did not run in an iteration, a nan after a blow-up): ParaView
# then draws no line at all. The chart reads a copy of the file with those
# rows left out, and the axis spans the values actually drawn.
def positiveCopy(name):
    """copy of a CSV file of dataDir without the rows that have a value
    that is not a finite number > 0 (iteration column excluded); returns
    its name and the min, max of the values kept"""
    with open(os.path.join(dataDir, name)) as f:
        rows = list(csv.reader(f))
    keep, lo, hi = [rows[0]], None, None
    for r in rows[1:]:
        try:
            v = [float(x) for x in r[1:]]
        except ValueError:
            continue
        if all(math.isfinite(x) and x > 0 for x in v):
            keep.append(r)
            lo = min(v) if lo is None else min(lo, min(v))
            hi = max(v) if hi is None else max(hi, max(v))
    out = name.replace('.csv', '_log.csv')
    with open(os.path.join(dataDir, out), 'w', newline='') as f:
        csv.writer(f).writerows(keep)
    return out, lo, hi


head, _ = columns('residuals.csv')
name, lo, hi = positiveCopy('residuals.csv')
v = chart('Iteration', 'Initial residual', log=True)
if lo is not None:
    v.LeftAxisUseCustomRange = 1
    v.LeftAxisRangeMinimum = 10 ** math.floor(math.log10(lo))
    v.LeftAxisRangeMaximum = 10 ** math.ceil(math.log10(hi))
series(v, name, head[0], head[1:])
save(v, 'residuals')

# ------------------------------------------- history of Cd and Cl
# the first tenth of the iterations is left out: it is the start from the
# uniform initial field, with values far from the final ones that would
# squeeze the rest of the history. The y axis contains all the values shown
# (and the experimental Cd), with a margin of a tenth of their span.
head, c = columns('coeffs.csv')
it = head[0]
first = c[it][0] + 0.1 * (max(c[it]) - c[it][0])
shown = [y for name in ('Cd', 'Cl') for x, y in zip(c[it], c[name])
         if x >= first] + ([cdRef] if cdRef is not None else [])
lo, hi = min(shown), max(shown)
pad = 0.1 * (hi - lo or abs(hi) or 1.0)
v = chart('Iteration', 'Coefficient')
v.BottomAxisUseCustomRange = 1
v.BottomAxisRangeMinimum, v.BottomAxisRangeMaximum = first, max(c[it])
v.LeftAxisUseCustomRange = 1
v.LeftAxisRangeMinimum, v.LeftAxisRangeMaximum = lo - pad, hi + pad
series(v, 'coeffs.csv', it, ['Cd', 'CdMean'], ['Cd', 'Cd mean'],
       colour=colours[0], dashed=('CdMean',))
series(v, 'coeffs.csv', it, ['Cl', 'ClMean'], ['Cl', 'Cl mean'],
       colour=colours[1], dashed=('ClMean',))
if cdRef is not None:
    series(v, 'coeffs.csv', it, ['CdReference'], ['Cd experimental'],
           colour=[0, 0, 0], dashed=('CdReference',))
save(v, 'coeffs')

# ------------------------------------- Cp on the wing section
# mainplane and flap, suction and pressure sides: one set of points (the
# two sides of a wing cannot be split by height as the roof and underbody
# of a car body)
if os.path.exists(os.path.join(dataDir, 'cp_section.csv')):
    v = chart('x / c', '-Cp')
    series(v, 'cp_section.csv', 'xc', ['minusCp'],
           ['wing section'], points=True)
    save(v, 'cp_centreline')


# ---------------------- history at the inlet
# from the end of the first tenth, as the coefficients


def history(csvName, yTitle, ys, labels, name):
    if not os.path.exists(os.path.join(dataDir, csvName)):
        return
    h, d = columns(csvName)
    v = chart('Iteration', yTitle)
    v.BottomAxisUseCustomRange = 1
    v.BottomAxisRangeMinimum, v.BottomAxisRangeMaximum = first, max(d[h[0]])
    # vertical range from the iterations shown only: the first ones (the
    # start from the uniform field, or a blow-up) would flatten the curves
    with open(os.path.join(dataDir, csvName)) as f:
        rows = list(csv.reader(f))
    ix = [rows[0].index(y) for y in ys if y in rows[0]]
    vals = []
    for r in rows[1:]:
        try:
            if float(r[0]) >= first:
                vals += [float(r[i]) for i in ix if r[i] not in ('', 'N/A')]
        except ValueError:
            pass
    vals = [x for x in vals if math.isfinite(x)]
    if vals:
        lo, hi = min(vals), max(vals)
        pad = 0.05 * (hi - lo) if hi > lo else 0.05 * (abs(hi) or 1)
        v.LeftAxisUseCustomRange = 1
        v.LeftAxisRangeMinimum, v.LeftAxisRangeMaximum = lo - pad, hi + pad
    series(v, csvName, h[0], ys, labels)
    save(v, name)


history('inlet.csv', 'Cp on the inlet', ['CpInlet'], ['Cp inlet'], 'inlet')

# ------------------- profiles behind the wing (x/H from the trailing point)
# z/H against the velocity: one line per position x/H
positions = ['0.25', '0.5', '1', '1.5', '2']
for comp, title, name in (('uU', 'u / Uinf', 'profiles_u'),
                          ('wU', 'w / Uinf', 'profiles_w')):
    files = [('profile%d.csv' % (i + 1), 'x/H = ' + x)
             for i, x in enumerate(positions)
             if os.path.exists(os.path.join(dataDir, 'profile%d.csv' % (i + 1)))]
    if not files:
        continue
    v = chart(title, 'z / H')
    for f, lab in files:
        series(v, f, comp, ['zH'], [lab])
    save(v, name)

# ------------------------------------------------- wake centre line
if os.path.exists(os.path.join(dataDir, 'centreline.csv')):
    v = chart('x / H (from the base)', 'u / Uinf, Cp')
    series(v, 'centreline.csv', 'xH', ['uU', 'Cp'], ['u / Uinf', 'Cp'])
    save(v, 'centreline')
