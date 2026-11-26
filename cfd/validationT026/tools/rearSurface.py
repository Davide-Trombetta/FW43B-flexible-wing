#!/usr/bin/env python3
# rearSurface.py -- extracts from the STL the rear upper surface that slopes
#                   down (the slant of the Ahmed body) and the end of the roof
#                   before it, for a distance refinement in snappyHexMesh
#
# Usage (scripts/rearRefinement):
#     python3 tools/rearSurface.py <stl> <region> <upstream [m]> <out.stl>
# Prints shell assignments: REAR_FACES, REAR_AREA, REAR_XMIN, REAR_XMAX,
# REAR_ANGLE (mean inclination of the slant from the horizontal, degrees).
#
# Faces of <region>, with the outward normal n (made outward from the body if
# the STL is oriented inwards):
#   slant  n_x > 0.1 and n_z > 0.1 (facing downstream and up), in the rear
#          half of the region: the flow separates at its upstream edge
#   roof   n_z > 0.95, from <upstream> before the slant edge to the edge:
#          the boundary layer that reaches the edge decides the separation
# ASCII STL only (as written by CAD tools for OpenFOAM: "solid <name>").

import math
import sys

stl, region, upstream, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


# triangles of the region
tris, solid, verts = [], None, []
with open(stl, errors='replace') as f:
    for line in f:
        w = line.split()
        if not w:
            continue
        if w[0] == 'solid':
            solid = w[1] if len(w) > 1 else ''
        elif w[0] == 'vertex' and solid == region:
            verts.append([float(w[1]), float(w[2]), float(w[3])])
            if len(verts) == 3:
                tris.append(verts)
                verts = []
if not tris:
    sys.exit('no ASCII triangles of region %s in %s' % (region, stl))

# area vectors (right-hand rule) and centres
data = []
for t in tris:
    a = cross(sub(t[1], t[0]), sub(t[2], t[0]))
    area = 0.5 * math.sqrt(sum(c * c for c in a))
    if area == 0.0:
        continue
    n = [c / (2 * area) for c in a]
    c = [(t[0][i] + t[1][i] + t[2][i]) / 3 for i in range(3)]
    data.append((t, n, c, area))
# outward: normals pointing away from the centre of the region on average
cen = [sum(d[2][i] * d[3] for d in data) / sum(d[3] for d in data) for i in range(3)]
s = sum(d[3] * sum(d[1][i] * (d[2][i] - cen[i]) for i in range(3)) for d in data)
sign = 1.0 if s > 0 else -1.0
xmin = min(d[2][0] for d in data)
xmax = max(d[2][0] for d in data)
xmid = 0.5 * (xmin + xmax)

slant = [d for d in data
         if sign * d[1][0] > 0.1 and sign * d[1][2] > 0.1 and d[2][0] > xmid]
if not slant:
    print('REAR_FACES=0')
    sys.exit(0)
sx = min(min(v[0] for v in d[0]) for d in slant)          # slant edge
roof = [d for d in data
        if sign * d[1][2] > 0.95 and sx - upstream <= d[2][0] <= sx]
faces = slant + roof

with open(out, 'w') as f:
    f.write('solid rear\n')
    for t, n, c, area in faces:
        f.write('  facet normal %g %g %g\n    outer loop\n'
                % tuple(sign * x for x in n))
        for v in t:
            f.write('      vertex %.9g %.9g %.9g\n' % tuple(v))
        f.write('    endloop\n  endfacet\n')
    f.write('endsolid rear\n')

aSlant = sum(d[3] for d in slant)
nx = sum(sign * d[1][0] * d[3] for d in slant) / aSlant
nz = sum(sign * d[1][2] * d[3] for d in slant) / aSlant
print('REAR_FACES=%d' % len(faces))
print('REAR_AREA=%.6g' % sum(d[3] for d in faces))
print('REAR_XMIN=%.6g' % min(min(v[0] for v in d[0]) for d in faces))
print('REAR_XMAX=%.6g' % max(max(v[0] for v in d[0]) for d in faces))
print('REAR_ANGLE=%.3g' % math.degrees(math.atan2(nx, nz)))
