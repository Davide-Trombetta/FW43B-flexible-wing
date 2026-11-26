#!/usr/bin/env python3
"""
t026STL.py -- STL of the Tyrrell 026 single-element wing with endplates
(Zerihan 2001, PhD thesis, University of Southampton; Zerihan & Zhang 2000).
Standard library only (no numpy/scipy): it runs with the python3 of the
cluster inside Allrun.pre.

Usage:  python3 t026STL.py <output.stl> [--alpha <deg>] [--coords t026_table2.csv]

Geometry (thesis, section 2.2.1 and Table 2):
  profile     Table 2, "single element wing coordinates at reference
              incidence of 1 deg" (x/c, y/c of suction and pressure surface),
              i.e. already rotated to the true incidence of 3.45 deg
  chord       c = 223.4 mm, constant; finite trailing edge 0.007c (1.65 mm)
  span        b = 1100 mm (between the inner faces of the endplates)
  endplates   250 x 100 x 4 mm, rectangular, parallel to the ground at the
              reference incidence. Their position relative to the wing is
              not given in numbers: it is measured on Figure 1a of the thesis
              (leading edge 2 mm behind the front of the endplate, bottom edge
              20.8 mm below the leading edge, i.e. about 3 mm below the lowest
              point of the suction surface).
  incidence   quoted alpha (thesis convention) = true incidence - 2.45 deg;
              --alpha rotates wing and endplates together about the quarter
              chord by (alpha - 1) deg; a larger alpha gives more downforce
              (leading edge down).

Output frame (as the FW43B case): x along the flow, y along the span (centred
on y = 0), z up; the suction surface faces down (-z): downforce wing, lift
coefficient of forceCoeffs (liftDir 0 0 1) negative. Two solids:
  wing       extruded profile, closed at the trailing edge and at the tips;
             it goes 2 mm into each endplate (closed, overlapping solids)
  endplates  the two plates
Units: metres.
"""
import argparse
import math
import os

C = 0.2234            # chord [m]
B = 1.100             # span [m]
EP_L, EP_H, EP_T = 0.250, 0.100, 0.004      # endplate length, height, thickness [m]
EP_X0 = -0.002        # front of the endplate from the leading edge [m] (Figure 1a)
EP_Z0 = -0.0208       # bottom of the endplate from the leading edge [m] (Figure 1a)
TIP_IN = 0.002        # the wing goes this far into each endplate [m]


# ---------------------------------------------------------------- vectors
def sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# ---------------------------------------------------------------- PCHIP
def pchip(x, y):
    """monotone piecewise cubic Hermite interpolant (Fritsch-Carlson, as
    scipy.interpolate.PchipInterpolator, endpoint rule included)"""
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    m = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    d = [0.0] * n
    for k in range(1, n - 1):
        if m[k - 1] * m[k] <= 0:
            d[k] = 0.0
        else:
            w1, w2 = 2 * h[k] + h[k - 1], h[k] + 2 * h[k - 1]
            d[k] = (w1 + w2) / (w1 / m[k - 1] + w2 / m[k])

    def edge(h0, h1, m0, m1):
        dd = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
        if math.copysign(1, dd) != math.copysign(1, m0):
            return 0.0
        if math.copysign(1, m0) != math.copysign(1, m1) and abs(dd) > abs(3 * m0):
            return 3 * m0
        return dd
    d[0] = edge(h[0], h[1], m[0], m[1])
    d[-1] = edge(h[-1], h[-2], m[-1], m[-2])

    def f(t):
        i = 0
        while i < n - 2 and t > x[i + 1]:
            i += 1
        s = (t - x[i]) / h[i]
        h00 = (1 + 2 * s) * (1 - s) ** 2
        h10 = s * (1 - s) ** 2
        h01 = s * s * (3 - 2 * s)
        h11 = s * s * (s - 1)
        return h00 * y[i] + h10 * h[i] * d[i] + h01 * y[i + 1] + h11 * h[i] * d[i + 1]
    return f


def profile(csv, n=160):
    """closed profile in the x-z plane: suction surface from the trailing
    edge to the leading edge, pressure surface back to the trailing edge;
    points clustered at both edges (cosine spacing)"""
    rows = []
    with open(csv) as fh:
        next(fh)
        for line in fh:
            if line.strip():
                rows.append([float(v) for v in line.split(',')])
    xs = [r[0] for r in rows]; ys = [r[1] for r in rows]
    xp = [r[2] for r in rows]; yp = [r[3] for r in rows]
    # interpolation in sqrt(x): smooth round the leading edge
    fs = pchip([math.sqrt(v) for v in xs], ys)
    fp = pchip([math.sqrt(v) for v in xp], yp)
    xn = [0.5 * (1 - math.cos(math.pi * k / (n - 1))) for k in range(n)]
    lower = [(x, fs(math.sqrt(x))) for x in xn][::-1]     # TE -> LE, suction (lower)
    upper = [(x, fp(math.sqrt(x))) for x in xn][1:]       # LE -> TE, pressure (upper)
    return [(p[0] * C, p[1] * C) for p in lower + upper]


def area2(P):
    return 0.5 * sum(P[i][0] * P[i - len(P) + 1][1] - P[i - len(P) + 1][0] * P[i][1] for i in range(len(P)))


def earclip(P):
    """triangulation of a simple polygon (counter-clockwise), ear clipping"""
    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def inside(p, a, b, c):
        return cr(a, b, p) >= 0 and cr(b, c, p) >= 0 and cr(c, a, p) >= 0
    idx = list(range(len(P))); tris = []
    while len(idx) > 3:
        for k in range(len(idx)):
            i0, i1, i2 = idx[k - 1], idx[k], idx[(k + 1) % len(idx)]
            a, b, c = P[i0], P[i1], P[i2]
            if cr(a, b, c) <= 1e-14:
                continue
            if any(inside(P[j], a, b, c) for j in idx if j not in (i0, i1, i2)):
                continue
            tris.append((i0, i1, i2)); idx.pop(k); break
        else:
            raise SystemExit('triangulation of the profile failed')
    tris.append(tuple(idx))
    return tris


def extrude(P, y0, y1):
    """closed prism of the polygon P (x, z) between y0 and y1"""
    n = len(P); T = []
    A = [(p[0], y0, p[1]) for p in P]; Bv = [(p[0], y1, p[1]) for p in P]
    for i in range(n):
        j = (i + 1) % n
        T.append((A[i], A[j], Bv[j])); T.append((A[i], Bv[j], Bv[i]))
    for a, b, c in earclip(P):
        T.append((Bv[a], Bv[b], Bv[c])); T.append((A[a], A[c], A[b]))
    return T


def orient_outward(T):
    vol = sum(dot(t[0], cross(t[1], t[2])) for t in T) / 6
    return ([t[::-1] for t in T] if vol < 0 else T), abs(vol)


def box(x0, x1, y0, y1, z0, z1):
    return extrude([(x0, z0), (x1, z0), (x1, z1), (x0, z1)], y0, y1)


def rotate(T, dalpha, xc, zc):
    """rotation about the y axis through (xc, zc): dalpha > 0 lowers the
    leading edge (more incidence for the inverted wing)"""
    ca, sa = math.cos(math.radians(dalpha)), math.sin(math.radians(dalpha))

    def r(p):
        X, Z = p[0] - xc, p[2] - zc
        return (xc + X * ca - Z * sa, p[1], zc + X * sa + Z * ca)
    return [tuple(r(p) for p in t) for t in T]


def write(f, name, T):
    f.write(f'solid {name}\n')
    for t in T:
        nv = cross(sub(t[1], t[0]), sub(t[2], t[0]))
        ln = math.sqrt(dot(nv, nv)) or 1e-30
        f.write(f'  facet normal {nv[0] / ln:.7e} {nv[1] / ln:.7e} {nv[2] / ln:.7e}\n    outer loop\n')
        for p in t:
            f.write(f'      vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n')
        f.write('    endloop\n  endfacet\n')
    f.write(f'endsolid {name}\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('output')
    ap.add_argument('--alpha', type=float, default=1.0, help='quoted incidence [deg] (thesis convention, reference 1)')
    ap.add_argument('--coords', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 't026_table2.csv'))
    a = ap.parse_args()

    P = profile(a.coords)
    if area2(P) < 0:
        P = P[::-1]
    h = B / 2
    wing, vw = orient_outward(extrude(P, -h - TIP_IN, h + TIP_IN))
    plates = []
    for s in (-1, 1):
        y0, y1 = (h, h + EP_T) if s > 0 else (-h - EP_T, -h)
        plates += orient_outward(box(EP_X0, EP_X0 + EP_L, y0, y1, EP_Z0, EP_Z0 + EP_H))[0]

    # quarter chord on the chord line (leading edge (0, 0) - mid trailing edge)
    te = (0.5 * (P[0][0] + P[-1][0]), 0.5 * (P[0][1] + P[-1][1]))
    xc, zc = 0.25 * te[0], 0.25 * te[1]
    d = a.alpha - 1.0
    wing, plates = rotate(wing, d, xc, zc), rotate(plates, d, xc, zc)

    with open(a.output, 'w') as f:
        write(f, 'wing', wing)
        write(f, 'endplates', plates)
    pts = [p for t in wing + plates for p in t]
    lo = [round(min(p[i] for p in pts), 4) for i in range(3)]
    hi = [round(max(p[i] for p in pts), 4) for i in range(3)]
    print(f'{a.output}: wing {len(wing)} triangles (volume {vw * 1e6:.1f} cm3), endplates {len(plates)} triangles')
    print(f'   quoted alpha {a.alpha:g} deg = true incidence {a.alpha + 2.45:g} deg')
    print(f'   bounding box {lo} - {hi} m')
    print(f'   reference: c = {C} m, b = {B} m, S = c b = {C * B:.5f} m2, '
          f'lowest point of the wing z = {min(p[2] for t in wing for p in t):.4f} m')


if __name__ == '__main__':
    main()
