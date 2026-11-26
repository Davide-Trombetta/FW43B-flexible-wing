#!/usr/bin/env python3
"""
t026STL.py -- STL of the Tyrrell 026 single-element wing with endplates
(Zerihan 2001, PhD thesis, University of Southampton; Zerihan & Zhang 2000).

Usage:  python3 t026STL.py <output.stl> [--alpha <deg>] [--coords table2.csv]

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
import os

import numpy as np
from scipy.interpolate import PchipInterpolator

C = 0.2234            # chord [m]
B = 1.100             # span [m]
EP_L, EP_H, EP_T = 0.250, 0.100, 0.004      # endplate length, height, thickness [m]
EP_X0 = -0.002        # front of the endplate from the leading edge [m] (Figure 1a)
EP_Z0 = -0.0208       # bottom of the endplate from the leading edge [m] (Figure 1a)
TIP_IN = 0.002        # the wing goes this far into each endplate [m]


def profile(csv, n=160):
    """closed profile (counter-clockwise in x-z seen from +y): suction surface
    from the trailing edge to the leading edge, pressure surface back to the
    trailing edge; points clustered at both edges (cosine spacing)"""
    t = np.loadtxt(csv, delimiter=',', skiprows=1)
    xs, ys, xp, yp = t.T
    beta = np.linspace(0, np.pi, n)
    xn = 0.5 * (1 - np.cos(beta))
    # interpolation in sqrt(x): smooth round the leading edge
    fs = PchipInterpolator(np.sqrt(xs), ys)
    fp = PchipInterpolator(np.sqrt(xp), yp)
    zs, zp = fs(np.sqrt(xn)), fp(np.sqrt(xn))
    lower = np.column_stack([xn, zs])[::-1]           # TE -> LE, suction (lower)
    upper = np.column_stack([xn, zp])[1:]             # LE -> TE, pressure (upper)
    P = np.vstack([lower, upper]) * C
    return P            # first point: TE suction, last point: TE pressure (open)


def area2(P):
    x, z = P[:, 0], P[:, 1]
    return 0.5 * np.sum(x * np.roll(z, -1) - np.roll(x, -1) * z)


def earclip(P):
    """triangulation of a simple polygon (counter-clockwise), ear clipping"""
    idx = list(range(len(P))); tris = []
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    def inside(p, a, b, c):
        return cross(a, b, p) >= 0 and cross(b, c, p) >= 0 and cross(c, a, p) >= 0
    guard = 0
    while len(idx) > 3 and guard < 100000:
        guard += 1
        for k in range(len(idx)):
            i0, i1, i2 = idx[k - 1], idx[k], idx[(k + 1) % len(idx)]
            a, b, c = P[i0], P[i1], P[i2]
            if cross(a, b, c) <= 1e-14:
                continue
            if any(inside(P[j], a, b, c) for j in idx if j not in (i0, i1, i2)):
                continue
            tris.append((i0, i1, i2)); idx.pop(k); break
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
    return np.array(T)


def orient_outward(T):
    vol = np.einsum('ij,ij->i', T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6
    return (T[:, ::-1] if vol < 0 else T), abs(vol)


def box(x0, x1, y0, y1, z0, z1):
    P = np.array([(x0, z0), (x1, z0), (x1, z1), (x0, z1)])
    return extrude(P, y0, y1)


def rotate(T, dalpha, xc, zc):
    """rotation about the y axis through (xc, zc): dalpha > 0 lowers the
    leading edge (more incidence for the inverted wing)"""
    a = np.radians(dalpha)
    X, Z = T[..., 0] - xc, T[..., 2] - zc
    T = T.copy()
    T[..., 0] = xc + X * np.cos(a) - Z * np.sin(a)
    T[..., 2] = zc + X * np.sin(a) + Z * np.cos(a)
    return T


def write(f, name, T):
    f.write(f'solid {name}\n')
    for t in T:
        n = np.cross(t[1] - t[0], t[2] - t[0]); n /= max(np.linalg.norm(n), 1e-30)
        f.write(f'  facet normal {n[0]:.7e} {n[1]:.7e} {n[2]:.7e}\n    outer loop\n')
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
        b, _ = orient_outward(box(EP_X0, EP_X0 + EP_L, y0, y1, EP_Z0, EP_Z0 + EP_H))
        plates.append(b)
    plates = np.vstack(plates)

    # quarter chord on the chord line (leading edge (0, 0) - mid trailing edge)
    te = 0.5 * (P[0] + P[-1])
    xc, zc = 0.25 * te[0], 0.25 * te[1]
    d = a.alpha - 1.0
    wing, plates = rotate(wing, d, xc, zc), rotate(plates, d, xc, zc)

    with open(a.output, 'w') as f:
        write(f, 'wing', wing)
        write(f, 'endplates', plates)
    allp = np.vstack([wing.reshape(-1, 3), plates.reshape(-1, 3)])
    print(f'{a.output}: wing {len(wing)} triangles (volume {vw * 1e6:.1f} cm3), endplates {len(plates)} triangles')
    print(f'   quoted alpha {a.alpha:g} deg = true incidence {a.alpha + 2.45:g} deg')
    print(f'   bounding box {np.round(allp.min(0), 4)} - {np.round(allp.max(0), 4)} m')
    print(f'   reference: c = {C} m, b = {B} m, S = c b = {C * B:.5f} m2, '
          f'lowest point of the wing z = {wing[..., 2].min():.4f} m')


if __name__ == '__main__':
    main()
