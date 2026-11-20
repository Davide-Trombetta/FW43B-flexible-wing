#!/usr/bin/env python3
# spectrum.py -- frequency content of Cd and Cl in the averaging window of a
#                URANS run (run by scripts/report)
#
# Usage, from the case directory:
#     python3 tools/report/spectrum.py <coeffs.dat> <from> <uInf> <H> <out.csv>
#
# <coeffs.dat>: the history of scripts/report (coeffHistory), one row per
#     time step: Time Cm Cd Cl Cl(f) Cl(r)
# <from>: start of the averaging window [s] (timeStart of system/fieldAverage)
# <uInf>, <H>: free-stream velocity and wing height, for the Strouhal number
#     St = f H / uInf
#
# Prints one line (for the shell "read"):
#     F_CL ST_CL AMP_CL F_CD ST_CD AMP_CD NPERIODS
# (F: frequency of the highest peak [Hz], AMP: its amplitude, NPERIODS:
#  periods of the Cl peak in the window; "-" for each value if the window
#  is too short) and writes <out.csv>: f, St, ampCd, ampCl up to five times
# the Cl peak frequency (or the Nyquist frequency).
#
# Method: the window (fixed time step, as in system/controlDict) minus its
# mean, times a Hann window, discrete Fourier transform on a frequency grid
# eight times finer than 1/T (zero padding: T = window length), single-sided
# amplitude corrected for the Hann window. The peak is searched above 2/T
# (slower changes are drift, not shedding) and refined with a parabola
# through the three grid points around it. Frequency resolution about 1/T:
# with the window of 0.2 s at 70 m/s, 5 Hz on a shedding of the order of
# 100 Hz. Plain Python, no numpy: it runs wherever python3 does (cluster).

import math
import sys


def load(path, start):
    t, cd, cl = [], [], []
    for line in open(path):
        v = line.split()
        if len(v) < 4 or v[0].startswith('#'):
            continue
        try:
            x = [float(a) for a in v[:4]]
        except ValueError:
            continue
        if x[0] >= start:
            t.append(x[0]); cd.append(x[2]); cl.append(x[3])
    return t, cd, cl


def amplitudes(y, dt, freqs):
    """single-sided amplitude of y (mean removed, Hann window) at freqs"""
    n = len(y)
    m = sum(y) / n
    w = [0.5 - 0.5 * math.cos(2 * math.pi * i / (n - 1)) for i in range(n)]
    yw = [(a - m) * b for a, b in zip(y, w)]
    norm = 2.0 / sum(w)
    out = []
    for f in freqs:
        # rotation by a fixed angle per sample (no cos/sin per point)
        c, s = math.cos(2 * math.pi * f * dt), math.sin(2 * math.pi * f * dt)
        re = im = 0.0
        cr, ci = 1.0, 0.0
        for a in yw:
            re += a * cr
            im -= a * ci
            cr, ci = cr * c - ci * s, cr * s + ci * c
        out.append(norm * math.hypot(re, im))
    return out


def peak(freqs, amp, fmin):
    """frequency and amplitude of the highest peak above fmin"""
    idx = [i for i, f in enumerate(freqs) if f >= fmin and 0 < i < len(freqs) - 1]
    if not idx:
        return None, None
    k = max(idx, key=lambda i: amp[i])
    a, b, c = amp[k - 1], amp[k], amp[k + 1]
    den = a - 2 * b + c
    d = 0.5 * (a - c) / den if den != 0 else 0.0
    d = max(-0.5, min(0.5, d))
    df = freqs[1] - freqs[0]
    return freqs[k] + d * df, b - 0.25 * (a - c) * d


def main():
    path, start, uInf, H, out = sys.argv[1], float(sys.argv[2]), \
        float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    t, cd, cl = load(path, start)
    if len(t) < 16:
        print('- - - - - - -')
        return
    dt = (t[-1] - t[0]) / (len(t) - 1)          # fixed time step
    T = t[-1] - t[0] + dt
    df = 1.0 / (8 * T)
    fnyq = 0.5 / dt
    freqs = [i * df for i in range(int(fnyq / df) + 1)]
    acl = amplitudes(cl, dt, freqs)
    acd = amplitudes(cd, dt, freqs)
    fcl, pcl = peak(freqs, acl, 2.0 / T)
    fcd, pcd = peak(freqs, acd, 2.0 / T)
    if fcl is None:
        print('- - - - - - -')
        return
    fmax = min(fnyq, 5 * fcl)
    with open(out, 'w') as f:
        f.write('f,St,ampCd,ampCl\n')
        for fr, a, b in zip(freqs, acd, acl):
            if fr <= fmax:
                f.write('%g,%g,%g,%g\n' % (fr, fr * H / uInf, a, b))
    st = lambda fr: fr * H / uInf
    print('%.4g %.4g %.3g %.4g %.4g %.3g %.1f' % (
        fcl, st(fcl), pcl, fcd, st(fcd), pcd, fcl * T))


if __name__ == '__main__':
    main()
