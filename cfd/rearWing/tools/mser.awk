#------------------------------------------------------------------------------
#  mser.awk -- start of the stable part of a history (MSER rule) and the mean
#              of Cd and Cl over it with their standard errors (batch means)
#
#  Usage: awk -v cd=<column> -v cl=<column> [-v nb=10] [-v drift=1] [-v from=<it>]
#             -f tools/mser.awk <file>
#         <file>: one line per iteration (scripts/history forceCoeffs ...),
#         iteration in column 1; for forceCoeffs.dat cd=3, cl=4
#  Prints one line (for the shell "read"):
#     FROM FROMCD FROMCL UNSTABLE N CDMEAN CDSE CLMEAN CLSE
#  -v from=<iteration>: the stable part starts there (no MSER search)
#  with -v drift=1 four more fields at the end:
#     ... CDDRIFT CDSIG CLDRIFT CLSIG
#  FROM      start of the stable part (iteration): the later of FROMCD, FROMCL
#  UNSTABLE  1 if the start is at the end of the search range (see below)
#  N         iterations from FROM to the end
#  *MEAN     mean over them;  *SE  standard error of that mean
#
#  MSER (Marginal Standard Error Rule, White 1997): the first iteration d
#  that minimises  sum_{i>=d} (y_i - mean_d)^2 / (n - d)^2, i.e. the
#  uncertainty of the mean of what is left, searched in the first half of the
#  history as the rule prescribes; if the minimum is at the end of the first
#  half the transient is not over and the run is too short (UNSTABLE = 1).
#  Standard error: method of batch means. Consecutive iterations of SIMPLE
#  are strongly correlated, so sigma / sqrt(N) would be far too small; the
#  stable part is cut in nb equal batches and SE = std(batch means)/sqrt(nb)
#  (nb = 10, a usual choice: 10-30 batches).
#  Drift (-v drift=1): mean of the second half of the stable part minus mean
#  of its first half (from the batch means), and whether that difference is
#  significant: Welch t test on the two halves of the batch means, 95 %.
#  A part that is still moving slowly can look flat to MSER and give a small
#  standard error; the two halves of a stationary part (also one with slow
#  oscillations) have the same mean within their scatter, those of a drift
#  do not.
#  The single implementation used by scripts/report and
#  scripts/convergenceMonitor.
#------------------------------------------------------------------------------
function mser(a,   i, s, q, S, Q, d, m, v, best, k) {
    s = 0; q = 0
    for (i = n; i >= 1; i--) { s += a[i]; q += a[i] * a[i]; S[i] = s; Q[i] = q }
    best = -1
    for (d = 1; d <= int(n / 2) + 1; d++) {
        m = n - d + 1; v = (Q[d] - S[d] * S[d] / m) / (m * m)
        if (best < 0 || v < best) { best = v; k = d }
    }
    return k
}
function tquant(nu,   z) {     # Student t, two-sided 95 % (Cornish-Fisher expansion)
    z = 1.959964
    return z + (z^3 + z) / (4 * nu) + (5 * z^5 + 16 * z^3 + 3 * z) / (96 * nu^2) \
             + (3 * z^7 + 19 * z^5 + 17 * z^3 - 15 * z) / (384 * nu^3)
}
function halves(a, d,   b, h, k, i, bm, s1, s2, q1, q2, m1, m2, v1, v2) {
    # sets DRIFT = mean of the 2nd half of the batch means - mean of the 1st
    # half, and SIG = 1 if that difference is significant (t test, 95 %)
    DRIFT = 0; SIG = 0
    b = int((n - d + 1) / nb); h = int(nb / 2)
    if (b < 1 || h < 2) return
    s1 = s2 = q1 = q2 = 0
    for (k = 0; k < 2 * h; k++) {
        bm = 0
        for (i = 0; i < b; i++) bm += a[d + (k + nb - 2 * h) * b + i]
        bm /= b
        if (k < h) { s1 += bm; q1 += bm * bm } else { s2 += bm; q2 += bm * bm }
    }
    m1 = s1 / h; m2 = s2 / h
    v1 = (q1 - h * m1 * m1) / (h - 1); v2 = (q2 - h * m2 * m2) / (h - 1)
    if (v1 < 0) v1 = 0
    if (v2 < 0) v2 = 0
    DRIFT = m2 - m1
    SIG = (DRIFT * DRIFT > tquant(2 * h - 2) ^ 2 * (v1 + v2) / h) ? 1 : 0
}
function batch(a, d,   b, k, i, bm, B, mean, ss) {   # sets M and SE
    b = int((n - d + 1) / nb)
    if (b < 1) { M = 0; for (i = d; i <= n; i++) M += a[i]; M /= (n - d + 1); SE = -1; return }
    mean = 0
    for (k = 0; k < nb; k++) {
        bm = 0
        for (i = 0; i < b; i++) bm += a[d + k * b + i]
        B[k] = bm / b; mean += B[k]
    }
    mean /= nb; ss = 0
    for (k = 0; k < nb; k++) ss += (B[k] - mean) ^ 2
    M = mean; SE = sqrt(ss / (nb - 1)) / sqrt(nb)
}
BEGIN { if (nb == "") nb = 10 }
{ n++; it[n] = $1; y[n] = $cd; z[n] = $cl }
END {
    if (n < 2) { print "0 0 0 1 0 0 -1 0 -1" (drift ? " 0 1 0 1" : ""); exit }
    dd = mser(y); dl = mser(z); d = (dl > dd) ? dl : dd
    if (from != "") {                      # start given: no MSER
        for (d = 1; d < n && it[d] < from; d++) ; dd = dl = d
    }
    unstable = (d >= int(n / 2) + 1) ? 1 : 0
    batch(y, d); cdm = M; cdse = SE
    batch(z, d); clm = M; clse = SE
    printf "%s %s %s %d %d %.6g %.3g %.6g %.3g", it[d], it[dd], it[dl], unstable,
           n - d + 1, cdm, cdse, clm, clse
    if (drift) {
        halves(y, d); printf " %.3g %d", DRIFT, SIG
        halves(z, d); printf " %.3g %d", DRIFT, SIG
    }
    printf "\n"
}
