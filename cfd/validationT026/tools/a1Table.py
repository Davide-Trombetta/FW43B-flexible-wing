#!/usr/bin/env python3
"""
a1Table.py -- the a1 study of the validation case: collects the
report/results/validation.csv of every run and compares CFD and experiment.

Usage:  python3 tools/a1Table.py <folder with the runs> [-o a1_study]
        e.g. python3 cfd/validationT026/tools/a1Table.py results/validation/T026/alpha11

Input:  <folder>/a1_<value>/report/results/validation.csv (written by
        scripts/report results for each run)
Output: <name>.csv  one line per run, sorted by a1
        <name>.png  CL (downforce) and CD against a1, with the experiment
        and, on the screen, the a1 with the smallest CL error and the one with
        the smallest combined error sqrt(eCL^2 + eCD^2)
"""
import argparse
import csv
import glob
import math
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder')
    ap.add_argument('-o', '--out', default='a1_study')
    ap.add_argument('--yplus', action='store_true', help='include the y+ sensitivity runs (a1_<A1>_yp<YPLUS>)')
    a = ap.parse_args()

    rows = []
    for f in glob.glob(os.path.join(a.folder, 'a1_*', 'report', 'results', 'validation.csv')):
        if '_yp' in f and not a.yplus:      # y+ sensitivity runs: only with --yplus
            continue
        with open(f) as fh:
            for r in csv.DictReader(fh):
                r['run'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(f))))
                rows.append(r)
    if not rows:
        raise SystemExit(f'no validation.csv under {a.folder}/*/report/results/')
    rows.sort(key=lambda r: float(r['a1']))

    keys = ['run', 'alpha', 'a1', 'mesh', 'cells', 'iterations', 'status', 'Cl_downforce', 'ClStdError',
            'Cd', 'CdStdError', 'CL_exp', 'CD_exp', 'ClError_pct', 'CdError_pct']
    with open(a.out + '.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)

    print(f"{'a1':>6} {'CL CFD':>8} {'CL exp':>7} {'err %':>7} {'CD CFD':>8} {'CD exp':>7} {'err %':>7}  status")
    for r in rows:
        print(f"{float(r['a1']):6.3g} {float(r['Cl_downforce']):8.4f} {float(r['CL_exp']):7.3f} "
              f"{float(r['ClError_pct']):7.1f} {float(r['Cd']):8.4f} {float(r['CD_exp']):7.3f} "
              f"{float(r['CdError_pct']):7.1f}  {r['status']}")
    best = min(rows, key=lambda r: abs(float(r['ClError_pct'])))
    comb = min(rows, key=lambda r: math.hypot(float(r['ClError_pct']), float(r['CdError_pct'])))
    print(f"smallest CL error: a1 = {best['a1']} ({float(best['ClError_pct']):+.1f} %)")
    print(f"smallest combined error: a1 = {comb['a1']} "
          f"(CL {float(comb['ClError_pct']):+.1f} %, CD {float(comb['CdError_pct']):+.1f} %)")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not found: no plot'); return
    x = [float(r['a1']) for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for k, (cfd, se, exp, lab) in enumerate((('Cl_downforce', 'ClStdError', 'CL_exp', '$C_L$ (downforce)'),
                                            ('Cd', 'CdStdError', 'CD_exp', '$C_D$'))):
        ax[k].errorbar(x, [float(r[cfd]) for r in rows], yerr=[float(r[se] or 0) for r in rows],
                       marker='o', capsize=3, label='CFD')
        ax[k].axhline(float(rows[0][exp]), color='k', ls='--', label='experiment (Zerihan 2001)')
        ax[k].set_xlabel('$a_1$'); ax[k].set_ylabel(lab); ax[k].grid(alpha=0.3); ax[k].legend()
    fig.suptitle(f"Tyrrell 026 wing, quoted incidence {rows[0]['alpha']} deg")
    fig.tight_layout(); fig.savefig(a.out + '.png', dpi=150)
    print(f'written {a.out}.csv, {a.out}.png')


if __name__ == '__main__':
    main()
