#!/usr/bin/env python3
# readSets.py -- one sampled line of the function "profiles" (type sets,
#                format raw) as a single table, whatever the OpenFOAM version
#
# Usage, from the case directory:
#     python3 tools/readSets.py <folder> <set>
#     e.g. python3 tools/readSets.py postProcessing/profiles/300 wake1
#
# Prints one line per point:  coordinate  <field columns>, with a header line
#     # coord Cp0 k p Ux Uy Uz
# (scalars as they are, vectors as their x, y, z components).
#
# Why: the raw writer names its files after the set and the fields, and the
# versions differ in how they split them:
#   OpenFOAM.org          wake1.xy                    (all fields, fixed order)
#   ESI v2306 / v2412     wake1_Cp0_k_p_U.xy          (fields sorted by name)
#   older ESI             wake1_Cp0_k_p.xy, wake1_U.xy (scalars and vectors)
# Every file <set>_*.xy of the folder is read; the field names come from the
# file name, the number of components of each field from the number of
# columns (a vector has 3; U, the only vector sampled here, is assumed to be
# the vector when the columns leave the choice open).

import glob
import os
import sys

VECTORS = ('U', 'wallShearStress', 'vorticity')


def read(path, fields):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    if not rows:
        return [], []
    ncol = len(rows[0]) - 1                   # first column: the coordinate
    extra = ncol - len(fields)                # 2 per vector field
    names = []
    for f in fields:
        if extra >= 2 and f in VECTORS:
            names += [f + 'x', f + 'y', f + 'z']
            extra -= 2
        else:
            names.append(f)
    if len(names) != ncol:
        sys.exit('%s: %d columns for fields %s' % (path, ncol, fields))
    return names, rows


def main():
    folder, name = sys.argv[1], sys.argv[2]
    files = sorted(glob.glob(os.path.join(folder, name + '_*.xy')))
    if not files and os.path.exists(os.path.join(folder, name + '.xy')):
        files = [os.path.join(folder, name + '.xy')]
    if not files:
        sys.exit('no %s*.xy in %s' % (name, folder))
    coord, table, header = None, None, []
    for path in files:
        base = os.path.basename(path)[:-3]
        fields = base[len(name) + 1:].split('_') if base != name else \
            ['U', 'p', 'k', 'Cp0']                # OpenFOAM.org order
        names, rows = read(path, fields)
        if coord is None:
            coord = [r[0] for r in rows]
            table = [[] for _ in rows]
        for t, r in zip(table, rows):
            t += r[1:]
        header += names
    print('# coord ' + ' '.join(header))
    for c, t in zip(coord, table):
        print(c + ' ' + ' '.join(t))


if __name__ == '__main__':
    main()
