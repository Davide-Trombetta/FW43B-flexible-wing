#!/usr/bin/env python3
"""
map_pressure.py - one-way FSI load transfer: OpenFOAM wall pressure -> Abaqus (*DLOAD) or Ansys Mechanical (Imported Pressure).

Inputs
  cfd    one or more OpenFOAM surface files in legacy VTK (ASCII, POLYDATA) with a CELL_DATA field 'p'
         (kinematic pressure p/rho for incompressible solvers), written by the 'surfaces' function object
         on the rear-wing patches - see the pipeline document for the controlDict entry.
  mesh   (recommended) the shell mesh used in Ansys (rw_shell_*.inp, imported with External Model).

Abaqus (--abaqus FILE): *DLOAD data lines  elem, P, value  for every shell element (STRI65, S3R or S6
mesh), value = -rho*(p - pref): for Abaqus shell elements "the pressure is positive in the direction of the positive
element normal", and the normals of our meshes point out of the wing. The file is included under *DLOAD in run_wing_abaqus.inp.
Same-side inverse-distance mapping at the element centroids; the force check is printed.

Mechanical APDL (--apdl FILE): SFE commands on face 2 (top = outward face) of the SHELL281 elements, value
rho*(p - pref) (positive pressure on face 2 acts along -N). Read with /INPUT in run_wing_apdl.inp.

Ansys Workbench output: CSV  X [m], Y [m], Z [m], Pressure [Pa]   (gauge pressure rho*(p - pref), positive = pushing on the
surface, the Ansys convention) for Mechanical -> External Data -> Imported Pressure.

With the shell mesh (recommended)
  the pressure is mapped by this script onto the NODES of the Ansys shell mesh: for every node the k nearest
  CFD faces ON THE SAME SIDE of the skin (agreement with the node normal) are inverse-distance weighted. The CSV points coincide with the Ansys nodes, so the Ansys mapping
  is exact and the pressures of the upper and lower surfaces are never mixed at the thin trailing edges
  (which can happen if Ansys interpolates the raw CFD point cloud).
  The script prints the force applied to the shell next to the CFD force (target: difference < 1-2 %).
Without the mesh
  the CSV holds the CFD face centres and values; Ansys does the interpolation (check the thin edges).

Usage
  python3 map_pressure.py rw_shell_medium_STRI65.inp postProcessing/wingSurfaces/<time>/*.vtk --abaqus loads_pressure_abaqus.inp
  python3 map_pressure.py rw_shell_medium.inp postProcessing/wingSurfaces/<time>/*.vtk -o wing_pressure.csv
  python3 map_pressure.py rw_shell_medium_apdl.inp postProcessing/wingSurfaces/<time>/walls.vtk --field pMean --apdl loads_pressure_apdl.inp
  (the mesh can be the Abaqus-style .inp or the APDL file written by to_apdl.py)
  python3 map_pressure.py postProcessing/wingSurfaces/<time>/*.vtk -o wing_pressure_raw.csv
"""
import argparse, re
import numpy as np
from scipy.spatial import cKDTree

# ----------------------------------------------------------------------------------------------
def read_apdl(txt):
    """Mechanical APDL mesh written by to_apdl.py (rw_shell_*_apdl.inp): N,id,x,y,z and
    EN,e,I,J,K,K,M,N,K,P (SHELL281 as degenerate 6-node triangle), regions from the comment lines
    '! --- region NAME: ...' -> element sets E_NAME with rows (e, I, J, K, M, N, P) as in the .inp."""
    nodes = {}; elems = {}; reg = "E_ALL"
    for l in txt.splitlines():
        if l.startswith("N,"):
            a = l.split(","); nodes[int(a[1])] = [float(v) for v in a[2:5]]
        elif l.startswith("EN,"):
            a = [int(v) for v in l.split(",")[1:10]]
            elems.setdefault(reg, []).append([a[0], a[1], a[2], a[3], a[5], a[6], a[8]])
        else:
            m = re.match(r"!\s*---\s*region\s+(\w+)", l, re.I)
            if m: reg = "E_" + m.group(1).upper()
    return nodes, {r: np.array(v) for r, v in elems.items()}

def read_inp(path):
    """shell mesh: Abaqus-style .inp (*NODE, *ELEMENT TYPE=S6 ...) or APDL (rw_shell_*_apdl.inp)"""
    txt = open(path).read()
    if re.search(r"^EN,", txt, re.M):
        nodes, elems = read_apdl(txt)
        ids = np.array(sorted(nodes)); X = np.zeros((ids.max() + 1, 3)); X[ids] = [nodes[i] for i in ids]
        return X, elems
    nodes = {}
    for blk in re.finditer(r"\*NODE[^\n]*\n(.*?)(?=\n\*)", txt, re.S | re.I):
        for l in blk.group(1).splitlines():
            if l.strip():
                a = l.split(","); nodes[int(a[0])] = [float(v) for v in a[1:4]]
    elems = {}
    for blk in re.finditer(r"\*ELEMENT,\s*TYPE=(?:S6|STRI65|S3R|S3),\s*ELSET=(\w+)\n(.*?)(?=\n\*|\Z)", txt, re.S | re.I):
        rows = [list(map(int, l.split(","))) for l in blk.group(2).splitlines() if l.strip()]
        elems[blk.group(1).upper()] = np.array(rows)
    ids = np.array(sorted(nodes)); X = np.zeros((ids.max() + 1, 3)); X[ids] = [nodes[i] for i in ids]
    return X, elems

def tri_geom(X, conn):
    a, b, c = X[conn[:, 1]], X[conn[:, 2]], X[conn[:, 3]]
    cr = np.cross(b - a, c - a); A = 0.5 * np.linalg.norm(cr, axis=1)
    return (a + b + c) / 3.0, cr / (2 * A[:, None]), A

def read_vtk_polydata(path, field="p"):
    """Minimal reader for legacy ASCII VTK POLYDATA (OpenFOAM 'vtk' surface writer, format ascii)."""
    tok = open(path).read().split()
    i = 0; pts = None; polys = []; vals = None
    while i < len(tok):
        t = tok[i].upper()
        if t == "POINTS":
            n = int(tok[i + 1]); pts = np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3); i += 3 + 3 * n; continue
        if t == "POLYGONS":
            n = int(tok[i + 1]); size = int(tok[i + 2]); arr = list(map(int, tok[i + 3:i + 3 + size])); i += 3 + size
            j = 0
            while j < len(arr):
                k = arr[j]; polys.append(arr[j + 1:j + 1 + k]); j += 1 + k
            continue
        if t in ("SCALARS",) and tok[i + 1] == field:
            # SCALARS p float [1] \n LOOKUP_TABLE default \n values...
            j = i + 3
            if tok[j].isdigit(): j += 1
            if tok[j].upper() == "LOOKUP_TABLE": j += 2
            vals = np.array(tok[j:j + len(polys)], float); i = j + len(polys); continue
        if t == "FIELD":                                          # FIELD attributes n \n name ncomp ntuples type
            nf = int(tok[i + 2]); j = i + 3
            for _ in range(nf):
                name, nc, nt = tok[j], int(tok[j + 1]), int(tok[j + 2]); j += 4
                data = np.array(tok[j:j + nc * nt], float); j += nc * nt
                if name == field and nc == 1: vals = data
            i = j; continue
        i += 1
    if pts is None or vals is None: raise ValueError(f"{path}: POINTS/POLYGONS/{field} not found")
    C = np.zeros((len(polys), 3)); An = np.zeros((len(polys), 3))
    for k, poly in enumerate(polys):                              # fan triangulation for polygons
        P = pts[poly]; C[k] = P.mean(0)
        An[k] = 0.5 * np.cross(P[1:-1] - P[0], P[2:] - P[0]).sum(0)
    A = np.linalg.norm(An, axis=1)
    return C, An / A[:, None], A, vals

# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="[shell mesh .inp] + OpenFOAM surface files (legacy VTK, ascii)")
    ap.add_argument("--field", default="p"); ap.add_argument("--rho", type=float, default=1.225)
    ap.add_argument("--pref", type=float, default=0.0, help="reference kinematic pressure subtracted")
    ap.add_argument("--k", type=int, default=16); ap.add_argument("-o", "--out", default="wing_pressure.csv")
    ap.add_argument("--far", type=float, default=0.02,
                    help="CFD faces farther than this from the shell [m] (walls not modelled) are left out of the check")
    ap.add_argument("--abaqus", metavar="FILE", help="write Abaqus *DLOAD data lines (needs the shell mesh)")
    ap.add_argument("--apdl", metavar="FILE", help="write Mechanical APDL SFE commands (needs the shell mesh)")
    ap.add_argument("--normals-outward", action="store_true",
                    help="(no mesh given) CFD normals already point out of the body; default: OpenFOAM convention")
    a = ap.parse_args()
    meshes = [f for f in a.inputs if f.lower().endswith(".inp")]
    cfd = [f for f in a.inputs if not f.lower().endswith(".inp")]
    if len(meshes) > 1 or not cfd: ap.error("give at most one .inp mesh and at least one CFD surface file")

    # CFD side
    parts = [read_vtk_polydata(f, a.field) for f in cfd]
    Cc = np.vstack([p[0] for p in parts]); Nc = np.vstack([p[1] for p in parts])
    Ac = np.concatenate([p[2] for p in parts]); pc = np.concatenate([p[3] for p in parts])
    pe_c = a.rho * (pc - a.pref)
    print(f"CFD: {len(Cc)} faces from {len(cfd)} file(s), wetted area {Ac.sum():.4f} m2, "
          f"p from {pe_c.min():.1f} to {pe_c.max():.1f} Pa")

    # orientation: OpenFOAM boundary normals point out of the FLUID (into the body)
    if meshes:
        X, elems = read_inp(meshes[0])
        regs = [r for r in elems if r.startswith("E_") and r != "E_LOWER_TEST"]
        geo = {r: tri_geom(X, elems[r]) for r in regs}
        Ce = np.vstack([geo[r][0] for r in regs]); Ne = np.vstack([geo[r][1] for r in regs])
        _, j1 = cKDTree(Cc).query(Ce, k=1)
        if np.median((Nc[j1] * Ne).sum(1)) < 0:
            Nc = -Nc; print("CFD normals point into the body (OpenFOAM convention) -> flipped to outward")
    elif not a.normals_outward:
        Nc = -Nc
    F_cfd = -(pe_c[:, None] * Nc * Ac[:, None]).sum(0)           # pressure force on the body
    print(f"CFD pressure force on the wing   Fx={F_cfd[0]:9.2f} Fy={F_cfd[1]:9.2f} Fz={F_cfd[2]:9.2f} N")
    if meshes:
        # CFD walls that the shell model does not have (e.g. the mount block, rearWing_mount):
        # faces farther than --far from every shell node/centroid; the check uses the rest
        used = np.unique(np.concatenate([elems[r][:, 1:].ravel() for r in regs]))
        dfar, _ = cKDTree(np.vstack([X[used], Ce])).query(Cc)
        far = dfar > a.far
        if far.any():
            F_far = -(pe_c[far, None] * Nc[far] * Ac[far, None]).sum(0)
            F_cfd = F_cfd - F_far
            print(f"  of which on walls not in the shell model ({far.sum()} faces, {Ac[far].sum():.4f} m2, "
                  f"more than {1e3*a.far:g} mm from it): Fx={F_far[0]:.2f} Fy={F_far[1]:.2f} Fz={F_far[2]:.2f} N")
            print(f"CFD force on the walls of the shell model  Fx={F_cfd[0]:9.2f} Fy={F_cfd[1]:9.2f} Fz={F_cfd[2]:9.2f} N")

    if a.abaqus or a.apdl:
        if not meshes: ap.error("--abaqus / --apdl need the shell mesh (.inp)")
        tree = cKDTree(Cc); kk = min(a.k, len(Cc)); lines = []; Ftot = np.zeros(3); nfb = 0
        for r in regs:
            C, N, A = geo[r]; conn = elems[r]
            d, j = tree.query(C, k=kk)
            same = (Nc[j] * N[:, None, :]).sum(2) > 0.5
            w = np.where(same, 1.0 / np.maximum(d, 1e-9), 0.0)
            none = w.sum(1) == 0; w[none, 0] = 1.0; nfb += none.sum()
            w /= w.sum(1, keepdims=True)
            pe = (w * pe_c[j]).sum(1)
            Ftot += -(pe[:, None] * N * A[:, None]).sum(0)
            lines += [(e, v) for e, v in zip(conn[:, 0], pe)]
        if a.abaqus:
            with open(a.abaqus, "w") as o:
                o.write("** aerodynamic pressure on the shell elements, generated by map_pressure.py\n")
                o.write("** data lines for *DLOAD. Abaqus shells: positive P acts along the positive (outward) element normal,\n")
                o.write("** so P = -rho*(p - pref): a suction (p < pref) pulls the skin outward.\n")
                o.write("\n".join(f"{e}, P, {-v:.6g}" for e, v in lines) + "\n")
            print(f"written {a.abaqus}: {len(lines)} element loads (Abaqus)")
        if a.apdl:
            with open(a.apdl, "w") as o:
                o.write("! aerodynamic pressure on the SHELL281 elements, generated by map_pressure.py\n")
                o.write("! face 2 = top face (normal out of the wing); positive pressure on face 2 acts along -N,\n")
                o.write("! so the value is rho*(p - pref): a suction (p < pref) is negative and pulls the skin outward.\n")
                o.write("\n".join(f"SFE,{e},2,PRES,,{v:.6g}" for e, v in lines) + "\n")
            print(f"written {a.apdl}: {len(lines)} element loads (Mechanical APDL)")
        err = np.linalg.norm(Ftot - F_cfd) / np.linalg.norm(F_cfd)
        print(f"load applied to the shell        Fx={Ftot[0]:9.2f} Fy={Ftot[1]:9.2f} Fz={Ftot[2]:9.2f} N  "
              f"({nfb} elements on sharp edges used the nearest face)")
        print(f"relative difference {100*err:.2f} %" + ("   <-- WARNING > 2 %" if err > 0.02 else "   OK"))
        print(f"check (Abaqus .dat step AERO_PRESSURE / APDL load step 2): total reaction = -({Ftot[0]:.1f}, {Ftot[1]:.1f}, {Ftot[2]:.1f}) N")
        return

    if not meshes:
        with open(a.out, "w") as o:
            o.write("X [m],Y [m],Z [m],Pressure [Pa]\n")
            o.write("".join(f"{x:.7g},{y:.7g},{z:.7g},{v:.6g}\n" for (x, y, z), v in zip(Cc, pe_c)))
        print(f"written {a.out}: {len(Cc)} CFD points (Ansys interpolates; check the thin edges)")
        print(f"check in Mechanical: sum of the reactions = -({F_cfd[0]:.1f}, {F_cfd[1]:.1f}, {F_cfd[2]:.1f}) N")
        return

    # nodal values: k nearest CFD faces on the same side of the skin, inverse-distance weighted
    Nn = np.zeros((len(X), 3)); used = np.zeros(len(X), bool)
    for r in regs:
        C, N, A = geo[r]; conn = elems[r]
        for c in range(1, conn.shape[1]):
            np.add.at(Nn, conn[:, c], N * A[:, None]); used[conn[:, c]] = True
    ids = np.nonzero(used)[0]
    nn = Nn[ids] / np.maximum(np.linalg.norm(Nn[ids], axis=1, keepdims=True), 1e-30)
    tree = cKDTree(Cc); kk = min(a.k, len(Cc))
    d, j = tree.query(X[ids], k=kk)
    same = (Nc[j] * nn[:, None, :]).sum(2) > 0.5
    w = np.where(same, 1.0 / np.maximum(d, 1e-9), 0.0)
    none = w.sum(1) == 0; w[none, 0] = 1.0                          # sharp edges: nearest face
    w /= w.sum(1, keepdims=True)
    pn = np.zeros(len(X)); pn[ids] = (w * pe_c[j]).sum(1)
    print(f"shell nodes {len(ids)}, nodes on sharp edges (nearest CFD face used) {none.sum()}")
    # force check with the nodal values (what Ansys will apply): quadratic triangle, mid-side nodes weight 1/3
    Ftot = np.zeros(3)
    for r in regs:
        C, N, A = geo[r]; conn = elems[r]
        pel = pn[conn[:, 4:7]].mean(1) if conn.shape[1] == 7 else pn[conn[:, 1:4]].mean(1)   # exact for quadratic / linear
        Ftot += -(pel[:, None] * N * A[:, None]).sum(0)
    with open(a.out, "w") as o:
        o.write("X [m],Y [m],Z [m],Pressure [Pa]\n")
        o.write("".join("%.9g,%.9g,%.9g,%.6g\n" % (*X[i], pn[i]) for i in ids))
    err = np.linalg.norm(Ftot - F_cfd) / np.linalg.norm(F_cfd)
    print(f"\nCFD force, walls of the shell     Fx={F_cfd[0]:9.2f} Fy={F_cfd[1]:9.2f} Fz={F_cfd[2]:9.2f} N")
    print(f"load applied to the shell        Fx={Ftot[0]:9.2f} Fy={Ftot[1]:9.2f} Fz={Ftot[2]:9.2f} N")
    print(f"relative difference {100*err:.2f} %" + ("   <-- WARNING > 2 %" if err > 0.02 else "   OK"))
    print(f"written {a.out}: {len(ids)} shell nodes")
    print(f"check in Mechanical: sum of the reactions = -({Ftot[0]:.1f}, {Ftot[1]:.1f}, {Ftot[2]:.1f}) N")

if __name__ == "__main__":
    main()
