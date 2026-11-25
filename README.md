# Flexible F1 rear wing – CFD and structural analysis

Course project, *High Performance Simulation Lab for Mechanical Engineering* (Politecnico di Milano).
Williams FW43B(AR) 2021 rear wing: CFD with OpenFOAM, structural analysis of the composite shell with Ansys
Mechanical APDL, one-way fluid–structure coupling, and a flap-angle study on one circuit, the Circuit de
Barcelona-Catalunya (rigid vs flexible wing).


## Repository layout

| Folder | Content |
|---|---|
| `geometry/` | CFD surface of the isolated rear wing for five flap angles (−4°, −2°, 0°, +2°, +4°) |
| `cfd/` | OpenFOAM cases: `rearWing` (steady RANS) and `rearWingURANS` (unsteady RANS, production runs) |
| `structure/flap_<angle>deg/` | ready-to-run Ansys Mechanical APDL model for each flap angle |
| `fsi/` | load transfer CFD → structure and displacement transfer structure → CFD mesh |
| `results/` | results copied from the cluster, one folder per speed and angle, and the summary tables and plots |
| `docs/` | figures, slides, report |

All commands are run from the repository root.

## Objective

The project studies how the flexibility of a Formula 1 rear wing changes its aerodynamic performance, and whether
a flexible wing leads to a different choice of flap angle than a rigid one.

### Approach

- **One wing.** The Williams FW43B(AR) 2021 rear wing, modelled as a composite shell (carbon/epoxy laminates,
  titanium DRS link). Material data and lay-ups are assumptions based on datasheets and literature
  (see `structure/README.md`).
- **Five flap angles.** −4°, −2°, 0°, +2°, +4° (rigid rotation of the flap about its leading edge; positive = steeper
  flap, more downforce and more drag).
- **One circuit, two speeds.** The Circuit de Barcelona-Catalunya (Spanish Grand Prix 2021): a medium-downforce
  circuit, where the efficiency optimum falls in the middle of the flap-angle range and the standard rear wing
  is representative, with a 1 km main straight and both slow and fast corners. Two reference speeds: a corner
  speed (where downforce matters) and a straight-line speed (where drag matters, DRS closed); see
  [Reference speeds](#reference-speeds).
- **Rigid vs flexible.** For every angle and speed the wing is analysed rigid and deformed under its own
  aerodynamic loads (one-way fluid–structure coupling).

### Reference speeds

The speeds are **assumed values**, rounded from published data on the circuit, not taken from the telemetry of the
FW43B: since October 2026 the Formula 1 data APIs used by FastF1 require an F1TV account, so
`circuits/circuit_speeds.py` (which extracts them from the 2021 qualifying telemetry of a Williams driver) can only
be run with such an account.

| Speed | Value | Basis |
|---|---|---|
| Straight line, DRS closed | 310 km/h (86 m/s) | top speeds above 310 km/h at the end of the main straight; the DRS adds about 10–15 km/h, the CFD wing is DRS-closed |
| Corner | 160 km/h (44 m/s) | typical apex speed of the medium-speed corners (Turns 1, 4, 7–8, 12); the slow corners (Turns 5, 10: 70–110 km/h) and the fast ones (Turn 3 about 240 km/h, Turn 9 about 265 km/h in 2021) lie on either side |

Load ratio straight / corner: (310 / 160)² ≈ 3.8, so on the straight the wing carries about four times the
corner load and deforms about four times as much. 

### Simulation pipeline (one angle, one speed)

The same chain is run for each of the 2 × 5 = 10 combinations of speed and flap angle:

1. **Rigid CFD** – unsteady RANS (OpenFOAM `pimpleFoam`, k-ω SST) on the isolated wing. The forces and the
   wall pressure are time averages over the developed flow (`fieldAverage`: `UMean`, `pMean`), after the
   start-up transient; the same run gives the shedding frequency of the wake.
2. **Load transfer** – the mean wall pressure `pMean` is interpolated onto the structural shell mesh and the
   total force is checked against the CFD (`fsi/map_pressure.py`, difference < 2 %; see
   [Load transfer](#load-transfer-pressure-mapping)).
3. **FEM** – linear static analysis in Ansys Mechanical APDL: displacements of the wing.
4. **Mesh morphing** – the structural displacements are imposed on the wing walls of the rigid CFD mesh
   (`fsi/write_point_displacement.py`, `pointDisplacement`) and the volume mesh is moved once with
   `moveDynamicMesh` (see below).
5. **Deformed CFD** – the flow is solved again on the morphed, fixed mesh: downforce and drag of the flexible
   wing.

The 10 chains are independent of each other: on the cluster each one runs as its own sequence of jobs
(rigid CFD → FEM → mesh morphing → deformed CFD), and all of them run in parallel.

| | Runs |
|---|---|
| Rigid CFD | 10 (5 angles × 2 speeds) |
| FEM | 10 |
| Deformed CFD | 10 |
| **Total** | **20 CFD + 10 FEM** |

The rigid wing is run at both speeds even though its force coefficients hardly depend on speed. This keeps
every chain identical and self-contained, the FEM load is the actual load at that speed (no scaling of loads or
displacements with V²), and rigid and flexible wing are compared at the same speed and Reynolds number, so the
difference is due to the deformation only. The two rigid runs of each angle also check that the coefficients
are independent of speed.

### Load transfer (pressure mapping)

The CFD wall mesh and the structural shell mesh describe the same surface but do not match: the CFD walls are
the polygonal boundary faces of the snappyHexMesh mesh (about 5 mm on the wing), the shell is made of quadratic
triangles (SHELL281, about 14 mm). The wall pressure therefore has to be interpolated from one mesh to the other.
`fsi/map_pressure.py` does this in a few steps.

**Inputs**

- CFD: `postProcessing/wingSurfaces/<last time>/walls.vtk`, written by the `wingSurfaces` function object. For
  every wall face $j$ it gives the face centre $\mathbf{x}_j$, the area vector (unit normal $\mathbf{n}_j$, area $A_j$) and the time-averaged
  kinematic pressure `pMean` (m²/s², pressure divided by density, as in every incompressible OpenFOAM solver).
- Structure: `structure/flap_<angle>deg/ansys/rw_shell_<level>_apdl.inp` (`N` and `EN` commands). For every shell
  element $e$ the script computes the centroid $\mathbf{c}_e$, the unit outward normal $\mathbf{n}_e$ and the area $A_e$ from its three
  corner nodes.

**Algorithm**

1. **Gauge pressure.** $p_j = \rho\,(\bar{p}_j - p_\mathrm{ref})$, where $\bar{p}_j$ is `pMean`, $\rho = 1.225\ \mathrm{kg/m^3}$
   and $p_\mathrm{ref} = 0$ (`--rho`, `--pref`).
2. **Orientation.** OpenFOAM boundary normals point out of the fluid, i.e. into the wing. The script compares them
   with the shell normals (median of $\mathbf{n}_j \cdot \mathbf{n}_e$ over nearest pairs) and flips them so that both point out of
   the wing.
3. **Neighbour search (k-d tree).** The CFD face centres are stored in a k-d tree (`scipy.spatial.cKDTree`), a
   binary tree that splits space alternately along x, y and z. For each element centroid the $k = 16$ nearest CFD
   faces are found in $O(\log N)$ operations instead of $O(N)$ for a brute-force search over all $N$ faces (`--k`).
4. **Same-side filter.** Only the neighbours whose normal agrees with the element normal are kept:
   $\mathbf{n}_j \cdot \mathbf{n}_e > 0.5$ (angle below 60°). Near the thin trailing edges of mainplane and flap the closest CFD faces
   can belong to the opposite surface; without this filter the suction and pressure sides would be averaged
   together and the load at the edges would be wrong. If no neighbour passes the filter (sharp edges), the nearest
   face is used; the number of such elements is printed.
5. **Inverse-distance weighting (Shepard's method).** The element pressure is a weighted mean of the selected
   faces, with weights equal to the inverse of the distance $d_j = \lvert \mathbf{x}_j - \mathbf{c}_e \rvert$:

   $$p_e = \frac{\sum_{j} w_j\, p_j}{\sum_{j} w_j}, \qquad w_j = \frac{1}{d_j}$$

   This is the local form of Shepard's inverse-distance interpolation (power 1, only the $k$ nearest points).
   Closer faces count more, the result always lies between the neighbouring CFD values (no overshoot), and a
   uniform pressure is reproduced exactly.
6. **Ansys load.** One command `SFE,e,2,PRES,,`$p_e$ per element: uniform pressure on face 2 (the outward face) of
   SHELL281. A positive pressure on face 2 pushes along $-\mathbf{n}_e$, so a suction ($p_e < 0$) pulls the skin outward.
   The file is read with `/INPUT` by `run_wing_apdl.inp`.
7. **Force check.** The resultant pressure force of the CFD is compared with the force applied to the shell:

   $$\mathbf{F}_\mathrm{CFD} = -\sum_j p_j\, \mathbf{n}_j\, A_j, \qquad \mathbf{F}_\mathrm{FEM} = -\sum_e p_e\, \mathbf{n}_e\, A_e$$

   CFD walls that the shell model does not have (the
   mount block, `rearWing_mount`) are detected as the faces farther than 20 mm from the shell (`--far`) and are
   left out of $\mathbf{F}_\mathrm{CFD}$; their force is printed separately. The mapping is accepted if

   $$\frac{\lvert \mathbf{F}_\mathrm{FEM} - \mathbf{F}_\mathrm{CFD} \rvert}{\lvert \mathbf{F}_\mathrm{CFD} \rvert} < 2\ \%$$

   In Ansys the total reaction of the aerodynamic load step must then be $-\mathbf{F}_\mathrm{FEM}$.

**Type of mapping.** In the terms of the fluid–structure coupling literature this is a *consistent* mapping: an
intensive quantity (the pressure) is interpolated, and the total force is not imposed but checked afterwards. A
*conservative* mapping would instead transfer nodal forces so that their sum is preserved exactly. The consistent
approach is used here because the structural model is loaded with pressures (`SFE`) and because the CFD mesh is
finer than the structural one, so each element pressure is an average of several CFD faces and the force error
stays small (0.9 % at −2°, 86 m/s, mesh `medium`).

**Limits.** The force is conserved only approximately (the error is reduced by refining the CFD or the shell
mesh); the pressure is uniform over each element; the viscous wall shear stress is not transferred (it
contributes to the drag but is small compared with the pressure loads that bend the wing).

**References**

- D. Shepard, *A two-dimensional interpolation function for irregularly-spaced data*, Proc. 23rd ACM National
  Conference, 1968, pp. 517–524. [doi:10.1145/800186.810616](https://doi.org/10.1145/800186.810616)
- J. L. Bentley, *Multidimensional binary search trees used for associative searching*, Communications of the ACM
  18(9), 1975, pp. 509–517. [doi:10.1145/361002.361007](https://doi.org/10.1145/361002.361007)
- A. de Boer, A. H. van Zuijlen, H. Bijl, *Review of coupling methods for non-matching meshes*, Computer Methods
  in Applied Mechanics and Engineering 196(8), 2007, pp. 1515–1525.
- A. de Boer, A. H. van Zuijlen, H. Bijl, *Comparison of conservative and consistent approaches for the coupling
  of non-matching meshes*, Computer Methods in Applied Mechanics and Engineering 197(49), 2008, pp. 4284–4297.
  [doi:10.1016/j.cma.2008.05.001](https://doi.org/10.1016/j.cma.2008.05.001)
- C. Farhat, M. Lesoinne, P. Le Tallec, *Load and motion transfer algorithms for fluid/structure interaction
  problems with non-matching discrete interfaces*, Computer Methods in Applied Mechanics and Engineering 157,
  1998, pp. 95–114. [doi:10.1016/S0045-7825(97)00216-8](https://doi.org/10.1016/S0045-7825(97)00216-8)
- preCICE documentation, *Mapping configuration* (consistent and conservative mappings):
  <https://precice.org/configuration-mapping>

### Mesh morphing

The coupling is one-way, so the CFD mesh is deformed only once, before the deformed CFD, and that run uses a
fixed mesh (no dynamic mesh during the flow solution, no extra cost per time step). The OpenFOAM utility
`moveDynamicMesh` moves the mesh without solving the flow:

- motion solver `displacementLaplacian`: a Laplace equation ∇·(γ ∇u) = 0 for the mesh displacement, with the
  FEM displacements on the wing walls and zero displacement on the outer boundaries;
- diffusivity γ `quadratic inverseDistance` from the wing walls: the boundary-layer cells move almost rigidly
  with the wall and keep their thickness and aspect ratio, the deformation is absorbed by the larger cells away
  from the wing. With a uniform diffusivity the thin prism layers (first layer about 1 mm, wing displacements of
  a few millimetres) would be compressed or stretched;
- one motion step (a few steps for large displacements), then `checkMesh` on the morphed mesh: no inverted
  cells, non-orthogonality and skewness not worse than in the rigid mesh.

The deformed CFD then starts from the same initial field as the rigid one (`potentialFoam`), so that the two runs
differ only in the geometry.

This is the same mesh-motion approach as the previous year's project, where `displacementLaplacian` with a
quadratic inverse-distance diffusivity moved the fluid mesh during a two-way preCICE coupling (there, radial
basis functions were used by preCICE to map forces and displacements between the non-matching fluid and solid
meshes, not to move the mesh). An RBF mesh-motion solver is not among the standard motion solvers of OpenFOAM
v2306 and is not needed for displacements this small.

### Flap-angle selection

The results are copied from the cluster to the PC, one folder per speed and angle:

```
results/
  U<corner speed>/flap<angle>/rigid/        report, summary.csv, history of the coefficients
  U<corner speed>/flap<angle>/deformed/
  U<straight speed>/flap<angle>/rigid/
  U<straight speed>/flap<angle>/deformed/
```

A script in the repository root reads the 20 `summary.csv` files and computes, for each angle and case
(rigid / flexible):

- downforce **L** at the corner speed,
- drag **D** at the straight-line speed,
- efficiency **E = L / D**,

the best flap angle of the rigid and of the flexible wing, and the differences flexible − rigid in L, D and E
(tables and plots in `results/`).

L and D both increase with the flap angle, so "maximum downforce" or "minimum drag" alone would always select the
steepest or the flattest flap. The efficiency, instead, has a maximum at an intermediate angle and is the reference
for the choice. The flexible wing bends more on the straight (load ∝ V²), which lowers its drag there: its
efficiency curve, and so the best flap angle, can shift with respect to the rigid wing. The main results are this
shift and the gain or loss in downforce and drag.

### Numerical and HPC studies

- **Mesh independence** at 0°: three CFD meshes and three structural meshes (coarse / medium / fine). The chosen
  meshes are then used for all angles.
- **Time step and averaging window** of the URANS runs: Courant number, number of time steps per shedding period,
  and the end of the start-up transient (MSER check) before the averaging window.
- **Scaling** of the CFD on the cluster: strong scaling (fixed mesh, increasing processes) and weak scaling
  (fixed cells per process) to choose the number of processes for the production runs.

### Assumptions and limits

- Isolated rear wing (no car upstream), DRS closed.
- One circuit, with two assumed reference speeds (rounded published values, not team telemetry); other circuits
  would need their own speeds. The corner speed is a single representative value: the circuit has corners from
  about 70 to 265 km/h.
- Laminates and material data are assumed, not team data; the shell model is the outer skin only.
- One-way coupling: the deformed-wing loads are not fed back to the structure; the structure is loaded with the
  time-averaged pressure (static analysis, no dynamic response to the wake shedding).

## Workflow (one flap angle, one speed)

1. **CFD, rigid wing** – OpenFOAM `pimpleFoam`, k-ω SST, on `geometry/FW43B_rearWing_flap_<angle>deg.stl`, at the
   chosen speed (`cfd/rearWingURANS`: `./Allrun.pre -angle <angle>`, then `./Allrun -speed <U>`); the
   `wingSurfaces` function object writes the wall pressure, instantaneous and time-averaged (VTK).
2. **Loads** – `fsi/map_pressure.py` maps the mean pressure onto the shell mesh of the same angle:
   ```
   python3 fsi/map_pressure.py structure/flap_0deg/ansys/rw_shell_medium_apdl.inp \
           results/U<U>/flap0/rigid/postProcessing/wingSurfaces/<last time>/walls.vtk --field pMean \
           --apdl structure/flap_0deg/ansys/loads_pressure_apdl.inp
   ```
   It prints the CFD force and the applied force (difference must be < 2 %).
3. **Structure** – Mechanical APDL, `structure/flap_<angle>deg/ansys/run_wing_apdl.inp` (see `structure/README.md`).
   Checks: self weight; total reaction of the aerodynamic load step = − applied force.
   The run writes `wing_disp.csv` (nodal displacements).
4. **Mesh morphing** – the structural displacements (`wing_disp.csv`) are written on the wing walls of a copy of
   the rigid case, and the mesh is moved:
   ```
   python3 fsi/write_point_displacement.py <deformed case> --table wing_disp.csv
   moveDynamicMesh -parallel      # displacementLaplacian, quadratic inverseDistance (constant/dynamicMeshDict)
   checkMesh -parallel
   ```
5. **Deformed CFD** – a second `pimpleFoam` run at the same speed, on the morphed mesh.

## Requirements

- WSL Ubuntu, OpenFOAM v2306, Python 3 with `numpy scipy matplotlib` (`fastf1` for `circuits/circuit_speeds.py`)
- ParaView (`pvbatch`) and `pdflatex` for the reports of the CFD runs
- Ansys Mechanical APDL (the Student version is enough: ~128k node limit)


## Team

Marchettini Bianca Maria, Perlini Aurora, Trombetta Davide (HPC)
