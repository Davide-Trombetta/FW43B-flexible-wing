# Flexible F1 rear wing – CFD and structural analysis

Course project, *High Performance Simulation Lab for Mechanical Engineering* (Politecnico di Milano).
Williams FW43B(AR) 2021 rear wing: CFD with OpenFOAM, structural analysis of the composite shell with Ansys
Mechanical APDL, one-way fluid–structure coupling, and a flap-angle study on three circuits (rigid vs flexible wing).


## Repository layout

| Folder | Content |
|---|---|
| `geometry/` | CFD surface of the isolated rear wing for five flap angles (−4°, −2°, 0°, +2°, +4°) |
| `cfd/` | The full CFD case with OpenFOAM goes here |
| `structure/flap_<angle>deg/` | ready-to-run Ansys Mechanical APDL model for each flap angle |
| `fsi/` | load transfer CFD → structure and displacement transfer structure → CFD mesh |
| `results/` | summary results only (CSV tables, plots) |
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
- **Three circuits.** Low, medium and high downforce, e.g. Monza, Barcelona and the Hungaroring. Two reference
  speeds per circuit from the 2021 qualifying telemetry: a corner speed and a straight-line speed
  (DRS closed).
- **Rigid vs flexible.** For every angle and speed the wing is analysed rigid and deformed under its own
  aerodynamic loads (one-way fluid–structure coupling).

### Simulation chain (one angle, one speed)

1. **Rigid CFD** – steady RANS (OpenFOAM `simpleFoam`, k-ω SST) on the isolated wing: downforce, drag and the
   wall-pressure field.
2. **Load transfer** – the pressure is mapped onto the structural shell mesh with force conservation
   (`fsi/map_pressure.py`, difference < 2 %).
3. **FEM** – linear static analysis in Ansys Mechanical APDL: displacements of the wing.
4. **Deformed CFD** – the CFD mesh is morphed with the structural displacements (`fsi/write_point_displacement.py`)
   and the flow is solved again: downforce and drag of the flexible wing.

### Flap-angle selection

For each angle and case (rigid / flexible) we compute:

- downforce **L** at the corner speed,
- drag **D** at the straight-line speed,
- aerodynamic efficiency **E = L / D**.

L and D both increase with the flap angle, so "maximum downforce" or "minimum drag" alone would always select the
steepest or the flattest flap. The efficiency, instead, has a maximum at an intermediate angle and is the reference
for the choice:

| Circuit type | Priority | Chosen angle |
|---|---|---|
| Low downforce (e.g. Monza) | low drag on long straights | below the efficiency maximum |
| Medium downforce (e.g. Barcelona) | balance | close to the efficiency maximum |
| High downforce (e.g. Hungaroring) | grip in slow corners | above the efficiency maximum |

The flexible wing bends more on the straights (load ∝ V²), which lowers its drag there: its efficiency curve, and so
the best angle for each circuit, can shift with respect to the rigid wing. The main results are this shift and
the gain or loss in downforce and drag for each circuit.

### Reducing the number of simulations

With 5 angles and 3 circuits there are 6 speeds per angle (3 circuits × corner / straight). Running every case
would need 60 CFD and 30 FEM runs. Two physical properties reduce this:

- **Rigid CFD: one per angle.** The rigid geometry does not change with speed and the force coefficients are
  practically independent of it (incompressible flow, Mach < 0.3), so F = ½ ρ V² S C gives the forces at every
  speed. This is checked once, running the same angle at two speeds.
- **FEM: one per angle.** The structural model is linear and the loads scale with V², so the displacements do too:
  u(V) = u(V_ref) · (V / V_ref)² (`--scale` option of `fsi/write_point_displacement.py`).
- **Deformed CFD: one per angle and speed.** The deformed shape changes with speed, so these runs are needed.

| | Every case | With scaling |
|---|---|---|
| Rigid CFD | 30 | 5 (+1 check) |
| FEM | 30 | 5 |
| Deformed CFD | 30 | 30 |
| **Total** | **60 CFD + 30 FEM** | **36 CFD + 5 FEM** |

### Numerical and HPC studies

- **Mesh independence** at 0°: three CFD meshes and three structural meshes (coarse / medium / fine). The chosen meshes are then used for all angles.
- **Scaling** of the CFD on the cluster: strong scaling (fixed mesh, increasing processes) and weak scaling
  (fixed cells per process) to choose the number of processes for the production runs.

### Assumptions and limits

- Isolated rear wing (no car upstream), DRS closed.
- Laminates and material data are assumed, not team data; the shell model is the outer skin only.
- One-way coupling: the deformed-wing loads are not fed back to the structure.

## Workflow (one flap angle)

1. **CFD, rigid wing** – OpenFOAM `simpleFoam`, k-ω SST, on `geometry/FW43B_rearWing_flap_<angle>deg.stl`;
   the `wingSurfaces` function object writes the wall pressure (VTK).
2. **Loads** – `fsi/map_pressure.py` maps the pressure onto the shell mesh of the same angle:
   ```
   python3 fsi/map_pressure.py structure/flap_0deg/ansys/rw_shell_medium_apdl.inp \
           postProcessing/wingSurfaces/<time>/*.vtk -o structure/flap_0deg/ansys/loads_pressure_apdl.inp
   ```
   It prints the CFD force and the applied force (difference must be < 2 %).
3. **Structure** – Mechanical APDL, `structure/flap_<angle>deg/ansys/run_wing_apdl.inp` (see `structure/README.md`).
   Checks: self weight; total reaction of the aerodynamic load step = − applied force.
   The run writes `wing_disp.csv` (nodal displacements).
4. **Deformed CFD** – the structural displacements (`wing_disp.csv`) deform the existing CFD mesh:
   ```
   python3 fsi/write_point_displacement.py <deformed case> --table wing_disp.csv [--scale (V/Vref)^2]
   ```
   then `checkMesh` and a second `simpleFoam` run. The structure is linear, so one structural run per angle gives
   the deformation at any speed through `--scale (V/V_ref)^2`.

## Requirements

- WSL Ubuntu, OpenFOAM v2306, Python 3 with `numpy scipy matplotlib` (`fastf1` for `circuits/circuit_speeds.py`)
- Ansys Mechanical APDL (the Student version is enough: ~128k node limit)


## Team

Marchettini Bianca Maria, Perlini Aurora, Trombetta Davide (HPC)
