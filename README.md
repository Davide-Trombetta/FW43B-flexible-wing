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
