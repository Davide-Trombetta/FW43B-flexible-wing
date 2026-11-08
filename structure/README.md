# Structural models (Ansys Mechanical APDL)

One folder per flap angle, `flap_<angle>deg/ansys/` (−4, −2, 0, +2, +4):

| File | Content |
|---|---|
| `run_wing_apdl_selfweight.inp` | first check: self weight only (no CFD loads needed) |
| `run_wing_apdl.inp` | load step 1 self weight, load step 2 aerodynamic pressure; prints the reactions, writes `wing_disp.csv` |
| `materials_laminate_apdl.inp` | materials, element coordinate systems, layered sections |
| `rw_shell_<mesh>_apdl.inp` | mesh: nodes, SHELL281 elements (6-node triangles), components |
| `loads_pressure_apdl.inp` | *not versioned*: written by `fsi/map_pressure.py` from the CFD of the same angle |

Meshes: `medium` (74 314 nodes, 37 180 elements) for every angle; `flap_0deg` also has `coarse` (46 617 nodes)
and `fine` (125 018 nodes) for the mesh-independence study (change the `/INPUT,rw_shell_medium_apdl,inp` line).
Element normals point out of the wing.

## How to run

1. Mechanical APDL Product Launcher: *Working Directory* = `structure/flap_<angle>deg/ansys`, Simulation Environment
   *ANSYS*, *High Performance Computing Setup* → SMP with the available cores, *Run*.
2. *File → Read Input from…* → `run_wing_apdl_selfweight.inp` (first check) or `run_wing_apdl.inp`
   (needs `loads_pressure_apdl.inp`). To start again in the same session: `/CLEAR,START`.
3. Reactions are printed in the output window (line *TOTAL VALUES*); `wing_disp.csv` is written in the working
   directory (displacements of load step 2, `X, Y, Z, UX, UY, UZ` in metres).

**Checks**: self weight → total FZ ≈ 61.5 N (mass 6.27 kg); aerodynamic load → total reaction = − force printed by
`map_pressure.py` (verified at 0° with a test load: reactions equal and opposite to the load to 0.1 N).
Pressure convention: `SFE,elem,2,PRES` on face 2 (outer face), positive pushes inwards, suction is negative.


## Regions

| Component | Part | Section |
|---|---|---|
| `E_MAINPLANE` | main plane | [F/0/45/−45/0]s, 1.438 mm, 0° spanwise (y) |
| `E_FLAP` | flap | [F/0/45/−45/0]s, 1.438 mm, 0° spanwise (y) |
| `E_DRSPOD` | DRS pod fairing | [F/0/45/−45/0]s, 1.438 mm, 0° vertical (z) |
| `E_SWANNECK` | swan-neck arms | [F/0/0/45/−45/0/0]s, 1.962 mm, 0° vertical |
| `E_DRSLINK` | DRS link | titanium, 2 mm shell (hollow approximation of a solid part) |
| `E_ENDPLATE_YPOS`, `E_ENDPLATE_YNEG` | endplates | [F/0/45/−45/0]s, 1.438 mm, 0° vertical |
| `E_PYLONS` | pylon forks (attachment) | [F/0/0/45/−45/0/0]s, 1.962 mm, 0° vertical |

Constraints: translations fixed at the pylon base and at the endplate feet (node components `N_PYLON_BASE`,
`N_ENDPLATE_FOOT_YPOS`, `N_ENDPLATE_FOOT_YNEG`, together `N_SUPPORT`). Laminate mid-surface on the aerodynamic surface.

## Materials (assumptions)

Literature and datasheet values, **not team data**. In-plane properties of the UD ply from the manufacturer's
datasheet [1]; shear moduli, Poisson's ratios and out-of-plane properties from [2]; titanium from [3].

| Property | UD carbon/epoxy, IM7/8552 | Woven carbon/epoxy (5HS), AGP370-5H/3501-6S | Ti-6Al-4V, Grade 5 annealed |
|---|---|---|---|
| E1 | 164 GPa [1] | 77 GPa [2, A.6] | 113.8 GPa [3] (isotropic) |
| E2 | 12 GPa [1] | 75 GPa [2, A.6] | – |
| E3 | 12 GPa (= E2, transverse isotropy) | 13.8 GPa [2, A.6] | – |
| G12 | 7.8 GPa [2, A.4, IM7/977-3] | 6.5 GPa [2, A.6] | 44 GPa [3] (not an input: APDL derives it from E, ν) |
| G13 | 7.8 GPa (= G12) | 5.1 GPa [2, A.6] | – |
| G23 | 3.7 GPa [2, A.6, AS4/3501-6] | 4.1 GPa [2, A.6] | – |
| ν12 | 0.35 [2, A.4, IM7/977-3] | 0.06 [2, A.6] | 0.342 [3] |
| ν13 | 0.35 (= ν12) | 0.50 [2, A.6] | – |
| ν23 | 0.54 [2, A.6, AS4/3501-6] | 0.37 [2, A.6] | – |
| Density | 1570 kg/m³ [1] | 1600 kg/m³ [2, A.6] | 4430 kg/m³ [3] |
| Ply thickness | 0.131 mm [1] | 0.195 mm [1] (AGP193-PW cured ply) | 2 mm shell (DRS link) |

Axes: 1 along the fibres (warp for the fabric), 2 transverse in the ply plane (fill), 3 through the thickness.
The UD values from [2] for an IM7 system with toughened epoxy (977-3) and the out-of-plane values of AS4/3501-6
complete the datasheet, which gives only the in-plane moduli. The transverse shear modulus is consistent with
transverse isotropy: E2 / (2 (1 + ν23)) = 3.9 GPa vs 3.7 GPa.

**References**

1. Hexcel Corporation, *HexPly® 8552 – Epoxy matrix (180 °C curing matrix), Product Data Sheet (EU)*, 2023.
   https://hexcel.com/wp-content/uploads/2025/12/HexPly_8552_eu_DataSheet1.pdf
2. I. M. Daniel, O. Ishai, *Engineering Mechanics of Composite Materials*, 2nd ed., Oxford University Press, New York,
   2006, Appendix A, Tables A.4 and A.6.
3. Aerospace Metals (data from Allvac / ASM), *Titanium Ti-6Al-4V (Grade 5), Annealed – data sheet*.
   https://www.aerospacemetals.com/wp-content/uploads/2023/07/Titanium-Ti-6Al-4V-Grade-5-Annealed.pdf

Fibre reference (intermediate modulus class): Toray, *TORAYCA® T800S Data Sheet* (fibre modulus 294 GPa,
UD composite 163 GPa at 60 % fibre volume), https://www.toraycma.com/wp-content/uploads/T800S-Data-Sheet.pdf

To change a property edit the `MP` lines; to change a laminate edit the `SECDATA` lines (one per ply: thickness,
material, angle, integration points) of the section in `materials_laminate_apdl.inp` – in every angle folder.

The model is the outer skin only (real wings also have internal spars): with the real CFD loads check that the
deflections are realistic and, if needed, add 0° plies.

## Flap angles

The flap of each angle is rotated rigidly about the hinge (see `geometry/README.md`); the endplate, DRS link and pod
nodes near the flap follow smoothly within 40 mm, the main plane, swan necks and pylons are unchanged. Element quality
is the same as at 0°, and the flap nodes lie on the STL surface of the same angle (mean distance 0.08 mm).
