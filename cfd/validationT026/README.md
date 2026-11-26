# Validation case: Tyrrell 026 single-element wing (a1 of k-ω SST)

Before the FW43B study, the constant a1 of k-ω SST is calibrated on a wing with known experimental
data, simulated with **the same workflow and settings as the FW43B case** (`cfd/rearWing`: steady RANS,
`simpleFoam`, snappyHexMesh mesh S, wall functions, convergence monitor, PDF reports).

## Experiment

J. Zerihan, *An investigation into the aerodynamics of wings in ground effect*, PhD thesis, University
of Southampton, 2001 (chapter 3, single-element wing in freestream; see also Zerihan & Zhang,
J. Aircraft 37(6), 2000). An 80 % model of the main element of the Tyrrell 026 front wing (1998), with
endplates, tested in the 2.1 m × 1.7 m wind tunnel 750 mm above the belt (3.36 c, freestream).

| | |
|---|---|
| Profile | Table 2 of the thesis (`tools/t026_table2.csv`), modified NASA GA(W) LS(1)-0413 |
| Chord c, span b | 223.4 mm, 1100 mm (aspect ratio 4.92), trailing edge 1.65 mm |
| Endplates | 250 × 100 × 4 mm (position read on Figure 1a of the thesis) |
| Speed, Reynolds | about 30 m/s, Re = 0.43–0.46 × 10⁶, free transition |
| Incidence | quoted α; true incidence = α + 2.45° (reference α = 1°: endplates parallel to the ground) |
| Coefficients | C_L = L/(q c b) (positive = downforce), C_D = D/(q c b) |

`experiment.csv` (read from Figure 4 of the thesis; the C_L at 1° is also given in the text):

| α quoted | α true | C_L | C_D | flow |
|---|---|---|---|---|
| 1° | 3.45° | 0.69 | 0.020 | attached |
| 11° | 13.45° | 1.34 | 0.098 | maximum downforce, trailing-edge separation starts at the centre |
| 13° | 15.45° | 1.32 | 0.14 | separated at the centre |

The default run is α = 11°: there a1 matters (separation). α = 1° (attached) is a check of the setup.

## What differs from `cfd/rearWing`

Everything else is the FW43B case; the changes are marked `[validation]` in the files.

| File | Change |
|---|---|
| `Allrun.pre` | `-alpha <deg>` (default 11): `tools/t026STL.py` writes `constant/triSurface/wing.stl` (regions `wing`, `endplates`); `elements.stl` = region `wing` |
| `tools/t026STL.py`, `tools/t026_table2.csv` | STL from the coordinates of the thesis, in the frame of the FW43B case (x flow, y span, z up, suction side down) |
| `system/meshParameters` | geometry names; `chord 0.2234`, `span 1.1` (lRef = c, Aref = c b); `lengthScale`, `heightScale` (domain and boxes the same size in chords as for the FW43B); `uInf 30`, `nu 1.5e-5` (Re = 0.447e6); **`yPlus 26`** (see below); `sectionY 0` (section at the centre of the span) |
| `scripts/geometryInfo` | single element; lRef, Aref fixed; L, H in chords; `wingTop`, `wingHeight` for the post-processing |
| `scripts/postSetup`, `tools/viewMesh.py`, `tools/report/images.py` | sections, profiles and image frames on the wing (not on the scaled H of the mesh); new image `wssx_bottom`: skin friction along the flow on the suction side (separation line) |
| `scripts/report` | incidence and a1 in the tables; section *Comparison with the experiment*; `report/results/validation.csv` |
| `constant/caseSettings`, `constant/turbulenceProperties` | `a1` of k-ω SST as a parameter |
| `system/surfaceFeatureExtractDict` | output `wing.eMesh` |

**y+.** In the FW43B case the surface cell follows the layers (last layer / `finalLayerThickness`). At
Re = 0.45 × 10⁶ the design y+ = 100 would give 13.5 mm surface cells (c/16) and a layer stack thicker than
the boundary layer. With y+ = 26 the surface cell is 3.5 mm = c/64, the same as the FW43B mesh S in
chords (5.6 mm / 361 mm), with the same wall functions, 4 layers and ratio 1.2. About 2 million cells.

## Runs on the cluster

```
cd /global-scratch/$USER/FW43B-flexible-wing/cfd/validationT026
qsub -N T11_a031 -v A1=0.31 OpenFOAMJob_validation.sh        # one run: a1 = 0.31, alpha 11
./submit_a1_sweep.sh                                         # a1 = 0.31 0.4 0.5 0.7 1
./submit_a1_sweep.sh 0.31 0.6 0.8                            # other values
ALPHA=1 ./submit_a1_sweep.sh 0.31 1                          # check at alpha 1
```

Variables of `OpenFOAMJob_validation.sh`: `A1` (default 0.31), `ALPHA` (11), `MESH` (S), `SPEED` (30).
Each run: `results/validation/T026/alpha<ALPHA>/a1_<A1>/` and its `rigid.tgz`-style archive
(`pack_rigid.sh` of the URANS case).

## Results and choice of a1

Each run writes `report/results/validation.csv` (CFD against experiment). On the PC, with the runs
copied to `results/validation/T026/alpha11/`:

```
python3 cfd/validationT026/tools/a1Table.py results/validation/T026/alpha11 -o results/validation/T026/a1_study
```

It prints the table, the a1 with the smallest C_L error and the one with the smallest combined error,
and writes `a1_study.csv` and `a1_study.png` (C_L and C_D against a1, with the experiment).

Limits to report: free transition in the experiment (laminar bubble at x/c ≈ 0.3), fully turbulent CFD;
Re 0.45 × 10⁶ against 1–2 × 10⁶ for the FW43B; endplate position read from a figure; experimental values
read from a plot (about ±0.01 on C_L); wind-tunnel walls not modelled (no corrections in the thesis).
