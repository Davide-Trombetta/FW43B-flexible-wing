# CFD (OpenFOAM v2306)

- `rearWing/`: the CFD case of the isolated rear wing (steady RANS, k-ω SST, `simpleFoam`), automated from the STL
  of one flap angle to the PDF reports; see `rearWing/README.md`. In short:

  ```bash
  cd cfd/rearWing
  ./Allrun.pre -angle 0 -mesh S -cores 8     # mesh of the 0° STL of ../../geometry
  ./Allrun -speed 70                         # run at 70 m/s, until Cd and Cl are converged
  ```

  The wall pressure for the structural load transfer (`fsi/map_pressure.py`) is written in
  `rearWing/postProcessing/wingSurfaces/<time>/walls.vtk`; one line per run (flap angle, speed, mesh, Cd, Cl, drag,
  downforce, efficiency) in `rearWing/report/results/summary.csv`.

Meshes, time directories, `postProcessing/`, `logs/` and `report/` are not versioned (see `rearWing/.gitignore`).
Run the cases on the Linux filesystem of WSL (not under `/mnt/c`), e.g. copying the case to `$FOAM_RUN` and setting
`GEOMETRY_DIR` to the `geometry/` folder of the repository.
