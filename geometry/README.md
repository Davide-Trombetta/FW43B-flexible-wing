# Geometry

`FW43B_rearWing_flap_<angle>deg.stl` (angle = −4, −2, 0, +2, +4): CFD surface of the Williams FW43B(AR) 2021 rear
wing, ASCII multi-solid STL, **units in metres**. Source: model extracted from the Williams 2021 AR app.

**Flap angle**: rigid rotation of the flap about a spanwise axis through its leading edge (x = 4.0235 m,
z = 0.9270 m); positive = steeper flap (more downforce and drag). 0° is the original geometry (flap chord at 66.3°
to the horizontal). The central part of the flap and the DRS link follow the rotation smoothly. The hinge position
is an assumption.

**Regions** (solids → OpenFOAM patches `wing_<region>` with snappyHexMesh):
`rearWing_mainPlane`, `rearWing_flap`, `rearWing_drsPod` (DRS pod, swan-neck arms, DRS link),
`rearWing_endplate_yPos`, `rearWing_endplate_yNeg`, `rearWing_pylons`, `rearWing_mount` (attachment to the gearbox).

**Reference frame**: x downstream (flow along +x), y spanwise (symmetry plane y = 0), z up, ground at z = 0,
x = 0 at the front axle; rake (1.96°) already applied.

Every CFD case must use the STL of its angle with **identical** meshing settings.
