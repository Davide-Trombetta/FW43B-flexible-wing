#!/bin/bash
#------------------------------------------------------------------------------
#  pack_rigid.sh -- after a rigid job: the files to bring to the PC
#
#  Usage:  cluster/pack_rigid.sh [-mesh] <case folder>
#     e.g. cluster/pack_rigid.sh -mesh ~/runs/U86/flap0/rigid
#
#  Writes  <case folder>.tgz : everything the reports and the FEM need
#  (logs, postProcessing with walls.vtk and the extracted planes, settings,
#  partial report), WITHOUT the decomposed fields and mesh (processor*/), which
#  stay on the cluster for the deformed CFD.
#
#  -mesh   also the snappyHexMesh mesh, reconstructed in constant/polyMesh
#          (reconstructParMesh -constant; needs the OpenFOAM module loaded):
#          an exact copy of the mesh for the deformed CFD if the case on the
#          cluster is lost, and the mesh images of 'scripts/report mesh' on
#          the PC. Also the number of layers per face written by snappyHexMesh
#          (0/nSurfaceLayers, reconstructed with reconstructPar), for the layer
#          image. A few hundred MB more.
#          Without -mesh, constant/polyMesh (there only the blockMesh
#          background) is left out.
#------------------------------------------------------------------------------
set -e
MESH=no
[ "$1" = -mesh ] && { MESH=yes; shift; }
[ $# -eq 1 ] || { sed -n 4,6p "$0"; exit 1; }
C=${1%/}
[ -d "$C/postProcessing" ] || { echo "$C: no postProcessing/ (run not done?)"; exit 1; }
B=$(basename "$C")

EXCL=(--exclude='processor*' --exclude="$B/tools/*/build" --exclude="$B/tools/bin"
      --exclude="$B/[0-9]" --exclude="$B/[0-9].[0-9]*")

if [ $MESH = yes ]; then
    command -v reconstructParMesh > /dev/null ||
        { echo "reconstructParMesh not found: module load openfoam/2306-gcc-13.2.0-tnytlfv"; exit 1; }
    [ -d "$C/processor0/constant/polyMesh" ] || { echo "$C: no decomposed mesh"; exit 1; }
    echo "reconstructing the mesh in $C/constant/polyMesh ..."
    (cd "$C" && reconstructParMesh -constant > logs/log.reconstructParMesh 2>&1) ||
        { echo "reconstructParMesh failed: see $C/logs/log.reconstructParMesh"; exit 1; }
    tail -5 "$C/logs/log.reconstructParMesh" | grep -q "^End" ||
        { echo "reconstructParMesh did not end: see $C/logs/log.reconstructParMesh"; exit 1; }
    grep -m1 "cells:" "$C/logs/log.reconstructParMesh" || true
    # layers per face (snappyHexMesh writes nSurfaceLayers in processor*/0)
    (cd "$C" && reconstructPar -time 0 -fields '(nSurfaceLayers)' > logs/log.reconstructParLayers 2>&1) ||
        echo "warning: nSurfaceLayers not reconstructed (see $C/logs/log.reconstructParLayers): no layer image"
else
    EXCL+=(--exclude="$B/constant/polyMesh")
fi

tar cf "$C.tar" -C "$(dirname "$C")" "${EXCL[@]}" "$B"
# the time folders are left out; with -mesh only 0/nSurfaceLayers is added
[ $MESH = yes ] && [ -f "$C/0/nSurfaceLayers" ] &&
    tar rf "$C.tar" -C "$(dirname "$C")" "$B/0/nSurfaceLayers"
gzip -f "$C.tar" && mv -f "$C.tar.gz" "$C.tgz"
du -h "$C.tgz"
