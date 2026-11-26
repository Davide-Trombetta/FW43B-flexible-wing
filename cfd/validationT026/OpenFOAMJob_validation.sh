#!/bin/bash
#------------------------------------------------------------------------------
#  OpenFOAMJob_validation.sh -- Grid Engine job: one run of the validation
#  case (Tyrrell 026 single-element wing in freestream, Zerihan 2001), steady
#  RANS (simpleFoam, k-omega SST) with the workflow of cfd/rearWing, on CALIMERO.
#
#  Variables (qsub -v):
#    A1     k-omega SST a1 of this run           (default 0.31, the standard value)
#    ALPHA  quoted incidence [deg]               (default 11; true = ALPHA + 2.45)
#    MESH   mesh size XS, S, M or L              (default S)
#    SPEED  free-stream velocity [m/s]           (default: uInf of meshParameters, 30)
#  Results in /global-scratch/<user>/results/validation/T026/alpha<ALPHA>/a1_<A1>
#  (report/results/validation.csv: CFD against experiment for this run)
#
#  Usage (from the folder of this script):
#     qsub -N T11_a031 -v A1=0.31 OpenFOAMJob_validation.sh
#     qsub -N T11_a1   -v A1=1,ALPHA=11 OpenFOAMJob_validation.sh
#  several values of a1 at once: submit_a1_sweep.sh
#------------------------------------------------------------------------------
#$ -S /bin/bash
#$ -N T026_VALID
#$ -q all.q
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y

source /etc/profile
module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv
set -eo pipefail

A1="${A1:-0.31}"
ALPHA="${ALPHA:-11}"
MESH="${MESH:-S}"
echo "$A1" | grep -qE '^[0-9]+([.][0-9]*)?$' || { echo "a1 non valido: $A1"; exit 1; }
echo "$ALPHA" | grep -qE '^[+-]?[0-9]+([.][0-9]*)?$' || { echo "Angolo non valido: $ALPHA"; exit 1; }
case "$MESH" in XS|S|M|L) ;; *) echo "Mesh non valida: $MESH"; exit 1 ;; esac

BASE="/global-scratch/$(id -un)"
REPO="$BASE/FW43B-flexible-wing"
SOURCE="$REPO/cfd/validationT026"
CASE_DIR="$BASE/results/validation/T026/alpha${ALPHA}/a1_${A1}"

[ -d "$SOURCE" ] || { echo "Caso sorgente non trovato: $SOURCE"; exit 1; }
if ls -d "$SOURCE"/processor* "$SOURCE"/constant/polyMesh "$SOURCE"/tools/stlInfo/build > /dev/null 2>&1; then
    echo "Il caso sorgente non e' pulito (mesh o build presenti): $SOURCE"; exit 1
fi
mkdir -p "$(dirname "$CASE_DIR")"
mkdir "$CASE_DIR" || { echo "La cartella del caso esiste gia: $CASE_DIR"; exit 1; }
cp -a "$SOURCE/." "$CASE_DIR/"
cd "$CASE_DIR"

# a1 of this run: constant/caseSettings, read by constant/turbulenceProperties
sed -i "s|^\(a1[[:space:]]\{1,\}\)[^;]*;|\1$A1;|" constant/caseSettings
A1CHECK=$(foamDictionary constant/turbulenceProperties -entry RAS/kOmegaSSTCoeffs/a1 -value)
echo "a1 = $A1CHECK (constant/turbulenceProperties)"
awk -v a="$A1" -v b="$A1CHECK" 'BEGIN { exit !(a + 0 == b + 0) }' ||
    { echo "a1 non impostato correttamente"; exit 1; }

echo "Nodo: $(hostname)   Core: $NSLOTS   alpha: $ALPHA deg   a1: $A1   mesh: $MESH"
./Allrun.pre -alpha "$ALPHA" -mesh "$MESH" -cores "$NSLOTS"
if [ -n "${SPEED:-}" ]; then ./Allrun -speed "$SPEED"; else ./Allrun; fi
echo "Simulazione completata: $CASE_DIR"
cat report/results/validation.csv 2> /dev/null || true

bash "$REPO/cfd/rearWingURANS/pack_rigid.sh" -mesh "$CASE_DIR"
cp "$SGE_STDOUT_PATH" "$(dirname "$CASE_DIR")/" 2>/dev/null || true
