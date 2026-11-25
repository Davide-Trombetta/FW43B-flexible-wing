#!/bin/bash
#------------------------------------------------------------------------------
#  OpenFOAMJob_simple.sh -- Grid Engine job: rigid STEADY CFD (simpleFoam) of one
#  case, k-omega SST with a1 = 1, on CALIMERO.
#
#  Same as OpenFOAMJob.sh (URANS), but:
#    - source case  cfd/rearWing          (steady, simpleFoam, stopped by
#                                          scripts/convergenceMonitor)
#    - results in   results/simple/U<speed>/flap<angle>/rigid
#    - a1 = 1 in constant/turbulenceProperties of the copy (taken from the URANS
#      case if the steady case does not have it), checked before the run
#
#  Usage (from the folder of this script):
#     qsub -N S86_f0 -v ANGLE=0,SPEED=86 OpenFOAMJob_simple.sh
#------------------------------------------------------------------------------
#$ -S /bin/bash
#$ -N FW43B_SIMPLE
#$ -q all.q
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y

source /etc/profile
module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv
set -eo pipefail

ANGLE="${ANGLE:?Specificare ANGLE}"
SPEED="${SPEED:?Specificare SPEED}"
case "$ANGLE" in -4|-2|0|+2|+4) ;; *) echo "Angolo non valido: $ANGLE"; exit 1 ;; esac
case "$SPEED" in 44|86) ;; *) echo "Velocita non valida: $SPEED"; exit 1 ;; esac

BASE="/global-scratch/$(id -un)"
REPO="$BASE/FW43B-flexible-wing"
SOURCE="$REPO/cfd/rearWing"
URANS="$REPO/cfd/rearWingURANS"
export GEOMETRY_DIR="$REPO/geometry"
CASE_DIR="$BASE/results/simple/U${SPEED}/flap${ANGLE}/rigid"

[ -d "$SOURCE" ] || { echo "Caso sorgente non trovato: $SOURCE"; exit 1; }
[ -f "$GEOMETRY_DIR/FW43B_rearWing_flap_${ANGLE}deg.stl" ] || { echo "Geometria STL non trovata"; exit 1; }
if ls -d "$SOURCE"/processor* "$SOURCE"/constant/polyMesh "$SOURCE"/tools/stlInfo/build > /dev/null 2>&1; then
    echo "Il caso sorgente non e' pulito (mesh o build presenti): $SOURCE"; exit 1
fi

mkdir -p "$(dirname "$CASE_DIR")"
mkdir "$CASE_DIR" || { echo "La cartella del caso esiste gia: $CASE_DIR"; exit 1; }
cp -a "$SOURCE/." "$CASE_DIR/"
cd "$CASE_DIR"

# --- a1 = 1 (k-omega SST): from the URANS case if the steady one does not set it
TP=constant/turbulenceProperties
if ! grep -qE "^\s*a1\s" "$TP"; then
    echo "a1 non presente in $TP del caso stazionario: copio quello del caso URANS"
    cp "$URANS/$TP" "$TP"
fi
A1=$(foamDictionary "$TP" -entry RAS/kOmegaSSTCoeffs/a1 -value 2>/dev/null || true)
echo "a1 = ${A1:-(non trovato)}"
[ "$A1" = 1 ] || { echo "a1 non vale 1 in $TP: controllare il file"; exit 1; }

echo "Nodo: $(hostname)   Core assegnati: $NSLOTS   Angolo: $ANGLE   Velocita: $SPEED m/s"
echo "Solutore: $(foamDictionary system/controlDict -entry application -value)"

./Allrun.pre -angle "$ANGLE" -mesh S -cores "$NSLOTS"
./Allrun -speed "$SPEED"
echo "Simulazione completata: $CASE_DIR"

bash "$URANS/pack_rigid.sh" -mesh "$CASE_DIR"
cp "$SGE_STDOUT_PATH" "$(dirname "$CASE_DIR")/" 2>/dev/null || true
