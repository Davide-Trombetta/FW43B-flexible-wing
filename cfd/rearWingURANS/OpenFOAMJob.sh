#!/bin/bash
#$ -S /bin/bash
#$ -N FW43B_URANS
#$ -q all.q
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y

source /etc/profile
set -euo pipefail

module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv

ANGLE="${ANGLE:?Specificare ANGLE}"
SPEED="${SPEED:?Specificare SPEED}"

case "$ANGLE" in
    -4|-2|0|+2|+4) ;;
    *) echo "Angolo non valido: $ANGLE"; exit 1 ;;
esac

case "$SPEED" in
    44|86) ;;
    *) echo "Velocita non valida: $SPEED"; exit 1 ;;
esac

BASE="/global-scratch/$(id -un)"
SOURCE="$BASE/FW43B-flexible-wing/cfd/rearWingURANS"
export GEOMETRY_DIR="$BASE/FW43B-flexible-wing/geometry"
CASE_DIR="$BASE/results/U${SPEED}/flap${ANGLE}/rigid"

if [ ! -d "$SOURCE" ]; then
    echo "Caso sorgente non trovato: $SOURCE"
    exit 1
fi

if [ ! -f "$GEOMETRY_DIR/FW43B_rearWing_flap_${ANGLE}deg.stl" ]; then
    echo "Geometria STL non trovata"
    exit 1
fi

mkdir -p "$(dirname "$CASE_DIR")"

if ! mkdir "$CASE_DIR"; then
    echo "La cartella del caso esiste gia o non puo essere creata: $CASE_DIR"
    exit 1
fi

cp -a "$SOURCE/." "$CASE_DIR/"
cd "$CASE_DIR"

echo "Nodo: $(hostname)"
echo "Core assegnati: $NSLOTS"
echo "Angolo: $ANGLE"
echo "Velocita: $SPEED m/s"

./Allrun.pre -angle "$ANGLE" -mesh S -cores "$NSLOTS"
./Allrun -speed "$SPEED"

echo "Simulazione completata: $CASE_DIR"

"$BASE/FW43B-flexible-wing/cfd/rearWingURANS/pack_rigid.sh" -mesh "$CASE_DIR"
