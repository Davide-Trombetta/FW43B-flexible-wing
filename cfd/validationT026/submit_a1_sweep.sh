#!/bin/bash
#------------------------------------------------------------------------------
#  submit_a1_sweep.sh -- submits one OpenFOAMJob_validation.sh per value of a1
#
#  Usage (from the folder of this script, on the login node: it only calls qsub):
#     ./submit_a1_sweep.sh                        a1 = 0.31 0.4 0.5 0.7 1, alpha 11
#     ./submit_a1_sweep.sh 0.31 0.6 1             these values
#     ALPHA=1 MESH=S ./submit_a1_sweep.sh 0.31 1  another incidence or mesh
#  Each run has its own folder (results/validation/T026/alpha<ALPHA>/a1_<A1>);
#  a run whose folder already exists is skipped.
#------------------------------------------------------------------------------
cd "${0%/*}" || exit 1
ALPHA="${ALPHA:-11}"
MESH="${MESH:-S}"
VALUES="${*:-0.31 0.4 0.5 0.7 1}"
BASE="/global-scratch/$(id -un)/results/validation/T026/alpha${ALPHA}"
mkdir -p "$(dirname "$BASE")"          # folder of the job logs (-o)

for a in $VALUES; do
    if [ -d "$BASE/a1_$a" ]; then
        echo "a1 = $a: $BASE/a1_$a exists, skipped"
        continue
    fi
    NAME="T${ALPHA}_a${a//./}"
    qsub -N "$NAME" -o "$(dirname "$BASE")/" -v A1="$a",ALPHA="$ALPHA",MESH="$MESH" OpenFOAMJob_validation.sh
done
qstat
