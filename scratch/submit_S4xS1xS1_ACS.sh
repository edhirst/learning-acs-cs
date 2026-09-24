#!/bin/bash
#PBS -q serial
#PBS -l mem=24G
#PBS -l nodes=1:ppn=8
#PBS -N J_S4xS1xS1_acs
#PBS -m abe
#PBS -k oe
#PBS -j oe

# Submit from either the repository itself or its parent directory.
set -euo pipefail
cd "$PBS_O_WORKDIR" || exit 1
mkdir -p outputs logs

archive_pbs_log () {
    JOBNUM="${PBS_JOBID%%.*}"
    LIVE_FILE="$HOME/${PBS_JOBNAME}.o${JOBNUM}"
    if [[ -f "$LIVE_FILE" ]]; then
        mv "$LIVE_FILE" "$PBS_O_WORKDIR/outputs/"
    fi
}
trap archive_pbs_log EXIT

module load miniconda3/22.11.1-gcc-9.4.0
eval "$(conda shell.bash hook)"
conda activate nijenhuis

if [[ -f extra/campaign_unit.py ]]; then
    CODE_ROOT="$PBS_O_WORKDIR"
elif [[ -f github/extra/campaign_unit.py ]]; then
    CODE_ROOT="$PBS_O_WORKDIR/github"
else
    echo "Cannot find extra/campaign_unit.py under $PBS_O_WORKDIR" >&2
    exit 2
fi

SEED="${SEED:-7}"
OUT_DIR="${ACS_OUT_DIR:-$PBS_O_WORKDIR/outputs/s4xs1xs1_ambient/acs_s${SEED}}"
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export PYTHONUNBUFFERED=1

cd "$CODE_ROOT"
python3 extra/campaign_unit.py \
    --case s4xs1xs1 --formulation ambient --task acs --seed "$SEED" \
    --steps "${STEPS:-10000}" --hidden "${HIDDEN:-256,256,256,256}" \
    --poly-steps "${POLY_STEPS:-3000}" --lbfgs "${LBFGS:-200}" \
    --wall-minutes "${WALL_MINUTES:-4000}" --out "$OUT_DIR"

# The CS job expects: outputs/s4xs1xs1_ambient/acs_s${SEED}/final.pt
echo "ACS output: $OUT_DIR"
