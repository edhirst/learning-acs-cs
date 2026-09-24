#!/bin/bash
#PBS -q serial
#PBS -l mem=24G
#PBS -l nodes=1:ppn=8
#PBS -N J_S4xS1xS1_cs
#PBS -m abe
#PBS -k oe
#PBS -j oe

# Submit after ACS succeeds, for example:
#   ACS_JOB=$(qsub scratch/submit_S4xS1xS1_ACS.sh)
#   qsub -W depend=afterok:${ACS_JOB} scratch/submit_S4xS1xS1_CS.sh
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
ACS_DIR="${ACS_OUT_DIR:-$PBS_O_WORKDIR/outputs/s4xs1xs1_ambient/acs_s${SEED}}"
OUT_DIR="${CS_OUT_DIR:-$PBS_O_WORKDIR/outputs/s4xs1xs1_ambient/cs_s${SEED}}"
if [[ ! -f "$ACS_DIR/final.pt" ]]; then
    echo "Missing ACS checkpoint: $ACS_DIR/final.pt" >&2
    exit 3
fi
mkdir -p "$OUT_DIR"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export PYTHONUNBUFFERED=1

cd "$CODE_ROOT"
python3 extra/campaign_unit.py \
    --case s4xs1xs1 --formulation ambient --task cs --seed "$SEED" \
    --steps "${STEPS:-80000}" --hidden "${HIDDEN:-256,256,256,256}" \
    --stages "${STAGES:-2}" --tail --lbfgs "${LBFGS:-200}" \
    --base-from "$ACS_DIR" --wall-minutes "${WALL_MINUTES:-4000}" \
    --out "$OUT_DIR"

echo "CS output: $OUT_DIR"
