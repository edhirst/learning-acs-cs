#!/usr/bin/env bash
# Local (CPU) campaign: ACS unit for every case, then the CS unit from the
# unit's OWN certified ACS (poly mode, extra/poly_cs.py), then the table.
#
#   bash extra/local_campaign.sh [SEED] [OUT_ROOT] [NPAR]
#
# Defaults: seed 7, outputs/local_campaign/s<seed>, 8 parallel single-threaded
# units. Python: $PYTHON or .venv/bin/python. S4 has no ACS: its CS unit
# refuses (no certified base), which is the expected outcome.
set -euo pipefail
SEED=${1:-7}
OUT=${2:-outputs/local_campaign/s${SEED}}
NPAR=${3:-8}
PY=${PYTHON:-.venv/bin/python}
cd "$(dirname "$0")/.."
mkdir -p "$OUT/logs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1
export PY OUT SEED

CASES="s2 s4 s6 s1xs1 s2xs2 s3xs1 s3xs3 s4xs2 s5xs1 s4xs1xs1"

run_acs () {
  "$PY" extra/campaign_unit.py --case "$1" --task acs --seed "$SEED" \
      --out "$OUT/$1_acs" > "$OUT/logs/$1_acs.txt" 2>&1 || echo "ACS $1 failed"
}
run_cs () {          # degree ladder 2,3; three-factor products "2t" for S4xS1xS1
  local classes=2,3
  [ "$1" = s4xs1xs1 ] && classes=2t
  "$PY" extra/campaign_unit.py --case "$1" --task cs --seed "$SEED" --cs-mode poly \
      --poly-degrees "$classes" --base-from "$OUT/$1_acs" \
      --out "$OUT/$1_cs" > "$OUT/logs/$1_cs.txt" 2>&1 || echo "CS $1 failed"
}
export -f run_acs run_cs

# ---- phase A: ACS (harmonic ladder, attempts plain -> Dirichlet -> orientation)
printf '%s\n' $CASES | xargs -P "$NPAR" -n 1 bash -c 'run_acs "$1"' _
# ---- phase B: CS from the unit's own certified ACS (S4 refuses: no certified base)
printf '%s\n' $CASES | xargs -P "$NPAR" -n 1 bash -c 'run_cs "$1"' _

"$PY" extra/campaign_summary.py "$OUT" | tee "$OUT/summary.md"
