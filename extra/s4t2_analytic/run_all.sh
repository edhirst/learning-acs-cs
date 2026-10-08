#!/usr/bin/env bash
# Run the manuscript's verification scripts and regenerate the reference outputs in data/.
# The scripts themselves live in manuscript_draft/code/ (single source of truth).
# Usage, from anywhere:  extra/s4t2_analytic/run_all.sh    (set PYTHON to choose the interpreter)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CODE="$HERE/../../manuscript_draft/code"
PY="${PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1    # leave no __pycache__ in manuscript_draft/

echo "== manuscript_draft/code/verify_construction.py"
"$PY" "$CODE/verify_construction.py" | tee "$HERE/data/verify_construction.txt"

echo "== manuscript_draft/code/sphere_products.py"
(cd "$CODE" && "$PY" sphere_products.py) | tee "$HERE/data/sphere_products.txt"

echo "== extra/s4t2_analytic/export_table.py"
"$PY" "$HERE/export_table.py" "$HERE/data/sphere_products.csv"
