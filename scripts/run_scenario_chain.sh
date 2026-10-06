#!/usr/bin/env bash
# The corrected scenario runs (seasonal, scenario-dependent atmosphere): validation grid first (quick), then the dense Armenia-only grid.
# Both read only cached inputs and local ISIMIP files, so no Drive access is needed; per-member checkpoints make a restart resume.
set -o pipefail
cd "$(dirname "$0")/.."
WORKERS=${WORKERS:-6}
echo "=== scenario run, validation grid, $(date) ==="
python3 scripts/run_future_projections.py --workers "$WORKERS" 2>&1 | tee /tmp/future_coarse.log || { echo "=== coarse run FAILED $(date) ==="; exit 1; }
echo "=== scenario run, dense grid, $(date) ==="
python3 scripts/run_future_projections.py --dense --workers "$WORKERS" 2>&1 | tee /tmp/future_dense.log || { echo "=== dense run FAILED $(date) ==="; exit 1; }
echo "=== SCENARIO CHAIN DONE $(date) ==="
