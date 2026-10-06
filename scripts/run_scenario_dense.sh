#!/usr/bin/env bash
# The dense-grid scenario run on its own (the validation-grid run is done). Reads only cached inputs and local ISIMIP files.
set -o pipefail
cd "$(dirname "$0")/.."
python3 scripts/run_future_projections.py --dense --workers "${WORKERS:-5}" 2>&1 | tee /tmp/future_dense.log && echo "=== DENSE SCENARIO RUN DONE $(date) ==="
