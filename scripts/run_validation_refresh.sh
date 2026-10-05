#!/usr/bin/env bash
# Refreshes the validation-grid (80-point, 25 Armenian) results so they carry every field the UI maps need:
# per-species 2019 water stress and P[V >= V*] (refugium), and every scenario quantity (future projections,
# run in parallel across the 45 GCM x SSP x horizon members). Same seeds as the earlier runs, so the fields that
# already existed must come out identical -- check with tests/test_refresh_regression.py (see IMPLEMENTATION_LOG).
set -o pipefail
cd "$(dirname "$0")/.."
MAX_TRIES=${MAX_TRIES:-5}
WAIT_S=${WAIT_S:-600}
WORKERS=${WORKERS:-5}

run_step() {
  local name="$1"; shift
  for ((i = 1; i <= MAX_TRIES; i++)); do
    echo "=== [$name] attempt $i/$MAX_TRIES at $(date) ==="
    if "$@" 2>&1 | tee "/tmp/valrefresh_${name}_try${i}.log"; then
      echo "=== [$name] OK at $(date) ==="
      return 0
    fi
    echo "=== [$name] FAILED attempt $i; waiting ${WAIT_S}s ==="
    sleep "$WAIT_S"
  done
  echo "=== [$name] GAVE UP after $MAX_TRIES attempts ==="
  return 1
}

run_step refugium python3 scripts/fit_refugium_viability.py &&
run_step future python3 scripts/run_future_projections.py --workers "$WORKERS"
