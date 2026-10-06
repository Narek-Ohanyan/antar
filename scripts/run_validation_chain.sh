#!/usr/bin/env bash
# Baseline (2019) chain on the 80-point validation grid with the seasonal atmosphere: TOPOHYDRO -> XYLEM -> REFUGIUM.
# Each step reads the cached Drive inputs written by the first, so only that first step touches Drive. Same retry loop
# as run_dense_chain.sh. The scenario run (needs the ISIMIP3b atmosphere files) is a separate step.
set -o pipefail
cd "$(dirname "$0")/.."
MAX_TRIES=${MAX_TRIES:-10}
WAIT_S=${WAIT_S:-900}

run_step() {
  local name="$1"; shift
  for ((i = 1; i <= MAX_TRIES; i++)); do
    echo "=== [$name] attempt $i/$MAX_TRIES at $(date) ==="
    if "$@" 2>&1 | tee "/tmp/valchain_${name}_try${i}.log"; then
      echo "=== [$name] OK at $(date) ==="
      return 0
    fi
    echo "=== [$name] FAILED attempt $i; waiting ${WAIT_S}s ==="
    sleep "$WAIT_S"
  done
  echo "=== [$name] GAVE UP after $MAX_TRIES attempts ==="
  return 1
}

run_step topohydro python3 scripts/run_topohydro_grid.py &&
run_step xylem python3 scripts/fit_xylem_mechanistic_hazard.py &&
run_step refugium python3 scripts/fit_refugium_viability.py
