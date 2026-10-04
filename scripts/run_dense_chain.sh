#!/usr/bin/env bash
# Runs the dense-grid steps in order. A step that fails (typically DriveCoverageError during a
# Google Drive disruption) is retried after a long wait instead of killing the whole chain: the
# first unattended overnight attempt lost ~33 hours because one outage window ended the chain and
# nothing retried it.
set -o pipefail
cd "$(dirname "$0")/.."
MAX_TRIES=${MAX_TRIES:-10}
WAIT_S=${WAIT_S:-900}

run_step() {
  local name="$1"; shift
  for ((i = 1; i <= MAX_TRIES; i++)); do
    echo "=== [$name] attempt $i/$MAX_TRIES at $(date) ==="
    if "$@" 2>&1 | tee "/tmp/dense_${name}_try${i}.log"; then
      echo "=== [$name] OK at $(date) ==="
      return 0
    fi
    echo "=== [$name] FAILED attempt $i; waiting ${WAIT_S}s for Drive to recover ==="
    sleep "$WAIT_S"
  done
  echo "=== [$name] GAVE UP after $MAX_TRIES attempts ==="
  return 1
}

run_step topohydro python3 scripts/run_topohydro_grid.py --dense &&
run_step xylem python3 scripts/fit_xylem_mechanistic_hazard.py --dense &&
run_step refugium python3 scripts/fit_refugium_viability.py --dense &&
run_step mneme python3 scripts/fit_mneme_hazard_panel.py --dense
