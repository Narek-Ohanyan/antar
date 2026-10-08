#!/usr/bin/env bash
# MNEME's panel fit (climate of the forest cells streamed from Drive, then events, placebo and hazard) and the vitality-response analysis. The panel step checkpoints each year, so a retry after a
# Drive outage resumes from the last finished year.
set -o pipefail
cd "$(dirname "$0")/.."
MAX_TRIES=${MAX_TRIES:-6}
WAIT_S=${WAIT_S:-300}
for ((i = 1; i <= MAX_TRIES; i++)); do
  echo "=== [mneme panel] attempt $i/$MAX_TRIES at $(date) ==="
  if python3 scripts/fit_mneme_hazard_panel.py --dense "$@"; then
    python3 scripts/fit_mneme_vitality_response.py && echo "=== [mneme] OK at $(date) ===" && exit 0
    echo "=== vitality analysis failed ==="; exit 1
  fi
  echo "=== FAILED attempt $i; waiting ${WAIT_S}s ==="
  sleep "$WAIT_S"
done
echo "=== GAVE UP after $MAX_TRIES attempts ==="
exit 1
