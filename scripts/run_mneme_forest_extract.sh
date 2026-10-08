#!/usr/bin/env bash
# Forest-pixel extraction for MNEME, retried after a pause when Drive drops out. The script saves its progress every few nodes and resumes from it, so a retry costs only the node in flight.
set -o pipefail
cd "$(dirname "$0")/.."
MAX_TRIES=${MAX_TRIES:-6}
WAIT_S=${WAIT_S:-300}
for ((i = 1; i <= MAX_TRIES; i++)); do
  echo "=== [mneme forest extraction] attempt $i/$MAX_TRIES at $(date) ==="
  if python3 scripts/extract_mneme_forest_pixels.py --dense "$@"; then
    echo "=== [mneme forest extraction] OK at $(date) ==="
    exit 0
  fi
  echo "=== FAILED attempt $i; waiting ${WAIT_S}s ==="
  sleep "$WAIT_S"
done
echo "=== GAVE UP after $MAX_TRIES attempts ==="
exit 1
