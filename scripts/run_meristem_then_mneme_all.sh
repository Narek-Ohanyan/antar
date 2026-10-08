#!/usr/bin/env bash
# After the MERISTEM water-deficit extraction has finished: refit the niche models with it, apply them at the nodes, then rerun MNEME on every qualifying forest pixel (so that Drive is
# never read by two jobs at once). Each step stops the chain if it fails.
set -o pipefail
cd "$(dirname "$0")/.."
LOG_CWD=${LOG_CWD:-/tmp/meristem_cwd.log}
echo "CHAIN: waiting for the MERISTEM water-deficit extraction ($LOG_CWD)"
until grep -q "wrote _real_cwd_for_meristem.npz" "$LOG_CWD" 2>/dev/null || grep -q "Traceback" "$LOG_CWD" 2>/dev/null; do sleep 30; done
if ! grep -q "wrote _real_cwd_for_meristem.npz" "$LOG_CWD"; then echo "CHAIN FAILED: the water-deficit extraction did not finish"; exit 1; fi
echo "CHAIN: water deficit done"
mkdir -p configs/fitted/archive
cp configs/fitted/meristem_adult_niche.yaml configs/fitted/archive/meristem_adult_niche_proxy_cwd.yaml
python3 scripts/fit_meristem_adult_niche.py || { echo "CHAIN FAILED: niche refit"; exit 1; }
echo "CHAIN: niche refit done"
python3 scripts/apply_meristem_niche.py || { echo "CHAIN FAILED: applying the niche"; exit 1; }
echo "CHAIN: niche applied at the nodes"
cp configs/fitted/mneme_hazard_panel_2010_2019_dense.yaml configs/fitted/archive/mneme_hazard_panel_k25.yaml
cp configs/fitted/mneme_vitality_response_dense.yaml configs/fitted/archive/mneme_vitality_response_k25.yaml
./scripts/run_mneme_forest_extract.sh --pixels-per-node=100000 || { echo "CHAIN FAILED: MNEME extraction"; exit 1; }
echo "CHAIN: MNEME pixels extracted"
./scripts/run_mneme_forest_fit.sh --pixels-per-node=100000 || { echo "CHAIN FAILED: MNEME fit"; exit 1; }
echo "CHAIN: ALL DONE"
