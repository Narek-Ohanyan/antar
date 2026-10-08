"""MNEME's continuous vitality response, from the forest-pixel panel saved by ``fit_mneme_hazard_panel.py`` (no network, no new extraction).

Asks, with pixel fixed effects, how the standardised kNDVI anomaly of a forest pixel moves with the year-to-year departure of the climatic water deficit of its cell, whether the weather
predicts an unseen year, whether the response grows with canopy height, and, per species group, whether the response is steeper where XYLEM's mechanistic hazard is higher: an
empirical check of XYLEM's drought-stress ranking that needs no dieback events. It is a vitality response, not mortality: a dry year can lower kNDVI without killing a tree.
See ``antar.hazard.vitality.analyse_vitality`` for the models and what was fixed in advance.
"""
import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from antar.hazard.forest_panel import GROUP_OF_CLASS  # noqa: E402
from antar.hazard.vitality import analyse_vitality  # noqa: E402
from run_topohydro_grid import CACHE_DIR, DENSE_GRID_COLS, DENSE_GRID_ROWS, extract_static_grid_inputs  # noqa: E402

PANEL_PATH = CACHE_DIR / "mneme_forest_panel_dense.csv.gz"
XYLEM_PATH = ROOT / "configs" / "fitted" / "xylem_mechanistic_hazard_2019_dense.yaml"
OUT_PATH = ROOT / "configs" / "fitted" / "mneme_vitality_response_dense.yaml"
CLIMATE_COLS = ["cwd_mm", "vpd24_mean_kpa", "gdd_annual", "t_mean_c_annual"]


def xylem_by_node(static):
    """{group: {node index: mean mechanistic hazard}} for the groups XYLEM fitted, matching its cells to the grid nodes by position."""
    d = yaml.load(open(XYLEM_PATH), Loader=yaml.CSafeLoader)
    key = lambda lat, lon: (round(float(lat), 6), round(float(lon), 6))
    pos = {key(a, b): i for i, (a, b) in enumerate(zip(static["lats"], static["lons"]))}
    out = {}
    for g, v in d["groups"].items():
        if v.get("status") != "ok":
            continue
        out[g] = {pos[key(c["lat"], c["lon"])]: float(c["h_mech_mean"]) for c in v["cells"] if key(c["lat"], c["lon"]) in pos}
    return out


def main():
    if not PANEL_PATH.exists():
        sys.exit(f"{PANEL_PATH.name} not found: run scripts/fit_mneme_hazard_panel.py --dense first")
    panel = pd.read_csv(PANEL_PATH)
    static = extract_static_grid_inputs(grid_rows=DENSE_GRID_ROWS, grid_cols=DENSE_GRID_COLS)
    res = analyse_vitality(panel, CLIMATE_COLS, xylem=xylem_by_node(static), group_of_class=GROUP_OF_CLASS)
    res = {"run_date": datetime.date.today().isoformat(), "grid": "dense_armenia_stride7", "panel_years": [int(panel["year"].min()), int(panel["year"].max())],
           "what": "vitality response (standardised kNDVI anomaly), not mortality", **res}
    OUT_PATH.write_text(yaml.dump(res, sort_keys=False, default_flow_style=False))
    a = res["model_a"]["terms"][0]
    print(f"=== {res['n_pixel_years']} pixel-years, {res['n_pixels']} pixels, {res['n_cells']} cells ({res['n_excluded_recent_disturbance']} pixel-years left out for recent harvest or fire) ===")
    print(f"=== drought response: {a['coef']:+.3f} SD of kNDVI per SD of CWD anomaly (95% CI {a['ci95'][0]:+.3f} to {a['ci95'][1]:+.3f}, p = {a['p']:.3g}); "
          f"out-of-sample R2 for an unseen year {res['leave_one_year_out']['cwd_only']['r2_out_of_sample']:+.4f} ===")
    for g, x in res["xylem_check"].items():
        if x.get("status") == "ok":
            t = x["model"]["terms"][1]
            print(f"=== XYLEM ranking check, {g}: interaction {t['coef']:+.3f} (95% CI {t['ci95'][0]:+.3f} to {t['ci95'][1]:+.3f}) -> {x['verdict']} ===")
        else:
            print(f"=== XYLEM ranking check, {g}: {x['status']} ===")
    print(f"=== wrote {OUT_PATH.name} ===")


if __name__ == "__main__":
    sys.exit(main())
