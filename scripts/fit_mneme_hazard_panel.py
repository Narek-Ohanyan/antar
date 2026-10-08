"""MNEME's person-period panel of forest pixels: dieback labels from the satellite vitality record crossed with the climate of the pixel's cell, and the hazard fit when there are events enough.

Replaces the first panel, which followed one 30 m pixel at each grid node whatever its land cover and read a harvest/fire layer that said "disturbed" everywhere after 2000 (so the
dieback rule could never fire). The frame is now forest pixels only (``extract_mneme_forest_pixels.py``) and the harvest/fire layer is the corrected export.

Rows: one per forest pixel and panel year (2010-2019). The panel years are the overlap of two real limits: the dieback rule needs a 10-year trailing window and a 2-year recovery window
inside the 2000-2024 record (so 2010-2022), and CHELSA-daily precipitation, which TOPOHYDRO's water balance needs, ends on 2019-12-31. Climate covariates belong to the node (the
climate cell of 30 arcseconds); every forest pixel of a node shares them and has its own vitality series, canopy height and elevation.

Design: Mundlak within/between split of the node's climate (cwd, vpd, gdd, mean temperature) plus the pixel's 2019 canopy height. Stand age, a lagged-drought cross-basis and the full
Eq. 7.3 stand terms need data that do not exist for Armenian plots and are not included; misclassification correction (Se/Sp) needs a labelled sample and is not applied. The stacked
hazard learner is fitted only when there are at least MIN_EVENTS_PER_PREDICTOR events per column of the design (a standard rule of thumb for stable coefficients); otherwise the
panel and the event rate are reported and no hazard is claimed.

Also saves the pixel-year table (data/_cache/mneme_forest_panel_dense.csv.gz) for the vitality-response analysis (``fit_mneme_vitality_response.py``).
"""
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import brier_score_loss, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    CACHE_DIR, DENSE_GRID_COLS, DENSE_GRID_ROWS, compute_forcing_for_year, extract_static_grid_inputs,
)
from mneme_inputs import PANEL_YEARS, VITALITY_YEARS  # noqa: E402

from antar.hazard import forest_panel as fp  # noqa: E402
from antar.hazard.observation import event_rate_interval, zero_event_upper_bound  # noqa: E402
from antar.hazard.panel import mundlak_decompose  # noqa: E402
from antar.hazard.pipeline import fit_stacked_hazard  # noqa: E402
from antar.validation.splits import spatial_block_ids  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT_PATH_DENSE = ROOT / "configs" / "fitted" / "mneme_hazard_panel_2010_2019_dense.yaml"
CHECKPOINT_PATH_DENSE = ROOT / "data" / "_mneme_panel_checkpoint_dense.json"
def _int_arg(name, default):
    for a in sys.argv[1:]:
        if a.startswith(name + "="):
            return int(a.split("=", 1)[1])
    return default


# the same cap the extraction was run with (see extract_mneme_forest_pixels.py): another cap reads its own pixel file
PIXELS_PER_NODE = _int_arg("--pixels-per-node", 25)
PIXELS_PATH = CACHE_DIR / f"mneme_forest_pixels_dense{'' if PIXELS_PER_NODE == 25 else f'_k{PIXELS_PER_NODE}'}.npz"
PANEL_PATH = CACHE_DIR / "mneme_forest_panel_dense.csv.gz"

CLIMATE_COLS = ["cwd_mm", "vpd24_mean_kpa", "gdd_annual", "t_mean_c_annual"]
SPATIAL_BLOCK_SIZE_DEG = 1.0
MIN_EVENTS_PER_PREDICTOR = 10
MIN_DECLINE_TO_RISE_RATIO = 2.0     # declared before the large run: the rule, applied to the mirrored series, must find at least this many times more declines than rises


def node_climate(static, nodes):
    """One row per node and panel year: the TOPOHYDRO forcing summaries of that year, checkpointed after each year so a Drive outage costs one year at most.

    The checkpoint is kept after a successful run (the climate of a fixed set of cells and years does not change, and the streaming is the slow part); it is used again only when its key
    (the cells and the years) matches.
    """
    key = hashlib.sha1((",".join(map(str, nodes.tolist())) + "|" + ",".join(map(str, PANEL_YEARS))).encode()).hexdigest()[:12]
    sub = fp.subset_static(static, nodes)
    saved = json.loads(CHECKPOINT_PATH_DENSE.read_text()) if CHECKPOINT_PATH_DENSE.exists() else {}
    rows = saved["rows"] if isinstance(saved, dict) and saved.get("key") == key else []
    done = {r["year"] for r in rows}
    if done:
        print(f"  resuming from checkpoint: years {sorted(done)} already computed", flush=True)
    for year in PANEL_YEARS:
        if year in done:
            continue
        cells = compute_forcing_for_year(sub, year)
        for k, c in enumerate(cells):
            if c is None:
                continue
            rows.append({"node": int(nodes[k]), "year": year, "cwd_mm": float(c.cwd_mm["pm_fao56"]), "vpd24_mean_kpa": float(np.mean(c.vpd_24h_kpa)),
                         "gdd_annual": float(c.gdd_cumulative[-1]), "t_mean_c_annual": float(np.mean(c.t_mean_c))})
        CHECKPOINT_PATH_DENSE.parent.mkdir(parents=True, exist_ok=True)
        CHECKPOINT_PATH_DENSE.write_text(json.dumps({"key": key, "rows": rows}))
        print(f"  {year}: checkpointed ({len(rows)} node-years)", flush=True)
    return pd.DataFrame(rows)


def main():
    if "--dense" not in sys.argv:
        sys.exit("MNEME's forest panel is defined on the dense grid only: pass --dense")
    if not PIXELS_PATH.exists():
        sys.exit(f"{PIXELS_PATH.name} not found: run scripts/extract_mneme_forest_pixels.py --dense first")
    static = extract_static_grid_inputs(grid_rows=DENSE_GRID_ROWS, grid_cols=DENSE_GRID_COLS)
    z = np.load(PIXELS_PATH, allow_pickle=False)
    pix = {k: z[k] for k in z.files}
    years = [int(y) for y in pix["years"]]
    assert years == VITALITY_YEARS, "the pixel file was built for a different set of years"
    problem = fp.disturbance_layer_problem(pix["no_disturbance"], years)
    if problem:
        sys.exit(f"STOP: {problem}")
    nodes = np.unique(pix["node"])
    n_pix = len(pix["node"])
    print(f"=== {n_pix} forest pixels in {len(nodes)} climate cells; {int(np.isfinite(pix['kndvi']).all(axis=1).sum())} have a complete 2000-2024 kNDVI series ===", flush=True)

    print(f"=== Climate of the {len(nodes)} cells, {PANEL_YEARS[0]}-{PANEL_YEARS[-1]} (TOPOHYDRO forcing, per year) ===", flush=True)
    climate = node_climate(static, nodes)

    panel = fp.build_person_years(pix, years, PANEL_YEARS, climate, CLIMATE_COLS)
    PANEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(PANEL_PATH, index=False)
    risk = panel[panel["at_risk"]]
    n_py, n_ev = len(risk), int(risk["event"].sum())
    n_cand = int(panel[panel["candidate_at_risk"]]["candidate"].sum())
    rate = n_ev / n_py if n_py else float("nan")
    lo, hi = event_rate_interval(n_ev, n_py) if n_py else (float("nan"), float("nan"))
    print(f"=== Panel: {len(panel)} pixel-years ({len(risk)} at risk), {panel['pixel'].nunique()} pixels, {panel['node'].nunique()} cells; "
          f"{n_ev} dieback onsets (strict rule), {n_cand} if harvest/fire flags were ignored ===", flush=True)

    full = np.isfinite(pix["kndvi"]).all(axis=1)
    onsets_down = fp.count_onset_pixels(pix["kndvi"][full], years, PANEL_YEARS)
    onsets_up = fp.mirror_onset_count(pix["kndvi"][full], years, PANEL_YEARS)
    adj = fp.remove_year_effect(pix["kndvi"])
    onsets_down_adj = fp.count_onset_pixels(adj[full], years, PANEL_YEARS)
    onsets_up_adj = fp.mirror_onset_count(adj[full], years, PANEL_YEARS)
    risk_adj = panel[panel["at_risk_adj"]]
    n_ev_adj, n_py_adj = int(risk_adj["event_adj"].sum()), len(risk_adj)
    lo_adj, hi_adj = event_rate_interval(n_ev_adj, n_py_adj) if n_py_adj else (float("nan"), float("nan"))
    print(f"=== The rule on the vitality series alone, flags ignored: {onsets_down} pixels with an onset of decline, {onsets_up} with an onset of rise (mirrored series); "
          f"with the common year effect removed: {onsets_down_adj} and {onsets_up_adj}. Strict events with the year effect removed: {n_ev_adj} ===", flush=True)

    X_cols = [f"{c}_dev" for c in CLIMATE_COLS] + [f"{c}_bar" for c in CLIMATE_COLS] + ["height_m"]
    result = {
        "run_date": datetime.date.today().isoformat(), "grid": "dense_armenia_stride7", "design": "forest pixels (closed forest on the Ecosystem Map), " + ("every qualifying pixel of each climate cell" if int(pix.get("pixels_per_node", 25)) >= 1000 else f"up to {int(pix.get('pixels_per_node', 25))} per climate cell"),
        "pixels_per_cell_cap": int(pix.get("pixels_per_node", 25)), "n_qualifying_pixels": int(pix["n_qualifying_pixels"]) if "n_qualifying_pixels" in pix else None,
        "panel_years": PANEL_YEARS, "n_points": int(len(nodes)), "n_pixels": int(panel["pixel"].nunique()), "n_points_valid_kndvi": int(panel["pixel"].nunique()),
        "n_cells": int(panel["node"].nunique()), "n_person_years": int(n_py), "n_pixel_years_total": int(len(panel)),
        "n_events": n_ev, "n_events_if_flags_ignored": n_cand, "n_pixels_with_event": int(risk.loc[risk["event"], "pixel"].nunique()),
        "n_cells_with_event": int(risk.loc[risk["event"], "node"].nunique()),
        "event_rate_per_pixel_year": None if not n_py else float(rate),
        "event_rate_ci95": None if not n_py else [float(lo), float(hi)],
        "event_rate_upper_bound_95_one_sided": float(zero_event_upper_bound(n_py)) if n_py and n_ev == 0 else None,
        "year_effect_removed": {"n_events": n_ev_adj, "n_person_years": int(n_py_adj), "event_rate_per_pixel_year": None if not n_py_adj else float(n_ev_adj / n_py_adj),
                                "event_rate_ci95": None if not n_py_adj else [float(lo_adj), float(hi_adj)],
                                "note": "the strict rule on series from which the year effect common to the forest pixels was removed (sensor and composite changes step all pixels together); a country-wide drought is removed with it, so this is a check next to the unadjusted rule, not a replacement"},
        "mirror_check": {"n_pixels_with_full_series": int(full.sum()), "pixels_with_onset_of_decline": onsets_down, "pixels_with_onset_of_rise": onsets_up,
                         "year_effect_removed": {"pixels_with_onset_of_decline": onsets_down_adj, "pixels_with_onset_of_rise": onsets_up_adj},
                         "note": "the dieback rule on the vitality series alone (harvest/fire flags ignored), and on the mirrored series, which counts persistent rises; clearly more declines than rises would be the signature of real dieback. In the real record rises outnumber declines: the unadjusted series carries a sensor-driven upward step, so the rule is not read as a dieback detector without the year-effect-removed variant"},
        "forest_class_pixels": {int(k): int(v) for k, v in panel.drop_duplicates("pixel")["forest_class"].value_counts().sort_index().items()},
        "design_columns": X_cols, "min_events_per_predictor": MIN_EVENTS_PER_PREDICTOR,
        "rule": "z < -2 against the pixel's own previous 10 years, starting a decline (previous year not below -2), no harvest or fire flag in the year or the year before, and no recovery above -1 in the next 2 years; first onset only",
        "scope_note": ("Climate-only Mundlak design (cwd_mm, vpd24_mean_kpa, gdd_annual, t_mean_c_annual, each split within/between by climate cell) plus the 2019 canopy height of the pixel; "
                        "no stand age, lagged-drought cross-basis or stand-interaction terms, because no such data exist for Armenian plots. Misclassification correction (Se/Sp) is not applied: "
                        "no labelled sample exists. Pixels of one cell share their climate, so person-years are not independent and the event-rate interval is narrower than the data justify."),
    }

    need = MIN_EVENTS_PER_PREDICTOR * len(X_cols)
    ratio = onsets_down / onsets_up if onsets_up else float("inf")
    result["decline_to_rise_ratio"] = None if ratio == float("inf") else float(ratio)
    result["min_decline_to_rise_ratio"] = MIN_DECLINE_TO_RISE_RATIO
    if n_ev < need:
        result["status"] = "too_few_events_for_a_hazard_fit"
        result["status_detail"] = f"{n_ev} events; a stable fit of {len(X_cols)} columns needs at least {need} ({MIN_EVENTS_PER_PREDICTOR} per column)"
        print(f"=== {n_ev} events < {need} needed: reporting the panel and the event rate only ===", flush=True)
    elif ratio < MIN_DECLINE_TO_RISE_RATIO:
        result["status"] = "events_not_distinguishable_from_noise"
        result["status_detail"] = (f"{n_ev} events, enough for the fit, but the same rule finds {onsets_down} pixels with a decline and {onsets_up} with a rise on the mirrored series (ratio {ratio:.2f}; at least "
                                   f"{MIN_DECLINE_TO_RISE_RATIO:g} is required): the events cannot be told from variation in the record, so no hazard is fitted")
        print(f"=== {n_ev} events but declines/rises = {onsets_down}/{onsets_up}: not distinguishable from noise, no hazard fitted ===", flush=True)
    else:
        dec = mundlak_decompose(risk.reset_index(drop=True), "node", CLIMATE_COLS)
        Xm = dec[X_cols].to_numpy()
        ok = np.all(np.isfinite(Xm), axis=1)
        dec, Xm = dec[ok], Xm[ok]
        y = dec["event"].to_numpy().astype(float)
        blocks = spatial_block_ids(dec["lat"].to_numpy(), dec["lon"].to_numpy(), SPATIAL_BLOCK_SIZE_DEG)
        n_blocks = len(np.unique(blocks))
        result["n_spatial_blocks"] = int(n_blocks)
        if n_blocks < 2:
            result["status"] = "fewer_than_two_spatial_blocks"
        else:
            fit = fit_stacked_hazard(Xm, y, blocks, glm_kwargs={"l2": 1e-3}, seed=0)
            valid = ~np.isnan(fit["oof_pred_stat"])
            result["status"] = "fitted"
            result["glm_coefficients"] = dict(zip(X_cols, [float(v) for v in fit["glm"].coef_]))
            result["glm_intercept"] = float(fit["glm"].intercept_)
            result["stack_weights_glm_gbm"] = [float(w) for w in fit["stack_weights"]]
            result["n_oof_valid"] = int(valid.sum())
            if valid.sum() > 0 and len(np.unique(y[valid])) > 1:
                result["oof_auc_stacked"] = float(roc_auc_score(y[valid], fit["oof_pred_stat"][valid]))
                result["oof_brier_stacked"] = float(brier_score_loss(y[valid], fit["oof_pred_stat"][valid]))
                result["oof_auc_rf_benchmark"] = float(roc_auc_score(y[valid], fit["oof_pred_rf"][valid]))
            print(f"=== Fitted. OOF AUC (stacked): {result.get('oof_auc_stacked', 'n/a')} ===", flush=True)

    OUT_PATH_DENSE.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH_DENSE.write_text(yaml.dump(result, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH_DENSE} and {PANEL_PATH.name} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
