"""MNEME's first real person-period hazard panel: real dieback labels (from the real
Sentinel/HLS/Landsat kNDVI vitality composites + Hansen/MODIS disturbance ancillary
already exported to Drive) crossed with real climate covariates (this session's
gridded TOPOHYDRO forcing, extended here across multiple years), fit with the
already-tested Sec. 7.3 stacked-learner pipeline (`antar.hazard.pipeline`).

Reuses the same 80-point grid this session's TOPOHYDRO/XYLEM runs validated (not
a new sample) -- one consistent real-data validation grid across engines for this
first end-to-end pass.

Two real, hard data-availability constraints set the panel's year range,
determined empirically rather than assumed:

* `dieback_event`'s trailing_window=10 and recovery_seasons=2 mean an event can
  only be *determined* (not merely observed) for years satisfying
  2000+10 <= year <= 2024-2, i.e. 2010-2022.
* CHELSA-daily's real `pr` coverage stops at 2019-12-31 (the source archive's own
  gap, confirmed earlier this session) -- TOPOHYDRO's water-balance-driven CWD
  covariate needs real precipitation, so climate covariates are only real through
  2019.

The overlap -- the only years with BOTH a real, determinable dieback label AND a
real climate covariate -- is 2010-2019 (10 years). That intersection, not either
constraint alone, sets PANEL_YEARS below.

One substantial, explicitly stated gap: Eq. 7.3's full design also includes an
age spline, a DLNM drought-legacy cross-basis and stand/size interaction terms.
No real stand-age, forest-structure or lagged-drought-response data has been
pulled for any Armenian plot -- `antar.hazard.pipeline.build_hazard_design`
requires a real, *varying* age array (`spline_basis` needs actual spread to
build a valid B-spline knot vector; a constant/fabricated age would not "do
nothing," it would corrupt the basis). Rather than invent one, this fit bypasses
`build_hazard_design` entirely and builds the design directly from
`antar.hazard.panel.mundlak_decompose` on the real climate covariates alone,
reusing `antar.hazard.pipeline.fit_stacked_hazard` (which only consumes an
already-built X/y/blocks, agnostic to how X was assembled) unmodified. This is a
real, reduced-scope hazard fit -- climate signal only -- not the full Eq. 7.3
model; the missing terms are a stated follow-up, not silently dropped.

Misclassification correction (Eq. 7.1, Se/Sp) is also skipped: no stratified
sample of interpreted pixel-years exists to estimate Se/Sp from. The fitted
hazard is therefore the raw satellite-kNDVI-based rate, not misclassification-
corrected -- stated, not hidden.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window
from sklearn.metrics import brier_score_loss, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    DATA_DIR, drive_vsicurl_url, extract_static_grid_inputs, compute_forcing_for_year,
)

from antar.hazard.observation import dieback_event  # noqa: E402
from antar.hazard.panel import mundlak_decompose  # noqa: E402
from antar.hazard.pipeline import fit_stacked_hazard  # noqa: E402
from antar.validation.splits import spatial_block_ids  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "mneme_hazard_panel_2010_2019.yaml"
CHECKPOINT_PATH = Path(__file__).resolve().parent.parent / "data" / "_mneme_panel_checkpoint.json"

VITALITY_TILE_SIZE_PX = 4864
VITALITY_TILE_IDS = {
    (0, 0): "1X5VF3X8k9DvJzpxVbO6Ccjfhq9tSvHLW", (0, 4864): "1iC4l1_UoG-2ZpAKP90CeM3e-1wbpZUYA",
    (4864, 0): "1CVhnSKyGkMuLQUgA4mKNoLNPvzZypZvj", (4864, 4864): "1JyNnH34ZEGSYcB4qMrthliOd35ThOSEP",
}
DISTURBANCE_FILE_ID = "1-m_xKLBqBZYY8EFhFb9C_9iCzcJcB-jv"
VITALITY_YEARS = list(range(2000, 2025))  # full real coverage, needed for dieback_event's trailing/recovery windows
PANEL_YEARS = list(range(2010, 2020))     # real climate-covariate x real-determinable-event overlap (see docstring)

SPATIAL_BLOCK_SIZE_DEG = 1.0  # coarse blocks over ~80 points spanning ~3.3x2.6 degrees


def extract_vitality_and_disturbance(token, row_px, col_px):
    n = len(row_px)
    kndvi = np.full((n, len(VITALITY_YEARS)), np.nan)
    tile_row = (row_px // VITALITY_TILE_SIZE_PX) * VITALITY_TILE_SIZE_PX
    tile_col = (col_px // VITALITY_TILE_SIZE_PX) * VITALITY_TILE_SIZE_PX
    for tile_key in sorted(set(zip(tile_row.tolist(), tile_col.tolist()))):
        file_id = VITALITY_TILE_IDS[tile_key]
        sel = (tile_row == tile_key[0]) & (tile_col == tile_key[1])
        url = drive_vsicurl_url(file_id)
        with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES"):
            with rasterio.open(url) as src:
                bidx = {b: j + 1 for j, b in enumerate(src.descriptions)}
                for i in np.where(sel)[0]:
                    lr, lc = int(row_px[i] - tile_key[0]), int(col_px[i] - tile_key[1])
                    w = Window(lc, lr, 1, 1)
                    for yi, y in enumerate(VITALITY_YEARS):
                        kndvi[i, yi] = src.read(bidx[f"kndvi_{y}"], window=w)[0, 0]
        print(f"  vitality tile {tile_key}: {sel.sum()} points", flush=True)

    no_disturbance = np.full((n, len(VITALITY_YEARS)), True)
    url = drive_vsicurl_url(DISTURBANCE_FILE_ID)
    with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES"):
        with rasterio.open(url) as src:
            bidx = {b: j + 1 for j, b in enumerate(src.descriptions)}
            for i in range(n):
                w = Window(int(col_px[i]), int(row_px[i]), 1, 1)
                for yi, y in enumerate(VITALITY_YEARS):
                    no_disturbance[i, yi] = bool(src.read(bidx[f"no_disturbance_{y}"], window=w)[0, 0])
    print(f"  disturbance ancillary: {n} points", flush=True)
    return kndvi, no_disturbance


def main():
    print("=== Static grid inputs (terrain/soil, streamed, reused across years) ===", flush=True)
    static = extract_static_grid_inputs()
    lats, lons = static["lats"], static["lons"]
    row_px, col_px, token = static["row_px"], static["col_px"], static["token"]
    n = len(lats)

    print("=== Real vitality (kNDVI) + disturbance, 2000-2024 (streamed) ===", flush=True)
    kndvi, no_disturbance = extract_vitality_and_disturbance(token, row_px, col_px)
    n_valid_kndvi = int(np.all(~np.isnan(kndvi), axis=1).sum())
    print(f"  {n_valid_kndvi}/{n} points have a fully valid (unmasked) 2000-2024 kNDVI series", flush=True)

    print("=== Real dieback-event labels (antar.hazard.observation.dieback_event) ===", flush=True)
    events = np.zeros((n, len(VITALITY_YEARS)), dtype=bool)
    for i in range(n):
        events[i] = dieback_event(kndvi[i], no_disturbance[i])
    n_events_total = int(events[:, [VITALITY_YEARS.index(y) for y in PANEL_YEARS]].sum())
    print(f"  {n_events_total} real dieback-onset events across {n} points x {len(PANEL_YEARS)} panel years", flush=True)

    print(f"=== Real climate covariates, {PANEL_YEARS[0]}-{PANEL_YEARS[-1]} (gridded TOPOHYDRO, per year) ===", flush=True)
    rows = []
    years_done = set()
    if CHECKPOINT_PATH.exists():
        rows = json.loads(CHECKPOINT_PATH.read_text())
        years_done = {r["year"] for r in rows}
        print(f"  resuming from checkpoint: {len(rows)} rows already computed for years {sorted(years_done)}",
              flush=True)

    for year in PANEL_YEARS:
        if year in years_done:
            print(f"  {year}: skipped, already in checkpoint", flush=True)
            continue
        cells = compute_forcing_for_year(static, year)
        yi = VITALITY_YEARS.index(year)
        for i in range(n):
            if cells[i] is None or np.isnan(kndvi[i, yi]):
                continue
            c = cells[i]
            rows.append({
                "place": i, "year": year, "lat": float(lats[i]), "lon": float(lons[i]),
                "event": bool(events[i, yi]),
                "cwd_mm": float(c.cwd_mm["pm_fao56"]),
                "vpd24_mean_kpa": float(np.mean(c.vpd_24h_kpa)),
                "gdd_annual": float(c.gdd_cumulative[-1]),
                "t_mean_c_annual": float(np.mean(c.t_mean_c)),
            })
        CHECKPOINT_PATH.parent.mkdir(parents=True, exist_ok=True)
        CHECKPOINT_PATH.write_text(json.dumps(rows))
        print(f"  {year}: checkpointed ({len(rows)} rows total)", flush=True)

    panel = pd.DataFrame(rows)
    print(f"=== Panel: {len(panel)} real person-year rows, {panel['event'].sum()} real events "
          f"({panel['place'].nunique()} places) ===", flush=True)

    climate_cols = ["cwd_mm", "vpd24_mean_kpa", "gdd_annual", "t_mean_c_annual"]
    panel_dec = mundlak_decompose(panel, "place", climate_cols)
    X_cols = [f"{c}_dev" for c in climate_cols] + [f"{c}_bar" for c in climate_cols]
    X = panel_dec[X_cols].to_numpy()
    y = panel_dec["event"].to_numpy().astype(float)
    blocks = spatial_block_ids(panel_dec["lat"].to_numpy(), panel_dec["lon"].to_numpy(), SPATIAL_BLOCK_SIZE_DEG)
    n_blocks = len(np.unique(blocks))
    print(f"=== Design: {X.shape}, {n_blocks} spatial blocks for CV ===", flush=True)

    result = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "panel_years": PANEL_YEARS,
        "n_points": int(n), "n_points_valid_kndvi": n_valid_kndvi,
        "n_person_years": int(len(panel)), "n_events": int(panel["event"].sum()),
        "n_spatial_blocks": int(n_blocks),
        "design_columns": X_cols,
        "scope_note": ("Climate-only Mundlak design (cwd_mm, vpd24_mean_kpa, gdd_annual, "
                        "t_mean_c_annual, each split within/between by place) -- no real "
                        "stand-age, DLNM drought-legacy or stand-structure data exists yet, "
                        "so Eq. 7.3's age spline/cross-basis/stand terms are not included "
                        "(see module docstring). Misclassification correction (Se/Sp) also "
                        "not applied -- no stratified ground-truth sample exists."),
    }

    if n_blocks < 2 or panel["event"].sum() < 2:
        result["status"] = "insufficient_real_events_or_blocks_for_a_meaningful_fit"
        print("=== Too few real events/blocks for a meaningful stacked fit -- reporting panel only ===", flush=True)
    else:
        fit = fit_stacked_hazard(X, y, blocks, glm_kwargs={"l2": 1e-3}, seed=0)
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
        print(f"=== GLM coefficients: {result['glm_coefficients']} ===", flush=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.dump(result, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH} ===", flush=True)

    if CHECKPOINT_PATH.exists():
        CHECKPOINT_PATH.unlink()


if __name__ == "__main__":
    sys.exit(main())
