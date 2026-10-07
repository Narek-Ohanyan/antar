"""Real treeline CHANGE across the ISIMIP3b scenario ensemble -- the user's explicit requirement
(2026-09-30: "make sure that we model the Treeline change and have it on the UI"). The earlier
`compute_treeline_diagnostic.py` produced only today's (2019) climatic ceiling; this extends it to
every real (GCM x SSP x horizon) member of the same 45-member ensemble future-projections uses.

Method -- the same real delta/change-factor construction as `run_future_projections.py`, applied
only to the one quantity treeline needs, so no Monte Carlo and no soil/ERA5 streaming is involved:

1. Reproduce the pipeline's own growing-season temperature exactly. `topoclimate_forcing` computes
   t_mean_c = downscale_temperature(t_mean_ref, z_cell, z_ref, gamma_of_day) and then
   GST = mean(t_mean_c[t_mean_c >= 0.9]) (`indices.growing_season_mean_temperature`). Both are
   called here directly, on the same real CHELSA-daily 2019 reference series, with the same real
   fitted monthly lapse rate -- not a re-derivation. Only terrain (elevation, z_ref_m) is needed
   (read from the cached static inputs; streamed from Drive only the first time), because those two numbers are the only non-CHELSA inputs to that formula.
2. Baseline check, stated not assumed: the reproduced 2019 GST must match the value the full
   pipeline stored in `treeline_diagnostic_2019.yaml` for the same cell. The maximum absolute
   difference is recorded in the output and the run aborts if it exceeds 0.05 degC.
3. For each member, add the real monthly ISIMIP3b temperature delta (future 5-year-window mean minus
   2015-2019 mean, from `run_future_projections.compute_deltas`) to the daily reference series, re-run
   the same two steps, and invert with `antar.niche.growth.potential_treeline_elevation`.

GST is recomputed from the shifted daily series rather than approximated as GST_2019 + mean monthly
delta, because warming also changes *which days* clear the 0.9 degC growing-season threshold (the
season lengthens) -- that second effect is real and the shortcut would miss it.

Stated limits: this is a CLIMATIC CEILING ("realised treelines lag their climatic potential by
decades", Sec. 8.1), not a forecast of where forest will stand; the thermal threshold (6.45 degC) is
the global Koerner-Paulsen default, not an Armenia-specific calibration; one representative year per
horizon, with 2019's real weather sequence perturbed by the monthly delta.

`--dense` runs the Armenia-only dense grid and writes a separate file.
"""
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import datetime  # noqa: E402

from run_topohydro_grid import (  # noqa: E402
    grid_latlon, extract_static_grid_inputs, _load_chelsa_arrays,
    GRID_ROWS, GRID_COLS, DENSE_GRID_ROWS, DENSE_GRID_COLS,
)
from run_future_projections import compute_deltas, GCMS, SCENARIOS, HORIZONS  # noqa: E402

from antar.climate import downscale, indices  # noqa: E402
from antar.niche.growth import potential_treeline_elevation  # noqa: E402

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"
OUT_PATH = CONFIG_DIR / "fitted" / "treeline_change.yaml"
OUT_PATH_DENSE = CONFIG_DIR / "fitted" / "treeline_change_dense.yaml"
GROWING_SEASON_MONTHS = ["Apr", "May", "Jun", "Jul", "Aug", "Sep"]
MONTH_ORDER = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
GROWING_SEASON_T0_C = 0.9   # matches topoclimate_forcing's own default (growing_season_t0_c)
BASELINE_TOLERANCE_C = 0.05


def growing_season_temperature(t_ref_daily, month, z_cell_m, z_ref_m, gamma_k_per_m_by_month, delta_by_month=None):
    """GST for one cell: the pipeline's own two-step computation, with an optional additive
    monthly delta applied to the reference series first. `month` is 1..12 per day."""
    t_ref = np.asarray(t_ref_daily, dtype=float)
    if delta_by_month is not None:
        t_ref = t_ref + np.asarray(delta_by_month, dtype=float)[month - 1]
    gamma_of_day = np.asarray(gamma_k_per_m_by_month, dtype=float)[month - 1]
    t_mean_c = downscale.downscale_temperature(t_ref, z_cell_m, z_ref_m, gamma_of_day)
    return indices.growing_season_mean_temperature(t_mean_c, t0_c=GROWING_SEASON_T0_C)


def main():
    dense = "--dense" in sys.argv
    grid_rows, grid_cols = (DENSE_GRID_ROWS, DENSE_GRID_COLS) if dense else (GRID_ROWS, GRID_COLS)
    out_path = OUT_PATH_DENSE if dense else OUT_PATH

    lapse = yaml.safe_load(open(CONFIG_DIR / "fitted" / "topohydro_lapse_rate.yaml"))
    gamma_by_month = np.array([lapse["temperature"]["gamma_k_per_m_by_month"][m] for m in MONTH_ORDER])
    gamma_gs = float(np.mean([lapse["temperature"]["gamma_k_per_m_by_month"][m] for m in GROWING_SEASON_MONTHS]))

    lats, lons, chelsa_row, chelsa_col = grid_latlon(grid_rows, grid_cols)
    n = len(lats)
    print(f"=== {n} grid points ({'DENSE' if dense else 'validation'}); terrain from the cached static inputs ===", flush=True)
    static = extract_static_grid_inputs(grid_rows, grid_cols)         # streamed from Drive once and cached (data/_cache/static_*.npz); no Drive access when the cache exists
    elevation, z_ref_m = static["elevation"], static["z_ref_m"]
    ok = ~np.isnan(elevation)
    print(f"  {int(ok.sum())}/{n} points have real terrain", flush=True)

    arrs = _load_chelsa_arrays()
    day0 = datetime.date(1979, 1, 1)
    i0 = (datetime.date(2019, 1, 1) - day0).days
    i1 = (datetime.date(2019, 12, 31) - day0).days + 1
    month = np.array([(day0 + datetime.timedelta(days=d)).month for d in range(i0, i1)])
    t_ref_all = arrs["tas"][i0:i1, chelsa_row, chelsa_col] - 273.15  # (days, n)

    gst_2019 = np.full(n, np.nan)
    for i in np.where(ok)[0]:
        gst_2019[i] = growing_season_temperature(t_ref_all[:, i], month, elevation[i], z_ref_m[i], gamma_by_month)
    ztl_2019 = potential_treeline_elevation(gst_2019, elevation, gamma_gs)

    # Baseline check against the full pipeline's own stored value (matched by lat/lon).
    diag = yaml.safe_load(open(CONFIG_DIR / "fitted" / "treeline_diagnostic_2019.yaml"))
    stored = {(round(c["lat"], 6), round(c["lon"], 6)): c["growing_season_mean_t_c"] for c in diag["cells"]}
    diffs = [abs(gst_2019[i] - stored[(round(float(lats[i]), 6), round(float(lons[i]), 6))])
             for i in np.where(ok)[0] if (round(float(lats[i]), 6), round(float(lons[i]), 6)) in stored]
    if diffs:
        max_diff = float(np.max(diffs))
        print(f"=== Baseline check: reproduced 2019 GST vs full-pipeline stored value, {len(diffs)} "
              f"matched cells, max |diff| = {max_diff:.4f} degC ===", flush=True)
        if max_diff > BASELINE_TOLERANCE_C:
            raise SystemExit(f"baseline GST mismatch {max_diff:.3f} > {BASELINE_TOLERANCE_C} -- refusing to "
                             "report treeline change from a reconstruction that does not match the pipeline")
    else:
        max_diff = None
        print("=== Baseline check skipped: no cell overlaps the stored 78-point diagnostic "
              "(expected for --dense, whose cells differ) ===", flush=True)

    cells = [{"lat": float(lats[i]), "lon": float(lons[i]), "elevation_m": float(elevation[i]),
              "gst_2019_c": float(gst_2019[i]), "potential_treeline_2019_m": float(ztl_2019[i])}
             for i in np.where(ok)[0]]
    keep = np.where(ok)[0]

    members = {}
    for gcm in GCMS:
        for scenario in SCENARIOS:
            deltas = compute_deltas(lats, lons, gcm, scenario)
            for horizon in HORIZONS:
                dtas = deltas[(horizon, "tas")]  # (12 months, n)
                ztl = np.full(n, np.nan)
                gst = np.full(n, np.nan)
                for i in keep:
                    gst[i] = growing_season_temperature(t_ref_all[:, i], month, elevation[i], z_ref_m[i],
                                                        gamma_by_month, delta_by_month=dtas[:, i])
                ztl = potential_treeline_elevation(gst, elevation, gamma_gs)
                shift = ztl - ztl_2019
                members[f"{gcm}__{scenario}__{horizon}"] = {
                    "gcm": gcm, "scenario": scenario, "horizon": int(horizon),
                    "gst_c": [round(float(gst[i]), 3) for i in keep],
                    "potential_treeline_m": [round(float(ztl[i]), 1) for i in keep],
                    "shift_m": [round(float(shift[i]), 1) for i in keep],
                    "mean_shift_m": round(float(np.nanmean(shift[keep])), 1),
                }
            print(f"  {gcm} / {scenario}: done", flush=True)

    summary = {}
    for scenario in SCENARIOS:
        for horizon in HORIZONS:
            vals = [m["mean_shift_m"] for m in members.values()
                    if m["scenario"] == scenario and m["horizon"] == int(horizon)]
            summary[f"{scenario}__{horizon}"] = {
                "ensemble_mean_shift_m": round(float(np.mean(vals)), 1),
                "ensemble_min_shift_m": round(float(np.min(vals)), 1),
                "ensemble_max_shift_m": round(float(np.max(vals)), 1),
                "n_gcms": len(vals),
            }
            print(f"  {scenario} {horizon}: mean treeline shift {summary[f'{scenario}__{horizon}']['ensemble_mean_shift_m']:+.0f} m "
                  f"(GCM range {summary[f'{scenario}__{horizon}']['ensemble_min_shift_m']:+.0f} to "
                  f"{summary[f'{scenario}__{horizon}']['ensemble_max_shift_m']:+.0f} m)", flush=True)

    result = {
        "run_date": datetime.date.today().isoformat(),
        "grid": "dense_armenia_stride7" if dense else "validation_80pt_stride40",
        "thermal_threshold_c": 6.45,
        "gamma_growing_season_k_per_km": gamma_gs * 1000,
        "growing_season_t0_c": GROWING_SEASON_T0_C,
        "baseline_check_max_abs_diff_c": max_diff,
        "method": ("Real monthly ISIMIP3b temperature delta (5-year-window future minus 2015-2019) added "
                   "to the real 2019 CHELSA-daily reference series; growing-season temperature recomputed "
                   "from the shifted daily series with the pipeline's own functions (so the lengthening "
                   "season is captured, not approximated); inverted with potential_treeline_elevation."),
        "scope_note": ("Climatic ceiling only, not a forecast of realised forest position (realised treelines "
                       "lag by decades). Global Koerner-Paulsen 6.45 degC threshold, not an Armenia "
                       "calibration. One representative year per horizon."),
        "n_cells": len(cells),
        "cells": cells,
        "summary_by_scenario_horizon": summary,
        "members": members,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.dump(result, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {out_path} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
