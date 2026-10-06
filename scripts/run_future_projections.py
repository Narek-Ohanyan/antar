"""Real future climate-change projections for TOPOHYDRO/XYLEM/REFUGIUM, the real next step
unlocked by ISIMIP3b's completion this session (60/60, all 5 real GCMs x 3 real SSPs).

Design: the real "delta method" / change-factor downscaling (standard in the climate-impacts
literature, not invented for this project) -- NOT a fresh elevation-lapse downscaling of
ISIMIP3b's own coarse values. Reasoning, stated because it isn't obvious: ISIMIP3b's native
resolution is 0.5 degrees (~50km, confirmed against the real data shape (5,6) over this
project's bbox, matching the standard global ISIMIP3b/W5E5 grid cell-center convention,
-89.75..89.75 / -179.75..179.75), far coarser than CHELSA-daily's ~1km. Re-deriving a "coarse
reference elevation" (the way z_ref_m already works for CHELSA, S:func:`compute_forcing_for_year`)
at ISIMIP3b's ~50km footprint would need a window roughly 50x wider than the one already used
for CHELSA, with no real way to verify the right answer the way the CHELSA z_ref_m choice could
be sanity-checked against station data. The delta method sidesteps this entirely: it only ever
asks ISIMIP3b for a smooth *monthly climate-change signal* (future minus 2019, additive for
temperature, multiplicative for precipitation, matching the standard convention that avoids
negative precipitation), applied on top of the REAL 2019 CHELSA-daily reference series --
preserving 2019's real day-to-day weather sequence/variability, perturbed by a real projected
shift, then run through the exact same already-tested topoclimate_forcing pipeline (same real
lapse rate, same real z_ref_m, same real cold-air-pooling, same real soil bucket) unchanged.

Two real, stated scope reductions, consistent with every other real run this session:

* **Single representative year per horizon (2050, 2080, 2100), not the full 20-year
  climatological window** `configs/scenarios.yaml`'s own convention specifies. The monthly
  delta itself IS averaged over a real 5-year window around each target year (2048-2052,
  2078-2082, 2096-2100 -- the last one real-data-limited to end at 2100, not extending past it)
  to damp single-year weather noise out of the climate-change signal; but the perturbed "future"
  daily series applied is still one representative year (2019's real shape), not a full
  20-year-mean climatology. The full 20-year version is real, substantial follow-up work, not
  done here.
* **Wind/radiation/dewpoint/pressure (ERA5-Land-sourced) held at their real 2019 values.** No
  real future projection for these exists in anything pulled this session -- ERA5-Land is a
  reanalysis product (historical observations only). Only temperature and precipitation, which
  DO have a real ISIMIP3b future projection, actually change between the historical run and
  these future scenarios -- stated explicitly, not silently assumed unchanged.

Runs all 15 real GCM x SSP members x 3 horizons x 3 real XYLEM functional groups (juniper
skipped, same precedent as every other real run -- n=1 XFT record, no real variance estimate)
through XYLEM's real two-level Monte Carlo and REFUGIUM's real viability/robust-refugium.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import datetime
import json
import multiprocessing
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    DATA_DIR, extract_static_grid_inputs, extract_era5land, _load_chelsa_arrays,
    saturation_vapour_pressure, wind_speed_2m, ROOTING_DEPTH_MM_BY_GROUP,
    GRID_ROWS, GRID_COLS, DENSE_GRID_ROWS, DENSE_GRID_COLS, ROOTING_DEPTH_MM_PLACEHOLDER, constant_atmosphere_requested, atmosphere_label,
)
from fit_xylem_mechanistic_hazard import load_functional_groups, PET_FORMULATION, OUTER_DRAWS, INNER_DRAWS  # noqa: E402
from fit_refugium_viability import V_STAR, RHO, LAM  # noqa: E402

from antar.climate import downscale  # noqa: E402
from antar.climate.atmosphere import AtmosphereShape  # noqa: E402
from antar.climate.forcing import topoclimate_forcing  # noqa: E402
from antar.climate.radiation import net_radiation_from_era5  # noqa: E402
from antar.hydraulics.monte_carlo import two_level_failure_probability  # noqa: E402
from antar.hydraulics.twophase import simulate_two_phase  # noqa: E402
from antar.viability.cohort import viability  # noqa: E402
from antar.viability.refugia import refugium_score, robust_refugium  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "future_projections.yaml"
CHECKPOINT_PATH = DATA_DIR / "_future_projections_checkpoint.json"
# Real 2026-10-02 densification, explicitly requested (user: "re-run even if it will take a day
# or two... I want real and solid outputs"): --dense runs the full real 45-member ensemble against
# the Armenia-only DENSE_GRID instead of 78, with its own output/checkpoint paths. A real, large
# cost -- 45 members x ~925 Armenian cells x 3 groups x the same real Monte Carlo that already ran 45x78x3
# times before -- deliberately accepted, not cut down, per that explicit instruction.
OUT_PATH_DENSE = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "future_projections_dense.yaml"
CHECKPOINT_PATH_DENSE = DATA_DIR / "_future_projections_checkpoint_dense.json"

ISIMIP_DIR = DATA_DIR / "isimip3b"
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]
BASELINE_YEAR = 2019  # matches the real CHELSA reference year every other run this session used
BASELINE_WINDOW = (2015, 2019)  # real ISIMIP3b data starts 2015 -- earliest real window available
HORIZONS = {2050: (2048, 2052), 2080: (2078, 2082), 2100: (2096, 2100)}  # 5yr windows; 2100's is
                                                                          # real-data-limited, not centred

# Standard global ISIMIP3b/W5E5 0.5deg grid cell-centre convention -- verified this session
# against the real data shape (5,6) over this project's bbox before relying on it.
ISIMIP_LATS = np.arange(-89.75, 90, 0.5)
ISIMIP_LONS = np.arange(-179.75, 180, 0.5)
BBOX = (43.4, 38.8, 46.7, 41.4)


def isimip_cell_indices(lats, lons):
    """Nearest ISIMIP3b 0.5deg cell (row, col) for each real grid point."""
    lat_sel = ISIMIP_LATS[(ISIMIP_LATS >= BBOX[1]) & (ISIMIP_LATS <= BBOX[3])]
    lon_sel = ISIMIP_LONS[(ISIMIP_LONS >= BBOX[0]) & (ISIMIP_LONS <= BBOX[2])]
    row = np.array([np.argmin(np.abs(lat_sel - lat)) for lat in lats])
    col = np.array([np.argmin(np.abs(lon_sel - lon)) for lon in lons])
    return row, col


def load_isimip_var(var, gcm, scenario):
    path = ISIMIP_DIR / f"CHELSA_ISIMIP3b_{var}__{gcm}__{scenario}.npz"
    npz = np.load(path)
    dates = npz["dates"]
    years = np.array([int(d[:4]) for d in dates])
    months = np.array([int(d[5:7]) for d in dates])
    return npz["data"], years, months


def monthly_means_for_window(data, years, months, row, col, year_start, year_end, n_points):
    """Real monthly mean at each point's nearest ISIMIP3b cell, averaged over [year_start, year_end]."""
    sel = (years >= year_start) & (years <= year_end)
    out = np.zeros((12, n_points))
    for m in range(1, 13):
        msel = sel & (months == m)
        for i in range(n_points):
            out[m - 1, i] = np.nanmean(data[msel, row[i], col[i]])
    return out


def compute_deltas(lats, lons, gcm, scenario):
    """Real monthly delta (temperature, additive K) and ratio (precip, multiplicative) per
    point per horizon, from real local ISIMIP3b data -- no network access needed.
    """
    row, col = isimip_cell_indices(lats, lons)
    n = len(lats)
    deltas = {}
    for var in ["tas", "tasmax", "tasmin", "pr"]:
        data, years, months = load_isimip_var(var, gcm, scenario)
        baseline = monthly_means_for_window(data, years, months, row, col, *BASELINE_WINDOW, n)
        for horizon, (y0, y1) in HORIZONS.items():
            future = monthly_means_for_window(data, years, months, row, col, y0, y1, n)
            key = (horizon, var)
            if var == "pr":
                deltas[key] = np.where(baseline > 0.01, future / np.maximum(baseline, 0.01), 1.0)
            else:
                deltas[key] = future - baseline  # already Kelvin, so additive delta is degC-equivalent
    return deltas


def build_future_cell(static, deltas, horizon, i, elevation, z_ref_m, slope, aspect, concavity,
                       soil, w_max_mm, lapse, lats, era5_2019, chelsa_ref, atmos=None):
    """``atmos``: daily multiplicative factors {rs, rl, wind, ea} for this cell and scenario (antar.climate.atmosphere); None =
    the old behaviour (ERA5-Land annual means repeated every day, unchanged across scenarios)."""
    t_mean_ref_c, t_max_ref_c, t_min_ref_c, p_ref_mm, doy, month = chelsa_ref
    # deltas[(horizon, var)] is shaped (12 months, n_points) -- index the point axis with [:, i],
    # not [i] (which would index into the month axis and crash once i >= 12).
    d_tas = deltas[(horizon, "tas")][:, i]
    d_tasmax = deltas[(horizon, "tasmax")][:, i]
    d_tasmin = deltas[(horizon, "tasmin")][:, i]
    r_pr = deltas[(horizon, "pr")][:, i]

    month_idx = month - 1
    future_t_mean = t_mean_ref_c + d_tas[month_idx]
    future_t_max = t_max_ref_c + d_tasmax[month_idx]
    future_t_min = t_min_ref_c + d_tasmin[month_idx]
    future_p = p_ref_mm * r_pr[month_idx]

    gamma_k_per_m, precip_gradient_per_m = lapse
    gamma_of_day = gamma_k_per_m[month - 1]
    t_mean_c_for_rn = downscale.downscale_temperature(future_t_mean, elevation[i], z_ref_m[i], gamma_of_day)
    wind10_i, ssrd_i, strd_i, dewpoint_k_i, pressure_pa_i = era5_2019
    n_days = len(doy)
    f = atmos if atmos is not None else AtmosphereShape.constant_factors(n_days)
    rn_mj_m2 = net_radiation_from_era5(ssrd_i[i] * f["rs"], strd_i[i] * f["rl"], t_mean_c_for_rn)
    u2_m_s_i = wind_speed_2m(wind10_i[i], z_m=10.0)
    ea_ref_kpa_i = saturation_vapour_pressure(dewpoint_k_i[i] - 273.15)

    return topoclimate_forcing(
        doy=doy, month=month,
        t_mean_ref_c=future_t_mean, t_max_ref_c=future_t_max, t_min_ref_c=future_t_min,
        p_ref_mm=future_p, ea_ref_kpa=ea_ref_kpa_i * f["ea"],
        u2_m_s=u2_m_s_i * f["wind"], rn_mj_m2=rn_mj_m2,
        z_cell_m=elevation[i], z_ref_m=z_ref_m[i], lat_deg=lats[i],
        slope_deg=slope[i], aspect_deg=aspect[i],
        gamma_k_per_m=gamma_k_per_m, precip_gradient_per_m=precip_gradient_per_m,
        w_max_mm=w_max_mm[i], theta_sat=soil["theta_sat"][i], psi_sat_mpa=soil["psi_sat_mpa"][i],
        b_clapp_hornberger=soil["b_clapp_hornberger"][i], theta_fc=soil["theta_fc"][i],
        theta_lim=soil["theta_lim"][i], gdd_budburst=200.0,
        concavity_index=concavity[i], calm_clear_night_frac=0.3,
        pressure_kpa=pressure_pa_i[i] / 1000.0,
    )


# ---------------------------------------------------------------------------------------------
# One (GCM, SSP, horizon) member: every valid cell x every real group. Independent of every other member
# (cell seeds are the cell's grid index, nothing is shared), so members run in parallel worker processes
# and give bit-identical numbers to a serial run.
# ---------------------------------------------------------------------------------------------
SCHEMA = 2     # bumped when the per-member record gains fields; older checkpoint files are ignored
_CTX = {}


def _r(x, nd=6):
    return None if x is None else round(float(x), nd)


def run_member_horizon(c, deltas, gcm, scenario, horizon):
    """c: the shared context dict (see main). Returns the member record."""
    lats, lons, n = c["lats"], c["lons"], c["n"]
    cell_results, climate_rows, generic_rows = {}, [], []
    for i in range(n):
        if not c["valid_mask"][i]:
            continue
        chelsa_ref_i = (c["t_mean_ref_all"][:, i], c["t_max_ref_all"][:, i], c["t_min_ref_all"][:, i],
                        c["p_ref_all"][:, i], c["doy"], c["month"])
        atmos = (c["atmos"].scenario_factors(lats[i], c["lons"][i], c["doy"], gcm, scenario, horizon)
                 if c.get("atmos") is not None else None)
        cell = None
        for gname, g in c["groups"].items():
            # each group's own w_max_mm (its own rooting depth) gives its own CWD / WSI / soil-potential signal
            cell = build_future_cell({}, deltas, horizon, i, c["elevation"], c["z_ref_m"], c["slope"], c["aspect"],
                                     c["concavity"], c["soil"], c["w_max_mm_by_group"][gname], c["lapse"], lats,
                                     c["era5_2019"], chelsa_ref_i, atmos)
            psi_soil = cell.psi_soil_mpa[PET_FORMULATION]

            def simulate(traits, psi_soil=psi_soil, cell=cell):
                sim = simulate_two_phase(psi_soil, cell.t_max_c, cell.vpd_24h_kpa, cell.pressure_kpa, traits, 1.0)
                return sim["hfi_max"]

            h = two_level_failure_probability(simulate, g["base"], g["hyper_sd"], g["individual_sd"],
                                              outer_draws=OUTER_DRAWS, inner_draws=INNER_DRAWS, seed=i)
            v = np.array([viability([hh], p_height_ok=1.0) for hh in h])
            p_viable = float(np.mean(v >= V_STAR))      # share of the 50 trait-knowledge draws with V >= V*
            cell_results.setdefault(gname, []).append({
                "lat": float(lats[i]), "lon": float(lons[i]),
                "viability_mean": float(v.mean()), "h_mech_mean": float(h.mean()),
                "viability_p10": _r(np.percentile(v, 10)), "viability_p90": _r(np.percentile(v, 90)),
                "p_viable": _r(p_viable), "robust_criterion_a": bool(p_viable >= RHO),
                "h_mech_sd": _r(h.std()),                                              # outer-loop (trait-knowledge) spread
                "refugium_score": _r(refugium_score(v[:, None], lam=LAM)[0]),         # median - lam * IQR across the 50 draws
                "cwd_mm": _r(cell.cwd_mm[PET_FORMULATION], 3), "wsi": _r(cell.wsi[PET_FORMULATION], 4),
                "psi_min_mpa": _r(np.min(psi_soil), 4),
            })
        # generic (species-agnostic, 1000 mm rooting depth) water balance: the scenario counterpart of the 2019
        # TOPOHYDRO layers, with all three PET formulations
        gcell = build_future_cell({}, deltas, horizon, i, c["elevation"], c["z_ref_m"], c["slope"], c["aspect"],
                                  c["concavity"], c["soil"], c["w_max_mm_generic"], c["lapse"], lats,
                                  c["era5_2019"], chelsa_ref_i, atmos)
        generic_rows.append({"lat": float(lats[i]), "lon": float(lons[i]),
                             "cwd_mm_by_pet": {k: _r(v, 3) for k, v in gcell.cwd_mm.items()},
                             "wsi": _r(gcell.wsi[PET_FORMULATION], 4),
                             "psi_min_mpa": _r(np.min(gcell.psi_soil_mpa[PET_FORMULATION]), 4)})
        # group-independent climate of the scenario year (identical for every group's cell)
        climate_rows.append({
            "lat": float(lats[i]), "lon": float(lons[i]),
            "t_mean_c": _r(np.mean(gcell.t_mean_c), 4), "precip_mm": _r(np.sum(gcell.p_mm), 2),
            "gdd": _r(gcell.gdd_cumulative[-1], 2), "late_frost_days": int(gcell.late_frost_days),
            "gsl_days": int(gcell.growing_season_length_days), "gst_c": _r(gcell.growing_season_mean_t_c, 4),
        })
    group_summary = {}
    for gname, rows in cell_results.items():
        group_summary[gname] = {
            "n_cells": len(rows),
            "mean_viability": float(np.mean([r["viability_mean"] for r in rows])),
            "mean_h_mech": float(np.mean([r["h_mech_mean"] for r in rows])),
            "mean_refugium_score": float(np.mean([r["refugium_score"] for r in rows])),
            "cells": rows,
        }
    return {"schema": SCHEMA, "gcm": gcm, "scenario": scenario, "horizon": horizon, "groups": group_summary,
            "climate": climate_rows, "generic": generic_rows}


def _init_worker(ctx_path):
    global _CTX
    with open(ctx_path, "rb") as f:
        _CTX = pickle.load(f)
    _CTX["delta_cache"] = {}
    _CTX["atmos"] = AtmosphereShape(_CTX["atmos_dir"]) if _CTX.get("atmos_dir") else None


def _run_task(task):
    gcm, scenario, horizon = task
    key = (gcm, scenario)
    if key not in _CTX["delta_cache"]:                       # local ISIMIP3b files, no network
        _CTX["delta_cache"][key] = compute_deltas(_CTX["lats"], _CTX["lons"], gcm, scenario)
    t0 = datetime.datetime.now()
    rec = run_member_horizon(_CTX, _CTX["delta_cache"][key], gcm, scenario, horizon)
    return f"{gcm}__{scenario}__{horizon}", rec, (datetime.datetime.now() - t0).total_seconds()


def _arg(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    return default


def main():
    dense = "--dense" in sys.argv
    grid_rows, grid_cols = (DENSE_GRID_ROWS, DENSE_GRID_COLS) if dense else (GRID_ROWS, GRID_COLS)
    out_path = OUT_PATH_DENSE if dense else OUT_PATH
    checkpoint_path = CHECKPOINT_PATH_DENSE if dense else CHECKPOINT_PATH

    print(f"=== Static grid inputs ({'DENSE Armenia-only (stride 7)' if dense else '80-point validation'} grid, "
          f"terrain/soil streamed once) ===", flush=True)
    static = extract_static_grid_inputs(grid_rows=grid_rows, grid_cols=grid_cols)
    lats, lons = static["lats"], static["lons"]
    elevation, z_ref_m = static["elevation"], static["z_ref_m"]
    slope, aspect, concavity = static["slope"], static["aspect"], static["concavity"]
    soil, valid_soil = static["soil"], static["valid_soil"]
    chelsa_row, chelsa_col = static["chelsa_row"], static["chelsa_col"]
    n = len(lats)

    print("=== Real 2019 ERA5-Land (held constant for all future scenarios -- see docstring) ===", flush=True)
    wind10, ssrd, strd, dewpoint_k, pressure_pa = extract_era5land(static["token"], static["row_px"], static["col_px"], year=2019)
    era5_2019 = (wind10, ssrd, strd, dewpoint_k, pressure_pa)
    valid_era5 = ~np.isnan(wind10)
    valid_mask = valid_soil & valid_era5

    print("=== Real 2019 CHELSA-daily reference series (the baseline every delta perturbs) ===", flush=True)
    arrs = _load_chelsa_arrays()
    day0 = datetime.date(1979, 1, 1)
    i0 = (datetime.date(2019, 1, 1) - day0).days
    i1 = (datetime.date(2019, 12, 31) - day0).days + 1
    dates = [day0 + datetime.timedelta(days=d) for d in range(i0, i1)]
    doy = np.array([d.timetuple().tm_yday for d in dates])
    month = np.array([d.month for d in dates])
    t_mean_ref_all = arrs["tas"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_max_ref_all = arrs["tasmax"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_min_ref_all = arrs["tasmin"][i0:i1, chelsa_row, chelsa_col] - 273.15
    p_ref_all = arrs["pr"][i0:i1, chelsa_row, chelsa_col]

    lapse_cfg = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_lapse_rate.yaml"))
    month_order = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    gamma_k_per_m = np.array([lapse_cfg["temperature"]["gamma_k_per_m_by_month"][m] for m in month_order])
    precip_gradient_per_m = np.array([lapse_cfg["precipitation"]["gradient_per_m_by_month"][m] for m in month_order])
    lapse = (gamma_k_per_m, precip_gradient_per_m)

    groups = load_functional_groups()
    real_groups = {k: v for k, v in groups.items() if v["status"] == "ok"}
    print(f"=== Real functional groups: {list(real_groups)} ===", flush=True)

    # Real per-group w_max_mm: cheap, local (no network) -- each real functional group's own
    # rooting depth (ROOTING_DEPTH_MM_BY_GROUP) now gives its own real water-holding capacity,
    # rather than every species sharing one generic value. See run_topohydro_grid.py's module
    # docstring for why this matters and why it's cheap to do per group here specifically.
    w_max_mm_by_group = {g: (soil["theta_fc"] - soil["theta_lim"]) * ROOTING_DEPTH_MM_BY_GROUP[g]
                          for g in real_groups}

    w_max_mm_generic = (soil["theta_fc"] - soil["theta_lim"]) * ROOTING_DEPTH_MM_PLACEHOLDER

    gcms = _arg("--gcms", ",".join(GCMS)).split(",")
    scenarios = _arg("--scenarios", ",".join(SCENARIOS)).split(",")
    workers = int(_arg("--workers", "6"))
    out_path = Path(_arg("--out", str(out_path)))
    ckpt_dir = DATA_DIR / f"_future_ckpt_{out_path.stem}"
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    tasks = [(g, s, h) for g in gcms for s in scenarios for h in HORIZONS]
    results = {}
    for g, s, h in tasks:
        f = ckpt_dir / f"{g}__{s}__{h}.json"
        if f.exists():
            rec = json.loads(f.read_text())
            if rec.get("schema") == SCHEMA:
                results[f"{g}__{s}__{h}"] = rec
    todo = [t for t in tasks if f"{t[0]}__{t[1]}__{t[2]}" not in results]
    print(f"=== {len(tasks)} members; {len(results)} already in {ckpt_dir.name}; {len(todo)} to run on {workers} worker(s) ===", flush=True)

    ctx = {"lats": lats, "lons": lons, "n": n, "valid_mask": valid_mask, "elevation": elevation, "z_ref_m": z_ref_m,
           "slope": slope, "aspect": aspect, "concavity": concavity, "soil": soil, "lapse": lapse,
           "era5_2019": era5_2019, "t_mean_ref_all": t_mean_ref_all, "t_max_ref_all": t_max_ref_all,
           "t_min_ref_all": t_min_ref_all, "p_ref_all": p_ref_all, "doy": doy, "month": month,
           "groups": real_groups, "w_max_mm_by_group": w_max_mm_by_group, "w_max_mm_generic": w_max_mm_generic,
           "atmos_dir": None if constant_atmosphere_requested() else str(ISIMIP_DIR)}

    def keep(mkey, rec, secs):
        results[mkey] = rec
        (ckpt_dir / f"{mkey}.json").write_text(json.dumps(rec))
        print(f"  {mkey}: done in {secs / 60:.1f} min ({len(results)}/{len(tasks)})", flush=True)

    if todo and workers > 1:
        ctx_path = DATA_DIR / f"_future_ctx_{out_path.stem}.pkl"
        with open(ctx_path, "wb") as f:
            pickle.dump(ctx, f, protocol=pickle.HIGHEST_PROTOCOL)
        try:
            with multiprocessing.get_context("spawn").Pool(workers, initializer=_init_worker, initargs=(str(ctx_path),)) as pool:
                for mkey, rec, secs in pool.imap_unordered(_run_task, todo):
                    keep(mkey, rec, secs)
        finally:
            ctx_path.unlink(missing_ok=True)
    elif todo:
        _init_ctx = dict(ctx, delta_cache={}, atmos=AtmosphereShape(ctx["atmos_dir"]) if ctx["atmos_dir"] else None)
        for g, s, h in todo:
            if (g, s) not in _init_ctx["delta_cache"]:
                _init_ctx["delta_cache"][(g, s)] = compute_deltas(lats, lons, g, s)
            t0 = datetime.datetime.now()
            keep(f"{g}__{s}__{h}", run_member_horizon(_init_ctx, _init_ctx["delta_cache"][(g, s)], g, s, h),
                 (datetime.datetime.now() - t0).total_seconds())

    missing = [t for t in tasks if f"{t[0]}__{t[1]}__{t[2]}" not in results]
    if missing:
        raise SystemExit(f"{len(missing)} members missing -- not writing a partial ensemble")

    final = {
        "run_date": datetime.date.today().isoformat(),
        "schema": SCHEMA,
        "baseline_year": BASELINE_YEAR,
        "baseline_window": list(BASELINE_WINDOW),
        "horizons": {str(k): list(v) for k, v in HORIZONS.items()},
        "atmosphere": atmosphere_label(),
        "method": "delta/change-factor downscaling: real monthly ISIMIP3b anomaly (additive "
                  "temperature, multiplicative precipitation) applied to the real 2019 CHELSA-daily "
                  "reference series, run through the unchanged real topoclimate_forcing pipeline, "
                  "once per real functional group using that group's own real rooting depth "
                  f"({ROOTING_DEPTH_MM_BY_GROUP}, Canadell et al. 1996) rather than one generic "
                  "value shared across all species, plus once with the generic 1000 mm rooting depth "
                  "(the scenario counterpart of the species-agnostic 2019 TOPOHYDRO layers). wind/"
                  "radiation/dewpoint/pressure held at real 2019 ERA5-Land values -- no real future "
                  "projection exists for these. Per cell and group: 50 trait-knowledge draws x 200 "
                  "individual draws; viability mean, P10/P90, P[V >= V*], robust-refugium criterion (a), "
                  "risk-averse score, outer-loop hazard spread.",
        "members": {k: results[k] for k in sorted(results)},
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Dumper = getattr(yaml, "CSafeDumper", yaml.SafeDumper)
    out_path.write_text(yaml.dump(final, Dumper=Dumper, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {out_path} ===", flush=True)
    for f in ckpt_dir.glob("*.json"):
        f.unlink()
    ckpt_dir.rmdir()


if __name__ == "__main__":
    sys.exit(main())
