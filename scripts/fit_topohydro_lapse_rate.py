"""TOPOHYDRO's first real fit: monthly temperature lapse rate and precipitation-
elevation gradient from real GHCN-Daily Armenian station observations
(configs/manifests/ghcnd_armenia.yaml), via antar.climate.downscale's existing
fit_monthly_lapse_rate / fit_precip_elevation_gradient -- reused directly, no
fitting logic duplicated here.

Fit window: the concept note's own reference_period (1991-2020, configs/scenarios.yaml),
not the station data's full span -- consistent with every other anomaly/index in this
project being standardised against that one period.

TAVG is used directly where GHCN reports it; where a station/day has TMAX and TMIN but
no TAVG, TAVG = (TMAX+TMIN)/2 is used as a standard meteorological convention (not an
invented value -- this is literally how daily mean temperature is defined in the
absence of sub-daily data), applied uniformly and reported separately in the output
diagnostics so it's never hidden which days came from which source.
"""
import gzip
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.climate.downscale import fit_monthly_lapse_rate, fit_precip_elevation_gradient

import numpy as np
import yaml

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_lapse_rate.yaml"
REFERENCE_PERIOD = (1991, 2020)  # configs/scenarios.yaml's reference_period, not the station data's full span


def load_stations():
    stations = json.loads((DATA_DIR / "_tmp_ghcnd_armenia_stations.json").read_text())
    return {s["station_id"]: s for s in stations}


def load_observations():
    """Returns dict: station_id -> date -> {element: value}."""
    obs = defaultdict(lambda: defaultdict(dict))
    with gzip.open(DATA_DIR / "_tmp_ghcnd_armenia.csv.gz", "rt") as f:
        reader = csv.DictReader(f)
        for row in reader:
            year = int(row["date"][:4])
            if not (REFERENCE_PERIOD[0] <= year <= REFERENCE_PERIOD[1]):
                continue
            obs[row["station_id"]][row["date"]][row["element"]] = float(row["value"])
    return obs


def build_temperature_arrays(stations, obs):
    elevations, temps_c, months, tavg_direct_n, tavg_derived_n = [], [], [], 0, 0
    for station_id, by_date in obs.items():
        if station_id not in stations:
            continue
        elev = stations[station_id]["elevation_m"]
        for date, elements in by_date.items():
            month = int(date[5:7])
            if "TAVG" in elements:
                elevations.append(elev); temps_c.append(elements["TAVG"]); months.append(month)
                tavg_direct_n += 1
            elif "TMAX" in elements and "TMIN" in elements:
                elevations.append(elev); temps_c.append((elements["TMAX"] + elements["TMIN"]) / 2.0); months.append(month)
                tavg_derived_n += 1
    return np.array(elevations), np.array(temps_c), np.array(months), tavg_direct_n, tavg_derived_n


def build_precip_arrays(stations, obs):
    elevations, precip_mm, months = [], [], []
    for station_id, by_date in obs.items():
        if station_id not in stations:
            continue
        elev = stations[station_id]["elevation_m"]
        for date, elements in by_date.items():
            if "PRCP" in elements:
                elevations.append(elev); precip_mm.append(elements["PRCP"]); months.append(int(date[5:7]))
    return np.array(elevations), np.array(precip_mm), np.array(months)


def r_squared_per_month(elevation_m, y, month, gamma_or_grad, log_y=False):
    """Diagnostic only -- not used by the fit itself."""
    out = {}
    for mo in range(1, 13):
        sel = month == mo
        if sel.sum() < 3 or np.isnan(gamma_or_grad[mo - 1]):
            out[mo] = None
            continue
        z, yv = elevation_m[sel], y[sel]
        target = np.log(np.maximum(yv, 1e-6)) if log_y else yv
        pred_slope = gamma_or_grad[mo - 1]
        intercept = np.mean(target) - pred_slope * np.mean(z)
        pred = pred_slope * z + intercept
        ss_res = np.sum((target - pred) ** 2)
        ss_tot = np.sum((target - np.mean(target)) ** 2)
        out[mo] = round(float(1 - ss_res / ss_tot), 3) if ss_tot > 0 else None
    return out


def main():
    stations = load_stations()
    obs = load_observations()
    print(f"{len(stations)} stations, {len(obs)} with observations in {REFERENCE_PERIOD}", flush=True)

    elev_t, t_c, month_t, tavg_direct_n, tavg_derived_n = build_temperature_arrays(stations, obs)
    print(f"temperature: {len(elev_t)} station-days ({tavg_direct_n} direct TAVG, {tavg_derived_n} derived (TMAX+TMIN)/2)", flush=True)
    gamma = fit_monthly_lapse_rate(elev_t, t_c, month_t)
    gamma_r2 = r_squared_per_month(elev_t, t_c, month_t, gamma)

    elev_p, p_mm, month_p = build_precip_arrays(stations, obs)
    print(f"precipitation: {len(elev_p)} station-days", flush=True)
    grad = fit_precip_elevation_gradient(elev_p, p_mm, month_p)
    grad_r2 = r_squared_per_month(elev_p, p_mm, month_p, grad, log_y=True)

    print("\nMonthly lapse rate Gamma_m (K/m) and precip gradient (m^-1):", flush=True)
    month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    for i, name in enumerate(month_names):
        print(f"  {name}: gamma={gamma[i]:.5f} K/m (R2={gamma_r2[i+1]}), "
              f"grad={grad[i]:.6f} /m (R2={grad_r2[i+1]})", flush=True)

    # Physical sanity check, not just "the fit ran": real environmental lapse rates are
    # negative and typically -0.004 to -0.010 K/m (moist to dry adiabatic bounds).
    valid_gamma = gamma[~np.isnan(gamma)]
    implausible = valid_gamma[(valid_gamma > 0.002) | (valid_gamma < -0.015)]
    if len(implausible) > 0:
        print(f"\nWARNING: {len(implausible)} monthly lapse rates outside the physically "
              f"plausible -0.015 to +0.002 K/m range: {implausible}", flush=True)

    result = {
        "fit_date": "2026-09-30",
        "reference_period": list(REFERENCE_PERIOD),
        "source": "configs/manifests/ghcnd_armenia.yaml (53 real GHCN-Daily Armenian stations)",
        "n_stations_with_obs": len(obs),
        "temperature": {
            "n_station_days": int(len(elev_t)),
            "n_tavg_direct": tavg_direct_n,
            "n_tavg_derived_from_tmax_tmin": tavg_derived_n,
            "gamma_k_per_m_by_month": {month_names[i]: (None if np.isnan(gamma[i]) else round(float(gamma[i]), 6)) for i in range(12)},
            "r_squared_by_month": {month_names[i]: gamma_r2[i+1] for i in range(12)},
        },
        "precipitation": {
            "n_station_days": int(len(elev_p)),
            "gradient_per_m_by_month": {month_names[i]: (None if np.isnan(grad[i]) else round(float(grad[i]), 8)) for i in range(12)},
            "r_squared_log_by_month": {month_names[i]: grad_r2[i+1] for i in range(12)},
        },
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.safe_dump(result, sort_keys=False))
    print(f"\nWrote {OUT_PATH}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
