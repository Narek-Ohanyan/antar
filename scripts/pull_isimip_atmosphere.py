"""Seasonal and scenario-dependent atmosphere for the water balance: shortwave, longwave, wind and vapour pressure.

Why this exists: the model's wind, radiation and humidity were ERA5-Land annual means, one number per cell repeated for
all 365 days and unchanged in every scenario. `sensitivity_atmospheric_forcing.py` showed that this removes most of the
summer water deficit and flattens the scenario response. This pulls the real seasonal cycle and the real projected
change from ISIMIP, through the same Files API and Armenia bbox as `pull_isimip3b.py`:

* ISIMIP3a obsclim (GSWP3-W5E5, observation-based, daily, 2011-2019 file): the BASELINE seasonal shape, 2015-2019.
* ISIMIP3b bias-adjusted (5 GCMs x SSP1-2.6/3-7.0/5-8.5, 2015-2100): the projected monthly change.

Variables: rsds (shortwave down, W m-2), rlds (longwave down, W m-2), sfcwind (m s-1) and an actual vapour pressure `ea`
(kPa) built per day as hurs/100 * es(tas) (so the temperature-driven rise in humidity is kept) and then averaged.
Surface pressure is left out: its seasonal and scenario variation is ~1% and only enters the psychrometric constant.

Storage-minimal: each job's zip is read, aggregated to MONTHLY means (year, month) and deleted; the .npz holds
(n_months, lat, lon) plus a 'dates' array of month-start dates, so `monthly_means_for_window` works unchanged.
"""
import concurrent.futures
import json
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import pull_isimip3b as P  # noqa: E402
from antar.climate.vapour import saturation_vapour_pressure  # noqa: E402

OUT_DIR = P.OUT_DIR
STAGING = OUT_DIR / "_staging_atmos"
GCMS, SCENARIOS = P.GCMS, P.SCENARIOS
RAW_VARS_3B = ["rsds", "rlds", "sfcwind", "hurs"]        # + tas (already local) for ea
BASELINE_YEARS = (2015, 2019)


def dataset_paths(sim_round, forcing, scenario, var):
    params = {"simulation_round": sim_round, "product": "InputData", "climate_forcing": forcing,
              "climate_scenario": scenario, "climate_variable": var}
    if sim_round == "ISIMIP3a":
        params["time_step"] = "daily"
    resp = requests.get(f"{P.DATA_API_URL}/datasets/", params=params, timeout=30)
    resp.raise_for_status()
    res = resp.json().get("results", [])
    if not res:
        return None
    paths = [f["path"] for f in res[0]["files"]]
    if sim_round == "ISIMIP3a":                           # only the file that holds 2015-2019
        paths = [p for p in paths if p.endswith("_2011_2019.nc")]
    else:                                                 # only the decade files that overlap a window we use (6 of 9 files)
        from antar.climate.atmosphere import BASELINE_WINDOW, HORIZON_WINDOWS
        windows = [BASELINE_WINDOW] + list(HORIZON_WINDOWS.values())
        keep = []
        for p in paths:
            m = re.search(r"_(\d{4})_(\d{4})\.nc$", p)
            if m is None or any(int(m.group(1)) <= w1 and int(m.group(2)) >= w0 for w0, w1 in windows):
                keep.append(p)
        paths = keep
    return paths


def fetch_daily(paths, var, key):
    """Server-side Armenia cutout, then read to (dates, daily array). Staging is deleted before returning."""
    import netCDF4 as nc
    job = P.submit_job(paths)
    finished = P.poll_job(job["job_url"])
    work = STAGING / key
    work.mkdir(parents=True, exist_ok=True)
    zpath = work / "result.zip"
    with requests.get(finished["file_url"], stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(zpath, "wb") as f:
            shutil.copyfileobj(r.raw, f)
    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(work)
    dates, vals = [], []
    for ncf in sorted(work.glob("*.nc")):
        ds = nc.Dataset(ncf)
        tv = ds.variables["time"]
        dts = nc.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
        arr = np.array(ds.variables[var][:], dtype=np.float32)
        fill = getattr(ds.variables[var], "_FillValue", None)
        if fill is not None:
            arr[arr == fill] = np.nan
        dates.extend(d.isoformat()[:10] for d in dts)
        vals.append(arr)
        ds.close()
    shutil.rmtree(work, ignore_errors=True)
    return np.array(dates), np.concatenate(vals, axis=0)


def monthly(dates, values):
    """Daily (time, lat, lon) -> monthly means keyed by 'YYYY-MM-01'."""
    ym = np.array([d[:7] for d in dates])
    keys = sorted(set(ym))
    out = np.stack([np.nanmean(values[ym == k], axis=0) for k in keys])
    return np.array([k + "-01" for k in keys]), out.astype(np.float32)


def local_tas(forcing, scenario):
    """Daily tas (degC) already pulled for the temperature deltas (ISIMIP3b only)."""
    z = np.load(OUT_DIR / f"CHELSA_ISIMIP3b_tas__{forcing}__{scenario}.npz")
    return z["dates"], z["data"] - 273.15


def process(task):
    kind, var, gcm, scenario = task            # kind: '3b' | '3a'
    key = f"{kind}_{var}__{gcm}__{scenario}"
    last = None
    for attempt in range(1, P.MAX_RETRIES + 1):
        try:
            if kind == "3b":
                paths = dataset_paths("ISIMIP3b", gcm, scenario, "hurs" if var == "ea" else var)
            else:
                paths = dataset_paths("ISIMIP3a", "gswp3-w5e5", "obsclim", "hurs" if var == "ea" else var)
            if not paths:
                return key, None, "no_dataset_found"
            dates, vals = fetch_daily(paths, "hurs" if var == "ea" else var, key)
            if var == "ea":                    # actual vapour pressure = RH x saturation at the same day's mean temperature
                if kind == "3b":
                    tdates, tas = local_tas(gcm, scenario)
                else:
                    tpaths = dataset_paths("ISIMIP3a", "gswp3-w5e5", "obsclim", "tas")
                    tdates, tas = fetch_daily(tpaths, "tas", key + "_tas")
                    tas = tas - 273.15
                index = {d: j for j, d in enumerate(np.asarray(tdates).tolist())}
                missing = [d for d in np.asarray(dates).tolist() if d not in index]
                if missing:
                    raise ValueError(f"{key}: {len(missing)} hurs days have no tas day (first {missing[0]}): calendars differ")
                tas = tas[[index[d] for d in np.asarray(dates).tolist()]]
                vals = np.clip(vals, 0, 100) / 100.0 * saturation_vapour_pressure(tas)       # kPa
            mdates, mvals = monthly(dates, vals)
            return key, (mdates, mvals), None
        except Exception as e:                  # noqa: BLE001 -- retried, then reported
            last = str(e)
            import time
            time.sleep(min(4 * attempt, 30))
    return key, None, last


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STAGING.mkdir(parents=True, exist_ok=True)
    VARS = ["rsds", "rlds", "sfcwind", "ea"]
    tasks = [("3a", v, "gswp3-w5e5", "obsclim") for v in VARS] + [("3b", v, g, s) for v in VARS for g in GCMS for s in SCENARIOS]
    prog_path = OUT_DIR / "atmos_pull_progress.json"
    done = json.loads(prog_path.read_text()) if prog_path.exists() else {}
    todo = [t for t in tasks if f"{t[0]}_{t[1]}__{t[2]}__{t[3]}" not in done or done[f"{t[0]}_{t[1]}__{t[2]}__{t[3]}"].get("status") != "done"]
    print(f"=== {len(tasks)} jobs, {len(tasks) - len(todo)} done, {len(todo)} to fetch ===", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(process, t): t for t in todo}
        for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            key, res, err = fut.result()
            if res is None:
                print(f"[{i}/{len(todo)}] {key}: FAILED ({err})", flush=True)
                done[key] = {"status": "failed", "error": err}
            else:
                mdates, mvals = res
                out = OUT_DIR / f"ATMOS_{key}.npz"
                tmp = out.with_name(out.stem + ".tmp.npz")
                np.savez_compressed(tmp, data=mvals, dates=mdates)
                os.replace(tmp, out)
                print(f"[{i}/{len(todo)}] {key}: DONE {mvals.shape} range=[{np.nanmin(mvals):.4g}, {np.nanmax(mvals):.4g}]", flush=True)
                done[key] = {"status": "done", "file": out.name, "n_months": len(mdates)}
            prog_path.write_text(json.dumps(done, indent=1))
    shutil.rmtree(STAGING, ignore_errors=True)
    failed = [k for k, v in done.items() if v["status"] != "done"]
    print("ALL DONE" if not failed else f"{len(failed)} FAILED: {failed}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
