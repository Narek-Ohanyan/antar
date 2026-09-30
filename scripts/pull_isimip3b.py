"""Pull Armenia's window from ISIMIP3b future climate projections (2015-2100): the
concept note's original 5-GCM ensemble under SSP1-2.6/SSP3-7.0/SSP5-8.5, bias-adjusted
against W5E5 v1.0. Primary CMIP6 future source per the 2026-09-30 decision to run this
alongside (not instead of) the CORDEX ensemble on hyperion -- see configs/scenarios.yaml.

Uses the ISIMIP Files API v2 directly (POST to files.isimip.org/api/v2), NOT the
isimip-client package's own cutout() method -- that method sends a request shape
{'task': 'cutout_bbox', ...} the live API rejects; the real schema is
{'paths': [...], 'operations': [{'operation': 'cutout_bbox', 'bbox': [...]}]},
confirmed by reading the API's own root endpoint. bbox order is
[lon_min, lon_max, lat_min, lat_max] -- confirmed by checking a real cutout's actual
output coordinates, not assumed (see IMPLEMENTATION_LOG.md 2026-09-30: the first,
wrong-order attempt "succeeded" and silently produced data for northern Georgia).

One job per (variable, GCM, scenario) bundles all 9 of that combination's ~decade-long
NetCDF files -- confirmed the API accepts multiple paths per job, cutting 540
individual cutouts down to 60. Storage-minimal by design: each job's zip is downloaded,
its 9 NetCDFs read into one concatenated daily array, then deleted immediately --
nothing but the final compact .npz per (variable, GCM, scenario) stays on disk.
"""
import concurrent.futures
import json
import os
import shutil
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import requests

FILES_API_URL = "https://files.isimip.org/api/v2"
DATA_API_URL = "https://data.isimip.org/api/v1"
BBOX = (43.4, 46.7, 38.8, 41.4)  # lon_min, lon_max, lat_min, lat_max -- confirmed order, do not reorder

VARIABLES = ["pr", "tas", "tasmax", "tasmin"]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENARIOS = ["ssp126", "ssp370", "ssp585"]

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "isimip3b"
STAGING_DIR = OUT_DIR / "_staging"
MAX_RETRIES = 4
POLL_SECONDS = 8
MAX_WAIT_SECONDS = 900  # a 9-file job observed to take ~3-4 min; allow generous headroom


def combo_key(var, gcm, scenario):
    return f"{var}__{gcm}__{scenario}"


def get_dataset_paths(var, gcm, scenario):
    resp = requests.get(f"{DATA_API_URL}/datasets/", params={
        "simulation_round": "ISIMIP3b", "product": "InputData",
        "climate_forcing": gcm, "climate_scenario": scenario, "climate_variable": var,
    }, timeout=30)
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        return None
    return [f["path"] for f in results[0]["files"]]


def submit_job(paths):
    payload = {"paths": paths, "operations": [{"operation": "cutout_bbox", "bbox": list(BBOX)}]}
    resp = requests.post(FILES_API_URL, json=payload, timeout=60)
    resp.raise_for_status()
    return resp.json()


def poll_job(job_url):
    waited = 0
    while waited < MAX_WAIT_SECONDS:
        resp = requests.get(job_url, timeout=30)
        resp.raise_for_status()
        job = resp.json()
        if job["status"] == "finished":
            return job
        if job["status"] == "failed":
            raise RuntimeError(f"job failed: {job}")
        time.sleep(POLL_SECONDS)
        waited += POLL_SECONDS
    raise TimeoutError(f"job did not finish within {MAX_WAIT_SECONDS}s: {job_url}")


def process_combo(var, gcm, scenario):
    """Returns (key, dates, values) on success, or (key, None, error_str) on failure."""
    key = combo_key(var, gcm, scenario)
    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            paths = get_dataset_paths(var, gcm, scenario)
            if not paths:
                return key, None, "no_dataset_found"
            job = submit_job(paths)
            finished = poll_job(job["job_url"])

            work_dir = STAGING_DIR / key
            work_dir.mkdir(parents=True, exist_ok=True)
            zip_path = work_dir / "result.zip"
            with requests.get(finished["file_url"], stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(zip_path, "wb") as f:
                    shutil.copyfileobj(r.raw, f)
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(work_dir)

            import netCDF4 as nc
            nc_files = sorted(work_dir.glob("*.nc"))
            all_dates, all_values = [], []
            for ncf in nc_files:
                ds = nc.Dataset(ncf)
                time_var = ds.variables["time"]
                dates = nc.num2date(time_var[:], time_var.units, only_use_cftime_datetimes=False)
                arr = np.array(ds.variables[var][:], dtype=np.float32)  # (time, lat, lon)
                fill = getattr(ds.variables[var], "_FillValue", None)
                if fill is not None:
                    arr[arr == fill] = np.nan
                all_dates.extend(d.isoformat()[:10] for d in dates)
                all_values.append(arr)
                ds.close()

            values = np.concatenate(all_values, axis=0)
            shutil.rmtree(work_dir)  # storage-minimal: raw NetCDFs discarded immediately after reading
            return key, (all_dates, values), None
        except Exception as e:
            last_error = str(e)
            time.sleep(min(4 * attempt, 30))
    return key, None, last_error


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STAGING_DIR.mkdir(parents=True, exist_ok=True)

    combos = [(v, g, s) for v in VARIABLES for g in GCMS for s in SCENARIOS]
    progress_path = OUT_DIR / "pull_progress.json"
    done = {}
    if progress_path.exists():
        done = json.loads(progress_path.read_text())

    todo = [c for c in combos if combo_key(*c) not in done]
    print(f"=== {len(combos)} combos total, {len(done)} already done, {len(todo)} to fetch ===", flush=True)

    # Modest concurrency: each job is a genuine multi-GB server-side operation, not a cheap
    # per-file HTTP fetch like the CHELSA pulls -- 16-way parallelism here would be disrespectful.
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        futures = {ex.submit(process_combo, *c): c for c in todo}
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            var, gcm, scenario = futures[fut]
            key = combo_key(var, gcm, scenario)
            _, result, error = fut.result()
            if result is None:
                print(f"[{i}/{len(todo)}] {key}: FAILED ({error})", flush=True)
                done[key] = {"status": "failed", "error": error}
            else:
                dates, values = result
                out_path = OUT_DIR / f"CHELSA_ISIMIP3b_{key}.npz"
                tmp_path = out_path.with_name(out_path.stem + ".tmp.npz")
                np.savez_compressed(tmp_path, data=values, dates=np.array(dates))
                os.replace(tmp_path, out_path)
                print(f"[{i}/{len(todo)}] {key}: DONE, shape={values.shape}, "
                      f"range=[{np.nanmin(values):.2f}, {np.nanmax(values):.2f}]", flush=True)
                done[key] = {"status": "done", "local_path": str(out_path.relative_to(OUT_DIR.parent.parent)),
                             "n_days": len(dates)}
            progress_path.write_text(json.dumps(done, indent=2))

    shutil.rmtree(STAGING_DIR, ignore_errors=True)
    print("ALL COMBOS DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
