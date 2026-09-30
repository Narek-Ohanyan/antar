"""Pull Armenia's window from CHELSA-BIOCLIM+ (Karger et al. 2022, ESSD) -- historical
(1981-2010) and future (2011-2040/2041-2070/2071-2100, 5 GCMs x 3 SSPs) bioclimatic
variables, via the same HTTP-range-request pattern already proven for CHELSA-daily
(scripts/pull_chelsa_daily.py): windowed /vsicurl/ reads, no full-file download,
incremental checkpointing with atomic writes.

Of the 74 variables the archive actually has (confirmed by listing the bucket, not
assumed from the paper), only the ~37 that map to something ANTAR's engines actually
use are pulled -- a deliberate scope decision, documented here and in
IMPLEMENTATION_LOG.md, not an attempt at completeness for its own sake:

  bio01-bio19        standard bioclim predictors -- MERISTEM's adult-niche Boyce index
  gdd0/5/10          growing-degree-days -- MERISTEM's GDD modifier cross-check
  gddlgd0/5/10        GDD to last growing degree day -- growing-season timing
  gsl/gsp/gst        TREELIM growing-season length/precip/temp -- MERISTEM's treeline ceiling
  fcf, fgd, lgd      frost-change frequency, first/last growing degree day -- late-frost modifier
  vpdmean, vpdmax    vapour-pressure deficit -- XYLEM's water-stress terms
  petmean, petmax    potential evapotranspiration -- TOPOHYDRO's water balance
  sfcWindmean        wind speed -- ERA5-Land cross-check
  rsdsmean           shortwave radiation -- slope-radiation energy term cross-check

Dropped (not pulled): the min/range variants of the above (redundant with mean/max for
this project's purposes), the Koppen-Geiger class bins (kg0-5, not used by any
mechanistic term here), npp/swe/swb/cltmax-range/cmimax-range/hursmax-range/ngd*
(genuinely out of scope for the engines as currently specified -- can be added later
against a real, specific need rather than pulled speculatively now).
"""
import concurrent.futures
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.errors import RasterioIOError

BASE = "/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim"
BBOX = (43.4, 38.8, 46.7, 41.4)
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "chelsa_bioclim"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VARIABLES = [
    *[f"bio{n:02d}" for n in range(1, 20)],
    "gdd0", "gdd5", "gdd10",
    "gddlgd0", "gddlgd5", "gddlgd10",
    "gsl", "gsp", "gst",
    "fcf", "fgd", "lgd",
    "vpdmean", "vpdmax",
    "petmean", "petmax",
    "sfcWindmean",
    "rsdsmean",
]

HISTORICAL_PERIOD = "1981-2010"
FUTURE_PERIODS = ["2011-2040", "2041-2070", "2071-2100"]
GCMS = ["GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"]
SSPS = ["ssp126", "ssp370", "ssp585"]

MAX_RETRIES = 5
CHECKPOINT_EVERY = 20


def scenario_key(period, gcm=None, ssp=None):
    return "historical" if gcm is None else f"{period}|{gcm}|{ssp}"


def all_scenarios():
    yield (HISTORICAL_PERIOD, None, None)
    for period in FUTURE_PERIODS:
        for gcm in GCMS:
            for ssp in SSPS:
                yield (period, gcm, ssp)


def url_for(var, period, gcm, ssp):
    if gcm is None:
        return f"{BASE}/{var}/{period}/CHELSA_{var}_{period}_V.2.1.tif"
    return f"{BASE}/{var}/{period}/{gcm}/{ssp}/CHELSA_{gcm.lower()}_{ssp}_{var}_{period}_V.2.1.tif"


def fetch_one(var, period, gcm, ssp):
    url = url_for(var, period, gcm, ssp)
    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with rasterio.Env(GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10, CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
                with rasterio.open(url) as src:
                    window = from_bounds(*BBOX, transform=src.transform)
                    data = src.read(1, window=window)
                    scale = src.scales[0] or 1.0
                    offset = src.offsets[0] or 0.0
                    # These BIOCLIM+ rasters are int32 with nodata=2147483647 (INT32_MAX) --
                    # unlike CHELSA-daily's uint16 files, reading raw and scaling without masking
                    # first turns nodata into a huge bogus value (2147483647 * 0.1 = 214748364.7),
                    # caught here by checking the actual fetched max, not assuming scale/offset
                    # alone was enough. Masked to NaN before scaling.
                    nodata = src.nodatavals[0]
                    arr = data.astype(np.float32)
                    if nodata is not None:
                        arr[data == nodata] = np.nan
                    return scenario_key(period, gcm, ssp), arr * scale + offset, None
        except RasterioIOError as e:
            if "404" in str(e):
                return scenario_key(period, gcm, ssp), None, "404"
            last_error = e
            time.sleep(min(2 ** attempt, 30))
    print(f"{var} {scenario_key(period, gcm, ssp)}: giving up after {MAX_RETRIES} attempts ({last_error})", flush=True)
    return scenario_key(period, gcm, ssp), None, "read_error"


def _checkpoint(var, done, missing_404, failed_read, out_path, progress_path):
    keys = sorted(done.keys())
    # A checkpoint can legitimately fire with done still empty -- e.g. a variable with no
    # future projections at all (confirmed for vpdmean: only 1981-2010 exists in the archive)
    # has 45 fast 404s racing ahead of the one slow historical-file read across 8 threads, so
    # the first checkpoint boundary can land before that single success has arrived. np.stack
    # needs at least one array; write the progress JSON either way, but only touch the .npz
    # once there is something real to stack.
    if keys:
        stacked = np.stack([done[k] for k in keys]).astype(np.float32)
        tmp_path = out_path.with_name(out_path.stem + ".tmp.npz")
        np.savez_compressed(tmp_path, data=stacked)
        os.replace(tmp_path, out_path)
    progress_path.write_text(json.dumps(
        {"keys": keys, "missing_404": missing_404, "failed_read": failed_read}
    ))
    return keys


def pull_variable(var):
    out_path = OUT_DIR / f"CHELSA_BIOCLIM_{var}_armenia.npz"
    progress_path = OUT_DIR / f".{var}_progress.json"

    done = {}
    missing_404 = []
    if progress_path.exists() and out_path.exists():
        prior = json.loads(progress_path.read_text())
        missing_404 = prior.get("missing_404", [])
        with np.load(out_path) as saved:
            stacked = saved["data"]
            for i, k in enumerate(prior["keys"]):
                done[k] = stacked[i]
    failed_read = []

    scenarios = list(all_scenarios())
    todo = [
        (period, gcm, ssp) for (period, gcm, ssp) in scenarios
        if scenario_key(period, gcm, ssp) not in done
        and scenario_key(period, gcm, ssp) not in missing_404
    ]
    print(f"=== {var}: {len(scenarios)} scenarios total, {len(done)} already done, {len(todo)} to fetch ===", flush=True)

    if not todo:
        print(f"{var}: DONE (already complete). {len(done)} scenarios.", flush=True)
        return sorted(done.keys()), missing_404, failed_read

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futures = {ex.submit(fetch_one, var, period, gcm, ssp): (period, gcm, ssp) for (period, gcm, ssp) in todo}
        n_done = 0
        for fut in concurrent.futures.as_completed(futures):
            key, arr, reason = fut.result()
            n_done += 1
            if arr is not None:
                done[key] = arr
            elif reason == "404":
                missing_404.append(key)
            else:
                failed_read.append(key)
            if n_done % CHECKPOINT_EVERY == 0:
                _checkpoint(var, done, missing_404, failed_read, out_path, progress_path)
                print(f"{var}: checkpointed at {n_done}/{len(todo)}", flush=True)

    keys = _checkpoint(var, done, missing_404, failed_read, out_path, progress_path)
    print(f"{var}: DONE. {len(keys)} scenarios, {len(missing_404)} confirmed-missing (404), "
          f"{len(failed_read)} failed-read (will retry next run)", flush=True)
    return keys, missing_404, failed_read


def main():
    manifest_rows = []
    for var in VARIABLES:
        out_path = OUT_DIR / f"CHELSA_BIOCLIM_{var}_armenia.npz"
        keys, missing_404, failed_read = pull_variable(var)
        import hashlib
        checksum = hashlib.sha256(out_path.read_bytes()).hexdigest() if out_path.exists() else None
        manifest_rows.append({
            "variable": var, "n_scenarios": len(keys),
            "n_missing_404": len(missing_404), "missing_404": missing_404,
            "n_failed_read": len(failed_read),
            "local_path": str(out_path.relative_to(OUT_DIR.parent.parent)) if out_path.exists() else None,
            "sha256": checksum,
        })
    Path(OUT_DIR / "chelsa_bioclim_manifest_rows.json").write_text(json.dumps(manifest_rows, indent=2))
    print("ALL VARIABLES DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
