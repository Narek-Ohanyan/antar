"""Pull Armenia's window from CHELSA-daily (native CHELSA V2.1, historical reference)
directly via HTTP range requests against the public COGs -- no full-file download.

tas/tasmax/tasmin: 1979-2024 (confirmed present). pr: 1979-2019 only (confirmed gap
2020-2024 in the public archive's own "active (incomplete)" status). Requires the
``geo`` extra (rasterio) -- see pyproject.toml.

Filenames are CHELSA_{var}_{DD}_{MM}_{YYYY} (day then month) -- verified empirically
against the actual seasonal temperature signal, not assumed; see IMPLEMENTATION_LOG.md.
"""
import concurrent.futures
import datetime
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.errors import RasterioIOError

BASE = "/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/daily"
BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "chelsa"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VARIABLES = {
    # tas/tasmax/tasmin: full available range. pr: available range ends 2019 (confirmed
    # gap 2020-2024 in the public archive's own "active (incomplete)" status).
    # 1979 start matches the reanalysis-era convention (ERA5 also starts 1979) and, more
    # importantly, actually covers the concept note's stated reference_period (1991-2020,
    # scenarios.yaml) -- the earlier 2000-2024 window did not.
    "tas": (datetime.date(1979, 1, 1), datetime.date(2024, 12, 31)),
    "tasmax": (datetime.date(1979, 1, 1), datetime.date(2024, 12, 31)),
    "tasmin": (datetime.date(1979, 1, 1), datetime.date(2024, 12, 31)),
    "pr": (datetime.date(1979, 1, 1), datetime.date(2019, 12, 31)),
}

MAX_RETRIES = 5


def date_range(start, end):
    d = start
    while d <= end:
        yield d
        d += datetime.timedelta(days=1)


def fetch_one(var, d):
    # Filenames are CHELSA_{var}_{DD}_{MM}_{YYYY} -- day *then* month, confirmed empirically
    # (the naive month_day guess returned real files for day<=12, silently reading the wrong
    # calendar date every time: e.g. "01_10_2000" is Oct 1, not Jan 10 -- caught by checking
    # the actual mean temperature, ~10 degC, unambiguously October not January in Armenia).
    url = f"{BASE}/{var}/{d.year}/CHELSA_{var}_{d.day:02d}_{d.month:02d}_{d.year}_V.2.1.tif"
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with rasterio.Env(GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10, CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
                with rasterio.open(url) as src:
                    window = from_bounds(*BBOX, transform=src.transform)
                    data = src.read(1, window=window)
                    scale = src.scales[0] or 1.0
                    offset = src.offsets[0] or 0.0
                    return d, data.astype(np.float32) * scale + offset, data.shape
        except RasterioIOError as e:
            if "404" in str(e):
                return d, None, None  # genuinely missing day, not a transient error
            if attempt == MAX_RETRIES:
                raise
            time.sleep(min(2 ** attempt, 30))
    return d, None, None


def pull_variable(var, start, end):
    out_path = OUT_DIR / f"CHELSA_{var}_2000_{end.year}_armenia_daily.npz"
    progress_path = OUT_DIR / f".{var}_progress.json"

    done = {}
    if progress_path.exists():
        saved = np.load(out_path, allow_pickle=True) if out_path.exists() else None
        prior_dates = json.loads(progress_path.read_text())
        if saved is not None:
            for i, ds in enumerate(prior_dates["dates"]):
                done[ds] = saved["data"][i]

    dates = list(date_range(start, end))
    results = {}
    missing = []
    shape_ref = None

    todo = [d for d in dates if d.isoformat() not in done]
    print(f"=== {var}: {len(dates)} days total, {len(done)} already done, {len(todo)} to fetch ===", flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(fetch_one, var, d): d for d in todo}
        n_done = 0
        for fut in concurrent.futures.as_completed(futures):
            d, arr, shape = fut.result()
            n_done += 1
            if arr is None:
                missing.append(d.isoformat())
            else:
                done[d.isoformat()] = arr
                shape_ref = shape
            if n_done % 500 == 0:
                print(f"{var}: {n_done}/{len(todo)} fetched this run ({len(done)} total, {len(missing)} missing)", flush=True)

    ordered_dates = [d.isoformat() for d in dates if d.isoformat() in done]
    stacked = np.stack([done[ds] for ds in ordered_dates]).astype(np.float32)
    np.savez_compressed(out_path, data=stacked)
    progress_path.write_text(json.dumps({"dates": ordered_dates, "missing": missing}))
    print(f"{var}: DONE. shape={stacked.shape}, {len(missing)} missing days, saved to {out_path}", flush=True)
    return out_path, ordered_dates, missing


if __name__ == "__main__":
    manifest_rows = []
    for var, (start, end) in VARIABLES.items():
        out_path, dates, missing = pull_variable(var, start, end)
        checksum = hashlib.sha256(out_path.read_bytes()).hexdigest()
        manifest_rows.append({
            "variable": var, "n_days": len(dates), "n_missing": len(missing),
            "first_missing": missing[:5], "local_path": str(out_path.relative_to(OUT_DIR.parent.parent)),
            "sha256": checksum,
        })
    Path(OUT_DIR / "chelsa_daily_manifest_rows.json").write_text(json.dumps(manifest_rows, indent=2))
    print("ALL VARIABLES DONE", flush=True)
