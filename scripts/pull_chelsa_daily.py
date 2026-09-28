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
import os
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
    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with rasterio.Env(GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10, CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
                with rasterio.open(url) as src:
                    window = from_bounds(*BBOX, transform=src.transform)
                    data = src.read(1, window=window)
                    scale = src.scales[0] or 1.0
                    offset = src.offsets[0] or 0.0
                    return d, data.astype(np.float32) * scale + offset, None
        except RasterioIOError as e:
            if "404" in str(e):
                return d, None, "404"  # genuinely missing day, not a transient error
            last_error = e
            time.sleep(min(2 ** attempt, 30))
    # Retries exhausted on something other than a 404 (e.g. a corrupt remote COG strile, as hit
    # for CHELSA_tas_18_02_2020) -- a real, distinct-from-absent failure. Recorded separately so
    # it isn't silently conflated with a confirmed-missing day, and retried on the next run rather
    # than raising and killing the whole multi-day pull over one bad file.
    print(f"{var} {d.isoformat()}: giving up after {MAX_RETRIES} attempts ({last_error})", flush=True)
    return d, None, "read_error"


CHECKPOINT_EVERY = 2000  # a kill mid-run (disk full, network drop, OOM) only loses work since the
                         # last checkpoint, not the whole variable -- lost a full 16.5k-day run to
                         # exactly this before checkpointing existed; see IMPLEMENTATION_LOG.md.


def _checkpoint(var, all_dates, done, missing_404, failed_read, out_path, progress_path):
    # Write to a temp file and rename into place -- np.savez_compressed writing straight to
    # out_path left a truncated, unreadable .npz (BadZipFile) when a checkpoint was interrupted
    # mid-write, destroying the previously-good, fully-fetched file along with it. os.replace
    # is atomic on the same filesystem, so an interruption now only loses the temp file, never
    # the last good checkpoint.
    ordered_dates = [d.isoformat() for d in all_dates if d.isoformat() in done]
    stacked = np.stack([done[ds] for ds in ordered_dates]).astype(np.float32)
    # Must itself end in ".npz" -- np.savez_compressed silently appends ".npz" to any filename
    # that doesn't already end in it, so a ".npz.tmp" path actually gets written to
    # ".npz.tmp.npz", and the rename below then can't find what it just wrote.
    tmp_path = out_path.with_name(out_path.stem + ".tmp.npz")
    np.savez_compressed(tmp_path, data=stacked)
    os.replace(tmp_path, out_path)
    progress_path.write_text(json.dumps(
        {"dates": ordered_dates, "missing_404": missing_404, "failed_read": failed_read}
    ))
    return ordered_dates


def pull_variable(var, start, end):
    out_path = OUT_DIR / f"CHELSA_{var}_2000_{end.year}_armenia_daily.npz"
    progress_path = OUT_DIR / f".{var}_progress.json"

    done = {}
    missing_404 = []
    if progress_path.exists():
        prior = json.loads(progress_path.read_text())
        missing_404 = prior.get("missing_404", prior.get("missing", []))
        if out_path.exists():
            # Decompress the stacked array once -- NpzFile.__getitem__ re-decompresses the whole
            # array on every call, so indexing it inside the loop (the original form here) turned
            # a 14k-day resume into ~14k redundant full-array decompressions and never progressed.
            with np.load(out_path, allow_pickle=True) as saved:
                stacked = saved["data"]
                for i, ds in enumerate(prior["dates"]):
                    done[ds] = stacked[i]
    # failed_read days (retries exhausted on something other than a confirmed 404) are always
    # retried on resume rather than persisted as permanent -- the failure mode (e.g. a corrupt
    # COG strile read) is not known to be permanent the way a 404 is.
    failed_read = []

    dates = list(date_range(start, end))
    todo = [d for d in dates if d.isoformat() not in done and d.isoformat() not in missing_404]
    print(f"=== {var}: {len(dates)} days total, {len(done)} already done, {len(todo)} to fetch ===", flush=True)

    if not todo:
        # Nothing new to fetch -- skip re-stacking and re-compressing the whole (possibly
        # multi-GB) array just to write back the file unchanged.
        ordered_dates = [d.isoformat() for d in dates if d.isoformat() in done]
        print(f"{var}: DONE (already complete). {len(ordered_dates)} days, "
              f"{len(missing_404)} confirmed-missing (404), saved to {out_path}", flush=True)
        return out_path, ordered_dates, missing_404, failed_read

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(fetch_one, var, d): d for d in todo}
        n_done = 0
        for fut in concurrent.futures.as_completed(futures):
            d, arr, reason = fut.result()
            n_done += 1
            if arr is not None:
                done[d.isoformat()] = arr
            elif reason == "404":
                missing_404.append(d.isoformat())
            else:
                failed_read.append(d.isoformat())
            if n_done % 500 == 0:
                print(f"{var}: {n_done}/{len(todo)} fetched this run ({len(done)} total, {len(missing_404)} missing, {len(failed_read)} failed-read)", flush=True)
            if n_done % CHECKPOINT_EVERY == 0:
                _checkpoint(var, dates, done, missing_404, failed_read, out_path, progress_path)
                print(f"{var}: checkpointed at {n_done}/{len(todo)}", flush=True)

    ordered_dates = _checkpoint(var, dates, done, missing_404, failed_read, out_path, progress_path)
    print(f"{var}: DONE. {len(ordered_dates)} days, {len(missing_404)} confirmed-missing (404), "
          f"{len(failed_read)} failed-read (will retry next run), saved to {out_path}", flush=True)
    return out_path, ordered_dates, missing_404, failed_read


if __name__ == "__main__":
    manifest_rows = []
    for var, (start, end) in VARIABLES.items():
        out_path, dates, missing_404, failed_read = pull_variable(var, start, end)
        checksum = hashlib.sha256(out_path.read_bytes()).hexdigest()
        manifest_rows.append({
            "variable": var, "n_days": len(dates),
            "n_missing_404": len(missing_404), "first_missing_404": missing_404[:5],
            "n_failed_read": len(failed_read), "failed_read_dates": failed_read,
            "local_path": str(out_path.relative_to(OUT_DIR.parent.parent)),
            "sha256": checksum,
        })
    Path(OUT_DIR / "chelsa_daily_manifest_rows.json").write_text(json.dumps(manifest_rows, indent=2))
    print("ALL VARIABLES DONE", flush=True)
