"""The first real, gridded TOPOHYDRO forcing run: ties every piece built this
session (the real monthly lapse-rate/precip-gradient fit, the Saxton & Rawls
pedotransfer, the terrain concavity index, the ERA5-derived net radiation and
wind-height correction) together through antar.climate.forcing.topoclimate_forcing,
against entirely real SRTM/SoilGrids/ERA5-Land/CHELSA-daily data, for the first
time end to end.

Deliberately reduced scope, stated plainly rather than presented as a final
production run:

* 80 grid points -- an 8x10 regular subsample of CHELSA-daily's own 312x396
  grid (the same subsample verified earlier this session), not the full 30m
  master grid (~92M cells is a separate, much larger computation).
* One year, 2019 -- the most recent year common to all four CHELSA-daily
  variables (pr stops at 2019-12-31; tas/tasmax/tasmin run to 2024), not the
  full 1991-2020 reference period.

Two inputs have no real data source pulled for this run and are explicit,
documented placeholders (not silently defaulted):

* ``calm_clear_night_frac`` -- no cloud-cover product was pulled (ERA5-Land's
  own export here is wind/radiation/dewpoint/pressure only). Fixed at 0.3,
  matching k_cap's own "to be calibrated" status in antar.climate.downscale.
* ``gdd_budburst`` -- species/functional-group specific per the concept note;
  no species-specific value has been derived yet. 200 GDD-days (base 5 degC)
  stands in as a generic literature-typical temperate-tree order of magnitude.
* ``rooting_depth_mm`` (1000 mm here) -- feeds w_max_mm = (theta_fc - theta_lim)
  * rooting_depth. This standalone run stays species-agnostic (one generic
  value for the whole grid) by design; real per-group values now exist
  (``ROOTING_DEPTH_MM_BY_GROUP``, Canadell et al. 1996, 2.9-9.5 m depending on
  functional group) and callers that know their real species (XYLEM, REFUGIUM,
  future-projections, MERISTEM's CWD extraction) pass them through
  ``extract_static_grid_inputs``/``compute_grid_forcing``'s ``rooting_depth_mm``
  parameter instead -- each calling this module's forcing pipeline once per
  real functional group rather than once per cell, since w_max_mm (and
  therefore the whole real water-balance signal: CWD, WSI, soil psi) was
  previously shared identically across all 4 species, which is itself a real
  simplification worth having fixed, not just the constant's value.

Two real, defensible modelling choices specific to this run (not framework
defaults, documented here rather than in library code):

* ``z_ref_m`` (the elevation the lapse-rate downscaling treats CHELSA's ~1km
  pixel as representing) is the SRTM elevation averaged over a ~930 m window
  (31x31 30m pixels) centred on each point, i.e. an explicit coarse-scale
  proxy -- not the exact 30m cell elevation used elsewhere as z_cell_m, and
  not a single fixed reference station. This lets the lapse-rate correction
  do real work (correcting for real within-CHELSA-pixel elevation variation)
  rather than degenerately downscaling a point to itself. The same z_ref_m is
  reused for the dew-point elevation adjustment (ea_ref_kpa): ERA5-Land's own
  native resolution (~9km) is coarser than CHELSA's, so treating it as
  representing the same coarse elevation is, if anything, an under-correction.
* Net radiation (rn_mj_m2, an input this function requires from upstream) is
  built from ERA5-Land's annual-mean ssrd/strd broadcast as a daily constant
  across 2019, using an emitted-longwave term evaluated at each day's
  downscaled t_mean_c. Since t_mean_c is itself computed inside
  topoclimate_forcing, a lapse-downscaled t_mean_c is computed once, upfront,
  with the same antar.climate.downscale.downscale_temperature call
  topoclimate_forcing makes internally, purely to drive this Rn term -- a
  harmless duplicate of a cheap, pure computation, not a divergent estimate.

terrain.tif, soils.tif and the relevant ERA5-Land tiles are all read via GDAL's
/vsicurl/ HTTP range-request streaming (Drive's download endpoint, authenticated
with the existing Earth Engine OAuth refresh token) rather than downloaded whole:
terrain.tif is 421 MB, soils.tif is 416 MB and the 16 ERA5-Land tiles are 15+ GB
combined, and only a few dozen small windows/pixels are actually needed --
consistent with this project's minimize-local-footprint convention. Nothing is
left on disk by this script.
"""
import datetime
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.climate import downscale
from antar.climate.atmosphere import AtmosphereShape
from antar.io.armenia_mask import load_mask as load_armenia_mask
from antar.climate.forcing import topoclimate_forcing
from antar.climate.radiation import net_radiation_from_era5
from antar.climate.soil_pedotransfer import soil_hydraulic_parameters
from antar.climate.terrain import concavity_index
from antar.climate.vapour import saturation_vapour_pressure, wind_speed_2m

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_grid_run_2019.yaml"

BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max
CHELSA_GRID_SHAPE = (312, 396)
GRID_ROWS = np.arange(0, 312, 40)   # 8 rows -- the original 80-point validation grid
GRID_COLS = np.arange(0, 396, 40)   # 10 cols
# The dense grid (redefined 2026-10-05). The first version (stride 11 over the whole rectangle, 1044
# points) was found -- while building the UI -- to be about two-thirds foreign territory: Armenia
# fills only 36.7% of the rectangle (45,322 of 123,552 CHELSA cells), so only ~330 of those 1044
# cells were in Armenia and ~2/3 of the compute was spent on Georgia, Azerbaijan, Turkey and Iran.
# The dense grid is now a stride-7 subsample of CHELSA's native cells RESTRICTED TO CELLS INSIDE
# ARMENIA (src/antar/io/armenia_mask.py): ~925 cells, ~5 x 6.5 km apart, all Armenian. Stride 7
# keeps the dominant cost (future-projections' 45 members x cells x 3 groups x ~1.3 s Monte Carlo)
# near ~45 h while giving ~2.8x more Armenian cells than the old dense grid.
class _ArmeniaOnlyAxis(np.ndarray):
    """A grid axis whose grid is restricted to cells inside Armenia (see grid_latlon)."""
    inside_armenia_only = True


DENSE_STRIDE = 7
DENSE_GRID_ROWS = np.arange(0, 312, DENSE_STRIDE).view(_ArmeniaOnlyAxis)
DENSE_GRID_COLS = np.arange(0, 396, DENSE_STRIDE).view(_ArmeniaOnlyAxis)
DENSE_GRID_ID = f"dense_armenia_stride{DENSE_STRIDE}"
DENSE_GRID_LABEL = f"Armenia-only dense grid (stride {DENSE_STRIDE})"
YEAR = 2019

TERRAIN_DRIVE_ID = "17zOkIKhKZiDQaRtF2SwPhk3XOa37hpbe"
SOILS_DRIVE_ID = "18LRnI4Nsnlj6ks1UClGDe4BdufqUG24n"
CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "_cache"


def _points_key(a, b):
    """Stable key for a set of grid points (rounded so float noise cannot split a cache)."""
    h = hashlib.sha1(np.round(np.concatenate([np.asarray(a, float), np.asarray(b, float)]), 6).tobytes()).hexdigest()[:14]
    return f"{len(a)}_{h}"


def _save_npz(path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


ERA5LAND_TILE_SIZE_PX = 3072  # GEE shard size confirmed against terrain.tif's own pixel grid

CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER = 0.3
GDD_BUDBURST_PLACEHOLDER = 200.0
ROOTING_DEPTH_MM_PLACEHOLDER = 1000.0  # generic default -- this script's own standalone grid run
# is deliberately species-agnostic (see module docstring); callers that know their real
# functional group should pass ROOTING_DEPTH_MM_BY_GROUP[group] to extract_static_grid_inputs/
# compute_grid_forcing instead. Real values: Canadell et al. (1996), Oecologia 108, 583-595 --
# biome-level (not genus-level: the paper reports by biome/functional-type, not individually for
# Fagus/Quercus/Pinus/Juniperus), still a real, large improvement over one flat guess shared by
# every species. mesic_diffuse_porous_broadleaf and ring_porous_oak share the same real value --
# the paper's temperate-deciduous-forest figure covers both (Fagus orientalis/Carpinus betulus AND
# Quercus macranthera/iberica fall in the same biome class, not two different real numbers).
ROOTING_DEPTH_MM_BY_GROUP = {
    "mesic_diffuse_porous_broadleaf": 2900.0,  # temperate deciduous forest
    "ring_porous_oak": 2900.0,                 # temperate deciduous forest (same biome class)
    "pine": 3900.0,                            # temperate coniferous forest
    "juniper_arid_conifer": 9500.0,            # desert/arid shrubland (arid-adapted, deep-rooted)
}

COARSE_REF_WINDOW_PX = 15  # +-15 px at 30m = 930m, matching CHELSA-daily's own ~927m pixel
CONCAVITY_WINDOW_PX = 3    # 7x7 (210m) local neighbourhood for cold-air-pooling concavity

ERA5LAND_TILE_IDS = {
    (0, 0): "1Osy7hh3XbhDJgOTuXnIKGvMeA1RLIG_H", (0, 3072): "1CJB9LZtDaIfVEXzbGV8N-skFElU17v7e",
    (0, 6144): "1XDzUgFxdulZWTvN90-x5fmhp018wKclz", (0, 9216): "1Ordcx93EtwZTRvaTpnp9w2aZwlWt5u68",
    (3072, 0): "1DFpW0KFisx2C9778QH66MHnibEuVLgY_", (3072, 3072): "1QrtWlHR19aGLl5p7D5BJ3EKTx4SzJK7J",
    (3072, 6144): "1BLG95JTJZ1soePO1bureX5UtXgmHizlY", (3072, 9216): "1_9Ajxw3rySxh9fOzZJMKfxWy3nJQ3fiv",
    (6144, 0): "1UT6WqhnlNftbhagsZ__ZlLOQMF1nHzd5", (6144, 3072): "113jX1fWiW0F4dawd9tC7KWSpByQ2PGWn",
    (6144, 6144): "1aYXhxaWjYQhl0qVZQHvMcOJAy_-TU6bQ", (6144, 9216): "1KiQzUMLM01Q5E-y5n-Qgi3krOXEOAR5W",
    (9216, 0): "1zIiisaynM-CWe_QreVrCT9qIRQLgOqGt", (9216, 3072): "13hrIrruW5H-Ux_CzbWy8RNYkLQH6gFjJ",
    (9216, 6144): "16V2QKatPrNfiNe1m5whc42L8vKHjvGN2", (9216, 9216): "1DvIQ6aV-rVge00BH4eXZPNrT0u7i-c9v",
}
# Bilinear-resampled, re-exported 2026-10-01 (configs/resampling_policy.yaml's fix for the real
# ~300x ERA5-Land-to-master-grid upsample). Every real run through REFUGIUM/future-projections/
# AEGIS/MERISTEM-in-progress used the ORIGINAL nearest-neighbor tiles, kept here for any re-run
# that needs to reproduce those exact numbers -- see configs/manifests/gee_exports.yaml for the
# full real id/size/MD5 record of both versions.
ERA5LAND_TILE_IDS_NEAREST_NEIGHBOR_ORIGINAL = {
    (0, 0): "1kQxh_EDcCn6Ih5xyfvfIZuQE2mWNPqaS", (0, 3072): "1ZCNxrs6P3-_Y8DYg-MetdslHjQnV0phA",
    (0, 6144): "13YdynAAjD0KwDYzpMt3C5GmnH3eKEDHP", (0, 9216): "1rm4-5kX8tO5i855b3BGEP7cXqbb8htzT",
    (3072, 0): "1psHQxLuxdRcN1Oczv7VDh1mcB8a7N38s", (3072, 3072): "1LkBgsLjeZ1_uRlO1c6FHmsZ8pxpIPjnT",
    (3072, 6144): "1FE5mcNk6xONJ4sChXZ7qfmgWyS_45Uk1", (3072, 9216): "1C5RBC59v1ytFwsTWcCme0VM-VrG8wORe",
    (6144, 0): "1wtVOpllg9Q-PmQDaJeetU8ZIbKI1pjMr", (6144, 3072): "1zP8XI026RlIaFCR-3LnF1qxdbQSORnDf",
    (6144, 6144): "1_8tMT-0SnBrRUepjOcWzIo4hjcQrAtcH", (6144, 9216): "1dp-96zRTfqwPDTTdn95tv1-yjmjHJDUp",
    (9216, 0): "1t-hg8gk3sKbWCSvFi7uom7PjGI5SNCwX", (9216, 3072): "1q6KXWJQAIOQgP-OGLtyTS7L3G7YIOvpk",
    (9216, 6144): "1qxYjAt8IwXexaStDAJCphgYXSCDxUzH_", (9216, 9216): "1sJHL1cLIxn9mSwpL9FrPXQBxl0LO1yFT",
}


def grid_latlon(grid_rows=GRID_ROWS, grid_cols=GRID_COLS):
    """Cell-centre lat/lon and CHELSA (row, col) for every grid cell. If either axis is an
    ``_ArmeniaOnlyAxis`` (the dense grid), cells outside Armenia are dropped."""
    lats = BBOX[3] - (np.asarray(grid_rows) + 0.5) * (BBOX[3] - BBOX[1]) / CHELSA_GRID_SHAPE[0]
    lons = BBOX[0] + (np.asarray(grid_cols) + 0.5) * (BBOX[2] - BBOX[0]) / CHELSA_GRID_SHAPE[1]
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    chelsa_row, chelsa_col = np.meshgrid(np.asarray(grid_rows), np.asarray(grid_cols), indexing="ij")
    lat_f, lon_f, row_f, col_f = lat_grid.ravel(), lon_grid.ravel(), chelsa_row.ravel(), chelsa_col.ravel()
    if getattr(grid_rows, "inside_armenia_only", False) or getattr(grid_cols, "inside_armenia_only", False):
        keep = load_armenia_mask()[row_f, col_f]
        lat_f, lon_f, row_f, col_f = lat_f[keep], lon_f[keep], row_f[keep], col_f[keep]
    return lat_f, lon_f, row_f, col_f


def get_access_token():
    """Real OAuth token refresh, retried on transient transport failure.

    A real crash hit mid-REFUGIUM-dense-run (2026-10-02): the token endpoint itself
    (oauth2.googleapis.com) dropped the connection (google.auth.exceptions.TransportError /
    requests.exceptions.ConnectionError), a failure mode entirely outside rasterio's scope -- every
    retry wrapper added tonight only caught rasterio.errors.RasterioIOError around the *raster*
    read, not a failure in the token refresh those retries call fresh each attempt. Fixing it once
    here, at the source, protects every caller automatically rather than teaching each retry loop
    about a second exception type.
    """
    import ee
    from google.oauth2.credentials import Credentials
    import google.auth.transport.requests as gareq
    import google.auth.exceptions
    import requests.exceptions

    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    last_err = None
    for attempt in range(3):
        try:
            creds = Credentials(
                None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
                client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"],
            )
            creds.refresh(gareq.Request())
            return creds.token
        except (google.auth.exceptions.TransportError, requests.exceptions.ConnectionError) as e:
            last_err = e
            print(f"  get_access_token: refresh failed (attempt {attempt + 1}/3): {e}", flush=True)
            time.sleep(5)
    raise last_err


class DriveCoverageError(RuntimeError):
    """Raised when a Drive-streamed extraction loses too many points to read failures. A real
    failure mode (2026-10-02/03): during a multi-hour Drive disruption, per-tile retries "succeeded"
    at degrading gracefully -- leaving the points NaN -- so the dense XYLEM/REFUGIUM runs finished
    "cleanly" with only 212/1044 and 418/1044 cells and wrote output files that looked like valid
    dense results but were badly biased subsamples. Graceful degradation is right for a stray bad
    tile; it is wrong when most of the grid is lost, because the result then silently misrepresents
    its own coverage. Fail loudly instead and let the caller retry once Drive recovers."""


MAX_LOST_FRACTION = 0.03  # >3% of points lost to read failures (beyond legitimately-missing soil) = abort
RETRY_SLEEPS_S = [5, 30, 120]  # exponential-ish backoff: a Drive disruption can last minutes, not seconds


def drive_vsicurl_url(file_id: str) -> str:
    return f"/vsicurl/https://www.googleapis.com/drive/v3/files/{file_id}?alt=media"


def extract_terrain(token, lats, lons):
    """Real per-point terrain windows (concavity + coarse z_ref), retried per-point with a fresh
    connection/token on transient failure. One connection serving all N points was fine at 80
    points but a real risk at the 1044-point dense grid: any single transient vsicurl failure
    partway through would previously crash the whole function and lose all prior points' work --
    a real concern for an unattended multi-hour overnight run. Re-opening per point costs a little
    overhead but bounds the real damage of one bad point to that point alone (left NaN after 3
    failed attempts), not a full restart."""
    xs, ys = None, None
    n = len(lats)
    elevation, slope, aspect, concavity, z_ref = (np.full(n, np.nan) for _ in range(5))
    row_px, col_px = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    failed_points = []

    for i in range(n):
        last_err = None
        for attempt in range(3):
            try:
                fresh_token = get_access_token()
                url = drive_vsicurl_url(TERRAIN_DRIVE_ID)
                with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {fresh_token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                    with rasterio.open(url) as src:
                        x, y = warp_transform("EPSG:4326", src.crs, [lons[i]], [lats[i]])
                        x, y = x[0], y[0]
                        row, col = src.index(x, y)
                        row_px[i], col_px[i] = row, col

                        cw = CONCAVITY_WINDOW_PX
                        r0, r1 = max(0, row - cw), min(src.height, row + cw + 1)
                        c0, c1 = max(0, col - cw), min(src.width, col + cw + 1)
                        elev_win = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0))
                        slope_win = src.read(2, window=Window(c0, r0, c1 - c0, r1 - r0))
                        aspect_win = src.read(3, window=Window(c0, r0, c1 - c0, r1 - r0))
                        lr, lc = row - r0, col - c0
                        elevation[i] = elev_win[lr, lc]
                        slope[i] = slope_win[lr, lc]
                        aspect[i] = aspect_win[lr, lc]
                        concavity[i] = concavity_index(elev_win)[lr, lc]

                        rw = COARSE_REF_WINDOW_PX
                        r0c, r1c = max(0, row - rw), min(src.height, row + rw + 1)
                        c0c, c1c = max(0, col - rw), min(src.width, col + rw + 1)
                        coarse_win = src.read(1, window=Window(c0c, r0c, c1c - c0c, r1c - r0c))
                        z_ref[i] = float(np.nanmean(coarse_win))
                last_err = None
                break
            except rasterio.errors.RasterioIOError as e:
                last_err = e
                time.sleep(2)
        if last_err is not None:
            failed_points.append(i)
        if (i + 1) % 200 == 0:
            print(f"    terrain: {i + 1}/{n}", flush=True)
    if failed_points:
        print(f"  terrain: {len(failed_points)} point(s) failed all retries, left NaN", flush=True)
        if len(failed_points) / n > MAX_LOST_FRACTION:
            raise DriveCoverageError(f"terrain lost {len(failed_points)}/{n} points to read failures")
    return elevation, slope, aspect, concavity, z_ref, row_px, col_px


def extract_soils(token, lats, lons):
    """Real soils sample, retried with a fresh token on failure. A real crash hit this function
    with a bare HTTP 401 mid-REFUGIUM-dense-run (2026-10-02): it used the single `token` captured
    once at the top of extract_static_grid_inputs rather than refreshing its own, unlike
    extract_terrain/extract_era5land which already learned this lesson earlier tonight -- a real
    gap, not a one-off, now closed here too."""
    last_err = None
    for attempt in range(3):
        try:
            fresh_token = get_access_token()
            url = drive_vsicurl_url(SOILS_DRIVE_ID)
            with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {fresh_token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                with rasterio.open(url) as src:
                    band_names = list(src.descriptions)
                    xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
                    vals = np.array(list(src.sample(zip(xs, ys))))
            last_err = None
            break
        except rasterio.errors.RasterioIOError as e:
            last_err = e
            print(f"  soils: read failed (attempt {attempt + 1}/3): {e}", flush=True)
            time.sleep(5)
    if last_err is not None:
        print(f"  soils: FAILED after 3 attempts: {last_err}", flush=True)
        raise DriveCoverageError(f"soils unreadable after retries: {last_err}")
    idx = {b: i for i, b in enumerate(band_names)}
    # SoilGrids' GEE-mapped units are per-mille (g/kg, 0-1000): verified empirically this
    # session -- clay+sand+silt sum to ~1000 at every real sample point -- so /10 gives the
    # percent (clay/sand) and g/kg (soc) units antar.climate.soil_pedotransfer expects.
    clay_pct = vals[:, idx["clay_0_30cm_mean"]] / 10.0
    sand_pct = vals[:, idx["sand_0_30cm_mean"]] / 10.0
    soc_g_kg = vals[:, idx["soc_0_30cm_mean"]] / 10.0
    return clay_pct, sand_pct, soc_g_kg


def extract_era5land(token, row_px, col_px, year=YEAR):
    """Per-tile ERA5-Land read with retry-on-exception (fresh token each attempt). A real
    transient vsicurl read failure crashed this function's pre-retry version mid-XYLEM-run
    (2026-10-01, RasterioIOError inside a plain per-point loop with no retry) -- the GDAL
    HTTP timeout alone bounds how long a single stalled request hangs, but does nothing for a
    request that fails fast and needs a retry; see compute_real_cwd_for_meristem.py's
    extract_era5land_vectorized, which already had this pattern, for the precedent."""
    n = len(row_px)
    cache = CACHE_DIR / f"era5_{year}_{_points_key(row_px, col_px)}.npz"
    if cache.exists():
        z = np.load(cache)
        print(f"=== ERA5-Land {year} for {n} points loaded from {cache.name} (no Drive access) ===", flush=True)
        return z["wind10"], z["ssrd"], z["strd"], z["dewpoint_k"], z["pressure_pa"]
    if token is None:
        token = get_access_token()
    wind10, ssrd, strd, dewpoint_k, pressure_pa = (np.full(n, np.nan) for _ in range(5))
    tile_row = (row_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX
    tile_col = (col_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX

    failed_tiles = []
    for tile_key in sorted(set(zip(tile_row.tolist(), tile_col.tolist()))):
        file_id = ERA5LAND_TILE_IDS[tile_key]
        sel = (tile_row == tile_key[0]) & (tile_col == tile_key[1])

        last_err = None
        for attempt in range(len(RETRY_SLEEPS_S) + 1):
            try:
                fresh_token = get_access_token()
                url = drive_vsicurl_url(file_id)
                with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {fresh_token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                    with rasterio.open(url) as src:
                        bidx = {b: j + 1 for j, b in enumerate(src.descriptions)}
                        if f"wind_speed_{year}" not in bidx:
                            # Drive occasionally serves a non-GeoTIFF body (HTML quota/error page)
                            # that GDAL opens as an empty dataset with no band names -- a real
                            # failure that surfaced as an uncaught KeyError, not a RasterioIOError.
                            raise rasterio.errors.RasterioIOError(
                                f"tile {tile_key} opened without band '{'wind_speed'}_{year}' "
                                f"(got {len(bidx)} named bands) -- not a valid ERA5 tile response")
                        tile_vals = {k: np.full(n, np.nan) for k in
                                     ["wind_speed", "ssrd", "strd", "dewpoint", "surface_pressure"]}
                        for i in np.where(sel)[0]:
                            lr, lc = int(row_px[i] - tile_key[0]), int(col_px[i] - tile_key[1])
                            w = Window(lc, lr, 1, 1)
                            for bandname in tile_vals:
                                tile_vals[bandname][i] = src.read(bidx[f"{bandname}_{year}"], window=w)[0, 0]
                for arr, key in [(wind10, "wind_speed"), (ssrd, "ssrd"), (strd, "strd"),
                                  (dewpoint_k, "dewpoint"), (pressure_pa, "surface_pressure")]:
                    arr[sel] = tile_vals[key][sel]
                last_err = None
                break
            except rasterio.errors.RasterioIOError as e:
                last_err = e
                print(f"  tile {tile_key}: read failed (attempt {attempt + 1}/{len(RETRY_SLEEPS_S) + 1}): {e}", flush=True)
                if attempt < len(RETRY_SLEEPS_S):
                    time.sleep(RETRY_SLEEPS_S[attempt])
        if last_err is not None:
            print(f"  tile {tile_key}: FAILED after {len(RETRY_SLEEPS_S) + 1} attempts, {sel.sum()} points left NaN", flush=True)
            failed_tiles.append(tile_key)
            continue
        print(f"  tile {tile_key}: {sel.sum()} points", flush=True)
    lost = int(np.isnan(wind10).sum())
    if lost / n > MAX_LOST_FRACTION:
        raise DriveCoverageError(
            f"ERA5-Land lost {lost}/{n} points ({len(failed_tiles)} tile(s) failed all retries: "
            f"{failed_tiles}) -- refusing to continue with a badly degraded sample")
    if failed_tiles:
        print(f"  === {len(failed_tiles)} tile(s) failed all retries: {failed_tiles} ===", flush=True)
    else:                                   # only a complete read is cached; a degraded one is re-read next time
        _save_npz(cache, {"wind10": wind10, "ssrd": ssrd, "strd": strd, "dewpoint_k": dewpoint_k, "pressure_pa": pressure_pa})
        print(f"=== cached ERA5-Land {year} -> {cache.name} ===", flush=True)
    return wind10, ssrd, strd, dewpoint_k, pressure_pa


def extract_static_grid_inputs(grid_rows=GRID_ROWS, grid_cols=GRID_COLS):
    """Terrain and soils: year-independent AND species-independent, extracted once and reused
    across years and across real functional groups.

    Factored out so a multi-year caller (MNEME's hazard-panel build) pays for terrain/soil
    extraction once, not once per year. Does NOT bake in rooting depth -- clay/sand/soc ->
    theta_fc/theta_lim/etc. (``soil_hydraulic_parameters``) are genuinely species-independent soil
    properties; only w_max_mm = (theta_fc - theta_lim) * rooting_depth varies by real functional
    group, and that multiply is cheap and local (no network), so it happens in
    :func:`compute_forcing_for_year` instead, where a per-group caller can repeat it 4 times for
    ~free rather than re-streaming terrain/soils/ERA5-Land 4 times for no reason.

    ``grid_rows``/``grid_cols`` default to the original 80-point validation grid; pass
    ``DENSE_GRID_ROWS``/``DENSE_GRID_COLS`` for the Armenia-only dense grid.
    """
    lats, lons, chelsa_row, chelsa_col = grid_latlon(grid_rows, grid_cols)
    n = len(lats)
    print(f"=== {n} grid points ({len(grid_rows)}x{len(grid_cols)}) ===", flush=True)

    # Terrain, soils and the ERA5 pixel indices take ~1-2 h to stream for the dense grid and never change, so they are
    # cached once per grid; every later step (XYLEM, REFUGIUM, scenarios, MNEME) reads the cache instead of Drive.
    cache = CACHE_DIR / f"static_{_points_key(lats, lons)}.npz"
    if cache.exists():
        z = np.load(cache)
        print(f"=== static inputs for {n} points loaded from {cache.name} (no Drive access) ===", flush=True)
        return {
            "lats": lats, "lons": lons, "chelsa_row": chelsa_row, "chelsa_col": chelsa_col, "token": None,
            "elevation": z["elevation"], "slope": z["slope"], "aspect": z["aspect"], "concavity": z["concavity"],
            "z_ref_m": z["z_ref_m"], "row_px": z["row_px"], "col_px": z["col_px"],
            "soil": {k[6:]: z[k] for k in z.files if k.startswith("soil__")}, "valid_soil": z["valid_soil"],
        }

    token = get_access_token()

    print("=== Terrain: elevation/slope/aspect/concavity/coarse-ref-elevation (streamed) ===", flush=True)
    elevation, slope, aspect, concavity, z_ref_m, row_px, col_px = extract_terrain(token, lats, lons)
    print(f"  elevation {elevation.min():.0f}-{elevation.max():.0f} m, "
          f"z_ref-z_cell spread {np.abs(z_ref_m - elevation).max():.0f} m max", flush=True)

    print("=== Soils (streamed, no full download) ===", flush=True)
    clay_pct, sand_pct, soc_g_kg = extract_soils(token, lats, lons)
    valid_soil = ~np.isnan(clay_pct)
    print(f"  {valid_soil.sum()}/{n} points have real soil data", flush=True)
    soil = soil_hydraulic_parameters(sand_pct, clay_pct, soc_g_kg)

    out = {
        "lats": lats, "lons": lons, "chelsa_row": chelsa_row, "chelsa_col": chelsa_col,
        "token": token, "elevation": elevation, "slope": slope, "aspect": aspect,
        "concavity": concavity, "z_ref_m": z_ref_m, "row_px": row_px, "col_px": col_px,
        "soil": soil, "valid_soil": valid_soil,
    }
    _save_npz(cache, {"elevation": elevation, "slope": slope, "aspect": aspect, "concavity": concavity, "z_ref_m": z_ref_m,
                      "row_px": row_px, "col_px": col_px, "valid_soil": valid_soil, **{f"soil__{k}": v for k, v in soil.items()}})
    print(f"=== cached static inputs -> {cache.name} ===", flush=True)
    return out


_CHELSA_CACHE: dict = {}


def _load_chelsa_arrays():
    if not _CHELSA_CACHE:
        _CHELSA_CACHE["tas"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tas_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["tasmax"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tasmax_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["tasmin"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tasmin_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["pr"] = np.load(DATA_DIR / "chelsa" / "CHELSA_pr_2000_2019_armenia_daily.npz")["data"]
    return _CHELSA_CACHE


_ATMOS = {}


def atmosphere_shape():
    """The monthly ISIMIP files (scripts/pull_isimip_atmosphere.py), loaded once per process."""
    if "a" not in _ATMOS:
        _ATMOS["a"] = AtmosphereShape(DATA_DIR / "isimip3b")
    return _ATMOS["a"]


def constant_atmosphere_requested():
    """ANTAR_CONSTANT_ATMOSPHERE=1 reproduces the earlier behaviour (annual-mean wind, radiation and humidity repeated every day)."""
    return os.environ.get("ANTAR_CONSTANT_ATMOSPHERE", "") == "1"


def extract_year_climate_inputs(static, year):
    """ERA5-Land + CHELSA-daily reference + lapse config for one year -- expensive (streamed) but
    genuinely species-independent, so it's extracted exactly once per (static, year) regardless of
    how many real functional groups later consume it. Separated out of
    :func:`compute_forcing_for_year` so a multi-group caller doesn't re-stream ERA5-Land once per
    group for no reason (only the cheap local w_max_mm multiply actually varies by group).

    Fetches a fresh access token rather than reusing ``static["token"]`` -- see
    :func:`compute_forcing_for_year`'s original docstring note on the real OAuth-expiry crash
    this guards against.
    """
    lats, lons, chelsa_row, chelsa_col = static["lats"], static["lons"], static["chelsa_row"], static["chelsa_col"]
    token = get_access_token()
    row_px, col_px = static["row_px"], static["col_px"]

    print(f"=== ERA5-Land {year} (streamed per-tile) ===", flush=True)
    wind10, ssrd, strd, dewpoint_k, pressure_pa = extract_era5land(token, row_px, col_px, year=year)
    valid_era5 = ~np.isnan(wind10)
    u2_m_s = wind_speed_2m(wind10, z_m=10.0)
    dewpoint_c = dewpoint_k - 273.15
    ea_ref_kpa = saturation_vapour_pressure(dewpoint_c)
    pressure_kpa_era5 = pressure_pa / 1000.0

    arrs = _load_chelsa_arrays()
    day0 = datetime.date(1979, 1, 1)
    i0 = (datetime.date(year, 1, 1) - day0).days
    i1 = (datetime.date(year, 12, 31) - day0).days + 1
    dates = [day0 + datetime.timedelta(days=d) for d in range(i0, i1)]
    doy = np.array([d.timetuple().tm_yday for d in dates])
    month = np.array([d.month for d in dates])

    pr = arrs["pr"]
    pr_available = i1 <= pr.shape[0]
    t_mean_ref_all = arrs["tas"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_max_ref_all = arrs["tasmax"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_min_ref_all = arrs["tasmin"][i0:i1, chelsa_row, chelsa_col] - 273.15
    p_ref_all = pr[i0:i1, chelsa_row, chelsa_col] if pr_available else None

    lapse = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_lapse_rate.yaml"))
    month_order = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    gamma_k_per_m = np.array([lapse["temperature"]["gamma_k_per_m_by_month"][m] for m in month_order])
    precip_gradient_per_m = np.array([lapse["precipitation"]["gradient_per_m_by_month"][m] for m in month_order])

    valid_mask = static["valid_soil"] & valid_era5 & pr_available

    # seasonal shape of shortwave, longwave, humidity and wind (annual mean of each cell unchanged); see antar.climate.atmosphere
    if constant_atmosphere_requested():
        print("=== ANTAR_CONSTANT_ATMOSPHERE=1: annual-mean wind/radiation/humidity repeated every day (old behaviour) ===", flush=True)
        atmos_factors = None
    else:
        atm = atmosphere_shape()
        atmos_factors = [atm.baseline_factors(lats[i], lons[i], doy) if valid_mask[i] else None for i in range(len(lats))]

    return {
        "atmos_factors": atmos_factors,
        "doy": doy, "month": month, "valid_mask": valid_mask,
        "ea_ref_kpa": ea_ref_kpa, "u2_m_s": u2_m_s, "ssrd": ssrd, "strd": strd,
        "pressure_kpa_era5": pressure_kpa_era5,
        "t_mean_ref_all": t_mean_ref_all, "t_max_ref_all": t_max_ref_all, "t_min_ref_all": t_min_ref_all,
        "p_ref_all": p_ref_all,
        "gamma_k_per_m": gamma_k_per_m, "precip_gradient_per_m": precip_gradient_per_m,
    }


def _forcing_cell(static, climate, w_max_mm, i):
    """One real cell's CellTopoclimate, given already-extracted static+climate inputs and a
    (possibly group-specific) w_max_mm array. The actual per-cell, per-group computation --
    cheap and local, no network -- shared by every caller below."""
    elevation, slope, aspect = static["elevation"], static["slope"], static["aspect"]
    concavity, z_ref_m, lats = static["concavity"], static["z_ref_m"], static["lats"]
    soil = static["soil"]
    doy, month = climate["doy"], climate["month"]
    n_days = len(doy)
    t_mean_ref_c = climate["t_mean_ref_all"][:, i]
    t_max_ref_c = climate["t_max_ref_all"][:, i]
    t_min_ref_c = climate["t_min_ref_all"][:, i]
    p_ref_mm = climate["p_ref_all"][:, i]

    gamma_of_day = climate["gamma_k_per_m"][month - 1]
    t_mean_c_for_rn = downscale.downscale_temperature(t_mean_ref_c, elevation[i], z_ref_m[i], gamma_of_day)
    af = climate.get("atmos_factors")
    f = af[i] if af is not None else AtmosphereShape.constant_factors(n_days)
    rn_mj_m2 = net_radiation_from_era5(climate["ssrd"][i] * f["rs"], climate["strd"][i] * f["rl"], t_mean_c_for_rn)

    return topoclimate_forcing(
        doy=doy, month=month,
        t_mean_ref_c=t_mean_ref_c, t_max_ref_c=t_max_ref_c, t_min_ref_c=t_min_ref_c,
        p_ref_mm=p_ref_mm, ea_ref_kpa=climate["ea_ref_kpa"][i] * f["ea"],
        u2_m_s=climate["u2_m_s"][i] * f["wind"], rn_mj_m2=rn_mj_m2,
        z_cell_m=elevation[i], z_ref_m=z_ref_m[i], lat_deg=lats[i],
        slope_deg=slope[i], aspect_deg=aspect[i],
        gamma_k_per_m=climate["gamma_k_per_m"], precip_gradient_per_m=climate["precip_gradient_per_m"],
        w_max_mm=w_max_mm[i], theta_sat=soil["theta_sat"][i], psi_sat_mpa=soil["psi_sat_mpa"][i],
        b_clapp_hornberger=soil["b_clapp_hornberger"][i], theta_fc=soil["theta_fc"][i],
        theta_lim=soil["theta_lim"][i], gdd_budburst=GDD_BUDBURST_PLACEHOLDER,
        concavity_index=concavity[i], calm_clear_night_frac=CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER,
        pressure_kpa=climate["pressure_kpa_era5"][i],
    )


def compute_forcing_for_year(static, year, rooting_depth_mm=ROOTING_DEPTH_MM_PLACEHOLDER):
    """Real per-cell CellTopoclimate for one year, given :func:`extract_static_grid_inputs`'s
    output. Separated from terrain/soil extraction so a multi-year caller only pays the
    (cheap) ERA5-Land/CHELSA/topoclimate_forcing cost once per year, not the (streamed,
    slower) terrain/soil cost.

    Single-rooting-depth convenience wrapper (one real value for the whole grid) -- for a real
    per-functional-group run, use :func:`compute_forcing_for_year_multi_group` instead, which
    shares this same expensive ERA5-Land/CHELSA extraction across all groups rather than
    repeating it.
    """
    climate = extract_year_climate_inputs(static, year)
    soil = static["soil"]
    w_max_mm = (soil["theta_fc"] - soil["theta_lim"]) * rooting_depth_mm
    n = len(static["lats"])

    results = []
    for i in range(n):
        if not climate["valid_mask"][i]:
            results.append(None)
            continue
        results.append(_forcing_cell(static, climate, w_max_mm, i))
    n_ok = sum(1 for r in results if r is not None)
    print(f"  {year}: {n_ok}/{n} points produced real forcing output", flush=True)
    return results


def compute_forcing_for_year_multi_group(static, year, rooting_depth_by_group):
    """Real per-cell CellTopoclimate for one year, once per real functional group -- the group-
    aware counterpart to :func:`compute_forcing_for_year`. Streams ERA5-Land and loads CHELSA
    exactly once (shared, species-independent); only the cheap local w_max_mm multiply and the
    final topoclimate_forcing call repeat per group, since w_max_mm = (theta_fc - theta_lim) *
    rooting_depth is the one real quantity in this whole pipeline that actually depends on species
    (via rooting_depth_by_group -- see ROOTING_DEPTH_MM_BY_GROUP).

    Returns ``{group_name: [CellTopoclimate or None, ...]}``.
    """
    climate = extract_year_climate_inputs(static, year)
    soil = static["soil"]
    n = len(static["lats"])

    results_by_group = {}
    for group, rooting_depth_mm in rooting_depth_by_group.items():
        w_max_mm = (soil["theta_fc"] - soil["theta_lim"]) * rooting_depth_mm
        group_results = []
        for i in range(n):
            if not climate["valid_mask"][i]:
                group_results.append(None)
                continue
            group_results.append(_forcing_cell(static, climate, w_max_mm, i))
        n_ok = sum(1 for r in group_results if r is not None)
        print(f"  {year} / {group} (rooting_depth_mm={rooting_depth_mm}): {n_ok}/{n} points produced real forcing output", flush=True)
        results_by_group[group] = group_results
    return results_by_group


def compute_grid_forcing(year=YEAR, rooting_depth_mm=ROOTING_DEPTH_MM_PLACEHOLDER,
                          grid_rows=GRID_ROWS, grid_cols=GRID_COLS):
    """Build the real CellTopoclimate for every valid grid point, for one year, one (generic or
    single real) rooting depth.

    Factored out of :func:`main` so other scripts can get the same real per-cell forcing without
    duplicating the extraction/orchestration logic. ``rooting_depth_mm`` defaults to the generic
    placeholder (this script's own standalone run is species-agnostic by design). For a real
    multi-group run, use :func:`compute_grid_forcing_multi_group` instead -- it shares the
    expensive streamed extraction across all real functional groups rather than repeating it once
    per group the way calling this function in a loop would.

    ``grid_rows``/``grid_cols`` default to the original 80-point validation grid; pass
    ``DENSE_GRID_ROWS``/``DENSE_GRID_COLS`` for the Armenia-only dense grid.

    Returns ``(lats, lons, elevation, results)``; ``results[i]`` is a
    :class:`antar.climate.forcing.CellTopoclimate` or ``None`` where soil/
    ERA5-Land/CHELSA-pr data was unavailable.
    """
    static = extract_static_grid_inputs(grid_rows=grid_rows, grid_cols=grid_cols)
    results = compute_forcing_for_year(static, year, rooting_depth_mm=rooting_depth_mm)
    return static["lats"], static["lons"], static["elevation"], results


def compute_grid_forcing_multi_group(year=YEAR, rooting_depth_by_group=ROOTING_DEPTH_MM_BY_GROUP,
                                      grid_rows=GRID_ROWS, grid_cols=GRID_COLS):
    """Real per-group counterpart to :func:`compute_grid_forcing`: one real, efficient call gets
    every real functional group's own CellTopoclimate (via its own real rooting depth) while
    streaming terrain/soils/ERA5-Land exactly once, not once per group.

    ``grid_rows``/``grid_cols`` default to the original 80-point validation grid; pass
    ``DENSE_GRID_ROWS``/``DENSE_GRID_COLS`` for the Armenia-only dense grid.

    Returns ``(lats, lons, elevation, results_by_group)``; ``results_by_group[group][i]`` is a
    :class:`antar.climate.forcing.CellTopoclimate` or ``None``.
    """
    static = extract_static_grid_inputs(grid_rows=grid_rows, grid_cols=grid_cols)
    results_by_group = compute_forcing_for_year_multi_group(static, year, rooting_depth_by_group)
    return static["lats"], static["lons"], static["elevation"], results_by_group


def main():
    dense = "--dense" in sys.argv
    grid_rows, grid_cols = (DENSE_GRID_ROWS, DENSE_GRID_COLS) if dense else (GRID_ROWS, GRID_COLS)
    out_path = (Path(__file__).resolve().parent.parent / "configs" / "fitted"
                / "topohydro_grid_run_2019_dense.yaml") if dense else OUT_PATH

    lats, lons, elevation, results = compute_grid_forcing(grid_rows=grid_rows, grid_cols=grid_cols)
    n = len(lats)
    n_ok = sum(1 for r in results if r is not None)

    summary = {
        "run_date": datetime.date.today().isoformat(),
        "year": YEAR,
        "grid": DENSE_GRID_ID if dense else "validation_80pt_stride40",
        "n_grid_points": int(n),
        "n_with_real_output": int(n_ok),
        "grid_definition": ("29x36 subsample (every 11th pixel) of CHELSA-daily's 312x396 armenia "
                             "grid" if dense else
                             "8x10 subsample (every 40th pixel) of CHELSA-daily's 312x396 armenia grid"),
        "placeholders": {
            "calm_clear_night_frac": CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER,
            "gdd_budburst": GDD_BUDBURST_PLACEHOLDER,
            "rooting_depth_mm": (f"{ROOTING_DEPTH_MM_PLACEHOLDER} (generic default -- this "
                                  "standalone run is species-agnostic by design; real per-group "
                                  f"values used elsewhere: {ROOTING_DEPTH_MM_BY_GROUP})"),
        },
        "points": [],
    }
    for i, out in enumerate(results):
        row = {"lat": float(lats[i]), "lon": float(lons[i]), "elevation_m": float(elevation[i])}
        if out is None:
            row["status"] = "skipped_missing_soil_or_era5land"
        else:
            row["status"] = "ok"
            row["t_mean_c_annual_mean"] = float(np.mean(out.t_mean_c))
            row["p_mm_annual_sum"] = float(np.sum(out.p_mm))
            row["gdd_cumulative_annual"] = float(out.gdd_cumulative[-1])
            row["late_frost_days"] = int(out.late_frost_days)
            row["growing_season_length_days"] = int(out.growing_season_length_days)
            row["cwd_mm_by_pet_formulation"] = {k: float(v) for k, v in out.cwd_mm.items()}
            row["wsi_by_pet_formulation"] = {k: float(v) for k, v in out.wsi.items()}
            row["psi_soil_mpa_annual_min_by_pet_formulation"] = {
                k: float(np.min(v)) for k, v in out.psi_soil_mpa.items()
            }
        summary["points"].append(row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.dump(summary, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {out_path} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
