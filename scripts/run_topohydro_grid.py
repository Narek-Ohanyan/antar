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
* ``rooting_depth_mm`` (1000 mm) -- feeds w_max_mm = (theta_fc - theta_lim) *
  rooting_depth; no XYLEM rooting-depth trait exists yet, so a generic
  temperate-tree order-of-magnitude value stands in (Canadell et al. 1996
  put median tree rooting depth around 1-2 m).

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
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.climate import downscale
from antar.climate.forcing import topoclimate_forcing
from antar.climate.radiation import net_radiation_from_era5
from antar.climate.soil_pedotransfer import soil_hydraulic_parameters
from antar.climate.terrain import concavity_index
from antar.climate.vapour import saturation_vapour_pressure, wind_speed_2m

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_grid_run_2019.yaml"

BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max
CHELSA_GRID_SHAPE = (312, 396)
GRID_ROWS = np.arange(0, 312, 40)   # 8 rows
GRID_COLS = np.arange(0, 396, 40)   # 10 cols
YEAR = 2019

TERRAIN_DRIVE_ID = "17zOkIKhKZiDQaRtF2SwPhk3XOa37hpbe"
SOILS_DRIVE_ID = "18LRnI4Nsnlj6ks1UClGDe4BdufqUG24n"
ERA5LAND_TILE_SIZE_PX = 3072  # GEE shard size confirmed against terrain.tif's own pixel grid

CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER = 0.3
GDD_BUDBURST_PLACEHOLDER = 200.0
ROOTING_DEPTH_MM_PLACEHOLDER = 1000.0

COARSE_REF_WINDOW_PX = 15  # +-15 px at 30m = 930m, matching CHELSA-daily's own ~927m pixel
CONCAVITY_WINDOW_PX = 3    # 7x7 (210m) local neighbourhood for cold-air-pooling concavity

ERA5LAND_TILE_IDS = {
    (0, 0): "1kQxh_EDcCn6Ih5xyfvfIZuQE2mWNPqaS", (0, 3072): "1ZCNxrs6P3-_Y8DYg-MetdslHjQnV0phA",
    (0, 6144): "13YdynAAjD0KwDYzpMt3C5GmnH3eKEDHP", (0, 9216): "1rm4-5kX8tO5i855b3BGEP7cXqbb8htzT",
    (3072, 0): "1psHQxLuxdRcN1Oczv7VDh1mcB8a7N38s", (3072, 3072): "1LkBgsLjeZ1_uRlO1c6FHmsZ8pxpIPjnT",
    (3072, 6144): "1FE5mcNk6xONJ4sChXZ7qfmgWyS_45Uk1", (3072, 9216): "1C5RBC59v1ytFwsTWcCme0VM-VrG8wORe",
    (6144, 0): "1wtVOpllg9Q-PmQDaJeetU8ZIbKI1pjMr", (6144, 3072): "1zP8XI026RlIaFCR-3LnF1qxdbQSORnDf",
    (6144, 6144): "1_8tMT-0SnBrRUepjOcWzIo4hjcQrAtcH", (6144, 9216): "1dp-96zRTfqwPDTTdn95tv1-yjmjHJDUp",
    (9216, 0): "1t-hg8gk3sKbWCSvFi7uom7PjGI5SNCwX", (9216, 3072): "1q6KXWJQAIOQgP-OGLtyTS7L3G7YIOvpk",
    (9216, 6144): "1qxYjAt8IwXexaStDAJCphgYXSCDxUzH_", (9216, 9216): "1sJHL1cLIxn9mSwpL9FrPXQBxl0LO1yFT",
}


def grid_latlon():
    lats = BBOX[3] - (GRID_ROWS + 0.5) * (BBOX[3] - BBOX[1]) / CHELSA_GRID_SHAPE[0]
    lons = BBOX[0] + (GRID_COLS + 0.5) * (BBOX[2] - BBOX[0]) / CHELSA_GRID_SHAPE[1]
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    chelsa_row, chelsa_col = np.meshgrid(GRID_ROWS, GRID_COLS, indexing="ij")
    return lat_grid.ravel(), lon_grid.ravel(), chelsa_row.ravel(), chelsa_col.ravel()


def get_access_token():
    import ee
    from google.oauth2.credentials import Credentials
    import google.auth.transport.requests as gareq

    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    creds = Credentials(
        None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
        client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"],
    )
    creds.refresh(gareq.Request())
    return creds.token


def drive_vsicurl_url(file_id: str) -> str:
    return f"/vsicurl/https://www.googleapis.com/drive/v3/files/{file_id}?alt=media"


def extract_terrain(token, lats, lons):
    url = drive_vsicurl_url(TERRAIN_DRIVE_ID)
    with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
        with rasterio.open(url) as src:
            xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
            n = len(xs)
            elevation, slope, aspect, concavity, z_ref = (np.zeros(n) for _ in range(5))
            row_px, col_px = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
            for i, (x, y) in enumerate(zip(xs, ys)):
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
    return elevation, slope, aspect, concavity, z_ref, row_px, col_px


def extract_soils(token, lats, lons):
    url = drive_vsicurl_url(SOILS_DRIVE_ID)
    with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
        with rasterio.open(url) as src:
            band_names = list(src.descriptions)
            xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
            vals = np.array(list(src.sample(zip(xs, ys))))
    idx = {b: i for i, b in enumerate(band_names)}
    # SoilGrids' GEE-mapped units are per-mille (g/kg, 0-1000): verified empirically this
    # session -- clay+sand+silt sum to ~1000 at every real sample point -- so /10 gives the
    # percent (clay/sand) and g/kg (soc) units antar.climate.soil_pedotransfer expects.
    clay_pct = vals[:, idx["clay_0_30cm_mean"]] / 10.0
    sand_pct = vals[:, idx["sand_0_30cm_mean"]] / 10.0
    soc_g_kg = vals[:, idx["soc_0_30cm_mean"]] / 10.0
    return clay_pct, sand_pct, soc_g_kg


def extract_era5land(token, row_px, col_px, year=YEAR):
    n = len(row_px)
    wind10, ssrd, strd, dewpoint_k, pressure_pa = (np.full(n, np.nan) for _ in range(5))
    tile_row = (row_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX
    tile_col = (col_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX

    for tile_key in sorted(set(zip(tile_row.tolist(), tile_col.tolist()))):
        file_id = ERA5LAND_TILE_IDS[tile_key]
        sel = (tile_row == tile_key[0]) & (tile_col == tile_key[1])
        url = drive_vsicurl_url(file_id)
        with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
            with rasterio.open(url) as src:
                bidx = {b: j + 1 for j, b in enumerate(src.descriptions)}
                for i in np.where(sel)[0]:
                    lr, lc = int(row_px[i] - tile_key[0]), int(col_px[i] - tile_key[1])
                    w = Window(lc, lr, 1, 1)
                    wind10[i] = src.read(bidx[f"wind_speed_{year}"], window=w)[0, 0]
                    ssrd[i] = src.read(bidx[f"ssrd_{year}"], window=w)[0, 0]
                    strd[i] = src.read(bidx[f"strd_{year}"], window=w)[0, 0]
                    dewpoint_k[i] = src.read(bidx[f"dewpoint_{year}"], window=w)[0, 0]
                    pressure_pa[i] = src.read(bidx[f"surface_pressure_{year}"], window=w)[0, 0]
        print(f"  tile {tile_key}: {sel.sum()} points", flush=True)
    return wind10, ssrd, strd, dewpoint_k, pressure_pa


def extract_static_grid_inputs():
    """Terrain and soils: year-independent, extracted once and reused across years.

    Factored out so a multi-year caller (MNEME's hazard-panel build) pays for
    terrain/soil extraction once, not once per year.
    """
    lats, lons, chelsa_row, chelsa_col = grid_latlon()
    n = len(lats)
    print(f"=== {n} grid points ({len(GRID_ROWS)}x{len(GRID_COLS)}) ===", flush=True)

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
    w_max_mm = (soil["theta_fc"] - soil["theta_lim"]) * ROOTING_DEPTH_MM_PLACEHOLDER

    return {
        "lats": lats, "lons": lons, "chelsa_row": chelsa_row, "chelsa_col": chelsa_col,
        "token": token, "elevation": elevation, "slope": slope, "aspect": aspect,
        "concavity": concavity, "z_ref_m": z_ref_m, "row_px": row_px, "col_px": col_px,
        "soil": soil, "w_max_mm": w_max_mm, "valid_soil": valid_soil,
    }


_CHELSA_CACHE: dict = {}


def _load_chelsa_arrays():
    if not _CHELSA_CACHE:
        _CHELSA_CACHE["tas"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tas_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["tasmax"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tasmax_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["tasmin"] = np.load(DATA_DIR / "chelsa" / "CHELSA_tasmin_2000_2024_armenia_daily.npz")["data"]
        _CHELSA_CACHE["pr"] = np.load(DATA_DIR / "chelsa" / "CHELSA_pr_2000_2019_armenia_daily.npz")["data"]
    return _CHELSA_CACHE


def compute_forcing_for_year(static, year):
    """Real per-cell CellTopoclimate for one year, given :func:`extract_static_grid_inputs`'s
    output. Separated from terrain/soil extraction so a multi-year caller only pays the
    (cheap) ERA5-Land/CHELSA/topoclimate_forcing cost once per year, not the (streamed,
    slower) terrain/soil cost.

    Fetches a fresh access token on every call rather than reusing ``static["token"]`` --
    Google OAuth access tokens expire after ~1 hour, and a multi-year caller (MNEME's hazard
    panel) runs for several hours; reusing one token caused a real crash (HTTP 400 on year 8
    of 10, `scripts/fit_mneme_hazard_panel.py`, 2026-09-30) that lost the whole run since it
    checkpointed nothing. Refreshing per year is cheap (one token exchange) next to the
    minutes of streaming each year already costs.
    """
    lats, lons, chelsa_row, chelsa_col = static["lats"], static["lons"], static["chelsa_row"], static["chelsa_col"]
    token = get_access_token()
    elevation, slope, aspect = static["elevation"], static["slope"], static["aspect"]
    concavity, z_ref_m, row_px, col_px = static["concavity"], static["z_ref_m"], static["row_px"], static["col_px"]
    soil, w_max_mm, valid_soil = static["soil"], static["w_max_mm"], static["valid_soil"]
    n = len(lats)

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
    n_days = len(dates)

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

    results = []
    valid_mask = valid_soil & valid_era5 & pr_available
    for i in range(n):
        if not valid_mask[i]:
            results.append(None)
            continue
        t_mean_ref_c, t_max_ref_c, t_min_ref_c = t_mean_ref_all[:, i], t_max_ref_all[:, i], t_min_ref_all[:, i]
        p_ref_mm = p_ref_all[:, i]

        gamma_of_day = gamma_k_per_m[month - 1]
        t_mean_c_for_rn = downscale.downscale_temperature(t_mean_ref_c, elevation[i], z_ref_m[i], gamma_of_day)
        rn_mj_m2 = net_radiation_from_era5(np.full(n_days, ssrd[i]), np.full(n_days, strd[i]), t_mean_c_for_rn)

        out = topoclimate_forcing(
            doy=doy, month=month,
            t_mean_ref_c=t_mean_ref_c, t_max_ref_c=t_max_ref_c, t_min_ref_c=t_min_ref_c,
            p_ref_mm=p_ref_mm, ea_ref_kpa=np.full(n_days, ea_ref_kpa[i]),
            u2_m_s=np.full(n_days, u2_m_s[i]), rn_mj_m2=rn_mj_m2,
            z_cell_m=elevation[i], z_ref_m=z_ref_m[i], lat_deg=lats[i],
            slope_deg=slope[i], aspect_deg=aspect[i],
            gamma_k_per_m=gamma_k_per_m, precip_gradient_per_m=precip_gradient_per_m,
            w_max_mm=w_max_mm[i], theta_sat=soil["theta_sat"][i], psi_sat_mpa=soil["psi_sat_mpa"][i],
            b_clapp_hornberger=soil["b_clapp_hornberger"][i], theta_fc=soil["theta_fc"][i],
            theta_lim=soil["theta_lim"][i], gdd_budburst=GDD_BUDBURST_PLACEHOLDER,
            concavity_index=concavity[i], calm_clear_night_frac=CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER,
            pressure_kpa=pressure_kpa_era5[i],
        )
        results.append(out)
    n_ok = sum(1 for r in results if r is not None)
    print(f"  {year}: {n_ok}/{n} points produced real forcing output", flush=True)
    return results


def compute_grid_forcing(year=YEAR):
    """Build the real CellTopoclimate for every valid grid point, for one year.

    Factored out of :func:`main` so other scripts (XYLEM's mechanistic-hazard
    re-fit, MNEME's hazard-panel build) can get the same real per-cell forcing
    without duplicating the extraction/orchestration logic -- re-running this
    (a few minutes, all transient/streamed, nothing extra cached) rather than
    trying to serialize full daily CellTopoclimate objects to disk.

    Returns ``(lats, lons, elevation, results)``; ``results[i]`` is a
    :class:`antar.climate.forcing.CellTopoclimate` or ``None`` where soil/
    ERA5-Land/CHELSA-pr data was unavailable.
    """
    static = extract_static_grid_inputs()
    results = compute_forcing_for_year(static, year)
    return static["lats"], static["lons"], static["elevation"], results


def main():
    lats, lons, elevation, results = compute_grid_forcing()
    n = len(lats)
    n_ok = sum(1 for r in results if r is not None)

    summary = {
        "run_date": datetime.date.today().isoformat(),
        "year": YEAR,
        "n_grid_points": int(n),
        "n_with_real_output": int(n_ok),
        "grid_definition": "8x10 subsample (every 40th pixel) of CHELSA-daily's 312x396 armenia grid",
        "placeholders": {
            "calm_clear_night_frac": CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER,
            "gdd_budburst": GDD_BUDBURST_PLACEHOLDER,
            "rooting_depth_mm": ROOTING_DEPTH_MM_PLACEHOLDER,
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

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.dump(summary, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
