"""Replace MERISTEM's `cwd_approx_mm = petmean - bio12` proxy with real TOPOHYDRO water-balance
CWD, at the exact 3524 real points (524 presence + 3000 shared background) the adult-niche fit
actually trains on -- not just the 78-point validation grid used everywhere else tonight.

Scale forced a real engineering choice, stated because it isn't obvious: these 3524 real points
(cached in data/_cache_niche_features.npz from the original fit) span almost the full country
(lat 38.89-41.21, lon 43.56-46.53 -- essentially the whole real BBOX), so there is no smaller
sub-region to shrink the problem to. Terrain and soils are downloaded fully, locally, once
(~840MB transient, matching this session's original pre-streaming convention -- defensible here
specifically because thousands of scattered local-window reads (concavity, the coarse z_ref_m
window) via individual network range requests would be far slower and noisier than one bulk
download followed by fast local numpy operations, deleted immediately after). ERA5-Land stays
streamed (16 tiles, 15+ GB combined -- a full download would be a real, unjustified storage
regression) but reads are batched per tile via vectorized rasterio .sample() calls, not a
per-point loop -- the same fix that made the 80-point grid runs fast, now necessary rather than
optional at 44x the point count.

**Real per-functional-group CWD, not one flat value for all 7 species.** MERISTEM fits a
presence-background model per species (not pooled), sharing the same 3000 background points across
every species' fit. w_max_mm (soil water-holding capacity) depends on rooting depth, which is real
and species-specific (``ROOTING_DEPTH_MM_BY_GROUP``, Canadell et al. 1996) -- so CWD is computed
once per real functional group, for every point (including background, since a background point
needs each group's own real CWD when it's used in that group's own fit), not once globally. Only
the cheap local `topoclimate_forcing` loop repeats per group; the expensive terrain/soils download
and ERA5-Land streaming above are each done exactly once, shared across all 4 groups -- the same
efficient split used in run_topohydro_grid.py/fit_xylem_mechanistic_hazard.py/
fit_refugium_viability.py/run_future_projections.py.
"""
import datetime
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    DATA_DIR, get_access_token, drive_vsicurl_url, _load_chelsa_arrays,
    saturation_vapour_pressure, wind_speed_2m, TERRAIN_DRIVE_ID, SOILS_DRIVE_ID,
    ERA5LAND_TILE_IDS, ERA5LAND_TILE_SIZE_PX, CONCAVITY_WINDOW_PX, COARSE_REF_WINDOW_PX,
    ROOTING_DEPTH_MM_BY_GROUP,
)
from antar.climate import downscale  # noqa: E402
from antar.climate.forcing import topoclimate_forcing  # noqa: E402
from antar.climate.radiation import net_radiation_from_era5  # noqa: E402
from antar.climate.soil_pedotransfer import soil_hydraulic_parameters  # noqa: E402
from antar.climate.terrain import concavity_index  # noqa: E402

CACHE_PATH = DATA_DIR / "_cache_niche_features.npz"
OUT_PATH = DATA_DIR / "_real_cwd_for_meristem.npz"
YEAR = 2019
GDD_BUDBURST_PLACEHOLDER = 200.0
CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER = 0.3
BBOX = (43.4, 38.8, 46.7, 41.4)

# Real species-to-functional-group mapping, from configs/species_traits.csv's own `example_taxa`
# column (not re-derived/guessed here). "background" points are shared across every species'
# presence-background fit (MERISTEM fits per species, not pooled), so for each real functional
# group a background point needs that group's own CWD too, computed uniformly for every point --
# see main()'s per-group loop below.
SPECIES_TO_FUNCTIONAL_GROUP = {
    "Fagus orientalis": "mesic_diffuse_porous_broadleaf",
    "Carpinus betulus": "mesic_diffuse_porous_broadleaf",
    "Quercus macranthera": "ring_porous_oak",
    "Quercus iberica": "ring_porous_oak",
    "Pinus kochiana": "pine",
    "Juniperus polycarpos": "juniper_arid_conifer",
    "Juniperus excelsa": "juniper_arid_conifer",
}


def load_all_points():
    d = np.load(CACHE_PATH, allow_pickle=True)
    groups, all_ll = [], []
    for k in d["keys"]:
        ll = d[f"{k}_ll"]
        groups.extend([k] * len(ll))
        all_ll.append(ll)
    all_ll = np.concatenate(all_ll, axis=0)
    return np.array(groups), all_ll[:, 0], all_ll[:, 1]  # group labels, lats, lons


def _download_full_raster(file_id, local_path, label):
    """Full local download, retried with a fresh token on failure. A real bug found during a
    broader audit (2026-10-02): this used a single token passed in once, never refreshed --
    the exact class of failure already fixed in extract_terrain/extract_soils/extract_era5land
    in run_topohydro_grid.py, just not yet ported here. Returns (data, profile, band_descriptions)."""
    last_err = None
    for attempt in range(3):
        try:
            print(f"  downloading {label} (transient, full, attempt {attempt + 1}/3)...", flush=True)
            fresh_token = get_access_token()
            url = drive_vsicurl_url(file_id)
            with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {fresh_token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                with rasterio.open(url) as src:
                    data = src.read()
                    profile = src.profile
                    band_descriptions = src.descriptions
            return data, profile, band_descriptions
        except rasterio.errors.RasterioIOError as e:
            last_err = e
            print(f"  {label}: download failed (attempt {attempt + 1}/3): {e}", flush=True)
            time.sleep(5)
    raise last_err


def extract_terrain_soils_local(lats, lons, token):
    """Full local download of terrain.tif + soils.tif (transient), fast local numpy extraction
    for all points -- see module docstring for why, at this point count."""
    n = len(lats)
    terrain_path = DATA_DIR / "_tmp_terrain_meristem.tif"
    soils_path = DATA_DIR / "_tmp_soils_meristem.tif"

    data, profile, _ = _download_full_raster(TERRAIN_DRIVE_ID, terrain_path, "terrain.tif")
    with rasterio.open(terrain_path, "w", **profile) as dst:
        dst.write(data)
    del data

    data, profile, band_descriptions = _download_full_raster(SOILS_DRIVE_ID, soils_path, "soils.tif")
    # .profile does NOT include band descriptions (names) -- a real bug this hit: soils needs
    # band names (clay_0_30cm_mean etc.) for extraction, unlike terrain which reads by index.
    # Copy descriptions onto the local copy explicitly, don't assume .profile carries them.
    with rasterio.open(soils_path, "w", **profile) as dst:
        dst.write(data)
        dst.descriptions = band_descriptions
    del data

    elevation, slope, aspect, concavity, z_ref = (np.zeros(n) for _ in range(5))
    row_px, col_px = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
    with rasterio.open(terrain_path) as src:
        xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
        for i, (x, y) in enumerate(zip(xs, ys)):
            row, col = src.index(x, y)
            row_px[i], col_px[i] = row, col
            cw = CONCAVITY_WINDOW_PX
            r0, r1 = max(0, row - cw), min(src.height, row + cw + 1)
            c0, c1 = max(0, col - cw), min(src.width, col + cw + 1)
            w = Window(c0, r0, c1 - c0, r1 - r0)
            elev_win, slope_win, aspect_win = src.read(1, window=w), src.read(2, window=w), src.read(3, window=w)
            lr, lc = row - r0, col - c0
            elevation[i], slope[i], aspect[i] = elev_win[lr, lc], slope_win[lr, lc], aspect_win[lr, lc]
            concavity[i] = concavity_index(elev_win)[lr, lc]
            rw = COARSE_REF_WINDOW_PX
            r0c, r1c = max(0, row - rw), min(src.height, row + rw + 1)
            c0c, c1c = max(0, col - rw), min(src.width, col + rw + 1)
            z_ref[i] = float(np.nanmean(src.read(1, window=Window(c0c, r0c, c1c - c0c, r1c - r0c))))
            if (i + 1) % 500 == 0:
                print(f"    terrain: {i + 1}/{n}", flush=True)

    with rasterio.open(soils_path) as src:
        band_names = list(src.descriptions)
        xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
        soil_vals = np.array(list(src.sample(zip(xs, ys))))
    idx = {b: j for j, b in enumerate(band_names)}
    clay_pct = soil_vals[:, idx["clay_0_30cm_mean"]] / 10.0
    sand_pct = soil_vals[:, idx["sand_0_30cm_mean"]] / 10.0
    soc_g_kg = soil_vals[:, idx["soc_0_30cm_mean"]] / 10.0

    terrain_path.unlink()
    soils_path.unlink()
    return elevation, slope, aspect, concavity, z_ref, row_px, col_px, clay_pct, sand_pct, soc_g_kg


def extract_era5land_vectorized(token, row_px, col_px, year=YEAR):
    n = len(row_px)
    wind10, ssrd, strd, dewpoint_k, pressure_pa = (np.full(n, np.nan) for _ in range(5))
    tile_row = (row_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX
    tile_col = (col_px // ERA5LAND_TILE_SIZE_PX) * ERA5LAND_TILE_SIZE_PX

    failed_tiles = []
    for tile_key in sorted(set(zip(tile_row.tolist(), tile_col.tolist()))):
        if tile_key not in ERA5LAND_TILE_IDS:
            continue
        file_id = ERA5LAND_TILE_IDS[tile_key]
        sel = np.where((tile_row == tile_key[0]) & (tile_col == tile_key[1]))[0]
        local_rows = row_px[sel] - tile_key[0]
        local_cols = col_px[sel] - tile_key[1]

        # Streaming individual Drive-hosted tiles over vsicurl occasionally hits a transient
        # read failure under sustained multi-tile load (real fault hit mid-run: TIFFReadEncodedTile
        # "got 0 bytes, expected N" on one tile after 3 others streamed fine). Retry with a fresh
        # token rather than letting one flaky tile abort all 16; only skip (leave NaN) after 3
        # failures, which the existing valid_era5 mask already handles honestly downstream.
        last_err = None
        for attempt in range(3):
            try:
                fresh_token = get_access_token()
                url = drive_vsicurl_url(file_id)
                with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {fresh_token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                    with rasterio.open(url) as src:
                        bidx = {b: j + 1 for j, b in enumerate(src.descriptions)}
                        xs_geo, ys_geo = src.xy(local_rows, local_cols)
                        coords = list(zip(xs_geo, ys_geo))
                        tile_vals = {}
                        for bandname in ["wind_speed", "ssrd", "strd", "dewpoint", "surface_pressure"]:
                            vals = np.array(list(src.sample(coords, indexes=bidx[f"{bandname}_{year}"])))[:, 0]
                            tile_vals[bandname] = vals
                for bandname, out in [("wind_speed", wind10), ("ssrd", ssrd), ("strd", strd),
                                       ("dewpoint", dewpoint_k), ("surface_pressure", pressure_pa)]:
                    out[sel] = tile_vals[bandname]
                last_err = None
                break
            except rasterio.errors.RasterioIOError as e:
                last_err = e
                print(f"  era5land tile {tile_key}: read failed (attempt {attempt + 1}/3): {e}", flush=True)
                time.sleep(5)
        if last_err is not None:
            print(f"  era5land tile {tile_key}: FAILED after 3 attempts, {len(sel)} points left NaN", flush=True)
            failed_tiles.append(tile_key)
            continue
        print(f"  era5land tile {tile_key}: {len(sel)} points", flush=True)
    if failed_tiles:
        print(f"  === {len(failed_tiles)} tile(s) failed all retries: {failed_tiles} ===", flush=True)
    return wind10, ssrd, strd, dewpoint_k, pressure_pa


def main():
    print("=== Loading the real 3524 MERISTEM points (presence + background) ===", flush=True)
    groups, lats, lons = load_all_points()
    n = len(lats)
    print(f"  {n} real points, lat {lats.min():.2f}-{lats.max():.2f}, lon {lons.min():.2f}-{lons.max():.2f}", flush=True)

    token = get_access_token()
    print("=== Terrain + soils (local, transient full download, fast local extraction) ===", flush=True)
    elevation, slope, aspect, concavity, z_ref_m, row_px, col_px, clay_pct, sand_pct, soc_g_kg = \
        extract_terrain_soils_local(lats, lons, token)
    valid_soil = ~np.isnan(clay_pct)
    print(f"  {valid_soil.sum()}/{n} points have real soil data", flush=True)
    soil = soil_hydraulic_parameters(sand_pct, clay_pct, soc_g_kg)
    # Real per-group w_max_mm (cheap, local): each real functional group's own rooting depth,
    # not one flat value shared by every species -- see run_topohydro_grid.py's module docstring.
    w_max_mm_by_group = {g: (soil["theta_fc"] - soil["theta_lim"]) * rd
                          for g, rd in ROOTING_DEPTH_MM_BY_GROUP.items()}

    print("=== ERA5-Land 2019 (streamed, vectorized per-tile) ===", flush=True)
    token = get_access_token()  # refresh -- the terrain/soils download above can take a while
    wind10, ssrd, strd, dewpoint_k, pressure_pa = extract_era5land_vectorized(token, row_px, col_px)
    valid_era5 = np.isfinite(wind10) & np.isfinite(ssrd) & np.isfinite(strd) & np.isfinite(dewpoint_k) & np.isfinite(pressure_pa)
    print(f"  {valid_era5.sum()}/{n} points have real ERA5-Land values", flush=True)
    u2_m_s = wind_speed_2m(wind10, z_m=10.0)
    ea_ref_kpa = saturation_vapour_pressure(dewpoint_k - 273.15)
    pressure_kpa_era5 = pressure_pa / 1000.0

    print("=== CHELSA-daily 2019 reference series (local, fast) ===", flush=True)
    arrs = _load_chelsa_arrays()
    day0 = datetime.date(1979, 1, 1)
    i0 = (datetime.date(YEAR, 1, 1) - day0).days
    i1 = (datetime.date(YEAR, 12, 31) - day0).days + 1
    dates = [day0 + datetime.timedelta(days=d) for d in range(i0, i1)]
    doy = np.array([d.timetuple().tm_yday for d in dates])
    month = np.array([d.month for d in dates])

    pixel_size = (BBOX[3] - BBOX[1]) / 312
    chelsa_row = np.clip(np.floor((BBOX[3] - lats) / pixel_size).astype(int), 0, 311)
    chelsa_col = np.clip(np.floor((lons - BBOX[0]) / pixel_size).astype(int), 0, 395)
    t_mean_ref_all = arrs["tas"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_max_ref_all = arrs["tasmax"][i0:i1, chelsa_row, chelsa_col] - 273.15
    t_min_ref_all = arrs["tasmin"][i0:i1, chelsa_row, chelsa_col] - 273.15
    p_ref_all = arrs["pr"][i0:i1, chelsa_row, chelsa_col]

    lapse_cfg = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_lapse_rate.yaml"))
    month_order = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    gamma_k_per_m = np.array([lapse_cfg["temperature"]["gamma_k_per_m_by_month"][m] for m in month_order])
    precip_gradient_per_m = np.array([lapse_cfg["precipitation"]["gradient_per_m_by_month"][m] for m in month_order])

    print("=== Running topoclimate_forcing per point, once per real functional group (real CWD) ===", flush=True)
    valid_mask = valid_soil & valid_era5
    real_cwd_mm_by_group = {}
    for gname, w_max_mm in w_max_mm_by_group.items():
        print(f"  === {gname} (rooting_depth_mm={ROOTING_DEPTH_MM_BY_GROUP[gname]}) ===", flush=True)
        real_cwd_mm = np.full(n, np.nan)
        for i in range(n):
            if not valid_mask[i]:
                continue
            gamma_of_day = gamma_k_per_m[month - 1]
            t_mean_c_for_rn = downscale.downscale_temperature(t_mean_ref_all[:, i], elevation[i], z_ref_m[i], gamma_of_day)
            rn_mj_m2 = net_radiation_from_era5(np.full(len(doy), ssrd[i]), np.full(len(doy), strd[i]), t_mean_c_for_rn)
            out = topoclimate_forcing(
                doy=doy, month=month,
                t_mean_ref_c=t_mean_ref_all[:, i], t_max_ref_c=t_max_ref_all[:, i], t_min_ref_c=t_min_ref_all[:, i],
                p_ref_mm=p_ref_all[:, i], ea_ref_kpa=np.full(len(doy), ea_ref_kpa[i]),
                u2_m_s=np.full(len(doy), u2_m_s[i]), rn_mj_m2=rn_mj_m2,
                z_cell_m=elevation[i], z_ref_m=z_ref_m[i], lat_deg=lats[i],
                slope_deg=slope[i], aspect_deg=aspect[i],
                gamma_k_per_m=gamma_k_per_m, precip_gradient_per_m=precip_gradient_per_m,
                w_max_mm=w_max_mm[i], theta_sat=soil["theta_sat"][i], psi_sat_mpa=soil["psi_sat_mpa"][i],
                b_clapp_hornberger=soil["b_clapp_hornberger"][i], theta_fc=soil["theta_fc"][i],
                theta_lim=soil["theta_lim"][i], gdd_budburst=GDD_BUDBURST_PLACEHOLDER,
                concavity_index=concavity[i], calm_clear_night_frac=CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER,
                pressure_kpa=pressure_kpa_era5[i],
            )
            real_cwd_mm[i] = out.cwd_mm["pm_fao56"]
            if (i + 1) % 500 == 0:
                print(f"    [{i + 1}/{n}] done", flush=True)
        n_ok = int(np.sum(~np.isnan(real_cwd_mm)))
        print(f"  {gname}: {n_ok}/{n} real CWD values computed", flush=True)
        real_cwd_mm_by_group[gname] = real_cwd_mm

    np.savez_compressed(
        OUT_PATH, lats=lats, lons=lons, groups=groups,
        species_to_functional_group=json.dumps(SPECIES_TO_FUNCTIONAL_GROUP),
        **{f"real_cwd_mm__{g}": arr for g, arr in real_cwd_mm_by_group.items()},
    )
    print(f"=== Wrote {OUT_PATH} (real per-group CWD: {list(real_cwd_mm_by_group)}) ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
