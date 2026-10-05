"""Elevation on the map's display grid, so interpolated surfaces can follow terrain.

Streams the project's 30 m terrain raster (band 1, elevation) from Drive and averages it onto the same
Web-Mercator grid as ui/assets/map (grid.json). Nothing is downloaded to disk but the small derived
raster: ui/assets/map/elevation.png, an 8-bit RGB PNG with elevation in metres = 256*R + G (lossless;
a canvas cannot read 16-bit PNGs), and elevation.json (source, resampling, and a check against the
elevations the model was actually run at).

Why: plain distance-based interpolation cannot reproduce quantities that are driven by elevation
(leave-one-out R^2 for annual temperature on the 25-node grid was -0.23). With elevation as a covariate,
temperature follows the fitted lapse rate between nodes.
"""
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import run_topohydro_grid as R  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ROOT / "ui" / "assets" / "map"


def main():
    grid = json.loads((ASSETS / "grid.json").read_text())
    w, h, px = grid["width"], grid["height"], grid["px_m"]
    dst_tf = from_origin(grid["x0"], grid["y_top"], px, px)
    region = np.array(Image.open(ASSETS / "region.png"))
    dem = np.full((h, w), np.nan, dtype=np.float32)

    token = R.get_access_token()
    url = R.drive_vsicurl_url(R.TERRAIN_DRIVE_ID)
    with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES",
                      GDAL_HTTP_TIMEOUT=60, GDAL_HTTP_CONNECTTIMEOUT=10, GDAL_CACHEMAX=512):
        with rasterio.open(url) as src:
            print(f"source {src.shape} {src.res} {src.crs}; averaging onto {w}x{h} at {px:.0f} m", flush=True)
            reproject(rasterio.band(src, 1), dem, dst_transform=dst_tf, dst_crs="EPSG:3857",
                      resampling=Resampling.average, src_nodata=src.nodata, dst_nodata=np.nan)

    inside = region > 0
    covered = np.isfinite(dem) & inside
    miss = float((inside & ~covered).sum() / inside.sum())
    print(f"inside-country pixels without elevation: {100 * miss:.2f}%", flush=True)
    if miss > 0.02:
        raise SystemExit("more than 2% of Armenian pixels have no elevation -- refusing to write a degraded raster")
    # the few uncovered pixels (edge of the source) take the nearest covered value
    from scipy import ndimage
    idx = ndimage.distance_transform_edt(~covered, return_distances=False, return_indices=True)
    filled = dem[idx[0], idx[1]]
    z = np.clip(np.round(np.where(inside, filled, 0)), 0, 65535).astype(np.uint16)
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[..., 0], rgb[..., 1] = z >> 8, z & 255
    Image.fromarray(rgb, mode="RGB").save(ASSETS / "elevation.png", optimize=True)

    # check against the elevations the model itself used (30 m point values at the node locations)
    checks = {}
    for f in sorted((ROOT / "ui" / "data").glob("grid_*.json")):
        g = json.loads(f.read_text())
        errs = []
        for la, lo, e in zip(g["lat"], g["lon"], g["elev"]):
            if e is None:
                continue
            x, y = R_merc(lo, lat=la)
            col, row = int((x - grid["x0"]) // px), int((grid["y_top"] - y) // px)
            if 0 <= row < h and 0 <= col < w and inside[row, col]:
                errs.append(float(z[row, col]) - e)
        if errs:
            checks[g["grid"]] = {"n": len(errs), "mean_diff_m": round(float(np.mean(errs)), 1), "rmse_m": round(float(np.sqrt(np.mean(np.square(errs)))), 1)}
    info = {"built": datetime.date.today().isoformat(), "source": "project terrain raster (SRTM-derived, 30 m, EPSG:32638), band 1",
            "resampling": f"average over each {px:.0f} m Mercator pixel (~{px * 0.766:.0f} m on the ground)",
            "encoding": "elevation_m = 256*R + G", "min_m": int(z[inside].min()), "max_m": int(z[inside].max()),
            "mean_m": round(float(z[inside].mean()), 1), "node_elevation_check": checks}
    (ASSETS / "elevation.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


def R_merc(lon, lat):
    r = 6378137.0
    return r * np.radians(lon), r * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


if __name__ == "__main__":
    sys.exit(main())
