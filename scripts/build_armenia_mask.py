"""Build the Armenia mask on the CHELSA-daily native grid (312 x 396 cells over BBOX).

Why this exists (found 2026-10-05, while building the UI): the project's study grid is a plain
rectangle (lon 43.4-46.7, lat 38.8-41.4) laid over Armenia, and Armenia fills only about a third of
it. Measured against the national Ecosystem Map, just 25 of the 80 stride-40 "validation grid" cells
were inside Armenia; 17 fell outside the map's extent and 38 more sat inside its rectangle over
Georgia, Azerbaijan, Turkey, Iran or Nakhchivan. Every earlier mean, validation correlation and the
AEGIS portfolio included that foreign land. This mask is the corrected sampling frame.

Definition: a CHELSA cell is "inside Armenia" if the Ecosystem Map of Armenia (10 m national
classification, EPSG:32638) has a mapped class (> 0) at the cell centre. Class 0 is the map's own
"outside mapped classes" background (its README says so); water bodies (class 1001, e.g. Lake Sevan)
count as inside the national territory. The raster is read at 1/16 resolution (160 m, nearest) --
ample for ~700-900 m cells, and it keeps the build to seconds.

Output: configs/armenia_mask_chelsa_312x396.npz (committed: tiny, and it lets the pipeline run
without the 46 MB raster present).
"""
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import transform

ROOT = Path(__file__).resolve().parent.parent
RASTER = ROOT / "data" / "ecosystem_map" / "Ecosystem_Map_of_Armenia.tif"
OUT = ROOT / "configs" / "armenia_mask_chelsa_312x396.npz"
BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max
H, W = 312, 396
DECIMATE = 16


def main():
    rows, cols = np.mgrid[0:H, 0:W]
    lat = BBOX[3] - (rows + 0.5) * (BBOX[3] - BBOX[1]) / H
    lon = BBOX[0] + (cols + 0.5) * (BBOX[2] - BBOX[0]) / W
    with rasterio.open(RASTER) as src:
        xs, ys = transform("EPSG:4326", src.crs, lon.ravel().tolist(), lat.ravel().tolist())
        small = src.read(1, out_shape=(src.height // DECIMATE, src.width // DECIMATE), resampling=Resampling.nearest)
        res = src.res[0] * DECIMATE
        left, top = src.bounds.left, src.bounds.top
    c = np.floor((np.array(xs) - left) / res).astype(int)
    r = np.floor((top - np.array(ys)) / res).astype(int)
    in_extent = (r >= 0) & (r < small.shape[0]) & (c >= 0) & (c < small.shape[1])
    cls = np.zeros(H * W, dtype=np.int32)
    cls[in_extent] = small[r[in_extent], c[in_extent]]
    inside = (cls > 0).reshape(H, W)
    cell_km2 = (111.195 * (BBOX[3] - BBOX[1]) / H) * (111.195 * np.cos(np.radians(lat)) * (BBOX[2] - BBOX[0]) / W)
    area = float(cell_km2[inside].sum())
    print(f"inside Armenia: {int(inside.sum())}/{H * W} cells ({100 * inside.mean():.1f}% of the rectangle); "
          f"area ≈ {area:,.0f} km² (Armenia's official area is 29,743 km²)")
    np.savez_compressed(OUT, inside=inside, bbox=np.array(BBOX), shape=np.array([H, W]),
                        source=np.array("Ecosystem Map of Armenia (class > 0 at cell centre), 160 m nearest"))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
