"""Static map layers for the UI: current forest cover, water, marz ids and borders, on one
Web-Mercator grid so Leaflet can place them exactly without a basemap.

Sources (see configs/manifests/):
* Ecosystem Map of Armenia (10 m, EPSG:32638): the national outline (class > 0), forest and woodland
  cover, water. It is block-averaged to 500 m cells in its own projection, then reprojected to the
  display grid with bilinear resampling.
* geoBoundaries gbOpen ARM ADM1: marz polygons, used only to give each in-country pixel a marz id.

Outputs (ui/assets/map/): grid.json, region.png, forest.png, woodland.png, forest_<group>.png, water.png
(8-bit, fraction x 255), borders.geojson (smoothed marz borders and country outline).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.features import rasterize, shapes
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
RASTER = ROOT / "data" / "ecosystem_map" / "Ecosystem_Map_of_Armenia.tif"
BOUNDARIES = ROOT / "data" / "boundaries" / "geoBoundaries-ARM-ADM1.geojson"
OUT = ROOT / "ui" / "assets" / "map"

R_EARTH = 6378137.0
PX_M = 650.0                      # Mercator metres per display pixel (~500 m on the ground at 40 N)
LON_W, LON_E, LAT_S, LAT_N = 43.40, 46.70, 38.80, 41.45
BLOCK = 50                        # 10 m source pixels per 500 m block

FOREST_CODES = [31, 32, 33, 34, 35, 36, 37]       # closed forests and closed forest plantations
WOODLAND_CODES = [39, 41, 43, 44]                 # open woodland plantations, subalpine, mixed and juniper woodlands
GROUP_CODES = {"broadleaf": [31], "oak": [32, 33, 34, 35], "pine": [36], "juniper": [44]}
WATER_CODE = 1001
OFFICIAL_AREA_KM2 = 29743.0       # national land area, used as a sanity check only


def merc_x(lon):
    return R_EARTH * np.radians(lon)


def merc_y(lat):
    return R_EARTH * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def inv_merc(x, y):
    return np.degrees(x / R_EARTH), np.degrees(2 * np.arctan(np.exp(y / R_EARTH)) - np.pi / 2)


def chaikin(ring: np.ndarray, iters: int = 2, closed: bool = True) -> np.ndarray:
    """Corner-cutting smoothing; keeps the curve inside the original polygon's hull."""
    pts = np.asarray(ring, dtype=float)
    if closed and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    for _ in range(iters):
        nxt = np.roll(pts, -1, axis=0) if closed else pts[1:]
        cur = pts if closed else pts[:-1]
        q, r = 0.75 * cur + 0.25 * nxt, 0.25 * cur + 0.75 * nxt
        out = np.empty((2 * len(cur), 2))
        out[0::2], out[1::2] = q, r
        pts = out
    return np.vstack([pts, pts[:1]]) if closed else pts


def rdp(pts: np.ndarray, tol: float) -> np.ndarray:
    """Ramer-Douglas-Peucker simplification of an open polyline (iterative, numpy)."""
    pts = np.asarray(pts, dtype=float)
    keep = np.zeros(len(pts), dtype=bool)
    keep[[0, -1]] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        d = pts[b] - pts[a]
        seg = pts[a + 1:b] - pts[a]
        norm = np.hypot(*d)
        dist = np.hypot(*seg.T) if norm == 0 else np.abs(d[0] * seg[:, 1] - d[1] * seg[:, 0]) / norm
        k = int(np.argmax(dist))
        if dist[k] > tol:
            keep[a + 1 + k] = True
            stack += [(a, a + 1 + k), (a + 1 + k, b)]
    return pts[keep]


def block_fractions(src, codes_by_name):
    """Per-500 m block fraction of each named code set, read strip by strip (the raster is ~1.5 GB)."""
    h, w = src.height, src.width
    nh, nw = h // BLOCK, w // BLOCK
    out = {k: np.zeros((nh, nw), dtype=np.float32) for k in list(codes_by_name) + ["inside"]}
    for br in range(nh):
        strip = src.read(1, window=((br * BLOCK, (br + 1) * BLOCK), (0, nw * BLOCK)))
        blocks = strip.reshape(BLOCK, nw, BLOCK).transpose(1, 0, 2).reshape(nw, -1)
        out["inside"][br] = (blocks > 0).mean(axis=1)
        for name, codes in codes_by_name.items():
            out[name][br] = np.isin(blocks, codes).mean(axis=1)
    return out, from_origin(src.bounds.left, src.bounds.top, BLOCK * 10.0, BLOCK * 10.0), src.crs


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    x0, x1 = merc_x(LON_W), merc_x(LON_E)
    y_top, y_bot = merc_y(LAT_N), merc_y(LAT_S)
    width, height = int(np.ceil((x1 - x0) / PX_M)), int(np.ceil((y_top - y_bot) / PX_M))
    dst_tf = from_origin(x0, y_top, PX_M, PX_M)

    codes = {"forest": FOREST_CODES, "woodland": WOODLAND_CODES, "water": [WATER_CODE], **{f"forest_{k}": v for k, v in GROUP_CODES.items()}}
    with rasterio.open(RASTER) as src:
        blk, src_tf, src_crs = block_fractions(src, codes)
    fr = {}
    for name, arr in blk.items():
        dst = np.zeros((height, width), dtype=np.float32)
        reproject(arr, dst, src_transform=src_tf, src_crs=src_crs, dst_transform=dst_tf, dst_crs="EPSG:3857",
                  resampling=Resampling.bilinear)
        fr[name] = np.clip(dst, 0, 1)

    inside = fr["inside"] >= 0.5
    assert not (inside[0].any() or inside[-1].any() or inside[:, 0].any() or inside[:, -1].any()), "grid clips the country"

    # marz ids: rasterize geoBoundaries polygons in Mercator, then give leftover in-country pixels the nearest marz.
    gj = json.loads(BOUNDARIES.read_text())
    names, shapes_in = [], []
    for k, f in enumerate(sorted(gj["features"], key=lambda f: f["properties"]["shapeName"]), start=1):
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        mp = [[np.column_stack([merc_x(np.array(r)[:, 0]), merc_y(np.array(r)[:, 1])]).tolist() for r in poly] for poly in polys]
        shapes_in.append(({"type": "MultiPolygon", "coordinates": mp}, k))
        names.append({"id": k, "name": f["properties"]["shapeName"], "iso": f["properties"].get("shapeISO")})
    region = rasterize(shapes_in, out_shape=(height, width), transform=dst_tf, fill=0, dtype="uint8")
    nearest = ndimage.distance_transform_edt(region == 0, return_distances=False, return_indices=True)
    filled = region[nearest[0], nearest[1]]
    region = np.where(inside, filled, 0).astype(np.uint8)
    n_reassigned = int(((region > 0) & (rasterize(shapes_in, out_shape=(height, width), transform=dst_tf, fill=0, dtype="uint8") == 0)).sum())

    rows = np.arange(height)
    lat_row = inv_merc(np.zeros(height), y_top - (rows + 0.5) * PX_M)[1]
    area_km2_px = (PX_M * np.cos(np.radians(lat_row))) ** 2 / 1e6      # ground area of one pixel, per row
    area_map = np.broadcast_to(area_km2_px[:, None], (height, width))
    country_km2 = float(area_map[region > 0].sum())
    for r in names:
        r["area_km2"] = round(float(area_map[region == r["id"]].sum()), 1)
    print(f"country area from raster: {country_km2:,.0f} km2 (official {OFFICIAL_AREA_KM2:,.0f}); "
          f"{n_reassigned} px reassigned to nearest marz")
    assert abs(country_km2 / OFFICIAL_AREA_KM2 - 1) < 0.03, "country area off by more than 3%"

    def save(name, arr, scale=255):
        a = np.where(region > 0, np.round(arr * scale), 0).astype(np.uint8) if arr.dtype != np.uint8 else arr
        Image.fromarray(a, mode="L").save(OUT / f"{name}.png", optimize=True)

    save("region", region)
    for name in codes:
        save(name, fr[name])

    # borders: smoothed outlines of each marz (shared edges appear twice; harmless at 1 px) + the country outline
    feats = []

    def ring_to_lonlat(ring):
        pts = rdp(np.array(ring), 0.7 * PX_M)                     # collapse the pixel staircase to a few vertices
        pts = chaikin(pts, iters=3)                                # then round the corners
        pts = np.vstack([rdp(pts[:-1], 0.05 * PX_M), pts[:1]])
        lon, lat = inv_merc(pts[:, 0], pts[:, 1])
        return [[round(float(a), 4), round(float(b), 4)] for a, b in zip(lon, lat)]

    for geom, val in shapes(region, mask=region > 0, transform=dst_tf, connectivity=4):
        for ring in geom["coordinates"][:1]:           # exterior rings only
            if len(ring) < 12:
                continue
            feats.append({"type": "Feature", "properties": {"kind": "marz", "id": int(val)},
                          "geometry": {"type": "LineString", "coordinates": ring_to_lonlat(ring)}})
    for geom, _ in shapes((region > 0).astype(np.uint8), mask=region > 0, transform=dst_tf, connectivity=4):
        for ring in geom["coordinates"][:1]:
            if len(ring) >= 40:
                feats.append({"type": "Feature", "properties": {"kind": "country"},
                              "geometry": {"type": "LineString", "coordinates": ring_to_lonlat(ring)}})
    (OUT / "borders.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}, separators=(",", ":")))

    lon_w, lat_n = inv_merc(x0, y_top)
    lon_e, lat_s = inv_merc(x0 + width * PX_M, y_top - height * PX_M)   # true raster edges, not the requested ones
    grid = {
        "crs": "EPSG:3857", "x0": x0, "y_top": float(y_top), "px_m": PX_M, "width": width, "height": height,
        "bounds_latlon": [[round(float(lat_s), 6), round(float(lon_w), 6)], [round(float(lat_n), 6), round(float(lon_e), 6)]],
        "earth_radius_m": R_EARTH, "regions": names, "country_area_km2": round(country_km2, 1),
        "official_area_km2": OFFICIAL_AREA_KM2,
        "forest_classes": {"forest": FOREST_CODES, "woodland": WOODLAND_CODES, **GROUP_CODES, "water": [WATER_CODE]},
        "class_names": {"31": "Fagus orientalis and other deciduous", "32": "Quercus macranthera and other deciduous",
                        "33": "Quercus macranthera", "34": "Quercus iberica and other deciduous", "35": "Quercus iberica",
                        "36": "Pinus kochiana", "37": "Closed forest plantations", "39": "Open woodland plantations",
                        "41": "Subalpine woodlands", "43": "Mixed woodlands", "44": "Juniper woodlands", "1001": "Water bodies"},
        "cover_area_km2": {k: round(float((fr[k] * area_map)[region > 0].sum()), 1) for k in codes},
        "sources": {"outline_and_cover": "Ecosystem Map of Armenia (BCC Armenia / Institute of Botany NAS RA / IOER, 2026), CC BY 4.0, 10 m, block-averaged to 500 m",
                    "marz_borders": "geoBoundaries gbOpen ARM ADM1 (Wikimedia-derived, 2005), CC BY 2.5; approximate"},
    }
    (OUT / "grid.json").write_text(json.dumps(grid, indent=1))
    print(f"grid {width}x{height}, forest {grid['cover_area_km2']['forest']:,.0f} km2, woodland {grid['cover_area_km2']['woodland']:,.0f} km2, "
          f"{len(feats)} border features -> {OUT}")


if __name__ == "__main__":
    sys.exit(main())
