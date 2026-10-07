"""Ground area of the grid cell that a model node stands for.

A node of the dense grid is one CHELSA pixel out of every ``stride`` (30 arcseconds, about 0.93 km, so about 6.5 km for stride 7). The node's results (viability,
eligibility) are taken to represent the whole cell around it, so a planting unit is that cell and not one hectare. The area is the exact area of a longitude-
latitude rectangle on a sphere: R^2 * dlon * (sin(lat + dlat/2) - sin(lat - dlat/2)).
"""
from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0088        # IUGG mean radius


def cell_area_ha(lat_deg, dlat_deg: float, dlon_deg: float):
    """Area in hectares of the cell centred on ``lat_deg`` that spans ``dlat_deg`` of latitude and ``dlon_deg`` of longitude (broadcasts over lat_deg)."""
    lat = np.asarray(lat_deg, dtype=float)
    south, north = np.radians(lat - dlat_deg / 2.0), np.radians(lat + dlat_deg / 2.0)
    km2 = EARTH_RADIUS_KM ** 2 * np.radians(dlon_deg) * (np.sin(north) - np.sin(south))
    return km2 * 100.0             # 1 km2 = 100 ha


def cell_land_shares(lat, lon, dlat_deg: float, dlon_deg: float, region, layers: dict, grid: dict) -> dict:
    """Mean of each fractional land-cover layer over the in-country raster pixels that a cell covers.

    ``region`` is the raster of region ids (> 0 inside the country), ``layers`` maps a name to a raster of fractions in [0, 1] on the same Web-Mercator grid
    (``grid``: x0, y_top, px_m, earth_radius_m; column = floor((x - x0) / px_m), row = floor((y_top - y) / px_m), as on the map page). Each cell is the
    rectangle of ``dlat_deg`` x ``dlon_deg`` around its node; a pixel counts whole if its index lies within the rectangle's pixel range. Returns one array per layer
    (NaN where the cell covers no in-country pixel) and ``coverage``: the share of the pixels in the rectangle that are inside the country.
    """
    lat, lon = np.atleast_1d(np.asarray(lat, dtype=float)), np.atleast_1d(np.asarray(lon, dtype=float))
    region = np.asarray(region)
    H, W = region.shape
    R, x0, y_top, px = grid["earth_radius_m"], grid["x0"], grid["y_top"], grid["px_m"]
    mx = lambda lo: R * np.radians(lo)
    my = lambda la: R * np.log(np.tan(np.pi / 4 + np.radians(la) / 2))
    out = {k: np.full(lat.shape, np.nan) for k in layers}
    coverage = np.zeros(lat.shape)
    for i in range(lat.size):
        c0, c1 = int(np.floor((mx(lon[i] - dlon_deg / 2) - x0) / px)), int(np.floor((mx(lon[i] + dlon_deg / 2) - x0) / px))
        r0, r1 = int(np.floor((y_top - my(lat[i] + dlat_deg / 2)) / px)), int(np.floor((y_top - my(lat[i] - dlat_deg / 2)) / px))
        r0, r1, c0, c1 = max(r0, 0), min(r1, H - 1), max(c0, 0), min(c1, W - 1)
        if r1 < r0 or c1 < c0:
            continue
        inside = region[r0:r1 + 1, c0:c1 + 1] > 0
        coverage[i] = inside.mean()
        if inside.any():
            for k, layer in layers.items():
                out[k][i] = float(np.asarray(layer)[r0:r1 + 1, c0:c1 + 1][inside].mean())
    out["coverage"] = coverage
    return out
