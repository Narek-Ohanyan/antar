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
