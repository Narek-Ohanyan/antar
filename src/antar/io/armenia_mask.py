"""Which CHELSA cells / lat-lon points lie inside Armenia.

The study grid is a rectangle over Armenia that Armenia fills only ~37% of; the rest is Georgia,
Azerbaijan, Turkey, Iran and Nakhchivan. Anything that reports on, averages over, or recommends
action in "Armenia" must be restricted to cells for which :func:`inside_armenia` is True. The mask is
built by ``scripts/build_armenia_mask.py`` from the national Ecosystem Map and committed as a small
.npz so this works without the large raster.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max -- the CHELSA-daily subset
SHAPE = (312, 396)
_PATH = Path(__file__).resolve().parents[3] / "configs" / "armenia_mask_chelsa_312x396.npz"


@lru_cache(maxsize=1)
def load_mask() -> np.ndarray:
    """Boolean (312, 396) array, True where the CHELSA cell centre is inside Armenia."""
    with np.load(_PATH) as f:
        return f["inside"].astype(bool)


def cell_index(lat, lon):
    """CHELSA (row, col) containing each point; may fall outside [0, shape) for points off the grid."""
    lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    row = np.floor((BBOX[3] - lat) / ((BBOX[3] - BBOX[1]) / SHAPE[0])).astype(int)
    col = np.floor((lon - BBOX[0]) / ((BBOX[2] - BBOX[0]) / SHAPE[1])).astype(int)
    return row, col


def inside_armenia(lat, lon) -> np.ndarray:
    """True for points whose CHELSA cell is inside Armenia; False off-grid or in a foreign cell."""
    row, col = cell_index(lat, lon)
    ok = (row >= 0) & (row < SHAPE[0]) & (col >= 0) & (col < SHAPE[1])
    out = np.zeros(np.shape(row), dtype=bool)
    out[ok] = load_mask()[row[ok], col[ok]]
    return out
