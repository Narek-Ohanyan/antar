"""Neighbourhood terrain operations on a gridded DEM -- flagged as missing since
``antar.io.gee_export.export_terrain`` was first built ("TWI, TPI and the cold-
air-pooling concavity index are not computed here: they are neighbourhood/flow-
routing operations GEE's per-pixel ``ee.Terrain`` functions do not provide, and
are left to local post-processing on the downloaded DEM"). DEM, slope and aspect
themselves come directly from that export (via ``ee.Terrain.products``); this
module is only the neighbourhood piece.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter


def concavity_index(elevation_m, radius_px: int = 1):
    """Topographic Position Index, sign-flipped so basins/valleys (cold-air-
    pooling sites) are positive and ridges/convex terrain are negative -- the
    convention ``antar.climate.downscale.cold_air_pooling_index`` expects
    (it clamps negative values to zero itself).

    concavity(x) = mean(elevation in a (2*radius_px+1)-square neighbourhood,
    excluding the centre pixel) - elevation(x), metres. Positive where a cell
    sits below its surroundings (a basin), negative on a local high point
    (a ridge) -- a standard, but not the only possible, TPI-style definition;
    documented as a choice, not derived from the concept note (which specifies
    the qualitative behaviour of ``cold_air_pooling_index`` but not how the
    concavity term itself should be computed from a DEM).

    Edge pixels (within ``radius_px`` of the array boundary) use whatever
    neighbourhood actually exists inside the array rather than wrapping or
    padding with a fabricated elevation value -- implemented via an explicit
    valid-pixel count per window (``uniform_filter`` on a ones-array), not by
    letting scipy silently pad with zeros or mirror values into the mean.
    """
    z = np.asarray(elevation_m, dtype=float)
    if z.ndim != 2:
        raise ValueError(f"elevation_m must be a 2D array, got shape {z.shape}")

    size = 2 * radius_px + 1
    window_sum = uniform_filter(z, size=size, mode="constant", cval=0.0) * size * size
    valid_count = uniform_filter(np.ones_like(z), size=size, mode="constant", cval=0.0) * size * size
    # Exclude the centre pixel itself from its own neighbourhood mean.
    neighbourhood_mean = (window_sum - z) / (valid_count - 1)

    return neighbourhood_mean - z
