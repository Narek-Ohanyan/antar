"""Degree-day snow module (mm w.e.).  Calibrated against MODIS/Sentinel-2 snow cover."""
from __future__ import annotations

import numpy as np


def partition_precipitation(p_mm, t_mean_c, t_snow_c: float = 1.0):
    """Split precipitation into liquid and solid parts by a temperature threshold."""
    p = np.asarray(p_mm, dtype=float)
    snow = np.where(np.asarray(t_mean_c, dtype=float) < t_snow_c, p, 0.0)
    return p - snow, snow


def simulate_snow(p_mm, t_mean_c, ddf: float = 4.0, t_melt_c: float = 0.0, t_snow_c: float = 1.0, swe0: float = 0.0):
    """Daily degree-day snowpack.

    Returns dict(swe, melt, liquid_input) where ``liquid_input`` = rain + melt,
    i.e. the water that reaches the soil-water balance.
    """
    p = np.asarray(p_mm, dtype=float)
    t = np.asarray(t_mean_c, dtype=float)
    rain, snow = partition_precipitation(p, t, t_snow_c)
    swe = np.zeros_like(p)
    melt = np.zeros_like(p)
    s = np.full(p.shape[1:], swe0, dtype=float) if p.ndim > 1 else float(swe0)
    for d in range(p.shape[0]):
        s = s + snow[d]
        m = np.minimum(s, ddf * np.maximum(t[d] - t_melt_c, 0.0))
        s = s - m
        swe[d] = s
        melt[d] = m
    return {"swe": swe, "melt": melt, "liquid_input": rain + melt}
