"""Station-based topoclimate downscaling (concept note Sec. 5.1, "Topoclimate" paragraph).

Daily temperature at cell x is

    T(x, d) = T_ref(d) + Gamma_m [z(x) - z_ref] + delta_CAP(x, d)

with the monthly lapse rate Gamma_m regressed from stations (not fixed at
-6.5 K/km) and delta_CAP <= 0 a cold-air-pooling term for T_min. Precipitation
follows a monthly elevation gradient with an exposure factor calibrated on
CHELSA and the stations.

The fitting functions here (``fit_monthly_lapse_rate``,
``fit_precip_elevation_gradient``) take station observations and must only
ever be called inside a training fold; the ``downscale_*`` functions are pure
transforms of already-fitted parameters and carry no fitting logic of their
own.
"""
from __future__ import annotations

import numpy as np


def fit_monthly_lapse_rate(elevation_m, t_obs_c, month, min_stations: int = 3):
    """Ordinary least squares slope of T on elevation, one fit per calendar month.

    ``month`` is a 1-12 integer label aligned with ``elevation_m``/``t_obs_c``
    (multiple stations and/or multiple days per station may repeat the same
    month). Returns a length-12 array (index 0 = January) of Gamma_m, K m-1;
    a month with fewer than ``min_stations`` distinct elevations is NaN.
    """
    z = np.asarray(elevation_m, dtype=float)
    t = np.asarray(t_obs_c, dtype=float)
    m = np.asarray(month, dtype=int)
    gamma = np.full(12, np.nan)
    for mo in range(1, 13):
        sel = m == mo
        if np.sum(sel) < min_stations or np.unique(z[sel]).size < 2:
            continue
        slope, _ = np.polyfit(z[sel], t[sel], 1)
        gamma[mo - 1] = slope
    return gamma


def cold_air_pooling_index(concavity_index, calm_clear_night_frac, k_cap: float = 1.0):
    """delta_CAP <= 0: scales terrain concavity by the share of clear, calm nights.

    The concept note specifies the qualitative behaviour (a cold-air-pooling
    term for T_min that grows with basin-like concavity and with radiative-
    cooling opportunity) but not a closed functional form. This implements the
    simplest form with that behaviour: convex/ridge terrain (concavity <= 0)
    contributes nothing, and the term is linear in concavity and in the
    calm/clear-night fraction, scaled by a fittable coefficient ``k_cap`` >= 0
    (neutral default 1.0, to be calibrated against station T_min in valleys
    vs. slopes). Documented assumption -- see IMPLEMENTATION_LOG.md.
    """
    concavity = np.maximum(np.asarray(concavity_index, dtype=float), 0.0)
    frac = np.clip(np.asarray(calm_clear_night_frac, dtype=float), 0.0, 1.0)
    return -abs(k_cap) * concavity * frac


def downscale_temperature(t_ref_c, z_cell_m, z_ref_m, gamma_k_per_m, delta_cap_k=0.0):
    """T(x, d) = T_ref(d) + Gamma_m [z(x) - z_ref] + delta_CAP(x, d)  (Sec. 5.1)."""
    t_ref = np.asarray(t_ref_c, dtype=float)
    return t_ref + gamma_k_per_m * (np.asarray(z_cell_m, dtype=float) - z_ref_m) + delta_cap_k


def fit_precip_elevation_gradient(elevation_m, p_obs_mm, month, min_stations: int = 3, eps: float = 1e-6):
    """Monthly log-linear precipitation-elevation gradient, one fit per calendar month.

    Fits ln(P) = ln(P0) + grad * z by OLS so the recovered gradient is a
    fractional (multiplicative) change per metre -- precipitation cannot be
    negative, an additive gradient can imply that it is. Returns a length-12
    array (index 0 = January) of grad, m-1; NaN where under-determined.
    """
    z = np.asarray(elevation_m, dtype=float)
    p = np.asarray(p_obs_mm, dtype=float)
    m = np.asarray(month, dtype=int)
    grad = np.full(12, np.nan)
    for mo in range(1, 13):
        sel = m == mo
        if np.sum(sel) < min_stations or np.unique(z[sel]).size < 2:
            continue
        slope, _ = np.polyfit(z[sel], np.log(np.maximum(p[sel], eps)), 1)
        grad[mo - 1] = slope
    return grad


def downscale_precipitation(p_ref_mm, z_cell_m, z_ref_m, gradient_per_m, exposure_factor: float = 1.0):
    """P(x) = P_ref * exp(grad * [z(x) - z_ref]) * exposure_factor, clipped at 0.

    ``exposure_factor`` is a multiplicative, station/CHELSA-calibrated
    correction for wind exposure and rain-shadow effects (neutral default
    1.0 = no correction) -- this is where 1 km CHELSA deltas erase kilometre-
    scale refugia and the exposure factor restores them.
    """
    p_ref = np.asarray(p_ref_mm, dtype=float)
    p = p_ref * np.exp(gradient_per_m * (np.asarray(z_cell_m, dtype=float) - z_ref_m)) * exposure_factor
    return np.maximum(p, 0.0)
