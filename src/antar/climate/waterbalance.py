"""Daily soil-water bucket with relative extractable water (REW) stress.

Following Granier et al. (1999): transpiration is demand-limited while
REW >= REW_crit (0.4) and declines linearly below it.  Storage W (mm) is the
water held between field capacity and the *species-specific extraction limit*;
``w_max`` = plant-available capacity x rooting depth (RZ-PAW).

Outputs feed XYLEM (soil water potential) and the climatic-water-deficit
indices used as hazard covariates.
"""
from __future__ import annotations

import numpy as np


def simulate_bucket(pet_mm, water_input_mm, w_max_mm, rew_crit: float = 0.4, w0_frac: float = 1.0):
    """Daily bucket.  Inputs may be shape (T,) or (T, N); ``w_max_mm`` scalar or (N,).

    Returns dict with ``w`` (end-of-day storage), ``rew``, ``aet``, ``deficit``
    (PET - AET), ``drainage`` and the closure error ``balance_error`` (should be ~0).
    """
    pet = np.asarray(pet_mm, dtype=float)
    inp = np.asarray(water_input_mm, dtype=float)
    w_max = np.asarray(w_max_mm, dtype=float)
    shape = pet.shape
    w = np.broadcast_to(w_max * w0_frac, shape[1:] if pet.ndim > 1 else ()).astype(float).copy()
    out_w = np.zeros(shape)
    out_aet = np.zeros(shape)
    out_dr = np.zeros(shape)
    out_rew = np.zeros(shape)
    w_start = w.copy()
    for d in range(shape[0]):
        w_star = w + inp[d]
        drain = np.maximum(w_star - w_max, 0.0)
        w_star = np.minimum(w_star, w_max)
        rew = np.where(w_max > 0, w_star / np.where(w_max > 0, w_max, 1.0), 0.0)
        supply_factor = np.minimum(1.0, rew / rew_crit)
        aet = np.minimum(w_star, pet[d] * supply_factor)
        w = w_star - aet
        out_w[d], out_aet[d], out_dr[d], out_rew[d] = w, aet, drain, rew
    balance_error = (w_start + inp.sum(axis=0) - out_aet.sum(axis=0) - out_dr.sum(axis=0) - w)
    return {
        "w": out_w,
        "rew": out_rew,
        "aet": out_aet,
        "deficit": pet - out_aet,
        "drainage": out_dr,
        "balance_error": balance_error,
    }


def climatic_water_deficit(pet_mm, aet_mm, axis: int = 0):
    """CWD = sum(PET - AET) (Stephenson 1998; Lutz et al. 2010)."""
    return np.sum(np.asarray(pet_mm, dtype=float) - np.asarray(aet_mm, dtype=float), axis=axis)


def water_stress_integral(rew, rew_crit: float = 0.4, axis: int = 0):
    """WSI = sum max(0, (REW_crit - REW)/REW_crit): a duration x intensity index (Granier 1999)."""
    rew = np.asarray(rew, dtype=float)
    return np.sum(np.maximum(0.0, (rew_crit - rew) / rew_crit), axis=axis)


def theta_from_storage(w_mm, w_max_mm, theta_fc, theta_lim):
    """Volumetric water content from storage between the extraction limit and field capacity."""
    frac = np.asarray(w_mm, dtype=float) / np.asarray(w_max_mm, dtype=float)
    return theta_lim + (theta_fc - theta_lim) * frac


def psi_clapp_hornberger(theta, theta_sat, psi_sat_mpa, b):
    """Soil water potential (MPa, negative): psi_sat * (theta/theta_sat)^(-b)."""
    theta = np.maximum(np.asarray(theta, dtype=float), 1e-6)
    return psi_sat_mpa * (theta / theta_sat) ** (-b)
