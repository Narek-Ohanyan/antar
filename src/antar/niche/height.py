"""Attainable (potential) canopy height by quantile regression.

A single-epoch canopy-height map mixes climate, soil *and* stand age/disturbance
history.  Modelling the conditional mean confuses these; modelling an upper
quantile (tau ~ 0.9) targets the height a site *can* reach - the analogue of
site index (Cade & Noon 2003; Skovsgaard & Vanclay 2008).
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import GradientBoostingRegressor


def fit_attainable_height(X, height, tau: float = 0.9, **kw):
    params = dict(loss="quantile", alpha=tau, n_estimators=300, learning_rate=0.05, max_depth=3, subsample=0.8, random_state=0)
    params.update(kw)
    return GradientBoostingRegressor(**params).fit(np.asarray(X, dtype=float), np.asarray(height, dtype=float))


def heteroscedastic_measurement_var(h_obs, a: float = 0.5, b: float = 0.05):
    """Simple variance model sigma_e^2 = (a + b H)^2 for map-product error; replace with validation-derived fit."""
    return (a + b * np.asarray(h_obs, dtype=float)) ** 2


def water_limited_height_fallback(h_star_reference, aet_mm, pet_mm):
    """Outside the area of applicability, H* falls back on a water-limited scaling by
    the transpiration ratio AET/PET from TOPOHYDRO (Sec. 8.1), rather than trusting
    the quantile-regression extrapolation.

    The note names the mechanism but not its exact form; implemented as
    ``H*_fallback = h_star_reference * clip(AET/PET, 0, 1)`` -- capped at the
    reference value, since a wetter-than-potential ratio cannot exceed the
    unconstrained attainable height. A documented simplification, not a quoted
    formula. Use with :func:`antar.validation.aoa.AreaOfApplicability.inside`
    to decide which cells need it.
    """
    aet = np.asarray(aet_mm, dtype=float)
    pet = np.asarray(pet_mm, dtype=float)
    ratio = np.clip(np.divide(aet, pet, out=np.zeros_like(aet, dtype=float), where=pet > 0), 0.0, 1.0)
    return np.asarray(h_star_reference, dtype=float) * ratio
