"""Quantile delta mapping (Cannon, Sobie & Murdock 2015).

Trend-preserving: the *model-projected change in every quantile* is transferred
onto the observed baseline distribution, which is what we need when the
quantity of interest is a tail (hot-dry extremes), not a mean.
"""
from __future__ import annotations

import numpy as np


def _empirical_cdf_of_sample(x: np.ndarray) -> np.ndarray:
    ranks = np.argsort(np.argsort(x))
    return (ranks + 0.5) / x.size


def quantile_delta_mapping(obs_hist, mod_hist, mod_proj, kind: str = "additive", eps: float = 1e-6):
    """Return bias-adjusted projected series with the same length as ``mod_proj``.

    kind='additive'        temperature-like variables
    kind='multiplicative'  precipitation-like variables (ratio change preserved)
    """
    obs_hist = np.asarray(obs_hist, dtype=float)
    mod_hist = np.asarray(mod_hist, dtype=float)
    mod_proj = np.asarray(mod_proj, dtype=float)
    tau = _empirical_cdf_of_sample(mod_proj)
    q_mod_hist = np.quantile(mod_hist, tau)
    q_obs_hist = np.quantile(obs_hist, tau)
    if kind == "additive":
        return q_obs_hist + (mod_proj - q_mod_hist)
    if kind == "multiplicative":
        ratio = mod_proj / np.maximum(q_mod_hist, eps)
        return q_obs_hist * ratio
    raise ValueError("kind must be 'additive' or 'multiplicative'")
