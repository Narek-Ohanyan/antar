"""Cohort viability: competing risks and survival to a functional canopy.

    h_tot(tau) = 1 - prod_k (1 - h_k(tau)),        k in {drought, frost, biotic, fire, background}
    V          = prod_{tau=0}^{T-1} (1 - h_tot(tau))  *  P[ H(T) >= H_min ]

Anthropogenic removal (logging, grazing) is an exogenous *scenario input*, not a
climate hazard, and enters as an additional cause k when a management scenario is specified.
"""
from __future__ import annotations

import numpy as np


def combine_competing_hazards(hazards):
    """``hazards``: (K, T) or list of K arrays of per-year cause-specific hazards -> (T,) total hazard."""
    h = np.asarray(hazards, dtype=float)
    return 1.0 - np.prod(1.0 - h, axis=0)


def survival_curve(h_total):
    return np.cumprod(1.0 - np.asarray(h_total, dtype=float))


def viability(h_total, p_height_ok: float = 1.0):
    """Probability that a planted cohort survives to horizon T *and* reaches the functional-canopy height."""
    return float(np.prod(1.0 - np.asarray(h_total, dtype=float)) * p_height_ok)


def cumulative_hazard_from_cloglog(eta):
    """Lambda(T) = sum exp(eta_y); S(T) = exp(-Lambda)."""
    return float(np.sum(np.exp(np.asarray(eta, dtype=float))))
