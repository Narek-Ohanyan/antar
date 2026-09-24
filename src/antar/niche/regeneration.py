"""Establishment (regeneration) as a discrete-time hazard over the first seasons after planting.

The regeneration niche is narrower than the adult niche (Jackson et al. 2009;
Davis et al. 2019; Dobrowski et al. 2015).  Modelling it with the same hazard
family as later life stages keeps the whole cohort model coherent.
"""
from __future__ import annotations

import numpy as np

from ..hazard.models import cloglog_inverse


def establishment_hazard(beta0, beta, x_by_season, age_effect=None):
    """h_tau = cloglog^-1(beta0 + x_tau . beta + age_effect_tau); x_by_season: (n_seasons, p)."""
    x = np.asarray(x_by_season, dtype=float)
    eta = beta0 + x @ np.asarray(beta, dtype=float)
    if age_effect is not None:
        eta = eta + np.asarray(age_effect, dtype=float)
    return cloglog_inverse(eta)


def establishment_survival(hazards):
    """S_est = prod (1 - h_tau) over the establishment window."""
    return float(np.prod(1.0 - np.asarray(hazards, dtype=float)))
