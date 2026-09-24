"""Extrapolation gating: fall back on physics where the data cannot speak.

    h_final = w(x) * h_stat + (1 - w(x)) * h_mech
    w(x)    = exp( - kappa * max(0, DI(x)/DI_thr - 1) )

w = 1 inside the area of applicability (AOA) and decays smoothly outside it.
"""
from __future__ import annotations

import numpy as np


def aoa_weight(di, di_threshold: float, kappa: float = 2.0):
    di = np.asarray(di, dtype=float)
    return np.exp(-kappa * np.maximum(0.0, di / di_threshold - 1.0))


def blend_hazards(h_stat, h_mech, w):
    """Blend on the complementary-log-log scale so that cumulative hazards, not probabilities, mix."""
    from .models import cloglog, cloglog_inverse

    eta = w * cloglog(h_stat) + (1.0 - w) * cloglog(h_mech)
    return cloglog_inverse(eta)
