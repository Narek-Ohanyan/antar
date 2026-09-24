"""Adult-niche diagnostics for presence-only occurrence data.

Presence-background models give *relative* suitability, not the probability of
occurrence (prevalence is unknown), and AUC from pseudo-absences is a poor
yardstick.  The continuous Boyce index (Hirzel et al. 2006) compares the
predicted-to-expected frequency ratio along the suitability gradient.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def boyce_index(pred_at_presence, pred_all, n_windows: int = 100, window_width: float = 0.1):
    pp = np.asarray(pred_at_presence, dtype=float)
    pa = np.asarray(pred_all, dtype=float)
    lo, hi = np.min(pa), np.max(pa)
    rng = hi - lo
    if rng <= 0:
        return np.nan
    w = window_width * rng
    starts = np.linspace(lo, hi - w, n_windows)
    pe, mids = [], []
    for s in starts:
        e = s + w
        p = np.mean((pp >= s) & (pp <= e))
        a = np.mean((pa >= s) & (pa <= e))
        if a > 0:
            pe.append(p / a)
            mids.append(0.5 * (s + e))
    if len(pe) < 3:
        return np.nan
    return float(spearmanr(mids, pe).statistic)
