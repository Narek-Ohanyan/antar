"""Conformal prediction utilities.

Split conformal gives finite-sample marginal coverage under *exchangeability*.  Spatial
data are not exchangeable, so calibration sets must be spatial blocks held out from
training, and under covariate shift (e.g. future climates) the weighted variant
(Tibshirani et al. 2019) with density-ratio weights is required.
"""
from __future__ import annotations

import numpy as np


def conformal_quantile(scores, alpha: float):
    """Finite-sample corrected (1 - alpha) quantile of nonconformity scores."""
    s = np.sort(np.asarray(scores, dtype=float))
    n = s.size
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return float(s[min(k, n) - 1]) if k <= n else float("inf")


def split_conformal_interval(pred_cal, y_cal, pred_new, alpha: float = 0.1):
    q = conformal_quantile(np.abs(np.asarray(y_cal) - np.asarray(pred_cal)), alpha)
    pred_new = np.asarray(pred_new, dtype=float)
    return pred_new - q, pred_new + q


def cqr_interval(lo_cal, hi_cal, y_cal, lo_new, hi_new, alpha: float = 0.1):
    """Conformalised quantile regression (Romano et al. 2019): adaptive-width intervals."""
    y = np.asarray(y_cal, dtype=float)
    scores = np.maximum(np.asarray(lo_cal) - y, y - np.asarray(hi_cal))
    q = conformal_quantile(scores, alpha)
    return np.asarray(lo_new) - q, np.asarray(hi_new) + q


def weighted_conformal_quantile(scores, weights_cal, weight_new: float, alpha: float = 0.1):
    """Weighted quantile of scores with the test point's mass placed at +inf (Tibshirani et al. 2019)."""
    s = np.asarray(scores, dtype=float)
    w = np.asarray(weights_cal, dtype=float)
    order = np.argsort(s)
    s, w = s[order], w[order]
    p = w / (w.sum() + weight_new)
    cum = np.cumsum(p)
    hit = np.where(cum >= 1 - alpha)[0]
    return float(s[hit[0]]) if hit.size else float("inf")
