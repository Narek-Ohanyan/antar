"""Proper scoring rules and calibration diagnostics.

Primary metrics for probabilistic outputs are *strictly proper* scores (Brier, log-loss)
reported as skill against an explicit reference, with the Murphy decomposition
BS = REL - RES + UNC.  Cohen's kappa is deliberately not provided: it is prevalence-
dependent and conflates quantity with allocation error (Pontius & Millones 2011).
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score


def brier(p, y):
    p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
    return float(np.mean((p - y) ** 2))


def brier_skill(p, y, p_ref=None):
    """BSS = 1 - BS / BS_ref; the default reference is the climatological event rate."""
    y = np.asarray(y, dtype=float)
    ref = np.full_like(y, y.mean()) if p_ref is None else np.asarray(p_ref, dtype=float)
    return float(1.0 - brier(p, y) / brier(ref, y))


def murphy_decomposition(p, y, n_bins: int = 10):
    """Return dict(brier_binned, reliability, resolution, uncertainty).

    Uses the binned forecast values so that BS_binned = REL - RES + UNC holds exactly.
    """
    p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
    edges = np.quantile(p, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    edges = np.unique(edges)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, len(edges) - 2)
    ybar = y.mean()
    rel = res = 0.0
    p_bin = np.empty_like(p)
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        pk, ok, nk = p[m].mean(), y[m].mean(), m.sum()
        p_bin[m] = pk
        rel += nk * (pk - ok) ** 2
        res += nk * (ok - ybar) ** 2
    n = len(y)
    unc = ybar * (1 - ybar)
    return {"brier_binned": brier(p_bin, y), "reliability": rel / n, "resolution": res / n, "uncertainty": unc}


def log_loss(p, y, eps: float = 1e-12):
    p = np.clip(np.asarray(p, dtype=float), eps, 1 - eps)
    y = np.asarray(y, dtype=float)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


def expected_calibration_error(p, y, n_bins: int = 10):
    p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
    edges = np.quantile(p, np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, len(edges) - 2)
    ece = 0.0
    for b in range(len(edges) - 1):
        m = idx == b
        if m.any():
            ece += m.mean() * abs(p[m].mean() - y[m].mean())
    return float(ece)


def calibration_slope_intercept(p, y, eps: float = 1e-6):
    """Cox calibration: logistic regression of y on logit(p); ideal = slope 1, intercept 0."""
    p = np.clip(np.asarray(p, dtype=float), eps, 1 - eps)
    z = np.log(p / (1 - p)).reshape(-1, 1)
    m = LogisticRegression(C=1e6, max_iter=1000).fit(z, np.asarray(y, dtype=int))
    return float(m.coef_[0, 0]), float(m.intercept_[0])


def pr_auc(p, y):
    return float(average_precision_score(np.asarray(y, dtype=int), np.asarray(p, dtype=float)))


def rmse(pred, obs):
    return float(np.sqrt(np.mean((np.asarray(pred, dtype=float) - np.asarray(obs, dtype=float)) ** 2)))


def mae(pred, obs):
    return float(np.mean(np.abs(np.asarray(pred, dtype=float) - np.asarray(obs, dtype=float))))


def r2_vs_mean(pred, obs, train_mean=None):
    """R^2 relative to a *training* mean (can be negative; that is information, not an error)."""
    obs = np.asarray(obs, dtype=float)
    mu = obs.mean() if train_mean is None else train_mean
    return float(1.0 - np.sum((obs - np.asarray(pred, dtype=float)) ** 2) / np.sum((obs - mu) ** 2))


def interval_coverage(lo, hi, y):
    y = np.asarray(y, dtype=float)
    return float(np.mean((y >= np.asarray(lo)) & (y <= np.asarray(hi))))


def pinball_loss(pred_q, y, tau):
    d = np.asarray(y, dtype=float) - np.asarray(pred_q, dtype=float)
    return float(np.mean(np.maximum(tau * d, (tau - 1) * d)))
