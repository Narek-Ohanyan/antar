"""Design-based accuracy and area estimation from a stratified probability sample.

Implements the estimators of Olofsson et al. (2014) for an error matrix of sample counts
n_hj (stratum h = map class, column j = reference class) with stratum area weights W_h:

    p_hat_j       = sum_h W_h * n_hj / n_h
    Var(p_hat_j)  = sum_h W_h^2 * [ (n_hj/n_h)(1 - n_hj/n_h) / (n_h - 1) ]
    A_hat_j       = A_tot * p_hat_j        (area estimate; SE scales identically)
    overall acc.  = sum_h W_h * n_hh / n_h

These estimates are what a map's *stated* accuracy and area should be based on; blocked CV
is for model selection, this is for the map product (cf. Wadoux et al. 2021; Stehman & Foody 2019).
"""
from __future__ import annotations

import numpy as np


def stratified_estimates(n_hj, stratum_area, total_area=None):
    n_hj = np.asarray(n_hj, dtype=float)
    a_h = np.asarray(stratum_area, dtype=float)
    a_tot = a_h.sum() if total_area is None else total_area
    w_h = a_h / a_h.sum()
    n_h = n_hj.sum(axis=1)
    frac = n_hj / n_h[:, None]
    p_j = (w_h[:, None] * frac).sum(axis=0)
    var_j = ((w_h[:, None] ** 2) * frac * (1 - frac) / (n_h[:, None] - 1)).sum(axis=0)
    overall = float((w_h * np.diag(frac)).sum())
    var_overall = float(((w_h**2) * np.diag(frac) * (1 - np.diag(frac)) / (n_h - 1)).sum())
    return {
        "p_hat": p_j,
        "se_p": np.sqrt(var_j),
        "area_hat": a_tot * p_j,
        "area_ci95": 1.96 * a_tot * np.sqrt(var_j),
        "overall_accuracy": overall,
        "se_overall": float(np.sqrt(var_overall)),
    }
