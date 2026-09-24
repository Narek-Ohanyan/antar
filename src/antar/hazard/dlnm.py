"""Distributed-lag non-linear cross-basis (Gasparrini et al. 2010), minimal version.

For exposure x observed at lags 0..L the linear predictor contribution is

    eta_lag(i) = sum_{l=0}^{L} f(x_{i,l}, l),    f(x, l) = sum_j sum_k theta_jk b_j(x) c_k(l)

The design columns returned here are  sum_l b_j(x_{i,l}) c_k(l)  for each (j, k).
Drought legacies (Anderegg et al. 2015) are exactly this kind of lagged, non-linear effect.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import BSpline


def _bspline_basis(x, knots_inner, degree, lo, hi):
    x = np.clip(np.asarray(x, dtype=float), lo, hi)
    t = np.r_[[lo] * (degree + 1), np.asarray(knots_inner, dtype=float), [hi] * (degree + 1)]
    n = len(t) - degree - 1
    return BSpline.design_matrix(x.ravel(), t, degree).toarray().reshape(*x.shape, n)


def spline_basis(x, knots_inner, degree: int = 3, x_range=None):
    """A plain (non-cross) B-spline basis for a single smooth term, e.g. f(a_i,t) in Eq. 7.3
    (time since disturbance/planting). Reuses the same B-spline machinery as
    :func:`crossbasis`'s exposure axis; f(a) = spline_basis(a, knots) @ coefficients,
    fitted as ordinary columns of the hazard model's design matrix.
    """
    x = np.asarray(x, dtype=float)
    lo, hi = (float(np.min(x)), float(np.max(x))) if x_range is None else x_range
    return _bspline_basis(x, knots_inner, degree, lo, hi)


def crossbasis(x_lagged, exposure_knots, lag_knots, degree_x: int = 2, degree_lag: int = 2, x_range=None):
    """Cross-basis matrix.

    x_lagged : (n, L+1) array; column l holds the exposure l years before the outcome year.
    Returns  : (n, J*K) design matrix.
    """
    x_lagged = np.asarray(x_lagged, dtype=float)
    n, nl = x_lagged.shape
    lo, hi = (np.min(x_lagged), np.max(x_lagged)) if x_range is None else x_range
    bx = _bspline_basis(x_lagged, exposure_knots, degree_x, lo, hi)          # (n, L+1, J)
    lags = np.arange(nl, dtype=float)
    bl = _bspline_basis(lags, lag_knots, degree_lag, 0.0, float(nl - 1))     # (L+1, K)
    cb = np.einsum("nlj,lk->njk", bx, bl)
    return cb.reshape(n, -1)


def lag_matrix(series, max_lag: int):
    """Build (n - max_lag, max_lag + 1) lagged exposure matrix from a 1-D annual series."""
    s = np.asarray(series, dtype=float)
    n = s.size
    cols = [s[max_lag - l: n - l] for l in range(max_lag + 1)]
    return np.column_stack(cols)
