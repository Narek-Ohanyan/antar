"""Ensemble uncertainty attribution.

* :func:`anova_fractions`  main-effect variance fractions for a full-factorial ensemble
  (Hawkins & Sutton 2009; Lehner et al. 2020) - scenario, GCM, downscaling, PET, statistical model, ...
* :func:`sobol_indices`    first-order and total Sobol' indices with the Saltelli (2010)
  design (Saltelli/Jansen estimators) for the continuous parameter inputs.
"""
from __future__ import annotations

from typing import Callable

import numpy as np


def anova_fractions(Y, factor_names):
    """``Y`` has one axis per factor (full factorial, no replication).  Returns dict name -> variance fraction,
    plus 'interaction' (everything not explained by main effects)."""
    Y = np.asarray(Y, dtype=float)
    total = Y.var()
    out = {}
    explained = 0.0
    for ax, name in enumerate(factor_names):
        other = tuple(i for i in range(Y.ndim) if i != ax)
        v = Y.mean(axis=other).var()
        out[name] = float(v / total) if total > 0 else 0.0
        explained += out[name]
    out["interaction"] = float(max(0.0, 1.0 - explained))
    return out


def sobol_indices(f: Callable[[np.ndarray], np.ndarray], bounds, n: int = 20000, seed: int = 0):
    """First-order (S1) and total (ST) indices; ``f`` maps (m, d) -> (m,)."""
    bounds = np.asarray(bounds, dtype=float)
    d = bounds.shape[0]
    rng = np.random.default_rng(seed)
    A = rng.random((n, d))
    B = rng.random((n, d))
    scale = lambda U: bounds[:, 0] + U * (bounds[:, 1] - bounds[:, 0])
    fA, fB = f(scale(A)), f(scale(B))
    var = np.var(np.concatenate([fA, fB]))
    s1, st = np.empty(d), np.empty(d)
    for i in range(d):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        fABi = f(scale(ABi))
        s1[i] = np.mean(fB * (fABi - fA)) / var            # Saltelli et al. 2010
        st[i] = 0.5 * np.mean((fA - fABi) ** 2) / var      # Jansen estimator
    return s1, st
