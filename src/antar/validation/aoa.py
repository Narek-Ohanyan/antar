"""Area of applicability (Meyer & Pebesma 2021), NumPy implementation.

DI(x) = min_j ||x - x_j|| / dbar, computed in a standardised, importance-weighted feature
space; dbar is the mean pairwise distance among training points.  The threshold is the
upper whisker (Q3 + 1.5 IQR) of the training points' DI values, each computed against
training points *outside its own CV fold*, so it reflects the distances the model was
actually validated over.  Cells with DI above the threshold are outside the AOA:
predictions there are not supported by the training data and are gated to the
mechanistic fallback (see hazard.gating).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


class AreaOfApplicability:
    def __init__(self, weights=None):
        self.weights = weights

    def fit(self, X_train, folds=None, rng_seed: int = 0, max_pairs: int = 200000):
        X = np.asarray(X_train, dtype=float)
        self.mu_ = X.mean(0)
        self.sd_ = X.std(0) + 1e-12
        w = np.ones(X.shape[1]) if self.weights is None else np.asarray(self.weights, dtype=float)
        self.w_ = w
        Z = self._transform(X)
        self.Z_ = Z
        rng = np.random.default_rng(rng_seed)
        n = len(Z)
        i = rng.integers(0, n, max_pairs)
        j = rng.integers(0, n, max_pairs)
        keep = i != j
        self.dbar_ = float(np.mean(np.linalg.norm(Z[i[keep]] - Z[j[keep]], axis=1)))
        folds = np.arange(n) if folds is None else np.asarray(folds)
        di = np.empty(n)
        for f in np.unique(folds):
            te = folds == f
            tr = ~te
            d, _ = cKDTree(Z[tr]).query(Z[te], k=1)
            di[te] = d / self.dbar_
        q1, q3 = np.percentile(di, [25, 75])
        self.threshold_ = float(q3 + 1.5 * (q3 - q1))
        self.tree_ = cKDTree(Z)
        return self

    def _transform(self, X):
        return (np.asarray(X, dtype=float) - self.mu_) / self.sd_ * self.w_

    def dissimilarity_index(self, X_new):
        d, _ = self.tree_.query(self._transform(X_new), k=1)
        return d / self.dbar_

    def inside(self, X_new):
        return self.dissimilarity_index(X_new) <= self.threshold_
