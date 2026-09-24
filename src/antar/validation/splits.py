"""Cross-validation designs that respect spatial, temporal and hierarchical structure.

* :func:`block_kfold`         spatial blocks (+ optional exclusion buffer around each test block)
* :func:`forward_chaining`    rolling-origin temporal splits
* :func:`leave_group_out`     e.g. leave-one-province-out, leave-one-extreme-year-out
* :func:`empirical_variogram` / :func:`range_from_variogram` choose the block size from the data
* :func:`nn_distance_summary` compare CV distances with prediction-to-training distances (NNDM logic)
"""
from __future__ import annotations

from typing import Iterator

import numpy as np
from scipy.spatial import cKDTree


def spatial_block_ids(x, y, block_size):
    bx = np.floor((np.asarray(x) - np.min(x)) / block_size).astype(int)
    by = np.floor((np.asarray(y) - np.min(y)) / block_size).astype(int)
    return bx * 100000 + by


def block_kfold(x, y, block_size, k: int = 5, buffer: float = 0.0, seed: int = 0) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield (train_idx, test_idx).  Coordinates and ``block_size``/``buffer`` share units.

    Test blocks are assigned to folds at random; training points closer than
    ``buffer`` to any test point are dropped (removes residual spatial leakage).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ids = spatial_block_ids(x, y, block_size)
    ub = np.unique(ids)
    rng = np.random.default_rng(seed)
    rng.shuffle(ub)
    fold_of = {b: i % k for i, b in enumerate(ub)}
    folds = np.array([fold_of[b] for b in ids])
    xy = np.column_stack([x, y])
    for f in range(k):
        te = np.where(folds == f)[0]
        tr = np.where(folds != f)[0]
        if te.size == 0:
            continue
        if buffer > 0 and tr.size:
            d, _ = cKDTree(xy[te]).query(xy[tr], k=1)
            tr = tr[d > buffer]
        yield tr, te


def forward_chaining(years, min_train_years: int, horizon: int = 1):
    """Rolling origin: train on years <= t, test on t+1 .. t+horizon."""
    years = np.asarray(years)
    uy = np.unique(years)
    for i in range(min_train_years, len(uy) - horizon + 1):
        tr = np.where(years <= uy[i - 1])[0]
        te = np.where((years > uy[i - 1]) & (years <= uy[i - 1 + horizon]))[0]
        yield tr, te


def leave_group_out(groups):
    groups = np.asarray(groups)
    for g in np.unique(groups):
        te = np.where(groups == g)[0]
        tr = np.where(groups != g)[0]
        yield g, tr, te


def empirical_variogram(x, y, z, n_bins: int = 20, max_dist: float | None = None):
    """Classical (Matheron) semivariogram: returns bin centres and semivariances."""
    xy = np.column_stack([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    z = np.asarray(z, dtype=float)
    n = len(z)
    iu = np.triu_indices(n, k=1)
    d = np.sqrt(((xy[iu[0]] - xy[iu[1]]) ** 2).sum(axis=1))
    g = 0.5 * (z[iu[0]] - z[iu[1]]) ** 2
    max_dist = max_dist or np.percentile(d, 50)
    edges = np.linspace(0, max_dist, n_bins + 1)
    idx = np.digitize(d, edges) - 1
    centres, gam = [], []
    for b in range(n_bins):
        m = idx == b
        if m.sum() > 5:
            centres.append(0.5 * (edges[b] + edges[b + 1]))
            gam.append(g[m].mean())
    return np.array(centres), np.array(gam)


def range_from_variogram(centres, gamma, frac: float = 0.95):
    """Practical range: first lag at which semivariance reaches ``frac`` of the plateau (median of last third)."""
    if len(gamma) < 4:
        return np.nan
    sill = np.median(gamma[-max(len(gamma) // 3, 1):])
    hit = np.where(gamma >= frac * sill)[0]
    return float(centres[hit[0]]) if hit.size else float(centres[-1])


def nn_distance_summary(train_xy, target_xy, cv_pairs):
    """Median nearest-neighbour distance for (a) target->train and (b) test->train within CV.

    If (b) is much smaller than (a), the CV is optimistic for the prediction task
    (Milà et al. 2022 make the two ECDFs match through NNDM LOO CV).
    """
    train_xy = np.asarray(train_xy, dtype=float)
    d_target, _ = cKDTree(train_xy).query(np.asarray(target_xy, dtype=float), k=1)
    d_cv = []
    for tr, te in cv_pairs:
        d, _ = cKDTree(train_xy[tr]).query(train_xy[te], k=1)
        d_cv.append(d)
    d_cv = np.concatenate(d_cv)
    return float(np.median(d_target)), float(np.median(d_cv))
