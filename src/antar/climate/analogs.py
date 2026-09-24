"""Climate analogues in interannual-variability units.

For a target cell/time we measure the Mahalanobis distance to every reference
cell *using the reference cell's own interannual covariance* - in the spirit of
the sigma-dissimilarity of Mahony et al. (2017).  Small distance = plausible
analogue (a seed-source candidate); large minimum distance = novel climate.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def mahalanobis_to_reference(x_target, ref_means, ref_cov):
    """Distance of ``x_target`` (k,) to each reference mean (m, k), using ref_cov (k, k) or (m, k, k)."""
    x = np.asarray(x_target, dtype=float)
    mu = np.asarray(ref_means, dtype=float)
    cov = np.asarray(ref_cov, dtype=float)
    diff = x[None, :] - mu
    if cov.ndim == 2:
        inv = np.linalg.inv(cov)
        d2 = np.einsum("mi,ij,mj->m", diff, inv, diff)
    else:
        inv = np.linalg.inv(cov)
        d2 = np.einsum("mi,mij,mj->m", diff, inv, diff)
    return np.sqrt(np.maximum(d2, 0.0))


def novelty_percentile(distance, k: int):
    """Chi-square percentile of a squared Mahalanobis distance with k variables."""
    return stats.chi2.cdf(np.asarray(distance, dtype=float) ** 2, df=k)


def best_analogue(x_target, ref_means, ref_cov):
    d = mahalanobis_to_reference(x_target, ref_means, ref_cov)
    j = int(np.argmin(d))
    return j, float(d[j])
