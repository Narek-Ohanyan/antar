"""Agro-climatic and drought indices used as hazard covariates."""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import gamma as _gamma


def growing_degree_days(t_mean_c, base_c: float = 5.0):
    return np.maximum(np.asarray(t_mean_c, dtype=float) - base_c, 0.0)


def late_frost_days(t_min_c, gdd_cumulative, gdd_budburst: float, t_crit_c: float = -2.0):
    """Number of days with frost after (thermal-time) budburst."""
    tmin = np.asarray(t_min_c, dtype=float)
    return int(np.sum((gdd_cumulative >= gdd_budburst) & (tmin < t_crit_c)))


def growing_season_length(t_mean_c, t0_c: float = 0.9) -> int:
    """LGS: number of days in the year with T_mean >= T0 (Koerner & Paulsen 2004; concept note Eq. 8.3).

    T0 = 0.9 degC is the pan-biome tree-form growing-season threshold; the
    framework's own operational definition of the growing season is exactly
    this day count, not a contiguous-run rule.
    """
    return int(np.sum(np.asarray(t_mean_c, dtype=float) >= t0_c))


def growing_season_mean_temperature(t_mean_c, t0_c: float = 0.9) -> float:
    """GST: mean T_mean over the days with T_mean >= T0 (concept note Eq. 8.3). NaN if LGS = 0."""
    t = np.asarray(t_mean_c, dtype=float)
    season = t[t >= t0_c]
    return float(np.mean(season)) if season.size else float("nan")


def _loglogistic_params_pwm(x: np.ndarray):
    """3-parameter log-logistic via probability-weighted moments (Vicente-Serrano et al. 2010)."""
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    i = np.arange(1, n + 1)
    f = (i - 0.35) / n
    w0 = np.mean((1 - f) ** 0 * x)
    w1 = np.mean((1 - f) ** 1 * x)
    w2 = np.mean((1 - f) ** 2 * x)
    beta = (2 * w1 - w0) / (6 * w1 - w0 - 6 * w2)
    g1 = _gamma(1 + 1 / beta)
    g2 = _gamma(1 - 1 / beta)
    alpha = (w0 - 2 * w1) * beta / (g1 * g2)
    gam = w0 - alpha * g1 * g2
    return alpha, beta, gam


def spei_from_series(d_series, scale: int = 3, period: int = 12, ref_slice: slice | None = None):
    """SPEI from a monthly climatic water balance series ``D = P - PET``.

    ``ref_slice`` (indices into the *accumulated* series) fixes the reference
    period, so that future series are standardised against the *historical*
    distribution.  Returns an array aligned with the input (NaN for the first
    ``scale-1`` months).
    """
    d = np.asarray(d_series, dtype=float)
    acc = np.full_like(d, np.nan)
    c = np.cumsum(np.insert(d, 0, 0.0))
    acc[scale - 1:] = c[scale:] - c[:-scale]
    out = np.full_like(d, np.nan)
    idx = np.arange(d.size)
    for m in range(period):
        sel = (idx % period == m) & ~np.isnan(acc)
        ref = sel.copy()
        if ref_slice is not None:
            mask = np.zeros(d.size, dtype=bool)
            mask[ref_slice] = True
            ref &= mask
        a, b, g = _loglogistic_params_pwm(acc[ref])
        with np.errstate(invalid="ignore", divide="ignore"):
            cdf = 1.0 / (1.0 + (a / (acc[sel] - g)) ** b)
        cdf = np.clip(cdf, 1e-6, 1 - 1e-6)
        out[sel] = stats.norm.ppf(cdf)
    return out
