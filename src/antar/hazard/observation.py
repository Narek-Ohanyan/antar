"""Turning raw remote-sensing/dendro observations into the person-period response (Sec. 7.1).

* :func:`dieback_event` -- the standardised-anomaly threshold rule that defines
  a dieback event from an annual-max-kNDVI series. The upstream temporal
  *segmentation* that would produce a cleaned kNDVI series in the first place
  (LandTrendr or CCDC; Kennedy et al. 2010; Zhu & Woodcock 2014) is a
  substantial, data-dependent algorithm in its own right and is not
  implemented here -- this function takes an already-segmented series (or, for
  now, the raw annual maxima) and applies the note's own numeric rule to it.
* :func:`correct_for_misclassification` -- Eq. 7.1, given Se/Sp. Estimating
  Se/Sp themselves needs a stratified sample of interpreted pixel-years
  (real data); only the correction formula is implemented here.
"""
from __future__ import annotations

import numpy as np


def trailing_z(kndvi_annual_max, trailing_window: int = 10):
    """z_t = (x_t - median of the previous ``trailing_window`` years) / (sample standard deviation of those years).

    NaN for the first ``trailing_window`` years, for any year whose window holds a missing value, and for a window whose spread is exactly zero (undetermined, never zero).
    The dieback rule and the continuous vitality response use this one definition.
    """
    x = np.asarray(kndvi_annual_max, dtype=float)
    z = np.full(x.size, np.nan)
    for t in range(trailing_window, x.size):
        window = x[t - trailing_window:t]
        med = np.median(window)
        sd = np.std(window, ddof=1)
        if sd > 0:
            z[t] = (x[t] - med) / sd
    return z


def dieback_event(kndvi_annual_max, no_disturbance, trailing_window: int = 10, onset_z: float = -2.0,
                   recovery_z: float = -1.0, recovery_seasons: int = 2, disturbance_lookback: int = 1, require_run_start: bool = True):
    """Standardised annual-max-kNDVI anomaly rule (Sec. 7.1).

    For year t, z_t = (kndvi_t - median(trailing `trailing_window` years)) /
    std(trailing `trailing_window` years) (the note specifies the thresholds
    numerically but not the spread estimator; trailing-window standard
    deviation is the direct reading of "standardised anomaly relative to its
    trailing 10-year median" and is the choice made here). An event fires at
    t when z_t < onset_z, no_disturbance[t] is True, and z has not risen above
    recovery_z within the following `recovery_seasons` years -- matching
    Table 3's caveat that dieback "must persist >= 2 growing seasons," not a
    single-year dip. Years too close to either end of the record to evaluate
    the trailing window or the recovery window are never flagged (undetermined,
    not asserted).

    ``no_disturbance`` (bool array, same length): True where harvest/fire flags
    are absent -- computed upstream (Table 4: Hansen GFC, MODIS burned area),
    not here.

    Two guards (added 2026-10-07) keep a harvest or fire from being read as dieback one year late. (1) The flag is checked over ``disturbance_lookback`` years before the onset as
    well as in it: a loss recorded late in year t-1 only shows in the Jul-Aug composite of year t. (2) With ``require_run_start`` an onset must begin a decline, that is z_{t-1} was not
    already below ``onset_z``: after a flagged harvest the collapse persists, the trailing window now holds the outlier, z stays below the threshold for years, and without this
    guard every one of those years would be an "event" the flag no longer covers.
    """
    nd = np.asarray(no_disturbance, dtype=bool)
    z = trailing_z(kndvi_annual_max, trailing_window)
    n = z.size
    onset = z < onset_z          # NaN comparisons are False: undetermined years never onset
    event = np.zeros(n, dtype=bool)
    for t in range(n):
        if not onset[t]:
            continue
        future = z[t + 1: t + 1 + recovery_seasons]
        if future.size < recovery_seasons:
            continue             # can't confirm persistence before the record ends: not flagged
        if require_run_start and t > 0 and onset[t - 1]:
            continue             # the decline began earlier: a continuation, not an onset
        if not nd[max(0, t - disturbance_lookback): t + 1].all():
            continue
        event[t] = not np.any(future > recovery_z)
    return event


def correct_for_misclassification(h_obs, se: float, sp: float):
    """Rogan & Gladen (1978), Eq. 7.1: h = (h_obs - (1 - Sp)) / (Se + Sp - 1).

    ``h_obs`` is the satellite-kNDVI-based event rate; ``se``/``sp`` are the
    sensitivity/specificity of that proxy against interpreted ground truth.
    Requires Se + Sp > 1 (the classifier must beat chance) for the correction
    to be identified. The corrected hazard is clipped to [0, 1]: sampling
    noise in Se/Sp can otherwise push it outside that range.
    """
    if se + sp <= 1.0:
        raise ValueError("Se + Sp must exceed 1 (classifier must beat chance) for Eq. 7.1 to be identified")
    h_obs = np.asarray(h_obs, dtype=float)
    h = (h_obs - (1.0 - sp)) / (se + sp - 1.0)
    return np.clip(h, 0.0, 1.0)


def event_rate_interval(n_events: int, n_person_years: int, level: float = 0.95):
    """Exact (Clopper-Pearson) two-sided interval for an event rate k / n, in events per person-year.

    With no events the lower limit is 0 and the upper limit is 1 - (alpha/2)^(1/n) for the two-sided interval; the one-sided 95 % upper bound quoted for "no events seen"
    is :func:`zero_event_upper_bound`. Person-years of the same pixel and of pixels that share a climate cell are not independent, so the interval is narrower than the data justify.
    """
    from scipy.stats import beta
    k, n = int(n_events), int(n_person_years)
    if n <= 0 or k < 0 or k > n:
        raise ValueError("need 0 <= events <= person-years and at least one person-year")
    a = (1.0 - level) / 2.0
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1.0 - a, k + 1, n - k))
    return lo, hi


def zero_event_upper_bound(n_person_years: int, level: float = 0.95) -> float:
    """One-sided upper bound on the rate when no event was seen in ``n_person_years``: 1 - (1 - level)^(1/n), about 3 / n for a 95 % bound ("rule of three")."""
    n = int(n_person_years)
    if n <= 0:
        raise ValueError("need at least one person-year")
    return float(1.0 - (1.0 - level) ** (1.0 / n))
