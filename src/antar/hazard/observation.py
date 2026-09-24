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


def dieback_event(kndvi_annual_max, no_disturbance, trailing_window: int = 10, onset_z: float = -2.0,
                   recovery_z: float = -1.0, recovery_seasons: int = 2):
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
    """
    x = np.asarray(kndvi_annual_max, dtype=float)
    nd = np.asarray(no_disturbance, dtype=bool)
    n = x.size
    z = np.full(n, np.nan)
    for t in range(trailing_window, n):
        window = x[t - trailing_window:t]
        med = np.median(window)
        sd = np.std(window, ddof=1)
        if sd > 0:
            z[t] = (x[t] - med) / sd
    onset = z < onset_z          # NaN comparisons are False: undetermined years never onset
    event = np.zeros(n, dtype=bool)
    for t in range(n):
        if not onset[t]:
            continue
        future = z[t + 1: t + 1 + recovery_seasons]
        if future.size < recovery_seasons:
            continue             # can't confirm persistence before the record ends: not flagged
        event[t] = bool(nd[t]) and not np.any(future > recovery_z)
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
