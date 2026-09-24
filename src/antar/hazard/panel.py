"""Panel utilities: within-between (Mundlak) decomposition and person-period expansion.

Why the decomposition matters
-----------------------------
With x_it = xbar_i + (x_it - xbar_i), the coefficient on the *between* part xbar_i
captures long-run spatial differences (adaptation, site quality, confounded by
everything that co-varies in space); the coefficient on the *within* part
(x_it - xbar_i) is identified from year-to-year anomalies at the same place and
is the quantity that transfers to a warming climate (Mundlak 1978).  A model that
pools both (as v1 did) cannot tell them apart.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def mundlak_decompose(df: pd.DataFrame, group: str, cols: list[str], suffix_bar: str = "_bar", suffix_dev: str = "_dev"):
    out = df.copy()
    for c in cols:
        bar = out.groupby(group)[c].transform("mean")
        out[c + suffix_bar] = bar
        out[c + suffix_dev] = out[c] - bar
    return out


def rare_event_weights(event, pi0: float):
    """1 for kept events, 1/pi0 for kept non-events sampled at a known rate pi0 (Sec. 7.2).

    All events are kept; non-events are sampled within blocks at a known
    fraction pi0 and weighted by 1/pi0, a consistent estimator of the
    full-data model under choice-based sampling (King & Zeng 2001).
    """
    if not (0.0 < pi0 <= 1.0):
        raise ValueError("pi0 (the non-event sampling fraction) must be in (0, 1]")
    event = np.asarray(event, dtype=bool)
    return np.where(event, 1.0, 1.0 / pi0)


def no_adapt_bracket(x_proj, xbar_hist):
    """eta^no_adapt (Eq. 7.5): x~ = x_proj - xbar_hist, xbar = xbar_hist.

    The whole projected climatic shift is treated as a year-to-year departure
    at a fixed historical place mean, so beta_W (identified from weather, not
    confounded by adaptation) applies to the entire shift.
    """
    x_proj = np.asarray(x_proj, dtype=float)
    xbar_hist = np.asarray(xbar_hist, dtype=float)
    return x_proj - xbar_hist, xbar_hist


def full_acclimation_bracket(x_proj, xbar_proj):
    """eta^full_accl (Eq. 7.5): x~ = x_proj - xbar_proj, xbar = xbar_proj.

    Trees are assumed to fully track their own projected place mean, so only
    beta_B (typically weaker, confounded with adaptation/site quality) applies
    to the shift of the mean itself.
    """
    x_proj = np.asarray(x_proj, dtype=float)
    xbar_proj = np.asarray(xbar_proj, dtype=float)
    return x_proj - xbar_proj, xbar_proj


def space_for_time_bracket(x_proj, xbar_hist, xbar_proj) -> dict:
    """Both Eq. 7.5 brackets at once: "every projection is reported as this bracket,
    and the width between the bounds is itself a result" (Sec. 7.4).

    Returns ``{"no_adapt": (x_tilde, xbar), "full_accl": (x_tilde, xbar)}``.
    """
    return {
        "no_adapt": no_adapt_bracket(x_proj, xbar_hist),
        "full_accl": full_acclimation_bracket(x_proj, xbar_proj),
    }


def person_period(ids, first_year, last_year, event_year=None) -> pd.DataFrame:
    """Expand units into unit-year rows with a discrete-time event indicator.

    ``event_year`` may be NaN/None for censored units.  Rows after the event are
    dropped (first-event analysis).  Columns: id, year, event.
    """
    ids = np.asarray(ids)
    first_year = np.broadcast_to(np.asarray(first_year), ids.shape)
    last_year = np.broadcast_to(np.asarray(last_year), ids.shape)
    ev = np.full(ids.shape, np.nan) if event_year is None else np.asarray(event_year, dtype=float)
    rows = []
    for i, f, l, e in zip(ids, first_year, last_year, ev):
        stop = l if np.isnan(e) else min(l, int(e))
        for y in range(int(f), int(stop) + 1):
            rows.append((i, y, int((not np.isnan(e)) and y == int(e))))
    return pd.DataFrame(rows, columns=["id", "year", "event"])
