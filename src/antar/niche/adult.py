"""Adult-niche model and diagnostics for presence-background occurrence data
(concept note Sec 7, Module D: "a penalised presence-background model on
physically meaningful predictors ... with the background drawn from a
target-group sample to absorb collection bias, and score it under spatial-
block CV with the continuous Boyce index").

Presence-background models give *relative* suitability, not the probability of
occurrence (prevalence is unknown), and AUC from pseudo-absences is a poor
yardstick.  The continuous Boyce index (Hirzel et al. 2006) compares the
predicted-to-expected frequency ratio along the suitability gradient.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression


def fit_presence_background(x_presence, x_background, C: float = 1.0):
    """Penalised (L2/ridge) logistic regression distinguishing presence from
    background points -- the standard practical equivalent of a presence-only
    point-process / MaxEnt-style model (infinitely-weighted logistic regression
    converges to the same likelihood as a Poisson point process; Fithian &
    Hastie 2013), used here at finite background weight since that is already
    the standard MaxEnt-equivalent practice at this sample size.

    ``x_presence``/``x_background`` are (n, p) feature matrices on the same
    physically meaningful predictors (e.g. CWD, GDD, winter minimum, VPD,
    soil -- see the concept note). ``class_weight='balanced'`` corrects for
    presence points typically being far outnumbered by background points,
    which would otherwise bias the penalised fit toward "predict background."
    Returns the fitted ``sklearn.linear_model.LogisticRegression``; predicted
    suitability is ``model.predict_proba(x)[:, 1]``.
    """
    x_presence = np.asarray(x_presence, dtype=float)
    x_background = np.asarray(x_background, dtype=float)
    x = np.vstack([x_presence, x_background])
    y = np.concatenate([np.ones(len(x_presence)), np.zeros(len(x_background))])
    # l1_ratio=0 is sklearn's current (>=1.8) spelling of a pure L2/ridge penalty;
    # penalty="l2" still works but is deprecated and warns.
    model = LogisticRegression(l1_ratio=0, C=C, class_weight="balanced", max_iter=1000)
    model.fit(x, y)
    return model


def boyce_index(pred_at_presence, pred_all, n_windows: int = 100, window_width: float = 0.1):
    pp = np.asarray(pred_at_presence, dtype=float)
    pa = np.asarray(pred_all, dtype=float)
    lo, hi = np.min(pa), np.max(pa)
    rng = hi - lo
    if rng <= 0:
        return np.nan
    w = window_width * rng
    starts = np.linspace(lo, hi - w, n_windows)
    pe, mids = [], []
    for s in starts:
        e = s + w
        p = np.mean((pp >= s) & (pp <= e))
        a = np.mean((pa >= s) & (pa <= e))
        if a > 0:
            pe.append(p / a)
            mids.append(0.5 * (s + e))
    if len(pe) < 3:
        return np.nan
    return float(spearmanr(mids, pe).statistic)
