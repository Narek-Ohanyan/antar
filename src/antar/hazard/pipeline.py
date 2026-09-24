"""MNEME's Sec. 7.3 assembly: design matrix, three learners, nested-CV stacking, and the gate.

:func:`build_hazard_design` assembles Eq. 7.3's linear predictor from labelled
pieces -- age spline, Mundlak within/between climate columns, DLNM drought-
legacy cross-basis, stand terms -- as plain column concatenation (no hidden
decision beyond what Eq. 7.3 itself specifies). Coordinates and raw elevation
are deliberately not accepted: "anything geographic must reach the model
through a physical driver, or a place label that cannot be projected."

:func:`fit_stacked_hazard` is Sec. 7.3's learner recipe: an elastic-net cloglog
GLM and a monotone-constrained GBM, combined by NNLS on out-of-fold link-scale
predictions inside nested blocked CV (no block ever contributes to its own
out-of-fold prediction, so the stacking weights are never fitted on data they
are scored on); a random forest is fit alongside for benchmark comparison only
-- Sec. 7.3 is explicit that RF "cannot extrapolate" and it is never part of
the stack.

:func:`gate_hazard` composes the already-existing, already-tested
``antar.validation.aoa.AreaOfApplicability`` and ``antar.hazard.gating`` into
Eq. 7.4's final blend in one call.

What this module does **not** do: fit any of this against real data. Every
piece here is exercised in its own tests against synthetic panels, proving the
wiring is correct; real skill (and therefore real use) waits on the real
dieback-labelled panel this session is currently pulling via
``antar.io.gee_export`` and ``antar.hazard.observation``.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier

from .dlnm import crossbasis, spline_basis
from .models import CloglogGLM, cloglog, cloglog_inverse, monotone_hgb, stack_nnls
from .panel import mundlak_decompose


def build_hazard_design(
    age,
    age_knots,
    climate_covariates: dict,
    place_id,
    z_lagged: dict | None = None,
    stand_terms=None,
    stand_term_names: list[str] | None = None,
    age_degree: int = 3,
):
    """Assemble Eq. 7.3's linear predictor design matrix.

    ``climate_covariates``: ``{name: values}``, each split into Mundlak
    within/between columns by ``place_id`` (Sec. 7.2). ``z_lagged`` (optional):
    ``{name: (lagged_array, exposure_knots, lag_knots)}``, each turned into a
    DLNM cross-basis (drought legacies, Sec. 7.2). ``stand_terms`` (optional):
    an already-assembled (n, k) array of stand/size x stress-interaction
    columns (gamma^T u in Eq. 7.3) -- signed a priori, not derived here.

    Returns ``(X, names)``.
    """
    import pandas as pd

    n = len(np.asarray(age))
    blocks = [spline_basis(age, age_knots, degree=age_degree)]
    names = [f"age_spline_{j}" for j in range(blocks[0].shape[1])]

    if climate_covariates:
        df = pd.DataFrame({"place": place_id, **{k: np.asarray(v, dtype=float) for k, v in climate_covariates.items()}})
        df = mundlak_decompose(df, "place", list(climate_covariates))
        for k in climate_covariates:
            blocks.append(df[[f"{k}_dev", f"{k}_bar"]].to_numpy())
            names += [f"{k}_within", f"{k}_between"]

    for zname, (z_lag, exposure_knots, lag_knots) in (z_lagged or {}).items():
        cb = crossbasis(z_lag, exposure_knots, lag_knots)
        blocks.append(cb)
        names += [f"{zname}_dlnm_{j}" for j in range(cb.shape[1])]

    if stand_terms is not None:
        stand_terms = np.asarray(stand_terms, dtype=float)
        if stand_terms.ndim == 1:
            stand_terms = stand_terms.reshape(-1, 1)
        blocks.append(stand_terms)
        names += stand_term_names or [f"stand_{j}" for j in range(stand_terms.shape[1])]

    X = np.column_stack(blocks)
    assert X.shape[0] == n
    return X, names


def _binned_cloglog_rate(link_for_binning, y, n_bins: int = 20):
    """A smoothed empirical link target for NNLS stacking (same quantile-binning
    idiom as antar.validation.metrics.murphy_decomposition): the cloglog of the
    empirical event rate within each bin of a prediction, not the (undefined at
    y in {0,1}) per-observation cloglog of the raw label.
    """
    link_for_binning = np.asarray(link_for_binning, dtype=float)
    y = np.asarray(y, dtype=float)
    edges = np.quantile(link_for_binning, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    edges = np.unique(edges)
    idx = np.clip(np.digitize(link_for_binning, edges[1:-1]), 0, len(edges) - 2)
    target = np.empty_like(y)
    for b in range(len(edges) - 1):
        m = idx == b
        if m.any():
            rate = np.clip(y[m].mean(), 1e-4, 1 - 1e-4)
            target[m] = cloglog(rate)
    return target


def fit_stacked_hazard(
    X,
    y,
    blocks,
    offset=None,
    sample_weight=None,
    glm_kwargs: dict | None = None,
    gbm_kwargs: dict | None = None,
    gbm_increasing=(),
    gbm_decreasing=(),
    rf_kwargs: dict | None = None,
    seed: int = 0,
) -> dict:
    """Sec. 7.3's three learners, stacked. ``blocks`` assigns each row to a CV block
    (e.g. from :func:`antar.validation.splits.block_kfold`'s fold ids). Returns a
    dict with the refit-on-all-data ``glm``/``gbm``/``rf`` models, the NNLS
    ``stack_weights``, and out-of-fold link/probability predictions
    (``oof_link_glm``, ``oof_link_gbm``, ``oof_link_stat``, ``oof_pred_stat``,
    ``oof_pred_rf``) -- the OOF stacked prediction is what Sec. 9's proper
    scoring rules should be evaluated against, not an in-sample refit.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    off = np.zeros(n) if offset is None else np.asarray(offset, dtype=float)
    sw = np.ones(n) if sample_weight is None else np.asarray(sample_weight, dtype=float)
    blocks = np.asarray(blocks)
    glm_kwargs = glm_kwargs or {}
    gbm_kwargs = gbm_kwargs or {}
    rf_kwargs = dict(n_estimators=500, random_state=seed) | (rf_kwargs or {})

    oof_link_glm = np.full(n, np.nan)
    oof_link_gbm = np.full(n, np.nan)
    oof_pred_rf = np.full(n, np.nan)
    for b in np.unique(blocks):
        tr, te = blocks != b, blocks == b
        if not tr.any() or not te.any():
            continue
        glm_f = CloglogGLM(**glm_kwargs).fit(X[tr], y[tr], offset=off[tr], sample_weight=sw[tr])
        oof_link_glm[te] = glm_f.linear_predictor(X[te], offset=off[te])

        gbm_f = monotone_hgb(X.shape[1], increasing=gbm_increasing, decreasing=gbm_decreasing, random_state=seed, **gbm_kwargs)
        gbm_f.fit(X[tr], y[tr], sample_weight=sw[tr])
        p_te = np.clip(gbm_f.predict_proba(X[te])[:, 1], 1e-10, 1 - 1e-10)
        oof_link_gbm[te] = cloglog(p_te)

        rf_f = RandomForestClassifier(**rf_kwargs).fit(X[tr], y[tr], sample_weight=sw[tr])
        oof_pred_rf[te] = rf_f.predict_proba(X[te])[:, 1]

    valid = ~np.isnan(oof_link_glm) & ~np.isnan(oof_link_gbm)
    target = _binned_cloglog_rate(0.5 * (oof_link_glm[valid] + oof_link_gbm[valid]), y[valid])
    stack_weights = stack_nnls(np.column_stack([oof_link_glm[valid], oof_link_gbm[valid]]), target)

    oof_link_stat = np.full(n, np.nan)
    oof_link_stat[valid] = stack_weights[0] * oof_link_glm[valid] + stack_weights[1] * oof_link_gbm[valid]
    oof_pred_stat = np.full(n, np.nan)
    oof_pred_stat[valid] = cloglog_inverse(oof_link_stat[valid])

    glm_full = CloglogGLM(**glm_kwargs).fit(X, y, offset=off, sample_weight=sw)
    gbm_full = monotone_hgb(X.shape[1], increasing=gbm_increasing, decreasing=gbm_decreasing, random_state=seed, **gbm_kwargs)
    gbm_full.fit(X, y, sample_weight=sw)
    rf_full = RandomForestClassifier(**rf_kwargs).fit(X, y, sample_weight=sw)

    return {
        "glm": glm_full,
        "gbm": gbm_full,
        "rf": rf_full,
        "stack_weights": stack_weights,
        "oof_link_glm": oof_link_glm,
        "oof_link_gbm": oof_link_gbm,
        "oof_link_stat": oof_link_stat,
        "oof_pred_stat": oof_pred_stat,
        "oof_pred_rf": oof_pred_rf,
    }


def gate_hazard(aoa, X_new, h_stat, h_mech, kappa: float = 2.0):
    """Eq. 7.4 end to end: DI from the fitted AOA, gate weight, cloglog-scale blend.

    ``aoa`` is a fitted :class:`antar.validation.aoa.AreaOfApplicability`.
    Thin composition of already-tested pieces -- adds no new hazard logic.
    """
    from .gating import aoa_weight, blend_hazards

    di = aoa.dissimilarity_index(X_new)
    w = aoa_weight(di, aoa.threshold_, kappa=kappa)
    return blend_hazards(h_stat, h_mech, w)
