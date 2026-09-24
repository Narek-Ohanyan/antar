"""Monotone surrogate of the hydraulic engine.

The mechanistic engine is run on Latin-hypercube samples that deliberately
extend *beyond* the observed climate range; a monotone gradient-boosting model
is then trained on those simulator runs.  Because the training domain was
designed to cover future, non-analogue conditions, the surrogate does not
extrapolate - it interpolates inside a domain we chose.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import qmc
from sklearn.ensemble import HistGradientBoostingRegressor


def latin_hypercube(bounds: dict[str, tuple[float, float]], n: int, seed: int = 0):
    names = list(bounds)
    lo = np.array([bounds[k][0] for k in names])
    hi = np.array([bounds[k][1] for k in names])
    sample = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    return names, qmc.scale(sample, lo, hi)


def fit_monotone_emulator(X, y, increasing: list[int] | None = None, decreasing: list[int] | None = None, **kw):
    """HistGradientBoostingRegressor with monotone constraints on the listed feature indices."""
    X = np.asarray(X, dtype=float)
    cst = np.zeros(X.shape[1], dtype=int)
    for j in increasing or []:
        cst[j] = 1
    for j in decreasing or []:
        cst[j] = -1
    params = dict(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, random_state=0)
    params.update(kw)
    return HistGradientBoostingRegressor(monotonic_cst=cst.tolist(), **params).fit(X, y)


def build_hazard_emulator(
    bounds: dict[str, tuple[float, float]],
    simulate_hazard,
    n_samples: int,
    increasing: list[str],
    decreasing: list[str],
    seed: int = 0,
    **kw,
):
    """The Sec. 6.3 "Emulation" recipe: LHS over named axes, run the full engine, fit monotone.

    ``bounds`` names both the stress-spell-summary axes (e.g.
    ``longest_closed_spell_days``, ``mean_vpd24_during_spell``,
    ``mean_t_during_spell`` -- see :mod:`antar.hydraulics.spells`) and the
    trait-group descriptor axes (e.g. ``capacitance``, ``buffer``).
    ``simulate_hazard`` receives one row as ``{name: value}`` and returns
    h_mech (or a hazard proxy) for it. ``increasing``/``decreasing`` name the
    axes with a known monotonic direction; unlisted axes are unconstrained.
    Returns ``(model, names)`` -- ``names`` gives the column order ``model``
    expects, needed to build matching prediction rows later.
    """
    names, X = latin_hypercube(bounds, n_samples, seed=seed)
    y = np.array([simulate_hazard(dict(zip(names, row))) for row in X])
    inc_idx = [names.index(k) for k in increasing]
    dec_idx = [names.index(k) for k in decreasing]
    model = fit_monotone_emulator(X, y, increasing=inc_idx, decreasing=dec_idx, **kw)
    return model, names


def emulator_max_abs_error(model, names, simulate_hazard, bounds, n_holdout: int = 500, seed: int = 1) -> float:
    """Sec. 6.3's own release gate: max absolute hazard error on a held-out hypercube.

    The concept note requires this below 0.02 before the emulator replaces the
    full engine anywhere.
    """
    _, X = latin_hypercube(bounds, n_holdout, seed=seed)
    y_true = np.array([simulate_hazard(dict(zip(names, row))) for row in X])
    y_hat = model.predict(X)
    return float(np.max(np.abs(y_hat - y_true)))
