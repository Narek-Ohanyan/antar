"""Where the spread of the scenario ensemble comes from.

The ensemble is a full factorial of climate model (GCM) x emissions path (SSP) x horizon. For every grid node and every metric (treeline shift, warming, precipitation change, ...)
the variance across the 45 members is split into the main effects of the three factors and a remainder (:func:`antar.uncertainty.variance.anova_fractions`; Hawkins & Sutton 2009,
Lehner et al. 2020). Two readings are reported: the pooled three-factor split, in which the horizon is itself a source of variance (the signal grows with time), and, at each fixed
horizon, the two-factor split between GCM and SSP, which is the split that matters for a decision at that date.

The ensemble holds one run per GCM and path, so internal variability is not separated from model response: it sits in the GCM main effect and the remainder. Fractions are computed
per node and summarised by their median and interquartile range over nodes; nodes where the metric does not vary at all are left out and counted.
"""
from __future__ import annotations

import numpy as np

from antar.uncertainty.variance import anova_fractions

FACTORS = ("gcm", "ssp", "horizon")
INTERACTIONS = ("gcm_x_ssp", "gcm_x_horizon", "ssp_x_horizon")


def cube_from_members(values: dict):
    """Arrange ``{(gcm, ssp, horizon): array over nodes}`` into a cube (nodes, gcm, ssp, horizon) and return it with the sorted factor levels.

    Raises if the ensemble is not a full factorial, because the variance split is only defined for a complete, balanced design.
    """
    gcms = sorted({k[0] for k in values})
    ssps = sorted({k[1] for k in values})
    hzs = sorted({k[2] for k in values})
    missing = [(g, s, h) for g in gcms for s in ssps for h in hzs if (g, s, h) not in values]
    if missing:
        raise ValueError(f"the ensemble is not a full factorial: {len(missing)} member(s) missing, e.g. {missing[0]}")
    n = len(next(iter(values.values())))
    cube = np.empty((n, len(gcms), len(ssps), len(hzs)))
    for gi, g in enumerate(gcms):
        for si, s in enumerate(ssps):
            for hi, h in enumerate(hzs):
                v = np.asarray(values[(g, s, h)], dtype=float)
                if v.shape != (n,):
                    raise ValueError(f"member {(g, s, h)} has shape {v.shape}, expected ({n},)")
                cube[:, gi, si, hi] = v
    return cube, gcms, ssps, hzs


def three_way_fractions(y):
    """Variance fractions of a balanced full-factorial 3-D array: the three main effects, the three two-way interactions and the three-way remainder.

    With population variances, V_A = Var(mean of y over the other axes), V_AB = Var(mean over the third axis) - V_A - V_B, and the remainder is what is left of Var(y). The nine
    fractions sum to one. A constant array gives all zeros.
    """
    y = np.asarray(y, dtype=float)
    total = y.var()
    if total <= 0:
        return {k: 0.0 for k in (*FACTORS, *INTERACTIONS, "remainder")}
    main = [y.mean(axis=tuple(j for j in range(3) if j != i)).var() for i in range(3)]
    two = {}
    for (i, j), name in (((0, 1), "gcm_x_ssp"), ((0, 2), "gcm_x_horizon"), ((1, 2), "ssp_x_horizon")):
        k = 3 - i - j
        two[name] = y.mean(axis=k).var() - main[i] - main[j]
    out = {n: m / total for n, m in zip(FACTORS, main)}
    out.update({n: max(0.0, v) / total for n, v in two.items()})
    out["remainder"] = max(0.0, 1.0 - sum(out.values()))
    return out


def partition_nodes(cube, min_total_variance: float = 1e-12):
    """Variance fractions per node. ``cube`` is (nodes, gcm, ssp, horizon).

    Returns ``{"pooled": {gcm, ssp, horizon, interaction (all interactions), gcm_x_ssp, gcm_x_horizon, ssp_x_horizon, remainder (three-way): (nodes,)}, "by_horizon": {index: {gcm, ssp, interaction: (nodes,)}}, "used": bool (nodes,)}``. Nodes with a non-finite
    value, or whose total variance is below ``min_total_variance``, have ``used`` False and NaN fractions.
    """
    cube = np.asarray(cube, dtype=float)
    n, G, S, H = cube.shape
    used = np.all(np.isfinite(cube.reshape(n, -1)), axis=1)
    pooled = {k: np.full(n, np.nan) for k in (*FACTORS, "interaction", *INTERACTIONS, "remainder")}
    by_h = {h: {k: np.full(n, np.nan) for k in ("gcm", "ssp", "interaction")} for h in range(H)}
    for i in range(n):
        if not used[i] or cube[i].var() < min_total_variance:
            used[i] = False
            continue
        f = anova_fractions(cube[i], list(FACTORS))
        f.update(three_way_fractions(cube[i]))
        f["interaction"] = 1.0 - f["gcm"] - f["ssp"] - f["horizon"]       # all interactions together, as in anova_fractions
        for k in pooled:
            pooled[k][i] = f[k]
        for h in range(H):
            y = cube[i, :, :, h]
            if y.var() < min_total_variance:
                continue
            fh = anova_fractions(y, ["gcm", "ssp"])
            for k in by_h[h]:
                by_h[h][k][i] = fh[k]
    return {"pooled": pooled, "by_horizon": by_h, "used": used}


def summarise(fractions: dict) -> dict:
    """Median and interquartile range over the nodes for which the fraction is defined."""
    out = {}
    for k, v in fractions.items():
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v)]
        out[k] = {"median": float(np.median(v)), "q25": float(np.percentile(v, 25)), "q75": float(np.percentile(v, 75))} if v.size else None
    return out


def partition_metric(values: dict, min_total_variance: float = 1e-12) -> dict:
    """Full partition of one metric: ``values`` is ``{(gcm, ssp, horizon): array over nodes}``."""
    cube, gcms, ssps, hzs = cube_from_members(values)
    p = partition_nodes(cube, min_total_variance)
    return {
        "levels": {"gcm": gcms, "ssp": ssps, "horizon": hzs},
        "n_nodes": int(cube.shape[0]), "n_nodes_used": int(p["used"].sum()),
        "pooled": summarise(p["pooled"]),
        "by_horizon": {int(hzs[h]): summarise(v) for h, v in p["by_horizon"].items()},
        "ensemble_sd": float(np.nanstd(cube[p["used"]])) if p["used"].any() else None,
    }
