"""Continuous vitality response: how a forest pixel's satellite vitality anomaly moves with the weather of the year.

The strict dieback rule needs a persistent collapse, so a short record can hold no event at all. The vitality anomaly (``observation.trailing_z``: kNDVI against the pixel's own previous ten
years) is defined for every pixel-year, so the same panel can still test the direction and size of the response to drought without waiting for events. It is a vitality response, not
mortality: a dry year can lower kNDVI without killing a tree.

Model, per pixel p in year t, with pixel fixed effects (everything that is constant at a pixel is removed):

    z_pt = a_p + sum_k b_k * x_k,nt + (interaction terms) + e_pt

where x_k,nt is the node's within-node climate anomaly (the year's value minus the node's mean over the panel years). Estimation is ordinary least squares on the pixel-demeaned data with
cluster-robust (CR1) standard errors. Climate is shared by every pixel of a node and, in a dry year, by every node, so errors are correlated within a node and within a year; both
clusterings are reported and the larger standard error decides.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def demean_by(values, groups):
    """``values`` (n,) or (n, k) minus the mean of its group; ``groups`` any hashable labels."""
    v = np.asarray(values, dtype=float)
    _, inv = np.unique(np.asarray(groups), return_inverse=True)
    counts = np.bincount(inv)
    if v.ndim == 1:
        means = np.bincount(inv, weights=v) / counts
        return v - means[inv]
    out = np.empty_like(v)
    for j in range(v.shape[1]):
        out[:, j] = v[:, j] - (np.bincount(inv, weights=v[:, j]) / counts)[inv]
    return out


def demean_two_way(values, g1, g2, tol: float = 1e-10, max_iter: int = 500):
    """Remove the means of two sets of groups at once (alternating projections); exact for a balanced panel and convergent for an unbalanced one."""
    v = np.asarray(values, dtype=float).copy()
    for _ in range(max_iter):
        before = v.copy()
        v = demean_by(demean_by(v, g1), g2)
        if np.max(np.abs(v - before)) < tol:
            break
    return v


def _cluster_covariance(Xd, resid, bread, clusters):
    """CR1 sandwich: G/(G-1) * (N-1)/(N-k) * bread * sum_g (X_g' e_g)(X_g' e_g)' * bread."""
    n, k = Xd.shape
    _, inv = np.unique(np.asarray(clusters), return_inverse=True)
    g = inv.max() + 1
    scores = np.zeros((g, k))
    np.add.at(scores, inv, Xd * resid[:, None])
    meat = scores.T @ scores
    c = (g / (g - 1.0)) * ((n - 1.0) / (n - k)) if g > 1 else np.nan
    return c * bread @ meat @ bread, g


def within_ols(y, X, pixel, clusters: dict, names=None, second_fe=None):
    """Fixed-effects regression of ``y`` on ``X`` with pixel effects, and cluster-robust inference for each clustering in ``clusters`` (name -> labels per row).

    Returns coefficients, the standard error under each clustering, ``se_conservative`` (the larger of them), the t-based 95 % interval and two-sided p-value from the
    conservative standard error (degrees of freedom: smallest number of clusters - 1), and the within R^2. Rows with a non-finite value are dropped first.
    ``second_fe`` (labels per row, e.g. the year) adds a second set of fixed effects; the coefficients are then identified from differences between pixels within the same year only.
    """
    y = np.asarray(y, dtype=float)
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != y.size:
        X = X.T
    pixel = np.asarray(pixel)
    keep = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X, pixel = y[keep], X[keep], pixel[keep]
    cl = {k: np.asarray(v)[keep] for k, v in clusters.items()}
    if second_fe is None:
        yd, Xd = demean_by(y, pixel), demean_by(X, pixel)
    else:
        fe2 = np.asarray(second_fe)[keep]
        yd, Xd = demean_two_way(y, pixel, fe2), demean_two_way(X, pixel, fe2)
    xtx = Xd.T @ Xd
    if np.linalg.cond(xtx) > 1e10:
        raise ValueError("the regressors are collinear after removing the pixel means; drop or combine a column")
    bread = np.linalg.inv(xtx)
    beta = bread @ (Xd.T @ yd)
    resid = yd - Xd @ beta
    k = X.shape[1]
    names = list(names) if names is not None else [f"x{j}" for j in range(k)]
    se, n_clusters = {}, {}
    for label, c in cl.items():
        V, g = _cluster_covariance(Xd, resid, bread, c)
        se[label] = np.sqrt(np.diag(V))
        n_clusters[label] = int(g)
    se_cons = np.nanmax(np.vstack(list(se.values())), axis=0)
    dof = max(1, min(n_clusters.values()) - 1)
    tcrit = stats.t.ppf(0.975, dof)
    tval = beta / se_cons
    ss_tot = float(yd @ yd)
    return {
        "names": names, "coef": beta, "se": se, "se_conservative": se_cons, "n_clusters": n_clusters, "dof": dof,
        "ci95_low": beta - tcrit * se_cons, "ci95_high": beta + tcrit * se_cons,
        "p_value": 2.0 * stats.t.sf(np.abs(tval), dof),
        "n_obs": int(y.size), "n_pixels": int(np.unique(pixel).size),
        "r2_within": float(1.0 - (resid @ resid) / ss_tot) if ss_tot > 0 else float("nan"),
    }


def leave_one_group_out_r2(y, X, pixel, group):
    """Out-of-sample R^2 for predicting a held-out group (e.g. a year) of pixel-years: 1 - SSE(model) / SSE(pixel mean of the training years).

    The model prediction for pixel p is its training mean of y plus the fitted slopes times (x - the pixel's training mean of x); the null prediction is the training mean of y alone.
    A pixel absent from the training rows is skipped. Positive means the weather terms predict an unseen year better than the pixel's own average does.
    """
    y = np.asarray(y, dtype=float)
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != y.size:
        X = X.T
    pixel, group = np.asarray(pixel), np.asarray(group)
    keep = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X, pixel, group = y[keep], X[keep], pixel[keep], group[keep]
    sse_model = sse_null = 0.0
    n_test = 0
    for g in np.unique(group):
        te = group == g
        tr = ~te
        yd, Xd = demean_by(y[tr], pixel[tr]), demean_by(X[tr], pixel[tr])
        beta = np.linalg.solve(Xd.T @ Xd, Xd.T @ yd)
        pix_tr, inv = np.unique(pixel[tr], return_inverse=True)
        cnt = np.bincount(inv)
        ybar = np.bincount(inv, weights=y[tr]) / cnt
        xbar = np.column_stack([np.bincount(inv, weights=X[tr][:, j]) / cnt for j in range(X.shape[1])])
        pos = {p: i for i, p in enumerate(pix_tr)}
        idx = np.array([pos.get(p, -1) for p in pixel[te]])
        ok = idx >= 0
        yt, Xt = y[te][ok], X[te][ok]
        pred = ybar[idx[ok]] + (Xt - xbar[idx[ok]]) @ beta
        sse_model += float(np.sum((yt - pred) ** 2))
        sse_null += float(np.sum((yt - ybar[idx[ok]]) ** 2))
        n_test += int(ok.sum())
    return {"r2_out_of_sample": float(1.0 - sse_model / sse_null) if sse_null > 0 else float("nan"), "n_test": n_test}


def interaction_verdict(coef: float, ci_low: float, ci_high: float, expected_sign: int = -1) -> str:
    """Fixed reading of an interaction estimate, decided before looking at it: 'supports' when the 95 % interval lies entirely on the expected side of zero, 'contradicts' when
    entirely on the other side, 'inconclusive' when it contains zero."""
    if not (np.isfinite(coef) and np.isfinite(ci_low) and np.isfinite(ci_high)):
        return "inconclusive"
    if ci_low > 0 or ci_high < 0:
        return "supports" if np.sign(coef) == np.sign(expected_sign) else "contradicts"
    return "inconclusive"


def rank_within(values):
    """Ranks scaled to [0, 1] (average ranks for ties); a constant input maps to 0.5 everywhere."""
    v = np.asarray(values, dtype=float)
    if np.nanmax(v) == np.nanmin(v):
        return np.full(v.shape, 0.5)
    r = stats.rankdata(v, method="average")
    return (r - 1.0) / (len(v) - 1.0)


# ---- the analysis ---------------------------------------------------------------------------------------------------------------------------------------------

def add_climate_anomalies(df, cols, node_col: str = "node", year_col: str = "year"):
    """Add ``<col>_anom``: the node's value in a year minus its mean over the panel years, divided by the standard deviation of those departures over all node-years.

    The mean and the standard deviation come from the unique node-year table, so a node counts once per year however many forest pixels it has; the coefficients then read per
    one standard deviation of the year-to-year departure of the climate of a cell.
    """
    nodeyear = df.drop_duplicates([node_col, year_col])[[node_col, year_col, *cols]]
    out = df.copy()
    scales = {}
    for c in cols:
        dev = nodeyear[c] - nodeyear.groupby(node_col)[c].transform("mean")
        sd = float(np.sqrt(np.mean(dev.to_numpy() ** 2)))
        scales[c] = sd
        mean = nodeyear.groupby(node_col)[c].mean()
        out[f"{c}_anom"] = (out[c] - out[node_col].map(mean)) / sd if sd > 0 else 0.0
    return out, scales


def _coef_table(r):
    rows = []
    for j, name in enumerate(r["names"]):
        rows.append({"term": name, "coef": float(r["coef"][j]), "se_node": float(r["se"]["node"][j]), "se_year": float(r["se"]["year"][j]), "se_used": float(r["se_conservative"][j]),
                     "ci95": [float(r["ci95_low"][j]), float(r["ci95_high"][j])], "p": float(r["p_value"][j])})
    return rows


def _fit(d, y_col, x_cols, names=None):
    r = within_ols(d[y_col].to_numpy(), d[x_cols].to_numpy(), d["pixel"].to_numpy(), {"node": d["node"].to_numpy(), "year": d["year"].to_numpy()}, names=names or x_cols)
    return r, {"n_obs": r["n_obs"], "n_pixels": r["n_pixels"], "n_cells": int(d["node"].nunique()), "n_clusters": r["n_clusters"], "r2_within": r["r2_within"], "terms": _coef_table(r)}


def analyse_vitality(df, climate_cols, xylem=None, group_of_class=None, min_nodes_for_ranking: int = 8, min_pixels_per_group: int = 30):
    """The vitality-response analysis on a pixel-year table with columns pixel, node, year, z, height_m, forest_class, disturbed_recent and the climate columns.

    Models, all with pixel fixed effects and the larger of the node-clustered and year-clustered standard errors:
    A (primary, fixed in advance): z on the climatic-water-deficit anomaly alone; the question is the sign and size of the drought response.
    B: z on all four climate anomalies, with leave-one-year-out R^2 for A and B (does the weather predict an unseen year?).
    C: A plus the interaction of the anomaly with the pixel's canopy height (taller stands more exposed?), exploratory.
    D: per species group (``group_of_class`` maps a forest class to a group): A plus the interaction with the rank, within the group, of XYLEM's mechanistic hazard at the pixel's cell.
       The expectation fixed in advance is a negative interaction (a steeper fall in vitality in a drier year where XYLEM says the hazard is higher); ``interaction_verdict`` reads it.
    Pixel-years with a harvest or fire flagged in the year or in the 10 years before are left out: their baseline is not weather. ``xylem`` is {group: {node: hazard}}.
    """
    anom_cols = [f"{c}_anom" for c in climate_cols]
    full, scales = add_climate_anomalies(df, climate_cols)
    d = full[(~full["disturbed_recent"]) & np.isfinite(full["z"]) & np.all(np.isfinite(full[anom_cols].to_numpy()), axis=1)].copy()
    out = {"n_pixel_years": int(len(d)), "n_pixels": int(d["pixel"].nunique()), "n_cells": int(d["node"].nunique()), "n_excluded_recent_disturbance": int(full["disturbed_recent"].sum()),
           "anomaly_sd": {c: float(v) for c, v in scales.items()}, "mean_z": float(d["z"].mean()), "sd_z": float(d["z"].std())}
    cwd = "cwd_mm_anom"
    rA, out["model_a"] = _fit(d, "z", [cwd], ["cwd anomaly"])
    out["leave_one_year_out"] = {"cwd_only": leave_one_group_out_r2(d["z"].to_numpy(), d[[cwd]].to_numpy(), d["pixel"].to_numpy(), d["year"].to_numpy())}
    try:
        _, out["model_b"] = _fit(d, "z", anom_cols, [c.replace("_anom", " anomaly") for c in anom_cols])
        out["leave_one_year_out"]["all_four"] = leave_one_group_out_r2(d["z"].to_numpy(), d[anom_cols].to_numpy(), d["pixel"].to_numpy(), d["year"].to_numpy())
    except ValueError as e:
        out["model_b"] = {"status": "collinear", "detail": str(e)}

    def brief(sub, y_col="z", second_fe=False):
        if sub["pixel"].nunique() < 30 or sub["node"].nunique() < 8:
            return {"status": "too_few_pixels_or_cells"}
        r = within_ols(sub[y_col].to_numpy(), sub[[cwd]].to_numpy(), sub["pixel"].to_numpy(), {"node": sub["node"].to_numpy(), "year": sub["year"].to_numpy()},
                       second_fe=sub["year"].to_numpy() if second_fe else None)
        return {"coef": float(r["coef"][0]), "ci95": [float(r["ci95_low"][0]), float(r["ci95_high"][0])], "p": float(r["p_value"][0]), "n_obs": r["n_obs"], "n_pixels": r["n_pixels"]}
    sens = {"pixel_and_year_fixed_effects": brief(d, second_fe=True)}
    if "valid_count" in d:
        sens["composites_with_at_least_5_valid_observations"] = brief(d[d["valid_count"] >= 5])
    sens["years_2014_2019"] = brief(d[d["year"] >= 2014])
    if "z_adj" in d:
        sens["year_effect_removed_from_the_series"] = brief(d[np.isfinite(d["z_adj"])], y_col="z_adj")
    out["sensitivity"] = sens

    h = d.drop_duplicates("pixel")["height_m"].to_numpy(dtype=float)
    hs = (np.nanmean(h), np.nanstd(h))
    dh = d[np.isfinite(d["height_m"])].copy()
    if hs[1] > 0 and dh["pixel"].nunique() > 50:
        dh["height_std"] = (dh["height_m"] - hs[0]) / hs[1]
        dh["cwd_x_height"] = dh[cwd] * dh["height_std"]
        _, out["model_c_height"] = _fit(dh, "z", [cwd, "cwd_x_height"], ["cwd anomaly", "cwd anomaly x canopy height (SD)"])
        out["model_c_height"]["height_mean_sd_m"] = [float(hs[0]), float(hs[1])]
    else:
        out["model_c_height"] = {"status": "no_height_spread"}

    # binned response for the chart: node-year means of z by quintile of the anomaly
    ny = d.groupby(["node", "year"]).agg(z=("z", "mean"), a=(cwd, "first")).reset_index()
    edges = np.unique(np.quantile(ny["a"], np.linspace(0, 1, 6)))
    if len(edges) > 2:
        ny["bin"] = np.clip(np.searchsorted(edges, ny["a"], side="right") - 1, 0, len(edges) - 2)
        out["response_by_anomaly_bin"] = [{"anomaly_mean": float(g["a"].mean()), "z_mean": float(g["z"].mean()), "z_se": float(g["z"].std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else None,
                                           "n_cell_years": int(len(g))} for _, g in ny.groupby("bin")]

    xy = {}
    if xylem and group_of_class:
        d["group"] = d["forest_class"].map(group_of_class)
        for g, h_by_node in xylem.items():
            sub = d[d["group"] == g].copy()
            entry = {"n_pixels": int(sub["pixel"].nunique()), "n_cells": int(sub["node"].nunique())}
            if entry["n_pixels"] < min_pixels_per_group or entry["n_cells"] < min_nodes_for_ranking:
                xy[g] = {**entry, "status": "too_few_pixels_or_cells"}
                continue
            nodes = np.array(sorted(sub["node"].unique()))
            hv = np.array([h_by_node.get(int(n), np.nan) for n in nodes], dtype=float)
            ok = np.isfinite(hv)
            nodes, hv = nodes[ok], hv[ok]
            if len(nodes) < min_nodes_for_ranking or len(np.unique(hv)) < 3:
                xy[g] = {**entry, "status": "no_spread_in_xylem_hazard", "n_distinct_hazard_values": int(len(np.unique(hv)))}
                continue
            rank = dict(zip(nodes.tolist(), rank_within(hv).tolist()))
            sub = sub[sub["node"].isin(rank)].copy()
            sub["xylem_rank"] = sub["node"].map(rank)
            sub["cwd_x_rank"] = sub[cwd] * sub["xylem_rank"]
            r, tab = _fit(sub, "z", [cwd, "cwd_x_rank"], ["cwd anomaly", "cwd anomaly x XYLEM hazard rank"])
            inter = tab["terms"][1]
            terciles = []
            cuts = np.quantile(hv, [1 / 3, 2 / 3])
            tercile_of = dict(zip(nodes.tolist(), np.searchsorted(cuts, hv, side="right").tolist()))
            sub["tercile"] = sub["node"].map(tercile_of)
            for t in range(3):
                st = sub[sub["tercile"] == t]
                if st["node"].nunique() >= 3 and st["pixel"].nunique() >= 10:
                    _, tt = _fit(st, "z", [cwd], ["cwd anomaly"])
                    terciles.append({"tercile": t, "hazard_range": [float(hv[[tercile_of[n] == t for n in nodes]].min()), float(hv[[tercile_of[n] == t for n in nodes]].max())],
                                     "n_pixels": tt["n_pixels"], "n_cells": tt["n_cells"], "slope": tt["terms"][0]["coef"], "slope_ci95": tt["terms"][0]["ci95"]})
            xy[g] = {**entry, "status": "ok", "model": tab, "expected_sign": "negative",
                     "verdict": interaction_verdict(inter["coef"], inter["ci95"][0], inter["ci95"][1], expected_sign=-1), "by_hazard_tercile": terciles,
                     "hazard_mean": float(np.mean(hv)), "hazard_sd": float(np.std(hv))}
    out["xylem_check"] = xy
    return out
