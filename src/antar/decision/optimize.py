"""Robust site-species portfolio via mixed-integer linear programming.

Decision x_{uj} in {0,1}: plant option j (species x provenance x method) in planning unit u.
Benefit in scenario c:  B_c = sum_{u,j} A_u b_{ujc} x_{uj}     (b = viability x value weight).

    maximise   (1 - lam) * sum_c w_c B_c   +   lam * CVaR_alpha(B)
    CVaR_alpha(B) = max_zeta  zeta - 1/(1-alpha) * sum_c w_c * s_c,    s_c >= zeta - B_c,  s_c >= 0
                    (mean of the worst (1-alpha) share of scenarios; Rockafellar & Uryasev 2000)
    s.t.       sum_j x_{uj} <= 1                     one option per unit
               sum_{u,j} A_u cost_{uj} x_{uj} <= budget
               area(option j) <= max_share * total planted area   (response diversity)
               area(options of group g) <= group_max_share * total planted area   (species-group diversity)
               x_{uj} <= e_{uj}                      eligibility (outside applicability region / legally excluded -> 0)
               optional per-basin water-use cap
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix



def non_dominated_options(benefit, cost, eligible=None, group=None, area=None):
    """Indices of the options worth keeping: drops every option that another option beats or ties.

    Option k makes option j unnecessary when, in every unit where j is allowed, k is allowed too, costs no more per hectare and gives at least the same benefit in
    every scenario (a tie in everything keeps the lower index). Any plan that uses j in a unit can use k there instead without costing more or benefiting less, so
    the optimum is unchanged. Only valid when no constraint couples units to specific options (no option-share cap, no water caps).

    With a species-group cap, pass ``group`` (J,) and ``area`` (U, J): option k then makes option j unnecessary only when both belong to the same group and act on the same area in
    every unit where j is allowed, so the swap leaves every group's planted area, and hence the cap, unchanged.
    """
    benefit = np.asarray(benefit, dtype=float)
    cost = np.asarray(cost, dtype=float)
    U, J, C = benefit.shape
    elig = np.ones((U, J), dtype=bool) if eligible is None else np.asarray(eligible, dtype=bool)
    keep = np.ones(J, dtype=bool)
    grp = None if group is None else np.asarray(group)
    ar = None if area is None else np.asarray(area, dtype=float)
    for j in range(J):
        for k in range(J):
            if k == j or not keep[k]:
                continue
            rows = elig[:, j]
            if not rows.any():
                keep[j] = False
                break
            if grp is not None and (grp[k] != grp[j] or (ar is not None and (ar[rows, k] != ar[rows, j]).any())):
                continue
            if not elig[rows, k].all() or (cost[rows, k] > cost[rows, j]).any() or (benefit[rows, k, :] < benefit[rows, j, :]).any():
                continue
            identical = (cost[rows, k] == cost[rows, j]).all() and (benefit[rows, k, :] == benefit[rows, j, :]).all() and (elig[:, k] == elig[:, j]).all()
            if not identical or k < j:
                keep[j] = False
                break
    return np.flatnonzero(keep)


def robust_portfolio(benefit, area, cost, budget, scenario_weights=None, lam: float = 0.5, alpha: float = 0.8,
                     max_share: float = 1.0, min_area: float = 0.0, water_use=None, water_caps=None, basin=None,
                     eligible=None, presolve: bool = True, time_limit: float | None = None, mip_rel_gap: float | None = None,
                     group=None, group_max_share: float = 1.0):
    """Solve the portfolio problem.

    benefit : (U, J, C) per-hectare benefit of option j in unit u under scenario c
    area    : (U,) hectares of each unit, or (U, J) hectares that option j can act on in unit u (an option may apply to only part of a unit);
              cost : (U, J) per-hectare cost;  budget : total cost cap
    water_use : (U, J) per-hectare extra water use; water_caps : {basin_id: cap}; basin : (U,) basin id
    eligible : (U, J) boolean mask; options with False can never be selected
    group, group_max_share : ``group`` (J,) gives each option a species-group id; the options of one group together may take at most ``group_max_share`` of the planted area (1 = no cap).
    presolve : drop options another option beats or ties (see :func:`non_dominated_options`); exact, and skipped when an option-share cap or water caps couple the options. With a group cap only
               options of the same group on identical areas are compared.
    time_limit : seconds for the MILP solver (None = no limit); mip_rel_gap : stop when the proven bound is within this relative gap of the best plan found (None = solver default).
    The result says whether optimality was proved and the gap that remained.
    Returns dict(x (U,J) binary, expected, cvar, scenario_benefit (C,), optimal (bool), mip_gap, n_options_solved).
    """
    capped = group is not None and group_max_share < 1.0
    if presolve and max_share >= 1.0 and not water_caps:
        a2 = np.asarray(area, dtype=float)
        a2 = a2[:, None] if a2.ndim == 1 else a2
        a2 = a2 * np.ones((1, np.asarray(benefit).shape[1]))
        keep = non_dominated_options(np.asarray(benefit, dtype=float) * a2[:, :, None], np.asarray(cost, dtype=float) * a2, eligible,
                                     group=np.asarray(group) if capped else None, area=a2 if capped else None)       # dominance on what a plan gets and pays in total
        if len(keep) < np.asarray(benefit).shape[1]:
            sub = robust_portfolio(np.asarray(benefit)[:, keep, :], np.asarray(area, dtype=float) if np.ndim(area) == 1 else np.asarray(area)[:, keep], np.asarray(cost)[:, keep], budget, scenario_weights=scenario_weights, lam=lam, alpha=alpha,
                                   max_share=max_share, min_area=min_area, eligible=None if eligible is None else np.asarray(eligible)[:, keep], presolve=False,
                                   time_limit=time_limit, mip_rel_gap=mip_rel_gap,
                                   group=None if group is None else np.asarray(group)[keep], group_max_share=group_max_share)
            x = np.zeros(np.asarray(benefit).shape[:2])
            x[:, keep] = sub["x"]
            return {**sub, "x": x}
    benefit = np.asarray(benefit, dtype=float)
    U, J, C = benefit.shape
    area = np.asarray(area, dtype=float)
    A = area[:, None] * np.ones((1, J)) if area.ndim == 1 else area          # (U, J) hectares each option acts on
    cost = np.asarray(cost, dtype=float)
    w = np.full(C, 1.0 / C) if scenario_weights is None else np.asarray(scenario_weights, dtype=float) / np.sum(scenario_weights)

    nx = U * J
    n = nx + 1 + C            # x, zeta, s_c
    iz = nx
    is_ = nx + 1
    coef = np.zeros((nx, C))
    for c in range(C):
        coef[:, c] = (benefit[:, :, c] * A).reshape(-1)
    exp_coef = coef @ w
    obj = np.zeros(n)
    obj[:nx] = -(1 - lam) * exp_coef
    obj[iz] = -lam
    obj[is_:] = lam * w / (1 - alpha)

    rows, lo, hi = [], [], []
    group_ids = np.unique(np.asarray(group)) if capped else np.array([], dtype=int)
    A_mat = lil_matrix((U + C + 1 + J + len(group_ids) + (len(water_caps) if water_caps else 0), n))
    r = 0
    for u in range(U):                                   # one option per unit
        A_mat[r, u * J:(u + 1) * J] = 1
        lo.append(-np.inf); hi.append(1.0); r += 1
    for c in range(C):                                   # s_c + B_c >= zeta  ->  zeta - B_c - s_c <= 0
        A_mat[r, :nx] = -coef[:, c]
        A_mat[r, iz] = 1.0
        A_mat[r, is_ + c] = -1.0
        lo.append(-np.inf); hi.append(0.0); r += 1
    A_mat[r, :nx] = (cost * A).reshape(-1)               # budget
    lo.append(-np.inf); hi.append(budget); r += 1
    total_area_expr = A.reshape(-1)
    for j in range(J):                                   # diversity: option share cap (linearised vs total area cap)
        col = np.zeros(nx)
        col[j::J] = A[:, j]
        A_mat[r, :nx] = col - max_share * total_area_expr
        lo.append(-np.inf); hi.append(0.0); r += 1
    for gid in group_ids:                                # species-group diversity: the group's area share of everything planted
        col = np.zeros(nx)
        for j in np.flatnonzero(np.asarray(group) == gid):
            col[j::J] = A[:, j]
        A_mat[r, :nx] = col - group_max_share * total_area_expr
        lo.append(-np.inf); hi.append(0.0); r += 1
    if water_caps:
        for b, cap in water_caps.items():
            col = np.zeros(nx)
            for u in range(U):
                if basin[u] == b:
                    col[u * J:(u + 1) * J] = np.asarray(water_use)[u] * A[u]
            A_mat[r, :nx] = col
            lo.append(-np.inf); hi.append(cap); r += 1
    A_mat = A_mat[:r].tocsr()
    cons = [LinearConstraint(A_mat, lo, hi)]
    if min_area > 0:
        cons.append(LinearConstraint(total_area_expr.reshape(1, -1) @ np.eye(nx, n), min_area, np.inf))
    integrality = np.r_[np.ones(nx), np.zeros(1 + C)]
    lb = np.r_[np.zeros(nx), -np.inf, np.zeros(C)]
    xub = np.ones(nx) if eligible is None else np.asarray(eligible, dtype=float).reshape(-1)
    ub = np.r_[xub, np.inf, np.full(C, np.inf)]
    options = {}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    if mip_rel_gap is not None:
        options["mip_rel_gap"] = float(mip_rel_gap)
    res = milp(obj, constraints=cons, integrality=integrality, bounds=Bounds(lb, ub), options=options)
    if res.x is None:
        raise RuntimeError(f"MILP failed: {res.message}")
    x = np.round(res.x[:nx]).reshape(U, J)
    bc = coef.T @ x.reshape(-1)
    srt = np.sort(bc)
    k = max(int(np.ceil((1 - alpha) * C)), 1)
    return {"x": x, "expected": float(bc @ w), "cvar": float(srt[:k].mean()), "scenario_benefit": bc, "optimal": bool(res.status == 0),
            "mip_gap": None if getattr(res, "mip_gap", None) is None else float(res.mip_gap), "n_options_solved": int(J)}


def efficient_frontier(benefit, area, cost, budget, lambdas=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                        scenario_weights=None, alpha: float = 0.8, **kw):
    """Sec. 10.3: sweep lambda from 0 (pure mean) to 1 (pure CVaR) to trace the
    mean-CVaR frontier. Solves :func:`robust_portfolio` once per lambda.

    "The distance between its ends is the price of robustness, expressed in
    units of benefit": ``price_of_robustness`` is the drop in expected benefit
    between the first and last lambda solved -- exactly that distance when the
    default grid (0 to 1) is used, and well-defined but not the note's literal
    quantity for any other grid endpoints.
    """
    lambdas = np.asarray(lambdas, dtype=float)
    plans = [
        robust_portfolio(benefit, area, cost, budget, scenario_weights=scenario_weights, lam=float(lam), alpha=alpha, **kw)
        for lam in lambdas
    ]
    expected = np.array([p["expected"] for p in plans])
    cvar = np.array([p["cvar"] for p in plans])
    return {
        "lambdas": lambdas,
        "expected": expected,
        "cvar": cvar,
        "price_of_robustness": float(expected[0] - expected[-1]),
        "plans": plans,
    }


def extrapolation_footprint(x, benefit, area, inside_aoa, scenario_weights=None):
    """Sec. 10.3: "areas with a large extrapolation footprint are not silently
    planted or excluded: they are flagged, and the plan states how much of its
    expected benefit comes from cells outside the area of applicability."

    A reporting diagnostic on an already-solved plan ``x`` -- it does not
    change the plan. ``inside_aoa``: (U,) or (U, J) boolean.
    """
    x = np.asarray(x, dtype=float)
    benefit = np.asarray(benefit, dtype=float)
    area = np.asarray(area, dtype=float)
    U, J, C = benefit.shape
    w = np.full(C, 1.0 / C) if scenario_weights is None else np.asarray(scenario_weights, dtype=float) / np.sum(scenario_weights)
    expected_per_unit_option = np.einsum("uj,ujc,c->uj", x, benefit, w) * area[:, None]

    inside = np.asarray(inside_aoa, dtype=bool)
    if inside.ndim == 1:
        inside = np.broadcast_to(inside[:, None], (U, J))

    total = float(expected_per_unit_option.sum())
    outside = float(expected_per_unit_option[~inside].sum())
    return {
        "total_expected_benefit": total,
        "outside_aoa_expected_benefit": outside,
        "outside_aoa_share": outside / total if total > 0 else 0.0,
    }


def evaluate_portfolio(x, benefit, area, alpha: float = 0.8, scenario_weights=None):
    """Score a fixed plan on an (independent) scenario set: expected benefit and lower-tail CVaR.

    Use on members that were *not* in the optimisation ensemble; a large gap to the in-sample values signals
    that the plan is over-fitted to its scenario sample (optimiser's curse).
    """
    x = np.asarray(x, dtype=float)
    benefit = np.asarray(benefit, dtype=float)
    area = np.asarray(area, dtype=float)
    C = benefit.shape[2]
    w = np.full(C, 1.0 / C) if scenario_weights is None else np.asarray(scenario_weights, dtype=float) / np.sum(scenario_weights)
    bc = np.einsum("uj,ujc,u->c", x, benefit, area)
    srt = np.argsort(bc)
    cum = np.cumsum(w[srt])
    k = int(np.searchsorted(cum, 1 - alpha - 1e-12) + 1)
    return {"expected": float(bc @ w), "cvar": float(np.average(bc[srt][:k], weights=w[srt][:k])), "scenario_benefit": bc}
