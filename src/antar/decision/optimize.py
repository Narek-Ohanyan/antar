"""Robust site-species portfolio via mixed-integer linear programming.

Decision x_{uj} in {0,1}: plant option j (species x provenance x method) in planning unit u.
Benefit in scenario c:  B_c = sum_{u,j} A_u b_{ujc} x_{uj}     (b = viability x value weight).

    maximise   (1 - lam) * sum_c w_c B_c   +   lam * CVaR_alpha(B)
    CVaR_alpha(B) = max_zeta  zeta - 1/(1-alpha) * sum_c w_c * s_c,    s_c >= zeta - B_c,  s_c >= 0
                    (mean of the worst (1-alpha) share of scenarios; Rockafellar & Uryasev 2000)
    s.t.       sum_j x_{uj} <= 1                     one option per unit
               sum_{u,j} A_u cost_{uj} x_{uj} <= budget
               area(option j) <= max_share * total planted area   (response diversity)
               x_{uj} <= e_{uj}                      eligibility (outside applicability region / legally excluded -> 0)
               optional per-basin water-use cap
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


def robust_portfolio(benefit, area, cost, budget, scenario_weights=None, lam: float = 0.5, alpha: float = 0.8,
                     max_share: float = 1.0, min_area: float = 0.0, water_use=None, water_caps=None, basin=None,
                     eligible=None):
    """Solve the portfolio problem.

    benefit : (U, J, C) per-hectare benefit of option j in unit u under scenario c
    area    : (U,) hectares;  cost : (U, J) per-hectare cost;  budget : total cost cap
    water_use : (U, J) per-hectare extra water use; water_caps : {basin_id: cap}; basin : (U,) basin id
    eligible : (U, J) boolean mask; options with False can never be selected
    Returns dict(x (U,J) binary, expected, cvar, scenario_benefit (C,)).
    """
    benefit = np.asarray(benefit, dtype=float)
    U, J, C = benefit.shape
    area = np.asarray(area, dtype=float)
    cost = np.asarray(cost, dtype=float)
    w = np.full(C, 1.0 / C) if scenario_weights is None else np.asarray(scenario_weights, dtype=float) / np.sum(scenario_weights)

    nx = U * J
    n = nx + 1 + C            # x, zeta, s_c
    iz = nx
    is_ = nx + 1
    coef = np.zeros((nx, C))
    for c in range(C):
        coef[:, c] = (benefit[:, :, c] * area[:, None]).reshape(-1)
    exp_coef = coef @ w
    obj = np.zeros(n)
    obj[:nx] = -(1 - lam) * exp_coef
    obj[iz] = -lam
    obj[is_:] = lam * w / (1 - alpha)

    rows, lo, hi = [], [], []
    A_mat = lil_matrix((U + C + 1 + J + (len(water_caps) if water_caps else 0), n))
    r = 0
    for u in range(U):                                   # one option per unit
        A_mat[r, u * J:(u + 1) * J] = 1
        lo.append(-np.inf); hi.append(1.0); r += 1
    for c in range(C):                                   # s_c + B_c >= zeta  ->  zeta - B_c - s_c <= 0
        A_mat[r, :nx] = -coef[:, c]
        A_mat[r, iz] = 1.0
        A_mat[r, is_ + c] = -1.0
        lo.append(-np.inf); hi.append(0.0); r += 1
    A_mat[r, :nx] = (cost * area[:, None]).reshape(-1)   # budget
    lo.append(-np.inf); hi.append(budget); r += 1
    total_area_expr = np.tile(area[:, None], (1, J)).reshape(-1)
    for j in range(J):                                   # diversity: option share cap (linearised vs total area cap)
        col = np.zeros(nx)
        col[j::J] = area
        A_mat[r, :nx] = col - max_share * total_area_expr
        lo.append(-np.inf); hi.append(0.0); r += 1
    if water_caps:
        for b, cap in water_caps.items():
            col = np.zeros(nx)
            for u in range(U):
                if basin[u] == b:
                    col[u * J:(u + 1) * J] = np.asarray(water_use)[u] * area[u]
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
    res = milp(obj, constraints=cons, integrality=integrality, bounds=Bounds(lb, ub))
    if not res.success:
        raise RuntimeError(f"MILP failed: {res.message}")
    x = np.round(res.x[:nx]).reshape(U, J)
    bc = coef.T @ x.reshape(-1)
    srt = np.sort(bc)
    k = max(int(np.ceil((1 - alpha) * C)), 1)
    return {"x": x, "expected": float(bc @ w), "cvar": float(srt[:k].mean()), "scenario_benefit": bc}


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
