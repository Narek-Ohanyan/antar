import numpy as np
import pytest

from antar.decision.optimize import evaluate_portfolio, robust_portfolio
from antar.io.grid import audit_dvpd_dtmax_ratio
from antar.uncertainty import variance


def test_anova_fractions_additive_model():
    a = np.array([0.0, 2.0, 4.0])[:, None, None]
    b = np.array([0.0, 1.0])[None, :, None]
    c = np.array([0.0, 0.5, 1.0, 1.5])[None, None, :]
    Y = a + b + c
    f = variance.anova_fractions(Y, ["ssp", "gcm", "downscaling"])
    assert sum(f[k] for k in ("ssp", "gcm", "downscaling")) == pytest.approx(1.0, abs=1e-9)
    assert f["interaction"] == pytest.approx(0.0, abs=1e-9)
    assert f["ssp"] > f["gcm"]


def test_sobol_ishigami():
    def f(X):
        return np.sin(X[:, 0]) + 7 * np.sin(X[:, 1]) ** 2 + 0.1 * X[:, 2] ** 4 * np.sin(X[:, 0])

    s1, st = variance.sobol_indices(f, [(-np.pi, np.pi)] * 3, n=60000, seed=0)
    assert s1 == pytest.approx([0.314, 0.442, 0.0], abs=0.05)
    assert st == pytest.approx([0.558, 0.442, 0.244], abs=0.05)


def _toy(U=3):
    # option 0: safe (10/ha in every scenario); option 1: higher mean but collapses in one of four scenarios
    b = np.zeros((U, 2, 4))
    b[:, 0, :] = 10.0
    b[:, 1, :] = [14.0, 14.0, 14.0, 0.0]
    return b


def test_portfolio_risk_neutral_takes_higher_mean_and_cvar_takes_safe_option():
    b = _toy()
    area, cost = np.ones(3), np.ones((3, 2))
    neutral = robust_portfolio(b, area, cost, budget=10, lam=0.0, alpha=0.75)
    robust = robust_portfolio(b, area, cost, budget=10, lam=1.0, alpha=0.75)
    assert neutral["x"][:, 1].sum() == 3
    assert robust["x"][:, 0].sum() == 3
    assert robust["cvar"] > neutral["cvar"]
    assert neutral["expected"] > robust["expected"]


def test_portfolio_diversification_constraint_and_budget():
    b = _toy(4)
    b[:, 1, :] = 11.0                                  # option 1 now dominates in every scenario
    area, cost = np.ones(4), np.ones((4, 2))
    r = robust_portfolio(b, area, cost, budget=10, lam=0.0, max_share=0.5)
    assert r["x"][:, 1].sum() <= 2 and r["x"][:, 0].sum() >= 1
    tight = robust_portfolio(b, area, cost, budget=2, lam=0.0)
    assert tight["x"].sum() <= 2


def test_v1_unit_audit_flags_the_training_table():
    ok = audit_dvpd_dtmax_ratio(0.057, 0.0686)
    assert not ok.ok
    assert audit_dvpd_dtmax_ratio(0.12, 1.0).ok


def test_portfolio_eligibility_mask_is_respected():
    b = _toy()
    b[:, 1, :] = 20.0                                   # option 1 dominates ...
    area, cost = np.ones(3), np.ones((3, 2))
    elig = np.ones((3, 2), dtype=bool)
    elig[0, 1] = False                                  # ... but is forbidden in unit 0
    r = robust_portfolio(b, area, cost, budget=10, lam=0.0, eligible=elig)
    assert r["x"][0, 1] == 0 and r["x"][0, 0] == 1
    assert r["x"][1:, 1].sum() == 2


def test_evaluate_portfolio_matches_optimiser_and_detects_overfit():
    b = _toy()
    area, cost = np.ones(3), np.ones((3, 2))
    r = robust_portfolio(b, area, cost, budget=10, lam=1.0, alpha=0.75)
    ev = evaluate_portfolio(r["x"], b, area, alpha=0.75)
    assert ev["expected"] == pytest.approx(r["expected"]) and ev["cvar"] == pytest.approx(r["cvar"])
    worse = b.copy()
    worse[:, 0, :] *= 0.5                                # a harsher, independent scenario set
    assert evaluate_portfolio(r["x"], worse, area, alpha=0.75)["cvar"] < ev["cvar"]
