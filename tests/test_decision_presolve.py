"""Dropping options that another option beats or ties changes nothing about the optimum, and says what it dropped."""
import numpy as np
import pytest

from antar.decision.optimize import non_dominated_options, robust_portfolio


def problem(seed, U=12, J=6, C=5, duplicate=True):
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.5, 1.0, size=(U, 2, C))                               # two genuinely different groups
    benefit = np.concatenate([base[:, [0], :]] * 3 + [base[:, [1], :]] * 3, axis=1) if duplicate else rng.uniform(0.5, 1.0, size=(U, J, C))
    cost = np.tile(rng.choice([1.0, 2.0, 5.0, 9.0], size=J), (U, 1))
    area = rng.uniform(1.0, 4.0, size=U)
    return benefit, area, cost


def test_a_costlier_copy_of_the_same_benefit_is_dropped_and_the_cheapest_is_kept():
    benefit = np.ones((4, 3, 2))
    cost = np.array([[5.0, 1.0, 9.0]] * 4)
    assert list(non_dominated_options(benefit, cost)) == [1]


def test_an_option_that_is_better_anywhere_is_kept():
    benefit = np.ones((3, 2, 2))
    benefit[1, 0, 1] = 2.0                                  # option 0 is better in one unit and one scenario
    cost = np.array([[5.0, 1.0]] * 3)
    assert list(non_dominated_options(benefit, cost)) == [0, 1]


def test_exact_ties_keep_the_lowest_index_and_eligibility_is_respected():
    benefit, cost = np.ones((3, 3, 2)), np.ones((3, 3))
    assert list(non_dominated_options(benefit, cost)) == [0]
    eligible = np.array([[1, 1, 1], [1, 1, 1], [0, 1, 1]], dtype=bool)         # option 0 is not allowed in unit 2, so it cannot stand in for the others there
    assert list(non_dominated_options(benefit, cost, eligible)) == [1]  # option 1 allowed everywhere option 0 and 2 are, equal in everything: lowest allowed index wins


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("lam", [0.0, 0.5, 1.0])
def test_the_optimum_is_the_same_with_and_without_the_presolve(seed, lam):
    benefit, area, cost = problem(seed)
    budget = 0.4 * float((area[:, None] * cost).max(axis=1).sum())             # a budget that binds
    a = robust_portfolio(benefit, area, cost, budget, lam=lam, presolve=False)
    b = robust_portfolio(benefit, area, cost, budget, lam=lam, presolve=True)
    assert b["n_options_solved"] < benefit.shape[1] and a["n_options_solved"] == benefit.shape[1]
    obj = lambda p: (1 - lam) * p["expected"] + lam * p["cvar"]
    assert obj(b) == pytest.approx(obj(a), rel=1e-7, abs=1e-7)
    assert b["x"].shape == a["x"].shape and (b["x"].sum(axis=1) <= 1).all()
    assert float((b["x"] * area[:, None] * cost).sum()) <= budget + 1e-6        # the plan the presolved problem returns is feasible in the full one


def test_the_presolve_is_skipped_when_an_option_share_cap_couples_the_options():
    benefit, area, cost = problem(1)
    p = robust_portfolio(benefit, area, cost, 1e9, max_share=0.5)
    assert p["n_options_solved"] == benefit.shape[1]


def test_the_result_says_whether_optimality_was_proved():
    benefit, area, cost = problem(2)
    p = robust_portfolio(benefit, area, cost, 1e9)
    assert p["optimal"] is True and (p["mip_gap"] is None or p["mip_gap"] < 1e-6)


# ---- an option may act on only part of a unit
def test_a_per_option_area_matrix_that_repeats_the_unit_area_changes_nothing():
    benefit, area, cost = problem(3)
    budget = 0.4 * float((area[:, None] * cost).max(axis=1).sum())
    a = robust_portfolio(benefit, area, cost, budget, lam=0.5, presolve=False)
    b = robust_portfolio(benefit, np.tile(area[:, None], (1, benefit.shape[1])), cost, budget, lam=0.5, presolve=False)
    assert b["expected"] == pytest.approx(a["expected"], rel=1e-9) and b["cvar"] == pytest.approx(a["cvar"], rel=1e-9)


def test_each_option_is_paid_and_rewarded_for_its_own_area_only():
    benefit = np.ones((1, 2, 2)) * 10.0
    cost = np.array([[1.0, 1.0]])
    area = np.array([[100.0, 40.0]])                        # option 0 acts on 100 ha, option 1 on 40 ha of the same unit
    p = robust_portfolio(benefit, area, cost, 1e9, lam=0.0, presolve=False)
    assert p["x"].tolist() == [[1.0, 0.0]] and p["expected"] == pytest.approx(1000.0)               # benefit = 10 x 100 ha
    tight = robust_portfolio(benefit, area, cost, 50.0, lam=0.0, presolve=False)                    # only option 1's 40 ha can be paid for
    assert tight["x"].tolist() == [[0.0, 1.0]] and tight["expected"] == pytest.approx(400.0)
    assert float((tight["x"] * area * cost).sum()) <= 50.0


@pytest.mark.parametrize("seed", range(4))
def test_the_presolve_is_still_exact_when_options_have_different_areas(seed):
    benefit, area, cost = problem(seed)
    rng = np.random.default_rng(100 + seed)
    A = area[:, None] * rng.uniform(0.2, 1.0, size=(1, benefit.shape[1])) * np.ones((benefit.shape[0], 1))
    A[:, 1] = A[:, 0]                                                           # keep some options comparable so that something can be dominated
    budget = 0.3 * float((A * cost).max(axis=1).sum())
    a = robust_portfolio(benefit, A, cost, budget, lam=0.5, presolve=False)
    b = robust_portfolio(benefit, A, cost, budget, lam=0.5, presolve=True)
    assert (0.5 * b["expected"] + 0.5 * b["cvar"]) == pytest.approx(0.5 * a["expected"] + 0.5 * a["cvar"], rel=1e-7, abs=1e-7)
    assert float((b["x"] * A * cost).sum()) <= budget + 1e-6


# ---- species-group share cap -------------------------------------------------------------------------------------------------------------------

def group_problem(seed, U=14, C=4):
    """Three options of group 0 (the same benefit at different costs) and two of group 1 (a worse benefit): without a cap the plan uses group 0 only."""
    rng = np.random.default_rng(seed)
    b0, b1 = rng.uniform(0.7, 1.0, (U, 1, C)), rng.uniform(0.2, 0.6, (U, 1, C))
    benefit = np.concatenate([b0, b0, b0, b1, b1], axis=1)
    cost = np.tile(np.array([3.0, 1.0, 2.0, 2.0, 1.0]), (U, 1))
    area = rng.uniform(1.0, 3.0, U)
    return benefit, area, cost, np.array([0, 0, 0, 1, 1])


def group_shares(p, area, group):
    planted = (p["x"] * area[:, None]).sum()
    return {int(g): float((p["x"][:, group == g] * area[:, None]).sum() / planted) for g in np.unique(group)}


def test_without_a_cap_the_better_group_takes_everything_and_with_one_it_does_not():
    benefit, area, cost, group = group_problem(0)
    free = robust_portfolio(benefit, area, cost, 1e9, group=group)
    assert group_shares(free, area, group)[0] == pytest.approx(1.0)
    capped = robust_portfolio(benefit, area, cost, 1e9, group=group, group_max_share=0.6)
    sh = group_shares(capped, area, group)
    assert sh[0] <= 0.6 + 1e-9 and sh[1] >= 0.4 - 1e-9
    assert capped["expected"] < free["expected"]                               # diversity costs benefit


def test_a_cap_below_one_over_the_number_of_groups_is_infeasible_unless_nothing_is_planted():
    benefit, area, cost, group = group_problem(1)
    p = robust_portfolio(benefit, area, cost, 1e9, group=group, group_max_share=0.4)       # two groups cannot each stay under 40 % of a plan
    assert p["x"].sum() == 0


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("cap", [0.6, 0.75])               # not 0.5: with two groups that would force their areas to be exactly equal, a subset-sum problem
def test_the_group_cap_optimum_is_the_same_with_and_without_the_presolve(seed, cap):
    benefit, area, cost, group = group_problem(seed)
    budget = 0.5 * float((area[:, None] * cost).max(axis=1).sum())
    a = robust_portfolio(benefit, area, cost, budget, group=group, group_max_share=cap, presolve=False)
    b = robust_portfolio(benefit, area, cost, budget, group=group, group_max_share=cap, presolve=True)
    assert b["n_options_solved"] < a["n_options_solved"] == benefit.shape[1]            # the costlier copies inside each group are dropped
    obj = lambda q: 0.5 * q["expected"] + 0.5 * q["cvar"]
    assert obj(b) == pytest.approx(obj(a), rel=1e-7, abs=1e-7)
    assert max(group_shares(b, area, group).values()) <= cap + 1e-9


def test_dominance_never_crosses_groups_or_unequal_areas_when_a_group_cap_is_active():
    benefit = np.ones((3, 2, 2))
    cost = np.array([[1.0, 5.0]] * 3)
    assert list(non_dominated_options(benefit, cost)) == [0]                                                   # no cap: the dearer copy goes
    assert list(non_dominated_options(benefit, cost, group=np.array([0, 1]), area=np.ones((3, 2)))) == [0, 1]  # different groups: both stay
    assert list(non_dominated_options(benefit, cost, group=np.array([0, 0]), area=np.ones((3, 2)))) == [0]     # same group, same area: dropped
    areas = np.array([[2.0, 1.0]] * 3)
    assert list(non_dominated_options(benefit, cost, group=np.array([0, 0]), area=areas)) == [0, 1]             # the copy acts on a different area, so the swap would move the group total
