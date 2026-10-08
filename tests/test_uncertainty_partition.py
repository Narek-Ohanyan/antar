import numpy as np
import pytest

from antar.uncertainty import partition as P


def _members(f, n_nodes=6, gcms=("a", "b", "c", "d", "e"), ssps=("ssp126", "ssp370", "ssp585"), hzs=(2050, 2080, 2100)):
    """An ensemble in which every member has the same value at all nodes, f(gcm index, ssp index, horizon index)."""
    return {(g, s, h): f(gi, si, hi) * np.ones(n_nodes) for gi, g in enumerate(gcms) for si, s in enumerate(ssps) for hi, h in enumerate(hzs)}


def test_additive_ensemble_splits_exactly_into_its_three_effects():
    a = np.array([0.0, 1.0, -1.0, 2.0, -2.0])            # GCM effect
    b = np.array([0.0, 3.0, 6.0])                        # SSP effect
    c = np.array([0.0, 5.0, 10.0])                       # horizon effect
    vals = _members(lambda g, s, h: a[g] + b[s] + c[h])
    # give every node its own scale so the node dimension is not degenerate
    vals = {k: v * np.linspace(1, 2, 6) for k, v in vals.items()}
    r = P.partition_metric(vals)
    tot = a.var() + b.var() + c.var()
    assert r["pooled"]["gcm"]["median"] == pytest.approx(a.var() / tot)
    assert r["pooled"]["ssp"]["median"] == pytest.approx(b.var() / tot)
    assert r["pooled"]["horizon"]["median"] == pytest.approx(c.var() / tot)
    assert r["pooled"]["interaction"]["median"] == pytest.approx(0.0, abs=1e-12)
    assert r["levels"]["horizon"] == [2050, 2080, 2100] and r["n_nodes_used"] == 6
    for h in (2050, 2080, 2100):                                       # at a fixed horizon only GCM and SSP remain
        bh = r["by_horizon"][h]
        assert bh["gcm"]["median"] == pytest.approx(a.var() / (a.var() + b.var()))
        assert bh["ssp"]["median"] == pytest.approx(b.var() / (a.var() + b.var()))
        assert bh["interaction"]["median"] == pytest.approx(0.0, abs=1e-12)


def test_a_pure_interaction_ends_up_in_the_remainder():
    a = np.array([1.0, -1.0, 1.0, -1.0, 0.0])
    b = np.array([1.0, -1.0, 0.0])
    vals = {k: v * np.linspace(1, 2, 4) for k, v in _members(lambda g, s, h: a[g] * b[s], n_nodes=4).items()}
    r = P.partition_metric(vals)
    assert r["pooled"]["interaction"]["median"] > 0.9
    assert r["pooled"]["gcm"]["median"] < 0.05


def test_fractions_sum_to_one_for_a_random_ensemble():
    rng = np.random.default_rng(3)
    cube = rng.normal(size=(20, 5, 3, 3))
    p = P.partition_nodes(cube)
    total = sum(p["pooled"][k] for k in ("gcm", "ssp", "horizon", "interaction"))
    assert total == pytest.approx(np.ones(20))
    assert p["used"].all()


def test_nodes_without_variation_or_with_missing_values_are_left_out_and_counted():
    cube = np.random.default_rng(4).normal(size=(5, 5, 3, 3))
    cube[1] = 7.0                         # the metric does not vary at this node
    cube[3, 2, 1, 0] = np.nan             # a member is missing here
    p = P.partition_nodes(cube)
    assert p["used"].tolist() == [True, False, True, False, True]
    assert np.isnan(p["pooled"]["gcm"][1]) and np.isnan(p["pooled"]["gcm"][3])
    vals = {(g, s, h): cube[:, gi, si, hi] for gi, g in enumerate("abcde") for si, s in enumerate(("ssp126", "ssp370", "ssp585")) for hi, h in enumerate((2050, 2080, 2100))}
    assert P.partition_metric(vals)["n_nodes_used"] == 3


def test_a_missing_member_is_refused():
    vals = {(g, s, h): np.zeros(3) for g in "ab" for s in ("x", "y") for h in (1, 2)}
    del vals[("b", "y", 2)]
    with pytest.raises(ValueError, match="full factorial"):
        P.cube_from_members(vals)


def test_cube_places_each_member_at_its_factor_levels():
    vals = {(g, s, h): np.full(2, 100 * gi + 10 * si + hi) for gi, g in enumerate("ab") for si, s in enumerate(("p", "q")) for hi, h in enumerate((10, 20))}
    cube, gcms, ssps, hzs = P.cube_from_members(vals)
    assert (gcms, ssps, hzs) == (["a", "b"], ["p", "q"], [10, 20])
    assert cube[0, 1, 0, 1] == 100 + 0 + 1 and cube.shape == (2, 2, 2, 2)


def test_two_way_interactions_are_told_apart_and_everything_sums_to_one():
    rng = np.random.default_rng(5)
    G, S, H = 5, 3, 3
    base = rng.normal(size=(G, S, H))
    f = P.three_way_fractions(base)
    assert sum(f.values()) == pytest.approx(1.0)
    ssp_h = np.array([[0.0, 0.5, 1.0], [0.0, 2.0, 4.0], [0.0, 4.0, 8.0]])         # paths diverge with time: a pure SSP x horizon pattern, centred
    y = np.broadcast_to((ssp_h - ssp_h.mean(0) - ssp_h.mean(1, keepdims=True) + ssp_h.mean())[None], (G, S, H))
    fy = P.three_way_fractions(y)
    assert fy["ssp_x_horizon"] == pytest.approx(1.0) and fy["ssp"] == pytest.approx(0.0, abs=1e-12) and fy["gcm"] == 0.0
    add = np.add.outer(np.add.outer(rng.normal(size=G), rng.normal(size=S)), rng.normal(size=H))
    fa = P.three_way_fractions(add)
    assert fa["remainder"] == pytest.approx(0.0, abs=1e-12) and all(fa[k] == pytest.approx(0.0, abs=1e-12) for k in P.INTERACTIONS)
    assert all(v == 0.0 for v in P.three_way_fractions(np.ones((G, S, H))).values())


def test_pooled_interaction_equals_the_sum_of_its_parts():
    cube = np.random.default_rng(6).normal(size=(10, 5, 3, 3))
    p = P.partition_nodes(cube)["pooled"]
    parts = p["gcm_x_ssp"] + p["gcm_x_horizon"] + p["ssp_x_horizon"] + p["remainder"]
    assert parts == pytest.approx(p["interaction"])
