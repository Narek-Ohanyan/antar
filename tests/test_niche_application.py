import numpy as np
import pytest

from antar.niche import application as A


def test_chelsa_keys_follow_the_period_model_and_path_of_a_member():
    assert A.chelsa_key(2050, "gfdl-esm4", "ssp126") == "2041-2070|GFDL-ESM4|ssp126"
    assert A.chelsa_key(2080, "ukesm1-0-ll", "ssp585") == A.chelsa_key(2100, "ukesm1-0-ll", "ssp585") == "2071-2100|UKESM1-0-LL|ssp585"
    with pytest.raises(KeyError):
        A.chelsa_key(2030, "gfdl-esm4", "ssp126")


def test_saturation_vapour_pressure_matches_reference_values():
    assert A.saturation_vapour_pressure_kpa(20.0) == pytest.approx(2.338, abs=0.005)
    assert A.saturation_vapour_pressure_kpa(0.0) == pytest.approx(0.6108)


def test_vpd_rises_about_seven_percent_per_degree_at_constant_humidity_and_not_at_all_without_warming():
    assert A.vpd_after_warming(500.0, 10.0, 0.0) == pytest.approx(500.0)
    up = A.vpd_after_warming(500.0, 10.0, 1.0) / 500.0
    assert 1.06 < up < 1.075
    assert A.vpd_after_warming(500.0, 10.0, 4.0) == pytest.approx(500.0 * A.saturation_vapour_pressure_kpa(14.0) / A.saturation_vapour_pressure_kpa(10.0))
    assert A.vpd_after_warming(np.array([300.0, 600.0]), np.array([5.0, 12.0]), 2.0).shape == (2,)


def test_logistic_score_and_missing_predictors():
    s = A.logistic_score([1.0, -2.0], 0.5, [[0.0, 0.0], [1.0, 1.0], [np.nan, 0.0]])
    assert s[0] == pytest.approx(1 / (1 + np.exp(-0.5))) and s[1] == pytest.approx(1 / (1 + np.exp(0.5)))
    assert np.isnan(s[2])


def test_presence_threshold_leaves_the_stated_share_of_records_below_it():
    scores = np.linspace(0.0, 1.0, 101)
    t = A.presence_threshold(scores, 0.10)
    assert t == pytest.approx(0.10) and (scores < t).mean() == pytest.approx(0.10, abs=0.011)
    assert A.presence_threshold(np.array([0.2, np.nan, 0.4]), 0.0) == pytest.approx(0.2)
    with pytest.raises(ValueError):
        A.presence_threshold(np.array([np.nan]))


def test_the_quality_gate_needs_a_fitted_model_with_enough_skill():
    assert A.passes_gate("fitted", 0.43, 0.3) and not A.passes_gate("fitted", 0.25, 0.3)
    assert not A.passes_gate("skipped_insufficient_data", None, 0.3) and not A.passes_gate("fitted", None, 0.3) and not A.passes_gate("fitted", float("nan"), 0.3)
    assert A.passes_gate("fitted", 0.3, 0.3)


def test_a_group_is_supported_where_any_passing_species_is_and_unassessed_when_none_passes():
    by = {"a": [True, False, False], "b": [False, False, True], "c": [True, True, True]}
    assert A.group_support(by, ["a", "b"]).tolist() == [True, False, True]
    assert A.group_support(by, ["b"]).tolist() == [False, False, True]
    assert A.group_support(by, []) is None and A.group_support(by, ["zzz"]) is None
