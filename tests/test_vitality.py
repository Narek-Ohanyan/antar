import numpy as np
import pytest

from antar.hazard import vitality as V


def _panel(rng, n_nodes=40, pix_per_node=5, n_years=10, beta=(-0.5, 0.2), noise=1.0):
    node = np.repeat(np.arange(n_nodes), pix_per_node * n_years)
    pixel = np.repeat(np.arange(n_nodes * pix_per_node), n_years)
    year = np.tile(np.arange(n_years), n_nodes * pix_per_node)
    a = rng.normal(0, 1, (n_nodes, n_years))                      # node-year climate anomaly 1
    b = rng.normal(0, 1, (n_nodes, n_years))
    x = np.column_stack([a[node, year], b[node, year]])
    fe = rng.normal(0, 3, n_nodes * pix_per_node)[pixel]          # large pixel effects the model must remove
    y = fe + x @ np.array(beta) + rng.normal(0, noise, pixel.size)
    return y, x, pixel, node, year


def test_demean_by_removes_group_means():
    out = V.demean_by(np.array([1.0, 3.0, 10.0, 14.0]), ["a", "a", "b", "b"])
    assert out.tolist() == [-1.0, 1.0, -2.0, 2.0]
    m = V.demean_by(np.array([[1.0, 5.0], [3.0, 9.0]]), [1, 1])
    assert m.tolist() == [[-1.0, -2.0], [1.0, 2.0]]


def test_within_ols_recovers_slopes_despite_large_pixel_effects():
    rng = np.random.default_rng(0)
    y, x, pixel, node, year = _panel(rng)
    r = V.within_ols(y, x, pixel, {"node": node, "year": year}, names=["a", "b"])
    truth = np.array([-0.5, 0.2])
    assert np.all(np.abs(r["coef"] - truth) < 3.5 * r["se_conservative"])     # coverage is tested separately, over many draws
    assert r["coef"][0] < -0.3 and r["coef"][1] > 0.05
    assert r["p_value"][0] < 1e-6 and r["n_pixels"] == 200 and r["n_obs"] == 2000
    assert r["se_conservative"][0] == max(r["se"]["node"][0], r["se"]["year"][0])


def test_pixel_constants_do_not_change_the_slopes():
    rng = np.random.default_rng(1)
    y, x, pixel, node, year = _panel(rng)
    base = V.within_ols(y, x, pixel, {"node": node})["coef"]
    shifted = V.within_ols(y + 100.0 * (pixel % 7), x, pixel, {"node": node})["coef"]
    assert shifted == pytest.approx(base, abs=1e-9)


def test_cluster_standard_error_is_close_to_the_classical_one_for_independent_errors():
    rng = np.random.default_rng(2)
    y, x, pixel, node, year = _panel(rng, n_nodes=80)
    r = V.within_ols(y, x, pixel, {"node": node})
    xd = V.demean_by(x, pixel)
    classical = np.sqrt(np.diag(np.linalg.inv(xd.T @ xd)) * (r["n_obs"] - 2 - r["n_pixels"]) ** -1 *
                        np.sum((V.demean_by(y, pixel) - xd @ r["coef"]) ** 2))
    assert r["se"]["node"] == pytest.approx(classical, rel=0.2)


def test_confidence_intervals_cover_the_truth_about_as_often_as_claimed():
    rng = np.random.default_rng(3)
    hits = 0
    reps = 120
    for _ in range(reps):
        y, x, pixel, node, year = _panel(rng, n_nodes=25, pix_per_node=2, noise=2.0)
        r = V.within_ols(y, x, pixel, {"node": node})
        hits += bool(r["ci95_low"][0] <= -0.5 <= r["ci95_high"][0])
    assert 0.88 <= hits / reps <= 0.995


def test_rows_with_missing_values_are_dropped():
    rng = np.random.default_rng(4)
    y, x, pixel, node, year = _panel(rng)
    y2 = y.copy()
    y2[::17] = np.nan
    r = V.within_ols(y2, x, pixel, {"node": node})
    assert r["n_obs"] == int(np.isfinite(y2).sum())


def test_out_of_sample_r2_is_positive_for_a_real_effect_and_not_for_none():
    rng = np.random.default_rng(5)
    y, x, pixel, node, year = _panel(rng, beta=(-0.8, 0.0))
    assert V.leave_one_group_out_r2(y, x, pixel, year)["r2_out_of_sample"] > 0.1
    y0, x0, pixel0, node0, year0 = _panel(rng, beta=(0.0, 0.0))
    assert V.leave_one_group_out_r2(y0, x0, pixel0, year0)["r2_out_of_sample"] < 0.02


def test_interaction_with_a_node_trait_is_detected_with_the_right_sign():
    rng = np.random.default_rng(6)
    y, x, pixel, node, year = _panel(rng, beta=(-0.3, 0.0), noise=0.7)
    stress = V.rank_within(rng.normal(size=40))[node]              # node trait in [0, 1]
    y = y + (-0.6) * stress * x[:, 0]                              # steeper decline where stress is high
    X = np.column_stack([x[:, 0], x[:, 0] * stress])
    r = V.within_ols(y, X, pixel, {"node": node, "year": year}, names=["a", "a_x_stress"])
    assert r["coef"][1] == pytest.approx(-0.6, abs=0.15)
    assert V.interaction_verdict(r["coef"][1], r["ci95_low"][1], r["ci95_high"][1]) == "supports"


def test_interaction_verdict_reads_the_interval_not_the_point_estimate():
    assert V.interaction_verdict(-0.2, -0.4, -0.01) == "supports"
    assert V.interaction_verdict(0.2, 0.01, 0.4) == "contradicts"
    assert V.interaction_verdict(-0.2, -0.5, 0.1) == "inconclusive"
    assert V.interaction_verdict(float("nan"), -1, 1) == "inconclusive"
    assert V.interaction_verdict(0.2, 0.01, 0.4, expected_sign=1) == "supports"


def test_rank_within_is_scaled_to_the_unit_interval_and_handles_a_constant():
    r = V.rank_within([10.0, 30.0, 20.0, 20.0])
    assert r.min() == 0.0 and r.max() == 1.0 and r[2] == r[3]
    assert V.rank_within([5.0, 5.0, 5.0]).tolist() == [0.5, 0.5, 0.5]


# ---- the analysis on a synthetic panel -----------------------------------------------------------------------------------------------------------------

def _table(seed=0, n_nodes=40, per_node=6, years=range(2010, 2020), slope_low=-0.2, slope_high=-0.6, disturbed_share=0.0):
    import pandas as pd
    rng = np.random.default_rng(seed)
    yrs = list(years)
    node_h = rng.uniform(0, 0.5, n_nodes)                                   # XYLEM hazard of each cell
    clim = rng.normal(size=(n_nodes, len(yrs)))                             # a node-year cwd departure
    clim2 = 0.6 * clim + 0.8 * rng.normal(size=(n_nodes, len(yrs)))         # a correlated second climate variable (vpd), not a copy of the first
    rows = []
    pixel = 0
    for n in range(n_nodes):
        rank = (np.argsort(np.argsort(node_h))[n]) / (n_nodes - 1)
        for _ in range(per_node):
            fe = rng.normal(0, 2)
            h = rng.uniform(5, 25)
            for t, y in enumerate(yrs):
                slope = slope_low + (slope_high - slope_low) * rank
                rows.append({"pixel": pixel, "node": n, "year": y, "z": fe + slope * clim[n, t] + rng.normal(0, 0.5), "cwd_mm": 100 + 20 * clim[n, t], "vpd": 1 + 0.1 * clim2[n, t],
                             "height_m": h, "forest_class": 32, "disturbed_recent": rng.random() < disturbed_share})
            pixel += 1
    return pd.DataFrame(rows), {n: float(node_h[n]) for n in range(n_nodes)}


def test_climate_anomalies_are_centred_per_cell_and_scaled_by_the_pooled_departure():
    df, _ = _table()
    out, scales = V.add_climate_anomalies(df, ["cwd_mm", "vpd"])
    assert out.groupby("node")["cwd_mm_anom"].mean().abs().max() < 1e-9
    ny = out.drop_duplicates(["node", "year"])
    assert np.sqrt(np.mean(ny["cwd_mm_anom"] ** 2)) == pytest.approx(1.0)
    assert scales["cwd_mm"] == pytest.approx(20.0, rel=0.15)


def test_the_analysis_recovers_the_drought_response_and_the_xylem_interaction():
    df, hazard = _table(seed=1)
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"], xylem={"ring_porous_oak": hazard}, group_of_class={32: "ring_porous_oak"})
    a = r["model_a"]["terms"][0]
    assert -0.55 < a["coef"] < -0.25 and a["ci95"][1] < 0 and a["p"] < 0.01                  # mean slope is about -0.4
    assert r["leave_one_year_out"]["cwd_only"]["r2_out_of_sample"] > 0.01
    x = r["xylem_check"]["ring_porous_oak"]
    assert x["status"] == "ok" and x["expected_sign"] == "negative"
    assert -0.7 < x["model"]["terms"][1]["coef"] < -0.1 and x["verdict"] == "supports"
    sl = [t["slope"] for t in x["by_hazard_tercile"]]
    assert len(sl) == 3 and sl[0] > sl[2]                                                     # least negative where XYLEM sees the least hazard
    assert len(r["response_by_anomaly_bin"]) == 5
    assert r["response_by_anomaly_bin"][0]["z_mean"] > r["response_by_anomaly_bin"][-1]["z_mean"]      # drier years, lower vitality


def test_no_response_gives_an_inconclusive_check():
    df, hazard = _table(seed=2, slope_low=0.0, slope_high=0.0)
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"], xylem={"ring_porous_oak": hazard}, group_of_class={32: "ring_porous_oak"})
    assert r["xylem_check"]["ring_porous_oak"]["verdict"] == "inconclusive"
    assert r["model_a"]["terms"][0]["ci95"][0] < 0 < r["model_a"]["terms"][0]["ci95"][1]


def test_recently_disturbed_pixel_years_are_left_out():
    df, _ = _table(seed=3, disturbed_share=0.3)
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"])
    assert r["n_excluded_recent_disturbance"] == int(df["disturbed_recent"].sum())
    assert r["n_pixel_years"] == int((~df["disturbed_recent"]).sum())
    assert r["xylem_check"] == {}


def test_a_group_with_too_few_cells_or_no_spread_is_reported_not_forced():
    df, hazard = _table(seed=4, n_nodes=6)
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"], xylem={"ring_porous_oak": hazard}, group_of_class={32: "ring_porous_oak"})
    assert r["xylem_check"]["ring_porous_oak"]["status"] == "too_few_pixels_or_cells"
    df, hazard = _table(seed=5)
    flat = {n: 0.01 for n in hazard}
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"], xylem={"ring_porous_oak": flat}, group_of_class={32: "ring_porous_oak"})
    assert r["xylem_check"]["ring_porous_oak"]["status"] == "no_spread_in_xylem_hazard"


def test_collinear_regressors_are_refused_by_the_regression_and_recorded_by_the_analysis():
    rng = np.random.default_rng(0)
    pixel = np.repeat(np.arange(50), 10)
    x = rng.normal(size=(500, 1))
    with pytest.raises(ValueError, match="collinear"):
        V.within_ols(rng.normal(size=500), np.hstack([x, 2 * x]), pixel, {"node": pixel})
    df, hazard = _table(seed=6)
    df["vpd"] = 1 + 0.1 * (df["cwd_mm"] - 100) / 20                                   # an exact copy of the water-deficit anomaly
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"])
    assert r["model_b"]["status"] == "collinear" and "all_four" not in r["leave_one_year_out"] and r["model_a"]["terms"][0]["coef"] < 0


def test_two_way_demeaning_removes_both_sets_of_means():
    rng = np.random.default_rng(0)
    pixel = np.repeat(np.arange(30), 8)
    year = np.tile(np.arange(8), 30)
    v = rng.normal(size=pixel.size) + 5 * (pixel % 7) + 3 * (year % 3)
    keep = rng.random(pixel.size) > 0.2                          # unbalanced
    out = V.demean_two_way(v[keep], pixel[keep], year[keep])
    assert np.abs(np.bincount(pixel[keep], weights=out)).max() < 1e-6 and np.abs(np.bincount(year[keep], weights=out)).max() < 1e-6


def test_a_year_effect_that_tracks_the_regressor_biases_the_pixel_only_fit_but_not_the_two_way_fit():
    rng = np.random.default_rng(1)
    n_nodes, n_years = 60, 10
    node = np.repeat(np.arange(n_nodes), n_years)
    year = np.tile(np.arange(n_years), n_nodes)
    common = rng.normal(size=n_years)                                          # a country-wide year pattern (a sensor step, a national drought)
    x = common[year] + 0.7 * rng.normal(size=node.size)                        # the regressor follows it, with local differences
    y = 2.0 * common[year] + 0.0 * x + rng.normal(0, 0.3, node.size)           # the true local response is zero; the year pattern is not a response to x
    pixel = node
    one = V.within_ols(y, x[:, None], pixel, {"node": node, "year": year})["coef"][0]
    two = V.within_ols(y, x[:, None], pixel, {"node": node, "year": year}, second_fe=year)["coef"][0]
    assert one > 0.5 and abs(two) < 0.12


def test_the_analysis_reports_how_the_response_holds_up_to_the_year_effect_and_to_poor_composites():
    df, hazard = _table(seed=7)
    df["valid_count"] = np.where(df["year"] < 2014, 3.0, 20.0)
    df["z_adj"] = df["z"] - df.groupby("year")["z"].transform("median")
    r = V.analyse_vitality(df, ["cwd_mm", "vpd"])
    s = r["sensitivity"]
    assert set(s) == {"pixel_and_year_fixed_effects", "composites_with_at_least_5_valid_observations", "years_2014_2019", "year_effect_removed_from_the_series"}
    assert s["years_2014_2019"]["coef"] < 0 and s["composites_with_at_least_5_valid_observations"]["n_obs"] < r["n_pixel_years"]
    assert s["pixel_and_year_fixed_effects"]["coef"] < 0 and s["year_effect_removed_from_the_series"]["coef"] < 0
