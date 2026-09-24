import numpy as np
import pytest
from scipy.ndimage import gaussian_filter
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score

from antar.validation import aoa, conformal, design_based, metrics, splits


def _field(n_side=60, corr_cells=6, seed=0):
    rng = np.random.default_rng(seed)
    f = gaussian_filter(rng.normal(size=(n_side, n_side)), corr_cells, mode="wrap")
    f = (f - f.mean()) / f.std()
    xx, yy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    return xx.ravel().astype(float), yy.ravel().astype(float), f.ravel()


# ----------------------------------------------------------------- CV designs
def test_block_kfold_no_block_overlap_and_buffer():
    x, y, _ = _field(40)
    ids = splits.spatial_block_ids(x, y, 8)
    for tr, te in splits.block_kfold(x, y, 8, k=4, buffer=2.0, seed=1):
        assert set(ids[tr]).isdisjoint(set(ids[te]))
        d = np.min(np.hypot(x[tr][:, None] - x[te][None, :], y[tr][:, None] - y[te][None, :]), axis=1)
        assert d.min() > 2.0


def test_forward_chaining_is_causal():
    years = np.repeat(np.arange(2000, 2015), 3)
    for tr, te in splits.forward_chaining(years, min_train_years=8, horizon=2):
        assert years[tr].max() < years[te].min()


def test_leave_group_out_partitions():
    g = np.array(list("aabbccc"))
    seen = set()
    for grp, tr, te in splits.leave_group_out(g):
        assert set(g[tr]).isdisjoint(set(g[te])) and len(tr) + len(te) == len(g)
        seen.add(grp)
    assert seen == {"a", "b", "c"}


def test_variogram_range_grows_with_correlation_length():
    r = []
    for cc in (3, 8):
        x, y, z = _field(50, cc, seed=2)
        idx = np.random.default_rng(0).choice(len(x), 1200, replace=False)
        c, g = splits.empirical_variogram(x[idx], y[idx], z[idx], n_bins=25, max_dist=40)
        r.append(splits.range_from_variogram(c, g))
    assert r[1] > r[0]


def test_random_cv_inflates_skill_on_spatially_autocorrelated_data():
    """The v1 failure mode in miniature: coordinates-only model, random vs blocked CV."""
    x, y, z = _field(50, 7, seed=3)
    idx = np.random.default_rng(1).choice(len(x), 1500, replace=False)
    x, y, z = x[idx], y[idx], z[idx]
    XY = np.column_stack([x, y])

    def cv_r2(splitter):
        pred = np.full(len(z), np.nan)
        for tr, te in splitter:
            m = RandomForestRegressor(60, min_samples_leaf=3, random_state=0, n_jobs=-1).fit(XY[tr], z[tr])
            pred[te] = m.predict(XY[te])
        ok = ~np.isnan(pred)
        return r2_score(z[ok], pred[ok])

    rng = np.random.default_rng(0)
    perm = rng.permutation(len(z))
    rand_folds = ((np.arange(len(z)) % 5) for _ in [0])
    fold = np.arange(len(z)) % 5
    rand = cv_r2((np.where(fold[perm] != f)[0], np.where(fold[perm] == f)[0]) for f in range(5))
    blocked = cv_r2(splits.block_kfold(x, y, 15, k=5, buffer=3.0, seed=0))
    assert rand > 0.6 and blocked < rand - 0.3


# ----------------------------------------------------------------- metrics
def test_murphy_decomposition_identity():
    rng = np.random.default_rng(0)
    p = rng.beta(1, 8, 20000)
    y = (rng.random(20000) < p).astype(float)
    d = metrics.murphy_decomposition(p, y, n_bins=12)
    assert d["brier_binned"] == pytest.approx(d["reliability"] - d["resolution"] + d["uncertainty"], abs=1e-9)


def test_brier_skill_zero_for_climatology_and_positive_for_informative():
    rng = np.random.default_rng(1)
    p = rng.beta(1, 8, 20000)
    y = (rng.random(20000) < p).astype(float)
    assert metrics.brier_skill(np.full_like(p, y.mean()), y) == pytest.approx(0.0, abs=1e-12)
    assert metrics.brier_skill(p, y) > 0.05


def test_calibration_diagnostics_detect_overconfidence():
    rng = np.random.default_rng(2)
    p = rng.beta(1, 8, 40000)
    y = (rng.random(40000) < p).astype(float)
    s, i = metrics.calibration_slope_intercept(p, y)
    assert s == pytest.approx(1.0, abs=0.08) and abs(i) < 0.1
    logit = np.log(p / (1 - p))
    p_over = 1 / (1 + np.exp(-2.0 * logit))
    s2, _ = metrics.calibration_slope_intercept(p_over, y)
    assert s2 < 0.7
    assert metrics.expected_calibration_error(p, y) < metrics.expected_calibration_error(p_over, y)


def test_r2_can_be_negative_and_is_reported_as_such():
    obs = np.array([1.0, 2.0, 3.0, 4.0])
    assert metrics.r2_vs_mean(np.array([4.0, 3.0, 2.0, 1.0]), obs) < 0


# ----------------------------------------------------------------- AOA
def test_aoa_flags_novel_conditions():
    rng = np.random.default_rng(0)
    Xtr = rng.normal(size=(600, 3))
    folds = rng.integers(0, 5, 600)
    a = aoa.AreaOfApplicability().fit(Xtr, folds)
    assert a.inside(rng.normal(size=(200, 3))).mean() > 0.9
    assert a.inside(rng.normal(size=(200, 3)) + 6.0).mean() < 0.05
    assert np.all(a.dissimilarity_index(np.full((3, 3), 8.0)) > a.threshold_)


# ----------------------------------------------------------------- conformal
def test_split_conformal_coverage():
    rng = np.random.default_rng(0)
    covs = []
    for t in range(200):
        y_cal = rng.normal(0, 1, 400)
        pred_cal = np.zeros(400)
        y_new = rng.normal(0, 1, 500)
        lo, hi = conformal.split_conformal_interval(pred_cal, y_cal, np.zeros(500), alpha=0.1)
        covs.append(metrics.interval_coverage(lo, hi, y_new))
    assert np.mean(covs) >= 0.895


def test_weighted_conformal_reduces_to_unweighted_with_equal_weights():
    rng = np.random.default_rng(1)
    s = np.abs(rng.normal(size=300))
    q_w = conformal.weighted_conformal_quantile(s, np.ones(300), 1.0, alpha=0.1)
    q_u = conformal.conformal_quantile(s, alpha=0.1)
    assert q_w == pytest.approx(q_u, abs=1e-12)


def test_cqr_widens_when_undercovering():
    rng = np.random.default_rng(2)
    y = rng.normal(0, 1, 500)
    lo, hi = np.full(500, -0.3), np.full(500, 0.3)          # far too narrow
    lo2, hi2 = conformal.cqr_interval(lo, hi, y, np.array([-0.3]), np.array([0.3]), alpha=0.1)
    assert hi2[0] > 1.0


# ----------------------------------------------------------------- design-based estimation
def test_design_based_equals_sample_proportion_under_proportional_allocation():
    n_hj = np.array([[45, 5], [10, 40]])       # rows: map class; cols: reference class
    area = np.array([500.0, 500.0])
    r = design_based.stratified_estimates(n_hj, area)
    assert r["p_hat"] == pytest.approx([(45 + 10) / 100, (5 + 40) / 100])
    assert r["overall_accuracy"] == pytest.approx(0.85)
    assert r["area_hat"].sum() == pytest.approx(1000.0)


def test_design_based_corrects_for_unequal_strata():
    # rare class mapped in a tiny stratum but with substantial omission error in the big stratum
    n_hj = np.array([[40, 10], [5, 95]])
    area = np.array([100.0, 9900.0])
    r = design_based.stratified_estimates(n_hj, area)
    assert r["p_hat"][0] == pytest.approx(0.01 * 0.8 + 0.99 * 0.05)
