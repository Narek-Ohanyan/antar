import numpy as np
import pandas as pd
import pytest

from antar.hazard import dlnm, gating, models, observation, panel, pipeline
from antar.validation.aoa import AreaOfApplicability


def test_mundlak_decomposition():
    df = pd.DataFrame({"id": [1, 1, 1, 2, 2, 2], "x": [1.0, 2, 3, 10, 12, 14]})
    out = panel.mundlak_decompose(df, "id", ["x"])
    assert out.groupby("id")["x_dev"].mean().abs().max() < 1e-12
    assert out.loc[out.id == 2, "x_bar"].iloc[0] == pytest.approx(12.0)
    assert (out["x_bar"] + out["x_dev"]).equals(out["x"])


def test_person_period_censoring_and_events():
    pp = panel.person_period([1, 2, 3], [2000, 2000, 2000], [2004, 2004, 2004], [np.nan, 2002, 2004])
    g = pp.groupby("id")
    assert g.size().to_dict() == {1: 5, 2: 3, 3: 5}
    assert g["event"].sum().to_dict() == {1: 0, 2: 1, 3: 1}
    assert pp[pp.id == 2].event.tolist() == [0, 0, 1]


def test_crossbasis_partition_of_unity_and_shape():
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, (50, 5))
    cb = dlnm.crossbasis(x, exposure_knots=[-0.5, 0.5], lag_knots=[2.0], degree_x=2, degree_lag=2)
    assert cb.shape == (50, (2 + 2 + 1) * (1 + 2 + 1))
    assert np.allclose(cb.sum(axis=1), 5.0)       # sum over (j,k) of b_j c_k = number of lags


def test_lag_matrix():
    s = np.arange(10.0)
    m = dlnm.lag_matrix(s, 2)
    assert m.shape == (8, 3) and m[0].tolist() == [2.0, 1.0, 0.0]


def test_cloglog_roundtrip():
    p = np.array([0.001, 0.05, 0.5, 0.9])
    assert models.cloglog_inverse(models.cloglog(p)) == pytest.approx(p)


def test_cloglog_glm_recovers_parameters_with_offset():
    rng = np.random.default_rng(42)
    n = 200_000
    X = rng.normal(size=(n, 2))
    exposure = rng.choice([1.0, 2.0], n)
    beta0, beta = -4.0, np.array([1.0, -0.5])
    eta = beta0 + X @ beta + np.log(exposure)
    y = rng.random(n) < models.cloglog_inverse(eta)
    m = models.CloglogGLM(l2=1e-8).fit(X, y, offset=np.log(exposure))
    assert m.intercept_ == pytest.approx(beta0, abs=0.1)
    assert m.coef_ == pytest.approx(beta, abs=0.1)


def test_weighted_penalized_wls_matches_ridge_closed_form_at_l1_zero():
    rng = np.random.default_rng(0)
    n, p = 200, 4
    A = np.column_stack([np.ones(n), rng.normal(size=(n, p - 1))])
    z = rng.normal(size=n)
    w = rng.uniform(0.5, 2.0, n)
    beta_cd = models._weighted_penalized_wls(A, z, w, l1=0.0, l2=0.7, fit_intercept=True, beta_init=np.zeros(p))
    pen = 0.7 * np.eye(p)
    pen[0, 0] = 0.0
    H = A.T @ (A * w[:, None]) + pen
    beta_direct = np.linalg.solve(H, A.T @ (w * z))
    assert beta_cd == pytest.approx(beta_direct, abs=1e-8)


def test_weighted_penalized_wls_satisfies_kkt_conditions():
    rng = np.random.default_rng(1)
    n, p = 300, 6
    A = np.column_stack([np.ones(n), rng.normal(size=(n, p - 1))])
    true_beta = np.array([0.5, 2.0, 0.0, -1.5, 0.0, 0.0])
    z = A @ true_beta + rng.normal(0, 0.3, n)
    w = np.ones(n)
    l1, l2 = 8.0, 0.5
    beta = models._weighted_penalized_wls(
        A, z, w, l1=l1, l2=l2, fit_intercept=True, beta_init=np.zeros(p), cd_max_iter=2000, cd_tol=1e-12
    )
    resid = z - A @ beta
    grad = -2.0 * (A.T @ (w * resid))
    assert abs(grad[0]) < 1e-4                          # intercept: unpenalised, exact stationarity
    for j in range(1, p):
        gj = grad[j] + 2 * l2 * beta[j]
        if abs(beta[j]) > 1e-6:
            assert gj == pytest.approx(-l1 * np.sign(beta[j]), abs=1e-3)   # active coefficient: exact KKT
        else:
            assert abs(gj) <= l1 + 1e-3                  # zeroed coefficient: inside the subgradient box


def test_cloglog_glm_elastic_net_shrinks_irrelevant_coefficients_more_than_ridge():
    rng = np.random.default_rng(7)
    n = 20_000
    X = rng.normal(size=(n, 5))
    true_coef = np.array([1.2, -0.8, 0.0, 0.0, 0.0])    # only the first two features matter
    eta = -3.0 + X @ true_coef
    y = rng.random(n) < models.cloglog_inverse(eta)
    ridge = models.CloglogGLM(l1=0.0, l2=1e-6).fit(X, y)
    enet = models.CloglogGLM(l1=0.15, l2=1e-6).fit(X, y)
    assert np.sum(np.abs(enet.coef_[2:])) < np.sum(np.abs(ridge.coef_[2:]))
    assert enet.coef_[0] > 0 and enet.coef_[1] < 0       # the real signals still recovered in sign


def test_monotone_hgb_respects_constraint_even_with_noise():
    rng = np.random.default_rng(1)
    n = 6000
    X = rng.normal(size=(n, 3))
    p = models.cloglog_inverse(-3 + 1.2 * X[:, 0] + rng.normal(0, 0.5, n))
    y = rng.random(n) < p
    clf = models.monotone_hgb(3, increasing=[0]).fit(X, y)
    grid = np.zeros((80, 3))
    grid[:, 0] = np.linspace(-3, 3, 80)
    assert np.all(np.diff(clf.predict_proba(grid)[:, 1]) >= -1e-9)


def test_spline_basis_partition_of_unity():
    x = np.linspace(-1.0, 1.0, 25)
    b = dlnm.spline_basis(x, knots_inner=[-0.3, 0.3], degree=3)
    assert np.allclose(b.sum(axis=1), 1.0)


def test_rare_event_weights_and_validation():
    event = np.array([True, False, False, True, False])
    w = panel.rare_event_weights(event, pi0=0.25)
    assert w.tolist() == [1.0, 4.0, 4.0, 1.0, 4.0]
    with pytest.raises(ValueError):
        panel.rare_event_weights(event, pi0=0.0)
    with pytest.raises(ValueError):
        panel.rare_event_weights(event, pi0=1.5)


def test_space_for_time_bracket_matches_eq_7_5():
    x_proj = np.array([12.0, 14.0])
    xbar_hist = 10.0
    xbar_proj = 13.0
    out = panel.space_for_time_bracket(x_proj, xbar_hist, xbar_proj)
    tilde_na, bar_na = out["no_adapt"]
    tilde_fa, bar_fa = out["full_accl"]
    assert tilde_na.tolist() == [2.0, 4.0] and bar_na == 10.0
    assert tilde_fa.tolist() == [-1.0, 1.0] and bar_fa == 13.0
    # the bracket has nonzero width whenever historical and projected place means differ
    assert not np.allclose(tilde_na, tilde_fa)


def test_correct_for_misclassification_recovers_true_hazard():
    h_true = 0.05
    se, sp = 0.9, 0.95
    h_obs = se * h_true + (1 - sp) * (1 - h_true)
    h_hat = observation.correct_for_misclassification(h_obs, se, sp)
    assert h_hat == pytest.approx(h_true, abs=1e-9)


def test_correct_for_misclassification_requires_better_than_chance():
    with pytest.raises(ValueError):
        observation.correct_for_misclassification(0.1, se=0.5, sp=0.4)


def test_correct_for_misclassification_clips_to_unit_interval():
    assert observation.correct_for_misclassification(0.0, se=0.9, sp=0.95) == pytest.approx(0.0)
    assert observation.correct_for_misclassification(1.0, se=0.9, sp=0.95) == pytest.approx(1.0)


def test_dieback_event_fires_on_persistent_decline_not_transient_dip():
    n = 20
    baseline = np.full(n, 0.6)
    persistent = baseline.copy()
    persistent[12:] = 0.05          # crashes at year 12 (0-indexed) and stays down
    nd = np.ones(n, dtype=bool)
    ev = observation.dieback_event(persistent, nd, trailing_window=10, recovery_seasons=2)
    assert ev[12]

    transient = baseline.copy()
    transient[12] = 0.05            # one bad year
    transient[13] = 0.62            # recovers immediately
    ev2 = observation.dieback_event(transient, nd, trailing_window=10, recovery_seasons=2)
    assert not ev2[12]


def test_dieback_event_suppressed_by_disturbance_flag():
    n = 20
    x = np.full(n, 0.6)
    x[12:] = 0.05
    nd_clean = np.ones(n, dtype=bool)
    nd_disturbed = np.ones(n, dtype=bool)
    nd_disturbed[12] = False        # harvest/fire flagged in the crash year
    assert observation.dieback_event(x, nd_clean, trailing_window=10)[12]
    assert not observation.dieback_event(x, nd_disturbed, trailing_window=10)[12]


def test_dieback_event_undetermined_near_record_end_is_not_flagged():
    n = 20
    x = np.full(n, 0.6)
    x[-1] = 0.05                    # crash in the very last year: no future data to confirm persistence
    nd = np.ones(n, dtype=bool)
    ev = observation.dieback_event(x, nd, trailing_window=10, recovery_seasons=2)
    assert not ev[-1]


def test_build_hazard_design_shapes_and_mundlak_columns():
    rng = np.random.default_rng(0)
    n = 60
    age = rng.uniform(0, 40, n)
    place = np.repeat(np.arange(6), 10)
    cwd = rng.normal(500, 100, n)
    X, names = pipeline.build_hazard_design(age, age_knots=[10, 20, 30], climate_covariates={"cwd": cwd}, place_id=place)
    n_age = dlnm.spline_basis(age, [10, 20, 30]).shape[1]
    assert X.shape == (n, n_age + 2)
    assert names[-2:] == ["cwd_within", "cwd_between"]
    within = X[:, names.index("cwd_within")]
    between = X[:, names.index("cwd_between")]
    assert (within + between) == pytest.approx(cwd)


def test_build_hazard_design_with_dlnm_and_stand_terms():
    rng = np.random.default_rng(1)
    n = 40
    age = rng.uniform(0, 30, n)
    place = np.repeat(np.arange(4), 10)
    cwd = rng.normal(400, 80, n)
    z_lag = rng.normal(0, 1, (n, 5))
    stand = rng.normal(0, 1, (n, 2))
    X, names = pipeline.build_hazard_design(
        age, age_knots=[10, 20], climate_covariates={"cwd": cwd}, place_id=place,
        z_lagged={"wsi": (z_lag, [-0.5, 0.5], [2.0])}, stand_terms=stand, stand_term_names=["ht_class", "density"],
    )
    assert X.shape[0] == n
    assert "ht_class" in names and "density" in names
    assert any("wsi_dlnm" in nm for nm in names)


def test_fit_stacked_hazard_recovers_signal_and_stack_weights_valid():
    rng = np.random.default_rng(3)
    n = 6000
    X = rng.normal(size=(n, 3))
    true_eta = -3.0 + 1.2 * X[:, 0] - 0.6 * X[:, 1]
    y = rng.random(n) < models.cloglog_inverse(true_eta)
    blocks = rng.integers(0, 5, n)

    out = pipeline.fit_stacked_hazard(
        X, y, blocks, glm_kwargs={"l2": 1e-4}, gbm_kwargs={"max_iter": 150},
        gbm_increasing=[0], gbm_decreasing=[1], seed=0,
    )
    w = out["stack_weights"]
    assert np.all(w >= -1e-9) and w.sum() == pytest.approx(1.0, abs=1e-6)
    valid = ~np.isnan(out["oof_pred_stat"])
    assert valid.sum() > n * 0.9
    assert out["oof_pred_stat"][valid][y[valid] == 1].mean() > out["oof_pred_stat"][valid][y[valid] == 0].mean()
    assert out["glm"].coef_[0] > 0 and out["glm"].coef_[1] < 0


def test_gate_hazard_favours_stat_inside_aoa_and_mech_outside():
    rng = np.random.default_rng(4)
    X_train = rng.normal(size=(500, 2))
    aoa = AreaOfApplicability().fit(X_train, folds=rng.integers(0, 5, 500))
    inside = X_train[:1]
    outside = np.array([[50.0, 50.0]])
    h_stat, h_mech = np.array([0.05]), np.array([0.4])
    g_in = pipeline.gate_hazard(aoa, inside, h_stat, h_mech)
    g_out = pipeline.gate_hazard(aoa, outside, h_stat, h_mech)
    assert g_in == pytest.approx(h_stat, abs=1e-6)
    assert g_out == pytest.approx(h_mech, abs=1e-3)


def test_gating_weights_and_blend():
    w = gating.aoa_weight(np.array([0.5, 1.0, 2.0, 5.0]), di_threshold=1.0, kappa=2.0)
    assert w[0] == 1.0 and w[1] == 1.0 and w[2] < 0.2 and w[3] < w[2]
    h = gating.blend_hazards(np.array([0.01]), np.array([0.2]), np.array([0.0]))
    assert h[0] == pytest.approx(0.2)
    h = gating.blend_hazards(np.array([0.01]), np.array([0.2]), np.array([1.0]))
    assert h[0] == pytest.approx(0.01)
