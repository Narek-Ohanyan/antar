import numpy as np
import pytest

from antar.climate import downscale
from antar.niche import adult, growth, height, regeneration
from antar.viability import cohort, refugia


def test_chapman_richards_asymptote_and_inverse():
    assert growth.chapman_richards(500, 25.0, 0.04, 1.5) == pytest.approx(25.0, rel=1e-6)
    a = growth.years_to_height(5.0, 25.0, 0.04, 1.5)
    assert growth.chapman_richards(a, 25.0, 0.04, 1.5) == pytest.approx(5.0, rel=1e-9)
    assert growth.years_to_height(30.0, 25.0, 0.04) == np.inf


def test_stress_years_delay_height():
    good = growth.height_trajectory(25.0, 0.04, 1.5, np.ones(40))
    stress = np.ones(40)
    stress[10:16] = 0.2
    bad = growth.height_trajectory(25.0, 0.04, 1.5, stress)
    assert np.all(bad <= good + 1e-12) and bad[-1] < good[-1]


def test_water_modifier_half_point():
    assert growth.water_modifier(300.0, 300.0) == pytest.approx(0.5)
    assert growth.water_modifier(0.0, 300.0) == pytest.approx(1.0)


def test_gdd_modifier_shape_and_half_point():
    assert growth.gdd_modifier(0.0, 800.0) == pytest.approx(0.0)
    assert growth.gdd_modifier(800.0, 800.0) == pytest.approx(0.5)
    gdd = np.linspace(0, 4000, 50)
    f = growth.gdd_modifier(gdd, 800.0)
    assert np.all(np.diff(f) > 0) and f[-1] < 1.0


def test_late_frost_modifier_shape_and_half_point():
    assert growth.late_frost_modifier(0.0, 5.0) == pytest.approx(1.0)
    assert growth.late_frost_modifier(5.0, 5.0) == pytest.approx(0.5)
    n = np.linspace(0, 20, 50)
    f = growth.late_frost_modifier(n, 5.0)
    assert np.all(np.diff(f) < 0)


def test_thermal_treeline_modifier_is_separate_from_gdd_modifier():
    # far above both thresholds -> near 1; far below either -> near 0
    assert growth.thermal_treeline_modifier(gst=15.0, lgs=150.0, s_t=0.5, s_l=5.0) > 0.99
    assert growth.thermal_treeline_modifier(gst=0.0, lgs=150.0, s_t=0.5, s_l=5.0) < 0.01
    assert growth.thermal_treeline_modifier(gst=15.0, lgs=10.0, s_t=0.5, s_l=5.0) < 0.01
    gst = np.linspace(-2, 15, 40)
    f = growth.thermal_treeline_modifier(gst, lgs=150.0, s_t=0.5, s_l=5.0)
    assert np.all(np.diff(f) >= -1e-12)
    # a plentiful heat sum (high f_T) says nothing about whether the season is long
    # enough to grow at all (f_tl): the two must disagree in at least one regime
    assert growth.gdd_modifier(3000.0, 800.0) > 0.7
    assert growth.thermal_treeline_modifier(gst=15.0, lgs=10.0, s_t=0.5, s_l=5.0) < 0.01


def test_potential_treeline_elevation_inverts_the_lapse_relation():
    gamma = -0.0060
    z_ref, gst_ref, t_tl = 1000.0, 10.0, growth._KORNER_PAULSEN_T_TL_C
    z_tl = growth.potential_treeline_elevation(gst_ref, z_ref, gamma)
    gst_at_z_tl = downscale.downscale_temperature(gst_ref, z_tl, z_ref, gamma)
    assert gst_at_z_tl == pytest.approx(t_tl, abs=1e-9)


def test_probability_reaches_height_ensemble_fraction():
    m, t = 200, 40
    h_min = 5.0
    phi = np.ones((m, t))
    phi[m // 2:] = 0.05                             # half the ensemble is heavily stressed
    h_star, k, p = np.full(m, 25.0), np.full(m, 0.04), np.full(m, 1.5)
    a_star = growth.years_to_height(h_min, 25.0, 0.04, 1.5)
    assert a_star < t                                # sanity: an unstressed member reaches it in time
    prob = growth.probability_reaches_height(phi, h_star, k, p, h_min, horizon_years=t)
    assert prob == pytest.approx(0.5)


def test_water_limited_height_fallback_scales_by_transpiration_ratio():
    assert height.water_limited_height_fallback(20.0, aet_mm=500.0, pet_mm=500.0) == pytest.approx(20.0)
    assert height.water_limited_height_fallback(20.0, aet_mm=0.0, pet_mm=500.0) == pytest.approx(0.0)
    assert height.water_limited_height_fallback(20.0, aet_mm=600.0, pet_mm=500.0) == pytest.approx(20.0)  # clipped at 1


def test_attainable_height_is_an_upper_quantile():
    rng = np.random.default_rng(0)
    X = rng.uniform(0, 1, (4000, 2))
    h = 8 + 14 * X[:, 0] - rng.exponential(4.0, 4000)     # age/disturbance-driven downward noise
    m = height.fit_attainable_height(X, h, tau=0.9)
    cover = np.mean(h <= m.predict(X))
    assert 0.85 < cover < 0.95
    assert m.predict(X).mean() > h.mean()


def test_boyce_index_informative_vs_random():
    rng = np.random.default_rng(1)
    pa = rng.beta(2, 2, 20000)
    presences_good = rng.choice(pa, 400, p=(pa**4) / (pa**4).sum())
    presences_random = rng.choice(pa, 400)
    assert adult.boyce_index(presences_good, pa) > 0.8
    assert abs(adult.boyce_index(presences_random, pa)) < 0.6


def test_establishment_survival_product():
    h = regeneration.establishment_hazard(-2.0, [0.5], np.array([[0.0], [1.0], [2.0]]))
    assert np.all(np.diff(h) > 0)
    assert regeneration.establishment_survival(h) == pytest.approx(np.prod(1 - h))


def test_competing_risks_and_viability():
    h = np.array([[0.01, 0.02, 0.03], [0.00, 0.10, 0.00], [0.005, 0.005, 0.005]])
    tot = cohort.combine_competing_hazards(h)
    assert tot == pytest.approx(1 - np.prod(1 - h, axis=0))
    assert cohort.viability(tot, 0.9) == pytest.approx(0.9 * np.prod(1 - tot))
    assert np.all(np.diff(cohort.survival_curve(tot)) < 0)


def test_survival_exp_cumhazard_identity_cloglog():
    eta = np.log(np.array([0.01, 0.02, 0.05]))
    from antar.hazard.models import cloglog_inverse
    h = cloglog_inverse(eta)
    assert np.prod(1 - h) == pytest.approx(np.exp(-cohort.cumulative_hazard_from_cloglog(eta)))


def test_buffer_index_positive_in_cool_pocket():
    g = np.full((101, 101), 1.0)
    g[48:53, 48:53] = -1.0            # cool pocket inside a hot region
    b = refugia.buffer_index(g, radius_cells=20)
    assert b[50, 50] > 0.5 and b[10, 10] < 0.1


def test_robust_refugium_and_score():
    v = np.array([[0.9, 0.5], [0.8, 0.4], [0.85, 0.7], [0.7, 0.3], [0.95, 0.65]])
    mask = refugia.robust_refugium(v, v_star=0.6, rho=0.8)
    assert mask.tolist() == [True, False]
    s = refugia.refugium_score(v, lam=0.5)
    assert s[0] > s[1]


def test_viability_composes_with_meristem_height_probability():
    """Eq. 2.1/8.5: V = [prod (1 - h_tot)] * P[H_T >= H_min] -- the survival term from
    REFUGIUM's own cohort.viability and the height term from MERISTEM's
    growth.probability_reaches_height, composed exactly as Fig. 2 wires Module D into E.
    """
    m, t = 300, 30
    h_min = 5.0
    phi = np.full((m, t), 0.9)                      # a mildly stressed but viable ensemble
    h_star, k, p = np.full(m, 22.0), np.full(m, 0.05), np.full(m, 1.5)
    p_height_ok = growth.probability_reaches_height(phi, h_star, k, p, h_min, horizon_years=t)
    assert 0.0 < p_height_ok <= 1.0

    hazards = np.full((3, t), 0.01)                 # three low, steady cause-specific hazards
    h_tot = cohort.combine_competing_hazards(hazards)
    v = cohort.viability(h_tot, p_height_ok)
    assert v == pytest.approx(np.prod(1 - h_tot) * p_height_ok)
    assert 0.0 < v < 1.0


def test_robust_refugium_full_eq_8_6_conjunction():
    v = np.array([[0.9, 0.9], [0.8, 0.8], [0.85, 0.85], [0.7, 0.7], [0.95, 0.95]])
    assert refugia.robust_refugium(v, v_star=0.6, rho=0.8).tolist() == [True, True]     # (a) alone: both pass

    buffered = np.array([1.0, -0.5])
    assert refugia.robust_refugium(v, v_star=0.6, rho=0.8, buffer_index_grid=buffered).tolist() == [True, False]

    inside_aoa = np.array([True, False])
    assert refugia.robust_refugium(v, v_star=0.6, rho=0.8, inside_aoa=inside_aoa).tolist() == [True, False]

    # conjunctive: buffered everywhere but outside AOA at cell 1 still fails there
    combined = refugia.robust_refugium(
        v, v_star=0.6, rho=0.8, buffer_index_grid=np.array([1.0, 1.0]), inside_aoa=inside_aoa
    )
    assert combined.tolist() == [True, False]
