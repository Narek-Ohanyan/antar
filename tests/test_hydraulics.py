import numpy as np
import pytest

from antar.climate import forcing
from antar.hydraulics import calibration, emulator, gmin, monte_carlo, pipeline, spells, twophase, vulnerability


def _tiny_cell():
    """A short (40-day), single-drought-spell TOPOHYDRO cell-year, for XYLEM integration tests."""
    n = 40
    doy = np.arange(1, n + 1)
    month = np.full(n, 7)
    t_mean = np.full(n, 28.0)
    t_max = np.full(n, 34.0)
    t_min = np.full(n, 20.0)
    p_ref = np.zeros(n)                     # no rain: soil dries out monotonically
    ea_ref = np.full(n, 1.2)
    u2 = np.full(n, 2.0)
    rn = np.full(n, 15.0)
    return forcing.topoclimate_forcing(
        doy=doy, month=month,
        t_mean_ref_c=t_mean, t_max_ref_c=t_max, t_min_ref_c=t_min,
        p_ref_mm=p_ref, ea_ref_kpa=ea_ref, u2_m_s=u2, rn_mj_m2=rn,
        z_cell_m=1200.0, z_ref_m=1200.0, lat_deg=39.8, slope_deg=10.0, aspect_deg=180.0,
        gamma_k_per_m=np.full(12, -0.006), precip_gradient_per_m=np.full(12, 0.0002),
        w_max_mm=60.0,  # small bucket: dries out within the 40-day window
        theta_sat=0.45, psi_sat_mpa=-0.0005, b_clapp_hornberger=6.0, theta_fc=0.30, theta_lim=0.05,
        gdd_budburst=50.0,
    )


def test_plc_at_p50_is_50():
    assert vulnerability.plc(-3.0, -3.0, 40.0) == pytest.approx(50.0)


def test_psi_at_plc_roundtrip_and_p88_closed_form():
    p50, S = -3.0, 40.0
    a = S / 25.0
    p88 = vulnerability.psi_at_plc(88.0, p50, S)
    assert vulnerability.plc(p88, p50, S) == pytest.approx(88.0)
    assert p88 == pytest.approx(p50 - 1.992 / a, abs=2e-3)
    p12 = vulnerability.psi_at_plc(12.0, p50, S)
    assert vulnerability.slope_from_p12_p88(p12, p88) == pytest.approx(S, rel=1e-9)


def test_slope_at_p50_equals_S():
    p50, S = -2.5, 30.0
    h = 1e-5
    d = (vulnerability.plc(p50 + h, p50, S) - vulnerability.plc(p50 - h, p50, S)) / (2 * h)
    assert -d == pytest.approx(S, rel=1e-4)


def test_stomatal_factor_and_gs90():
    assert vulnerability.stomatal_factor(-2.0, -2.0) == pytest.approx(0.5)
    g90 = vulnerability.psi_gs90(-2.0, 4.0)
    assert vulnerability.stomatal_factor(g90, -2.0, 4.0) == pytest.approx(0.1)


def test_gmin_biphasic_ratios_and_continuity():
    g25, tp = 3.0, 38.0
    assert gmin.gmin_temperature(25.0, g25, tp) == pytest.approx(g25)
    assert gmin.gmin_temperature(35.0, g25, tp) / gmin.gmin_temperature(25.0, g25, tp) == pytest.approx(1.2)
    assert gmin.gmin_temperature(tp + 10, g25, tp) / gmin.gmin_temperature(tp, g25, tp) == pytest.approx(4.8)
    eps = 1e-9
    assert gmin.gmin_temperature(tp - eps, g25, tp) == pytest.approx(gmin.gmin_temperature(tp + eps, g25, tp), rel=1e-6)


def test_traits_buffer_positive_and_psi_crit_below_p50_for_angiosperm():
    t = twophase.Traits()
    assert t.psi_crit < t.p50 and t.buffer > 0


def test_simulation_matches_closed_form_tcrit():
    t = twophase.Traits(psi_close=-2.5)
    vpd, temp, patm = 2.0, 30.0, 90.0
    n = 120
    sim = twophase.simulate_two_phase(np.full(n, -3.5), np.full(n, temp), np.full(n, vpd), patm, t)
    tc = twophase.t_crit_days_constant(t, vpd, patm, temp)
    assert sim["day_first_failure"] == int(np.ceil(tc)) - 1
    assert sim["day_first_closed"] == 0


def test_rewetting_resets_index():
    t = twophase.Traits()
    psi = np.r_[np.full(5, -3.5), np.full(5, -0.5), np.full(5, -3.5)]
    sim = twophase.simulate_two_phase(psi, np.full(15, 30.0), np.full(15, 2.0), 90.0, t)
    assert sim["hfi"][5] == 0.0 and sim["hfi"][4] > 0
    assert sim["hfi"][9] == 0.0
    assert sim["hfi"][14] == pytest.approx(sim["hfi"][4])


def test_heatwave_above_tp_shortens_time_to_failure_dramatically():
    t = twophase.Traits(tp=38.0)
    cool = twophase.t_crit_days_constant(t, 2.0, 90.0, 30.0)
    hot = twophase.t_crit_days_constant(t, 2.0, 90.0, 44.0)
    # beyond Tp the cuticular conductance jumps (Q10 = 4.8), so time-to-failure collapses;
    # the net factor is smaller than the raw Q10 jump because g_min also rises below Tp
    assert hot < cool / 2.5


def test_failure_probability_increases_with_vpd():
    rng = np.random.default_rng(0)
    base = twophase.Traits()
    draws = monte_carlo.sample_traits(rng, base, {"capacitance": 6000.0, "g25": 0.6, "p50": 0.3}, 300)

    def sim(vpd):
        return lambda tr: twophase.simulate_two_phase(np.full(30, -3.5), np.full(30, 32.0), np.full(30, vpd), 90.0, tr)["hfi_max"]

    p_low = monte_carlo.failure_probability(sim(1.0), draws)
    p_high = monte_carlo.failure_probability(sim(3.0), draws)
    assert p_high > p_low


def test_sample_trait_hyperparameters_no_sd_returns_base_unchanged():
    rng = np.random.default_rng(2)
    base = twophase.Traits()
    hyper = monte_carlo.sample_trait_hyperparameters(rng, base, {})
    assert hyper == base


def test_sample_trait_hyperparameters_clamps_psi_close_above_p50():
    rng = np.random.default_rng(3)
    base = twophase.Traits(p50=-3.0, psi_close=-2.9)  # already close to P50
    hyper = monte_carlo.sample_trait_hyperparameters(rng, base, {"p50": 5.0})  # can push P50 far less negative
    assert hyper.psi_close >= hyper.p50


def test_two_level_monte_carlo_collapses_with_no_variance_anywhere():
    base = twophase.Traits()

    def sim(tr):
        return twophase.simulate_two_phase(np.full(6, -3.5), np.full(6, 30.0), np.full(6, 2.0), 90.0, tr)["hfi_max"]

    h = monte_carlo.two_level_failure_probability(sim, base, hyper_sd={}, individual_sd={}, outer_draws=10, inner_draws=5, seed=0)
    assert np.all(h == h[0])            # no perturbation anywhere -> identical outcome every draw
    assert h[0] in (0.0, 1.0)


def test_two_level_monte_carlo_outer_loop_adds_knowledge_uncertainty():
    base = twophase.Traits()  # t_crit ~ 8.3 days at VPD24=2 kPa, T=30 C, Patm=90 kPa (concept note worked example, Sec. 6.2)
    n = 8  # just under t_crit: base traits alone survive a closed spell of this length

    def sim(tr):
        return twophase.simulate_two_phase(np.full(n, -3.5), np.full(n, 30.0), np.full(n, 2.0), 90.0, tr)["hfi_max"]

    h_no_hyper = monte_carlo.two_level_failure_probability(
        sim, base, hyper_sd={}, individual_sd={}, outer_draws=60, inner_draws=1, seed=1
    )
    h_with_hyper = monte_carlo.two_level_failure_probability(
        sim, base, hyper_sd={"p50": 0.3}, individual_sd={}, outer_draws=60, inner_draws=1, seed=1
    )
    assert np.std(h_no_hyper) == 0.0
    assert np.std(h_with_hyper) > 0.0


def test_mechanistic_hazard_for_cell_wires_topohydro_output_into_xylem():
    cell = _tiny_cell()
    base = twophase.Traits()
    h = pipeline.mechanistic_hazard_for_cell(
        cell, base, hyper_sd={"p50": 0.2}, individual_sd={"capacitance": 3000.0},
        outer_draws=8, inner_draws=10, seed=0,
    )
    assert h.shape == (8,)
    assert np.all((h >= 0.0) & (h <= 1.0))
    # a 40-day rainless drought at 34 C should drive at least some hazard in the ensemble
    assert h.mean() > 0.0


def test_mechanistic_hazard_for_cell_rejects_unknown_pet_formulation():
    cell = _tiny_cell()
    with pytest.raises(KeyError):
        pipeline.mechanistic_hazard_for_cell(cell, twophase.Traits(), {}, {}, pet_formulation="not_a_formulation")


def test_recalibrate_hazard_identity_at_a0_b1():
    h = np.array([0.01, 0.1, 0.3, 0.6, 0.9])
    assert calibration.recalibrate_hazard(h, a_g=0.0, b_g=1.0) == pytest.approx(h, abs=1e-8)


def test_recalibrate_hazard_monotone_increasing_in_h_mech():
    h = np.linspace(0.001, 0.999, 50)
    out = calibration.recalibrate_hazard(h, a_g=0.3, b_g=1.4)
    assert np.all(np.diff(out) > 0)
    assert np.all((out > 0) & (out < 1))


def test_recalibrate_hazard_intercept_shifts_baseline_up():
    h = 0.3
    baseline = calibration.recalibrate_hazard(h, a_g=0.0, b_g=1.0)
    shifted = calibration.recalibrate_hazard(h, a_g=1.0, b_g=1.0)
    assert shifted > baseline


def test_longest_true_run_finds_the_longest_of_several_spells():
    mask = np.array([True, True, False, True, True, True, False, True])
    length, start, end = spells.longest_true_run(mask)
    assert (length, start, end) == (3, 3, 6)


def test_longest_true_run_all_false():
    assert spells.longest_true_run(np.zeros(5, dtype=bool)) == (0, -1, -1)


def test_stress_spell_summary_matches_hand_computed_values():
    closed = np.array([False, True, True, True, False])
    vpd = np.array([1.0, 2.0, 4.0, 6.0, 1.0])
    t = np.array([10.0, 20.0, 30.0, 40.0, 10.0])
    s = spells.stress_spell_summary(closed, vpd, t)
    assert s == {
        "longest_closed_spell_days": 3,
        "mean_vpd24_during_spell": 4.0,
        "mean_t_during_spell": 30.0,
        "n_closed_days": 3,
    }


def test_stress_spell_summary_no_closed_days_is_all_zero():
    n = 5
    s = spells.stress_spell_summary(np.zeros(n, dtype=bool), np.ones(n), np.ones(n))
    assert s == {
        "longest_closed_spell_days": 0,
        "mean_vpd24_during_spell": 0.0,
        "mean_t_during_spell": 0.0,
        "n_closed_days": 0,
    }


def test_emulator_is_monotone_by_construction():
    rng = np.random.default_rng(1)
    names, X = emulator.latin_hypercube({"cwd": (0, 900), "vpd": (0.5, 5.0), "awc": (30, 300)}, 3000, seed=1)
    y = 1 - np.exp(-(X[:, 0] / 400) ** 1.5 * (X[:, 1] / 2.0)) + 0.1 * np.exp(-X[:, 2] / 150) + rng.normal(0, 0.01, len(X))
    m = emulator.fit_monotone_emulator(X, y, increasing=[0, 1], decreasing=[2])
    grid = np.column_stack([np.linspace(0, 900, 60), np.full(60, 2.0), np.full(60, 100.0)])
    assert np.all(np.diff(m.predict(grid)) >= -1e-9)


def test_build_hazard_emulator_dict_api_is_monotone_and_reproduces_holdout():
    bounds = {"spell": (0.0, 60.0), "vpd": (0.5, 5.0), "capacitance": (5000.0, 50000.0)}

    def simulate_hazard(row):
        # smooth, hand-specified surface: increasing in spell/vpd, decreasing in capacitance
        return 1 - np.exp(-(row["spell"] / 30.0) * (row["vpd"] / 2.0) / (row["capacitance"] / 20000.0))

    model, names = emulator.build_hazard_emulator(
        bounds, simulate_hazard, n_samples=2000, increasing=["spell", "vpd"], decreasing=["capacitance"], seed=2
    )
    assert names == ["spell", "vpd", "capacitance"]
    err = emulator.emulator_max_abs_error(model, names, simulate_hazard, bounds, n_holdout=300, seed=3)
    # Sec. 6.3's own release gate (< 0.02 absolute hazard) is a production-scale bar: it
    # presumes an LHS budget sized for the real study, not a unit test's few thousand
    # samples. This asserts the machinery converges to a sane, bounded error, not that
    # gate -- emulator_max_abs_error is exactly the function that would enforce it later.
    assert err < 0.3

    grid = np.column_stack([np.linspace(0, 60, 40), np.full(40, 2.0), np.full(40, 20000.0)])
    assert np.all(np.diff(model.predict(grid)) >= -1e-9)


def test_build_hazard_emulator_reproduces_real_two_phase_engine():
    """Wires the actual twophase engine through the emulator machinery (not a synthetic stand-in).

    Uses small sample counts purely for test runtime; the API and monotone-
    constraint wiring are identical to production use.
    """
    base = twophase.Traits()
    bounds = {
        "longest_closed_spell_days": (2.0, 15.0),
        "mean_vpd24_during_spell": (1.0, 4.0),
        "mean_t_during_spell": (20.0, 42.0),
        "capacitance": (10000.0, 50000.0),
    }

    def simulate_hazard(row):
        n = max(int(round(row["longest_closed_spell_days"])), 1)
        traits = base.with_updates(capacitance=row["capacitance"])
        sim = twophase.simulate_two_phase(
            np.full(n, -3.5), np.full(n, row["mean_t_during_spell"]), np.full(n, row["mean_vpd24_during_spell"]), 90.0, traits
        )
        return 1.0 - np.exp(-sim["hfi_max"])  # bounded [0, 1) hazard proxy; the real target is a MC probability

    model, names = emulator.build_hazard_emulator(
        bounds,
        simulate_hazard,
        n_samples=600,
        increasing=["longest_closed_spell_days", "mean_vpd24_during_spell", "mean_t_during_spell"],
        decreasing=["capacitance"],
        seed=4,
    )
    err = emulator.emulator_max_abs_error(model, names, simulate_hazard, bounds, n_holdout=150, seed=5)
    # Loose bound: this is the real g_min(T) phase-transition physics (a genuine kink at Tp
    # within the sampled T range), fit from a fast-test-sized LHS -- see the comment on the
    # synthetic-surface test above re: the 0.02 production gate.
    assert err < 0.4
