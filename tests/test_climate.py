import numpy as np
import pytest

from antar.climate import (
    analogs,
    bias,
    downscale,
    forcing,
    indices,
    pet,
    radiation,
    snow,
    stand_coupling,
    vapour,
    waterbalance,
)


# ---------------------------------------------------------------- vapour pressure
@pytest.mark.parametrize("t,expected", [(20.0, 2.338), (25.0, 3.168), (30.0, 4.243)])
def test_saturation_vapour_pressure_matches_fao56_table(t, expected):
    assert vapour.saturation_vapour_pressure(t) == pytest.approx(expected, abs=0.003)


def test_svp_slope_fao56_value():
    assert vapour.svp_slope(20.0) == pytest.approx(0.1447, abs=5e-4)


def test_dewpoint_roundtrip():
    ea = vapour.saturation_vapour_pressure(11.3)
    assert vapour.dewpoint_from_vapour_pressure(ea) == pytest.approx(11.3, abs=1e-9)


def test_pressure_fao56_example():
    assert vapour.pressure_from_elevation(1800.0) == pytest.approx(81.8, abs=0.1)


def test_v1_vpd_formula_understates_afternoon_vpd():
    """Worked example: Tmax 30, Tmean 22, daily-mean RH 50 %. Physical e_a = RH * es(Tmean)."""
    ea = vapour.actual_vapour_pressure_from_rh(50.0, 22.0)
    correct = float(vapour.vpd_max(30.0, ea))
    v1 = float(vapour.vpd_v1_formula(30.0, 50.0))
    assert correct == pytest.approx(2.92, abs=0.02)
    assert v1 == pytest.approx(2.12, abs=0.02)
    assert (correct - v1) / correct > 0.25          # >25 % understatement


def test_vpd_ordering_and_positivity():
    ea = 1.2
    assert vapour.vpd_max(30, ea) > vapour.vpd_daytime(22, 30, ea) > 0
    assert vapour.vpd_24h(30, 14, ea) >= 0


def test_dewpoint_downscaling_cools_with_elevation():
    ea_low = 1.5
    ea_high = vapour.downscale_vapour_pressure_via_dewpoint(ea_low, 1000.0, 2000.0)
    assert ea_high < ea_low


# ---------------------------------------------------------------- PET
def test_pm_fao56_example_18_brussels():
    et0 = pet.pm_fao56(rn=13.28, g=0.0, t_mean_c=16.9, u2=2.078, es_kpa=1.997, ea_kpa=1.409, pressure_kpa=100.1)
    assert et0 == pytest.approx(3.88, abs=0.05)


def test_canopy_pm_reduces_to_fao56_for_reference_grass():
    t, u2, es, ea, P = 25.0, 2.0, 3.17, 1.6, 95.0
    rn = 15.0
    ref = pet.pm_fao56(rn, 0.0, t, u2, es, ea, P)
    can = pet.canopy_pm(rn, 0.0, t, es - ea, g_a=u2 / 208.0, g_c=1.0 / 70.0, pressure_kpa=P)
    assert can == pytest.approx(ref, rel=0.04)


def test_priestley_taylor_and_energy_only_ordering():
    rn = 14.0
    pt = pet.priestley_taylor(rn, 0.0, 20.0, 95.0)
    eo = pet.energy_only(rn, 0.0)
    assert 0 < pt < eo * 1.26


def test_stomatal_vpd_response_bounds_and_monotonicity():
    assert pet.stomatal_vpd_response(0.0, k_vpd=0.6) == pytest.approx(1.0)
    vpd = np.array([0.5, 1.0, 2.0, 4.0])
    f = pet.stomatal_vpd_response(vpd, k_vpd=0.6)
    assert np.all((f > 0) & (f <= 1.0))
    assert np.all(np.diff(f) < 0)


def test_canopy_conductance_scales_linearly_with_lai():
    g1 = pet.canopy_conductance(g_s_max=0.01, vpd_kpa=1.5, lai=2.0, k_vpd=0.6)
    g2 = pet.canopy_conductance(g_s_max=0.01, vpd_kpa=1.5, lai=4.0, k_vpd=0.6)
    assert g2 == pytest.approx(2.0 * g1)


# ---------------------------------------------------------------- snow
def test_snow_mass_conservation():
    rng = np.random.default_rng(1)
    t = rng.normal(0, 6, 400)
    p = rng.gamma(0.4, 6, 400)
    out = snow.simulate_snow(p, t)
    rain, snowfall = snow.partition_precipitation(p, t)
    assert snowfall.sum() == pytest.approx(out["melt"].sum() + out["swe"][-1], rel=1e-9)
    assert np.all(out["swe"] >= 0)


# ---------------------------------------------------------------- water balance
def test_bucket_closes_and_bounds():
    rng = np.random.default_rng(2)
    T, N = 365, 4
    pet_ = np.abs(rng.normal(3, 1.5, (T, N)))
    inp = rng.gamma(0.3, 8, (T, N)) * (rng.random((T, N)) < 0.3)
    w_max = np.array([50.0, 100.0, 150.0, 250.0])
    out = waterbalance.simulate_bucket(pet_, inp, w_max)
    assert np.max(np.abs(out["balance_error"])) < 1e-9
    assert np.all((out["rew"] >= 0) & (out["rew"] <= 1 + 1e-12))
    assert np.all(out["aet"] <= pet_ + 1e-12)
    # deeper soils buffer drought: less deficit
    cwd = waterbalance.climatic_water_deficit(pet_, out["aet"])
    assert cwd[0] >= cwd[-1]


def test_stress_integral_zero_when_wet():
    assert waterbalance.water_stress_integral(np.full(100, 0.8)) == 0.0
    assert waterbalance.water_stress_integral(np.full(10, 0.2)) == pytest.approx(5.0)


def test_clapp_hornberger_at_saturation():
    assert waterbalance.psi_clapp_hornberger(0.45, 0.45, -0.0005, 6.0) == pytest.approx(-0.0005)


# ---------------------------------------------------------------- SPEI / indices
def test_spei_is_standard_normal_on_reference():
    rng = np.random.default_rng(3)
    d = rng.normal(0, 40, 12 * 60)
    s = indices.spei_from_series(d, scale=3)
    s = s[~np.isnan(s)]
    assert abs(s.mean()) < 0.1 and abs(s.std() - 1.0) < 0.15


def test_spei_reference_fixing_detects_future_drying():
    rng = np.random.default_rng(4)
    hist = rng.normal(0, 40, 12 * 30)
    fut = rng.normal(-60, 40, 12 * 30)          # persistently drier climate
    d = np.r_[hist, fut]
    s = indices.spei_from_series(d, scale=3, ref_slice=slice(0, 12 * 30))
    assert np.nanmean(s[12 * 30:]) < -0.8 and abs(np.nanmean(s[:12 * 30])) < 0.15


def test_gdd_and_frost():
    t = np.array([0, 5, 10, 15.0])
    assert indices.growing_degree_days(t).tolist() == [0, 0, 5, 10]


def test_growing_season_length_and_mean_temperature():
    # 100 days at 5 degC (>= T0), 265 days at -5 degC (< T0)
    t = np.r_[np.full(100, 5.0), np.full(265, -5.0)]
    assert indices.growing_season_length(t, t0_c=0.9) == 100
    assert indices.growing_season_mean_temperature(t, t0_c=0.9) == pytest.approx(5.0)


def test_growing_season_length_zero_when_never_reached():
    t = np.full(365, -1.0)
    assert indices.growing_season_length(t) == 0
    assert np.isnan(indices.growing_season_mean_temperature(t))


# ---------------------------------------------------------------- topoclimate downscaling
def test_lapse_rate_recovers_known_slope():
    rng = np.random.default_rng(7)
    z = rng.uniform(600, 2600, 240)
    month = rng.integers(1, 13, 240)
    true_gamma = -0.0065 + 0.001 * np.sin(2 * np.pi * (month - 1) / 12)   # seasonally varying, K/m
    t = 15.0 + true_gamma * z
    gamma_hat = downscale.fit_monthly_lapse_rate(z, t, month)
    assert np.all(np.isfinite(gamma_hat))
    for mo in range(1, 13):
        expected = -0.0065 + 0.001 * np.sin(2 * np.pi * (mo - 1) / 12)
        assert gamma_hat[mo - 1] == pytest.approx(expected, abs=1e-6)


def test_lapse_rate_nan_when_underdetermined():
    gamma = downscale.fit_monthly_lapse_rate([500.0], [10.0], [1])
    assert np.isnan(gamma[0])


def test_downscale_temperature_matches_standard_lapse_rate():
    t_ref = 20.0
    t_cell = downscale.downscale_temperature(t_ref, z_cell_m=2000.0, z_ref_m=1000.0, gamma_k_per_m=-0.0065)
    assert t_cell == pytest.approx(20.0 - 6.5, abs=1e-9)


def test_cold_air_pooling_is_never_positive_and_zero_on_ridges():
    idx = downscale.cold_air_pooling_index(concavity_index=np.array([-1.0, 0.0, 2.0]), calm_clear_night_frac=0.8)
    assert idx[0] == 0.0 and idx[1] == 0.0
    assert idx[2] < 0.0


def test_cold_air_pooling_scales_with_calm_clear_fraction():
    still = downscale.cold_air_pooling_index(2.0, 1.0, k_cap=0.5)
    windy = downscale.cold_air_pooling_index(2.0, 0.1, k_cap=0.5)
    assert still < windy <= 0.0


def test_precip_gradient_recovers_known_slope_and_stays_nonnegative():
    rng = np.random.default_rng(8)
    z = rng.uniform(600, 2600, 240)
    month = rng.integers(1, 13, 240)
    true_grad = 0.0003
    p = 400.0 * np.exp(true_grad * z)
    grad_hat = downscale.fit_precip_elevation_gradient(z, p, month)
    assert np.all(np.isfinite(grad_hat))
    assert grad_hat == pytest.approx(np.full(12, true_grad), abs=1e-8)

    p_cell = downscale.downscale_precipitation(400.0, z_cell_m=2000.0, z_ref_m=1000.0, gradient_per_m=true_grad)
    assert p_cell == pytest.approx(400.0 * np.exp(true_grad * 1000.0), rel=1e-9)
    assert downscale.downscale_precipitation(0.0, 2000.0, 1000.0, true_grad) == 0.0


# ---------------------------------------------------------------- QDM
def test_qdm_preserves_additive_change():
    rng = np.random.default_rng(5)
    obs = rng.normal(15, 5, 3000)
    mh = rng.normal(13, 6, 3000)
    mp = rng.normal(13, 6, 3000) + 2.5
    out = bias.quantile_delta_mapping(obs, mh, mp, "additive")
    assert out.mean() - obs.mean() == pytest.approx(mp.mean() - mh.mean(), abs=0.3)


def test_qdm_preserves_multiplicative_change():
    rng = np.random.default_rng(6)
    obs = rng.gamma(2, 3, 5000)
    mh = rng.gamma(2, 2, 5000)
    mp = mh * 0.8
    out = bias.quantile_delta_mapping(obs, mh, mp, "multiplicative")
    assert out.mean() / obs.mean() == pytest.approx(0.8, rel=0.05)


# ---------------------------------------------------------------- analogues
def test_analogue_finds_nearest_and_novelty():
    means = np.array([[10.0, 500.0], [14.0, 450.0], [18.0, 300.0]])
    cov = np.diag([1.0, 40.0**2])
    j, d = analogs.best_analogue(np.array([14.2, 445.0]), means, cov)
    assert j == 1 and d < 0.5
    far = analogs.mahalanobis_to_reference(np.array([30.0, 100.0]), means, cov).min()
    assert analogs.novelty_percentile(far, k=2) > 0.999


# ---------------------------------------------------------------- stand coupling (placeholder LAI)
def test_load_placeholder_lai_reads_labelled_config(tmp_path):
    cfg = tmp_path / "stand_defaults.yaml"
    cfg.write_text("lai_by_functional_group:\n  pine: 3.0\n  oak: 3.5\nstatus: ILLUSTRATIVE_PLACEHOLDER\n")
    lai = stand_coupling.load_placeholder_lai(cfg)
    assert lai == {"pine": 3.0, "oak": 3.5}


def test_load_placeholder_lai_rejects_unlabelled_config(tmp_path):
    cfg = tmp_path / "stand_defaults.yaml"
    cfg.write_text("lai_by_functional_group:\n  pine: 3.0\nstatus: CALIBRATED\n")
    with pytest.raises(ValueError):
        stand_coupling.load_placeholder_lai(cfg)


# ---------------------------------------------------------------- radiation geometry
def test_flat_surface_ratio_is_one():
    assert radiation.slope_radiation_ratio(40.0, 0.0, 0.0, doy=355) == pytest.approx(1.0, abs=1e-6)


def test_south_facing_slope_gains_in_winter_north_loses():
    south = radiation.slope_radiation_ratio(40.0, 30.0, 180.0, doy=355)
    north = radiation.slope_radiation_ratio(40.0, 30.0, 0.0, doy=355)
    assert south > 1.3 and north < 0.3


def test_svf_bounds():
    assert radiation.sky_view_factor_from_horizon(np.zeros(16)) == pytest.approx(1.0)
    assert radiation.sky_view_factor_from_horizon(np.full(16, np.pi / 2)) == pytest.approx(0.0, abs=1e-12)
    assert radiation.sky_view_factor_planar_slope(0.0) == 1.0


def test_extraterrestrial_radiation_fao56_example():
    # FAO-56 Example 8: latitude 20 S (-0.35 rad), 3 September (doy 246) -> Ra = 32.2 MJ m-2 d-1
    assert radiation.extraterrestrial_radiation(-0.35, 246) == pytest.approx(32.2, abs=0.2)


# ---------------------------------------------------------------- forcing pipeline (TOPOHYDRO orchestration)
def _synthetic_year(seed=0):
    rng = np.random.default_rng(seed)
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month = np.repeat(np.arange(1, 13), days_in_month)
    doy = np.arange(1, 366)
    t_mean = 10.0 + 12.0 * np.sin(2 * np.pi * (doy - 80) / 365.0)
    t_max = t_mean + 8.0 + rng.normal(0, 0.5, 365)
    t_min = t_mean - 6.0 + rng.normal(0, 0.5, 365)
    p_ref = rng.gamma(0.4, 6.0, 365) * (rng.random(365) < 0.3)
    ea_ref = vapour.actual_vapour_pressure_from_rh(60.0, t_mean)
    u2 = np.full(365, 2.0) + rng.normal(0, 0.2, 365)
    rn = 5.0 + 10.0 * np.sin(2 * np.pi * (doy - 80) / 365.0)
    return dict(doy=doy, month=month, t_mean=t_mean, t_max=t_max, t_min=t_min, p_ref=p_ref, ea_ref=ea_ref, u2=u2, rn=rn)


def _call_forcing(y, **overrides):
    kwargs = dict(
        doy=y["doy"], month=y["month"],
        t_mean_ref_c=y["t_mean"], t_max_ref_c=y["t_max"], t_min_ref_c=y["t_min"],
        p_ref_mm=y["p_ref"], ea_ref_kpa=y["ea_ref"], u2_m_s=y["u2"], rn_mj_m2=y["rn"],
        z_cell_m=1500.0, z_ref_m=1000.0, lat_deg=40.2, slope_deg=15.0, aspect_deg=180.0,
        gamma_k_per_m=np.full(12, -0.0060), precip_gradient_per_m=np.full(12, 0.0002),
        w_max_mm=120.0,
        theta_sat=0.45, psi_sat_mpa=-0.0005, b_clapp_hornberger=6.0, theta_fc=0.30, theta_lim=0.10,
        gdd_budburst=100.0,
        concavity_index=1.5, calm_clear_night_frac=0.3, k_cap=1.0,
    )
    kwargs.update(overrides)
    return forcing.topoclimate_forcing(**kwargs)


def test_forcing_pipeline_shapes_and_pet_ensemble_without_canopy():
    y = _synthetic_year()
    out = _call_forcing(y)
    assert out.t_mean_c.shape == (365,)
    assert set(out.pet_mm) == {"pm_fao56", "priestley_taylor", "energy_only"}
    assert set(out.water_balance) == set(out.pet_mm)
    assert out.rs_slope_mj_m2 is None
    assert out.pressure_kpa == pytest.approx(vapour.pressure_from_elevation(1500.0))


def test_forcing_pipeline_adds_canopy_pm_only_when_conductances_given():
    y = _synthetic_year()
    g_c = pet.canopy_conductance(g_s_max=0.01, vpd_kpa=1.2, lai=3.0, k_vpd=0.6)
    out = _call_forcing(y, g_a=0.01, g_c=g_c)
    assert "canopy_pm" in out.pet_mm and "canopy_pm" in out.water_balance


def test_forcing_pipeline_water_balance_closes_for_every_pet_formulation():
    y = _synthetic_year()
    out = _call_forcing(y)
    for name, wb in out.water_balance.items():
        assert np.max(np.abs(wb["balance_error"])) < 1e-6, name
        assert np.all((wb["rew"] >= -1e-12) & (wb["rew"] <= 1 + 1e-12)), name


def test_forcing_pipeline_soil_water_potential_matches_water_balance_keys():
    y = _synthetic_year()
    out = _call_forcing(y)
    assert set(out.psi_soil_mpa) == set(out.water_balance)
    for name, psi in out.psi_soil_mpa.items():
        assert psi.shape == (365,)
        assert np.all(psi <= 0.0), name          # water potential is never positive
    # a wetter bucket (larger w_max, same inputs) should sit closer to saturation (less negative psi)
    wet = _call_forcing(y, w_max_mm=400.0)
    assert np.mean(wet.psi_soil_mpa["pm_fao56"]) >= np.mean(out.psi_soil_mpa["pm_fao56"])


def test_forcing_pipeline_vpd_ordering_and_gdd_monotonicity():
    y = _synthetic_year()
    out = _call_forcing(y)
    assert np.all(out.vpd_max_kpa + 1e-9 >= out.vpd_day_kpa)
    assert np.all(np.diff(out.gdd_cumulative) >= -1e-9)
    assert 0 < out.growing_season_length_days < 365


def test_forcing_pipeline_reports_slope_radiation_only_when_requested():
    y = _synthetic_year()
    rs_flat = 10.0 + 10.0 * np.sin(2 * np.pi * (y["doy"] - 80) / 365.0)
    out = _call_forcing(y, rs_flat_mj_m2=rs_flat)
    assert out.rs_slope_mj_m2 is not None and out.rs_slope_mj_m2.shape == (365,)


def test_forcing_pipeline_requires_gdd_budburst():
    y = _synthetic_year()
    with pytest.raises(TypeError):
        forcing.topoclimate_forcing(
            doy=y["doy"], month=y["month"],
            t_mean_ref_c=y["t_mean"], t_max_ref_c=y["t_max"], t_min_ref_c=y["t_min"],
            p_ref_mm=y["p_ref"], ea_ref_kpa=y["ea_ref"], u2_m_s=y["u2"], rn_mj_m2=y["rn"],
            z_cell_m=1500.0, z_ref_m=1000.0, lat_deg=40.2, slope_deg=15.0, aspect_deg=180.0,
            gamma_k_per_m=np.full(12, -0.006), precip_gradient_per_m=np.full(12, 0.0002),
            w_max_mm=120.0,
            theta_sat=0.45, psi_sat_mpa=-0.0005, b_clapp_hornberger=6.0, theta_fc=0.30, theta_lim=0.10,
        )
