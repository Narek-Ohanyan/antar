"""TOPOHYDRO: the per-cell daily topoclimate forcing pipeline (concept note Sec. 5).

:func:`topoclimate_forcing` composes the tested kernels elsewhere in
``antar.climate`` into one daily time series for a single 30 m cell, in the
order the concept note specifies:

  1. downscaled T_mean/T_max/T_min (Sec. 5.1) from a pre-fitted monthly lapse
     rate and, for T_min only, a cold-air-pooling term
  2. downscaled precipitation (Sec. 5.1) from a pre-fitted monthly elevation
     gradient and exposure factor
  3. dew-point-downscaled actual vapour pressure and the three VPD
     constructions (Sec. 5.1)
  4. PET, carried as an ensemble rather than collapsed to one number
     (Sec. 5.2): ``pm_fao56``, ``priestley_taylor`` and ``energy_only`` are
     always computed (fully specified, closed-form in the inputs taken here);
     ``canopy_pm`` is computed only if the caller supplies an aerodynamic
     conductance ``g_a`` and canopy conductance ``g_c``, because the concept
     note does not specify a canopy-scale g_a formula -- flagged rather than
     approximated (see IMPLEMENTATION_LOG.md)
  5. snow partition/melt and the soil bucket, run once per PET formulation so
     the structural PET ensemble survives into CWD/WSI rather than being
     collapsed early (Sec. 5.3)
  6. soil water potential psi_s (Eq. 5.8, Clapp-Hornberger), again once per PET
     formulation -- this is TOPOHYDRO's half of the Sec. 6.2 handoff ("Phase 1
     ... is delivered by Module A"): XYLEM's stomatal-closure trigger is
     psi_s,d <= psi_close
  7. GDD, late-frost days, growing-season length and mean temperature
     (Secs. 5.3, 8.1)

Two things are deliberately upstream of this function, matching the concept
note's own module boundaries:

* Bias adjustment (QDM, ``antar.climate.bias``) is a per-month, distributional
  correction against a historical baseline; it is applied to the reference
  forcing *before* it reaches this function, not inside it.
* Absolute net radiation ``rn_mj_m2`` is taken as a direct input. If
  ``rs_flat_mj_m2`` is also supplied, a shortwave-only, geometry-only slope
  ratio (``antar.climate.radiation.slope_radiation_ratio``) is reported as a
  diagnostic, but it is not multiplied into ``rn_mj_m2``: net radiation also
  carries a longwave term that a beam-geometry ratio does not apply to, and
  the full beam/diffuse/terrain-shading/reflected decomposition of Eq. 5.4 is
  not implemented here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import downscale, indices, radiation, snow, vapour, waterbalance
from . import pet as pet_mod


@dataclass
class CellTopoclimate:
    """Daily outputs for one cell; array fields are aligned with the input day axis.

    ``pet_mm``, ``water_balance``, ``cwd_mm`` and ``wsi`` are keyed by PET
    formulation name so the structural PET ensemble is preserved end to end.
    """

    t_mean_c: np.ndarray
    t_max_c: np.ndarray
    t_min_c: np.ndarray
    p_mm: np.ndarray
    pressure_kpa: float
    rain_mm: np.ndarray
    snowfall_mm: np.ndarray
    swe_mm: np.ndarray
    liquid_input_mm: np.ndarray
    ea_kpa: np.ndarray
    vpd_day_kpa: np.ndarray
    vpd_max_kpa: np.ndarray
    vpd_24h_kpa: np.ndarray
    rs_slope_mj_m2: np.ndarray | None
    pet_mm: dict
    water_balance: dict
    cwd_mm: dict
    wsi: dict
    psi_soil_mpa: dict
    gdd_cumulative: np.ndarray
    late_frost_days: int
    growing_season_length_days: int
    growing_season_mean_t_c: float


def topoclimate_forcing(
    *,
    doy,
    month,
    t_mean_ref_c,
    t_max_ref_c,
    t_min_ref_c,
    p_ref_mm,
    ea_ref_kpa,
    u2_m_s,
    rn_mj_m2,
    z_cell_m,
    z_ref_m,
    lat_deg,
    slope_deg,
    aspect_deg,
    gamma_k_per_m,
    precip_gradient_per_m,
    w_max_mm,
    theta_sat,
    psi_sat_mpa,
    b_clapp_hornberger,
    theta_fc,
    theta_lim,
    gdd_budburst,
    concavity_index=0.0,
    calm_clear_night_frac=0.0,
    k_cap: float = 1.0,
    exposure_factor: float = 1.0,
    pressure_kpa=None,
    g_mj_m2=0.0,
    rs_flat_mj_m2=None,
    g_a=None,
    g_c=None,
    gdd_base_c: float = 5.0,
    late_frost_t_crit_c: float = -2.0,
    growing_season_t0_c: float = 0.9,
    rew_crit: float = 0.4,
) -> CellTopoclimate:
    """Run the Sec. 5 pipeline for one cell over an arbitrary-length daily record.

    ``month`` is a 1-12 integer array aligned with the day axis, used to index
    the pre-fitted ``gamma_k_per_m`` / ``precip_gradient_per_m`` (each length
    12, from :func:`antar.climate.downscale.fit_monthly_lapse_rate` /
    :func:`antar.climate.downscale.fit_precip_elevation_gradient`, fit inside
    the training fold only). ``theta_sat``/``psi_sat_mpa``/``b_clapp_hornberger``/
    ``theta_fc``/``theta_lim`` are the per-cell Clapp & Hornberger (1978) soil
    pedotransfer parameters (SoilGrids 2.0 via Saxton & Rawls 2006, Table 3);
    they drive ``psi_soil_mpa``, XYLEM's Phase-1 (stomatal-closure) trigger.
    ``gdd_budburst`` (GDD units) has no framework-wide default: it is
    species/functional-group specific and must be supplied by the caller.
    """
    doy = np.asarray(doy)
    month = np.asarray(month, dtype=int)
    t_mean_ref_c = np.asarray(t_mean_ref_c, dtype=float)
    t_max_ref_c = np.asarray(t_max_ref_c, dtype=float)
    t_min_ref_c = np.asarray(t_min_ref_c, dtype=float)
    p_ref_mm = np.asarray(p_ref_mm, dtype=float)
    ea_ref_kpa = np.asarray(ea_ref_kpa, dtype=float)
    rn = np.asarray(rn_mj_m2, dtype=float)

    if pressure_kpa is None:
        pressure_kpa = vapour.pressure_from_elevation(z_cell_m)

    # --- Sec. 5.1: topoclimate downscaling ----------------------------------
    gamma_of_day = np.asarray(gamma_k_per_m, dtype=float)[month - 1]
    delta_cap = downscale.cold_air_pooling_index(concavity_index, calm_clear_night_frac, k_cap=k_cap)
    t_mean_c = downscale.downscale_temperature(t_mean_ref_c, z_cell_m, z_ref_m, gamma_of_day)
    t_max_c = downscale.downscale_temperature(t_max_ref_c, z_cell_m, z_ref_m, gamma_of_day)
    # delta_CAP is specified as a T_min-only nocturnal cold-air-drainage correction.
    t_min_c = downscale.downscale_temperature(t_min_ref_c, z_cell_m, z_ref_m, gamma_of_day, delta_cap_k=delta_cap)

    grad_of_day = np.asarray(precip_gradient_per_m, dtype=float)[month - 1]
    p_mm = downscale.downscale_precipitation(p_ref_mm, z_cell_m, z_ref_m, grad_of_day, exposure_factor=exposure_factor)

    # --- Sec. 5.1: humidity travels as vapour pressure, VPD computed last ---
    ea_kpa = vapour.downscale_vapour_pressure_via_dewpoint(ea_ref_kpa, z_ref_m, z_cell_m)
    vpd_day = vapour.vpd_daytime(t_mean_c, t_max_c, ea_kpa)
    vpd_mx = vapour.vpd_max(t_max_c, ea_kpa)
    vpd_24 = vapour.vpd_24h(t_max_c, t_min_c, ea_kpa)

    # --- Sec. 5.2: shortwave slope-geometry diagnostic (not folded into Rn) -
    rs_slope = None
    if rs_flat_mj_m2 is not None:
        ratio = np.array([radiation.slope_radiation_ratio(lat_deg, slope_deg, aspect_deg, d) for d in doy])
        rs_slope = np.asarray(rs_flat_mj_m2, dtype=float) * ratio

    # --- Sec. 5.2: PET ensemble ----------------------------------------------
    es_kpa = vapour.saturation_vapour_pressure(t_mean_c)
    pet_mm = {
        "pm_fao56": pet_mod.pm_fao56(rn, g_mj_m2, t_mean_c, u2_m_s, es_kpa, ea_kpa, pressure_kpa),
        "priestley_taylor": pet_mod.priestley_taylor(rn, g_mj_m2, t_mean_c, pressure_kpa),
        "energy_only": pet_mod.energy_only(rn, g_mj_m2),
    }
    if g_a is not None and g_c is not None:
        # Stomatal regulation is governed by daytime VPD (design rule, Sec. 5.1).
        pet_mm["canopy_pm"] = pet_mod.canopy_pm(rn, g_mj_m2, t_mean_c, vpd_day, g_a, g_c, pressure_kpa)

    # --- Sec. 5.3: snow, then the soil bucket per PET formulation ------------
    rain_mm, snowfall_mm = snow.partition_precipitation(p_mm, t_mean_c)
    snow_out = snow.simulate_snow(p_mm, t_mean_c)

    water_balance: dict = {}
    cwd_mm: dict = {}
    wsi: dict = {}
    psi_soil_mpa: dict = {}
    for name, pet_series in pet_mm.items():
        wb = waterbalance.simulate_bucket(pet_series, snow_out["liquid_input"], w_max_mm, rew_crit=rew_crit)
        water_balance[name] = wb
        cwd_mm[name] = float(waterbalance.climatic_water_deficit(pet_series, wb["aet"]))
        wsi[name] = float(waterbalance.water_stress_integral(wb["rew"], rew_crit=rew_crit))
        theta = waterbalance.theta_from_storage(wb["w"], w_max_mm, theta_fc, theta_lim)
        psi_soil_mpa[name] = waterbalance.psi_clapp_hornberger(theta, theta_sat, psi_sat_mpa, b_clapp_hornberger)

    # --- Secs. 5.3 / 8.1: GDD, late frost, growing season --------------------
    gdd_cumulative = np.cumsum(indices.growing_degree_days(t_mean_c, base_c=gdd_base_c))
    n_late_frost = indices.late_frost_days(t_min_c, gdd_cumulative, gdd_budburst, t_crit_c=late_frost_t_crit_c)
    lgs = indices.growing_season_length(t_mean_c, t0_c=growing_season_t0_c)
    gst = indices.growing_season_mean_temperature(t_mean_c, t0_c=growing_season_t0_c)

    return CellTopoclimate(
        t_mean_c=t_mean_c,
        t_max_c=t_max_c,
        t_min_c=t_min_c,
        p_mm=p_mm,
        pressure_kpa=float(pressure_kpa) if np.ndim(pressure_kpa) == 0 else pressure_kpa,
        rain_mm=rain_mm,
        snowfall_mm=snowfall_mm,
        swe_mm=snow_out["swe"],
        liquid_input_mm=snow_out["liquid_input"],
        ea_kpa=ea_kpa,
        vpd_day_kpa=vpd_day,
        vpd_max_kpa=vpd_mx,
        vpd_24h_kpa=vpd_24,
        rs_slope_mj_m2=rs_slope,
        pet_mm=pet_mm,
        water_balance=water_balance,
        cwd_mm=cwd_mm,
        wsi=wsi,
        psi_soil_mpa=psi_soil_mpa,
        gdd_cumulative=gdd_cumulative,
        late_frost_days=n_late_frost,
        growing_season_length_days=lgs,
        growing_season_mean_t_c=gst,
    )
