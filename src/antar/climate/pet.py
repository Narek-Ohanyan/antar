"""Potential evapotranspiration - four formulations kept side by side.

PET is a *structural* uncertainty in climate-impact work: Penman-Monteith with
fixed surface resistance over-predicts future atmospheric drying (Milly & Dunne
2016; Yang et al. 2019).  v2 therefore carries an ensemble of formulations
instead of one number:

    pm_fao56            reference grass PM (FAO-56)
    canopy_pm           Monteith form with explicit canopy/aerodynamic conductance
                        (allows CO2-dependent surface resistance upstream)
    priestley_taylor    equilibrium evaporation x alpha
    energy_only         (Rn - G) / lambda        (Milly & Dunne 2016)

Units: energy MJ m-2 day-1, PET mm day-1, pressure kPa, wind m s-1.
"""
from __future__ import annotations

import numpy as np

from ..constants import CP_MJ_KG_K, LATENT_HEAT_MJ_KG, R_SPECIFIC_AIR, SECONDS_PER_DAY
from .vapour import psychrometric_constant, svp_slope


def pm_fao56(rn, g, t_mean_c, u2, es_kpa, ea_kpa, pressure_kpa):
    """FAO-56 reference evapotranspiration (mm day-1), eq. 6."""
    delta = svp_slope(t_mean_c)
    gamma = psychrometric_constant(pressure_kpa)
    num = 0.408 * delta * (rn - g) + gamma * (900.0 / (np.asarray(t_mean_c, dtype=float) + 273.0)) * u2 * (
        es_kpa - ea_kpa
    )
    den = delta + gamma * (1.0 + 0.34 * u2)
    return num / den


def air_density(t_c, pressure_kpa):
    """kg m-3, FAO-56 Box 6 (virtual temperature 1.01 T)."""
    return pressure_kpa / (1.01 * (np.asarray(t_c, dtype=float) + 273.0) * R_SPECIFIC_AIR)


def canopy_pm(rn, g, t_mean_c, vpd_kpa, g_a, g_c, pressure_kpa):
    """Penman-Monteith with explicit conductances (m s-1); returns mm day-1.

    lambda E = [ D (Rn-G) + rho cp VPD g_a ] / [ D + gamma (1 + g_a/g_c) ]
    """
    delta = svp_slope(t_mean_c)
    gamma = psychrometric_constant(pressure_kpa)
    rho = air_density(t_mean_c, pressure_kpa)
    aero = rho * CP_MJ_KG_K * vpd_kpa * g_a * SECONDS_PER_DAY
    le = (delta * (rn - g) + aero) / (delta + gamma * (1.0 + g_a / g_c))
    return le / LATENT_HEAT_MJ_KG


def stomatal_vpd_response(vpd_kpa, k_vpd: float):
    """f(VPD) in (0, 1]: multiplicative stomatal-closure response to atmospheric demand (Eq. 5.5).

    Exponential Jarvis-type form f(VPD) = exp(-k_vpd * VPD), f(0) = 1,
    monotonically decreasing. The concept note fixes the *role* of this term
    (g_c = g_s,max f(VPD) LAI) but not its functional form or coefficient;
    ``k_vpd`` (kPa-1) is therefore a required, species/scenario-calibrated
    input with no built-in default -- see IMPLEMENTATION_LOG.md.
    """
    return np.exp(-k_vpd * np.asarray(vpd_kpa, dtype=float))


def canopy_conductance(g_s_max, vpd_kpa, lai, k_vpd: float):
    """g_c = g_s,max * f(VPD) * LAI (Eq. 5.5); the CO2 scenario knob enters through g_s,max."""
    return g_s_max * stomatal_vpd_response(vpd_kpa, k_vpd) * np.asarray(lai, dtype=float)


def priestley_taylor(rn, g, t_mean_c, pressure_kpa, alpha: float = 1.26):
    delta = svp_slope(t_mean_c)
    gamma = psychrometric_constant(pressure_kpa)
    return alpha * delta / (delta + gamma) * (rn - g) / LATENT_HEAT_MJ_KG


def energy_only(rn, g):
    """Available-energy PET (Milly & Dunne 2016)."""
    return (np.asarray(rn, dtype=float) - np.asarray(g, dtype=float)) / LATENT_HEAT_MJ_KG
