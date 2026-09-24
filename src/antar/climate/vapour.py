"""Vapour pressure and vapour-pressure deficit (VPD).

Design rule of v2: *humidity is downscaled as actual vapour pressure (or dew
point), never as relative humidity*, and VPD is computed last, at daily
resolution, with the temperature that actually drives the process:

* ``vpd_daytime``  - stomatal regulation (daytime temperature)
* ``vpd_max``      - upper bound of afternoon demand
* ``vpd_24h``      - cuticular / residual water loss (24 h mean saturation)
"""
from __future__ import annotations

import numpy as np

DAYTIME_T_COEF = 0.45   # Tday = Tmean + 0.45 (Tmax - Tmean); MTCLIM / Biome-BGC convention


def saturation_vapour_pressure(t_c):
    """e_s(T) in kPa for T in degrees C (FAO-56 eq. 11; Tetens form)."""
    t_c = np.asarray(t_c, dtype=float)
    return 0.6108 * np.exp(17.27 * t_c / (t_c + 237.3))


def svp_slope(t_c):
    """Slope of the saturation vapour pressure curve, kPa K-1 (FAO-56 eq. 13)."""
    t_c = np.asarray(t_c, dtype=float)
    return 4098.0 * saturation_vapour_pressure(t_c) / (t_c + 237.3) ** 2


def dewpoint_from_vapour_pressure(ea_kpa):
    """Inverse of :func:`saturation_vapour_pressure`."""
    ea = np.asarray(ea_kpa, dtype=float)
    ln = np.log(ea / 0.6108)
    return 237.3 * ln / (17.27 - ln)


def actual_vapour_pressure_from_rh(rh_pct, t_mean_c):
    """e_a from *daily-mean* RH and *daily-mean* T (they belong together)."""
    return np.asarray(rh_pct, dtype=float) / 100.0 * saturation_vapour_pressure(t_mean_c)


def daytime_temperature(t_mean_c, t_max_c, coef: float = DAYTIME_T_COEF):
    t_mean_c = np.asarray(t_mean_c, dtype=float)
    return t_mean_c + coef * (np.asarray(t_max_c, dtype=float) - t_mean_c)


def vpd_daytime(t_mean_c, t_max_c, ea_kpa):
    """Daytime VPD (kPa): e_s(T_day) - e_a, floored at 0."""
    t_day = daytime_temperature(t_mean_c, t_max_c)
    return np.maximum(saturation_vapour_pressure(t_day) - np.asarray(ea_kpa, dtype=float), 0.0)


def vpd_max(t_max_c, ea_kpa):
    """Afternoon-peak VPD (kPa): e_s(T_max) - e_a, floored at 0."""
    return np.maximum(saturation_vapour_pressure(t_max_c) - np.asarray(ea_kpa, dtype=float), 0.0)


def vpd_24h(t_max_c, t_min_c, ea_kpa):
    """24 h VPD (kPa) using the FAO-56 mean saturation vapour pressure."""
    es = 0.5 * (saturation_vapour_pressure(t_max_c) + saturation_vapour_pressure(t_min_c))
    return np.maximum(es - np.asarray(ea_kpa, dtype=float), 0.0)


def vpd_v1_formula(t_max_c, rh_mean_pct):
    """The v1 README formula, e_s(Tmax) * (1 - RH_mean/100).

    Kept **only** for regression tests that quantify its bias: it pairs a daily
    *mean* humidity with the daily *maximum* temperature and therefore
    understates afternoon VPD.
    """
    return saturation_vapour_pressure(t_max_c) * (1.0 - np.asarray(rh_mean_pct, dtype=float) / 100.0)


def downscale_vapour_pressure_via_dewpoint(ea_cell_kpa, z_cell_m, z_target_m, dewpoint_lapse_k_per_m=-0.002):
    """Elevation-adjust e_a through the dew point (not through RH).

    The dew-point lapse rate is smaller in magnitude than the dry lapse rate and
    is seasonal; the default is a placeholder that must be re-estimated per month
    from ERA5 pressure levels / station regressions (see docs).
    """
    td = dewpoint_from_vapour_pressure(ea_cell_kpa) + dewpoint_lapse_k_per_m * (
        np.asarray(z_target_m, dtype=float) - np.asarray(z_cell_m, dtype=float)
    )
    return saturation_vapour_pressure(td)


def pressure_from_elevation(z_m):
    """Atmospheric pressure (kPa), FAO-56 eq. 7."""
    z_m = np.asarray(z_m, dtype=float)
    return 101.3 * ((293.0 - 0.0065 * z_m) / 293.0) ** 5.26


def psychrometric_constant(pressure_kpa):
    return 0.000665 * np.asarray(pressure_kpa, dtype=float)
