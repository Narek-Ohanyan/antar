"""Saxton & Rawls (2006) pedotransfer functions: soil texture (sand/clay/organic
matter fractions, e.g. from SoilGrids) -> the hydraulic parameters
``antar.climate.waterbalance.theta_from_storage``/``psi_clapp_hornberger`` and
``antar.climate.forcing.topoclimate_forcing`` need (theta_sat, theta_fc,
theta_lim, psi_sat_mpa, b_clapp_hornberger) -- referenced by
``antar.io.gee_export.export_soils`` and this project's README since TOPOHYDRO
was first built, but never implemented: ``waterbalance.py`` only ever consumed
these parameters, nothing derived them from soil texture.

Equations and every coefficient below are transcribed from Saxton and Rawls's
own Table 1 (as reproduced, with source attribution, in Aliku & Oshunsaya,
2016, "Establishing relationship between measured and predicted soil water
characteristics using SOILWAT model," Geosci. Model Dev. Discuss.,
doi:10.5194/gmd-2016-165 -- fetched and read directly rather than recalled
from memory, given how consequential a wrong pedotransfer coefficient would
be; the original Saxton & Rawls (2006) PDF was unreachable (server down) when
this was written).

Units, per the paper's own Table 2 (not assumed): sand and clay are DECIMAL
FRACTIONS BY WEIGHT (0-1, not 0-100, and not by volume); organic matter is a
DECIMAL FRACTION BY VOLUME (0-1). SoilGrids reports sand/clay as g/kg (0-1000)
and SOC as g/kg -- converted here, not left for the caller to get wrong.

Citation: Saxton, K.E., Rawls, W.J. (2006). Soil water characteristic
estimates by texture and organic matter for hydrologic solutions. Soil
Science Society of America Journal, 70, 1569-1578.
"""
from __future__ import annotations

import numpy as np


def _theta_1500(sand_frac, clay_frac, om_frac):
    """Wilting point moisture (1500 kPa), m3/m3 -- Saxton & Rawls Eq. 1."""
    s, c, om = sand_frac, clay_frac, om_frac
    theta_1500t = (-0.024 * s + 0.487 * c + 0.006 * om + 0.005 * (s * om)
                   - 0.013 * (c * om) + 0.068 * (s * c) + 0.031)
    return theta_1500t + (0.14 * theta_1500t - 0.02)


def _theta_33(sand_frac, clay_frac, om_frac):
    """Field capacity moisture (33 kPa), m3/m3 -- Saxton & Rawls Eq. 2."""
    s, c, om = sand_frac, clay_frac, om_frac
    theta_33t = (-0.251 * s + 0.195 * c + 0.011 * om + 0.006 * (s * om)
                 - 0.027 * (c * om) + 0.452 * (s * c) + 0.299)
    return theta_33t + (1.283 * theta_33t ** 2 - 0.374 * theta_33t - 0.015)


def _theta_s_minus_33(sand_frac, clay_frac, om_frac):
    """SAT-33kPa moisture, m3/m3 -- Saxton & Rawls Eq. 3."""
    s, c, om = sand_frac, clay_frac, om_frac
    t = (0.278 * s + 0.034 * c + 0.022 * om - 0.018 * (s * om)
         - 0.027 * (c * om) - 0.584 * (s * c) + 0.078)
    return t + (0.636 * t - 0.107)


def soil_hydraulic_parameters(sand_pct, clay_pct, soc_g_kg):
    """SoilGrids clay/sand/SOC (its own native units) -> the parameters
    ``topoclimate_forcing`` needs. Returns a dict with theta_sat, theta_fc,
    theta_lim (m3/m3), psi_sat_mpa (negative, MPa) and b_clapp_hornberger
    (dimensionless), each broadcastable over the input arrays' shape.

    ``sand_pct``/``clay_pct``: percent by weight, 0-100 (SoilGrids' own
    convention, e.g. from ``export_soils``'s clay/sand bands -- already percent,
    not the g/kg some other SoilGrids layers use; confirm against the actual
    band units before reusing this for a different soil source).
    ``soc_g_kg``: soil organic carbon, g/kg (SoilGrids' native SOC unit).
    Converted to the organic-matter volume fraction Saxton & Rawls's equations
    need via the standard van Bemmelen factor (OM% = SOC% x 1.724) and an
    assumed soil bulk density of 1.3 Mg/m3 (a common default for mineral
    topsoil where no measured bulk density is available -- not a measured
    value, flagged here rather than silently treated as one).
    """
    sand_pct = np.asarray(sand_pct, dtype=float)
    clay_pct = np.asarray(clay_pct, dtype=float)
    soc_g_kg = np.asarray(soc_g_kg, dtype=float)

    sand_frac = sand_pct / 100.0
    clay_frac = clay_pct / 100.0
    # SOC (g/kg dry soil) -> OM mass fraction (x1.724, van Bemmelen) -> OM
    # volume fraction (x assumed bulk density / water density) -- see docstring.
    assumed_bulk_density_mg_m3 = 1.3
    om_mass_frac = (soc_g_kg / 1000.0) * 1.724
    om_frac = om_mass_frac * assumed_bulk_density_mg_m3  # water density = 1 Mg/m3

    theta_1500 = _theta_1500(sand_frac, clay_frac, om_frac)
    theta_33 = _theta_33(sand_frac, clay_frac, om_frac)
    theta_s_minus_33 = _theta_s_minus_33(sand_frac, clay_frac, om_frac)
    theta_s = theta_33 + theta_s_minus_33 - 0.097 * sand_frac + 0.043  # Eq. 5

    # Eqs. 14-15: the same power-law form as Clapp & Hornberger (1978),
    # psi = A * theta^-B, fit through the paper's own 33kPa/1500kPa points --
    # A/B reused directly as antar's psi_sat_mpa/b_clapp_hornberger (see module
    # docstring) rather than introducing a second, redundant power-law fit.
    b_clapp_hornberger = (np.log(1500.0) - np.log(33.0)) / (np.log(theta_33) - np.log(theta_1500))
    a_kpa = np.exp(np.log(33.0) + b_clapp_hornberger * np.log(theta_33))
    # A*theta^-B is a positive tension (kPa) at theta=theta_s; antar's convention
    # is a negative matric potential in MPa (see waterbalance.psi_clapp_hornberger).
    psi_sat_mpa = -(a_kpa * theta_s ** (-b_clapp_hornberger)) / 1000.0

    return {
        "theta_sat": theta_s,
        "theta_fc": theta_33,      # field capacity, by definition theta_33 (33 kPa)
        "theta_lim": theta_1500,   # wilting point, by definition theta_1500 (1500 kPa)
        "psi_sat_mpa": psi_sat_mpa,
        "b_clapp_hornberger": b_clapp_hornberger,
    }
