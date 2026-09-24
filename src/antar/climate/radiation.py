"""Terrain-corrected radiation geometry.

Only *geometry* is implemented here (slope/aspect illumination ratio, sky-view
factor).  Absolute radiation comes from ERA5-Land / CHELSA ``rsds`` and is
partitioned into beam and diffuse parts upstream.
"""
from __future__ import annotations

import numpy as np


def solar_declination(doy):
    return 0.409 * np.sin(2.0 * np.pi / 365.0 * np.asarray(doy, dtype=float) - 1.39)


def inverse_relative_distance(doy):
    return 1.0 + 0.033 * np.cos(2.0 * np.pi / 365.0 * np.asarray(doy, dtype=float))


def sunset_hour_angle(lat_rad, doy):
    dec = solar_declination(doy)
    return np.arccos(np.clip(-np.tan(lat_rad) * np.tan(dec), -1.0, 1.0))


def extraterrestrial_radiation(lat_rad, doy):
    """Ra, MJ m-2 day-1 on a horizontal surface (FAO-56 eq. 21)."""
    dr = inverse_relative_distance(doy)
    dec = solar_declination(doy)
    ws = sunset_hour_angle(lat_rad, doy)
    gsc = 0.0820
    return (24 * 60 / np.pi) * gsc * dr * (
        ws * np.sin(lat_rad) * np.sin(dec) + np.cos(lat_rad) * np.cos(dec) * np.sin(ws)
    )


def cos_incidence(slope_rad, aspect_rad, zenith_rad, azimuth_rad):
    """cos(theta_i) between the sun and the slope normal.

    Aspect and sun azimuth are both measured clockwise from north.
    cos(theta) = cos(b) cos(Z) + sin(b) sin(Z) cos(phi_s - phi_a)
    """
    return np.cos(slope_rad) * np.cos(zenith_rad) + np.sin(slope_rad) * np.sin(zenith_rad) * np.cos(
        azimuth_rad - aspect_rad
    )


def _sun_position(lat_rad, dec, hour_angle):
    """Zenith and azimuth (clockwise from north) for hour angle(s) (rad, 0 = solar noon)."""
    cos_z = np.sin(lat_rad) * np.sin(dec) + np.cos(lat_rad) * np.cos(dec) * np.cos(hour_angle)
    cos_z = np.clip(cos_z, -1.0, 1.0)
    zen = np.arccos(cos_z)
    sin_z = np.sqrt(np.maximum(1.0 - cos_z**2, 1e-12))
    # azimuth measured from north, clockwise
    cos_az = (np.sin(dec) - cos_z * np.sin(lat_rad)) / (sin_z * np.cos(lat_rad))
    az = np.arccos(np.clip(cos_az, -1.0, 1.0))
    az = np.where(hour_angle > 0, 2.0 * np.pi - az, az)  # afternoon: west of south
    return zen, az


def slope_radiation_ratio(lat_deg, slope_deg, aspect_deg, doy, n_steps: int = 288, tau: float = 0.75):
    """Ratio of clear-sky *direct-beam* radiation on a slope to that on flat ground (daily integral).

    Simple Bouguer-Lambert attenuation ``tau ** air_mass`` is used so the ratio
    reflects both illumination angle and path length.  No cast shadows here
    (those come from a horizon model in the production code).
    """
    lat = np.deg2rad(lat_deg)
    slope = np.deg2rad(slope_deg)
    aspect = np.deg2rad(aspect_deg)
    dec = float(solar_declination(doy))
    ws = float(sunset_hour_angle(lat, doy))
    if ws <= 1e-6:
        return np.nan
    h = np.linspace(-ws, ws, n_steps)
    zen, az = _sun_position(lat, dec, h)
    cz = np.cos(zen)
    up = cz > 0.02
    air_mass = 1.0 / np.maximum(cz, 0.02)
    beam = tau**air_mass
    flat = np.where(up, beam * cz, 0.0).sum()
    ci = cos_incidence(slope, aspect, zen, az)
    sloped = np.where(up, beam * np.maximum(ci, 0.0), 0.0).sum()
    return float(sloped / flat) if flat > 0 else np.nan


def sky_view_factor_from_horizon(horizon_elev_rad):
    """Sky-view factor of a horizontal surface, (1/N) sum cos^2(H_k) (Dozier & Frew 1990)."""
    h = np.asarray(horizon_elev_rad, dtype=float)
    return float(np.mean(np.cos(h) ** 2))


def sky_view_factor_planar_slope(slope_rad):
    """Isotropic sky-view factor of an unobstructed planar slope, (1 + cos b)/2."""
    return 0.5 * (1.0 + np.cos(slope_rad))
