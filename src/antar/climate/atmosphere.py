"""Seasonal shape and scenario change of the atmospheric inputs to the water balance.

The pipeline's wind, shortwave, longwave and humidity come from ERA5-Land as ANNUAL means: one value per cell, which was
repeated for every day of the year and left unchanged in every scenario. That removes the summer peak of evaporative
demand (shortwave in July is ~1.6x its annual mean) and with it most of the water deficit, and it makes the scenario
response of the water balance almost purely a temperature-and-rain effect (``scripts/sensitivity_atmospheric_forcing.py``:
mean climatic water deficit of broadleaf cells 14 mm with constant forcing vs 153 mm with a seasonal cycle; its change under
SSP5-8.5 by 2100 +4 mm vs +212 mm).

This module keeps ERA5-Land's cell-level annual mean (it carries the local topography) and adds

* the SEASONAL SHAPE of each variable from an observation-based daily product (ISIMIP3a GSWP3-W5E5, 2015-2019
  climatology, 0.5 degree), as a daily factor whose mean over the year is exactly 1, so the annual mean is unchanged;
* the SCENARIO CHANGE from ISIMIP3b (same bias-adjustment family): the ratio of the future 5-year monthly climatology to
  the model's own 2015-2019 monthly climatology, multiplied onto the shape.

    x(day) = x_annual(cell) * shape(day) * ratio(day)          shape: mean 1;  ratio: 1 in the baseline year

Monthly values are joined by cyclic linear interpolation between mid-month days so there are no steps at month edges.
Variables (keys): ``rs`` shortwave down, ``rl`` longwave down, ``wind`` wind speed, ``ea`` actual vapour pressure
(hurs/100 x es(tas), built daily before averaging so the temperature-driven rise in humidity is kept).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

BBOX = (43.4, 38.8, 46.7, 41.4)                     # lon_min, lat_min, lon_max, lat_max: the ISIMIP cutout
_LATS = np.arange(-89.75, 90, 0.5)
_LONS = np.arange(-179.75, 180, 0.5)
LAT_SEL = _LATS[(_LATS >= BBOX[1]) & (_LATS <= BBOX[3])]
LON_SEL = _LONS[(_LONS >= BBOX[0]) & (_LONS <= BBOX[2])]

BASELINE_WINDOW = (2015, 2019)
HORIZON_WINDOWS = {2050: (2048, 2052), 2080: (2078, 2082), 2100: (2096, 2100)}   # must equal run_future_projections.HORIZONS
VARIABLES = {"rs": "rsds", "rl": "rlds", "wind": "sfcwind", "ea": "ea"}
_DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
MID_DOY = np.cumsum(_DAYS) - _DAYS / 2.0 + 0.5       # mid-month day of year (1-based, 365-day year)


def cell_index(lats, lons):
    """Nearest ISIMIP 0.5 degree cell (row, col) of the cutout for each point."""
    lats, lons = np.atleast_1d(lats), np.atleast_1d(lons)
    row = np.array([np.argmin(np.abs(LAT_SEL - la)) for la in lats])
    col = np.array([np.argmin(np.abs(LON_SEL - lo)) for lo in lons])
    return row, col


def climatology(monthly, dates, y0, y1):
    """(n_months, ...) series keyed by 'YYYY-MM-01' -> (12, ...) mean of each calendar month over y0..y1."""
    years = np.array([int(str(d)[:4]) for d in dates])
    months = np.array([int(str(d)[5:7]) for d in dates])
    out = []
    for m in range(1, 13):
        sel = (years >= y0) & (years <= y1) & (months == m)
        if not sel.any():
            raise ValueError(f"no data for month {m} in {y0}-{y1}")
        out.append(np.nanmean(monthly[sel], axis=0))
    return np.stack(out)


def cyclic_interp(values12, doy):
    """12 monthly values -> the given days of year, linear between mid-month points, wrapping December to January."""
    v = np.asarray(values12, dtype=float)
    x = np.concatenate([MID_DOY - 365.0, MID_DOY, MID_DOY + 365.0])
    y = np.concatenate([v, v, v])
    return np.interp(np.asarray(doy, dtype=float), x, y)


def seasonal_shape(baseline12, doy):
    """Daily factor with mean exactly 1 over ``doy`` (the baseline's seasonal cycle, annual level removed)."""
    d = cyclic_interp(baseline12, doy)
    return d / d.mean()


def scenario_factor(baseline12, member_baseline12, member_future12, doy):
    """Shape of the observation-based baseline times the model's monthly future/baseline ratio."""
    mb = np.asarray(member_baseline12, dtype=float)
    ratio12 = np.where(mb > 0, np.asarray(member_future12, dtype=float) / np.where(mb > 0, mb, 1.0), 1.0)
    return seasonal_shape(baseline12, doy) * cyclic_interp(ratio12, doy)


class AtmosphereShape:
    """Loads the monthly ISIMIP files written by ``scripts/pull_isimip_atmosphere.py`` and serves daily factors per point."""

    def __init__(self, data_dir):
        self.dir = Path(data_dir)
        self._base = {}
        self._member = {}

    def _load(self, name):
        z = np.load(self.dir / name)
        return z["data"], z["dates"]

    def _baseline12(self, key):
        if key not in self._base:
            data, dates = self._load(f"ATMOS_3a_{VARIABLES[key]}__gswp3-w5e5__obsclim.npz")
            self._base[key] = climatology(data, dates, *BASELINE_WINDOW)
        return self._base[key]

    def _member12(self, key, gcm, scenario, window):
        k = (key, gcm, scenario, window)
        if k not in self._member:
            data, dates = self._load(f"ATMOS_3b_{VARIABLES[key]}__{gcm}__{scenario}.npz")
            self._member[k] = climatology(data, dates, *window)
        return self._member[k]

    def baseline_factors(self, lat, lon, doy):
        r, c = (int(a[0]) for a in cell_index(lat, lon))
        return {k: seasonal_shape(self._baseline12(k)[:, r, c], doy) for k in VARIABLES}

    def scenario_factors(self, lat, lon, doy, gcm, scenario, horizon):
        r, c = (int(a[0]) for a in cell_index(lat, lon))
        win = HORIZON_WINDOWS[horizon]
        return {k: scenario_factor(self._baseline12(k)[:, r, c], self._member12(k, gcm, scenario, BASELINE_WINDOW)[:, r, c],
                                   self._member12(k, gcm, scenario, win)[:, r, c], doy) for k in VARIABLES}

    @staticmethod
    def constant_factors(n_days):
        """The old behaviour: every day the annual mean (kept so earlier results can be reproduced and compared)."""
        return {k: np.ones(n_days) for k in VARIABLES}
