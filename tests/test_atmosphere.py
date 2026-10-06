"""Seasonal shape and scenario change of the atmospheric inputs (src/antar/climate/atmosphere.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

from antar.climate import atmosphere as A

DOY = np.arange(1, 366)


def monthly_series(base12, years, scale_by_year=None):
    """(n_months, 5, 6) series and 'YYYY-MM-01' dates for the given years."""
    data, dates = [], []
    for y in years:
        for m in range(12):
            v = base12[m] * (1.0 if scale_by_year is None else scale_by_year(y))
            data.append(np.full((5, 6), v))
            dates.append(f"{y}-{m + 1:02d}-01")
    return np.array(data, dtype=np.float32), np.array(dates)


def test_mid_month_days_and_wraparound():
    assert A.MID_DOY[0] == pytest.approx(16.0) and A.MID_DOY[6] == pytest.approx(197.0)
    flat = A.cyclic_interp(np.full(12, 3.0), DOY)
    assert np.allclose(flat, 3.0)
    jan_dec = np.zeros(12)
    jan_dec[0], jan_dec[11] = 10.0, 10.0
    v = A.cyclic_interp(jan_dec, np.array([1, 365]))
    assert v[0] == pytest.approx(v[1], rel=0.05) and v.min() > 8          # continuous across new year, no cliff


def test_shape_has_mean_one_and_keeps_the_seasonal_cycle():
    rs = np.array([4, 7, 12, 17, 22, 26, 27, 23, 17, 10, 5, 3], dtype=float)      # a continental shortwave cycle
    s = A.seasonal_shape(rs, DOY)
    assert s.mean() == pytest.approx(1.0, abs=1e-12)
    assert s.max() / s.min() > 5                                                    # summer is far above winter
    assert DOY[np.argmax(s)] in range(180, 215)                                     # peaks in July


def test_annual_level_is_untouched_by_the_shape():
    """x(day) = annual * shape: whatever the shape, the year's mean is the ERA5 annual value."""
    annual = 16.7
    for base in (np.linspace(2, 30, 12), np.full(12, 5.0), np.array([9, 8, 7, 6, 5, 4, 3, 4, 5, 6, 7, 8], dtype=float)):
        assert (annual * A.seasonal_shape(base, DOY)).mean() == pytest.approx(annual, rel=1e-12)


def test_scenario_factor_baseline_is_identity_and_ratio_scales():
    base = np.linspace(3, 25, 12)
    same = A.scenario_factor(base, base, base, DOY)
    assert np.allclose(same, A.seasonal_shape(base, DOY))                           # future == baseline -> just the shape
    doubled = A.scenario_factor(base, base, 2 * base, DOY)
    assert np.allclose(doubled, 2 * A.seasonal_shape(base, DOY))
    summer_only = base.copy()
    summer_only[5:8] *= 1.5
    f = A.scenario_factor(base, base, summer_only, DOY) / A.seasonal_shape(base, DOY)
    assert f[170:230].mean() > 1.3 and f[:60].mean() == pytest.approx(1.0, abs=0.01)


def test_cell_index_matches_the_isimip_grid_over_armenia():
    assert len(A.LAT_SEL) == 5 and len(A.LON_SEL) == 6                              # the (5, 6) cutout confirmed in the log
    r, c = A.cell_index([40.18, 39.25], [44.51, 43.75])                              # Yerevan; the south-west corner cell
    assert (r[0], c[0]) == (2, 2) and (r[1], c[1]) == (0, 0)      # 44.51E is 0.24 deg from the 44.75 centre, 0.26 from 44.25


def test_horizons_agree_with_the_projection_script():
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import run_future_projections as F
    assert A.HORIZON_WINDOWS == {k: tuple(v) for k, v in F.HORIZONS.items()}
    assert A.BASELINE_WINDOW == tuple(F.BASELINE_WINDOW)


def test_atmosphere_shape_from_files(tmp_path):
    rs12 = np.array([4, 7, 12, 17, 22, 26, 27, 23, 17, 10, 5, 3], dtype=float)
    for var in ("rsds", "rlds", "sfcwind", "ea"):
        d, t = monthly_series(rs12, range(2011, 2020))
        np.savez(tmp_path / f"ATMOS_3a_{var}__gswp3-w5e5__obsclim.npz", data=d, dates=t)
        d, t = monthly_series(rs12, range(2015, 2101), scale_by_year=lambda y: 1.0 + 0.2 * (y - 2015) / 85)   # +20% by 2100
        np.savez(tmp_path / f"ATMOS_3b_{var}__gfdl-esm4__ssp585.npz", data=d, dates=t)
    atm = A.AtmosphereShape(tmp_path)
    base = atm.baseline_factors(40.2, 44.5, DOY)
    assert all(base[k].mean() == pytest.approx(1.0) for k in A.VARIABLES)
    fut = atm.scenario_factors(40.2, 44.5, DOY, "gfdl-esm4", "ssp585", 2100)
    ratio = fut["rs"].mean() / base["rs"].mean()
    assert 1.15 < ratio < 1.22                                                       # the +20% trend, read at 2096-2100 vs 2015-2019
    const = A.AtmosphereShape.constant_factors(365)
    assert all(np.all(v == 1.0) for v in const.values())
