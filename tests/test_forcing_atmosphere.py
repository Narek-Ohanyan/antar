"""The cell forcing path with seasonal atmosphere factors (scripts/run_topohydro_grid.py::_forcing_cell), on a synthetic cell."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import run_topohydro_grid as R  # noqa: E402
from antar.climate.atmosphere import AtmosphereShape, seasonal_shape  # noqa: E402

DOY = np.arange(1, 366)
MONTH = np.array([int(m) for m in np.repeat(np.arange(1, 13), [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])])


def synthetic(atmos_factors):
    t = 6.0 + 14.0 * np.sin(2 * np.pi * (DOY - 110) / 365)                 # continental annual cycle, mean 6 C
    p = np.where((MONTH >= 5) & (MONTH <= 9), 0.9, 2.2)                     # drier summer
    static = {"elevation": np.array([1500.0]), "slope": np.array([5.0]), "aspect": np.array([180.0]), "concavity": np.array([0.0]),
              "z_ref_m": np.array([1500.0]), "lats": np.array([40.2]),
              "soil": {"theta_sat": np.array([0.45]), "psi_sat_mpa": np.array([-0.0005]), "b_clapp_hornberger": np.array([6.0]),
                       "theta_fc": np.array([0.30]), "theta_lim": np.array([0.12])}}
    climate = {"doy": DOY, "month": MONTH, "t_mean_ref_all": t[:, None], "t_max_ref_all": (t + 6)[:, None], "t_min_ref_all": (t - 6)[:, None],
               "p_ref_all": p[:, None], "gamma_k_per_m": np.full(12, -0.0055), "precip_gradient_per_m": np.zeros(12),
               "ssrd": np.array([16.7e6]), "strd": np.array([24.7e6]), "ea_ref_kpa": np.array([0.7]), "u2_m_s": np.array([1.2]),
               "pressure_kpa_era5": np.array([85.0]), "atmos_factors": atmos_factors}
    return R._forcing_cell(static, climate, np.array([0.18 * 1500.0]), 0)


def shape_factors():
    rs = seasonal_shape(np.array([4, 7, 12, 17, 22, 26, 27, 23, 17, 10, 5, 3.0]), DOY)
    ea = seasonal_shape(np.array([0.3, 0.35, 0.5, 0.7, 0.9, 1.2, 1.6, 1.5, 1.1, 0.8, 0.5, 0.35]), DOY)
    ones = np.ones(365)
    return {"rs": rs, "rl": seasonal_shape(np.array([230, 235, 250, 270, 295, 320, 345, 340, 310, 280, 250, 235.0]), DOY), "wind": ones, "ea": ea}


def test_unit_factors_reproduce_the_constant_behaviour_exactly():
    a = synthetic(None)                                                     # no factors: constant annual means
    b = synthetic([AtmosphereShape.constant_factors(365)])                  # explicit all-ones factors
    assert a.cwd_mm["pm_fao56"] == b.cwd_mm["pm_fao56"]
    assert np.array_equal(a.psi_soil_mpa["pm_fao56"], b.psi_soil_mpa["pm_fao56"])


def test_seasonal_radiation_and_humidity_raise_the_summer_water_deficit():
    const = synthetic(None)
    seasonal = synthetic([shape_factors()])
    assert seasonal.cwd_mm["pm_fao56"] > 2.0 * const.cwd_mm["pm_fao56"] + 5          # summer demand no longer averaged away
    assert seasonal.psi_soil_mpa["pm_fao56"].min() < const.psi_soil_mpa["pm_fao56"].min()
    summer = (MONTH >= 6) & (MONTH <= 8)
    assert seasonal.pet_mm["pm_fao56"][summer].mean() > const.pet_mm["pm_fao56"][summer].mean()
    winter = (MONTH == 12) | (MONTH <= 2)
    assert seasonal.pet_mm["pm_fao56"][winter].mean() < const.pet_mm["pm_fao56"][winter].mean()   # and winter demand falls


def test_annual_mean_of_the_inputs_is_unchanged_by_the_shape():
    f = shape_factors()
    for k, v in f.items():
        assert v.mean() == pytest.approx(1.0, abs=1e-12), k


def test_cache_key_is_stable_and_round_trips(tmp_path):
    a, b = np.array([40.1, 40.2]), np.array([44.5, 44.6])
    assert R._points_key(a, b) == R._points_key(a + 1e-9, b - 1e-9)           # float noise does not split the cache
    assert R._points_key(a, b) != R._points_key(a, b + 0.01)
    R._save_npz(tmp_path / "x.npz", {"v": np.arange(5.0)})
    assert np.array_equal(np.load(tmp_path / "x.npz")["v"], np.arange(5.0))
