"""A node whose ERA5-Land wind is present but whose radiation or dew point is missing has a NaN water balance; because NaN >= 1 is False, the
hydraulic engine reported a hazard of exactly 0 (perfect viability) for three dense-grid nodes. Every band must be finite for a node to run."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("run_topohydro_grid_mask_test", SCRIPTS / "run_topohydro_grid.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_a_node_needs_every_era5_band():
    era5_complete = load().era5_complete
    ok = np.array([1.0, 1.0, 1.0, 1.0])
    wind = np.array([1.0, 1.0, 1.0, np.nan])
    ssrd = np.array([1.0, np.nan, 1.0, 1.0])
    strd, dew, pres = ok.copy(), ok.copy(), ok.copy()
    dew2 = np.array([1.0, 1.0, np.nan, 1.0])
    got = era5_complete(wind, ssrd, strd, dew2, pres)
    assert got.tolist() == [True, False, False, False]
    assert era5_complete(ok, ok, ok, ok, ok).all()


def test_nan_forcing_cannot_pass_as_zero_hazard():
    """The failure mechanism itself: a NaN index is never 'failed', so it must never be allowed to reach the Monte Carlo."""
    hfi = np.array([np.nan, np.nan])
    assert not np.any(hfi >= 1.0)            # reads as 'no failure' -> hazard 0: the reason the mask must be strict
