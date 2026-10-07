"""A node whose water balance is NaN must never run. NaN >= 1 is False, so the hydraulic engine reported a hazard of exactly 0 (perfect viability) for
three dense-grid nodes whose soil parameters were NaN although their clay value was present. The masks now require every input to be finite."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))


def load(name):
    spec = importlib.util.spec_from_file_location(name + "_masktest", SCRIPTS / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_a_node_needs_every_soil_parameter():
    soil_complete = load("run_topohydro_grid").soil_complete
    soil = {"theta_sat": np.array([0.4, 0.4, np.nan, 0.4]), "theta_fc": np.array([0.3, 0.3, 0.3, 0.3]), "b": np.array([5.0, np.nan, 5.0, 5.0])}
    clay = np.array([20.0, 20.0, 20.0, np.nan])
    assert soil_complete(soil, clay).tolist() == [True, False, False, False]


def test_a_node_needs_every_era5_band():
    era5_complete = load("run_topohydro_grid").era5_complete
    ok = np.ones(4)
    wind = np.array([1.0, 1.0, 1.0, np.nan])
    ssrd = np.array([1.0, np.nan, 1.0, 1.0])
    dew = np.array([1.0, 1.0, np.nan, 1.0])
    assert era5_complete(wind, ssrd, ok, dew, ok).tolist() == [True, False, False, False]
    assert era5_complete(ok, ok, ok, ok, ok).all()


def test_nan_cannot_pass_as_a_zero_hazard():
    """The failure mechanism: a NaN index is never 'failed', so NaN must be stopped before the Monte Carlo."""
    assert not np.any(np.array([np.nan, np.nan]) >= 1.0)


def test_dropping_invalid_nodes_recomputes_the_summaries():
    drop = load("drop_invalid_nodes")
    bad = {(40.0, 44.0)}
    topo = {"points": [{"lat": 40.0, "lon": 44.0, "status": "skipped_nonfinite_water_balance"}, {"lat": 40.1, "lon": 44.0, "status": "ok"}]}
    assert drop.invalid_nodes(topo) == bad
    xylem = {"groups": {"g": {"cells": [{"lat": 40.0, "lon": 44.0, "h_mech_mean": 0.0}, {"lat": 40.1, "lon": 44.0, "h_mech_mean": 0.2}]}}}
    assert drop.clean_xylem(xylem, bad) == 1
    g = xylem["groups"]["g"]
    assert g["n_cells"] == 1 and abs(g["h_mech_mean_across_cells"] - 0.2) < 1e-12
    assert drop.clean_xylem(xylem, bad) == 0                                   # idempotent
    scen = {"members": {"m": {"groups": {"g": {"cells": [{"lat": 40.0, "lon": 44.0, "viability_mean": 1.0, "h_mech_mean": 0.0, "refugium_score": 1.0},
                                                         {"lat": 40.1, "lon": 44.0, "viability_mean": 0.9, "h_mech_mean": 0.1, "refugium_score": 0.8}]}},
                          "climate": [{"lat": 40.0, "lon": 44.0}, {"lat": 40.1, "lon": 44.0}], "generic": [{"lat": 40.0, "lon": 44.0}]}}}
    assert drop.clean_scenarios(scen, bad) == 1
    m = scen["members"]["m"]
    assert m["groups"]["g"]["mean_viability"] == 0.9 and len(m["climate"]) == 1 and m["generic"] == []
