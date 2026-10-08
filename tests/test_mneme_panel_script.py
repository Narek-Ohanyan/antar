"""The MNEME panel script end to end on synthetic inputs: nothing here touches Drive or the real rasters."""
import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

pytest.importorskip("rasterio")
SCRIPTS = str(Path(__file__).resolve().parent.parent / "scripts")


@pytest.fixture
def script(tmp_path, monkeypatch):
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    mod = importlib.import_module("fit_mneme_hazard_panel")
    monkeypatch.setattr(mod, "PIXELS_PATH", tmp_path / "pixels.npz")
    monkeypatch.setattr(mod, "PANEL_PATH", tmp_path / "panel.csv.gz")
    monkeypatch.setattr(mod, "OUT_PATH_DENSE", tmp_path / "out.yaml")
    monkeypatch.setattr(mod, "CHECKPOINT_PATH_DENSE", tmp_path / "ckpt.json")
    monkeypatch.setattr(sys, "argv", ["fit_mneme_hazard_panel.py", "--dense"])
    return mod


def _install(script, monkeypatch, n_nodes, per_node, collapse_fraction, flag_all_disturbed=False, seed=0, noisy=True):
    rng = np.random.default_rng(seed)
    n = n_nodes * per_node
    years = list(range(2000, 2025))
    # noisy: random year-to-year variation, which on its own now and then satisfies the rule; otherwise a steady rise, which never does
    kn = 0.6 + 0.01 * rng.normal(size=(n, 25)) if noisy else 0.6 + 0.001 * np.arange(25)[None, :] + 0.01 * rng.random((n, 1))
    collapse = rng.random(n) < collapse_fraction
    when = rng.integers(12, 18, n)
    for i in np.flatnonzero(collapse):
        kn[i, when[i]:] = 0.1
    nd = np.ones((n, 25), dtype=bool)
    nd[:, 0] = True
    if flag_all_disturbed:
        nd[:, 1:] = False
    node = np.repeat(np.arange(n_nodes), per_node)
    np.savez_compressed(script.PIXELS_PATH, node=node, row=np.zeros(n, int), col=np.zeros(n, int), x=np.zeros(n), y=np.zeros(n),
                        lat=np.repeat(38.8 + 0.05 * np.arange(n_nodes), per_node), lon=np.repeat(43.5 + 0.05 * np.arange(n_nodes), per_node), forest_class=np.full(n, 31),
                        kndvi=kn, valid_count=np.full((n, 25), 20.0), no_disturbance=nd, elevation=rng.uniform(800, 2000, n), height_mean=rng.uniform(8, 25, n),
                        height_n_subpoints=np.full(n, 9), years=np.array(years), failed_nodes=np.array([], dtype=int))
    static = {"lats": np.zeros(n_nodes), "valid_soil": np.ones(n_nodes, bool), "soil": {"theta_fc": np.zeros(n_nodes)}, "token": None}
    monkeypatch.setattr(script, "extract_static_grid_inputs", lambda **kw: static)

    def fake_forcing(sub, year):
        k = len(sub["lats"])
        return [SimpleNamespace(cwd_mm={"pm_fao56": 100.0 + 5 * (year - 2010) + 3 * j}, vpd_24h_kpa=np.full(5, 1.0 + 0.01 * j), gdd_cumulative=np.array([0.0, 3000.0 + 10 * (year - 2010)]),
                                t_mean_c=np.full(5, 10.0 + 0.1 * (year - 2010))) for j in range(k)]
    monkeypatch.setattr(script, "compute_forcing_for_year", fake_forcing)
    return n


def test_a_panel_without_events_reports_the_rate_bound_and_no_hazard(script, monkeypatch):
    _install(script, monkeypatch, n_nodes=6, per_node=8, collapse_fraction=0.0, noisy=False)
    script.main()
    out = yaml.safe_load(script.OUT_PATH_DENSE.read_text())
    assert out["status"] == "too_few_events_for_a_hazard_fit" and out["n_events"] == 0 and out["n_pixels"] == 48
    assert out["n_person_years"] == 48 * 10 and out["event_rate_ci95"][0] == 0.0
    assert out["event_rate_upper_bound_95_one_sided"] == pytest.approx(1 - 0.05 ** (1 / 480))
    assert out["design"].startswith("forest pixels") and "glm_coefficients" not in out
    assert out["mirror_check"]["pixels_with_onset_of_decline"] == 0 and out["mirror_check"]["pixels_with_onset_of_rise"] == 0
    import pandas as pd
    panel = pd.read_csv(script.PANEL_PATH)
    assert len(panel) == 480 and {"z", "event", "candidate", "height_m", "cwd_mm", "at_risk"} <= set(panel.columns)
    assert script.CHECKPOINT_PATH_DENSE.exists()                                       # the climate checkpoint is kept for the next run
    assert out["year_effect_removed"]["n_events"] == 0


def test_every_collapse_is_counted_once_and_the_hazard_is_fitted_when_events_are_enough(script, monkeypatch):
    n = _install(script, monkeypatch, n_nodes=30, per_node=40, collapse_fraction=0.12, seed=1)
    script.main()
    out = yaml.safe_load(script.OUT_PATH_DENSE.read_text())
    assert out["n_events"] >= 90 and out["n_events"] <= n * 0.12 * 1.3                  # one onset per collapsed pixel, never one per year
    assert out["n_events_if_flags_ignored"] == out["n_events"]
    assert out["mirror_check"]["pixels_with_onset_of_decline"] > out["mirror_check"]["pixels_with_onset_of_rise"] + 100       # 12 % of 1200 pixels collapsed for real
    assert out["status"] == "fitted" and set(out["glm_coefficients"]) == set(out["design_columns"])
    assert out["n_person_years"] < n * 10                                               # years after an onset leave the risk set


def test_a_layer_that_flags_everything_as_disturbed_stops_the_script(script, monkeypatch):
    _install(script, monkeypatch, n_nodes=5, per_node=6, collapse_fraction=0.0, flag_all_disturbed=True, noisy=False)
    with pytest.raises(SystemExit, match="not a plausible record"):
        script.main()


def test_the_script_needs_the_dense_flag_and_the_extraction_first(script, monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["fit_mneme_hazard_panel.py"])
    with pytest.raises(SystemExit, match="dense"):
        script.main()
    monkeypatch.setattr(sys, "argv", ["fit_mneme_hazard_panel.py", "--dense"])
    with pytest.raises(SystemExit, match="extract_mneme_forest_pixels"):
        script.main()


def test_a_second_run_reuses_the_climate_checkpoint_and_a_changed_cell_set_does_not(script, monkeypatch):
    _install(script, monkeypatch, n_nodes=6, per_node=5, collapse_fraction=0.0, noisy=False)
    script.main()
    calls = []
    real = script.compute_forcing_for_year
    monkeypatch.setattr(script, "compute_forcing_for_year", lambda sub, year: calls.append(year) or real(sub, year))
    script.main()
    assert calls == []                                                                  # all ten years came from the checkpoint
    _install(script, monkeypatch, n_nodes=7, per_node=5, collapse_fraction=0.0, noisy=False)
    monkeypatch.setattr(script, "compute_forcing_for_year", lambda sub, year: calls.append(year) or real(sub, year))
    script.main()
    assert sorted(calls) == list(range(2010, 2020))                                     # a different set of cells: recomputed


def test_enough_events_that_are_not_clearly_declines_rather_than_rises_are_not_fitted(script, monkeypatch):
    _install(script, monkeypatch, n_nodes=30, per_node=60, collapse_fraction=0.0, seed=5)         # symmetric noise: persistent rises as often as declines
    monkeypatch.setattr(script, "MIN_EVENTS_PER_PREDICTOR", 1)                                  # so that the event count alone would allow a fit
    script.main()
    out = yaml.safe_load(script.OUT_PATH_DENSE.read_text())
    assert out["n_events"] >= 9
    assert out["status"] == "events_not_distinguishable_from_noise" and "glm_coefficients" not in out
    assert out["decline_to_rise_ratio"] < script.MIN_DECLINE_TO_RISE_RATIO
    assert "not distinguishable" in out["status_detail"] or "cannot be told" in out["status_detail"]
