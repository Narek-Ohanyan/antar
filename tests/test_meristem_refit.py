"""The MERISTEM refit with the real, group-specific water deficit: the right column is replaced for the right points, and the fitted model carries what its application needs."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = str(Path(__file__).resolve().parent.parent / "scripts")


@pytest.fixture
def mod():
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    return importlib.import_module("fit_meristem_adult_niche")


def _features(mod, rng, n_per=12, n_bg=30):
    out = {}
    for sp in mod.TARGET_SPECIES:
        out[sp] = (rng.uniform(39, 41, (n_per, 2)), rng.normal(size=(n_per, 8)))
    out["background"] = (rng.uniform(39, 41, (n_bg, 2)), rng.normal(size=(n_bg, 8)))
    return out


def _write_cwd(mod, path, features):
    labels = np.array([k for k in list(mod.TARGET_SPECIES) + ["background"] for _ in range(len(features[k][0]))])
    n = len(labels)
    cwd = {g: 1000.0 * (j + 1) + np.arange(n, dtype=float) for j, g in enumerate(["mesic_diffuse_porous_broadleaf", "ring_porous_oak", "pine", "juniper_arid_conifer"])}
    np.savez_compressed(path, lats=np.zeros(n), lons=np.zeros(n), groups=labels, **{f"real_cwd_mm__{g}": v for g, v in cwd.items()})
    return labels, cwd


def test_each_species_gets_the_water_deficit_of_its_own_group_at_its_records_and_at_the_background(mod, tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    feats = _features(mod, rng)
    monkeypatch.setattr(mod, "REAL_CWD_PATH", tmp_path / "cwd.npz")
    labels, cwd = _write_cwd(mod, tmp_path / "cwd.npz", feats)
    out = mod.with_real_cwd(feats)
    col = mod.CWD_COLUMN
    # Quercus iberica (4th species) belongs to the oak group; its records are the 4th block of 12
    ll, X, bll, bX = out["Quercus iberica"]
    assert X[:, col].tolist() == cwd["ring_porous_oak"][36:48].tolist()
    assert bX[:, col].tolist() == cwd["ring_porous_oak"][len(labels) - 30:].tolist()               # the shared background carries the same group's values
    assert np.array_equal(np.delete(X, col, axis=1), np.delete(feats["Quercus iberica"][1], col, axis=1))   # nothing else changed
    f = out["Fagus orientalis"]
    assert f[1][:, col].tolist() == cwd["mesic_diffuse_porous_broadleaf"][0:12].tolist()
    assert f[3][:, col].tolist() != out["Pinus kochiana"][3][:, col].tolist()                          # a different group sees different background values


def test_a_water_deficit_file_that_does_not_line_up_with_the_feature_cache_is_refused(mod, tmp_path, monkeypatch):
    rng = np.random.default_rng(1)
    feats = _features(mod, rng)
    monkeypatch.setattr(mod, "REAL_CWD_PATH", tmp_path / "cwd.npz")
    _write_cwd(mod, tmp_path / "cwd.npz", feats)
    shorter = dict(feats)
    shorter["background"] = (feats["background"][0][:-3], feats["background"][1][:-3])
    with pytest.raises(ValueError, match="points"):
        mod.with_real_cwd(shorter)
    swapped = dict(feats)
    swapped["Fagus orientalis"], swapped["Carpinus betulus"] = feats["Carpinus betulus"], feats["Fagus orientalis"]       # same counts, labels in the wrong blocks
    np.savez_compressed(tmp_path / "cwd2.npz", **{k: v for k, v in np.load(tmp_path / "cwd.npz", allow_pickle=True).items()})
    monkeypatch.setattr(mod, "REAL_CWD_PATH", tmp_path / "cwd2.npz")
    z = dict(np.load(tmp_path / "cwd2.npz", allow_pickle=True))
    z["groups"] = z["groups"][::-1].copy()
    np.savez_compressed(tmp_path / "cwd2.npz", **z)
    with pytest.raises(ValueError, match="records"):
        mod.with_real_cwd(feats)


def test_the_fitted_model_keeps_full_precision_coefficients_and_a_presence_threshold(mod):
    rng = np.random.default_rng(2)
    pres = rng.normal(1.0, 1.0, (80, 8))
    bg = rng.normal(0.0, 1.0, (300, 8))
    ll_p, ll_b = rng.uniform(39, 41, (80, 2)), rng.uniform(39, 41, (300, 2))
    r = mod.fit_one_species("x", ll_p, pres, ll_b, bg)
    assert r["status"] == "fitted" and set(r["coefficients"]) == set(mod.FEATURE_NAMES) and "cwd_mm" in r["coefficients"]
    assert any(len(repr(v).split(".")[-1]) > 7 for v in r["coefficients"].values())                         # not rounded to five decimals
    from antar.niche import application as A
    score = A.logistic_score(list(r["coefficients"].values()), r["intercept"], pres)
    thr = r["presence_threshold"]
    assert thr["omission"] == 0.10 and thr["score"] == pytest.approx(A.presence_threshold(score, 0.10))
    assert thr["n_records_below"] == int((score < thr["score"]).sum()) and 6 <= thr["n_records_below"] <= 10


def test_a_species_with_too_few_records_is_skipped_not_forced(mod):
    rng = np.random.default_rng(3)
    r = mod.fit_one_species("x", rng.uniform(39, 41, (8, 2)), rng.normal(size=(8, 8)), rng.uniform(39, 41, (50, 2)), rng.normal(size=(50, 8)))
    assert r["status"] == "skipped_insufficient_data" and "presence_threshold" not in r
