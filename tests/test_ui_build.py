"""The UI data bundle must be internally consistent whatever fitted outputs exist: every scenario series lines up
with its member list and cell list, nothing is NaN, and the map plan agrees with what was actually exported."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "ui" / "data"


@pytest.fixture(scope="module")
def built():
    spec = importlib.util.spec_from_file_location("ui_build_data", ROOT / "ui" / "build_data.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.build()
    manifest = json.loads((DATA / "manifest.json").read_text())
    grids = {g: json.loads((DATA / f["file"]).read_text()) for g, f in manifest["grids"].items()}
    return manifest, grids


def test_scenario_series_line_up(built):
    manifest, grids = built
    for gid, g in grids.items():
        n, sc = g["n_cells"], g["scenarios"]
        for qid, mat in sc.get("fp_series", {}).items():
            assert len(mat) == len(sc["members"]) == 45, (gid, qid)
            assert all(len(row) == n for row in mat), (gid, qid)
        for qid, mat in sc.get("tl_series", {}).items():
            assert len(mat) == len(sc["treeline_members"]) == 45, (gid, qid)
            assert all(len(row) == n for row in mat), (gid, qid)
        for lid, row in g["layers"].items():
            assert len(row) == n, (gid, lid)
        assert all(0 <= la <= 90 for la in g["lat"])


def test_no_nan_in_exported_values(built):
    _, grids = built
    for g in grids.values():
        rows = list(g["layers"].values()) + [r for m in g["scenarios"].get("fp_series", {}).values() for r in m] \
            + [r for m in g["scenarios"].get("tl_series", {}).values() for r in m]
        for row in rows:
            assert all(v is None or math.isfinite(v) for v in row)


def test_every_cell_is_inside_armenia(built):
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from antar.io.armenia_mask import inside_armenia
    _, grids = built
    for g in grids.values():
        assert inside_armenia(g["lat"], g["lon"]).all()


def test_map_plan_agrees_with_exports(built):
    manifest, grids = built
    plan = manifest["map_plan"]
    ids = [p["id"] for p in plan]
    assert len(ids) == len(set(ids))
    series = set()
    for g in grids.values():
        series |= set(g["scenarios"].get("fp_series", {})) | set(g["scenarios"].get("tl_series", {}))
    layers = set(manifest["layers"])
    for p in plan:
        assert p["scenario"] == (p["id"] in series), p["id"]
        if p["baseline"] is not None:
            assert p["baseline"] == (p["id"] in layers), p["id"]
        assert p["engine"] and p["scenario_from"]


def test_series_and_baseline_share_a_grid_when_both_exist(built):
    """A change map is scenario minus baseline: both must live on the same grid or the change cannot be shown."""
    _, grids = built
    for gid, g in grids.items():
        for qid in g["scenarios"].get("fp_series", {}):
            if qid in g["layers"]:
                assert any(v is not None for v in g["layers"][qid]), (gid, qid)
