"""Checks on the committed map layers (ui/assets/map): the geography has to be right, because the
map has no basemap to reveal a misplaced layer."""
import json
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

ASSETS = Path(__file__).resolve().parent.parent / "ui" / "assets" / "map"
R = 6378137.0
pytestmark = pytest.mark.skipif(not (ASSETS / "grid.json").exists(), reason="map assets not built")

# real towns and the marz they are in (Yerevan is its own unit in the boundary file)
TOWNS = {"Yerevan": (40.18, 44.51, "Yerevan"), "Gyumri": (40.79, 43.85, "Shirak"), "Goris": (39.51, 46.34, "Syunik"),
         "Vanadzor": (40.81, 44.49, "Lori"), "Ijevan": (40.88, 45.15, "Tavush"), "Jermuk": (39.84, 45.67, "Vayots Dzor"),
         "Artashat": (39.96, 44.55, "Ararat"), "Armavir": (40.15, 43.87, "Armavir"), "Aparan": (40.59, 44.36, "Aragatsotn"),
         "Hrazdan": (40.53, 44.77, "Kotayk"), "Gavar": (40.36, 45.13, "Gegharkunik")}


@pytest.fixture(scope="module")
def grid():
    return json.loads((ASSETS / "grid.json").read_text())


def load(name):
    return np.array(Image.open(ASSETS / f"{name}.png"))


def pixel(grid, lat, lon):
    x, y = R * np.radians(lon), R * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))
    return int(np.floor((grid["y_top"] - y) / grid["px_m"])), int(np.floor((x - grid["x0"]) / grid["px_m"]))


def test_towns_fall_in_their_marz(grid):
    region = load("region")
    names = {r["id"]: r["name"] for r in grid["regions"]}
    for town, (lat, lon, expected) in TOWNS.items():
        row, col = pixel(grid, lat, lon)
        assert names[int(region[row, col])] == expected, town


def test_lake_sevan_is_water_and_outside_points_are_outside(grid):
    row, col = pixel(grid, 40.35, 45.35)
    assert load("water")[row, col] > 250
    row, col = pixel(grid, 41.70, 44.80)                       # Tbilisi
    assert load("region")[row, col] == 0


def test_shapes_and_area(grid):
    for name in ["region", "forest", "woodland", "water", "forest_broadleaf", "forest_oak", "forest_pine", "forest_juniper"]:
        assert load(name).shape == (grid["height"], grid["width"]), name
    assert abs(grid["country_area_km2"] / grid["official_area_km2"] - 1) < 0.03
    assert len(grid["regions"]) == 11
    assert abs(sum(r["area_km2"] for r in grid["regions"]) - grid["country_area_km2"]) < 1.0
    inside = load("region") > 0
    assert not (inside[0].any() or inside[-1].any() or inside[:, 0].any() or inside[:, -1].any())


def test_cover_only_inside_the_country_and_totals_are_plausible(grid):
    inside = load("region") > 0
    for name in ["forest", "woodland", "water", "forest_oak"]:
        assert not load(name)[~inside].any(), name
    a = grid["cover_area_km2"]
    assert 2500 < a["forest"] < 4000 and 1500 < a["woodland"] < 3000       # national map: ~3.2k km2 forest, ~2.2k km2 woodland
    assert 1200 < a["water"] < 1600                                           # Lake Sevan alone is ~1240 km2
    assert a["forest_pine"] < 20                                              # class 36 is tiny (about 5 km2) in the map


def test_borders_are_valid_and_inside_the_grid(grid):
    gj = json.loads((ASSETS / "borders.geojson").read_text())
    (s, w), (n, e) = grid["bounds_latlon"]
    kinds = [f["properties"]["kind"] for f in gj["features"]]
    assert kinds.count("country") >= 1 and kinds.count("marz") >= 11
    for f in gj["features"]:
        pts = np.array(f["geometry"]["coordinates"])
        assert pts[:, 0].min() >= w - 0.01 and pts[:, 0].max() <= e + 0.01 and pts[:, 1].min() >= s - 0.01 and pts[:, 1].max() <= n + 0.01
        assert np.allclose(pts[0], pts[-1])


def test_geometry_helpers():
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    from build_map_assets import chaikin, inv_merc, merc_x, merc_y, rdp
    lon, lat = inv_merc(merc_x(44.5), merc_y(40.2))
    assert (lon, lat) == pytest.approx((44.5, 40.2))
    line = np.column_stack([np.arange(100.0), np.zeros(100)])
    assert len(rdp(line, 0.1)) == 2                                           # a straight line collapses to its ends
    zig = np.array([[0, 0], [1, 0.05], [2, 0], [3, 5.0], [4, 0]])
    assert len(rdp(zig, 0.5)) == 4                                            # but a real spike survives
    sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], dtype=float)
    sm = chaikin(sq, iters=2)
    assert np.allclose(sm[0], sm[-1]) and sm[:, 0].min() >= 0 and sm[:, 0].max() <= 1   # smoothing stays inside the hull
