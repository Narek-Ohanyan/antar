"""Cell areas: closed-form values, additivity, the cosine shrink with latitude, and the real dense-grid cell."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from antar.decision.units import EARTH_RADIUS_KM, cell_area_ha  # noqa: E402


def test_a_one_degree_cell_at_the_equator_is_12364_square_kilometres():
    km2 = cell_area_ha(0.0, 1.0, 1.0) / 100.0
    assert km2 == pytest.approx(12364.4, rel=1e-4)                    # R^2 * (pi/180) * 2 sin(0.5 deg)
    assert km2 == pytest.approx(EARTH_RADIUS_KM ** 2 * np.radians(1.0) * 2 * np.sin(np.radians(0.5)), rel=1e-12)


def test_stacked_cells_add_up_to_the_taller_cell():
    parts = sum(cell_area_ha(40.0 - 0.45 + 0.1 * k, 0.1, 0.5) for k in range(10))
    assert parts == pytest.approx(cell_area_ha(40.0, 1.0, 0.5), rel=1e-12)


def test_a_cell_shrinks_with_latitude_like_the_cosine():
    a0, a40 = cell_area_ha(0.0, 0.05, 0.05), cell_area_ha(40.5, 0.05, 0.05)
    assert a40 / a0 == pytest.approx(np.cos(np.radians(40.5)), rel=1e-3)
    lat = np.linspace(0, 80, 17)
    assert np.all(np.diff(cell_area_ha(lat, 0.05, 0.05)) < 0)         # strictly smaller towards the pole
    assert cell_area_ha(lat, 0.05, 0.05).shape == lat.shape


def test_a_node_of_the_dense_grid_stands_for_about_three_thousand_hectares():
    import run_topohydro_grid as g
    dlat = (g.BBOX[3] - g.BBOX[1]) / g.CHELSA_GRID_SHAPE[0] * g.DENSE_STRIDE
    dlon = (g.BBOX[2] - g.BBOX[0]) / g.CHELSA_GRID_SHAPE[1] * g.DENSE_STRIDE
    assert dlat == pytest.approx(7 * 30 / 3600)                       # seven 30-arcsecond pixels
    ha = cell_area_ha([39.0, 40.5, 41.2], dlat, dlon)
    assert np.all((ha > 2900) & (ha < 3700)) and ha[0] > ha[1] > ha[2]
    assert dlat * 111.195 == pytest.approx(6.49, abs=0.01)             # about 6.5 km north-south, as the grid's own description says


# ---- land cover per cell
from antar.decision.units import cell_land_shares  # noqa: E402

GRID = {"x0": 4_800_000.0, "y_top": 5_100_000.0, "px_m": 650.0, "earth_radius_m": 6378137.0}


def pixel_centre_latlon(row, col):
    x, y = GRID["x0"] + (col + 0.5) * GRID["px_m"], GRID["y_top"] - (row + 0.5) * GRID["px_m"]
    return np.degrees(2 * np.arctan(np.exp(y / GRID["earth_radius_m"])) - np.pi / 2), np.degrees(x / GRID["earth_radius_m"])


def cell_deg_for(rows, cols, lat):
    """Degrees that span just under `rows` x `cols` map pixels at this latitude (so the cell covers exactly those pixels). In Web Mercator a projected
    metre of northing is cos(lat) of a ground metre, so a pixel is px*cos(lat)/R radians of latitude and px/R of longitude."""
    R, px = GRID["earth_radius_m"], GRID["px_m"]
    return np.degrees((rows - 0.2) * px * np.cos(np.radians(lat)) / R), np.degrees((cols - 0.2) * px / R)


def test_land_shares_are_the_mean_over_the_pixels_a_cell_covers():
    region = np.ones((30, 30), dtype=np.uint8)
    forest = np.zeros((30, 30)); forest[10:15, 10:15] = 1.0                    # a 5 x 5 block of forest
    wood = np.full((30, 30), 0.25)
    lat, lon = pixel_centre_latlon(12, 12)                                       # centre of the block
    dlat, dlon = cell_deg_for(5, 5, lat)
    s = cell_land_shares([lat], [lon], dlat, dlon, region, {"forest": forest, "woodland": wood}, GRID)
    assert s["forest"][0] == pytest.approx(1.0) and s["woodland"][0] == pytest.approx(0.25) and s["coverage"][0] == pytest.approx(1.0)
    dlat, dlon = cell_deg_for(10, 10, lat)                                       # a cell twice as wide: a quarter of it is forest
    s2 = cell_land_shares([lat], [lon], dlat, dlon, region, {"forest": forest}, GRID)
    assert 0.2 <= s2["forest"][0] <= 0.3


def test_pixels_outside_the_country_do_not_count_and_a_cell_outside_has_no_value():
    region = np.zeros((30, 30), dtype=np.uint8); region[:, :15] = 1              # the left half is the country
    forest = np.where(np.arange(30)[None, :] < 15, 0.8, 0.0) * np.ones((30, 1))  # the right half has no forest, but it is not the country
    lat, lon = pixel_centre_latlon(15, 14)
    dlat, dlon = cell_deg_for(6, 6, lat)
    s = cell_land_shares([lat], [lon], dlat, dlon, region, {"forest": forest}, GRID)
    assert s["forest"][0] == pytest.approx(0.8) and 0.3 < s["coverage"][0] < 0.7        # the mean ignores the outside; the coverage reports how much of the cell it was
    far_lat, far_lon = pixel_centre_latlon(15, 25)
    s = cell_land_shares([far_lat], [far_lon], dlat, dlon, region, {"forest": forest}, GRID)
    assert np.isnan(s["forest"][0]) and s["coverage"][0] == 0.0
    off = cell_land_shares([10.0], [10.0], dlat, dlon, region, {"forest": forest}, GRID)             # far outside the raster altogether
    assert np.isnan(off["forest"][0])
