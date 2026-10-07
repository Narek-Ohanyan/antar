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
