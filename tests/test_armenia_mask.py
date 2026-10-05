"""The Armenia mask: known places, and that its area matches Armenia's official area."""
import numpy as np
import pytest

from antar.io import armenia_mask as am

pytestmark = pytest.mark.skipif(not am._PATH.exists(), reason="armenia mask not built")


def test_known_places():
    assert am.inside_armenia(40.18, 44.51)       # Yerevan
    assert am.inside_armenia(40.30, 45.35)       # Lake Sevan (water bodies count as national territory)
    assert am.inside_armenia(39.55, 46.30)       # Syunik highlands
    assert not am.inside_armenia(41.72, 44.79)   # Tbilisi, Georgia
    assert not am.inside_armenia(39.20, 45.40)   # Nakhchivan
    assert not am.inside_armenia(38.50, 43.80)   # eastern Turkey / Lake Van side


def test_off_grid_is_outside():
    assert not am.inside_armenia(10.0, 44.0)
    assert not am.inside_armenia(40.0, 80.0)


def test_area_matches_official_figure():
    lat = am.BBOX[3] - (np.arange(am.SHAPE[0]) + 0.5) * (am.BBOX[3] - am.BBOX[1]) / am.SHAPE[0]
    cell_km2 = (111.195 * (am.BBOX[3] - am.BBOX[1]) / am.SHAPE[0]) * (
        111.195 * np.cos(np.radians(lat)) * (am.BBOX[2] - am.BBOX[0]) / am.SHAPE[1])
    area = float((am.load_mask() * cell_km2[:, None]).sum())
    assert abs(area - 29743) / 29743 < 0.02  # within 2% of Armenia's official 29,743 km2


def test_vectorised_matches_scalar():
    lats, lons = np.array([40.18, 41.72, 39.55]), np.array([44.51, 44.79, 46.30])
    assert list(am.inside_armenia(lats, lons)) == [True, False, True]
