"""Treeline-change helper: the growing-season temperature recomputation the UI's treeline layer
rests on. Synthetic series with known answers -- no real data needed."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from compute_treeline_change import growing_season_temperature  # noqa: E402

from antar.niche.growth import potential_treeline_elevation  # noqa: E402

GAMMA = np.full(12, -0.0065)  # K/m, uniform across months for a clean hand check
DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MONTH = np.concatenate([np.full(d, m + 1) for m, d in enumerate(DAYS_PER_MONTH)])


def _series(summer_c, shoulder_c, winter_c):
    t = np.full(MONTH.size, winter_c, dtype=float)
    t[(MONTH >= 4) & (MONTH <= 9)] = summer_c
    t[(MONTH == 3) | (MONTH == 10)] = shoulder_c
    return t


def test_uniform_delta_shifts_treeline_by_delta_over_lapse_rate():
    # Season days are far above the 0.9 degC threshold and winter days far below it, so a +2 K
    # delta cannot change which days count: GST rises by exactly 2 K and the treeline by 2/0.0065 m.
    t = _series(summer_c=12.0, shoulder_c=-10.0, winter_c=-10.0)
    z = 1800.0
    gst0 = growing_season_temperature(t, MONTH, z, z, GAMMA)
    gst1 = growing_season_temperature(t, MONTH, z, z, GAMMA, delta_by_month=np.full(12, 2.0))
    assert abs((gst1 - gst0) - 2.0) < 1e-9
    shift = potential_treeline_elevation(gst1, z, -0.0065) - potential_treeline_elevation(gst0, z, -0.0065)
    assert abs(shift - 2.0 / 0.0065) < 1e-6


def test_season_lengthening_can_lower_gst_despite_warming():
    # March/October sit at 0 degC: excluded at baseline (< 0.9), included once +2 K lifts them to
    # 2 degC. GST is the MEAN over qualifying days (Eq. 8.3), and those newly counted shoulder days
    # are far colder than the 14 degC summer days, so GST FALLS (12.0 -> 10.96) even though every
    # single day warmed by 2 K. A real property of the framework's own definition, not a bug here:
    # it means treeline change is not guaranteed monotone in warming, and it is exactly what a
    # "GST_2019 + mean monthly delta" shortcut would get wrong (it would say 14.0).
    t = _series(summer_c=12.0, shoulder_c=0.0, winter_c=-10.0)
    z = 1800.0
    gst0 = growing_season_temperature(t, MONTH, z, z, GAMMA)
    gst1 = growing_season_temperature(t, MONTH, z, z, GAMMA, delta_by_month=np.full(12, 2.0))
    assert gst0 == 12.0  # baseline season is summer only (183 days)
    expected = (183 * 14.0 + 62 * 2.0) / (183 + 62)
    assert abs(gst1 - expected) < 1e-9
    assert gst1 < gst0


def test_no_delta_is_identity():
    t = _series(summer_c=11.0, shoulder_c=3.0, winter_c=-8.0)
    a = growing_season_temperature(t, MONTH, 2000.0, 1900.0, GAMMA)
    b = growing_season_temperature(t, MONTH, 2000.0, 1900.0, GAMMA, delta_by_month=np.zeros(12))
    assert a == b
