"""Late-frost days count the frosts after budburst and up to the warmest day of the year, and only those."""
import numpy as np
import pytest

from antar.climate import indices
from antar.climate.forcing import topoclimate_forcing


def _year(peak_day=200, amplitude=12.0, mean=8.0, n=365):
    d = np.arange(n)
    return mean + amplitude * np.cos(2 * np.pi * (d - peak_day) / n)


def test_warmest_day_is_the_peak_of_a_smooth_year():
    assert indices.warmest_day_index(_year(peak_day=200)) == 200
    assert indices.warmest_day_index(_year(peak_day=185)) == 185


def test_one_hot_day_does_not_move_the_warmest_day():
    t = _year(peak_day=200)
    t[60] += 40.0                                    # a single freak afternoon in early March
    assert indices.warmest_day_index(t, window_days=1) == 60
    assert indices.warmest_day_index(t, window_days=15) == 200


def test_ends_of_the_record_are_not_penalised():
    t = np.linspace(0.0, 10.0, 100)                  # warmest at the last day
    assert indices.warmest_day_index(t) == 99
    assert indices.warmest_day_index(t[::-1]) == 0


def test_window_must_be_odd():
    with pytest.raises(ValueError):
        indices.warmest_day_index(_year(), window_days=14)


def test_ties_go_to_the_earlier_day():
    assert indices.warmest_day_index(np.full(30, 5.0), window_days=1) == 0


def test_only_spring_frosts_are_counted_with_an_end_day():
    tmin = np.full(365, 3.0)
    tmin[[100, 110, 120]] = -4.0                     # spring frosts
    tmin[[300, 310]] = -6.0                          # autumn frosts
    gdd = np.where(np.arange(365) >= 90, 250.0, 0.0)   # budburst reached on day 90
    assert indices.late_frost_days(tmin, gdd, 200.0) == 5               # whole record: the old definition
    assert indices.late_frost_days(tmin, gdd, 200.0, end_day=200) == 3  # budburst to the warmest day


def test_frost_before_budburst_and_milder_than_critical_do_not_count():
    tmin = np.full(365, 3.0)
    tmin[[10, 50]] = -9.0                            # deep frost, but the buds are still closed
    tmin[[100, 101]] = -1.9                          # not below -2
    tmin[102] = -2.0                                 # equal to the threshold: not below
    gdd = np.where(np.arange(365) >= 90, 250.0, 0.0)
    assert indices.late_frost_days(tmin, gdd, 200.0, end_day=200) == 0


def test_end_day_is_inclusive():
    tmin = np.full(365, 3.0)
    tmin[150] = -4.0
    gdd = np.full(365, 300.0)
    assert indices.late_frost_days(tmin, gdd, 200.0, end_day=150) == 1
    assert indices.late_frost_days(tmin, gdd, 200.0, end_day=149) == 0


def test_budburst_after_the_warmest_day_gives_zero():
    tmin = np.full(365, 3.0)
    tmin[300] = -8.0
    gdd = np.where(np.arange(365) >= 280, 250.0, 0.0)  # a very cold site: buds open only after the summer peak
    assert indices.late_frost_days(tmin, gdd, 200.0, end_day=200) == 0


def test_forcing_counts_to_the_warmest_day_not_the_year_end():
    n = 365
    doy = np.arange(1, n + 1)
    month = np.minimum((doy - 1) // 31 + 1, 12)
    t_mean = _year(peak_day=200, amplitude=11.0, mean=6.0)
    t_min = t_mean - 6.0
    kw = dict(
        doy=doy, month=month, t_mean_ref_c=t_mean, t_max_ref_c=t_mean + 5.0, t_min_ref_c=t_min, p_ref_mm=np.full(n, 1.0),
        ea_ref_kpa=np.full(n, 0.8), u2_m_s=np.full(n, 2.0), rn_mj_m2=np.full(n, 8.0), z_cell_m=1500.0, z_ref_m=1500.0, lat_deg=40.0,
        slope_deg=0.0, aspect_deg=0.0, gamma_k_per_m=np.full(12, -0.006), precip_gradient_per_m=np.full(12, 0.0), w_max_mm=100.0,
        theta_sat=0.45, psi_sat_mpa=-0.002, b_clapp_hornberger=5.0, theta_fc=0.3, theta_lim=0.1, gdd_budburst=200.0,
    )
    cell = topoclimate_forcing(**kw)
    end = indices.warmest_day_index(cell.t_mean_c)
    assert cell.late_frost_days == indices.late_frost_days(cell.t_min_c, cell.gdd_cumulative, 200.0, end_day=end)
    whole_year = indices.late_frost_days(cell.t_min_c, cell.gdd_cumulative, 200.0)
    assert cell.late_frost_days < whole_year          # this year has autumn frosts that the old definition counted


def _cell_inputs(seed, elevation_offset=300.0, concavity=0.4):
    rng = np.random.default_rng(seed)
    n = 365
    doy = np.arange(1, n + 1)
    month = np.minimum((doy - 1) // 31 + 1, 12)
    t_mean = _year(peak_day=int(rng.integers(180, 215)), amplitude=float(rng.uniform(9, 13)), mean=float(rng.uniform(3, 9))) + rng.normal(0, 2.0, n)
    t_min = t_mean - rng.uniform(3, 8, n)
    return dict(doy=doy, month=month, t_mean_ref_c=t_mean, t_max_ref_c=t_mean + 5.0, t_min_ref_c=t_min, p_ref_mm=np.full(n, 1.0),
                ea_ref_kpa=np.full(n, 0.8), u2_m_s=np.full(n, 2.0), rn_mj_m2=np.full(n, 8.0), z_cell_m=1500.0 + elevation_offset,
                z_ref_m=1500.0, lat_deg=40.0, slope_deg=0.0, aspect_deg=0.0, gamma_k_per_m=np.full(12, -0.0062),
                precip_gradient_per_m=np.full(12, 0.0), w_max_mm=100.0, theta_sat=0.45, psi_sat_mpa=-0.002, b_clapp_hornberger=5.0,
                theta_fc=0.3, theta_lim=0.1, gdd_budburst=200.0, concavity_index=concavity, calm_clear_night_frac=0.3)


@pytest.mark.parametrize("seed", range(6))
def test_reference_shortcut_equals_the_full_pipeline(seed):
    from antar.climate.forcing import late_frost_days_from_reference
    kw = _cell_inputs(seed, elevation_offset=[-200, 0, 300, 700, 1000, 450][seed], concavity=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0][seed])
    cell = topoclimate_forcing(**kw)
    short = dict(t_mean_ref_c=kw["t_mean_ref_c"], t_min_ref_c=kw["t_min_ref_c"], month=kw["month"], z_cell_m=kw["z_cell_m"], z_ref_m=kw["z_ref_m"],
                 gamma_k_per_m=kw["gamma_k_per_m"], concavity_index=kw["concavity_index"], calm_clear_night_frac=kw["calm_clear_night_frac"],
                 gdd_budburst=kw["gdd_budburst"])
    assert late_frost_days_from_reference(**short) == cell.late_frost_days
    whole = indices.late_frost_days(cell.t_min_c, cell.gdd_cumulative, 200.0)
    assert late_frost_days_from_reference(**short, spring_only=False) == whole
    assert cell.late_frost_days <= whole
