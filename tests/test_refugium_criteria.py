import numpy as np
import pytest

from antar.viability import criteria as C

BBOX = (43.4, 38.8, 46.7, 41.4)
SHAPE = (312, 396)


def test_lattice_index_recovers_the_row_and_column_of_a_node():
    d_lat, d_lon = (BBOX[3] - BBOX[1]) / SHAPE[0], (BBOX[2] - BBOX[0]) / SHAPE[1]
    rows, cols = np.array([0, 7, 14, 308]), np.array([0, 21, 392, 7])
    lat = BBOX[3] - (rows + 0.5) * d_lat
    lon = BBOX[0] + (cols + 0.5) * d_lon
    ii, jj = C.lattice_index(lat, lon, BBOX, SHAPE, 7)
    assert ii.tolist() == [0, 1, 2, 44] and jj.tolist() == [0, 3, 56, 1]


def test_lattice_index_refuses_nodes_that_are_not_on_the_stride():
    d_lat, d_lon = (BBOX[3] - BBOX[1]) / SHAPE[0], (BBOX[2] - BBOX[0]) / SHAPE[1]
    with pytest.raises(ValueError):
        C.lattice_index([BBOX[3] - 3.5 * d_lat], [BBOX[0] + 0.5 * d_lon], BBOX, SHAPE, 7)


def test_buffer_index_is_positive_in_a_dip_and_negative_on_a_bump():
    n = 15
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()
    e = np.zeros(ii.size)
    e[(ii == 7) & (jj == 7)] = -3.0                  # a cool, moist pocket
    e[(ii == 3) & (jj == 11)] = +3.0                 # an exposed spot
    b = C.node_buffer_index(e, ii, jj, (n, n), radius_cells=3)
    assert b[(ii == 7) & (jj == 7)][0] > 0 > b[(ii == 3) & (jj == 11)][0]
    flat = C.node_buffer_index(np.full(ii.size, 5.0) + np.random.default_rng(0).normal(0, 1.0, ii.size), ii, jj, (n, n), 3)
    assert np.all(np.isfinite(flat)) and abs(flat.mean()) < 0.5       # pure noise has no systematic buffer


def test_buffer_index_ignores_lattice_gaps():
    ii, jj = np.meshgrid(np.arange(9), np.arange(9), indexing="ij")
    keep = ~((ii > 5) & (jj > 5))                    # a corner outside the country
    ii, jj = ii[keep], jj[keep]
    e = (ii + jj).astype(float)
    b = C.node_buffer_index(e, ii, jj, (9, 9), radius_cells=2)
    assert np.all(np.isfinite(b)) and b.shape == e.shape


def test_exposure_composite_averages_standardised_variables():
    out = C.exposure_composite([np.array([10.0, 20.0]), np.array([1.0, 3.0])], [(10.0, 5.0), (2.0, 1.0)])
    assert out.tolist() == pytest.approx([(0.0 + -1.0) / 2, (2.0 + 1.0) / 2])


def test_class_fractions():
    f = C.class_fractions(np.array([[31, 31, 32], [12, 36, 36]]), {"beech": [31], "oak": [32, 33], "pine": [36]})
    assert f == pytest.approx({"beech": 2 / 6, "oak": 1 / 6, "pine": 2 / 6})


def test_group_aoa_marks_a_distant_node_as_outside_and_a_typical_one_as_inside():
    rng = np.random.default_rng(1)
    X = rng.normal(0, 1, (120, 3))
    lat, lon = rng.uniform(39, 41, 120), rng.uniform(44, 46, 120)
    r = C.group_aoa(X, lat, lon, np.array([[0.0, 0.0, 0.0], [25.0, 25.0, 25.0], [np.nan, 0.0, 0.0]]))
    assert r["status"] == "ok" and r["n_train"] == 120
    assert r["inside"].tolist() == [True, False, False] and np.isnan(r["di"][2])


def test_group_aoa_is_not_established_from_a_tiny_sample():
    X = np.random.default_rng(2).normal(size=(8, 3))
    r = C.group_aoa(X, np.linspace(39, 41, 8), np.linspace(44, 46, 8), np.zeros((2, 3)))
    assert r["status"] == "insufficient_training_sample" and not r["inside"].any() and r["threshold"] is None


def test_conjunction_needs_all_three_criteria_and_a_finite_buffer():
    a = np.array([True, True, True, False, True])
    b = np.array([1.0, -1.0, np.nan, 1.0, 0.0])
    c = np.array([True, True, True, True, True])
    assert C.robust_conjunction(a, b, c).tolist() == [True, False, False, False, False]
    assert C.robust_conjunction(a, b, np.array([False, True, True, True, True])).tolist() == [False] * 5
