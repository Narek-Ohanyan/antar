import numpy as np
import pytest

from antar.hazard import forest_panel as fp
from antar.hazard import observation


def test_subpoint_offsets_are_the_centres_of_the_thirds():
    assert fp.subpoint_offsets(30.0, 3).tolist() == [5.0, 15.0, 25.0]


def test_master_subpoints_cover_the_pixel_and_run_south_with_the_row():
    x, y = fp.master_subpoints([0, 2], [0, 1], x0=1000.0, y0=2000.0)
    assert x.shape == y.shape == (2, 9)
    assert sorted(set(x[0])) == [1005.0, 1015.0, 1025.0]
    assert sorted(set(y[0])) == [1975.0, 1985.0, 1995.0]
    assert sorted(set(x[1])) == [1035.0, 1045.0, 1055.0]          # column 1 starts 30 m east
    assert sorted(set(y[1])) == [1915.0, 1925.0, 1935.0]          # row 2 starts 60 m south


def test_sample_raster_reads_pixels_and_returns_nan_outside():
    arr = np.arange(12, dtype=float).reshape(3, 4)                 # 3 rows x 4 cols of 10 m, upper-left (1000, 2000)
    v = fp.sample_raster(arr, [1005, 1035, 999, 1005], [1995, 1975, 1995, 2005], 1000.0, 2000.0, 10.0)
    assert v[0] == arr[0, 0] and v[1] == arr[2, 3]
    assert np.isnan(v[2]) and np.isnan(v[3])                       # west of the raster, north of the raster


def test_sample_raster_keeps_integer_dtype_with_an_integer_fill():
    arr = np.array([[31, 36]], dtype=np.int16)
    v = fp.sample_raster(arr, [1005, 1015, 1025], [1995, 1995, 1995], 1000.0, 2000.0, 10.0, fill=0)
    assert v.tolist() == [31, 36, 0] and v.dtype == np.int16


def test_a_ten_metre_raster_offset_by_ten_metres_is_sampled_without_resampling():
    # 10 m map whose upper-left corner is 10 m east of the master pixel's edge: the first sub-point (x = 1005) falls outside it
    classes = np.full((3, 3), 36, dtype=np.int16)
    x, y = fp.master_subpoints([0], [0], 1000.0, 2000.0)
    cls = fp.sample_raster(classes, x[0], y[0], 1010.0, 2000.0, 10.0, fill=0)
    assert (cls == 0).sum() == 3 and (cls == 36).sum() == 6        # one of the three columns of sub-points is missing
    share, modal = fp.forest_pixel_properties(cls[None, :])
    assert share[0] == pytest.approx(6 / 9) and modal[0] == 36


def test_forest_pixel_properties_share_and_modal_class():
    sub = np.array([[36] * 9, [31, 31, 31, 12, 12, 12, 13, 13, 13], [12] * 9, [31, 31, 32, 32, 32, 36, 36, 36, 36]])
    share, modal = fp.forest_pixel_properties(sub)
    assert share.tolist() == pytest.approx([1.0, 1 / 3, 0.0, 1.0])
    assert modal.tolist() == [36, 31, 0, 36]
    assert fp.forest_pixel_properties(np.array([[31, 31, 32, 32, 12, 12, 12, 12, 12]]))[1].tolist() == [31]   # tie: lower code


def test_plantations_are_forest_but_belong_to_no_species_group():
    assert 37 in fp.CLOSED_FOREST_CLASSES and 37 not in fp.GROUP_OF_CLASS
    assert fp.GROUP_OF_CLASS[36] == "pine" and fp.GROUP_OF_CLASS[31] == "mesic_diffuse_porous_broadleaf"
    assert {fp.GROUP_OF_CLASS[c] for c in (32, 33, 34, 35)} == {"ring_porous_oak"}


def test_in_cell_lower_edge_inclusive_upper_exclusive():
    h = 15 / 3600
    assert fp.in_cell(40.0 - h, 45.0, 40.0, 45.0, h)
    assert not fp.in_cell(40.0 + h, 45.0, 40.0, 45.0, h)
    assert not fp.in_cell(40.0, 45.0 + 2 * h, 40.0, 45.0, h)


def test_choose_pixels_is_reproducible_bounded_sorted_and_a_subset():
    cand = np.arange(100, 200)
    a, b = fp.choose_pixels(cand, 25, seed=7), fp.choose_pixels(cand, 25, seed=7)
    assert a.tolist() == b.tolist() and len(a) == 25 and len(set(a)) == 25
    assert np.all(np.diff(a) > 0) and set(a) <= set(cand)
    assert fp.choose_pixels(cand, 25, seed=8).tolist() != a.tolist()
    assert fp.choose_pixels(cand[:10], 25, seed=1).tolist() == cand[:10].tolist()


def test_node_seed_depends_on_the_node_and_not_on_the_others():
    assert fp.node_seed(5, 3) == fp.node_seed(5, 3) != fp.node_seed(5, 4)


def test_mean_over_subpoints_needs_enough_valid_values():
    v = np.array([[10.0] * 9, [10.0] * 5 + [np.nan] * 4, [10.0] * 4 + [np.nan] * 5])
    out = fp.mean_over_subpoints(v, min_valid=5)
    assert out[0] == 10.0 and out[1] == 10.0 and np.isnan(out[2])


def test_trailing_z_matches_the_definition_and_leaves_undetermined_years_nan():
    x = np.array([0.50, 0.52, 0.48, 0.51, 0.49, 0.50, 0.53, 0.47, 0.50, 0.51, 0.30, 0.50])
    z = observation.trailing_z(x, 10)
    assert np.all(np.isnan(z[:10]))
    w = x[:10]
    assert z[10] == pytest.approx((0.30 - np.median(w)) / np.std(w, ddof=1))
    assert np.isfinite(z[11])
    flat = observation.trailing_z(np.full(12, 0.5), 10)
    assert np.all(np.isnan(flat))                                  # exactly zero spread: undetermined, not zero
    gap = x.copy()
    gap[3] = np.nan
    assert np.isnan(observation.trailing_z(gap, 10)[10])           # a missing value in the window


# ---- frame and extraction on synthetic rasters ---------------------------------------------------------------------------------------------------------------

X0, Y0 = 500_000.0, 4_500_000.0           # master grid corner (UTM-like)
M_PER_DEG = 1000.0                         # toy conversion: 1 degree = 1000 m in both directions, lon 0 at x = X0, lat 0 at y = Y0 - 100000


def _to_xy(lon, lat):
    return X0 + np.asarray(lon, dtype=float) * M_PER_DEG, Y0 - 100_000.0 + np.asarray(lat, dtype=float) * M_PER_DEG


def _to_lonlat(x, y):
    return (np.asarray(x, dtype=float) - X0) / M_PER_DEG, (np.asarray(y, dtype=float) - Y0 + 100_000.0) / M_PER_DEG


def _class_reader_factory(forest_fn):
    """A 10 m class map whose class at (x, y) is forest_fn(x, y), exposed through the reader interface of build_frame."""
    def reader(xmin, ymin, xmax, ymax):
        ox, oy = np.floor((xmin - 20.0) / 10) * 10, np.ceil((ymax + 20.0) / 10) * 10
        w, h = int((xmax + 20 - ox) // 10) + 2, int((oy - (ymin - 20)) // 10) + 2
        xs = ox + (np.arange(w) + 0.5) * 10
        ys = oy - (np.arange(h) + 0.5) * 10
        XX, YY = np.meshgrid(xs, ys)
        return forest_fn(XX, YY).astype(np.int16), ox, oy, 10.0
    return reader


def test_frame_keeps_only_pure_forest_pixels_inside_the_cell():
    half = 0.15                              # 150 m, so a cell holds 10 x 10 master pixels
    lat0, lon0 = 0.5, 0.5                    # node centre at x = X0 + 500, y = Y0 - 100000 + 500
    cx, cy = _to_xy(lon0, lat0)
    # forest = the western half of the cell, class 36; elsewhere 12
    reader = _class_reader_factory(lambda X, Y: np.where(X < cx, 36, 12))
    f = fp.build_frame([lat0], [lon0], reader, _to_xy, _to_lonlat, half, X0, Y0, k=1000, seed=1, min_candidates=5)
    assert f["kept"].tolist() == [True]
    # the cell holds 10 x 10 master pixels; the forest edge is at the centre line, which cuts the fifth column of pixels, so 4 of the 10 columns are pure forest
    assert f["n_candidates"][0] == len(f["row"]) == 40
    assert np.all(f["x"] + 15.0 <= cx)                            # every sub-point of a chosen pixel lies west of the edge
    assert np.all(f["forest_class"] == 36)
    inside = fp.in_cell(f["lat"], f["lon"], lat0, lon0, half)
    assert inside.all()
    # the stored row/col reproduce the pixel centre
    assert f["x"] == pytest.approx(X0 + (f["col"] + 0.5) * 30.0) and f["y"] == pytest.approx(Y0 - (f["row"] + 0.5) * 30.0)


def test_frame_drops_nodes_with_too_few_candidates_and_caps_the_draw():
    half = 0.15
    reader_all = _class_reader_factory(lambda X, Y: np.full(X.shape, 31))
    reader_none = _class_reader_factory(lambda X, Y: np.full(X.shape, 12))
    f = fp.build_frame([0.5, 3.5], [0.5, 3.5], reader_all, _to_xy, _to_lonlat, half, X0, Y0, k=20, seed=3, min_candidates=5)
    assert f["kept"].tolist() == [True, True]
    assert np.bincount(f["node"]).tolist() == [20, 20]
    g = fp.build_frame([0.5], [0.5], reader_none, _to_xy, _to_lonlat, half, X0, Y0, k=20, seed=3)
    assert g["kept"].tolist() == [False] and g["n_candidates"].tolist() == [0] and len(g["node"]) == 0


def test_frame_draw_for_a_node_does_not_depend_on_the_other_nodes():
    half = 0.15
    reader = _class_reader_factory(lambda X, Y: np.full(X.shape, 36))
    both = fp.build_frame([0.5, 3.5], [0.5, 3.5], reader, _to_xy, _to_lonlat, half, X0, Y0, k=10, seed=9)
    second_only = fp.build_frame([3.5], [3.5], reader, _to_xy, _to_lonlat, half, X0, Y0, k=10, seed=9)
    # the seed is keyed by the node's id, so a node keeps its draw when the others are removed, and a different id gives a different draw
    solo_second = fp.build_frame([3.5], [3.5], reader, _to_xy, _to_lonlat, half, X0, Y0, k=10, seed=9, node_ids=[1])
    assert both["row"][both["node"] == 1].tolist() == solo_second["row"].tolist() and solo_second["node"].tolist() == [1] * 10
    assert second_only["row"].tolist() != solo_second["row"].tolist()      # same pixels available, but node id 0 draws differently from id 1


def test_mixed_pixels_are_not_candidates():
    half = 0.15
    # one 10 m stripe of non-forest every 30 m: every master pixel then holds a non-forest sub-point (at x = 5 mod 30)
    reader = _class_reader_factory(lambda X, Y: np.where(((X - X0) % 30.0) < 10.0, 12, 36))
    f = fp.build_frame([0.5], [0.5], reader, _to_xy, _to_lonlat, half, X0, Y0, k=100)
    assert not f["kept"][0] and f["n_candidates"][0] == 0


YEARS = [2000, 2001, 2002]
TILE = 8


class FakeBackend:
    """A tiny stand-in for the Drive rasters: every value is a known function of its master-grid position."""
    x0, y0 = X0, Y0
    height_origin = (X0 - 40.0, Y0 + 40.0, 10.0)           # a 10 m grid commensurate with the master grid, starting 40 m north-west of it
    height_shape = (60, 60)

    def __init__(self):
        self.calls = []

    @staticmethod
    def kn(row, col, y):
        return 0.001 * row + 0.0001 * col + 0.01 * (y - 2000)

    def vitality(self, tr0, tc0, window):
        self.calls.append(("vitality", tr0, tc0, window))
        r, c, h, w = window
        names = [n for y in YEARS for n in (f"kndvi_{y}", f"valid_count_{y}")]
        arr = np.zeros((len(names), h, w), dtype=np.float32)
        for j, n in enumerate(names):
            y = int(n[-4:])
            for a in range(h):
                for b in range(w):
                    arr[j, a, b] = self.kn(tr0 + r + a, tc0 + c + b, y) if n.startswith("kndvi") else 100 + y - 2000
        return arr, names

    def disturbance(self, window):
        r, c, h, w = window
        names = [f"no_disturbance_{y}" for y in YEARS]
        arr = np.ones((3, h, w), dtype=np.uint8)
        arr[1] = ((np.arange(r, r + h)[:, None] + np.arange(c, c + w)[None, :]) % 2)       # 2001: a checkerboard
        return arr, names

    def terrain(self, window):
        r, c, h, w = window
        rr, cc = np.meshgrid(np.arange(r, r + h), np.arange(c, c + w), indexing="ij")
        return np.stack([1000.0 + rr + cc, 5.0 * np.ones((h, w)), np.zeros((h, w))]).astype(np.float32), ["elevation", "slope", "aspect"]

    def height(self, window):
        r, c, h, w = window
        rr, cc = np.meshgrid(np.arange(r, r + h), np.arange(c, c + w), indexing="ij")
        return (rr + 0.01 * cc).astype(np.float32)


def test_extract_node_reads_each_pixel_across_a_tile_boundary_and_maps_the_bands():
    rows = np.array([5, 6, 9, 10])             # the tile edge is at row 8 and column 8
    cols = np.array([6, 9, 7, 12])
    bk = FakeBackend()
    x = X0 + (cols + 0.5) * 30.0
    y = Y0 - (rows + 0.5) * 30.0
    out = fp.extract_node(bk, rows, cols, x, y, YEARS, TILE)
    for i, (r, c) in enumerate(zip(rows, cols)):
        for j, yr in enumerate(YEARS):
            assert out["kndvi"][i, j] == pytest.approx(bk.kn(r, c, yr), abs=1e-5)
            assert out["valid_count"][i, j] == 100 + yr - 2000
        assert out["elevation"][i] == 1000.0 + r + c
        assert out["no_disturbance"][i, 0] and out["no_disturbance"][i, 2]
        assert out["no_disturbance"][i, 1] == bool((r + c) % 2)
    tiles = {(t[1], t[2]) for t in bk.calls}
    assert tiles == {(0, 0), (0, 8), (8, 0), (8, 8)}             # one read per tile that holds a pixel


def test_extract_node_height_is_the_mean_of_the_nine_sub_points():
    bk = FakeBackend()
    rows, cols = np.array([2]), np.array([3])
    out = fp.extract_node(bk, rows, cols, np.array([0.0]), np.array([0.0]), YEARS, TILE)
    # master pixel (2, 3) spans x in [X0+90, X0+120], y in [Y0-90, Y0-60]; sub-point k sits at 5, 15, 25 m from the west/north edge.
    # The height grid starts at (X0-40, Y0+40): column = floor((x - (X0-40)) / 10), row = floor(((Y0+40) - y) / 10)
    expect = []
    for dy in (5, 15, 25):
        for dx in (5, 15, 25):
            col = int(np.floor((90 + dx + 40) / 10))
            row = int(np.floor((40 + 60 + dy) / 10))
            expect.append(row + 0.01 * col)
    assert out["height_mean"][0] == pytest.approx(np.mean(expect), abs=1e-5)
    assert out["height_n_subpoints"][0] == 9


def test_extract_node_height_is_nan_when_the_pixel_is_off_the_height_grid():
    bk = FakeBackend()
    out = fp.extract_node(bk, np.array([200]), np.array([200]), np.array([0.0]), np.array([0.0]), YEARS, 512)
    assert np.isnan(out["height_mean"][0]) and out["height_n_subpoints"][0] == 0


def test_band_index_names_the_missing_band():
    with pytest.raises(KeyError, match="kndvi_2005"):
        fp.band_index(["kndvi_2000"], "kndvi_{year}", [2000, 2005])


def test_split_by_tile_and_window_helpers():
    g = fp.split_by_tile([1, 5, 9, 20], [1, 12, 3, 20], 8)
    assert {k: v.tolist() for k, v in g.items()} == {(0, 0): [0], (0, 8): [1], (8, 0): [2], (16, 16): [3]}
    assert fp.window_of([3, 5], [10, 7]) == (3, 7, 3, 4)
    win = np.arange(2 * 3 * 4).reshape(2, 3, 4)
    assert fp.values_in_window(win, 10, 20, [10, 12], [21, 23]).tolist() == [[1, 13], [11, 23]]


def _toy_frame():
    # two nodes (ids 3 and 7), three pixels each, in a 16 x 16 master-pixel toy world with 8-pixel tiles
    node = np.array([3, 3, 3, 7, 7, 7])
    row = np.array([1, 2, 2, 12, 13, 12])
    col = np.array([1, 1, 3, 12, 12, 14])
    return {"node": node, "row": row, "col": col, "x": X0 + (col + 0.5) * 30.0, "y": Y0 - (row + 0.5) * 30.0,
            "lat": np.zeros(6), "lon": np.zeros(6), "forest_class": np.full(6, 36)}


class FailingFor(FakeBackend):
    def __init__(self, bad_rows):
        super().__init__()
        self.bad_rows = bad_rows

    def disturbance(self, window):
        if window[0] in self.bad_rows:
            raise OSError("Drive said no")
        return super().disturbance(window)


def test_extract_all_lists_a_failed_node_and_keeps_the_others():
    frame = _toy_frame()
    done, failed = fp.extract_all(frame, FailingFor({12}), YEARS, TILE)         # node 7's window starts at row 12
    assert sorted(done) == [3] and [n for n, _ in failed] == [7] and "Drive said no" in failed[0][1]
    out = fp.assemble(frame, done)
    assert out["node"].tolist() == [3, 3, 3] and out["kndvi"].shape == (3, 3)


def test_extract_all_resumes_without_rereading_finished_nodes():
    frame = _toy_frame()
    bk = FakeBackend()
    first, _ = fp.extract_all(frame, FailingFor({12}), YEARS, TILE)
    seen = []
    done, failed = fp.extract_all(frame, bk, YEARS, TILE, done=first, on_node=lambda d, n: seen.append(n))
    assert failed == [] and sorted(done) == [3, 7] and seen == [7]
    assert all(c[1:3] != (0, 0) for c in bk.calls if c[0] == "vitality")      # node 3 sits in tile (0, 0) and was not read again
    out = fp.assemble(frame, done)
    assert out["node"].tolist() == frame["node"].tolist()
    assert out["kndvi"][3, 0] == pytest.approx(FakeBackend.kn(12, 12, 2000), abs=1e-5)       # node 7's first pixel, read in the second run
    assert out["kndvi"][0, 0] == pytest.approx(FakeBackend.kn(1, 1, 2000), abs=1e-5)         # node 3's first pixel, carried over from the first run


def test_other_exceptions_are_not_swallowed():
    class Broken(FakeBackend):
        def terrain(self, window):
            raise KeyError("elevation")
    with pytest.raises(KeyError):
        fp.extract_all(_toy_frame(), Broken(), YEARS, TILE)


def test_node_order_groups_nodes_by_export_block():
    frame = {"node": np.array([1, 2, 3]), "row": np.array([600, 10, 20]), "col": np.array([10, 700, 20])}
    assert fp.node_order(frame, 256) == [3, 2, 1]           # block (0,0) = node 3, then (0,2) = node 2, then (2,0) = node 1


def test_assemble_requires_pixels_grouped_by_node():
    frame = _toy_frame()
    frame["node"] = frame["node"][::-1]
    with pytest.raises(ValueError):
        fp.assemble(frame, {})


def test_a_layer_that_reads_disturbed_everywhere_is_rejected():
    years = list(range(2000, 2025))
    broken = np.zeros((40, 25), dtype=bool)
    broken[:, 0] = True                                   # the first export: 1 in 2000 (no loss data), 0 in every later year
    msg = fp.disturbance_layer_problem(broken, years)
    assert msg is not None and "2001" in msg and "not a plausible record" in msg


def test_a_layer_with_rare_disturbance_passes():
    rng = np.random.default_rng(0)
    years = list(range(2000, 2025))
    ok = rng.random((400, 25)) > 0.01                     # 1 % of pixel-years flagged
    assert fp.disturbance_layer_problem(ok, years) is None
    assert fp.disturbance_layer_problem(np.zeros((0, 25), dtype=bool), years) is not None


def test_event_rate_interval_matches_known_exact_values():
    lo, hi = observation.event_rate_interval(5, 100)
    assert lo == pytest.approx(0.01643, abs=2e-5) and hi == pytest.approx(0.11283, abs=2e-5)       # Clopper-Pearson 5/100
    lo0, hi0 = observation.event_rate_interval(0, 1000)
    assert lo0 == 0.0 and hi0 == pytest.approx(1 - 0.025 ** (1 / 1000))
    with pytest.raises(ValueError):
        observation.event_rate_interval(3, 2)


def test_zero_event_upper_bound_is_about_three_over_n():
    assert observation.zero_event_upper_bound(778) == pytest.approx(3.0 / 778, rel=0.02)
    assert observation.zero_event_upper_bound(778) == pytest.approx(1 - 0.05 ** (1 / 778))
    with pytest.raises(ValueError):
        observation.zero_event_upper_bound(0)


def _pixels(n=4, n_years=25):
    rng = np.random.default_rng(0)
    kn = 0.6 + 0.01 * rng.normal(size=(n, n_years))
    kn[0, 15:] = 0.1                       # pixel 0 collapses in 2015 and stays down: a dieback onset
    kn[1, 15:] = 0.1                       # pixel 1 does the same, but a harvest is flagged in 2015
    nd = np.ones((n, n_years), dtype=bool)
    nd[1, 15] = False
    kn[2, 5] = np.nan
    return {"node": np.array([10, 10, 20, 20]), "lat": np.array([40.0, 40.0, 41.0, 41.0]), "lon": np.array([44.0, 44.0, 45.0, 45.0]), "forest_class": np.array([31, 31, 33, 36]),
            "kndvi": kn, "valid_count": np.full((n, n_years), 20.0), "no_disturbance": nd, "height_mean": np.array([15.0, 16.0, 12.0, np.nan]), "elevation": np.array([1200.0, 1250.0, 900.0, 1500.0])}


def test_person_years_flag_the_event_and_the_harvest_suppressed_candidate_and_drop_unusable_rows():
    import pandas as pd
    years = list(range(2000, 2025))
    pix = _pixels()
    clim = pd.DataFrame([{"node": nd, "year": y, "cwd_mm": 100.0 + y - 2000 + nd} for nd in (10, 20) for y in range(2010, 2020)])
    df = fp.build_person_years(pix, years, list(range(2010, 2020)), clim, ["cwd_mm"])
    assert set(df.columns) >= {"pixel", "node", "year", "event", "candidate", "z", "kndvi", "height_m", "cwd_mm", "forest_class"}
    onset = df[(df.year == 2015)].set_index("pixel")
    assert bool(onset.loc[0, "event"]) and bool(onset.loc[0, "candidate"])
    assert not bool(onset.loc[1, "event"]) and bool(onset.loc[1, "candidate"])         # the harvest flag removes the event, the candidate stays
    assert df.event.sum() == 1 and df.candidate.sum() == 2                              # a persistent collapse counts once, at its first onset
    p0 = df[df.pixel == 0].set_index("year")
    assert p0.loc[2015, "at_risk"] and not p0.loc[2016, "at_risk"] and p0.loc[2010, "at_risk"]
    assert p0.loc[2016:, "event"].sum() == 0 and not p0.loc[2019, "candidate_at_risk"]
    p1 = df[df.pixel == 1].set_index("year")                                           # harvest-flagged: no strict event, so it stays at risk, but its candidate onset is counted
    assert p1["at_risk"].all() and p1.loc[2015, "candidate"] and not p1.loc[2016, "candidate_at_risk"]
    assert not p1.loc[2014, "disturbed_recent"] and p1.loc[2015, "disturbed_recent"] and p1.loc[2019, "disturbed_recent"]       # the flag stays in the baseline window for ten years
    assert not df[df.pixel == 0]["disturbed_recent"].any()
    assert df.loc[df.pixel == 0, "cwd_mm"].iloc[0] == 100.0 + 10 + 10 and df.loc[df.pixel == 2, "cwd_mm"].iloc[0] == 100.0 + 10 + 20
    assert np.isnan(df.loc[df.pixel == 3, "height_m"]).all()
    assert len(df) == 4 * 10


def test_person_years_drop_missing_vitality_and_missing_climate():
    import pandas as pd
    years = list(range(2000, 2025))
    pix = _pixels()
    pix["kndvi"][2, 12] = np.nan                                                    # 2012 missing for pixel 2
    clim = pd.DataFrame([{"node": nd, "year": y, "cwd_mm": 1.0} for nd in (10, 20) for y in range(2010, 2020) if not (nd == 20 and y == 2018)])
    df = fp.build_person_years(pix, years, list(range(2010, 2020)), clim, ["cwd_mm"])
    assert df[(df.pixel == 2) & (df.year == 2012)].empty
    assert df[(df.node == 20) & (df.year == 2018)].empty
    assert len(df) == 4 * 10 - 1 - 2


def test_subset_static_indexes_arrays_and_the_nested_soil_dict():
    static = {"lats": np.arange(5.0), "soil": {"theta_fc": np.arange(5.0) * 2, "psi": np.arange(5.0) * 3}, "token": "secret", "valid_soil": np.array([1, 1, 0, 1, 1], bool)}
    s = fp.subset_static(static, [4, 1])
    assert s["lats"].tolist() == [4.0, 1.0] and s["soil"]["theta_fc"].tolist() == [8.0, 2.0] and s["soil"]["psi"].tolist() == [12.0, 3.0]
    assert s["token"] is None and s["valid_soil"].tolist() == [True, True]


def test_mirror_count_is_symmetric_for_symmetric_noise_and_separates_real_declines():
    rng = np.random.default_rng(11)
    years, panel_years = list(range(2000, 2025)), list(range(2010, 2020))
    noise = 0.6 + 0.01 * rng.normal(size=(1500, 25))
    down, up = fp.count_onset_pixels(noise, years, panel_years), fp.mirror_onset_count(noise, years, panel_years)
    assert down > 0 and up > 0 and abs(down - up) < 4 * np.sqrt(down + up)             # symmetric noise: as many persistent rises as declines
    collapsed = noise.copy()
    collapsed[:120, 14:] = 0.1                                                          # 8 % of the pixels really die back in 2014
    d2, u2 = fp.count_onset_pixels(collapsed, years, panel_years), fp.mirror_onset_count(collapsed, years, panel_years)
    assert d2 - down > 100 and abs(u2 - up) < 4 * np.sqrt(up) + 5                       # declines rise by about the real events, rises do not


def test_mirror_of_a_series_with_one_rise_and_no_decline():
    x = np.full(25, 0.5) + 0.002 * np.sin(np.arange(25))
    x[14:] = 0.9
    assert fp.count_onset_pixels(x[None], list(range(2000, 2025)), list(range(2010, 2020))) == 0
    assert fp.mirror_onset_count(x[None], list(range(2000, 2025)), list(range(2010, 2020))) == 1


def test_removing_the_year_effect_takes_out_a_step_every_pixel_shares_but_not_a_local_collapse():
    rng = np.random.default_rng(21)
    kn = 0.5 + 0.01 * rng.normal(size=(300, 25))
    kn[:, 13:] += 0.15                                   # a sensor step in 2013 that every pixel shares
    kn[:20, 18:] -= 0.35                                 # a real, local collapse in 2018
    adj = fp.remove_year_effect(kn)
    assert abs(np.median(adj[:, 13:].mean(axis=1)) - np.median(adj[:, :13].mean(axis=1))) < 0.01          # the shared step is gone
    assert (adj[:20, 18:].mean() - adj[:20, :13].mean()) < -0.25                                          # the local collapse is still there
    assert fp.count_onset_pixels(adj, list(range(2000, 2025)), list(range(2010, 2020))) >= 18
    nan = kn.copy()
    nan[0, 5] = np.nan
    assert np.isnan(fp.remove_year_effect(nan)[0, 5]) and np.isfinite(fp.remove_year_effect(nan)[1, 5])


def test_person_years_carry_the_year_adjusted_columns():
    import pandas as pd
    years = list(range(2000, 2025))
    pix = _pixels()
    clim = pd.DataFrame([{"node": nd, "year": y, "cwd_mm": 100.0} for nd in (10, 20) for y in range(2010, 2020)])
    df = fp.build_person_years(pix, years, list(range(2010, 2020)), clim, ["cwd_mm"])
    assert {"z_adj", "event_adj", "candidate_adj", "at_risk_adj"} <= set(df.columns)
    assert df[df.pixel == 0].set_index("year").loc[2015, "event_adj"] in (True, False)
    assert np.isfinite(df["z_adj"]).sum() > 0
