"""Which 30 m pixels MNEME's panel follows: forest pixels, not grid nodes.

A dieback rule applied to the vitality of a pixel that is not forest (grassland, cropland, bare ground) says nothing about tree dieback. The panel therefore follows forest pixels
only. Each grid node represents one climate cell (the 30-arcsecond CHELSA pixel around it); the forest pixels of that cell are drawn at random, and they share the node's climate and
differ in their own vitality series, canopy height and elevation.

A master-grid pixel (30 m, UTM zone 38N) counts as forest when every one of the 3 x 3 sub-points at 10 m inside it falls on a closed-forest class of the national
Ecosystem Map. The Ecosystem Map and the height model are 10 m rasters whose pixel edges are offset from the master grid by a whole number of 10 m, so sub-point sampling needs no resampling.
"""
from __future__ import annotations

import numpy as np

# Ecosystem Map classes counted as closed forest: 31 Fagus orientalis and other deciduous, 32-35 Quercus, 36 Pinus kochiana, 37 closed forest plantations.
CLOSED_FOREST_CLASSES = (31, 32, 33, 34, 35, 36, 37)
GROUP_OF_CLASS = {31: "mesic_diffuse_porous_broadleaf", 32: "ring_porous_oak", 33: "ring_porous_oak", 34: "ring_porous_oak", 35: "ring_porous_oak",
                  36: "pine"}                  # 37 (plantations) has no species group


def subpoint_offsets(master_res_m: float = 30.0, n: int = 3) -> np.ndarray:
    """Offsets (m) of the n sub-point centres from a master pixel's edge: (k + 1/2) * res / n."""
    return (np.arange(n) + 0.5) * master_res_m / n


def master_subpoints(rows, cols, x0: float, y0: float, master_res_m: float = 30.0, n: int = 3):
    """Map coordinates of the n x n sub-points of master pixels (``rows``, ``cols``): arrays x, y of shape (len(rows), n * n).

    The master grid has its upper-left corner at (x0, y0); x grows with the column, y falls with the row.
    """
    rows = np.atleast_1d(np.asarray(rows, dtype=float))
    cols = np.atleast_1d(np.asarray(cols, dtype=float))
    off = subpoint_offsets(master_res_m, n)
    ox, oy = np.meshgrid(off, off)                      # (n, n): ox varies along columns, oy along rows
    x = x0 + cols[:, None] * master_res_m + ox.ravel()[None, :]
    y = y0 - rows[:, None] * master_res_m - oy.ravel()[None, :]
    return x, y


def sample_raster(arr, x, y, raster_x0: float, raster_y0: float, res_m: float, fill=np.nan):
    """Value of a north-up raster (upper-left corner raster_x0, raster_y0; square pixels of res_m) at map coordinates x, y; ``fill`` outside the array."""
    arr = np.asarray(arr)
    col = np.floor((np.asarray(x, dtype=float) - raster_x0) / res_m).astype(int)
    row = np.floor((raster_y0 - np.asarray(y, dtype=float)) / res_m).astype(int)
    ok = (row >= 0) & (row < arr.shape[0]) & (col >= 0) & (col < arr.shape[1])
    out = np.full(col.shape, fill, dtype=float if np.isnan(fill) else arr.dtype)
    out[ok] = arr[row[ok], col[ok]]
    return out


def forest_pixel_properties(class_sub, forest_classes=CLOSED_FOREST_CLASSES):
    """For each pixel (row of ``class_sub``, its sub-point class codes): the share of sub-points in closed forest and the most frequent forest class (0 when none)."""
    class_sub = np.asarray(class_sub)
    is_forest = np.isin(class_sub, forest_classes)
    share = is_forest.mean(axis=1)
    modal = np.zeros(class_sub.shape[0], dtype=int)
    for k, row in enumerate(class_sub):
        f = row[is_forest[k]]
        if f.size:
            vals, counts = np.unique(f, return_counts=True)
            modal[k] = int(vals[np.argmax(counts)])       # ties go to the lower class code
    return share, modal


def in_cell(lat, lon, lat0, lon0, half_deg: float):
    """True where (lat, lon) lies inside the square of half-width ``half_deg`` around (lat0, lon0); lower edges inclusive, upper edges exclusive."""
    lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    return (lat >= lat0 - half_deg) & (lat < lat0 + half_deg) & (lon >= lon0 - half_deg) & (lon < lon0 + half_deg)


def choose_pixels(candidate_index, k: int, seed: int):
    """Up to ``k`` of the candidate indices, drawn without replacement and returned in ascending order; the draw depends on ``seed`` only."""
    candidate_index = np.asarray(candidate_index)
    if candidate_index.size <= k:
        return np.sort(candidate_index)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(candidate_index, size=k, replace=False))


def mean_over_subpoints(values_sub, min_valid: int = 5):
    """Mean of a pixel's sub-point values, NaN when fewer than ``min_valid`` are finite."""
    v = np.asarray(values_sub, dtype=float)
    ok = np.isfinite(v)
    n = ok.sum(axis=1)
    s = np.where(ok, v, 0.0).sum(axis=1)
    return np.where(n >= min_valid, s / np.maximum(n, 1), np.nan)


def node_seed(base_seed: int, node: int) -> int:
    """A seed for one node's draw: reproducible, and unaffected by which other nodes are in the run."""
    return int(base_seed) * 1_000_003 + int(node)


# ---- frame: which pixels of a node's cell are followed ------------------------------------------------------------------------------------------------------

def build_frame(node_lat, node_lon, class_reader, to_xy, to_lonlat, half_deg: float, x0: float, y0: float, master_res_m: float = 30.0,
                k: int = 25, seed: int = 0, min_candidates: int = 5, forest_classes=CLOSED_FOREST_CLASSES, node_ids=None):
    """Choose the forest pixels that stand for each node.

    A master pixel is a candidate for a node when its centre lies in the node's climate cell (``half_deg`` around the node) and all nine of its sub-points are on a closed-forest
    class. Nodes with fewer than ``min_candidates`` candidates are left out. Up to ``k`` candidates per node are drawn at random (:func:`choose_pixels`, seeded per node).

    ``class_reader(xmin, ymin, xmax, ymax)`` returns ``(array, raster_x0, raster_y0, res)``: the Ecosystem Map window over that rectangle (class 0 where there is no data);
    ``to_xy(lon, lat)`` and ``to_lonlat(x, y)`` convert between longitude/latitude and map coordinates.

    ``node_ids`` names the nodes (default 0..n-1); the draw of a node depends on its id only, so it does not change when other nodes are added or removed.
    Returns arrays over the chosen pixels (node = id, row, col, x, y, lat, lon, forest_class) and, aligned with the input nodes, n_candidates and kept.
    """
    node_lat, node_lon = np.asarray(node_lat, dtype=float), np.asarray(node_lon, dtype=float)
    n_nodes = node_lat.size
    ids = np.arange(n_nodes) if node_ids is None else np.asarray(node_ids, dtype=int)
    n_candidates = np.zeros(n_nodes, dtype=int)
    out = {k_: [] for k_ in ("node", "row", "col", "x", "y", "lat", "lon", "forest_class")}
    for i in range(n_nodes):
        lon_c = node_lon[i] + half_deg * np.array([-1, 1, -1, 1])
        lat_c = node_lat[i] + half_deg * np.array([-1, -1, 1, 1])
        cx, cy = to_xy(lon_c, lat_c)
        c_lo, c_hi = int(np.floor((cx.min() - x0) / master_res_m)), int(np.floor((cx.max() - x0) / master_res_m))
        r_lo, r_hi = int(np.floor((y0 - cy.max()) / master_res_m)), int(np.floor((y0 - cy.min()) / master_res_m))
        rr, cc = np.meshgrid(np.arange(r_lo, r_hi + 1), np.arange(c_lo, c_hi + 1), indexing="ij")
        rr, cc = rr.ravel(), cc.ravel()
        px = x0 + (cc + 0.5) * master_res_m
        py = y0 - (rr + 0.5) * master_res_m
        lon_p, lat_p = to_lonlat(px, py)
        inside = in_cell(lat_p, lon_p, node_lat[i], node_lon[i], half_deg)
        rr, cc, px, py, lon_p, lat_p = rr[inside], cc[inside], px[inside], py[inside], lon_p[inside], lat_p[inside]
        if rr.size == 0:
            continue
        sx, sy = master_subpoints(rr, cc, x0, y0, master_res_m)
        arr, rx0, ry0, res = class_reader(sx.min(), sy.min(), sx.max(), sy.max())
        cls = sample_raster(arr, sx, sy, rx0, ry0, res, fill=0)
        share, modal = forest_pixel_properties(cls, forest_classes)
        cand = np.flatnonzero(share >= 1.0)
        n_candidates[i] = cand.size
        if cand.size < min_candidates:
            continue
        pick = choose_pixels(cand, k, node_seed(seed, int(ids[i])))
        out["node"].append(np.full(pick.size, int(ids[i])))
        for key, val in (("row", rr), ("col", cc), ("x", px), ("y", py), ("lat", lat_p), ("lon", lon_p), ("forest_class", modal)):
            out[key].append(val[pick])
    frame = {key: (np.concatenate(v) if v else np.zeros(0)) for key, v in out.items()}
    for key in ("node", "row", "col", "forest_class"):
        frame[key] = frame[key].astype(int)
    frame["n_candidates"] = n_candidates
    frame["kept"] = np.isin(ids, np.unique(frame["node"]))
    return frame


# ---- extraction: values of the chosen pixels from the exported rasters ----------------------------------------------------------------------------------------

def split_by_tile(rows, cols, tile_px: int):
    """Indices of master-grid pixels grouped by the export tile that holds them: {(tile_row0, tile_col0): indices}."""
    rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
    tr, tc = (rows // tile_px) * tile_px, (cols // tile_px) * tile_px
    groups = {}
    for key in sorted(set(zip(tr.tolist(), tc.tolist()))):
        groups[key] = np.flatnonzero((tr == key[0]) & (tc == key[1]))
    return groups


def window_of(rows, cols):
    """(row_off, col_off, height, width) of the smallest window holding every pixel."""
    rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
    return int(rows.min()), int(cols.min()), int(rows.max() - rows.min() + 1), int(cols.max() - cols.min() + 1)


def values_in_window(win, row_off, col_off, rows, cols):
    """``win`` has shape (bands, h, w) and starts at (row_off, col_off); returns (len(rows), bands) with the value of every band at each pixel."""
    win = np.asarray(win)
    return win[:, np.asarray(rows) - row_off, np.asarray(cols) - col_off].T


def band_index(names, template: str, years):
    """Positions of ``template.format(year=y)`` in ``names`` for every year; raises when a band is missing."""
    pos = {n: j for j, n in enumerate(names)}
    try:
        return [pos[template.format(year=y)] for y in years]
    except KeyError as e:
        raise KeyError(f"band {e} not found in the raster") from None


def extract_node(backend, rows, cols, x, y, years, tile_px: int, master_res_m: float = 30.0):
    """Values for one node's chosen pixels, read window by window from ``backend``.

    ``backend`` provides ``vitality(tile_row0, tile_col0, window)`` -> (array (bands, h, w) in tile-local coordinates, band names), ``disturbance(window)`` and ``terrain(window)``
    (master-grid windows, same return shape), and ``height(window)`` -> (h, w) array on the 10 m height grid with ``height_origin`` = (x, y, res) and ``height_shape``.
    A window is (row_off, col_off, height, width).
    Returns arrays kndvi, valid_count, no_disturbance (n, len(years)), elevation (n,), height_mean (n,) and the number of height sub-points used.
    """
    rows, cols = np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)
    n, ny = rows.size, len(years)
    kndvi = np.full((n, ny), np.nan)
    valid = np.full((n, ny), np.nan)
    for (tr0, tc0), idx in split_by_tile(rows, cols, tile_px).items():
        r_off, c_off, h, w = window_of(rows[idx] - tr0, cols[idx] - tc0)
        arr, names = backend.vitality(tr0, tc0, (r_off, c_off, h, w))
        vals = values_in_window(arr, r_off, c_off, rows[idx] - tr0, cols[idx] - tc0)
        kndvi[idx] = vals[:, band_index(names, "kndvi_{year}", years)]
        valid[idx] = vals[:, band_index(names, "valid_count_{year}", years)]
    win = window_of(rows, cols)
    arr, names = backend.disturbance(win)
    nodist = values_in_window(arr, win[0], win[1], rows, cols)[:, band_index(names, "no_disturbance_{year}", years)] > 0
    arr, names = backend.terrain(win)
    elevation = values_in_window(arr, win[0], win[1], rows, cols)[:, list(names).index("elevation")]

    sx, sy = master_subpoints(rows, cols, backend.x0, backend.y0, master_res_m)
    hx, hy, hres = backend.height_origin
    col0 = max(int(np.floor((sx.min() - hx) / hres)), 0)
    col1 = min(int(np.floor((sx.max() - hx) / hres)), backend.height_shape[1] - 1)
    row0 = max(int(np.floor((hy - sy.max()) / hres)), 0)
    row1 = min(int(np.floor((hy - sy.min()) / hres)), backend.height_shape[0] - 1)
    if col1 >= col0 and row1 >= row0:
        harr = backend.height((row0, col0, row1 - row0 + 1, col1 - col0 + 1))
        sub = sample_raster(harr, sx, sy, hx + col0 * hres, hy - row0 * hres, hres)
    else:
        sub = np.full(sx.shape, np.nan)
    return {"kndvi": kndvi, "valid_count": valid, "no_disturbance": nodist, "elevation": np.asarray(elevation, dtype=float),
            "height_mean": mean_over_subpoints(sub), "height_n_subpoints": np.isfinite(sub).sum(axis=1)}


def node_order(frame, block_px: int = 256):
    """Nodes sorted by the export block (block_px x block_px) of their first pixel, so consecutive nodes reuse the blocks already fetched."""
    nodes = np.unique(frame["node"])
    first = np.array([np.flatnonzero(frame["node"] == n)[0] for n in nodes])
    key = np.column_stack([frame["row"][first] // block_px, frame["col"][first] // block_px])
    return [int(nodes[i]) for i in np.lexsort((key[:, 1], key[:, 0]))]


def extract_all(frame, backend, years, tile_px: int, master_res_m: float = 30.0, done=None, retriable=(OSError,), on_node=None):
    """Extract every node of ``frame`` through ``backend``; returns ``(done, failed)``.

    ``done`` maps node -> :func:`extract_node` result and may hold nodes finished in an earlier run (they are not read again). A node whose read raises one of ``retriable`` after the
    backend's own retries is listed in ``failed`` and skipped; any other exception propagates. ``on_node(done, node)`` is called after each newly finished node.
    """
    done = {} if done is None else dict(done)
    failed = []
    for node in node_order(frame):
        if node in done:
            continue
        sel = np.flatnonzero(frame["node"] == node)
        try:
            done[node] = extract_node(backend, frame["row"][sel], frame["col"][sel], frame["x"][sel], frame["y"][sel], years, tile_px, master_res_m)
        except retriable as e:
            failed.append((int(node), str(e)))
            continue
        if on_node is not None:
            on_node(done, node)
    return done, failed


def assemble(frame, done, fields=("kndvi", "valid_count", "no_disturbance", "elevation", "height_mean", "height_n_subpoints")):
    """Frame pixels of the finished nodes with their extracted values, in frame order. Pixels of nodes not in ``done`` are dropped."""
    nodes = frame["node"]
    if nodes.size and np.any(np.diff(nodes) < 0):
        raise ValueError("frame pixels must be grouped by ascending node id")
    keep = np.flatnonzero(np.isin(nodes, list(done)))
    out = {k: np.asarray(frame[k])[keep] for k in ("node", "row", "col", "x", "y", "lat", "lon", "forest_class")}
    for field in fields:
        out[field] = np.concatenate([done[n][field] for n in np.unique(out["node"])]) if keep.size else np.zeros(0)
    return out


def disturbance_layer_problem(no_disturbance, years, first_year: int = 2001, last_year: int = 2023, min_share: float = 0.5):
    """None when a harvest/fire layer looks usable, else a sentence saying what is wrong.

    ``no_disturbance`` is (pixels, years) with True = no harvest or fire flagged. Harvest and fire touch a few percent of forest at most in a year, so in every year with loss data
    (Hansen covers 2001-2023) most pixels must be undisturbed. A layer that reads "disturbed" almost everywhere cannot be a record of events: it is a masking error (the first export
    of this layer was exactly that), and an event rule that needs "no disturbance" can then never fire.
    """
    nd = np.asarray(no_disturbance, dtype=bool)
    if nd.ndim != 2 or nd.shape[0] == 0:
        return "no pixels to check the disturbance layer on"
    years = list(years)
    share = nd.mean(axis=0)
    bad = [(y, float(share[years.index(y)])) for y in range(first_year, last_year + 1) if y in years and share[years.index(y)] < min_share]
    if bad:
        shown = ", ".join(f"{y}: {s:.0%}" for y, s in bad[:5])
        return (f"the disturbance layer flags more than {1 - min_share:.0%} of forest pixels as disturbed in {len(bad)} of the years with loss data ({shown}); "
                f"it is not a plausible record of harvest and fire")
    return None


# ---- panel assembly ---------------------------------------------------------------------------------------------------------------------------------------------

def subset_static(static: dict, idx):
    """The ``static`` grid inputs of the nodes ``idx`` only (arrays are indexed, the nested soil dict too; the token is dropped)."""
    idx = np.asarray(idx, dtype=int)
    out = {}
    for k, v in static.items():
        if k == "token":
            out[k] = None
        elif isinstance(v, dict):
            out[k] = {kk: np.asarray(vv)[idx] for kk, vv in v.items()}
        else:
            out[k] = np.asarray(v)[idx]
    return out


def remove_year_effect(kndvi):
    """Subtract the year effect common to the forest pixels from every series.

    The vitality composites merge sensors whose number and observations per composite change over the record (a median of 1-8 valid observations before 2014, 13-31 after), so every
    pixel steps together in some years. A pixel's deviation from its own mean is computed, the median of those deviations over all pixels gives the year effect, and it is removed. The
    median ignores a local collapse; a country-wide drought is removed with the artefacts, which is why this variant is reported next to the unadjusted rule and never in place of it.
    """
    k = np.asarray(kndvi, dtype=float)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            dev = k - np.nanmean(k, axis=1, keepdims=True)
            year_effect = np.nanmedian(dev, axis=0)
    return k - np.where(np.isfinite(year_effect), year_effect, 0.0)[None, :]


def _onsets(kndvi, no_disturbance, years, panel_years):
    """z (n, years) and, per pixel, the index of the first strict onset and of the first candidate onset (flag ignored) inside the panel years (len(years) = none)."""
    from antar.hazard.observation import dieback_event, trailing_z
    n, ny = kndvi.shape
    in_panel = np.isin(years, list(panel_years))
    all_nd = np.ones(ny, dtype=bool)
    ev = np.zeros((n, ny), dtype=bool)
    cand = np.zeros((n, ny), dtype=bool)
    z = np.full((n, ny), np.nan)
    for i in range(n):
        ev[i] = dieback_event(kndvi[i], no_disturbance[i]) & in_panel
        cand[i] = dieback_event(kndvi[i], all_nd) & in_panel
        z[i] = trailing_z(kndvi[i])
    first_ev = np.where(ev.any(axis=1), ev.argmax(axis=1), ny)
    first_cand = np.where(cand.any(axis=1), cand.argmax(axis=1), ny)
    return z, first_ev, first_cand


def build_person_years(pix: dict, years, panel_years, climate, climate_cols):
    """One row per forest pixel and panel year.

    ``pix`` holds arrays over pixels: node, lat, lon, forest_class, kndvi (n, years), valid_count, no_disturbance, height_mean, elevation. ``climate`` is a DataFrame with one row per
    node and year (columns node, year and ``climate_cols``). Events and the continuous anomaly are computed on each pixel's full 2000-2024 series: the strict dieback rule needs no
    harvest or fire flag, ``candidate`` is the same rule ignoring the flag, ``z`` is the standardised anomaly. A collapse that persists satisfies the rule in every following year, so
    ``event`` and ``candidate`` mark only the first onset of a pixel within the panel years, and ``at_risk`` is False for the years after it (first-event analysis: a pixel that has
    died back is no longer at risk of dying back). The columns ``z_adj``, ``event_adj`` and ``candidate_adj`` are the same on the series with the common year effect removed
    (:func:`remove_year_effect`). Pixel-years without a finite kNDVI or without climate for the node are dropped.
    """
    import pandas as pd
    years = list(years)
    nd_all = np.asarray(pix["no_disturbance"], dtype=bool)
    kn = np.asarray(pix["kndvi"], dtype=float)
    z, first_ev, first_cand = _onsets(kn, nd_all, years, panel_years)
    kn_adj = remove_year_effect(kn)
    z_adj, first_ev_adj, first_cand_adj = _onsets(kn_adj, nd_all, years, panel_years)
    clim = climate.set_index(["node", "year"])
    rows = []
    for yr in panel_years:
        yi = years.index(yr)
        # a harvest or fire flagged in the year or in the 10 years before it sits in the baseline window of z: the anomaly is then not a weather response
        disturbed_recent = ~nd_all[:, max(0, yi - 10): yi + 1].all(axis=1)
        sel = np.flatnonzero(np.isfinite(kn[:, yi]))
        keys = list(zip(np.asarray(pix["node"])[sel].tolist(), [yr] * len(sel)))
        have = np.array([k in clim.index for k in keys], dtype=bool)
        sel = sel[have]
        if not len(sel):
            continue
        c = clim.loc[[(int(pix["node"][i]), yr) for i in sel], list(climate_cols)].to_numpy()
        df = pd.DataFrame({
            "pixel": sel, "node": np.asarray(pix["node"])[sel], "year": yr, "lat": np.asarray(pix["lat"])[sel], "lon": np.asarray(pix["lon"])[sel],
            "forest_class": np.asarray(pix["forest_class"])[sel], "event": first_ev[sel] == yi, "candidate": first_cand[sel] == yi,
            "at_risk": yi <= first_ev[sel], "candidate_at_risk": yi <= first_cand[sel], "z": z[sel, yi],
            "event_adj": first_ev_adj[sel] == yi, "candidate_adj": first_cand_adj[sel] == yi, "at_risk_adj": yi <= first_ev_adj[sel], "z_adj": z_adj[sel, yi],
            "kndvi": kn[sel, yi], "valid_count": np.asarray(pix["valid_count"])[sel, yi],
            "no_disturbance": nd_all[sel, yi], "disturbed_recent": disturbed_recent[sel], "height_m": np.asarray(pix["height_mean"])[sel], "elevation_m": np.asarray(pix["elevation"])[sel],
        })
        for j, col in enumerate(climate_cols):
            df[col] = c[:, j]
        rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def count_onset_pixels(kndvi, years, panel_years):
    """Number of pixels with at least one onset of the dieback rule inside the panel years, the harvest/fire flag ignored (the rule on the vitality series alone)."""
    from antar.hazard.observation import dieback_event
    years = np.asarray(list(years))
    in_panel = np.isin(years, list(panel_years))
    all_nd = np.ones(len(years), dtype=bool)
    return int(sum(bool((dieback_event(series, all_nd) & in_panel).any()) for series in np.asarray(kndvi, dtype=float)))


def mirror_onset_count(kndvi, years, panel_years):
    """Pixels with a persistent *rise* by the same rule: the series is mirrored (x -> -x, so z -> -z) and the dieback rule is applied, the flag ignored.

    A decline that is real shows up as more pixels with an onset of decline than with an onset of rise; as many of one as of the other means the rule is reading variation, not
    dieback. (A shuffled record is not a fair null: kNDVI series are autocorrelated and carry sensor-driven level shifts, so shuffling manufactures far more "collapses" than exist.)
    """
    return count_onset_pixels(-np.asarray(kndvi, dtype=float), years, panel_years)
