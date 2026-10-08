"""Forest-pixel frame and Drive extraction for MNEME's panel.

MNEME's first panel followed one 30 m pixel at each grid node, whatever the land cover, so a "dieback" rule was applied to grassland and cropland as well as forest. This script
builds the frame the panel should have used and reads what the panel needs from the exported rasters:

1. Frame: for every node with usable soil data, up to K random 30 m pixels inside the node's climate cell (the 30-arcsecond CHELSA pixel around it) whose nine 10 m sub-points all
   lie on a closed-forest class of the national Ecosystem Map (``antar.hazard.forest_panel``). Local, no network. Cached with its parameters.
2. Extraction, one node at a time, in the order of the export blocks so neighbouring nodes reuse the blocks already fetched: kNDVI and valid-observation counts 2000-2024 (vitality
   tiles), the harvest/fire mask (disturbance layer), elevation (terrain layer) and the 2019 canopy height from the Sentinel-2 vegetation height model (WSL, Drive copy of the file
   registered in configs/manifests/s2_vhm.yaml), averaged over the pixel's nine sub-points. Progress is saved every few nodes and a rerun resumes; a node that fails every retry is
   left out and listed, and the run fails loudly if more than 3 % of nodes are lost.

Output: data/_cache/mneme_forest_pixels_dense.npz (arrays over pixels) and mneme_forest_frame_dense.npz.
"""
import hashlib
import os
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

os.environ.setdefault("GDAL_CACHEMAX", "1024")                      # MB; the blocks hold 50 bands of 256 x 256 floats, so keep many of them
os.environ.setdefault("CPL_VSIL_CURL_CACHE_SIZE", "268435456")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    CACHE_DIR, DENSE_GRID_COLS, DENSE_GRID_ROWS, DriveCoverageError, MAX_LOST_FRACTION, RETRY_SLEEPS_S, TERRAIN_DRIVE_ID,
    drive_vsicurl_url, extract_static_grid_inputs, get_access_token,
)
from mneme_inputs import DISTURBANCE_FILE_ID, VITALITY_TILE_IDS, VITALITY_TILE_SIZE_PX, VITALITY_YEARS  # noqa: E402
from antar.hazard import forest_panel as fp  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MAP_PATH = ROOT / "data" / "ecosystem_map" / "Ecosystem_Map_of_Armenia.tif"
MANIFEST = ROOT / "configs" / "manifests" / "s2_vhm.yaml"
HEIGHT_VARIABLE = "canopy_height_max_2019"

MASTER_CRS = "EPSG:32638"
MASTER_X0, MASTER_Y0, MASTER_RES_M = 361050.0, 4584570.0, 30.0
MASTER_SHAPE = (9626, 9553)                                          # rows, cols (terrain layer)
HEIGHT_X0, HEIGHT_Y0, HEIGHT_RES_M, HEIGHT_SHAPE = 366100.0, 4585600.0, 10.0, (28890, 28260)
HALF_CELL_DEG = 15.0 / 3600.0                                        # the CHELSA cell is 30 arcseconds wide
def _int_arg(name, default):
    for a in sys.argv[1:]:
        if a.startswith(name + "="):
            return int(a.split("=", 1)[1])
    return default


# Up to this many pixels are drawn per climate cell (pass --pixels-per-node=100000 for every qualifying pixel). Another cap writes its own files, so a trial run never overwrites a result.
PIXELS_PER_NODE, MIN_CANDIDATES, SEED = _int_arg("--pixels-per-node", 25), 5, 20261007
FILE_TAG = "" if PIXELS_PER_NODE == 25 else f"_k{PIXELS_PER_NODE}"
SAVE_EVERY_NODES = 10
TOKEN_MAX_AGE_S = 2400                                               # Google access tokens last an hour

FRAME_PATH = CACHE_DIR / f"mneme_forest_frame_dense{FILE_TAG}.npz"
PIXELS_PATH = CACHE_DIR / f"mneme_forest_pixels_dense{FILE_TAG}.npz"
PARTIAL_PATH = CACHE_DIR / f"mneme_forest_pixels_dense{FILE_TAG}_partial.npz"


def to_xy(lon, lat):
    x, y = warp_transform("EPSG:4326", MASTER_CRS, np.asarray(lon, float).tolist(), np.asarray(lat, float).tolist())
    return np.array(x), np.array(y)


def to_lonlat(x, y):
    lon, lat = warp_transform(MASTER_CRS, "EPSG:4326", np.asarray(x, float).tolist(), np.asarray(y, float).tolist())
    return np.array(lon), np.array(lat)


def make_class_reader(src):
    """Ecosystem Map window over a rectangle, class 0 where the rectangle leaves the map."""
    ox, oy, res = src.transform.c, src.transform.f, src.res[0]

    def reader(xmin, ymin, xmax, ymax):
        c0, c1 = int(np.floor((xmin - ox) / res)), int(np.floor((xmax - ox) / res))
        r0, r1 = int(np.floor((oy - ymax) / res)), int(np.floor((oy - ymin) / res))
        c0c, c1c, r0c, r1c = max(c0, 0), min(c1, src.width - 1), max(r0, 0), min(r1, src.height - 1)
        if c1c < c0c or r1c < r0c:
            return np.zeros((1, 1), dtype=np.int16), xmin, ymax, res
        arr = src.read(1, window=Window(c0c, r0c, c1c - c0c + 1, r1c - r0c + 1))
        return arr, ox + c0c * res, oy - r0c * res, res
    return reader


def frame_key(static):
    h = hashlib.sha1(np.round(np.concatenate([static["lats"], static["lons"]]), 6).tobytes()).hexdigest()[:10]
    return f"{len(static['lats'])}_{h}_k{PIXELS_PER_NODE}_m{MIN_CANDIDATES}_s{SEED}_c{''.join(map(str, fp.CLOSED_FOREST_CLASSES))}"


def build_or_load_frame(static):
    key = frame_key(static)
    if FRAME_PATH.exists():
        z = np.load(FRAME_PATH, allow_pickle=False)
        if str(z["key"]) == key:
            print(f"=== forest frame loaded from {FRAME_PATH.name} ===", flush=True)
            return {k: z[k] for k in z.files if k != "key"}
    usable = np.flatnonzero(static["valid_soil"])
    print(f"=== Forest frame: {len(usable)} nodes with soil data, up to {PIXELS_PER_NODE} pure closed-forest 30 m pixels per climate cell ===", flush=True)
    with rasterio.open(MAP_PATH) as src:
        assert src.crs.to_string() == MASTER_CRS, f"Ecosystem Map CRS {src.crs} is not {MASTER_CRS}"
        frame = fp.build_frame(static["lats"][usable], static["lons"][usable], make_class_reader(src), to_xy, to_lonlat, HALF_CELL_DEG,
                               MASTER_X0, MASTER_Y0, MASTER_RES_M, k=PIXELS_PER_NODE, seed=SEED, min_candidates=MIN_CANDIDATES, node_ids=usable)
    inside = (frame["row"] >= 0) & (frame["row"] < MASTER_SHAPE[0]) & (frame["col"] >= 0) & (frame["col"] < MASTER_SHAPE[1])
    if not inside.all():
        print(f"  {int((~inside).sum())} pixels fall off the master grid and are dropped", flush=True)
        for k in ("node", "row", "col", "x", "y", "lat", "lon", "forest_class"):
            frame[k] = frame[k][inside]
    frame["usable_nodes"] = usable
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FRAME_PATH, key=key, **frame)
    return frame


class DriveBackend:
    """The ``backend`` that ``forest_panel.extract_node`` reads through: windows of the Drive-hosted rasters, with retry, backoff and a token refreshed before it expires."""
    x0, y0 = MASTER_X0, MASTER_Y0
    height_origin = (HEIGHT_X0, HEIGHT_Y0, HEIGHT_RES_M)
    height_shape = HEIGHT_SHAPE

    def __init__(self, height_file_id):
        self.height_file_id = height_file_id
        self._tok, self._tok_t = None, 0.0
        self.height_nodata = None

    def _token(self, force=False):
        if force or self._tok is None or time.time() - self._tok_t > TOKEN_MAX_AGE_S:
            self._tok, self._tok_t = get_access_token(), time.time()
        return self._tok

    def _read(self, file_id, window):
        r, c, h, w = window
        last = None
        for attempt in range(len(RETRY_SLEEPS_S) + 1):
            try:
                opts = dict(GDAL_HTTP_HEADERS=f"Authorization: Bearer {self._token(force=attempt > 0)}", GDAL_DISABLE_READDIR_ON_OPEN="YES",
                            GDAL_HTTP_TIMEOUT=60, GDAL_HTTP_CONNECTTIMEOUT=15)
                with rasterio.Env(**opts):
                    with rasterio.open(drive_vsicurl_url(file_id)) as src:
                        arr = src.read(window=Window(c, r, w, h))
                        return arr, list(src.descriptions), src.nodata
            except rasterio.errors.RasterioIOError as e:
                last = e
                if attempt < len(RETRY_SLEEPS_S):
                    print(f"    read failed ({e}); retry {attempt + 1}/{len(RETRY_SLEEPS_S)} in {RETRY_SLEEPS_S[attempt]} s", flush=True)
                    time.sleep(RETRY_SLEEPS_S[attempt])
        raise last

    def vitality(self, tile_row0, tile_col0, window):
        arr, names, _ = self._read(VITALITY_TILE_IDS[(tile_row0, tile_col0)], window)
        return arr, names

    def disturbance(self, window):
        arr, names, _ = self._read(DISTURBANCE_FILE_ID, window)
        return arr, names

    def terrain(self, window):
        arr, names, _ = self._read(TERRAIN_DRIVE_ID, window)
        return arr, names

    def height(self, window):
        arr, _, nodata = self._read(self.height_file_id, window)
        arr = arr[0].astype(float)
        if nodata is not None and np.isfinite(nodata):
            arr[arr == nodata] = np.nan
        return arr

    def check_layout(self):
        """Fail before any extraction if a layer is not on the grid the frame assumes."""
        for label, fid, shape, origin, res in (("terrain", TERRAIN_DRIVE_ID, MASTER_SHAPE, (MASTER_X0, MASTER_Y0), MASTER_RES_M),
                                              ("height model", self.height_file_id, HEIGHT_SHAPE, (HEIGHT_X0, HEIGHT_Y0), HEIGHT_RES_M)):
            with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {self._token()}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=60):
                with rasterio.open(drive_vsicurl_url(fid)) as src:
                    ok = ((src.height, src.width) == tuple(shape) and src.crs.to_string() == MASTER_CRS and
                          (src.transform.c, src.transform.f) == origin and src.res[0] == res)
                    assert ok, f"{label} layer is not on the expected grid: {src.height}x{src.width}, {src.crs}, origin {(src.transform.c, src.transform.f)}, res {src.res}"
        print("=== layer grids checked: terrain (30 m) and height model (10 m) are where the frame expects them ===", flush=True)


def height_file_id():
    entries = yaml.safe_load(MANIFEST.read_text())
    e = next(e for e in entries if e["variable"] == HEIGHT_VARIABLE)
    import re
    return re.search(r"file id=([\w-]+)", e["notes"]).group(1)


def save_partial(done):
    arrays = {}
    for n, res in done.items():
        for k, v in res.items():
            arrays[f"n{n}__{k}"] = v
    arrays["done_nodes"] = np.array(sorted(done), dtype=int)
    tmp = PARTIAL_PATH.with_name(PARTIAL_PATH.stem + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.replace(PARTIAL_PATH)


def load_partial():
    if not PARTIAL_PATH.exists():
        return {}
    z = np.load(PARTIAL_PATH, allow_pickle=False)
    done = {}
    for n in z["done_nodes"].tolist():
        pre = f"n{n}__"
        done[n] = {k[len(pre):]: z[k] for k in z.files if k.startswith(pre)}
    return done


def main():
    if "--dense" not in sys.argv:
        sys.exit("MNEME's forest panel is defined on the dense grid only: pass --dense")
    static = extract_static_grid_inputs(grid_rows=DENSE_GRID_ROWS, grid_cols=DENSE_GRID_COLS)
    frame = build_or_load_frame(static)
    n_nodes = len(np.unique(frame["node"]))
    print(f"=== {len(frame['node'])} forest pixels in {n_nodes} nodes "
          f"({int((frame['n_candidates'] > 0).sum())} nodes have some pure forest pixel; {int((frame['n_candidates'] < MIN_CANDIDATES).sum())} have fewer than "
          f"{MIN_CANDIDATES} and are left out) ===", flush=True)
    if n_nodes == 0:
        sys.exit("no node has forest pixels: check the Ecosystem Map and the grid")

    backend = DriveBackend(height_file_id())
    backend.check_layout()
    done = load_partial()
    if done:
        print(f"=== resuming: {len(done)}/{n_nodes} nodes already extracted ===", flush=True)
    t0 = time.time()

    def on_node(finished, node):
        if len(finished) == 5:                                            # fail on the first few nodes, not after hours, if the harvest/fire layer is not usable
            problem = fp.disturbance_layer_problem(np.concatenate([r["no_disturbance"] for r in finished.values()]), VITALITY_YEARS)
            if problem:
                raise SystemExit(f"STOP: {problem}. Check DISTURBANCE_FILE_ID in mneme_inputs.py.")
        if len(finished) % SAVE_EVERY_NODES == 0:
            save_partial(finished)
            print(f"  {len(finished)}/{n_nodes} nodes ({(time.time() - t0) / 60:.1f} min elapsed)", flush=True)

    done, failed_info = fp.extract_all(frame, backend, VITALITY_YEARS, VITALITY_TILE_SIZE_PX, MASTER_RES_M, done=done, retriable=(OSError,), on_node=on_node)
    save_partial(done)
    failed = [n for n, _ in failed_info]
    for n, msg in failed_info:
        print(f"  node {n}: FAILED after all retries ({msg})", flush=True)
    if len(failed) / n_nodes > MAX_LOST_FRACTION:
        raise DriveCoverageError(f"{len(failed)}/{n_nodes} nodes could not be read ({failed[:10]}...): refusing to build a panel from a degraded sample; rerun to resume")

    out = fp.assemble(frame, done)
    out["pixels_per_node"] = np.array(PIXELS_PER_NODE)
    out["n_qualifying_pixels"] = np.array(int(frame["n_candidates"][frame["kept"]].sum()))
    out["years"] = np.array(VITALITY_YEARS)
    out["failed_nodes"] = np.array(failed, dtype=int)
    np.savez_compressed(PIXELS_PATH, **out)
    k19 = VITALITY_YEARS.index(2019)
    print(f"=== wrote {PIXELS_PATH.name}: {len(out['node'])} pixels, {len(np.unique(out['node']))} nodes, {len(failed)} nodes failed ===", flush=True)
    print(f"  kNDVI 2019: mean {np.nanmean(out['kndvi'][:, k19]):.3f}, missing {np.isnan(out['kndvi'][:, k19]).mean():.1%}; "
          f"canopy height: median {np.nanmedian(out['height_mean']):.1f} m, missing {np.isnan(out['height_mean']).mean():.1%}; "
          f"elevation {np.nanmin(out['elevation']):.0f}-{np.nanmax(out['elevation']):.0f} m", flush=True)
    if PARTIAL_PATH.exists() and not failed:
        PARTIAL_PATH.unlink()


if __name__ == "__main__":
    sys.exit(main())
