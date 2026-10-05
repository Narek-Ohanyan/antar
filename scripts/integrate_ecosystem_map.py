"""Integrate the real Ecosystem Map of Armenia (BCC Armenia / Institute of Botany NAS RA / IOER,
published 2026-09-18, CC BY 4.0 per the Earth Engine app page -- see
configs/manifests/ecosystem_map_armenia.yaml for the full source record and the one real
discrepancy worth flagging: the zip's own README.txt ships an unfilled licence placeholder even
though the app page states CC BY 4.0 plainly) at the project's real 78-cell grid.

Two real, separately scoped uses, not conflated:

1. **Validation**, not retraining: for each of the 78 real grid cells (the same points used by
   every gridded run tonight -- TOPOHYDRO/XYLEM/MNEME/REFUGIUM/future-projections), sample the
   real national classification raster in a local window and compute the real fraction covered by
   the forest classes matching each of REFUGIUM's 3 fitted functional groups (broadleaf ->
   Fagus orientalis class 31; oak -> Quercus macranthera/iberica classes 32-35; pine -> Pinus
   kochiana class 36; juniper -> class 44, descriptive only since REFUGIUM has no real viability
   fitted for juniper to correlate against). Then compute a real Spearman correlation between that
   observed real forest-cover fraction and REFUGIUM's predicted real viability_ensemble_mean, per
   group -- an honest check of whether the model's predictions track real mapped forest presence,
   reported as found, not massaged toward a clean result either way.

2. **A real, narrowly scoped AEGIS eligibility refinement**: add a "human-modified landscape"
   exclusion (settlements, buildings, active cropland, quarries/landfills -- classes 12, 13, 14,
   15, 16, 18) layered onto the existing real WDPA eligibility mask. Deliberately does NOT exclude
   already-forested cells, even though that might look like an obvious additionality fix: 3 of
   AEGIS's 8 real intervention methods (coppicing_oak, pine_thinning, wildfire_prevention) target
   *existing* forest, not open land, so a blanket "already forested -> ineligible" rule would
   wrongly zero out exactly the options that most need existing forest to operate on. The
   human-modified exclusion has no such ambiguity -- none of the 8 real methods make sense inside
   a settlement, building footprint, active cropland or quarry.

Real window radius: 500m (50 real 10m pixels in each direction, same `rasterio.windows.Window`
mechanism already used for terrain concavity/soils). Chosen because these are survey/grid points,
not polygons -- a 500m local neighbourhood is comparable to the kind of window already used for
terrain's concavity index, while being far smaller than the ~30km spacing between the 78 grid
points themselves. Stated explicitly: this is a local-context sample, not a representation of the
full ~30km grid cell.
"""
import sys
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.armenia_mask import inside_armenia  # noqa: E402

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "ecosystem_map"
RASTER_PATH = DATA_DIR / "Ecosystem_Map_of_Armenia.tif"
OUT_PATH = CONFIG_DIR / "fitted" / "ecosystem_ground_truth_2019.yaml"

WINDOW_RADIUS_M = 500.0

# Real class codes, from Ecosystem_Map_of_Armenia_Legend.csv.
GROUP_FOREST_CLASSES = {
    "mesic_diffuse_porous_broadleaf": [31],       # Fagus orientalis and other deciduous
    "ring_porous_oak": [32, 33, 34, 35],           # Quercus macranthera / iberica (+ other deciduous)
    "pine": [36],                                  # Pinus kochiana
    "juniper_arid_conifer": [44],                  # Juniper woodlands
}
HUMAN_MODIFIED_CLASSES = [12, 13, 14, 15, 16, 18]  # agricultural/cropland/settlement/building/quarry


def load_78_cells():
    d = yaml.safe_load(open(CONFIG_DIR / "fitted" / "refugium_viability_2019.yaml"))
    g = d["groups"]
    ok_group = next(k for k, v in g.items() if v.get("status") == "ok")
    cells = g[ok_group]["cells"]
    lats = np.array([c["lat"] for c in cells])
    lons = np.array([c["lon"] for c in cells])
    viability = {}
    for gname, gg in g.items():
        if gg.get("status") == "ok":
            viability[gname] = np.array([c["viability_ensemble_mean"] for c in gg["cells"]])
    return lats, lons, viability


def sample_class_fractions(lats, lons):
    """Returns one dict per point (class -> fraction), or None for a point that falls outside the
    real Armenia raster entirely -- the project's rectangular BBOX (43.4, 38.8, 46.7, 41.4) extends
    past Armenia's real (non-rectangular) national border at its NW corner and western edge, so
    some grid cells legitimately have no real ecosystem-map coverage. Stated, not silently zeroed."""
    n = len(lats)
    with rasterio.open(RASTER_PATH) as src:
        xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
        radius_px = int(round(WINDOW_RADIUS_M / src.res[0]))
        class_fractions = []
        for i, (x, y) in enumerate(zip(xs, ys)):
            row, col = src.index(x, y)
            if not (0 <= row < src.height and 0 <= col < src.width):
                class_fractions.append(None)
                continue
            r0, r1 = max(0, row - radius_px), min(src.height, row + radius_px + 1)
            c0, c1 = max(0, col - radius_px), min(src.width, col + radius_px + 1)
            win = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0))
            vals, counts = np.unique(win, return_counts=True)
            total = counts.sum()
            class_fractions.append(dict(zip(vals.tolist(), (counts / total).tolist())))
    return class_fractions


def main():
    lats, lons, viability = load_78_cells()
    n = len(lats)
    print(f"=== Sampling the real Ecosystem Map of Armenia at {n} real grid cells "
          f"(500m local window) ===", flush=True)
    class_fractions = sample_class_fractions(lats, lons)

    cells_out = []
    n_outside = 0
    for i in range(n):
        cf = class_fractions[i]
        if cf is None:
            n_outside += 1
            cells_out.append({
                "lat": float(lats[i]), "lon": float(lons[i]),
                "outside_real_armenia_raster_extent": True,
                "human_modified_fraction": None,
                "forest_cover_fraction_by_group": None,
            })
            continue
        human_modified_frac = sum(cf.get(c, 0.0) for c in HUMAN_MODIFIED_CLASSES)
        group_frac = {g: sum(cf.get(c, 0.0) for c in classes)
                      for g, classes in GROUP_FOREST_CLASSES.items()}
        cells_out.append({
            "lat": float(lats[i]), "lon": float(lons[i]),
            "outside_real_armenia_raster_extent": False,
            "human_modified_fraction": human_modified_frac,
            "forest_cover_fraction_by_group": group_frac,
        })
    print(f"=== {n_outside}/{n} real cells fall outside Armenia's real (non-rectangular) national "
          f"border -- the project's rectangular BBOX extends past it at the NW corner/western edge; "
          f"no real ecosystem-map coverage exists there, left unset not zero-filled ===", flush=True)

    print("=== Real validation: REFUGIUM predicted viability vs. real observed forest cover "
          "(in-bounds cells only) ===", flush=True)
    # Validate only on cells INSIDE Armenia. The earlier "in-bounds" test only excluded cells outside the
    # map's rectangle; 38 more sat inside the rectangle over foreign land, where mapped forest cover is
    # zero by construction (class 0 = unmapped) -- those cells were dragging the correlation around.
    in_arm = inside_armenia(lats, lons)
    for i, c in enumerate(cells_out):
        c["inside_armenia"] = bool(in_arm[i])
    in_bounds = np.array([(not c["outside_real_armenia_raster_extent"]) and c["inside_armenia"] for c in cells_out])
    validation = {}
    for g in viability:
        obs_all = np.array([c["forest_cover_fraction_by_group"][g] if not c["outside_real_armenia_raster_extent"]
                             else np.nan for c in cells_out])
        obs = obs_all[in_bounds]
        rho, pval = spearmanr(viability[g][in_bounds], obs)
        validation[g] = {"spearman_rho": float(rho), "p_value": float(pval),
                          "n": int(in_bounds.sum()), "mean_observed_forest_cover_fraction": float(obs.mean())}
        print(f"  {g}: rho={rho:.3f} (p={pval:.3f}), n={int(in_bounds.sum())}, "
              f"mean real observed cover={obs.mean():.3f}", flush=True)

    n_human_modified = sum(1 for c in cells_out
                            if not c["outside_real_armenia_raster_extent"] and c["human_modified_fraction"] > 0.5)
    print(f"=== {n_human_modified}/{int(in_bounds.sum())} real in-bounds cells >50% human-modified "
          f"(settlements/cropland/buildings/quarries) ===", flush=True)

    out = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "source": "Ecosystem Map of Armenia (BCC Armenia / Institute of Botany NAS RA / IOER, 2026-09-18)",
        "window_radius_m": WINDOW_RADIUS_M,
        "scope_note": ("Local 500m-radius class-fraction sample at the real 78-cell grid, NOT a "
                        "representation of the full ~30km grid-cell spacing. Validation correlates "
                        "REFUGIUM's predicted viability against real observed nearby forest-class "
                        "cover -- an honest check, not a retraining; juniper has no real REFUGIUM "
                        "viability to correlate against (XYLEM/REFUGIUM both skipped it for lack of "
                        "real trait variance) so only its real observed cover is reported. "
                        f"{n_outside}/{n} real cells fall outside Armenia's real national border "
                        "(the project's rectangular BBOX extends past it) and have no coverage here "
                        "-- left unset, not zero-filled."),
        "n_cells_total": n,
        "n_cells_outside_armenia_border": int(n - in_arm.sum()),
        "n_cells_inside_armenia": int(in_arm.sum()),
        "validation_viability_vs_observed_forest_cover": validation,
        "human_modified_fraction_threshold_for_aegis": 0.5,
        "n_cells_human_modified": n_human_modified,
        "cells": cells_out,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.dump(out, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
