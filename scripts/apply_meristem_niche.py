"""Apply the fitted MERISTEM niche models at the dense grid nodes, for the present climate and for each of the 45 scenario members (antar.niche.application).

For each species whose model passes the quality gate (mean spatial-block Boyce index >= MIN_BOYCE), a node is within the niche when its fitted score reaches the species' presence threshold
(the score below which a tenth of its own records fall). A species group is niche-supported at a node when any passing species of the group is. A group with no passing species gets no
support anywhere and is reported as not assessed, which is not the same as suitable everywhere.

Predictors at a node, in the order the models were fitted: winter minimum temperature, growing degree days (5 C) and mean vapour-pressure deficit from CHELSA-BIOCLIM+ (the member's own
period, model and path for the future; the vapour-pressure deficit has no future layer and is scaled with the member's warming at constant relative humidity), the water deficit of the
group (TOPOHYDRO, 2019 or the member's), and the SoilGrids clay, sand, silt and organic carbon of 0-30 cm.

Reads Drive once (SoilGrids at the nodes, cached); everything else is local.

    python3 scripts/apply_meristem_niche.py
"""
import datetime
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from antar.niche import application as A  # noqa: E402
from run_topohydro_grid import CACHE_DIR, DENSE_GRID_COLS, DENSE_GRID_ROWS, SOILS_DRIVE_ID, drive_vsicurl_url, extract_static_grid_inputs, get_access_token  # noqa: E402

FITTED = ROOT / "configs" / "fitted"
BIOCLIM = ROOT / "data" / "chelsa_bioclim"
NICHE_PATH = FITTED / "meristem_adult_niche.yaml"
OUT_PATH = FITTED / "meristem_niche_nodes_dense.yaml"
MIN_BOYCE = 0.30
OMISSION = 0.10
GROUP_SPECIES = {"mesic_diffuse_porous_broadleaf": ["Fagus orientalis", "Carpinus betulus"], "ring_porous_oak": ["Quercus macranthera", "Quercus iberica"], "pine": ["Pinus kochiana"]}
SOIL_BANDS = ["clay_0_30cm_mean", "sand_0_30cm_mean", "silt_0_30cm_mean", "soc_0_30cm_mean"]


def key(lat, lon):
    return (round(float(lat), 6), round(float(lon), 6))


def load(name):
    return yaml.load(open(FITTED / name), Loader=yaml.CSafeLoader)


def soils_at(lats, lons):
    """Raw SoilGrids values (clay, sand, silt, organic carbon, 0-30 cm, as the models were fitted on them) at the nodes; cached."""
    cache = CACHE_DIR / f"niche_soils_{len(lats)}_{abs(hash(tuple(np.round(lats, 5)))) % 10**8}.npz"
    if cache.exists():
        return np.load(cache)["soil"]
    last = None
    for attempt in range(4):
        try:
            with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {get_access_token()}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=60):
                with rasterio.open(drive_vsicurl_url(SOILS_DRIVE_ID)) as src:
                    assert list(src.descriptions) == SOIL_BANDS, f"soil bands are {src.descriptions}"
                    xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
                    soil = np.array(list(src.sample(zip(xs, ys))), dtype=float)
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cache, soil=soil)
            return soil
        except rasterio.errors.RasterioIOError as e:
            last = e
            time.sleep(10 * (attempt + 1))
    raise last


def bioclim_layers(variables):
    out = {}
    for v in variables:
        z = np.load(BIOCLIM / f"CHELSA_BIOCLIM_{v}_armenia.npz")["data"]
        keys = json.load(open(BIOCLIM / f".{v}_progress.json"))["keys"]
        out[v] = (z, {k: i for i, k in enumerate(keys)})
    return out


def main():
    niche = yaml.safe_load(open(NICHE_PATH))
    names = niche["features"]
    expected = ["bio06_winter_min_C", "gdd5", "vpdmean_Pa", "cwd_mm", "clay_0_30cm_mean", "sand_0_30cm_mean", "silt_0_30cm_mean", "soc_0_30cm_mean"]
    if names != expected:
        sys.exit(f"the niche file's features are {names}; refit MERISTEM with the real water deficit first (scripts/fit_meristem_adult_niche.py)")
    ref = load("refugium_viability_2019_dense.yaml")
    proj = load("future_projections_dense.yaml")
    topo = load("topohydro_grid_run_2019_dense.yaml")
    static = extract_static_grid_inputs(DENSE_GRID_ROWS, DENSE_GRID_COLS)
    srow = {key(a, b): (int(r), int(c)) for a, b, r, c in zip(static["lats"], static["lons"], static["chelsa_row"], static["chelsa_col"])}
    base_pt = {key(p["lat"], p["lon"]): p for p in topo["points"] if p.get("status") == "ok"}

    ref_cells = {g: ref["groups"][g]["cells"] for g in GROUP_SPECIES if ref["groups"][g].get("status") == "ok"}
    nodes = [(c["lat"], c["lon"]) for c in next(iter(ref_cells.values()))]
    keys = [key(a, b) for a, b in nodes]
    rows = np.array([srow[k][0] for k in keys])
    cols = np.array([srow[k][1] for k in keys])
    lats, lons = np.array([a for a, _ in nodes]), np.array([b for _, b in nodes])
    soil = soils_at(lats, lons)
    layers = bioclim_layers(["bio06", "gdd5", "vpdmean"])
    hist = {v: layers[v][0][layers[v][1]["historical"], rows, cols].astype(float) for v in layers}
    t_base = np.array([base_pt[k]["t_mean_c_annual_mean"] for k in keys])

    species_info, passing = {}, {}
    for sp, v in niche["species"].items():
        ok = A.passes_gate(v["status"], (v.get("spatial_block_cv") or {}).get("boyce_index_mean"), MIN_BOYCE)
        species_info[sp] = {"status": v["status"], "boyce_mean": (v.get("spatial_block_cv") or {}).get("boyce_index_mean"), "passes_gate": bool(ok),
                            "presence_threshold": (v.get("presence_threshold") or {}).get("score")}
    for g, spp in GROUP_SPECIES.items():
        passing[g] = [sp for sp in spp if species_info.get(sp, {}).get("passes_gate")]

    def suitability(group, cwd, bio06, gdd5, vpd):
        by_species = {}
        X = np.column_stack([bio06, gdd5, vpd, cwd, soil])
        for sp in passing[group]:
            c = niche["species"][sp]["coefficients"]
            score = A.logistic_score([c[n] for n in names], niche["species"][sp]["intercept"], X)
            by_species[sp] = np.isfinite(score) & (score >= niche["species"][sp]["presence_threshold"]["score"])
        return A.group_support(by_species, passing[group])

    out = {"run_date": datetime.date.today().isoformat(), "grid": "dense_armenia_stride7", "min_boyce_gate": MIN_BOYCE, "presence_omission": OMISSION,
           "cells": {"lat": [float(v) for v in lats], "lon": [float(v) for v in lons]}, "species": species_info,
           "method": ("A node is within a species' niche when the fitted score reaches the score below which a tenth of the species' own records fall. Only species whose mean spatial-block "
                      "Boyce index is at least the gate are used; a group is supported where any of its passing species is. Future predictors: CHELSA-BIOCLIM+ winter minimum and growing "
                      "degree days of the member's period (2041-2070 for 2050, 2071-2100 for 2080 and 2100), the vapour-pressure deficit scaled by the saturation vapour pressure at the member's "
                      "warming (constant relative humidity), the member's own water deficit, unchanged soils."), "groups": {}}
    for g, cells in ref_cells.items():
        entry = {"species_used": passing[g], "n_nodes": len(nodes)}
        if not passing[g]:
            entry.update(status="not_assessed", reason="no species model of this group passes the quality gate (mean Boyce index >= %.2f)" % MIN_BOYCE)
            out["groups"][g] = entry
            continue
        cwd0 = np.array([c["cwd_mm"] for c in cells])
        base = suitability(g, cwd0, hist["bio06"], hist["gdd5"], hist["vpdmean"])
        entry.update(status="applied", baseline=[int(v) for v in base], n_supported_baseline=int(base.sum()), members={})
        for mid, m in proj["members"].items():
            ck = A.chelsa_key(m["horizon"], m["gcm"], m["scenario"])
            fut = {v: layers[v][0][layers[v][1][ck], rows, cols].astype(float) for v in ("bio06", "gdd5")}
            clim = {key(c["lat"], c["lon"]): c for c in m["climate"]}
            gcells = {key(c["lat"], c["lon"]): c for c in m["groups"][g]["cells"]}
            cwd = np.array([gcells[k]["cwd_mm"] if k in gcells else np.nan for k in keys])
            t_mem = np.array([clim[k]["t_mean_c"] if k in clim else np.nan for k in keys])
            vpd = A.vpd_after_warming(hist["vpdmean"], t_base, t_mem - t_base)
            s = suitability(g, cwd, fut["bio06"], fut["gdd5"], vpd)
            entry["members"][mid] = [int(v) for v in s]
        counts = {}
        for mid, v in entry["members"].items():
            counts.setdefault(mid.split("__")[1] + "__" + mid.split("__")[2], []).append(sum(v))
        entry["n_supported_by_path_and_horizon"] = {k: [int(np.min(v)), float(np.mean(v)), int(np.max(v))] for k, v in sorted(counts.items())}
        out["groups"][g] = entry
        print(f"=== {g}: species used {passing[g]}; supported at {entry['n_supported_baseline']}/{len(nodes)} nodes now; by path and horizon (min, mean, max over models): "
              f"{ {k: (v[0], round(v[1]), v[2]) for k, v in entry['n_supported_by_path_and_horizon'].items()} } ===", flush=True)
    for g, e in out["groups"].items():
        if e["status"] != "applied":
            print(f"=== {g}: not assessed ({e['reason']}) ===", flush=True)
    OUT_PATH.write_text(yaml.dump(out, sort_keys=False, default_flow_style=False))
    print(f"=== wrote {OUT_PATH.name} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
