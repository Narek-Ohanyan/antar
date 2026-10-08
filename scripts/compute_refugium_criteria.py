"""Criteria (b) and (c) of the robust refugium for every dense node, group and scenario member (antar.viability.criteria).

Post-processes results that already exist (no new simulation): the 2019 REFUGIUM run, the 45 scenario members and the Ecosystem Map.

(a) is the stored criterion: viable in at least rho of the trait-knowledge ensemble.
(b) buffer index B > 0 on the node lattice: the node's composite exposure (standardised climatic water deficit, water-stress integral and growing degree days, each against its 2019
    national mean and sd) is below the mean of the nodes within RADIUS lattice cells. The water-stress integral stands in for growing-season vapour-pressure deficit, which is not stored per node.
(c) inside the group's area of applicability: built on the nodes where the group's forest covers at least MIN_FOREST_SHARE of the node's climate cell on the Ecosystem Map, in the space
    of annual mean temperature, annual precipitation, growing degree days, and the group's climatic water deficit and water-stress integral. Members whose climate lies outside it are
    flagged, not asserted. Pine has too little mapped forest to build one, so no pine node passes (c).
"""
import datetime
import sys
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from antar.viability import criteria as C  # noqa: E402

FITTED = ROOT / "configs" / "fitted"
OUT = FITTED / "refugium_criteria_dense.yaml"
MAP_PATH = ROOT / "data" / "ecosystem_map" / "Ecosystem_Map_of_Armenia.tif"
BBOX, GRID_SHAPE, STRIDE = (43.4, 38.8, 46.7, 41.4), (312, 396), 7
LATTICE_SHAPE = (45, 57)
RADIUS = 3                                   # lattice cells, about 20 km
RADIUS_CHECK = (2, 5)                        # sensitivity of the buffer criterion to the radius
MIN_FOREST_SHARE = 0.02
HALF_CELL_DEG = 15.0 / 3600.0
GROUP_CLASSES = {"mesic_diffuse_porous_broadleaf": [31], "ring_porous_oak": [32, 33, 34, 35], "pine": [36]}
CRITERIA_GROUPS = list(GROUP_CLASSES)


def key(lat, lon):
    return (round(float(lat), 6), round(float(lon), 6))


def load(name):
    return yaml.load(open(FITTED / name), Loader=yaml.CSafeLoader)


def forest_shares(lat, lon):
    """Per node and group: the share of the node's climate cell (30 arcseconds square) covered by the group's classes on the Ecosystem Map."""
    out = {g: np.zeros(len(lat)) for g in GROUP_CLASSES}
    with rasterio.open(MAP_PATH) as src:
        for i in range(len(lat)):
            xs, ys = warp_transform("EPSG:4326", src.crs, [lon[i] - HALF_CELL_DEG, lon[i] + HALF_CELL_DEG], [lat[i] - HALF_CELL_DEG, lat[i] + HALF_CELL_DEG])
            r0, c0 = src.index(min(xs), max(ys))
            r1, c1 = src.index(max(xs), min(ys))
            r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, src.height - 1), min(c1, src.width - 1)
            if r1 < r0 or c1 < c0:
                continue
            f = C.class_fractions(src.read(1, window=Window(c0, r0, c1 - c0 + 1, r1 - r0 + 1)), GROUP_CLASSES)
            for g, v in f.items():
                out[g][i] = v
    return out


def main():
    ref = load("refugium_viability_2019_dense.yaml")
    topo = load("topohydro_grid_run_2019_dense.yaml")
    proj = load("future_projections_dense.yaml")
    base_pt = {key(p["lat"], p["lon"]): p for p in topo["points"] if p.get("status") == "ok"}

    out = {"run_date": datetime.date.today().isoformat(), "grid": "dense_armenia_stride7", "radius_lattice_cells": RADIUS, "min_forest_share": MIN_FOREST_SHARE,
           "exposure_variables": ["climatic water deficit (group rooting depth)", "water-stress integral (group rooting depth)", "growing degree days"],
           "note": "Water-stress integral stands in for growing-season VPD, which is not stored per node. Pine has too little mapped forest for an area of applicability, so no pine node passes (c).",
           "groups": {}}
    for g in CRITERIA_GROUPS:
        cells = ref["groups"][g]["cells"]
        keys = [key(c["lat"], c["lon"]) for c in cells]
        lat, lon = np.array([c["lat"] for c in cells]), np.array([c["lon"] for c in cells])
        ii, jj = C.lattice_index(lat, lon, BBOX, GRID_SHAPE, STRIDE)
        cwd0 = np.array([c["cwd_mm"] for c in cells])
        wsi0 = np.array([c["wsi"] for c in cells])
        gdd0 = np.array([base_pt[k]["gdd_cumulative_annual"] for k in keys])
        t0 = np.array([base_pt[k]["t_mean_c_annual_mean"] for k in keys])
        p0 = np.array([base_pt[k]["p_mm_annual_sum"] for k in keys])
        a0 = np.array([c["robust_refugium_criterion_a"] for c in cells])
        reference = [(float(np.mean(v)), float(np.std(v))) for v in (cwd0, wsi0, gdd0)]
        shares = forest_shares(lat, lon)[g]
        train = shares >= MIN_FOREST_SHARE
        X0 = np.column_stack([t0, p0, gdd0, cwd0, wsi0])

        def criteria_for(cwd, wsi, gdd, t, p, a, radius=RADIUS):
            b = C.node_buffer_index(C.exposure_composite([cwd, wsi, gdd], reference), ii, jj, LATTICE_SHAPE, radius)
            return b

        members = {mid: m for mid, m in proj["members"].items()}
        mem_rows = []
        for mid, m in members.items():
            clim = {key(c["lat"], c["lon"]): c for c in m["climate"]}
            gcells = {key(c["lat"], c["lon"]): c for c in m["groups"][g]["cells"]}
            cwd = np.array([gcells[k]["cwd_mm"] if k in gcells else np.nan for k in keys])
            wsi = np.array([gcells[k]["wsi"] if k in gcells else np.nan for k in keys])
            gdd = np.array([clim[k]["gdd"] if k in clim else np.nan for k in keys])
            t = np.array([clim[k]["t_mean_c"] if k in clim else np.nan for k in keys])
            p = np.array([clim[k]["precip_mm"] if k in clim else np.nan for k in keys])
            a = np.array([gcells[k]["robust_criterion_a"] if k in gcells else False for k in keys])
            mem_rows.append((mid, m, a, np.column_stack([t, p, gdd, cwd, wsi]), criteria_for(cwd, wsi, gdd, t, p, a)))

        aoa = C.group_aoa(X0[train], lat[train], lon[train], np.vstack([X0] + [r[3] for r in mem_rows]))
        n = len(cells)
        inside_all = aoa["inside"].reshape(1 + len(mem_rows), n)
        b0 = criteria_for(cwd0, wsi0, gdd0, t0, p0, a0)
        full0 = C.robust_conjunction(a0, b0, inside_all[0])
        sens = {}
        for rad in RADIUS_CHECK:
            sens[str(rad)] = int(C.robust_conjunction(a0, criteria_for(cwd0, wsi0, gdd0, t0, p0, a0, rad), inside_all[0]).sum())
        entry = {"status": aoa["status"], "n_aoa_training_nodes": aoa["n_train"], "aoa_threshold": aoa["threshold"], "n_nodes": n,
                 "baseline": {"n_a": int(a0.sum()), "n_b": int((b0 > 0).sum()), "n_c": int(inside_all[0].sum()), "n_full": int(full0.sum()),
                              "n_full_at_other_radii": sens},
                 "cells": [{"lat": float(lat[i]), "lon": float(lon[i]), "buffer_index": None if not np.isfinite(b0[i]) else round(float(b0[i]), 3),
                            "inside_aoa": bool(inside_all[0][i]), "robust_full": bool(full0[i])} for i in range(n)],
                 "members": {}}
        for r, (mid, m, a, X, b) in enumerate(mem_rows, start=1):
            full = C.robust_conjunction(a, b, inside_all[r])
            entry["members"][mid] = {"n_a": int(a.sum()), "n_b": int((b > 0).sum()), "n_c": int(inside_all[r].sum()), "n_full": int(full.sum()),
                                     "robust_full": [int(v) for v in full], "inside_aoa": [int(v) for v in inside_all[r]]}
        out["groups"][g] = entry
        e = entry["baseline"]
        print(f"=== {g}: AOA {aoa['status']} on {aoa['n_train']} forest nodes | 2019: (a) {e['n_a']}, (b) {e['n_b']}, (c) {e['n_c']}, all three {e['n_full']} of {n} "
              f"(radius {RADIUS_CHECK[0]}/{RADIUS_CHECK[1]}: {sens}) ===", flush=True)
        for mid in ("gfdl-esm4__ssp126__2050", "gfdl-esm4__ssp585__2100", "mri-esm2-0__ssp585__2100"):
            if mid in entry["members"]:
                em = entry["members"][mid]
                print(f"      {mid}: (a) {em['n_a']}, (b) {em['n_b']}, (c) {em['n_c']}, all three {em['n_full']}", flush=True)
    OUT.write_text(yaml.dump(out, sort_keys=False, default_flow_style=False))
    print(f"=== wrote {OUT.name} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
