#!/usr/bin/env python3
"""Bundle every real fitted output under configs/fitted/ into compact JSON for the static UI.

Nothing here computes science: it only reshapes already-fitted results. Two rules keep the UI
honest about what it shows:

1. **Dense results are used only when complete.** For each dataset the Armenia-only dense file is
   preferred, but only if it covers at least DENSE_MIN_COVERAGE of the dense grid. A dense file
   that exists yet falls short (e.g. one written by a run that lost most of its Drive tiles) is
   *rejected and the rejection is recorded* in the manifest, so a partial subsample can never be
   presented as the dense result. Datasets can therefore sit on different grids; every layer
   carries its grid id and the UI labels it.
2. **Status and limitations travel with the data.** Each engine block in the manifest states what
   is fitted, what is reduced-scope, and what is still a placeholder -- the UI never shows a number
   without the caveats that came with it.

Run:  python3 ui/build_data.py        (also the `ui_data` rule in workflow/Snakefile)
"""
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.armenia_mask import inside_armenia, load_mask  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
FITTED = ROOT / "configs" / "fitted"
MANIFESTS = ROOT / "configs" / "manifests"
OUT = Path(__file__).resolve().parent / "data"

DENSE_STRIDE = 7
DENSE_N = int(load_mask()[np.ix_(np.arange(0, 312, DENSE_STRIDE), np.arange(0, 396, DENSE_STRIDE))].sum())  # Armenian cells on the dense grid
DENSE_MIN_COVERAGE = 0.90
GRIDS = {
    "validation": "stride-40 validation grid (Armenian cells only)",
    "dense": f"Armenia-only dense grid (stride {DENSE_STRIDE})",
}
GROUPS = {
    "mesic_diffuse_porous_broadleaf": {"label": "Broadleaf (Fagus, Carpinus)", "short": "broadleaf"},
    "ring_porous_oak": {"label": "Oak (Quercus macranthera, Q. iberica)", "short": "oak"},
    "pine": {"label": "Pine (Pinus kochiana)", "short": "pine"},
}
ROOTING_DEPTH_M = {"mesic_diffuse_porous_broadleaf": 2.9, "ring_porous_oak": 2.9, "pine": 3.9,
                   "juniper_arid_conifer": 9.5}

try:
    _Loader = yaml.CSafeLoader
except AttributeError:  # pragma: no cover
    _Loader = yaml.SafeLoader


def load(name):
    path = FITTED / name
    if not path.exists():
        return None
    with open(path) as f:
        return yaml.load(f, Loader=_Loader)


def r(x, nd=4):
    if x is None:
        return None
    x = float(x)
    if x != x:  # NaN
        return None
    return round(x, nd)


def clean(o):
    """Recursively turn NaN / +-inf into None. Python's json writes NaN as a bare token, which is
    not valid JSON and breaks the browser's parser -- a real bug hit on the pine Spearman rho,
    which is genuinely undefined (zero mapped pine cover in any sampled window)."""
    if isinstance(o, float):
        return o if np.isfinite(o) else None
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


def dump(o, path):
    path.write_text(json.dumps(clean(o), separators=(",", ":"), allow_nan=False))


def in_arm(cells):
    """Boolean mask over a list of {lat, lon} dicts: True where the cell is inside Armenia."""
    return np.array([bool(inside_armenia(c["lat"], c["lon"])) for c in cells], dtype=bool)


def key(lat, lon):
    return (round(float(lat), 5), round(float(lon), 5))


# --------------------------------------------------------------------------------------------
# Dataset selection: prefer dense only when complete
# --------------------------------------------------------------------------------------------
def choose(label, base_name, dense_name, coverage_fn):
    """Return (data, grid_id, info). `coverage_fn(data)` -> number of cells the dataset covers."""
    info = {"dataset": label, "files_considered": [], "rejected": []}
    dense = load(dense_name) if dense_name else None
    if dense is not None:
        cov = coverage_fn(dense)
        info["files_considered"].append({"file": dense_name, "n_cells": cov})
        if cov >= DENSE_MIN_COVERAGE * DENSE_N:
            info.update(file=dense_name, n_cells=cov, grid="dense")
            return dense, "dense", info
        info["rejected"].append(
            f"{dense_name}: covers {cov}/{DENSE_N} cells ({100 * cov / DENSE_N:.0f}%), below the "
            f"{int(100 * DENSE_MIN_COVERAGE)}% completeness bar -- not used")
    base = load(base_name)
    if base is None:
        return None, None, info
    cov = coverage_fn(base)
    info["files_considered"].append({"file": base_name, "n_cells": cov})
    info.update(file=base_name, n_cells=cov, grid="validation")
    return base, "validation", info


def _min_group_cells(d):
    ns = [g.get("n_cells", 0) for g in d.get("groups", {}).values() if g.get("status") == "ok"]
    return min(ns) if ns else 0


def _future_cells(d):
    members = d.get("members", {})
    if not members:
        return 0
    first = next(iter(members.values()))
    ns = [len(g.get("cells", [])) for g in first.get("groups", {}).values()]
    return min(ns) if ns else 0


# --------------------------------------------------------------------------------------------
# Grid assembly
# --------------------------------------------------------------------------------------------
class GridBuilder:
    def __init__(self, gid):
        self.gid = gid
        self.keys = {}      # (lat, lon) -> {"elev": ...}
        self.layers = {}    # layer_id -> {(lat,lon): value}
        self.scen = {}      # name -> structure
        self.layer_meta = {}

    def touch(self, lat, lon, elev=None):
        k = key(lat, lon)
        rec = self.keys.setdefault(k, {"elev": None})
        if elev is not None and rec["elev"] is None:
            rec["elev"] = float(elev)
        return k

    def set(self, layer_id, k, value, meta=None):
        self.layers.setdefault(layer_id, {})[k] = value
        if meta and layer_id not in self.layer_meta:
            self.layer_meta[layer_id] = meta

    def finalise(self):
        # The sampling frame is a rectangle that Armenia fills only ~37% of; cells in Georgia, Azerbaijan,
        # Turkey, Iran or Nakhchivan are dropped here so nothing shown, averaged or recommended is foreign.
        kept = [k for k in self.keys if bool(inside_armenia(k[0], k[1]))]
        self.n_dropped_foreign = len(self.keys) - len(kept)
        order = sorted(kept, key=lambda k: (-k[0], k[1]))
        index = {k: i for i, k in enumerate(order)}
        out = {
            "grid": self.gid, "label": GRIDS[self.gid], "n_cells": len(order),
            "lat": [k[0] for k in order], "lon": [k[1] for k in order],
            "elev": [r(self.keys[k]["elev"], 0) for k in order],
            "layers": {}, "scenarios": {}, "n_cells_dropped_outside_armenia": self.n_dropped_foreign,
        }
        for lid, vals in self.layers.items():
            arr = [None] * len(order)
            for k, v in vals.items():
                if k in index:
                    arr[index[k]] = v
            out["layers"][lid] = arr
        return out, index


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    grids = {g: GridBuilder(g) for g in GRIDS}
    provenance = {}
    engines = {}

    # ---- TOPOHYDRO ---------------------------------------------------------------------------
    topo, gid, info = choose("TOPOHYDRO grid run", "topohydro_grid_run_2019.yaml",
                             "topohydro_grid_run_2019_dense.yaml", lambda d: d.get("n_with_real_output", 0))
    provenance["topohydro"] = info
    pets = []
    if topo:
        G = grids[gid]
        for p in topo["points"]:
            if p.get("status") != "ok":
                continue
            k = G.touch(p["lat"], p["lon"], p["elevation_m"])
            G.set("t_mean", k, r(p["t_mean_c_annual_mean"], 2), {"label": "Mean annual temperature", "unit": "°C", "engine": "TOPOHYDRO"})
            G.set("precip", k, r(p["p_mm_annual_sum"], 0), {"label": "Annual precipitation", "unit": "mm", "engine": "TOPOHYDRO"})
            G.set("gdd", k, r(p["gdd_cumulative_annual"], 0), {"label": "Growing degree days (base 5 °C)", "unit": "°C·d", "engine": "TOPOHYDRO"})
            G.set("late_frost", k, p["late_frost_days"], {"label": "Late-frost days (after budburst GDD)", "unit": "days", "engine": "TOPOHYDRO", "placeholder": "budburst GDD = 200 (placeholder)"})
            G.set("gsl", k, p["growing_season_length_days"], {"label": "Growing-season length", "unit": "days", "engine": "TOPOHYDRO"})
            for pet, v in p["cwd_mm_by_pet_formulation"].items():
                if pet not in pets:
                    pets.append(pet)
                G.set(f"cwd_{pet}", k, r(v, 1), {"label": f"Climatic water deficit ({pet})", "unit": "mm", "engine": "TOPOHYDRO", "pet": pet})
            wsi = p["wsi_by_pet_formulation"].get("pm_fao56")
            G.set("wsi", k, r(wsi, 2), {"label": "Water-stress integral (PM-FAO56)", "unit": "", "engine": "TOPOHYDRO"})
            psi = p["psi_soil_mpa_annual_min_by_pet_formulation"].get("pm_fao56")
            G.set("psi_min", k, r(psi, 2), {"label": "Minimum soil water potential (PM-FAO56)", "unit": "MPa", "engine": "TOPOHYDRO"})
        engines["topohydro"] = {"grid": gid, "n_cells": info["n_cells"], "pet_formulations": pets,
                                "placeholders": topo.get("placeholders")}

    # ---- XYLEM -------------------------------------------------------------------------------
    xy, gid, info = choose("XYLEM mechanistic hazard", "xylem_mechanistic_hazard_2019.yaml",
                           "xylem_mechanistic_hazard_2019_dense.yaml", _min_group_cells)
    provenance["xylem"] = info
    xylem_summary = {}
    if xy:
        G = grids[gid]
        for gname, g in xy["groups"].items():
            if g.get("status") != "ok":
                xylem_summary[gname] = {"status": g.get("status")}
                continue
            for c in g["cells"]:
                k = G.touch(c["lat"], c["lon"], c.get("elevation_m"))
                G.set(f"hmech_{gname}", k, r(c["h_mech_mean"], 5),
                      {"label": f"Hydraulic-failure hazard — {GROUPS[gname]['label']}", "unit": "probability/yr", "engine": "XYLEM", "group": gname})
                G.set(f"hmech_sd_{gname}", k, r(c["h_mech_outer_spread_sd"], 5),
                      {"label": f"Hazard uncertainty (outer-loop sd) — {GROUPS[gname]['short']}", "unit": "", "engine": "XYLEM", "group": gname})
            ins = in_arm(g["cells"])
            xylem_summary[gname] = {"status": "ok", "n_cells": int(ins.sum()),
                                    "mean_h_mech": r(np.mean([c["h_mech_mean"] for c, i in zip(g["cells"], ins) if i]), 5) if ins.any() else None}
        engines["xylem"] = {"grid": gid, "groups": xylem_summary, "outer_draws": xy.get("outer_draws"),
                            "inner_draws": xy.get("inner_draws"),
                            "individual_sd_fraction": xy.get("individual_sd_fraction_of_hyper_sd_placeholder")}
        for gname, g in xy["groups"].items():
            if g.get("status") != "ok":
                engines["xylem"].setdefault("skipped", {})[gname] = g.get("status")

    # ---- REFUGIUM ----------------------------------------------------------------------------
    rf, gid, info = choose("REFUGIUM viability", "refugium_viability_2019.yaml",
                           "refugium_viability_2019_dense.yaml", _min_group_cells)
    provenance["refugium"] = info
    baseline_viab = {}
    refugium_arm = {}
    if rf:
        G = grids[gid]
        for gname, g in rf["groups"].items():
            if g.get("status") != "ok":
                continue
            for c in g["cells"]:
                k = G.touch(c["lat"], c["lon"], c.get("elevation_m"))
                G.set(f"viab_{gname}", k, r(c["viability_ensemble_mean"], 4),
                      {"label": f"Viability (2019, hydraulic survival) — {GROUPS[gname]['label']}", "unit": "probability", "engine": "REFUGIUM", "group": gname})
                G.set(f"refscore_{gname}", k, r(c["refugium_score"], 4),
                      {"label": f"Risk-averse refugium score — {GROUPS[gname]['short']}", "unit": "", "engine": "REFUGIUM", "group": gname})
                G.set(f"robust_{gname}", k, 1 if c["robust_refugium_criterion_a"] else 0,
                      {"label": f"Robust refugium, criterion (a) — {GROUPS[gname]['short']}", "unit": "0/1", "engine": "REFUGIUM", "group": gname, "binary": True})
            ins = in_arm(g["cells"])
            vals = [c["viability_ensemble_mean"] for c, i in zip(g["cells"], ins) if i]
            baseline_viab[gname] = r(np.mean(vals), 4) if vals else None
            refugium_arm[gname] = {"mean_viability": baseline_viab[gname], "n_cells": int(ins.sum()),
                                   "n_robust": int(sum(1 for c, i in zip(g["cells"], ins) if i and c["robust_refugium_criterion_a"]))}
        engines["refugium"] = {"grid": gid, "v_star": rf.get("v_star"), "rho": rf.get("rho"), "lam": rf.get("lam"),
                               "groups": refugium_arm,
                               "scope_note": rf.get("scope_note")}

    # ---- Treeline (today + change) -----------------------------------------------------------
    tc, gid, info = choose("Treeline change", "treeline_change.yaml", "treeline_change_dense.yaml",
                           lambda d: d.get("n_cells", 0))
    provenance["treeline_change"] = info
    treeline_summary = None
    if tc:
        G = grids[gid]
        cell_keys = []
        for c in tc["cells"]:
            k = G.touch(c["lat"], c["lon"], c["elevation_m"])
            cell_keys.append(k)
            z0 = c["potential_treeline_2019_m"]
            G.set("treeline_2019", k, r(z0, 0), {"label": "Potential treeline elevation (2019 climate)", "unit": "m", "engine": "Treeline"})
            G.set("treeline_margin", k, r(z0 - c["elevation_m"], 0), {"label": "Headroom below climatic treeline (treeline − cell elevation)", "unit": "m", "engine": "Treeline"})
            G.set("treeline_above", k, 1 if c["elevation_m"] > z0 else 0, {"label": "Cell above its own climatic treeline (2019)", "unit": "0/1", "engine": "Treeline", "binary": True})
        member_keys = sorted(tc["members"].keys())
        G.scen["treeline_members"] = [{"gcm": tc["members"][m]["gcm"], "ssp": tc["members"][m]["scenario"],
                                       "horizon": tc["members"][m]["horizon"]} for m in member_keys]
        G.scen["treeline_cell_keys"] = cell_keys
        G.scen["treeline_shift"] = [[r(v, 0) for v in tc["members"][m]["shift_m"]] for m in member_keys]
        tl_in = in_arm(tc["cells"])
        tl_summary = {}
        for ssp in ("ssp126", "ssp370", "ssp585"):
            for hz in (2050, 2080, 2100):
                per_gcm = [float(np.nanmean(np.array(m["shift_m"], dtype=float)[tl_in]))
                           for m in tc["members"].values() if m["scenario"] == ssp and int(m["horizon"]) == hz]
                tl_summary[f"{ssp}__{hz}"] = {"ensemble_mean_shift_m": round(float(np.mean(per_gcm)), 1),
                                              "ensemble_min_shift_m": round(float(np.min(per_gcm)), 1),
                                              "ensemble_max_shift_m": round(float(np.max(per_gcm)), 1), "n_gcms": len(per_gcm)}
        treeline_summary = {"grid": gid, "summary": tl_summary, "n_cells_armenia": int(tl_in.sum()),
                            "baseline_check_max_abs_diff_c": tc.get("baseline_check_max_abs_diff_c"),
                            "method": tc["method"], "scope_note": tc["scope_note"],
                            "gamma_k_per_km": r(tc["gamma_growing_season_k_per_km"], 2),
                            "threshold_c": tc["thermal_threshold_c"],
                            "frac_pairs_negative": r(float(np.mean([v < 0 for m in tc["members"].values() for v, i in zip(m["shift_m"], tl_in) if i])), 4)}

    # ---- Future projections --------------------------------------------------------------------
    fp, gid, info = choose("Future projections", "future_projections.yaml", "future_projections_dense.yaml", _future_cells)
    provenance["future_projections"] = info
    scen_summary = None
    if fp:
        G = grids[gid]
        member_keys = sorted(fp["members"].keys())
        first = fp["members"][member_keys[0]]
        gnames = list(first["groups"].keys())
        cell_keys = [G.touch(c["lat"], c["lon"]) for c in first["groups"][gnames[0]]["cells"]]
        G.scen["members"] = [{"gcm": fp["members"][m]["gcm"], "ssp": fp["members"][m]["scenario"],
                              "horizon": int(fp["members"][m]["horizon"])} for m in member_keys]
        G.scen["cell_keys"] = cell_keys
        G.scen["viab"] = {gn: [[r(c["viability_mean"], 4) for c in fp["members"][m]["groups"][gn]["cells"]]
                               for m in member_keys] for gn in gnames}
        G.scen["hmech"] = {gn: [[r(c["h_mech_mean"], 5) for c in fp["members"][m]["groups"][gn]["cells"]]
                                for m in member_keys] for gn in gnames}
        fp_in = in_arm(first["groups"][gnames[0]]["cells"])
        # Precomputed ensemble summary over ARMENIAN cells only: mean over cells, then mean / min / max across GCMs.
        scen_summary = {}
        for gn in gnames:
            scen_summary[gn] = {}
            for ssp in ("ssp126", "ssp370", "ssp585"):
                scen_summary[gn][ssp] = {}
                for hz in (2050, 2080, 2100):
                    vals = [float(np.mean([c["viability_mean"] for c, i in zip(fp["members"][m]["groups"][gn]["cells"], fp_in) if i]))
                            for m in member_keys
                            if fp["members"][m]["scenario"] == ssp and int(fp["members"][m]["horizon"]) == hz]
                    scen_summary[gn][ssp][str(hz)] = {"mean": r(np.mean(vals), 4), "min": r(np.min(vals), 4),
                                                       "max": r(np.max(vals), 4), "n_gcms": len(vals)}
        engines["future_projections"] = {"grid": gid, "n_members": len(member_keys), "method": fp.get("method"),
                                         "gcms": sorted({fp["members"][m]["gcm"] for m in member_keys})}

    # ---- Ecosystem-map validation (78-point grid only) ----------------------------------------
    eco = load("ecosystem_ground_truth_2019.yaml")
    ecosystem = None
    if eco:
        G = grids["validation"]
        for c in eco["cells"]:
            if c.get("outside_real_armenia_raster_extent"):
                continue
            k = G.touch(c["lat"], c["lon"])
            G.set("human_modified", k, r(c["human_modified_fraction"], 3),
                  {"label": "Human-modified land within 500 m (Ecosystem Map)", "unit": "fraction", "engine": "Ecosystem map"})
            for gn, frac in (c.get("forest_cover_fraction_by_group") or {}).items():
                if gn in GROUPS:
                    G.set(f"forest_cover_{gn}", k, r(frac, 3),
                          {"label": f"Mapped forest cover within 500 m — {GROUPS[gn]['short']} (Ecosystem Map)", "unit": "fraction", "engine": "Ecosystem map", "group": gn})
        ecosystem = {"n_cells_total": eco["n_cells_total"], "n_outside": eco["n_cells_outside_armenia_border"],
                     "validation": eco["validation_viability_vs_observed_forest_cover"],
                     "window_radius_m": eco["window_radius_m"], "scope_note": eco["scope_note"]}

    # ---- AEGIS --------------------------------------------------------------------------------
    ae, gid, info = choose("AEGIS portfolio", "aegis_portfolio.yaml", "aegis_portfolio_dense.yaml",
                           lambda d: d.get("n_units", 0))
    provenance["aegis"] = info
    aegis = None
    if ae:
        sys.path.insert(0, str(ROOT / "scripts"))
        try:
            from fit_aegis_portfolio import INTERVENTIONS
            cost = {i["name"]: i["cost_per_ha"] for i in INTERVENTIONS}
        except Exception:  # pragma: no cover
            cost = {}
        sweep = [{"budget_usd": float(k.replace("$", "").replace(",", "")), **{a: r(b, 2) if a != "n_units_planted" else b for a, b in v.items()}}
                 for k, v in ae["budget_sweep"].items()]
        fr = ae["lambda_frontier_at_representative_budget"]
        aegis = {"grid": gid, "scenario_source": ae["scenario_source"], "n_scenarios": ae["n_scenarios"],
                 "n_units": ae["n_units"], "value_per_ha_year_usd": ae["value_per_ha_year_usd"],
                 "scope_note": ae["scope_note"], "budget_sweep": sweep,
                 "frontier": {"lambdas": fr.get("lambdas"), "expected": [r(v, 2) for v in fr.get("expected", [])],
                              "cvar": [r(v, 2) for v in fr.get("cvar", [])],
                              "price_of_robustness": r(fr.get("price_of_robustness"), 4)},
                 "options": [{**o, "cost_per_ha_usd": cost.get(o["intervention"]),
                              "group_label": GROUPS.get(o["group"], {}).get("label", o["group"])}
                             for o in ae["option_labels"]]}

    # ---- MERISTEM / MNEME / variogram ---------------------------------------------------------
    mer = load("meristem_adult_niche.yaml")
    meristem = None
    if mer:
        meristem = {"features": mer["features"], "cwd_caveat": mer["cwd_caveat"], "pooling_note": mer["pooling_note"],
                    "background_type": mer.get("background_type"), "species": {
                        sp: {"status": v["status"], "n_presence": v.get("n_presence_after_nan_drop"),
                             "boyce_mean": (v.get("spatial_block_cv") or {}).get("boyce_index_mean"),
                             "boyce_folds": (v.get("spatial_block_cv") or {}).get("boyce_index_per_fold")}
                        for sp, v in mer["species"].items()}}
    mn, gid, info = choose("MNEME hazard panel", "mneme_hazard_panel_2010_2019.yaml",
                           "mneme_hazard_panel_2010_2019_dense.yaml", lambda d: d.get("n_points_valid_kndvi", 0))
    provenance["mneme"] = info
    mneme = {k: mn[k] for k in ("grid", "n_points", "n_points_valid_kndvi", "n_person_years", "n_events",
                                "n_spatial_blocks", "status", "scope_note") if k in mn} if mn else None
    vg, gid, info = choose("Variogram", "variogram.yaml", "variogram_dense.yaml", lambda d: d.get("n_points", 0))
    provenance["variogram"] = info
    variogram = vg

    # ---- write grid files ---------------------------------------------------------------------
    grid_files, layer_catalogue = {}, {}
    for gid, G in grids.items():
        if not G.keys:
            continue
        data, index = G.finalise()
        # remap scenario arrays (stored in their own cell order) onto this grid's cell order
        scen = {}
        if "members" in G.scen:
            n = len(index)
            pairs = [(index[k], j) for j, k in enumerate(G.scen["cell_keys"]) if k in index]

            def remap(mat):
                out = []
                for row in mat:
                    full = [None] * n
                    for pos, j in pairs:
                        full[pos] = row[j]
                    out.append(full)
                return out
            scen["members"] = G.scen["members"]
            scen["viab"] = {gn: remap(m) for gn, m in G.scen["viab"].items()}
            scen["hmech"] = {gn: remap(m) for gn, m in G.scen["hmech"].items()}
        if "treeline_members" in G.scen:
            n = len(index)
            pairs = [(index[k], j) for j, k in enumerate(G.scen["treeline_cell_keys"]) if k in index]
            full_rows = []
            for row in G.scen["treeline_shift"]:
                full = [None] * n
                for pos, j in pairs:
                    full[pos] = row[j]
                full_rows.append(full)
            scen["treeline_members"] = G.scen["treeline_members"]
            scen["treeline_shift"] = full_rows
        data["scenarios"] = scen
        fname = f"grid_{gid}.json"
        dump(data, OUT / fname)
        grid_files[gid] = {"file": fname, "n_cells": data["n_cells"], "label": GRIDS[gid], "n_dropped_outside_armenia": data["n_cells_dropped_outside_armenia"],
                           "bytes": (OUT / fname).stat().st_size,
                           "has_scenarios": bool(scen.get("members")), "has_treeline_change": bool(scen.get("treeline_members"))}
        for lid, meta in G.layer_meta.items():
            layer_catalogue[lid] = {**meta, "grid": gid, "grid_label": GRIDS[gid]}

    # ---- references from the dataset manifests ------------------------------------------------
    refs, seen = [], set()
    for mf in sorted(MANIFESTS.glob("*.yaml")):
        try:
            doc = yaml.safe_load(open(mf))
        except Exception:
            continue
        if not isinstance(doc, list):
            continue
        for e in doc:
            if not isinstance(e, dict) or not e.get("citation"):
                continue
            cite = str(e["citation"]).strip()
            if cite.startswith("See per-source"):
                continue
            sig = (e.get("variable"), cite)
            if sig in seen:
                continue
            seen.add(sig)
            refs.append({"dataset": e.get("variable"), "source": e.get("source"), "version": e.get("version"),
                         "citation": cite, "url": e.get("url"), "license": e.get("license"),
                         "accessed": str(e.get("access_date")) if e.get("access_date") else None,
                         "manifest": mf.name})
    eco_manifest = yaml.safe_load(open(MANIFESTS / "ecosystem_map_armenia.yaml")) if (MANIFESTS / "ecosystem_map_armenia.yaml").exists() else None
    if eco_manifest:
        s = eco_manifest["source"]
        refs.append({"dataset": s["title"], "source": "; ".join(s["publishers"]), "version": "2026-09-18",
                     "citation": f'{s["title"]} — {"; ".join(s["publishers"])}. {s["project"].strip()}',
                     "url": s["app_url"], "license": "CC BY 4.0 (per the app page; the zip's own README carries an unfilled licence placeholder)",
                     "accessed": s["downloaded"], "manifest": "ecosystem_map_armenia.yaml"})

    manifest = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "groups": {gn: {**meta, "rooting_depth_m": ROOTING_DEPTH_M[gn]} for gn, meta in GROUPS.items()},
        "grids": grid_files,
        "layers": layer_catalogue,
        "provenance": provenance,
        "engines": engines,
        "treeline": treeline_summary,
        "scenario_summary": scen_summary,
        "baseline_viability_2019": baseline_viab,
        "ecosystem_validation": ecosystem,
        "aegis": aegis,
        "meristem": meristem,
        "mneme": mneme,
        "variogram": variogram,
        "references": refs,
        "dense_min_coverage": DENSE_MIN_COVERAGE,
    }
    dump(manifest, OUT / "manifest.json")
    total = sum(f["bytes"] for f in grid_files.values()) + (OUT / "manifest.json").stat().st_size
    print(f"wrote {len(grid_files)} grid file(s) + manifest.json  ({total / 1e3:.0f} kB total); "
          f"{len(layer_catalogue)} layers; {len(refs)} references")
    for name, info in provenance.items():
        grid = info.get("grid", "-")
        print(f"  {name:20s} -> {info.get('file', 'MISSING'):52s} grid={grid}")
        for rej in info.get("rejected", []):
            print(f"      REJECTED {rej}")


if __name__ == "__main__":
    build()
