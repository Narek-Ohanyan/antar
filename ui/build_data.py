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
import re
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
# group-specific (rooting-depth-dependent) water-stress fields: (field in yaml, layer-id prefix, label, unit, decimals)
HYDRO_FIELDS = [("cwd_mm", "cwd", "Climatic water deficit", "mm", 1), ("wsi", "wsi", "Water-stress integral", "", 2),
                ("psi_min_mpa", "psimin", "Minimum soil water potential", "MPa", 2)]
# For colouring change maps: which direction is favourable for the forest. Matched on the layer-id prefix.
GOOD_DIRECTION = [("viab_", "high"), ("pviab_", "high"), ("robust_", "high"), ("robustfull_", "high"), ("aoa_", "high"), ("buffer_", "high"), ("refscore_", "high"), ("hmech_sd", None),
                  ("hmech_", "low"), ("cwd_", "low"), ("wsi", "low"), ("psimin_", "high"), ("psi_min", "high"),
                  ("treeline_margin", "high"), ("treeline_shift", "high"), ("treeline_2019", "high")]
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
SEASONAL_FORCING_FROM = "2026-10-06"      # date the seasonal ISIMIP atmosphere replaced the constant annual means


def forcing_current(doc, strict=False):
    """True if a result file was produced with the seasonal, scenario-dependent atmosphere. New files say so in `atmosphere`;
    older ones are judged by run_date (strict: a missing `atmosphere` key means superseded)."""
    if doc is None:
        return False
    if "atmosphere" in doc:
        return str(doc["atmosphere"]).startswith("seasonal")
    return (not strict) and str(doc.get("run_date", "")) >= SEASONAL_FORCING_FROM


def choose(label, base_name, dense_name, coverage_fn, gate=None):
    """Return (data, grid_id, info). `coverage_fn(data)` -> number of cells the dataset covers. `gate(doc)` -> bool: a file that
    fails it is WITHHELD (listed in info["withheld"], shown on the Status page, never drawn)."""
    info = {"dataset": label, "files_considered": [], "rejected": [], "withheld": []}
    dense = load(dense_name) if dense_name else None
    if dense is not None and gate is not None and not gate(dense):
        info["withheld"].append(f"{dense_name}: computed with the earlier constant annual-mean wind, radiation and humidity; withheld until recomputed")
        dense = None
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
    if base is not None and gate is not None and not gate(base):
        info["withheld"].append(f"{base_name}: computed with the earlier constant annual-mean wind, radiation and humidity; withheld until recomputed")
        base = None
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


def build_map_plan(layers, series_ids):
    """Every map the UI is meant to offer, with whether its 2019 and scenario results exist yet and which script
    produces what is missing. This is what the Status page shows, so 'no data' is never unexplained."""
    FP, XY, RF, TOPO, TL = ("run_future_projections.py", "fit_xylem_mechanistic_hazard.py", "fit_refugium_viability.py",
                            "run_topohydro_grid.py", "compute_treeline_change.py")
    CR = "compute_refugium_criteria.py"
    plan = []

    def engine_of(lid):
        pre = lid.split("_")[0]
        return {"hmech": "XYLEM", "viab": "REFUGIUM", "pviab": "REFUGIUM", "refscore": "REFUGIUM", "robust": "REFUGIUM", "robustfull": "REFUGIUM", "buffer": "REFUGIUM", "aoa": "REFUGIUM", "treeline": "Treeline"}.get(pre, "TOPOHYDRO")

    def add(lid, label, group, base_from, scen_from):
        good = next((d for pre, d in GOOD_DIRECTION if lid.startswith(pre)), None)
        plan.append({"id": lid, "label": label, "group": group, "engine": engine_of(lid), "good": good,
                     "baseline": lid in layers, "scenario": lid in series_ids,
                     "baseline_from": base_from, "scenario_from": scen_from})
    for lid, label in [("t_mean", "Mean annual temperature"), ("precip", "Annual precipitation"), ("gdd", "Growing degree days"),
                       ("gsl", "Growing-season length")]:
        add(lid, label, None, TOPO, FP)
    for lid, label in [("cwd_pm_fao56", "Climatic water deficit (pm_fao56, generic rooting depth)"),
                       ("cwd_priestley_taylor", "Climatic water deficit (priestley_taylor, generic rooting depth)"),
                       ("cwd_energy_only", "Climatic water deficit (energy_only, generic rooting depth)"),
                       ("wsi", "Water-stress integral (PM-FAO56, generic rooting depth)"),
                       ("psi_min", "Minimum soil water potential (PM-FAO56, generic rooting depth)")]:
        add(lid, label, None, TOPO, FP)
    for g, meta in GROUPS.items():
        short = meta["short"]
        for lid, label, bf in [("cwd", "Climatic water deficit", RF), ("wsi", "Water-stress integral", RF), ("psimin", "Minimum soil water potential", RF),
                               ("hmech", "Hydraulic-failure hazard", XY), ("hmech_sd", "Hazard uncertainty (trait-knowledge spread)", XY),
                               ("viab", "Viability", RF), ("pviab", "Probability that viability stays above the threshold", RF),
                               ("refscore", "Risk-averse refugium score", RF), ("robust", "Robust refugium (criterion a)", RF),
                               ("buffer", "Buffer index (criterion b)", CR), ("aoa", "Inside the area of applicability (criterion c)", CR),
                               ("robustfull", "Robust refugium, all three criteria", CR)]:
            add(f"{lid}_{g}", f"{label} — {short}", g, bf, FP)
    for lid, label in [("treeline_2019", "Potential treeline elevation"), ("treeline_margin", "Headroom below climatic treeline"),
                       ("treeline_above", "Above own climatic treeline")]:
        add(lid, label, None, TL, TL)
    plan.append({"id": "treeline_shift", "label": "Treeline shift vs 2019", "group": None, "engine": "Treeline", "good": "high", "baseline": None,
                 "scenario": "treeline_shift" in series_ids, "baseline_from": None, "scenario_from": TL})
    return plan


def add_water_layers(G, k, p, pets):
    """Climatic water deficit (3 PET formulations), water-stress integral and minimum soil potential for one TOPOHYDRO point."""
    for pet, v in p["cwd_mm_by_pet_formulation"].items():
        if pet not in pets:
            pets.append(pet)
        G.set(f"cwd_{pet}", k, r(v, 1), {"label": f"Climatic water deficit ({pet}, generic rooting depth)", "unit": "mm", "engine": "TOPOHYDRO", "pet": pet})
    G.set("wsi", k, r(p["wsi_by_pet_formulation"].get("pm_fao56"), 2), {"label": "Water-stress integral (PM-FAO56, generic rooting depth)", "unit": "", "engine": "TOPOHYDRO"})
    G.set("psi_min", k, r(p["psi_soil_mpa_annual_min_by_pet_formulation"].get("pm_fao56"), 2), {"label": "Minimum soil water potential (PM-FAO56, generic rooting depth)", "unit": "MPa", "engine": "TOPOHYDRO"})


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
        # temperature, rain, degree days, season length and frost do not depend on wind, radiation or humidity; the water balance does
        water_ok = forcing_current(topo)
        if not water_ok:
            info["withheld"].append(f"{info['file']}: climatic water deficit, water-stress integral and soil potential were computed with the earlier constant annual-mean wind, radiation and humidity; withheld until recomputed")
        G = grids[gid]
        for p in topo["points"]:
            if p.get("status") != "ok":
                continue
            k = G.touch(p["lat"], p["lon"], p["elevation_m"])
            G.set("t_mean", k, r(p["t_mean_c_annual_mean"], 2), {"label": "Mean annual temperature", "unit": "°C", "engine": "TOPOHYDRO"})
            G.set("precip", k, r(p["p_mm_annual_sum"], 0), {"label": "Annual precipitation", "unit": "mm", "engine": "TOPOHYDRO"})
            G.set("gdd", k, r(p["gdd_cumulative_annual"], 0), {"label": "Growing degree days (base 5 °C)", "unit": "°C·d", "engine": "TOPOHYDRO"})
            G.set("gsl", k, p["growing_season_length_days"], {"label": "Growing-season length", "unit": "days", "engine": "TOPOHYDRO"})
            if water_ok:
                add_water_layers(G, k, p, pets)
        if not water_ok:
            # the chosen (dense) file predates the seasonal forcing: take the water-balance layers from the other file if it is current
            alt = load("topohydro_grid_run_2019.yaml" if gid == "dense" else "topohydro_grid_run_2019_dense.yaml")
            alt_gid = "validation" if gid == "dense" else "dense"
            if alt is not None and forcing_current(alt):
                for p in alt["points"]:
                    if p.get("status") == "ok":
                        add_water_layers(grids[alt_gid], grids[alt_gid].touch(p["lat"], p["lon"], p["elevation_m"]), p, pets)
                info["withheld"].append(f"water-balance layers taken from the {alt_gid} grid file, which is current")
        engines["topohydro"] = {"grid": gid, "n_cells": info["n_cells"], "pet_formulations": pets,
                                "placeholders": topo.get("placeholders")}

    # ---- XYLEM -------------------------------------------------------------------------------
    xy, gid, info = choose("XYLEM mechanistic hazard", "xylem_mechanistic_hazard_2019.yaml",
                           "xylem_mechanistic_hazard_2019_dense.yaml", _min_group_cells, gate=forcing_current)
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
                           "refugium_viability_2019_dense.yaml", _min_group_cells, gate=forcing_current)
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
                if "p_viable" in c:
                    G.set(f"pviab_{gname}", k, r(c["p_viable"], 3),
                          {"label": f"Probability that viability ≥ {rf.get('v_star')} — {GROUPS[gname]['short']}", "unit": "probability", "engine": "REFUGIUM", "group": gname})
                for fld, lid, label, unit, nd in HYDRO_FIELDS:
                    if fld in c:
                        G.set(f"{lid}_{gname}", k, r(c[fld], nd),
                              {"label": f"{label} — {GROUPS[gname]['short']} rooting depth", "unit": unit, "engine": "REFUGIUM", "group": gname})
            ins = in_arm(g["cells"])
            vals = [c["viability_ensemble_mean"] for c, i in zip(g["cells"], ins) if i]
            baseline_viab[gname] = r(np.mean(vals), 4) if vals else None
            refugium_arm[gname] = {"mean_viability": baseline_viab[gname], "n_cells": int(ins.sum()),
                                   "n_robust": int(sum(1 for c, i in zip(g["cells"], ins) if i and c["robust_refugium_criterion_a"]))}
        crit = load("refugium_criteria_dense.yaml") if gid == "dense" else None
        if crit:
            for gname, cg in crit["groups"].items():
                if cg["status"] != "ok":
                    continue                       # no area of applicability can be built for this group, so there is no layer to draw (the Models and Method pages say so)
                short = GROUPS[gname]["short"]
                note = ""
                for c in cg["cells"]:
                    k = G.touch(c["lat"], c["lon"])
                    if c["buffer_index"] is not None:
                        G.set(f"buffer_{gname}", k, r(c["buffer_index"], 2), {"label": f"Buffer index (2019), criterion (b) — {short}", "unit": "SD", "engine": "REFUGIUM", "group": gname})
                    G.set(f"aoa_{gname}", k, 1 if c["inside_aoa"] else 0, {"label": f"Inside the area of applicability (2019), criterion (c) — {short}{note}", "unit": "0/1", "engine": "REFUGIUM", "group": gname, "binary": True})
                    G.set(f"robustfull_{gname}", k, 1 if c["robust_full"] else 0, {"label": f"Robust refugium, all three criteria (2019) — {short}{note}", "unit": "0/1", "engine": "REFUGIUM", "group": gname, "binary": True})
            refugium_arm["criteria"] = {gname: {"status": cg["status"], "n_aoa_training_cells": cg["n_aoa_training_nodes"], **cg["baseline"]} for gname, cg in crit["groups"].items()}
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
        G.scen["tl_series"] = {"treeline_shift": [[r(v, 0) for v in tc["members"][m]["shift_m"]] for m in member_keys],
                               "treeline_2019": [[r(v, 0) for v in tc["members"][m]["potential_treeline_m"]] for m in member_keys]}
        elev_tl = [c["elevation_m"] for c in tc["cells"]]
        G.scen["tl_series"]["treeline_margin"] = [[r(z - e, 0) for z, e in zip(tc["members"][m]["potential_treeline_m"], elev_tl)] for m in member_keys]
        G.scen["tl_series"]["treeline_above"] = [[1 if e > z else 0 for z, e in zip(tc["members"][m]["potential_treeline_m"], elev_tl)] for m in member_keys]
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
    fp, gid, info = choose("Future projections", "future_projections.yaml", "future_projections_dense.yaml", _future_cells,
                           gate=lambda d: forcing_current(d, strict=True))
    provenance["future_projections"] = info
    scen_summary = None
    fp_in = []
    if fp:
        G = grids[gid]
        member_keys = sorted(fp["members"].keys())
        first = fp["members"][member_keys[0]]
        gnames = list(first["groups"].keys())
        cell_keys = [G.touch(c["lat"], c["lon"]) for c in first["groups"][gnames[0]]["cells"]]
        G.scen["members"] = [{"gcm": fp["members"][m]["gcm"], "ssp": fp["members"][m]["scenario"],
                              "horizon": int(fp["members"][m]["horizon"])} for m in member_keys]
        G.scen["cell_keys"] = cell_keys
        series = {}
        # (yaml field, layer-id prefix, decimals); a series is exported only if EVERY member has the field
        for fld, lid, nd in [("viability_mean", "viab", 4), ("h_mech_mean", "hmech", 5), ("p_viable", "pviab", 3),
                             ("cwd_mm", "cwd", 1), ("wsi", "wsi", 2), ("psi_min_mpa", "psimin", 2),
                             ("h_mech_sd", "hmech_sd", 5), ("refugium_score", "refscore", 4), ("robust_criterion_a", "robust", 0)]:
            for gn in gnames:
                rows = [fp["members"][m]["groups"][gn]["cells"] for m in member_keys]
                if all(fld in c for cells in rows for c in cells):
                    conv = (lambda v: int(bool(v))) if fld == "robust_criterion_a" else (lambda v, nd=nd: r(v, nd))
                    series[f"{lid}_{gn}"] = [[conv(c[fld]) for c in cells] for cells in rows]
        crit_s = load("refugium_criteria_dense.yaml") if gid == "dense" else None
        if crit_s:
            for gname, cg in crit_s["groups"].items():
                if gname not in gnames or cg["status"] != "ok":
                    continue
                pos = {key(c["lat"], c["lon"]): i for i, c in enumerate(cg["cells"])}
                idx = [pos.get(ck) for ck in cell_keys]
                for lid, field in (("robustfull", "robust_full"), ("aoa", "inside_aoa")):
                    series[f"{lid}_{gname}"] = [[None if i is None else cg["members"][m][field][i] for i in idx] for m in member_keys]
        if all("climate" in fp["members"][m] for m in member_keys):
            for fld, lid, nd in [("t_mean_c", "t_mean", 2), ("precip_mm", "precip", 0), ("gdd", "gdd", 0),
                                 ("gsl_days", "gsl", 0)]:
                series[lid] = [[r(c[fld], nd) for c in fp["members"][m]["climate"]] for m in member_keys]
        if all("generic" in fp["members"][m] for m in member_keys):
            for pet in ("pm_fao56", "priestley_taylor", "energy_only"):
                series[f"cwd_{pet}"] = [[r(c["cwd_mm_by_pet"].get(pet), 1) for c in fp["members"][m]["generic"]] for m in member_keys]
            series["wsi"] = [[r(c["wsi"], 2) for c in fp["members"][m]["generic"]] for m in member_keys]
            series["psi_min"] = [[r(c["psi_min_mpa"], 2) for c in fp["members"][m]["generic"]] for m in member_keys]
        G.scen["fp_series"] = series
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
                           lambda d: d.get("n_units", 0), gate=forcing_current)       # built from the scenario viabilities, so it inherits their forcing
    provenance["aegis"] = info
    aegis = None
    if ae:
        sys.path.insert(0, str(ROOT / "scripts"))
        try:
            from fit_aegis_portfolio import INTERVENTIONS
            cost = {i["name"]: i["cost_per_ha"] for i in INTERVENTIONS}
        except Exception:  # pragma: no cover
            cost = {}
        def sweep_rows(src):
            rows = []
            for k, v in src.items():
                row = {"budget_usd": float(k.replace("$", "").replace(",", ""))}
                for a, val in v.items():
                    if a == "selected":
                        continue                       # the treated cells go to data/aegis_units.json, loaded only by the Decision page
                    row[a] = val if (isinstance(val, (dict, str, bool)) or a == "n_units_planted") else r(val, 6 if a == "mip_gap" else 2)
                rows.append(row)
            return rows
        sweep = sweep_rows(ae["budget_sweep"])
        diversity = {cap: sweep_rows(src) for cap, src in (ae.get("diversity_sweep") or {}).items()}
        fr = ae["lambda_frontier_at_representative_budget"]
        aegis = {"grid": gid, "scenario_source": ae["scenario_source"], "n_scenarios": ae["n_scenarios"],
                 "n_units": ae["n_units"], "n_eligible_units": ae.get("n_eligible_units"), "unit_area_ha": ae.get("unit_area_ha"),
                 "eligible_area_ha": ae.get("eligible_area_ha"), "representative_budget_usd": ae.get("representative_budget_usd", 45_000_000.0),
                 "value_per_ha_year_usd": ae["value_per_ha_year_usd"], "incremental_share": ae.get("incremental_share"),
                 "net_value_per_ha_year_usd": ae.get("net_value_per_ha_year_usd"),
                 "scope_note": ae["scope_note"], "budget_sweep": sweep,
                 "cost_table": [{k: c.get(k) for k in ("name", "cost_per_ha", "cost_low", "cost_high", "cost_basis", "source", "role")} for c in ae.get("cost_table", [])],
                 "not_ranked": ae.get("not_ranked", []), "niche": ae.get("niche"), "not_ranked_groups": ae.get("not_ranked_groups", []), "n_candidate_cells_with_a_supported_group": ae.get("n_candidate_cells_with_a_supported_group"), "group_caps": ae.get("group_caps", []),
                 "diversity_sweep": diversity, "diversity_note": ae.get("diversity_note"), "min_open_share": ae.get("min_open_share"), "n_eligible_before_land_cover": ae.get("n_eligible_before_land_cover"),
                 "planting_cost_sensitivity": ae.get("planting_cost_sensitivity"),
                 "frontier": {"lambdas": fr.get("lambdas"), "expected": [r(v, 2) for v in fr.get("expected", [])],
                              "cvar": [r(v, 2) for v in fr.get("cvar", [])],
                              "price_of_robustness": r(fr.get("price_of_robustness"), 4), "resolution_usd": r(fr.get("resolution_usd"), 0),
                              "solve_optimal": fr.get("solve_optimal"), "solve_mip_gap": fr.get("solve_mip_gap")},
                 "options": [{**o, "cost_per_ha_usd": cost.get(o["intervention"]),
                              "group_label": GROUPS.get(o["group"], {}).get("label", o["group"])}
                             for o in ae["option_labels"]]}

    units_path = OUT / "aegis_units.json"
    if aegis and ae.get("units"):
        dump({"units": ae["units"], "intervention_ids": ae["intervention_ids"],
              "budgets": [{"budget_usd": float(k.replace("$", "").replace(",", "")), "selected": v.get("selected", [])} for k, v in ae["budget_sweep"].items()],
              "capped": {cap: [{"budget_usd": float(k.replace("$", "").replace(",", "")), "selected": v.get("selected", [])} for k, v in src.items()]
                         for cap, src in (ae.get("diversity_sweep") or {}).items()}}, units_path)
        aegis["units_file"] = "aegis_units.json"
    elif units_path.exists():
        units_path.unlink()                        # never leave the cells of a withheld or older portfolio behind

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
        nn = load("meristem_niche_nodes_dense.yaml")
        if nn:
            meristem["min_boyce_gate"] = nn["min_boyce_gate"]
            meristem["presence_omission"] = nn["presence_omission"]
            for sp, v in nn["species"].items():
                if sp in meristem["species"]:
                    meristem["species"][sp]["passes_gate"] = v["passes_gate"]
            meristem["group_support"] = {g: {"status": e["status"], "species_used": e.get("species_used", []), "n_supported_baseline": e.get("n_supported_baseline"), "n_nodes": e["n_nodes"],
                                              "reason": e.get("reason"), "by_path_and_horizon": e.get("n_supported_by_path_and_horizon")} for g, e in nn["groups"].items()}
    mn, vit = load("mneme_hazard_panel_2010_2019_dense.yaml"), load("mneme_vitality_response_dense.yaml")
    provenance["mneme"] = {"dataset": "MNEME forest-pixel panel", "file": "mneme_hazard_panel_2010_2019_dense.yaml" if mn else None, "grid": "forest pixels (up to 25 per climate cell)", "rejected": [], "withheld": []}
    mneme = None
    if mn:
        mneme = {k: mn[k] for k in ("grid", "design", "panel_years", "n_points", "n_pixels", "n_cells", "n_person_years", "n_events", "n_events_if_flags_ignored", "n_pixels_with_event",
                                    "n_cells_with_event", "pixels_per_cell_cap", "n_qualifying_pixels", "event_rate_per_pixel_year", "event_rate_ci95", "year_effect_removed", "mirror_check", "forest_class_pixels", "min_events_per_predictor",
                                    "design_columns", "rule", "status", "status_detail", "decline_to_rise_ratio", "min_decline_to_rise_ratio", "scope_note", "glm_coefficients", "oof_auc_stacked") if k in mn}
        if vit:
            def rd(o):
                if isinstance(o, float):
                    return round(o, 4)
                if isinstance(o, dict):
                    return {k: rd(v) for k, v in o.items()}
                if isinstance(o, list):
                    return [rd(v) for v in o]
                return o
            mneme["vitality"] = rd({k: vit[k] for k in ("what", "n_pixel_years", "n_pixels", "n_cells", "n_excluded_recent_disturbance", "anomaly_sd", "model_a", "model_b", "model_c_height",
                                                        "leave_one_year_out", "sensitivity", "response_by_anomaly_bin", "xylem_check") if k in vit})
    up = load("uncertainty_partition_dense.yaml")
    uncertainty = None
    if up:
        def med(d):
            return {k: (None if v is None else r(v["median"], 3)) for k, v in d.items()}
        uncertainty = {"design": up["design"], "method": up["method"], "internal_variability_note": up["internal_variability_note"], "grid": "dense",
                       "metrics": {mid: {"label": m["label"], "unit": m["unit"], "n_nodes_used": m["n_nodes_used"], "n_nodes": m["n_nodes"], "pooled": med(m["pooled"]),
                                         "by_horizon": {str(h): med(v) for h, v in m["by_horizon"].items()}, "ensemble_sd": r(m["ensemble_sd"], 3)} for mid, m in up["metrics"].items()}}
    provenance["uncertainty_partition"] = {"dataset": "Ensemble uncertainty partition", "file": "uncertainty_partition_dense.yaml" if up else None, "grid": "dense", "rejected": [], "withheld": []}
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
            scen["fp_series"] = {qid: remap(m) for qid, m in G.scen["fp_series"].items()}
        if "treeline_members" in G.scen:
            n = len(index)
            pairs = [(index[k], j) for j, k in enumerate(G.scen["treeline_cell_keys"]) if k in index]
            def remap_tl(mat):
                rows = []
                for row in mat:
                    full = [None] * n
                    for pos, j in pairs:
                        full[pos] = row[j]
                    rows.append(full)
                return rows
            scen["treeline_members"] = G.scen["treeline_members"]
            scen["tl_series"] = {qid: remap_tl(m) for qid, m in G.scen["tl_series"].items()}
        data["scenarios"] = scen
        fname = f"grid_{gid}.json"
        dump(data, OUT / fname)
        grid_files[gid] = {"file": fname, "n_cells": data["n_cells"], "label": GRIDS[gid], "n_dropped_outside_armenia": data["n_cells_dropped_outside_armenia"],
                           "bytes": (OUT / fname).stat().st_size,
                           "has_scenarios": bool(scen.get("members")), "has_treeline_change": bool(scen.get("treeline_members"))}
        for lid, meta in G.layer_meta.items():
            good = next((d for pre, d in GOOD_DIRECTION if lid.startswith(pre)), None)
            layer_catalogue[lid] = {**meta, "grid": gid, "grid_label": GRIDS[gid], "good": good,
                                    # sampled-window layers are 500 m samples at the nodes, not fields: the real raster is the map
                                    "map": meta.get("engine") != "Ecosystem map"}

    series_ids = set()
    for gd in (json.loads((OUT / f["file"]).read_text()) for f in grid_files.values()):
        series_ids |= set(gd["scenarios"].get("fp_series", {})) | set(gd["scenarios"].get("tl_series", {}))
    map_plan = build_map_plan(layer_catalogue, series_ids)

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

    # photographs on the home page: every image file must have a credit entry (author, licence, source), and every entry a file
    photo_dir = Path(__file__).resolve().parent / "assets" / "photos"
    photos = json.loads((photo_dir / "credits.json").read_text()) if (photo_dir / "credits.json").exists() else []
    for ph in photos:
        for w in (800, 1600):
            if not (photo_dir / f"{ph['file']}-{w}.jpg").exists():
                raise SystemExit(f"photo {ph['file']}-{w}.jpg listed in credits.json is missing")

    manifest = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "groups": {gn: {**meta, "rooting_depth_m": ROOTING_DEPTH_M[gn]} for gn, meta in GROUPS.items()},
        "grids": grid_files,
        "layers": layer_catalogue,
        "map_plan": map_plan,
        "provenance": provenance,
        "engines": engines,
        "treeline": treeline_summary,
        "scenario_summary": scen_summary,
        "node_counts": {"viability": int(sum(fp_in)) if scen_summary else None, "treeline": treeline_summary["n_cells_armenia"] if treeline_summary else None,
                        "aegis": aegis["n_units"] if aegis else None},          # the number of model nodes each headline family is taken over (they differ: terrain alone reaches nodes the water balance cannot)
        "baseline_viability_2019": baseline_viab,
        "ecosystem_validation": ecosystem,
        "aegis": aegis,
        "meristem": meristem,
        "mneme": mneme,
        "uncertainty": uncertainty,
        "variogram": variogram,
        "references": refs,
        "photos": photos,
        "dense_min_coverage": DENSE_MIN_COVERAGE,
    }
    dump(manifest, OUT / "manifest.json")
    total = sum(f["bytes"] for f in grid_files.values()) + (OUT / "manifest.json").stat().st_size
    print(f"wrote {len(grid_files)} grid file(s) + manifest.json  ({total / 1e3:.0f} kB total); "
          f"{len(layer_catalogue)} layers; {len(refs)} references")
    print(f"index.html stamped with build {stamp_index()}")
    for name, info in provenance.items():
        grid = info.get("grid", "-")
        print(f"  {name:20s} -> {info.get('file', 'MISSING'):52s} grid={grid}")
        for rej in info.get("rejected", []):
            print(f"      REJECTED {rej}")



def stamp_index():
    """Give every asset URL in index.html a content hash, so a browser can never keep serving a stale script, style or data
    file after they change (a cached app.js showed an outdated warning after it had been corrected)."""
    import hashlib
    ui = Path(__file__).resolve().parent
    files = ([ui / n for n in ("styles.css", "charts.js", "raster.js", "interp.js", "map.js", "place.js", "math.js", "home.js", "decision.js", "palette.js", "app.js")] + sorted((ui / "data").glob("*.json"))
             + sorted((ui / "assets" / "map").glob("*")) + [ui / "assets" / "architecture.svg", ui / "assets" / "relief.svg"] + sorted((ui / "assets" / "photos").glob("*.*")) + sorted((ui / "assets" / "katex").glob("*.*")))
    h = hashlib.sha1()
    for f in files:
        if not f.is_file():
            continue
        h.update(f.name.encode())
        if f.name == "manifest.json":                       # drop the build timestamp so an unchanged result keeps its hash
            m = json.loads(f.read_text()); m.pop("generated", None)
            h.update(json.dumps(m, sort_keys=True).encode())
        else:
            h.update(f.read_bytes())
    build = h.hexdigest()[:10]
    idx = ui / "index.html"
    html = idx.read_text()
    new = re.sub(r"\?v=[A-Za-z0-9]+", f"?v={build}", html)
    new = re.sub(r'window\.ANTAR_BUILD = "[A-Za-z0-9]+"', f'window.ANTAR_BUILD = "{build}"', new)
    if new != html:
        idx.write_text(new)
    return build


if __name__ == "__main__":
    build()
