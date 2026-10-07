"""Removes from the dense-grid results every node that TOPOHYDRO did not accept (status other than ``ok``), and recomputes the summary fields.

Why: three dense nodes had clay data but NaN soil hydraulic parameters, so their water balance was NaN and the hydraulic engine, for which
``NaN >= 1`` is False, reported a hazard of exactly 0 (viability 1.0). The cause is fixed in run_topohydro_grid.py (``soil_complete``), and
TOPOHYDRO's own output already marks the nodes ``skipped_nonfinite_water_balance``; this script brings the XYLEM, REFUGIUM and scenario outputs, which
were computed with the old mask, into line without recomputing hours of Monte Carlo. Nodes are independent and every seed is the node's grid index, so
the remaining nodes are exactly what a rerun would give. Idempotent; the removed nodes are recorded in the file (``dropped_nodes``).

    python3 scripts/drop_invalid_nodes.py [--files xylem refugium scenarios]
"""
import sys
from pathlib import Path

import numpy as np
import yaml

FITTED = Path(__file__).resolve().parent.parent / "configs" / "fitted"
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
DUMPER = getattr(yaml, "CSafeDumper", yaml.SafeDumper)


def key(c):
    return (round(float(c["lat"]), 5), round(float(c["lon"]), 5))


def invalid_nodes(topohydro: dict) -> set:
    return {key(p) for p in topohydro["points"] if p.get("status") != "ok"}


def drop(cells: list, bad: set) -> list:
    return [c for c in cells if key(c) not in bad]


def clean_xylem(doc: dict, bad: set) -> int:
    n = 0
    for g in doc["groups"].values():
        if "cells" not in g:
            continue
        before = len(g["cells"])
        g["cells"] = drop(g["cells"], bad)
        n += before - len(g["cells"])
        g["n_cells"] = len(g["cells"])
        g["h_mech_mean_across_cells"] = float(np.mean([c["h_mech_mean"] for c in g["cells"]]))
    return n


def clean_refugium(doc: dict, bad: set) -> int:
    n = 0
    for g in doc["groups"].values():
        if "cells" not in g:
            continue
        before = len(g["cells"])
        g["cells"] = drop(g["cells"], bad)
        n += before - len(g["cells"])
        g["n_cells"] = len(g["cells"])
        g["n_robust_refugia_criterion_a"] = int(sum(1 for c in g["cells"] if c["robust_refugium_criterion_a"]))
        g["mean_viability_across_cells"] = float(np.mean([c["viability_ensemble_mean"] for c in g["cells"]]))
        g["mean_refugium_score_across_cells"] = float(np.mean([c["refugium_score"] for c in g["cells"]]))
    return n


def clean_scenarios(doc: dict, bad: set) -> int:
    n = 0
    for m in doc["members"].values():
        for g in m["groups"].values():
            before = len(g["cells"])
            g["cells"] = drop(g["cells"], bad)
            n += before - len(g["cells"])
            g["n_cells"] = len(g["cells"])
            g["mean_viability"] = float(np.mean([c["viability_mean"] for c in g["cells"]]))
            g["mean_h_mech"] = float(np.mean([c["h_mech_mean"] for c in g["cells"]]))
            g["mean_refugium_score"] = float(np.mean([c["refugium_score"] for c in g["cells"]]))
        for part in ("climate", "generic"):
            if part in m:
                m[part] = drop(m[part], bad)
    return n


TARGETS = {
    "xylem": ("xylem_mechanistic_hazard_2019_dense.yaml", clean_xylem),
    "refugium": ("refugium_viability_2019_dense.yaml", clean_refugium),
    "scenarios": ("future_projections_dense.yaml", clean_scenarios),
}


def main():
    wanted = sys.argv[sys.argv.index("--files") + 1:] if "--files" in sys.argv else list(TARGETS)
    topo = yaml.load(open(FITTED / "topohydro_grid_run_2019_dense.yaml"), Loader=LOADER)
    bad = invalid_nodes(topo)
    print(f"{len(bad)} nodes are not accepted by TOPOHYDRO's dense run")
    for name in wanted:
        fname, fn = TARGETS[name]
        path = FITTED / fname
        if not path.exists():
            print(f"  {fname}: not there yet, skipped")
            continue
        doc = yaml.load(open(path), Loader=LOADER)
        removed = fn(doc, bad)
        if removed:
            doc["dropped_nodes"] = f"{removed} cell records of nodes TOPOHYDRO rejected (non-finite soil parameters) removed by scripts/drop_invalid_nodes.py"
            path.write_text(yaml.dump(doc, Dumper=DUMPER, sort_keys=False, default_flow_style=False))
        print(f"  {fname}: {removed} cell records removed")


if __name__ == "__main__":
    main()
