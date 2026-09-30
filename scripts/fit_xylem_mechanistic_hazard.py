"""XYLEM's real mechanistic hydraulic-failure hazard, run for the first time against
real per-cell climate forcing rather than synthetic test data.

Reuses `run_topohydro_grid.compute_grid_forcing()` (this session's real gridded
TOPOHYDRO run -- SRTM/SoilGrids/ERA5-Land/CHELSA-daily, 78 valid points, 2019) to get
real `antar.climate.forcing.CellTopoclimate` objects, then runs
`antar.hydraulics.pipeline.mechanistic_hazard_for_cell`'s full Sec. 6.3 two-level
Monte Carlo (50 outer hyperparameter draws x 200 inner individual draws, the spec's
own numbers -- no reduction needed, this runs in ~5 minutes) for each of Armenia's
four defined functional groups against every real cell.

This produces the first real hazard surface anywhere in ANTAR: for each (cell,
functional group), the fraction of the outer-loop hyperparameter draws at which
mean individual-level failure probability reaches >=1 event -- i.e. `h_mech`'s own
spread across the outer loop IS the framework's trait-knowledge uncertainty
(Sec. 6.3), reported per cell rather than averaged away.

One real gap, stated plainly rather than silently worked around: Sec. 6.3 specifies
two variance components -- an outer-loop *hyperparameter* uncertainty (group-level
posterior spread) and an inner-loop *individual* uncertainty (tree-to-tree variation
within a population). `configs/fitted/xylem_trait_hyper_sd.yaml` only has the former
(the real between-record P50/slope spread XFT's pooled, multi-species, multi-study
records give us) -- XFT's compilation structure does not let the two be separated.
No individual-level variance estimate exists anywhere in the pulled data. Using
individual_sd = {} would make the inner loop degenerate (every individual identical,
collapsing h_mech to exactly 0 or 1 per outer draw instead of a continuous fraction),
so INDIVIDUAL_SD_FRACTION below scales the real hyper_sd down as an explicit,
documented placeholder standing in for the missing finer-grained estimate -- not
presented as measured.

juniper_arid_conifer has only n=1 XFT record (Juniperus thurifera) and therefore no
defined spread at all (`xylem_trait_hyper_sd.yaml`: p50_sd/slope_sd both null) --
skipped here rather than run with a fabricated hyper_sd, the same precedent
MERISTEM's adult-niche fit set for Pinus kochiana (n=8, skipped as too sparse,
not forced).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import compute_grid_forcing  # noqa: E402

from antar.hydraulics.pipeline import mechanistic_hazard_for_cell  # noqa: E402
from antar.hydraulics.twophase import Traits  # noqa: E402

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"
OUT_PATH = CONFIG_DIR / "fitted" / "xylem_mechanistic_hazard_2019.yaml"

PET_FORMULATION = "pm_fao56"
OUTER_DRAWS = 50   # Sec. 6.3's own stated outer-loop draw count
INNER_DRAWS = 200  # Sec. 6.3's own stated inner-loop draw count
INDIVIDUAL_SD_FRACTION = 0.5  # documented placeholder -- see module docstring


def load_functional_groups():
    df = pd.read_csv(CONFIG_DIR / "species_traits.csv", comment="#")
    hyper = yaml.safe_load(open(CONFIG_DIR / "fitted" / "xylem_trait_hyper_sd.yaml"))["groups"]
    groups = {}
    for _, row in df.iterrows():
        group = row["group"]
        base = Traits(
            p50=float(row["p50_mpa"]), slope=float(row["slope_pct_per_mpa"]),
            psi_close=float(row["psi_close_mpa"]), capacitance=float(row["capacitance_mmol_m2_mpa"]),
            g25=float(row["gmin25_mmol_m2_s"]), tp=float(row["tp_c"]), lethal_plc=float(row["lethal_plc"]),
        )
        h = hyper.get(group, {})
        if h.get("p50_sd") is None or h.get("slope_sd") is None:
            groups[group] = {"base": base, "status": "skipped_no_real_trait_variance_estimate",
                              "example_taxa": row["example_taxa"]}
            continue
        hyper_sd = {"p50": float(h["p50_sd"]), "slope": float(h["slope_sd"])}
        individual_sd = {k: v * INDIVIDUAL_SD_FRACTION for k, v in hyper_sd.items()}
        groups[group] = {
            "base": base, "status": "ok", "hyper_sd": hyper_sd, "individual_sd": individual_sd,
            "example_taxa": row["example_taxa"],
            "hyper_sd_source_n": {"p50": h["p50_n"], "slope": h["slope_n"]},
        }
    return groups


def main():
    groups = load_functional_groups()
    for name, g in groups.items():
        print(f"{name}: {g['status']}" + (f" (hyper_sd={g.get('hyper_sd')})" if g["status"] == "ok" else ""),
              flush=True)

    print("=== Building real 2019 gridded TOPOHYDRO forcing (reused from run_topohydro_grid) ===", flush=True)
    lats, lons, elevation, cells = compute_grid_forcing()
    n = len(lats)
    n_valid_cells = sum(1 for c in cells if c is not None)
    print(f"=== {n_valid_cells}/{n} cells have real forcing ===", flush=True)

    results = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "year": 2019,
        "pet_formulation": PET_FORMULATION,
        "outer_draws": OUTER_DRAWS,
        "inner_draws": INNER_DRAWS,
        "individual_sd_fraction_of_hyper_sd_placeholder": INDIVIDUAL_SD_FRACTION,
        "refill_fraction": 1.0,
        "groups": {},
    }

    for name, g in groups.items():
        print(f"=== {name} ===", flush=True)
        if g["status"] != "ok":
            results["groups"][name] = {"status": g["status"], "example_taxa": g["example_taxa"]}
            continue
        cell_results = []
        for i in range(n):
            if cells[i] is None:
                continue
            h = mechanistic_hazard_for_cell(
                cells[i], g["base"], g["hyper_sd"], g["individual_sd"],
                pet_formulation=PET_FORMULATION, outer_draws=OUTER_DRAWS, inner_draws=INNER_DRAWS, seed=i,
            )
            cell_results.append({
                "lat": float(lats[i]), "lon": float(lons[i]), "elevation_m": float(elevation[i]),
                "h_mech_mean": float(h.mean()), "h_mech_min": float(h.min()), "h_mech_max": float(h.max()),
                "h_mech_outer_spread_sd": float(h.std()),
            })
            if (len(cell_results)) % 20 == 0:
                print(f"  [{len(cell_results)}/{n_valid_cells}] done", flush=True)
        results["groups"][name] = {
            "status": "ok", "example_taxa": g["example_taxa"],
            "hyper_sd": g["hyper_sd"], "individual_sd": g["individual_sd"],
            "hyper_sd_source_n": g["hyper_sd_source_n"],
            "n_cells": len(cell_results),
            "h_mech_mean_across_cells": float(np.mean([c["h_mech_mean"] for c in cell_results])),
            "cells": cell_results,
        }
        print(f"  mean h_mech across cells: {results['groups'][name]['h_mech_mean_across_cells']:.4f}", flush=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.dump(results, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
