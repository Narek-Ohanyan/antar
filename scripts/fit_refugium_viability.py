"""REFUGIUM's first real pass: robust refugia and risk-averse viability scores,
computed against the real per-cell XYLEM mechanistic-hazard ensemble this session
just produced, for the same 78 real grid cells.

Deliberately reduced scope, stated plainly (the same pattern as every other real
run this session), not the full Eq. 8.6 conjunction:

* **Viability = hydraulic survival only, p_height_ok fixed at 1.0.** REFUGIUM's
  `antar.viability.cohort.viability` multiplies survival-to-horizon by
  `P[H_T >= H_min]`, MERISTEM's *growth/attainable-height* model output -- a
  different MERISTEM piece from the adult-niche presence-background model this
  session already fit. The growth/height model is NOT fit against real data
  anywhere in this project yet (ROADMAP.md: needs real GEDI canopy-height and
  tree-ring/plantation growth data, neither obtained). Rather than substitute
  the adult-niche probability (a different quantity MERISTEM computes, not
  what `viability()` calls for) or fabricate a height term, `p_height_ok=1.0`
  is left neutral and documented -- this run reports real one-year hydraulic
  survival, not the full multi-decade cohort-to-canopy viability Eq. 8.6 means.
* **One year (2019), not a multi-decade horizon.** `cohort.viability` is meant
  to compound survival across years to a horizon T; only one real year's
  mechanistic hazard exists, so "survival" here is a single-year probability,
  not a compounded trajectory to 2050/2080/2100.
* **`robust_refugium`'s criterion (a) only** (ensemble viability threshold).
  Criteria (b) (`buffer_index`, topographic exposure buffering) and (c) (area
  of applicability) are both skipped, not faked: `buffer_index` needs a dense
  raster (`scipy.ndimage.uniform_filter` over a regular grid) to compute a
  meaningful regional-minus-local exposure contrast, and this session's 80
  points are a sparse irregular subsample, not a raster -- applying it here
  would silently misuse the function outside its valid domain. No
  `AreaOfApplicability` has been fit at this scale either. `robust_refugium`
  already supports (a) alone by design (`buffer_index_grid`/`inside_aoa`
  default to `None`), so this is the function's own documented partial mode,
  not a workaround.

The real 50-draw outer-loop ensemble from XYLEM's Sec. 6.3 Monte Carlo (this
session's real trait-knowledge uncertainty) is what `robust_refugium`'s
criterion (a) is evaluated against -- not a fabricated ensemble.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    compute_grid_forcing_multi_group, ROOTING_DEPTH_MM_BY_GROUP,
    GRID_ROWS, GRID_COLS, DENSE_GRID_ROWS, DENSE_GRID_COLS,
)
from fit_xylem_mechanistic_hazard import load_functional_groups, PET_FORMULATION, OUTER_DRAWS, INNER_DRAWS  # noqa: E402

from antar.hydraulics.monte_carlo import two_level_failure_probability  # noqa: E402
from antar.hydraulics.twophase import simulate_two_phase  # noqa: E402
from antar.viability.cohort import viability  # noqa: E402
from antar.viability.refugia import refugium_score, robust_refugium  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "refugium_viability_2019.yaml"
# Real 2026-10-01 densification: --dense runs against the Armenia-only DENSE_GRID, writing to a
# separate file so the original validation-grid result (what AEGIS keys off of) stays intact.
OUT_PATH_DENSE = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "refugium_viability_2019_dense.yaml"
V_STAR = 0.6   # viability threshold (robust_refugium's own default)
RHO = 0.8      # required ensemble fraction (robust_refugium's own default)
LAM = 0.5      # refugium_score's own default risk-aversion weight


def main():
    dense = "--dense" in sys.argv
    grid_rows, grid_cols = (DENSE_GRID_ROWS, DENSE_GRID_COLS) if dense else (GRID_ROWS, GRID_COLS)
    out_path = OUT_PATH_DENSE if dense else OUT_PATH

    groups = load_functional_groups()

    # Real per-group forcing, one efficient call (shared streamed extraction, see
    # run_topohydro_grid.py's module docstring and fit_xylem_mechanistic_hazard.py's same pattern).
    print(f"=== Building real {'DENSE (Armenia-only, stride 7)' if dense else '2019'} gridded TOPOHYDRO "
          f"forcing per group (rooting_depth_mm={ROOTING_DEPTH_MM_BY_GROUP}) ===", flush=True)
    rooting_depth_for_fit_groups = {g: ROOTING_DEPTH_MM_BY_GROUP[g] for g in groups
                                     if groups[g]["status"] == "ok"}
    lats, lons, elevation, cells_by_group = compute_grid_forcing_multi_group(
        rooting_depth_by_group=rooting_depth_for_fit_groups, grid_rows=grid_rows, grid_cols=grid_cols)
    n = len(lats)

    results = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "year": 2019,
        "grid": "dense_armenia_stride7" if dense else "validation_80pt_stride40",
        "v_star": V_STAR, "rho": RHO, "lam": LAM,
        "scope_note": ("Viability = one-year hydraulic survival only (p_height_ok=1.0 -- "
                        "MERISTEM's growth/attainable-height model is not fit against real "
                        "data anywhere in this project yet). robust_refugium's criterion (a) "
                        "only -- buffer index and area-of-applicability both need a dense "
                        "raster / fitted domain this sparse 78-point sample does not provide."),
        "groups": {},
    }

    for name, g in groups.items():
        if g["status"] != "ok":
            results["groups"][name] = {"status": g["status"]}
            continue

        cells = cells_by_group[name]
        valid_idx = [i for i in range(n) if cells[i] is not None]
        print(f"=== {name}: {len(valid_idx)}/{n} cells have real forcing ===", flush=True)

        print(f"=== {name}: real viability ensemble per cell ===", flush=True)
        v_ensemble = np.full((OUTER_DRAWS, len(valid_idx)), np.nan)
        hydro = []   # group-specific (rooting-depth-dependent) water stress, kept so scenario change can be shown against it
        for col, i in enumerate(valid_idx):
            cell = cells[i]
            psi_soil = cell.psi_soil_mpa[PET_FORMULATION]
            hydro.append({"cwd_mm": float(cell.cwd_mm[PET_FORMULATION]), "wsi": float(cell.wsi[PET_FORMULATION]),
                          "psi_min_mpa": float(np.min(psi_soil))})

            def simulate(traits, psi_soil=psi_soil, cell=cell):
                sim = simulate_two_phase(psi_soil, cell.t_max_c, cell.vpd_24h_kpa, cell.pressure_kpa, traits, 1.0)
                return sim["hfi_max"]

            h = two_level_failure_probability(
                simulate, g["base"], g["hyper_sd"], g["individual_sd"],
                outer_draws=OUTER_DRAWS, inner_draws=INNER_DRAWS, seed=i,
            )
            v_ensemble[:, col] = np.array([viability([hh], p_height_ok=1.0) for hh in h])
            if (col + 1) % 20 == 0:
                print(f"  [{col + 1}/{len(valid_idx)}] done", flush=True)

        score = refugium_score(v_ensemble, lam=LAM)
        robust_mask = robust_refugium(v_ensemble, v_star=V_STAR, rho=RHO)

        cell_rows = []
        for col, i in enumerate(valid_idx):
            cell_rows.append({
                "lat": float(lats[i]), "lon": float(lons[i]), "elevation_m": float(elevation[i]),
                "viability_ensemble_mean": float(v_ensemble[:, col].mean()),
                "viability_ensemble_min": float(v_ensemble[:, col].min()),
                "viability_ensemble_max": float(v_ensemble[:, col].max()),
                "refugium_score": float(score[col]),
                "robust_refugium_criterion_a": bool(robust_mask[col]),
                "p_viable": float(np.mean(v_ensemble[:, col] >= V_STAR)),
                **hydro[col],
            })
        n_robust = int(robust_mask.sum())
        results["groups"][name] = {
            "status": "ok", "example_taxa": g["example_taxa"],
            "n_cells": len(valid_idx), "n_robust_refugia_criterion_a": n_robust,
            "mean_viability_across_cells": float(v_ensemble.mean()),
            "mean_refugium_score_across_cells": float(np.mean(score)),
            "cells": cell_rows,
        }
        print(f"  {name}: mean viability={v_ensemble.mean():.4f}, "
              f"{n_robust}/{len(valid_idx)} cells meet criterion (a)", flush=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.dump(results, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {out_path} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
