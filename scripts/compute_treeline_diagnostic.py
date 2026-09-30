"""The real potential-treeline-elevation diagnostic, run for the first time against real
gridded climate, explicitly required by the user (2026-09-30): "make sure that we model the
Treeline change."

`antar.niche.growth.potential_treeline_elevation` inverts the Sec. 5.1 lapse relation
T(z) = T_ref + Gamma*(z - z_ref) for the elevation at which growing-season mean temperature
equals the Koerner-Paulsen thermal treeline threshold (6.45 degC) -- a statement of each cell's
real climatic ceiling, not a forecast of where forest actually stands (the function's own
docstring: "realised treelines lag their climatic potential by decades").

Reuses the real 2019 gridded TOPOHYDRO run (78 cells) this session already produced. One real
simplification, noted because it isn't obvious: the function takes a single reference elevation/
temperature pair, but the lapse relation is linear, so inverting from the cell's own real
downscaled (z_cell_m, growing_season_mean_t_c) gives the exact same z_tl as inverting from the
true upstream reference (z_ref_m, gst_ref_c) -- algebraically identical, not an approximation,
so no extra data pull was needed to run this.

`gamma_k_per_m`: the real fitted monthly lapse rate (`topohydro_lapse_rate.yaml`) is monthly,
but this diagnostic wants one representative value. Using April-September's real fitted values
(the temperate growing-season months) rather than a flat 12-month average -- Armenia's real
fitted lapse rate is steeper in summer than winter (-7.6 vs -4.9 K/km, IMPLEMENTATION_LOG.md
2026-09-30), so a growing-season-specific average is the more defensible choice for a diagnostic
that is itself about growing-season temperature.
"""
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import compute_grid_forcing  # noqa: E402

from antar.niche.growth import potential_treeline_elevation  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "treeline_diagnostic_2019.yaml"
GROWING_SEASON_MONTHS = ["Apr", "May", "Jun", "Jul", "Aug", "Sep"]


def main():
    lapse = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "fitted" / "topohydro_lapse_rate.yaml"))
    gamma_growing_season = float(np.mean([lapse["temperature"]["gamma_k_per_m_by_month"][m] for m in GROWING_SEASON_MONTHS]))
    print(f"=== Growing-season (Apr-Sep) mean lapse rate: {gamma_growing_season*1000:.2f} K/km ===", flush=True)

    print("=== Building real 2019 gridded TOPOHYDRO forcing (reused) ===", flush=True)
    lats, lons, elevation, cells = compute_grid_forcing()
    n = len(lats)

    rows = []
    for i in range(n):
        if cells[i] is None:
            continue
        c = cells[i]
        z_tl = potential_treeline_elevation(c.growing_season_mean_t_c, elevation[i], gamma_growing_season)
        rows.append({
            "lat": float(lats[i]), "lon": float(lons[i]),
            "elevation_m": float(elevation[i]),
            "growing_season_mean_t_c": float(c.growing_season_mean_t_c),
            "potential_treeline_elevation_m": float(z_tl),
            "above_potential_treeline": bool(elevation[i] > z_tl),
        })

    n_above = sum(r["above_potential_treeline"] for r in rows)
    z_tls = [r["potential_treeline_elevation_m"] for r in rows]
    print(f"=== {len(rows)} real cells; potential treeline elevation range "
          f"{min(z_tls):.0f}-{max(z_tls):.0f} m; {n_above}/{len(rows)} cells sit above their "
          f"own real climatic ceiling ===", flush=True)

    result = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "year": 2019,
        "thermal_threshold_c": 6.45,
        "threshold_source": "Koerner & Paulsen (2004)-style global treeline GST threshold, "
                             "the value already hardcoded in antar.niche.growth as the default",
        "gamma_growing_season_k_per_km": gamma_growing_season * 1000,
        "gamma_note": "mean of the real fitted Apr-Sep monthly lapse rates, not a flat annual average",
        "scope_note": "Diagnostic climatic ceiling only, one real year (2019), not yet run across "
                       "the ISIMIP3b scenario ensemble/horizons -- that is the real next step for "
                       "showing treeline CHANGE rather than today's treeline position.",
        "n_cells": len(rows),
        "n_cells_above_potential_treeline": n_above,
        "cells": rows,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(yaml.dump(result, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {OUT_PATH} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
