"""How much do the two placeholder inputs of TOPOHYDRO matter?

* the budburst threshold (200 degC-days above 5 degC) that sets when late frost starts to count, and
* the fraction of clear, calm nights (0.3) that scales cold-air pooling of the minimum temperature.

Neither has a data source (a cited species threshold could not be verified; ERA5-Land has no cloud band). Instead of defending a number, this
script sweeps each over a wide range with the real seasonal, scenario-dependent atmosphere (everything else identical, same seeds) and reports
what moves: the late-frost day count, the summer 24 h VPD, and the XYLEM hydraulic hazard of each group. Run for the 2019 reference and for one
member of the ensemble (GFDL-ESM4, SSP5-8.5, 2100).

    python3 scripts/sensitivity_placeholders.py --ctx <pickle written by run_future_projections.py> --out configs/fitted/sensitivity_placeholders.yaml

The context pickle exists for the duration of a scenario run (data/_future_ctx_*.pkl); copy it before the run ends, or keep it with
``run_future_projections.py --keep-ctx``.
"""
import pickle
import sys
import datetime
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_future_projections import compute_deltas  # noqa: E402
from sensitivity_atmospheric_forcing import forcing  # noqa: E402
from fit_xylem_mechanistic_hazard import PET_FORMULATION, OUTER_DRAWS, INNER_DRAWS  # noqa: E402

from antar.hydraulics.monte_carlo import two_level_failure_probability  # noqa: E402
from antar.hydraulics.twophase import simulate_two_phase  # noqa: E402
from antar.io.armenia_mask import inside_armenia  # noqa: E402

BASE = {"gdd_budburst": 200.0, "calm_frac": 0.3}
SWEEPS = {
    "gdd_budburst": [100.0, 150.0, 200.0, 300.0, 400.0],
    "calm_frac": [0.0, 0.15, 0.3, 0.45, 0.6],
}


def evaluate(c, cells, deltas, horizon, params, with_hazard):
    out = {"late_frost_days": [], "vpd24_summer": [], "tmin_mean": [], "hazard": {g: [] for g in c["groups"]}}
    for i in cells:
        res = None
        for gname, g in c["groups"].items():
            cell = forcing(c, i, deltas, horizon, c["w_max_mm_by_group"][gname], "V2", **params)
            if res is None:
                res = cell
                out["late_frost_days"].append(float(cell.late_frost_days))
                out["vpd24_summer"].append(float(np.mean(cell.vpd_24h_kpa[(c["month"] >= 6) & (c["month"] <= 8)])))
                out["tmin_mean"].append(float(np.mean(cell.t_min_c)))
            if with_hazard:
                psi = cell.psi_soil_mpa[PET_FORMULATION]

                def simulate(traits, psi=psi, cell=cell):
                    return simulate_two_phase(psi, cell.t_max_c, cell.vpd_24h_kpa, cell.pressure_kpa, traits, 1.0)["hfi_max"]
                h = two_level_failure_probability(simulate, g["base"], g["hyper_sd"], g["individual_sd"], outer_draws=OUTER_DRAWS, inner_draws=INNER_DRAWS, seed=i)
                out["hazard"][gname].append(float(h.mean()))
    return {k: ({g: float(np.mean(v)) for g, v in val.items() if v} if k == "hazard" else float(np.mean(val))) for k, val in out.items()}


def main():
    ctx = sys.argv[sys.argv.index("--ctx") + 1]
    out_path = Path(sys.argv[sys.argv.index("--out") + 1])
    c = pickle.load(open(ctx, "rb"))
    cells = [i for i in range(c["n"]) if c["valid_mask"][i] and bool(inside_armenia(c["lats"][i], c["lons"][i]))]
    print(f"{len(cells)} Armenian cells", flush=True)
    deltas = compute_deltas(c["lats"], c["lons"], "gfdl-esm4", "ssp585")
    result = {"run_date": datetime.date.today().isoformat(), "n_cells": len(cells), "base": BASE,
              "note": "Mean over the Armenian cells of the grid; atmosphere = ISIMIP3a seasonal shape (+ ISIMIP3b change for the member).",
              "periods": {}}
    for period, d, hz in (("2019", None, None), ("gfdl-esm4_ssp585_2100", deltas, 2100)):
        rows = {}
        for name, values in SWEEPS.items():
            rows[name] = []
            for v in values:
                params = dict(BASE, **{name: v})
                # the budburst threshold only counts late-frost days: it cannot reach the hazard, so the Monte Carlo is skipped for it
                r = evaluate(c, cells, d, hz, params, with_hazard=(name != "gdd_budburst"))
                rows[name].append({"value": v, **r})
                print(f"  {period} {name}={v}: late frost {r['late_frost_days']:.1f} d, summer VPD24 {r['vpd24_summer']:.3f} kPa, hazard {r['hazard']}", flush=True)
        result["periods"][period] = rows
    out_path.write_text(yaml.safe_dump(result, sort_keys=False))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
