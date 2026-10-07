"""How much do the two TOPOHYDRO placeholders matter for late-frost days, now that they are counted from budburst to the warmest day?

Sweeps the budburst threshold (200 degC-days above 5 degC) and the calm-clear-night fraction (0.3, which scales cold-air pooling of the minimum
temperature) over every dense-grid node TOPOHYDRO accepted, for the 2019 reference and for GFDL-ESM4, SSP5-8.5, 2100, and records the mean count and
the number of nodes with at least one late-frost day. The count depends on the daily temperatures alone, so no water balance or hydraulic Monte
Carlo is involved (``antar.climate.forcing.late_frost_days_from_reference``, pinned to the full pipeline by a test).

Updates configs/fitted/sensitivity_placeholders.yaml: the ``late_frost_days`` entries, which were measured under the earlier whole-year definition,
are replaced by the ``late_frost_budburst_to_warmest_day`` section; the VPD and hazard sweeps, which late frost does not feed, are left as they were.

    python3 scripts/sensitivity_late_frost.py
"""
import datetime
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import recompute_late_frost as rl  # noqa: E402
from antar.climate.forcing import late_frost_days_from_reference  # noqa: E402

OUT = rl.FITTED / "sensitivity_placeholders.yaml"
BASE = {"gdd_budburst": 200.0, "calm_frac": 0.3}
SWEEPS = {
    "gdd_budburst": [25.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0],
    "calm_frac": [0.0, 0.15, 0.3, 0.45, 0.6],
}
MEMBER = ("gfdl-esm4", "ssp585", 2100)


def main():
    spec = rl.GRIDS["dense"]
    static = rl.extract_static_grid_inputs(spec["rows"], spec["cols"])
    t_mean, t_min, month, gamma = rl.reference_2019(static)
    base = yaml.load(open(rl.FITTED / spec["baseline"]), Loader=rl.LOADER)
    cells = [i for i, p in enumerate(base["points"]) if p.get("status") == "ok"]
    gcm, ssp, horizon = MEMBER
    deltas = rl.compute_deltas(static["lats"], static["lons"], gcm, ssp)

    def mean_and_share(params, member):
        counts = []
        for i in cells:
            n = rl.Node(static, t_mean, t_min, month, gamma, i)
            tm, tn = n.t_mean_ref, n.t_min_ref
            if member:
                tm, tn = tm + deltas[(horizon, "tas")][month - 1, i], tn + deltas[(horizon, "tasmin")][month - 1, i]
            counts.append(late_frost_days_from_reference(
                t_mean_ref_c=tm, t_min_ref_c=tn, month=month, z_cell_m=n.z, z_ref_m=n.z_ref, gamma_k_per_m=gamma, concavity_index=n.concavity,
                calm_clear_night_frac=params["calm_frac"], gdd_budburst=params["gdd_budburst"]))
        c = np.asarray(counts)
        return {"mean_days": float(c.mean()), "nodes_with_a_frost_day": int((c > 0).sum()), "max_days": int(c.max())}

    section = {"definition": "frost days (T_min < -2 degC) from budburst to the warmest day of the year", "n_cells": len(cells), "base": BASE, "periods": {}}
    for period, member in (("2019", False), (f"{gcm}_{ssp}_{horizon}", True)):
        section["periods"][period] = {}
        for name, values in SWEEPS.items():
            section["periods"][period][name] = []
            for v in values:
                r = mean_and_share(dict(BASE, **{name: v}), member)
                section["periods"][period][name].append({"value": v, **r})
                print(f"  {period} {name}={v}: mean {r['mean_days']:.2f} d, {r['nodes_with_a_frost_day']}/{len(cells)} nodes, max {r['max_days']}", flush=True)

    data = yaml.safe_load(open(OUT))
    for period in data["periods"].values():
        for rows in period.values():
            for row in rows:
                row.pop("late_frost_days", None)
    data["late_frost_budburst_to_warmest_day"] = {"run_date": datetime.date.today().isoformat(), **section}
    OUT.write_text(yaml.safe_dump(data, sort_keys=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
