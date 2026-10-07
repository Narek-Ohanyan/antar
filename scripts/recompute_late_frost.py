"""Recomputes late-frost days under the current definition (frosts from budburst to the warmest day of the year) in the stored results.

Late frost used to count every day below -2 degC from budburst to the end of the year, which includes the autumn frosts after growth has stopped. It
now counts only the frosts between budburst and the warmest day (``antar.climate.indices.warmest_day_index``). Only this one count changes, and it
depends on the daily temperatures alone, so the stored results are updated without redoing the water balance or the hydraulic Monte Carlo.

Before anything is written both definitions are recomputed for every node, and each stored count must equal one of them; if a single node matches
neither the script stops, because that would mean the shortcut is not the pipeline that produced the file. The new count is never larger than the old one.

    python3 scripts/recompute_late_frost.py                 # check only: reproduce the stored counts, report what would change
    python3 scripts/recompute_late_frost.py --apply         # write the new counts (idempotent)
    python3 scripts/recompute_late_frost.py --apply --only dense
"""
import datetime
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER, DENSE_GRID_COLS, DENSE_GRID_ROWS, GDD_BUDBURST_PLACEHOLDER, GRID_COLS, GRID_ROWS,
    _load_chelsa_arrays, extract_static_grid_inputs,
)
from run_future_projections import compute_deltas  # noqa: E402

from antar.climate.forcing import late_frost_days_from_reference  # noqa: E402

FITTED = Path(__file__).resolve().parent.parent / "configs" / "fitted"
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
DUMPER = getattr(yaml, "CSafeDumper", yaml.SafeDumper)
DEFINITION = "budburst_to_warmest_day"
YEAR = 2019
GRIDS = {
    "validation": {"rows": GRID_ROWS, "cols": GRID_COLS, "baseline": "topohydro_grid_run_2019.yaml", "members": "future_projections.yaml"},
    "dense": {"rows": DENSE_GRID_ROWS, "cols": DENSE_GRID_COLS, "baseline": "topohydro_grid_run_2019_dense.yaml", "members": "future_projections_dense.yaml"},
}


def key(lat, lon):
    return (round(float(lat), 5), round(float(lon), 5))


class Node:
    """One grid node's reference series and static inputs, with the count under either definition."""

    def __init__(self, static, t_mean_ref, t_min_ref, month, gamma, i):
        self.t_mean_ref, self.t_min_ref, self.month, self.gamma = t_mean_ref[:, i], t_min_ref[:, i], month, gamma
        self.z, self.z_ref, self.concavity = float(static["elevation"][i]), float(static["z_ref_m"][i]), float(static["concavity"][i])

    def count(self, spring_only, d_tas=None, d_tasmin=None):
        t_mean, t_min = self.t_mean_ref, self.t_min_ref
        if d_tas is not None:
            t_mean = t_mean + d_tas[self.month - 1]
            t_min = t_min + d_tasmin[self.month - 1]
        return late_frost_days_from_reference(
            t_mean_ref_c=t_mean, t_min_ref_c=t_min, month=self.month, z_cell_m=self.z, z_ref_m=self.z_ref, gamma_k_per_m=self.gamma,
            concavity_index=self.concavity, calm_clear_night_frac=CALM_CLEAR_NIGHT_FRAC_PLACEHOLDER, gdd_budburst=GDD_BUDBURST_PLACEHOLDER,
            spring_only=spring_only)


def reference_2019(static):
    arrs = _load_chelsa_arrays()
    day0 = datetime.date(1979, 1, 1)
    i0 = (datetime.date(YEAR, 1, 1) - day0).days
    i1 = (datetime.date(YEAR, 12, 31) - day0).days + 1
    month = np.array([(day0 + datetime.timedelta(days=d)).month for d in range(i0, i1)])
    r, c = static["chelsa_row"], static["chelsa_col"]
    lapse = yaml.safe_load(open(FITTED / "topohydro_lapse_rate.yaml"))["temperature"]["gamma_k_per_m_by_month"]
    gamma = np.array([lapse[m] for m in ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]])
    return arrs["tas"][i0:i1, r, c] - 273.15, arrs["tasmin"][i0:i1, r, c] - 273.15, month, gamma


def summarise(name, old_stored, old_re, new):
    """A stored count must equal one of the two recomputed counts (a run restarted after the change may hold some of each)."""
    old_stored, old_re, new = map(np.asarray, (old_stored, old_re, new))
    bad = int(np.sum((old_stored != old_re) & (old_stored != new)))
    print(f"  {name}: {len(new)} nodes; the stored count matches the old or the new definition at {len(new) - bad}/{len(new)}; "
          f"mean {old_stored.mean():.1f} -> {new.mean():.1f} days; new > old at {int(np.sum(new > old_re))} nodes", flush=True)
    return bad


def do_baseline(grid, spec, node_for, apply):
    path = FITTED / spec["baseline"]
    data = yaml.load(open(path), Loader=LOADER)
    stored, old_re, new, rows = [], [], [], []
    for i, p in enumerate(data["points"]):
        if p.get("status") != "ok":
            continue
        n = node_for(i)
        stored.append(p["late_frost_days"]); old_re.append(n.count(False)); new.append(n.count(True)); rows.append(p)
    bad = summarise(f"{grid} 2019 grid", stored, old_re, new)
    if bad:
        raise SystemExit(f"{grid}: {bad} stored counts match neither definition; refusing to write")
    if apply and data.get("late_frost_definition") != DEFINITION:
        for p, v in zip(rows, new):
            p["late_frost_days"] = int(v)
        data = {**{k: v for k, v in data.items() if k != "points"}, "late_frost_definition": DEFINITION, "points": data["points"]}
        path.write_text(yaml.dump(data, Dumper=DUMPER, sort_keys=False, default_flow_style=False))
        print(f"  wrote {path.name}", flush=True)


def do_members(grid, spec, static, node_for, apply):
    path = FITTED / spec["members"]
    if not path.exists():
        print(f"  {grid} scenarios: {path.name} does not exist yet; run this again once the run has finished", flush=True)
        return
    data = yaml.load(open(path), Loader=LOADER)
    index = {key(a, b): i for i, (a, b) in enumerate(zip(static["lats"], static["lons"]))}
    pairs = sorted({(m["gcm"], m["scenario"]) for m in data["members"].values()})
    stored, old_re, new, targets = [], [], [], []
    for gcm, ssp in pairs:
        deltas = compute_deltas(static["lats"], static["lons"], gcm, ssp)
        for m in data["members"].values():
            if (m["gcm"], m["scenario"]) != (gcm, ssp):
                continue
            h = m["horizon"]
            for row in m["climate"]:
                i = index[key(row["lat"], row["lon"])]
                d = (deltas[(h, "tas")][:, i], deltas[(h, "tasmin")][:, i])
                n = node_for(i)
                stored.append(row["late_frost_days"]); old_re.append(n.count(False, *d)); new.append(n.count(True, *d)); targets.append(row)
        print(f"  {gcm} {ssp} done", flush=True)
    bad = summarise(f"{grid} scenarios", stored, old_re, new)
    if bad:
        raise SystemExit(f"{grid}: {bad} stored counts match neither definition; refusing to write")
    if apply and data.get("late_frost_definition") != DEFINITION:
        for row, v in zip(targets, new):
            row["late_frost_days"] = int(v)
        data = {**{k: v for k, v in data.items() if k != "members"}, "late_frost_definition": DEFINITION, "members": data["members"]}
        path.write_text(yaml.dump(data, Dumper=DUMPER, sort_keys=False, default_flow_style=False))
        print(f"  wrote {path.name}", flush=True)


def main():
    apply = "--apply" in sys.argv
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    for grid, spec in GRIDS.items():
        if only and grid != only:
            continue
        print(f"=== {grid} grid ===", flush=True)
        static = extract_static_grid_inputs(spec["rows"], spec["cols"])
        t_mean, t_min, month, gamma = reference_2019(static)
        cache = {}

        def node_for(i, static=static, t_mean=t_mean, t_min=t_min, month=month, gamma=gamma, cache=cache):
            if i not in cache:
                cache[i] = Node(static, t_mean, t_min, month, gamma, i)
            return cache[i]

        do_baseline(grid, spec, node_for, apply)
        do_members(grid, spec, static, node_for, apply)
    if not apply:
        print("check only; nothing written (add --apply to write)")


if __name__ == "__main__":
    main()
