"""Split the spread of the 45-member scenario ensemble into climate model, emissions path and horizon (antar.uncertainty.partition).

Reads the dense scenario and treeline results already computed (no new simulation) and writes configs/fitted/uncertainty_partition_dense.yaml: for each of five metrics, the share of the
variance across members that belongs to each factor, pooled over horizons and at each horizon, as median and interquartile range over grid nodes.
"""
import datetime
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from antar.uncertainty.partition import partition_metric  # noqa: E402

FITTED = ROOT / "configs" / "fitted"
OUT = FITTED / "uncertainty_partition_dense.yaml"
OAK = "ring_porous_oak"                      # the group whose viability varies between nodes (broadleaf and pine sit near 1 almost everywhere)


def load(name):
    return yaml.load(open(FITTED / name), Loader=yaml.CSafeLoader)


def key(lat, lon):
    return (round(float(lat), 6), round(float(lon), 6))


def member_id(m):
    return (m["gcm"], m["scenario"], int(m["horizon"]))


def main():
    proj = load("future_projections_dense.yaml")
    tl = load("treeline_change_dense.yaml")
    topo = load("topohydro_grid_run_2019_dense.yaml")
    ref = load("refugium_viability_2019_dense.yaml")

    base_climate = {key(p["lat"], p["lon"]): p for p in topo["points"] if p.get("status") == "ok"}
    base_oak = {key(c["lat"], c["lon"]): c for c in ref["groups"][OAK]["cells"]}
    nodes = [k for k in base_oak if k in base_climate]                    # nodes with a valid water balance in every run
    members = proj["members"]
    print(f"=== {len(members)} scenario members, {len(nodes)} nodes with a valid water balance, {tl['n_cells']} treeline nodes ===", flush=True)

    def per_member(extract):
        out = {}
        for m in members.values():
            out[member_id(m)] = np.array(extract(m), dtype=float)
        return out

    def by_node(rows, field, baseline):
        d = {key(r["lat"], r["lon"]): r[field] for r in rows}
        return lambda: np.array([d.get(k, np.nan) - baseline(k) for k in nodes])

    metrics = {}
    # 1. treeline shift (all treeline nodes; the file stores it per member as a list over its cells)
    metrics["treeline_shift_m"] = {"label": "Treeline shift", "unit": "m", "values": {member_id(m): np.array(m["shift_m"], dtype=float) for m in tl["members"].values()}}
    # 2-3. warming and precipitation change against the 2019 climate of the same node
    metrics["warming_c"] = {"label": "Warming (annual mean temperature)", "unit": "°C",
                             "values": {member_id(m): np.array([{key(c["lat"], c["lon"]): c["t_mean_c"] for c in m["climate"]}.get(k, np.nan) - base_climate[k]["t_mean_c_annual_mean"] for k in nodes]) for m in members.values()}}
    metrics["precip_change_pct"] = {"label": "Precipitation change", "unit": "%",
                                     "values": {member_id(m): np.array([100.0 * ({key(c["lat"], c["lon"]): c["precip_mm"] for c in m["climate"]}.get(k, np.nan) / base_climate[k]["p_mm_annual_sum"] - 1.0) for k in nodes]) for m in members.values()}}
    # 4-5. oak water deficit and viability change against the 2019 run
    metrics["cwd_change_oak_mm"] = {"label": "Climatic water deficit change (oak rooting depth)", "unit": "mm",
                                     "values": {member_id(m): np.array([{key(c["lat"], c["lon"]): c["cwd_mm"] for c in m["groups"][OAK]["cells"]}.get(k, np.nan) - base_oak[k]["cwd_mm"] for k in nodes]) for m in members.values()}}
    metrics["viability_change_oak"] = {"label": "Oak viability change", "unit": "probability",
                                        "values": {member_id(m): np.array([{key(c["lat"], c["lon"]): c["viability_mean"] for c in m["groups"][OAK]["cells"]}.get(k, np.nan) - base_oak[k]["viability_ensemble_mean"] for k in nodes]) for m in members.values()}}

    out = {"run_date": datetime.date.today().isoformat(), "grid": proj.get("grid") or "dense_armenia_stride7",
           "design": "full factorial: 5 climate models x 3 emissions paths x 3 horizons = 45 members, one run each",
           "method": "main-effect variance fractions (Hawkins & Sutton 2009; Lehner et al. 2020) per node, summarised by median and interquartile range over nodes",
           "internal_variability_note": "One realisation per model and path, so internal variability is not separated from the model effect: it is inside the climate-model share and the remainder.",
           "metrics": {}}
    for mid, spec in metrics.items():
        r = partition_metric(spec["values"])
        r["label"], r["unit"] = spec["label"], spec["unit"]
        out["metrics"][mid] = r
        p = r["pooled"]
        print(f"  {mid}: used {r['n_nodes_used']}/{r['n_nodes']} nodes | pooled gcm {p['gcm']['median']:.2f} ssp {p['ssp']['median']:.2f} horizon {p['horizon']['median']:.2f} "
              f"remainder {p['interaction']['median']:.2f} | at 2100: gcm {r['by_horizon'][2100]['gcm']['median']:.2f} ssp {r['by_horizon'][2100]['ssp']['median']:.2f}", flush=True)
    OUT.write_text(yaml.dump(out, sort_keys=False, default_flow_style=False, allow_unicode=True))
    print(f"=== wrote {OUT.name} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
