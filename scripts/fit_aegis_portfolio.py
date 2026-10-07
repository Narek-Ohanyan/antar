"""AEGIS's first real portfolio pass: CVaR-robust site-option selection against real viability,
real cost, real value, real eligibility and (once `future_projections.yaml` exists) a real
multi-scenario ensemble -- the four scope decisions the user made explicitly this session, all
applied: (1) the real $417/ha/yr national ecosystem-services value as the value weight, not a
stakeholder-elicited figure (the concept note's own stated preference, but genuinely unavailable);
(2) options j = species x intervention-method combined; (3) budget swept across the real range,
not a single figure; (4) real water-basin figures extracted, but NOT wired into a constraint here
-- see the water section of the docstring below for why.

**Real, stated gap in the options dimension**: ANTAR's biophysical model (XYLEM/REFUGIUM)
predicts viability per *functional group*, never per *intervention method* -- there is no real
model anywhere in this project for how "natural regeneration" vs. "coppicing" vs. "windbreaks"
changes a tree's survival probability. So in the real 3 (group) x 8 (method) = 24 option grid
built here, **benefit varies only by group** (real, from REFUGIUM) while **cost varies only by
method** (real, from the World Bank/government source) -- the two dimensions are independent by
construction, not because ecology doesn't matter, but because no real data or model connects
method to viability. Stated plainly, not hidden in the output.

**Real gap in cost granularity**: the source gives only three real anchors across the 8 named
methods -- a real "cheapest mix" ($186/ha, covering natural regeneration/coppicing/wildfire
prevention), a real most-expensive single method (windbreaks, $13,260/ha), and a real blended
average across all 8 ($900/ha). The remaining 4 methods (degraded-forest planting, thinning,
anti-erosion plantations, mining-site reclamation) have no individually-sourced real cost --
the blended $900/ha stands in for them, flagged per-method in the output, not silently applied
to all 8.

**Water**: real water_use-per-option figures do not exist in anything pulled (the real RBMP
documents track irrigation/industry/municipal allocation, never a forestry sector, because
Armenian afforestation is documented as mostly rainfed). `water_use`/`water_caps`/`basin` are
left unset -- not faked, not silently defaulted to zero and presented as "no constraint found,"
genuinely omitted because no real number exists to populate them with.

**Value vs. cost units, stated not glossed over**: value ($417/ha) is a real *annual* ecosystem-
service flow; cost is a real *one-time* establishment expenditure. No real discount rate was
sourced, so no NPV conversion is applied -- `benefit` in the optimisation is annual value accrued,
`cost` is the one-time budget constraint. Mechanically consistent (cost only enters as a budget
constraint, never subtracted from the maximised objective), but the two are different kinds of
quantity, not a mistake to silently reconcile.

**Scenario ensemble**: uses `configs/fitted/future_projections.yaml`'s real per-(GCM, SSP,
horizon) viability if it exists (the real, intended multi-scenario ensemble this is built for);
falls back to `refugium_viability_2019.yaml`'s single real 2019 scenario if not yet available
(a real but degenerate C=1 ensemble -- CVaR has nothing to be robust across with one scenario,
stated explicitly in the output, not hidden behind a number that looks like a real frontier).
"""
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform as warp_transform

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import get_access_token, drive_vsicurl_url, BBOX, CHELSA_GRID_SHAPE, DENSE_STRIDE  # noqa: E402
from antar.decision.optimize import robust_portfolio, efficient_frontier  # noqa: E402
from antar.decision.units import cell_area_ha, cell_land_shares  # noqa: E402
from PIL import Image  # noqa: E402
from antar.io.armenia_mask import inside_armenia  # noqa: E402
import time  # noqa: E402

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"
OUT_PATH = CONFIG_DIR / "fitted" / "aegis_portfolio.yaml"
# Real 2026-10-02 densification, explicitly requested: --dense reads future_projections_dense.yaml
# (the real 1044-point, 45-member ensemble) instead of the 78-point one, writing its own output
# file so the original result stays intact for comparison.
OUT_PATH_DENSE = CONFIG_DIR / "fitted" / "aegis_portfolio_dense.yaml"

VALUE_PER_HA_YEAR = 417.0  # real, World Bank 2023 -- national-average ecosystem-services value
# The same note's own cost-benefit analysis assumes degraded land still delivers 25% of a fully restored forest's ecosystem services, so the INCREMENTAL
# benefit of restoring it is 75% of the 417 (and it phases in: 50% of the increment in years 1-5, 100% from year 6). Using the full 417 overstated every
# absolute benefit by a third; the ramp is not applied here (all benefits are steady-state, year 6 onward). A constant factor does not change which plan is
# optimal, only the benefit totals.
INCREMENTAL_SHARE = 0.75
# Solver settings. Cells are nearly identical in benefit (about 3,200 ha, ~99% survival), so proving a tight gap only separates cells that differ by far less than the
# uncertainty of the 417 USD/ha/yr value. The budget sweep is accepted within 0.5% of the proven bound; the frontier, whose two ends differ by very little, is solved much more
# tightly so that the price of robustness can be told from solver noise. Whether each solve reached its gap is recorded.
SWEEP_REL_GAP, SWEEP_TIME_S = 5e-3, 180
FRONTIER_REL_GAP, FRONTIER_TIME_S = 1e-4, 180
SENS_REL_GAP, SENS_TIME_S = 1e-2, 120
GRID_STRIDE = 40           # CHELSA pixels between the nodes of the coarse validation grid (the dense grid's is DENSE_STRIDE)
LAND_TENURE_FILE_ID = "1s9ENmIU67yyi9SJcXoOtrmERYTE-i9zN"

# Real cost per hectare by intervention method (World Bank 2023 / Armenia's Landscape Restoration
# Strategy 2022-2032, Box 3) -- "sourced" = one of the report's own real anchor figures;
# "blended" = the report's own all-8-methods average, standing in because no individual real
# figure exists for that specific method (see module docstring).
INTERVENTIONS = [
    # cost_per_ha is the establishment cost used in the budget constraint; cost_low / cost_high are the range found in the sources (US dollars per hectare).
    {"name": "degraded_forest_planting", "cost_per_ha": 6800.0, "cost_low": 3300.0, "cost_high": 12600.0, "cost_basis": "sourced_armenian_cost_data",
     "source": ("Central: 2.6 million AMD per hectare, the Ministry of Environment's average cost to establish a hectare of forest including five years of maintenance, "
                "mowing, fencing and nursery work (Hetq, 1 December 2025; about US$6,800 at that article's own conversion). Low: EUR 2,945 per hectare, planting with "
                "survival taken into account, KfW ex-post evaluation of the IKI Caucasus natural-forest programme (2017; about US$3,300, prices of 2013-2016). "
                "High: about US$12,600 per hectare, the items of an itemised intensive planting site in Armenia's Adaptation Fund proposal (2025), excluding irrigation "
                "and soil works.")},
    {"name": "natural_regeneration", "cost_per_ha": 186.0, "cost_basis": "sourced_cheapest_mix", "source": "World Bank (2023): the cheapest real option mix, US$9.3 million for 50,000 ha; the KfW evaluation also finds natural regeneration much cheaper than planting, fencing being the main cost."},
    {"name": "coppicing_oak", "cost_per_ha": 186.0, "cost_basis": "sourced_cheapest_mix", "source": "World Bank (2023): the cheapest real option mix (shared with natural regeneration and wildfire prevention)."},
    {"name": "pine_thinning", "cost_per_ha": 900.0, "cost_basis": "blended", "source": "World Bank (2023): the average over all eight options, US$45 million for 50,000 ha; no individual figure found."},
    {"name": "wildfire_prevention", "cost_per_ha": 186.0, "cost_basis": "sourced_cheapest_mix", "source": "World Bank (2023): the cheapest real option mix (shared with natural regeneration and coppicing)."},
    {"name": "anti_erosion_plantation", "cost_per_ha": 900.0, "cost_basis": "blended", "source": "World Bank (2023): the average over all eight options; no individual figure found."},
    {"name": "windbreaks_hedgerows", "cost_per_ha": 13260.0, "cost_basis": "sourced_most_expensive", "source": "World Bank (2023): the most expensive single option, US$663 million for 50,000 ha."},
    {"name": "mining_site_reclamation", "cost_per_ha": 900.0, "cost_basis": "blended", "source": "World Bank (2023): the average over all eight options; no individual figure found."},
]

# What each method acts on. The benefit figure (417 USD/ha/yr x 75%) is the value a RESTORED DEGRADED hectare adds, so only the restoration methods are ranked, and each acts on
# the open land of a cell (not forest, woodland or water). The three methods that act on existing forest have no benefit figure for maintaining forest, and mining reclamation has
# no data on where the mine sites are: they are listed but not ranked, which also stops the plan from "restoring" land that is already forest.
ROLE = {"degraded_forest_planting": "restoration", "natural_regeneration": "restoration", "anti_erosion_plantation": "restoration", "windbreaks_hedgerows": "restoration",
        "coppicing_oak": "forest_management", "pine_thinning": "forest_management", "wildfire_prevention": "forest_management", "mining_site_reclamation": "site_specific"}
NOT_RANKED_REASON = {"forest_management": "acts on existing forest, and no figure exists for the benefit of maintaining forest (the 417 USD/ha/yr x 75% value is that of a restored degraded hectare)",
                     "site_specific": "acts on mine sites, and there is no data on where they are"}
for _i in INTERVENTIONS:
    _i["role"] = ROLE[_i["name"]]
RANKED = [i for i in INTERVENTIONS if i["role"] == "restoration"]
MIN_OPEN_SHARE = 0.10      # a cell is a candidate only if at least this share of it is open land (a cell that is nearly all forest, woodland or water has nothing to restore)

# Real budget range (World Bank 2023): $9.3M (cheapest real option mix) to $663M (most expensive
# single real option) for the 50,000 ha NDC target; $45M is the real figure tied to the
# government's own actual recommended proportional mix -- swept per the user's explicit decision.
BUDGET_LEVELS_USD = [9.3e6, 45e6, 100e6, 300e6, 663e6]
REPRESENTATIVE_BUDGET_USD = 45e6


WATER_CLASS = 1001  # Ecosystem Map "Water bodies" (e.g. Lake Sevan): inside Armenia but not plantable


def load_eligibility(lats, lons):
    """Eligible = inside Armenia AND not protected (WDPA) AND not mostly human-modified or water.

    The first condition was missing until 2026-10-05: the study grid is a rectangle that Armenia
    fills only ~37% of, and a foreign cell has no protected-area flag and no mapped human-modified
    land, so it passed every other test -- the portfolio was recommending planting in Georgia,
    Azerbaijan, Turkey and Iran. (Found while building the UI.)

    The ecosystem-map exclusion is sampled LIVE here (not looked up from a file keyed to the 78-point
    grid's exact coordinates, which would silently miss at any other grid density)."""
    inside = inside_armenia(lats, lons)

    protected = None
    for attempt in range(3):
        try:
            token = get_access_token()
            url = drive_vsicurl_url(LAND_TENURE_FILE_ID)
            with rasterio.Env(GDAL_HTTP_HEADERS=f"Authorization: Bearer {token}", GDAL_DISABLE_READDIR_ON_OPEN="YES", GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10):
                with rasterio.open(url) as src:
                    xs, ys = warp_transform("EPSG:4326", src.crs, lons.tolist(), lats.tolist())
                    protected = np.array(list(src.sample(zip(xs, ys))))[:, 0]
            break
        except rasterio.errors.RasterioIOError as e:
            print(f"  WDPA read failed (attempt {attempt + 1}/3): {e}", flush=True)
            time.sleep(5 * (attempt + 1))
    if protected is None:
        raise RuntimeError("WDPA protected-area raster unreadable after retries -- refusing to assume 'not protected'")
    not_protected = protected < 0.5

    # Deliberately does NOT exclude already-forested cells -- 3 of the 8 real intervention methods
    # (coppicing_oak, pine_thinning, wildfire_prevention) target existing forest, so that would wrongly
    # zero out the options that most need it. See integrate_ecosystem_map.py's docstring.
    from integrate_ecosystem_map import sample_class_fractions, HUMAN_MODIFIED_CLASSES
    class_fractions = sample_class_fractions(lats, lons)
    plantable = np.ones(len(lats), dtype=bool)
    for i, cf in enumerate(class_fractions):
        if cf is not None:
            plantable[i] = sum(cf.get(c, 0.0) for c in HUMAN_MODIFIED_CLASSES + [WATER_CLASS]) <= 0.5
    print(f"  eligibility: {int(inside.sum())} inside Armenia | {int((inside & not_protected).sum())} also not protected "
          f"| {int((inside & not_protected & plantable).sum())} also plantable (not human-modified / water)", flush=True)
    return inside & not_protected & plantable


def load_refugium_only():
    """Fallback: real 2019 single-scenario viability (C=1, degenerate for CVaR -- stated in output)."""
    d = yaml.safe_load(open(CONFIG_DIR / "fitted" / "refugium_viability_2019.yaml"))
    groups = {k: v for k, v in d["groups"].items() if v.get("status") == "ok"}
    group_names = list(groups.keys())
    lats = np.array([c["lat"] for c in groups[group_names[0]]["cells"]])
    lons = np.array([c["lon"] for c in groups[group_names[0]]["cells"]])
    n = len(lats)
    # viability[group][cell] -> (n, 1) scenario axis
    viability = {g: np.array([c["viability_ensemble_mean"] for c in groups[g]["cells"]]).reshape(n, 1)
                 for g in group_names}
    return lats, lons, group_names, viability, ["2019_single_scenario"]


def load_future_projections(dense=False):
    """Real multi-scenario ensemble, if scripts/run_future_projections.py has finished."""
    path = CONFIG_DIR / "fitted" / ("future_projections_dense.yaml" if dense else "future_projections.yaml")
    if not path.exists():
        return None
    d = yaml.safe_load(open(path))
    members = d["members"]
    scenario_names = sorted(members.keys())
    first = members[scenario_names[0]]
    group_names = list(first["groups"].keys())
    cells0 = first["groups"][group_names[0]]["cells"]
    lats = np.array([c["lat"] for c in cells0])
    lons = np.array([c["lon"] for c in cells0])
    n = len(lats)
    viability = {}
    for g in group_names:
        arr = np.zeros((n, len(scenario_names)))
        for ci, sname in enumerate(scenario_names):
            cells = members[sname]["groups"][g]["cells"]
            arr[:, ci] = [c["viability_mean"] for c in cells]
        viability[g] = arr
    return lats, lons, group_names, viability, scenario_names



def allocation_summary(plan, area, cost, option_labels):
    """What a plan treats: units, hectares and cost in all, and the share by intervention method and by species group."""
    x = np.asarray(plan["x"]) > 0.5
    area = area[:, None] * np.ones((1, x.shape[1])) if np.ndim(area) == 1 else area
    ha, usd = x * area, x * area * cost
    out = {"n_units_planted": int(x.sum()), "area_ha": float(ha.sum()), "cost_usd": float(usd.sum()), "by_intervention": {}, "by_group": {}}
    groups = list(dict.fromkeys(o["group"] for o in option_labels))
    methods = list(dict.fromkeys(o["intervention"] for o in option_labels))
    # [unit index, index into units.group_ids, index into intervention_ids] for every treated cell; the map on the Decision page draws these
    out["selected"] = [[int(u), groups.index(option_labels[j]["group"]), methods.index(option_labels[j]["intervention"])] for u, j in zip(*np.nonzero(x))]
    for j, o in enumerate(option_labels):
        if x[:, j].any():
            for key, name in (("by_intervention", o["intervention"]), ("by_group", o["group"])):
                d = out[key].setdefault(name, {"units": 0, "area_ha": 0.0})
                d["units"] += int(x[:, j].sum()); d["area_ha"] += float(ha[:, j].sum())
    return out


def main():
    dense = "--dense" in sys.argv
    out_path = OUT_PATH_DENSE if dense else OUT_PATH

    future = load_future_projections(dense=dense)
    if future is not None:
        lats, lons, group_names, viability, scenario_names = future
        print(f"=== Real {'DENSE ' if dense else ''}multi-scenario ensemble: "
              f"{len(scenario_names)} real members ===", flush=True)
    elif dense:
        print("=== future_projections_dense.yaml not ready yet -- cannot run AEGIS --dense "
              "without it (no dense single-scenario fallback defined) ===", flush=True)
        return
    else:
        lats, lons, group_names, viability, scenario_names = load_refugium_only()
        print("=== future_projections.yaml not ready yet -- real single-scenario (2019) fallback, "
              "C=1, CVaR degenerate -- stated in output, not hidden ===", flush=True)

    # Foreign cells are not candidate units at all: drop them rather than merely mark them ineligible,
    # so the reported unit count and every benefit total refer to Armenia only.
    in_arm = inside_armenia(lats, lons)
    n_total = len(lats)
    lats, lons = lats[in_arm], lons[in_arm]
    viability = {g: v[in_arm] for g, v in viability.items()}
    print(f"=== {int(in_arm.sum())}/{n_total} sampled cells are inside Armenia; the other "
          f"{n_total - int(in_arm.sum())} (Georgia/Azerbaijan/Turkey/Iran/Nakhchivan) are dropped ===", flush=True)

    n = len(lats)
    C = len(scenario_names)
    J = len(group_names) * len(RANKED)
    print(f"=== {n} real units, {J} options ({len(group_names)} groups x {len(RANKED)} ranked restoration methods), "
          f"{C} scenarios ===", flush=True)

    print("=== Eligibility: inside Armenia + not protected (WDPA) + not human-modified/water ===", flush=True)
    eligible_mask = load_eligibility(lats, lons)
    print(f"  {eligible_mask.sum()}/{n} real units eligible (inside Armenia, not protected, plantable)", flush=True)

    # A unit is the grid cell its node stands for (stride x stride CHELSA pixels of 30 arcseconds), treated whole or not at all; the node's viability and eligibility are
    # taken to represent the whole cell. Restoration acts on the cell's OPEN land: its area minus forest, woodland and water (national Ecosystem Map, averaged over the cell).
    stride = DENSE_STRIDE if dense else GRID_STRIDE
    px_lat, px_lon = (BBOX[3] - BBOX[1]) / CHELSA_GRID_SHAPE[0], (BBOX[2] - BBOX[0]) / CHELSA_GRID_SHAPE[1]
    cell_ha = cell_area_ha(lats, stride * px_lat, stride * px_lon)
    mapdir = Path(__file__).resolve().parent.parent / "ui" / "assets" / "map"
    raster = lambda name: np.asarray(Image.open(mapdir / f"{name}.bin").convert("L"), dtype=float)
    shares = cell_land_shares(lats, lons, stride * px_lat, stride * px_lon, raster("region"),
                              {"forest": raster("forest") / 255.0, "woodland": raster("woodland") / 255.0, "water": raster("water") / 255.0}, json.load(open(mapdir / "grid.json")))
    no_cover = ~np.isfinite(shares["forest"])
    covered = {k: np.nan_to_num(shares[k], nan=1.0) for k in ("forest", "woodland", "water")}        # a cell the land-cover raster does not cover has no known open land
    open_share = np.clip(1.0 - covered["forest"] - covered["woodland"] - covered["water"], 0.0, 1.0)
    open_ha = cell_ha * open_share
    candidate = eligible_mask & (open_share >= MIN_OPEN_SHARE)
    print(f"  cell: {cell_ha.min():,.0f} to {cell_ha.max():,.0f} ha; open land {open_ha.sum():,.0f} ha of {cell_ha.sum():,.0f} ha in all cells; {int(no_cover.sum())} cells outside the land-cover raster", flush=True)
    print(f"  {int(candidate.sum())}/{int(eligible_mask.sum())} eligible cells have at least {MIN_OPEN_SHARE:.0%} open land; candidate open land {open_ha[candidate].sum():,.0f} ha", flush=True)

    benefit = np.zeros((n, J, C))
    cost = np.zeros((n, J))
    area = np.zeros((n, J))                                    # hectares each option acts on in each cell
    eligible = np.zeros((n, J), dtype=bool)
    option_labels = []
    j = 0
    for g in group_names:
        for interv in RANKED:
            benefit[:, j, :] = viability[g] * VALUE_PER_HA_YEAR * INCREMENTAL_SHARE
            cost[:, j] = interv["cost_per_ha"]
            area[:, j] = open_ha
            eligible[:, j] = candidate
            option_labels.append({"group": g, "intervention": interv["name"], "cost_basis": interv["cost_basis"]})
            j += 1

    results = {
        "run_date": __import__("datetime").date.today().isoformat(),
        "scenario_source": "future_projections" if future is not None else "refugium_2019_fallback",
        "n_scenarios": C,
        "scenario_names": scenario_names,
        "n_units": n,
        "n_eligible_units": int(candidate.sum()),
        "n_eligible_before_land_cover": int(eligible_mask.sum()),
        "min_open_share": MIN_OPEN_SHARE,
        "unit_area_ha": {"min": float(cell_ha.min()), "mean": float(cell_ha.mean()), "max": float(cell_ha.max()),
                         "basis": f"the grid cell a node stands for: {stride} x {stride} CHELSA pixels of 30 arcseconds, treated whole or not at all"},
        "eligible_area_ha": float(open_ha[candidate].sum()),
        "intervention_ids": [i["name"] for i in RANKED],
        "not_ranked": [{"name": i["name"], "role": i["role"], "reason": NOT_RANKED_REASON[i["role"]]} for i in INTERVENTIONS if i["role"] != "restoration"],
        "units": {
            "lat": [round(float(v), 5) for v in lats], "lon": [round(float(v), 5) for v in lons], "area_ha": [round(float(v), 1) for v in cell_ha],
            "open_ha": [round(float(v), 1) for v in open_ha], "forest_share": [round(float(v), 3) for v in covered["forest"]],
            "woodland_share": [round(float(v), 3) for v in covered["woodland"]], "water_share": [round(float(v), 3) for v in covered["water"]],
            "eligible": [bool(v) for v in candidate], "group_ids": list(group_names),
            "cell_deg": {"dlat": float(stride * px_lat), "dlon": float(stride * px_lon)},
            "viability_mean": {g: [round(float(v), 4) for v in viability[g].mean(axis=1)] for g in group_names},
            "viability_worst20": {g: [round(float(v), 4) for v in np.sort(viability[g], axis=1)[:, :max(int(np.ceil(0.2 * C)), 1)].mean(axis=1)] for g in group_names},
        },
        "representative_budget_usd": float(REPRESENTATIVE_BUDGET_USD),
        "solver": {"sweep": {"mip_rel_gap": SWEEP_REL_GAP, "time_limit_s": SWEEP_TIME_S}, "frontier": {"mip_rel_gap": FRONTIER_REL_GAP, "time_limit_s": FRONTIER_TIME_S}, "sensitivity": {"mip_rel_gap": SENS_REL_GAP, "time_limit_s": SENS_TIME_S}},
        "n_sampled_cells_total": n_total,
        "n_cells_dropped_outside_armenia": n_total - n,
        "value_per_ha_year_usd": VALUE_PER_HA_YEAR,
        "incremental_share": INCREMENTAL_SHARE,
        "net_value_per_ha_year_usd": VALUE_PER_HA_YEAR * INCREMENTAL_SHARE,
        "option_labels": option_labels,
        "cost_table": [{k: v for k, v in i.items()} for i in INTERVENTIONS],
        "scope_note": ("A planting unit is the grid cell around a model node (about 3,200 ha on the dense grid), treated whole or not at all, with the node's viability and eligibility "
                        "standing for the whole cell. Only the four restoration methods are ranked, and each acts on the cell's open land (its area minus forest, woodland and water, from the "
                        "national Ecosystem Map averaged over the cell); a cell needs at least 10% open land. Benefit varies only by functional group (REFUGIUM viability) times one national "
                        "value, 417 USD/ha/yr x 75% (the share a restored degraded hectare adds, World Bank 2023); cost varies only by intervention method -- the two are independent by "
                        "construction since no real data or model connects method to survival. Planting uses Armenian cost data (range shown); the other restoration methods use World Bank "
                        "2023 figures (cheapest mix, an all-method average, the most expensive option). Coppicing, pine thinning and wildfire prevention act on existing forest and are not "
                        "ranked: no figure exists for the benefit of maintaining forest. Mining reclamation is not ranked: there is no data on where the sites are. Water use, water caps and "
                        "basins are omitted: no real figure exists for any option's water use. Eligibility = real WDPA (not protected) AND real Ecosystem Map human-modified exclusion "
                        "(not >50% settlements/cropland/buildings/quarries in a local 500 m window)."),
        "budget_sweep": {}, "lambda_frontier_at_representative_budget": {},
    }

    print(f"=== Budget sweep (lambda=0.5, real levels ${BUDGET_LEVELS_USD}) ===", flush=True)
    for budget in BUDGET_LEVELS_USD:
        try:
            plan = robust_portfolio(benefit, area, cost, budget, lam=0.5, eligible=eligible, time_limit=SWEEP_TIME_S, mip_rel_gap=SWEEP_REL_GAP)
            alloc = allocation_summary(plan, area, cost, option_labels)
            results["budget_sweep"][f"${budget:,.0f}"] = {"expected": plan["expected"], "cvar": plan["cvar"], **alloc, "optimal": plan["optimal"], "mip_gap": plan["mip_gap"]}
            print(f"  ${budget:,.0f}: expected={plan['expected']:.1f}, cvar={plan['cvar']:.1f}, {alloc['n_units_planted']} units, "
                  f"{alloc['area_ha']:,.0f} ha, cost ${alloc['cost_usd']:,.0f}; methods {sorted(alloc['by_intervention'])}", flush=True)
        except Exception as e:
            results["budget_sweep"][f"${budget:,.0f}"] = {"error": str(e)}
            print(f"  ${budget:,.0f}: FAILED ({e})", flush=True)

    print(f"=== Mean-CVaR efficient frontier at the real representative budget (${REPRESENTATIVE_BUDGET_USD:,.0f}) ===", flush=True)
    try:
        frontier = efficient_frontier(benefit, area, cost, REPRESENTATIVE_BUDGET_USD, eligible=eligible, time_limit=FRONTIER_TIME_S, mip_rel_gap=FRONTIER_REL_GAP)
        results["lambda_frontier_at_representative_budget"] = {
            "lambdas": frontier["lambdas"].tolist(),
            "expected": frontier["expected"].tolist(),
            "cvar": frontier["cvar"].tolist(),
            "price_of_robustness": frontier["price_of_robustness"],
            "solve_optimal": [bool(p["optimal"]) for p in frontier["plans"]],
            "solve_mip_gap": [p["mip_gap"] for p in frontier["plans"]],
            # the smallest difference between two plans that the solves can tell from noise: the largest (gap x expected benefit) among them
            "resolution_usd": float(max((p["mip_gap"] or 0.0) * p["expected"] for p in frontier["plans"])),
        }
        print(f"  price of robustness: {frontier['price_of_robustness']:.2f}", flush=True)
    except Exception as e:
        results["lambda_frontier_at_representative_budget"] = {"error": str(e)}
        print(f"  FAILED ({e})", flush=True)

    # How much does the uncertain cost of planting matter? Same problem at the representative budget with planting at its low, central and high cost.
    print("=== Sensitivity to the cost of planting ===", flush=True)
    planting_cols = [jj for jj, o in enumerate(option_labels) if o["intervention"] == "degraded_forest_planting"]
    results["planting_cost_sensitivity"] = {}
    for label, key in (("low", "cost_low"), ("central", "cost_per_ha"), ("high", "cost_high")):
        c2 = cost.copy()
        c2[:, planting_cols] = next(i for i in INTERVENTIONS if i["name"] == "degraded_forest_planting")[key]
        try:
            plan = robust_portfolio(benefit, area, c2, REPRESENTATIVE_BUDGET_USD, lam=0.5, eligible=eligible, time_limit=SENS_TIME_S, mip_rel_gap=SENS_REL_GAP)
            results["planting_cost_sensitivity"][label] = {
                "planting_cost_usd_per_ha": float(c2[0, planting_cols[0]]), "expected": plan["expected"], "cvar": plan["cvar"],
                "n_units_planted": int(plan["x"].sum()), "n_units_planting_option": int(plan["x"][:, planting_cols].sum())}
            print(f"  planting at ${c2[0, planting_cols[0]]:,.0f}/ha: expected {plan['expected']:.1f}, units planted {int(plan['x'].sum())}, "
                  f"of which planting {int(plan['x'][:, planting_cols].sum())}", flush=True)
        except Exception as e:
            results["planting_cost_sensitivity"][label] = {"error": str(e)}

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.dump(results, sort_keys=False, default_flow_style=False))
    print(f"=== Wrote {out_path} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
