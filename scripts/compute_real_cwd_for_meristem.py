"""Real, water-balance climatic water deficit (CWD) at the 3,524 occurrence and background points the MERISTEM niche models are fitted on, once per functional group.

The niche models used ``cwd_approx = petmean - bio12`` (CHELSA-BIOCLIM+) as a stand-in for CWD. TOPOHYDRO's water balance now exists, so this computes the real thing at every point with
exactly the pipeline that produced the grid results: terrain and soils of the point, ERA5-Land and CHELSA-daily for 2019, the seasonal atmosphere (ISIMIP3a shape; an earlier version of
this script held the ERA5-Land annual means constant through the year, which removed most of the summer deficit and left the file invalid), and the group's own rooting depth
(Canadell et al. 1996), because the soil water store is the one thing in the water balance that depends on the species group.

Reads Drive (terrain, soils, ERA5-Land tiles) and the local CHELSA-daily series; writes data/_real_cwd_for_meristem.npz with one CWD array per group (NaN where a point has no
usable soil or ERA5-Land value).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_topohydro_grid import (  # noqa: E402
    DATA_DIR, ROOTING_DEPTH_MM_BY_GROUP, YEAR, atmosphere_label, compute_forcing_for_year_multi_group, extract_static_points,
)

CACHE_PATH = DATA_DIR / "_cache_niche_features.npz"
OUT_PATH = DATA_DIR / "_real_cwd_for_meristem.npz"

# Species to functional group, from configs/species_traits.csv's own `example_taxa` column.
SPECIES_TO_FUNCTIONAL_GROUP = {
    "Fagus orientalis": "mesic_diffuse_porous_broadleaf",
    "Carpinus betulus": "mesic_diffuse_porous_broadleaf",
    "Quercus macranthera": "ring_porous_oak",
    "Quercus iberica": "ring_porous_oak",
    "Pinus kochiana": "pine",
    "Juniperus polycarpos": "juniper_arid_conifer",
    "Juniperus excelsa": "juniper_arid_conifer",
}


def load_all_points():
    d = np.load(CACHE_PATH, allow_pickle=True)
    groups, all_ll = [], []
    for k in d["keys"]:
        ll = d[f"{k}_ll"]
        groups.extend([k] * len(ll))
        all_ll.append(ll)
    all_ll = np.concatenate(all_ll, axis=0)
    return np.array(groups), all_ll[:, 0], all_ll[:, 1]            # label of each point (a species or "background"), latitudes, longitudes


def main():
    labels, lats, lons = load_all_points()
    n = len(lats)
    print(f"=== {n} points (presence of 7 species + shared background): lat {lats.min():.2f}-{lats.max():.2f}, lon {lons.min():.2f}-{lons.max():.2f} ===", flush=True)
    static = extract_static_points(lats, lons)
    by_group = compute_forcing_for_year_multi_group(static, YEAR, ROOTING_DEPTH_MM_BY_GROUP)
    cwd = {g: np.array([np.nan if c is None else float(c.cwd_mm["pm_fao56"]) for c in cells]) for g, cells in by_group.items()}
    for g, v in cwd.items():
        ok = np.isfinite(v)
        print(f"  {g}: {int(ok.sum())}/{n} points with a real CWD (mean {np.nanmean(v):.0f} mm, 95th percentile {np.nanpercentile(v, 95):.0f} mm)", flush=True)
    np.savez_compressed(OUT_PATH, lats=lats, lons=lons, groups=labels, year=YEAR, atmosphere=atmosphere_label(),
                        species_to_functional_group=json.dumps(SPECIES_TO_FUNCTIONAL_GROUP), **{f"real_cwd_mm__{g}": v for g, v in cwd.items()})
    print(f"=== wrote {OUT_PATH.name} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
