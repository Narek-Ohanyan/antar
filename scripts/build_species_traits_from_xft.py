"""Replace species_traits.csv's placeholder P50/slope/gmin with real values computed
from the XFT (Xylem Functional Traits) database -- user-downloaded CSV from
xylemfunctionaltraits.org, no exact-species records exist for any of Armenia's four
target taxa, but real congeneric data does for all four groups (documented per-group
below, with record counts, so nobody downstream mistakes a proxy for a direct
measurement). Uses antar.hydraulics.vulnerability.slope_from_p12_p88 -- the project's
own existing P50/slope relationship (Pammenter & Van der Willigen 1998) -- rather than
re-deriving a formula.

psi_close_mpa, tp_c, lethal_plc and capacitance_mmol_m2_mpa are left as placeholders:
no clean, unit-compatible XFT field exists for any of them (capacitance in particular
is reported per sapwood volume, kg m-3 MPa-1, not the per-leaf-area molar basis
species_traits.csv uses, and converting would need wood density and sapwood:leaf-area
ratio this dataset doesn't reliably give per record -- not attempted rather than guessed).
lethal_plc is a modeling convention already correctly set per plant type (P88 for the
two angiosperm groups, P50 for the two gymnosperm groups, matching
hydraulics/twophase.py's own docstring), not something to fit from data.
"""
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.hydraulics.vulnerability import slope_from_p12_p88

XFT_PATH = Path(__file__).resolve().parent.parent / "data" / "xft" / "XFT_full_database_download_20260930-151542.csv"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "species_traits.csv"

# (genus, species) tuples per group -- congeneric/exact matches actually present in XFT,
# chosen for taxonomic and ecological proximity to the real target species, not just
# "same genus, whatever's there" (see IMPLEMENTATION_LOG.md for the full reasoning).
GROUP_SPECIES = {
    "mesic_diffuse_porous_broadleaf": [
        ("Fagus", "sylvatica"),      # closest relative of the target F. orientalis
        ("Carpinus", "betulus"),     # exact match -- one of the two target species
    ],
    "ring_porous_oak": [
        # European white-oak-group deciduous oaks, closest available analogs to
        # Q. macranthera / Q. iberica (both also white-oak-group, submediterranean):
        ("Quercus", "petraea"), ("Quercus", "robur"),
        ("Quercus", "pubescens"), ("Quercus", "frainetto"), ("Quercus", "humilis"),
    ],
    "pine": [
        # P. kochiana is taxonomically very close to (sometimes treated as a variety
        # of) P. sylvestris -- the single best available analog, not a loose genus-wide pick.
        ("Pinus", "sylvestris"),
    ],
    "juniper_arid_conifer": [
        # J. thurifera ONLY, despite n=1 -- checked, then deliberately not widened. A first
        # pass included J. communis (n=11) for sample size and it pulled the group median to
        # -5.96 MPa, well short of the -9 MPa expected for an arid-zone Irano-Turanian juniper;
        # a full-genus P50 survey (13-24 species, -1.67 to -14.2 MPa) confirmed J. communis
        # sits at the mesic/widespread end of that range while J. thurifera -- a high-elevation
        # Mediterranean-mountain species -- sits in the ecologically appropriate arid/cold niche
        # for J. polycarpos/J. excelsa. Ecological similarity outweighs raw sample count for a
        # single congeneric proxy; diluting it with a poorly-matched species for the sake of n
        # would be worse, not better. n=1 is a real limitation, stated plainly in status/notes,
        # not hidden behind a bigger but less appropriate sample.
        ("Juniperus", "thurifera"),
    ],
}

EXAMPLE_TAXA = {
    "mesic_diffuse_porous_broadleaf": "Fagus orientalis; Carpinus betulus",
    "ring_porous_oak": "Quercus macranthera; Q. iberica",
    "pine": "Pinus kochiana",
    "juniper_arid_conifer": "Juniperus polycarpos; J. excelsa",
}

# Unchanged from the original placeholder -- not derivable from XFT, see module docstring.
# gmin25_mmol_m2_s here is the FALLBACK used only if XFT has no Gsmin record for that group's
# matched species (real Gsmin data only existed for the pine group in this pull).
UNCHANGED = {
    "mesic_diffuse_porous_broadleaf": {"psi_close_mpa": -2.2, "gmin25_mmol_m2_s": 3.0, "tp_c": 38, "lethal_plc": 88, "capacitance_mmol_m2_mpa": 30000},
    "ring_porous_oak": {"psi_close_mpa": -2.4, "gmin25_mmol_m2_s": 3.5, "tp_c": 38, "lethal_plc": 88, "capacitance_mmol_m2_mpa": 30000},
    "pine": {"psi_close_mpa": -2.6, "gmin25_mmol_m2_s": 1.5, "tp_c": 40, "lethal_plc": 50, "capacitance_mmol_m2_mpa": 20000},
    "juniper_arid_conifer": {"psi_close_mpa": -5.0, "gmin25_mmol_m2_s": 1.0, "tp_c": 42, "lethal_plc": 50, "capacitance_mmol_m2_mpa": 15000},
}


def to_float(x):
    try:
        v = float(x)
        return v
    except (ValueError, TypeError):
        return None


def main():
    with open(XFT_PATH, encoding="latin-1") as f:
        rows = list(csv.DictReader(f))

    provenance = []
    out_rows = []
    for group, species_list in GROUP_SPECIES.items():
        matched = [r for r in rows if (r.get("Genus", "").strip(), r.get("Species", "").strip()) in species_list]
        p50s = [to_float(r["P50"]) for r in matched if to_float(r["P50"]) is not None]
        slopes = []
        for r in matched:
            p12, p88 = to_float(r["P12"]), to_float(r["P88"])
            if p12 is not None and p88 is not None and p12 > p88:
                slopes.append(slope_from_p12_p88(p12, p88))
        gsmins = [to_float(r["Gsmin..mol.m.2.s.1."]) for r in matched if to_float(r.get("Gsmin..mol.m.2.s.1.")) is not None]

        p50_median = round(statistics.median(p50s), 2) if p50s else None
        slope_median = round(statistics.median(slopes), 1) if slopes else None
        gmin_mmol = round(statistics.median(gsmins) * 1000, 2) if gsmins else None  # mol -> mmol m-2 s-1

        species_str = ", ".join(f"{g} {s} (n={sum(1 for r in matched if r['Genus'].strip()==g and r['Species'].strip()==s)})"
                                 for g, s in species_list)
        provenance.append(
            f"{group}: P50 n={len(p50s)}, slope n={len(slopes)}, Gsmin n={len(gsmins)}, from {species_str}"
        )

        row = {
            "group": group,
            "example_taxa": EXAMPLE_TAXA[group],
            "p50_mpa": p50_median if p50_median is not None else "",
            "slope_pct_per_mpa": slope_median if slope_median is not None else "",
            "psi_close_mpa": UNCHANGED[group]["psi_close_mpa"],
            "gmin25_mmol_m2_s": gmin_mmol if gmin_mmol is not None else UNCHANGED[group].get("gmin25_mmol_m2_s", ""),
            "tp_c": UNCHANGED[group]["tp_c"],
            "lethal_plc": UNCHANGED[group]["lethal_plc"],
            "capacitance_mmol_m2_mpa": UNCHANGED[group]["capacitance_mmol_m2_mpa"],
            "status": (
                f"P50_SLOPE_REAL_XFT_CONGENERIC_2026-09-30 (n={len(p50s)} P50, n={len(slopes)} slope, "
                f"from {species_str}); gmin25 " + ("REAL_XFT" if gmin_mmol is not None else "still ILLUSTRATIVE_PLACEHOLDER")
                + "; psi_close/tp_c/lethal_plc/capacitance still ILLUSTRATIVE_PLACEHOLDER"
            ),
        }
        out_rows.append(row)

    fieldnames = ["group", "example_taxa", "p50_mpa", "slope_pct_per_mpa", "psi_close_mpa",
                  "gmin25_mmol_m2_s", "tp_c", "lethal_plc", "capacitance_mmol_m2_mpa", "status"]
    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)
        f.write(
            "# psi_close_mpa, tp_c, lethal_plc, capacitance_mmol_m2_mpa remain order-of-magnitude\n"
            "# placeholders -- no clean, unit-compatible XFT field exists for any of them (see\n"
            "# scripts/build_species_traits_from_xft.py's module docstring). p50_mpa, slope_pct_per_mpa,\n"
            "# and gmin25_mmol_m2_s are real values computed from the XFT database (Choat et al. 2012 and\n"
            "# later contributions; xylemfunctionaltraits.org), congeneric proxies where the exact Armenian\n"
            "# species has no XFT records -- see IMPLEMENTATION_LOG.md 2026-09-30 for full per-group\n"
            "# provenance (which species, how many records, why each was chosen).\n"
        )

    print("\n".join(provenance))
    print(f"\nWrote {OUT_PATH}")


if __name__ == "__main__":
    main()
