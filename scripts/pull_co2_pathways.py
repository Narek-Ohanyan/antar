"""Pull the standard published CO2 concentration trajectories (RCMIP v5.1.0,
Meinshausen et al. 2020) for the scenarios this project actually uses: the CMIP5/RCP
pathways backing the CORDEX ensemble on hyperion (rcp26/45/60/85) and the CMIP6/SSP
pathways backing CHELSA-BIOCLIM+/NASA-GDDP-CMIP6/ISIMIP3b (ssp126/245/370/585), plus
historical, all in one small table -- so the RCP<->SSP correspondence this project's
ROADMAP.md (S:0) requires can be read directly off real forcing numbers instead of
assumed from matching pathway numbers (RCP2.6 and SSP1-2.6 are not the same thing).

Tiny output (9 rows) -- kept directly in configs/, not Drive; the raw 21 MB multi-
model/multi-scenario source file is discarded immediately after filtering.
"""
import csv
import io
import sys
from pathlib import Path

import requests

SOURCE_URL = "https://zenodo.org/record/4589756/files/rcmip-concentrations-annual-means-v5-1-0.csv"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "co2_concentration_pathways.csv"
KEEP_SCENARIOS = {"historical", "rcp26", "rcp45", "rcp60", "rcp85", "ssp126", "ssp245", "ssp370", "ssp585"}


def main():
    resp = requests.get(SOURCE_URL, timeout=60)
    resp.raise_for_status()
    reader = csv.DictReader(io.StringIO(resp.text))
    year_cols = [c for c in reader.fieldnames if c.isdigit() and 1979 <= int(c) <= 2100]
    rows = [
        row for row in reader
        if row["Variable"] == "Atmospheric Concentrations|CO2"
        and row["Region"] == "World"
        and row["Scenario"] in KEEP_SCENARIOS
    ]
    print(f"{len(rows)} scenario rows matched", flush=True)

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["scenario", "mip_era", "model", "unit"] + year_cols)
        for row in rows:
            writer.writerow(
                [row["Scenario"], row["Mip_Era"], row["Model"], row["Unit"]]
                + [row.get(y, "") for y in year_cols]
            )
    print(f"Wrote {OUT_PATH} ({OUT_PATH.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
