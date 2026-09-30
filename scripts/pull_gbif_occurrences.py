"""Pull real GBIF occurrence records for Armenia's target tree species (the taxa named
in configs/species_traits.csv's example_taxa column) -- MERISTEM's adult-niche Boyce
index needs real presence data, not a placeholder. Uses GBIF's public search API
directly (no registration needed for search; see api.gbif.org/v1/occurrence/search).

Small by nature (~700 records total across 7 species in this small a country) -- kept
directly in configs/, not Drive.
"""
import csv
import json
import sys
import time
import urllib.parse
from pathlib import Path

import requests

SPECIES = [
    "Fagus orientalis", "Carpinus betulus",
    "Quercus macranthera", "Quercus iberica",
    "Pinus kochiana",
    "Juniperus polycarpos", "Juniperus excelsa",
]
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "gbif_occurrences_armenia.csv"
API_BASE = "https://api.gbif.org/v1/occurrence/search"
PAGE_SIZE = 300
# Armenia's real bbox (configs/study_area.yaml) -- country=AM alone does not guarantee
# correct coordinates: a real "null island" (0,0) artifact was caught in this pull before
# this filter was added, a classic GBIF bad-georeferencing pattern, not assumed absent.
BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max


def in_bbox(lat: float, lon: float) -> bool:
    lon_min, lat_min, lon_max, lat_max = BBOX
    return lon_min <= lon <= lon_max and lat_min <= lat <= lat_max


def fetch_species(name: str) -> list[dict]:
    records = []
    offset = 0
    while True:
        params = {
            "country": "AM", "scientificName": name,
            "hasCoordinate": "true", "limit": PAGE_SIZE, "offset": offset,
        }
        resp = requests.get(API_BASE, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        for r in data.get("results", []):
            lat, lon = r.get("decimalLatitude"), r.get("decimalLongitude")
            if lat is None or lon is None or not in_bbox(lat, lon):
                continue  # bad georeferencing (e.g. (0,0) "null island") or genuinely outside Armenia
            records.append({
                "species": name,
                "decimalLatitude": lat,
                "decimalLongitude": lon,
                "eventDate": r.get("eventDate", ""),
                "year": r.get("year", ""),
                "basisOfRecord": r.get("basisOfRecord", ""),
                "coordinateUncertaintyInMeters": r.get("coordinateUncertaintyInMeters", ""),
                "gbifID": r.get("gbifID", ""),
                "datasetKey": r.get("datasetKey", ""),
            })
        if data.get("endOfRecords", True) or len(data.get("results", [])) == 0:
            break
        offset += PAGE_SIZE
        time.sleep(0.2)  # polite pacing, not required by GBIF but costs nothing
    return records


def main():
    all_records = []
    for sp in SPECIES:
        recs = fetch_species(sp)
        print(f"{sp}: {len(recs)} occurrences with coordinates", flush=True)
        all_records.extend(recs)

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "species", "decimalLatitude", "decimalLongitude", "eventDate", "year",
            "basisOfRecord", "coordinateUncertaintyInMeters", "gbifID", "datasetKey",
        ])
        writer.writeheader()
        writer.writerows(all_records)

    print(f"=== {len(all_records)} total occurrences written to {OUT_PATH} ===", flush=True)


if __name__ == "__main__":
    sys.exit(main())
