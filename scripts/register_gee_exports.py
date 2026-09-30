"""Register the completed Earth Engine -> Drive exports (terrain, soils, ERA5-Land
forcing, vegetation state, snow, land tenure, vitality composites, LandTrendr
segmentation, structure, disturbance ancillary) as manifest entries.

These were submitted as ee.batch export tasks in an earlier session and left
unverified -- this script reads their actual state from the Earth Engine API and
the actual files from Google Drive (never assumes a submitted task finished), and
writes one ManifestEntry per Drive file. Files are large (tens of MB to several GB
each, ~85 GB total) and stay on Drive by design (local_path=None) -- storage-minimal,
matching the project's cloud-first convention; nothing here is downloaded.

Checksums: Drive's own MD5 (returned by the Drive API on upload) is recorded in
``notes`` rather than a freshly computed SHA-256, since hashing would mean
downloading every multi-GB file for no integrity benefit Drive doesn't already
provide -- an explicit, documented choice, not an omission.
"""
import json
import sys
from pathlib import Path

import ee
import yaml
from google.oauth2.credentials import Credentials
import google.auth.transport.requests as gareq
from googleapiclient.discovery import build

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.manifest import ManifestEntry, write_manifest

GCP_PROJECT = "pure-highlander-495708-a9"
DRIVE_FOLDER_ID = "1tcE7Uuo0-ZyBqKBy1DUWxK7c3QkIQxFe"  # antar_gee_exports
ACCESS_DATE = "2026-09-30"

DESCRIPTIONS = {
    "antar_terrain": "SRTM 30 m DEM, slope, aspect (Copernicus GLO-30 substitute -- see IMPLEMENTATION_LOG.md for the confirmed Armenia coverage gap in EE's GLO-30 copy).",
    "antar_soils": "SoilGrids 2.0 clay/sand/silt fraction and organic carbon, 0-30cm mean, for Saxton & Rawls (2006) pedotransfer functions.",
    "antar_era5land_forcing": "Annual-mean wind speed, shortwave/longwave radiation, dewpoint temperature and surface pressure, ERA5-Land, 2000-2024. Tiled (multiple files) -- EE split the export by pixel-count limits, not by variable.",
    "antar_vegetation_state": "Growing-season-mean MODIS LAI, annual land-surface phenology (start/end of season, 2001-2024 only -- MCD12Q2 coverage), ESA WorldCover tree-cover class. Re-submitted after fixing a null-image bug for 2000 (MCD12Q2 has no data before 2001).",
    "antar_snow": "Annual MODIS MOD10A1 snow-cover fraction and duration, 2000-2024. Tiled.",
    "antar_land_tenure": "WDPA protected-area boundaries rasterised to a boolean 'protected' band, feeding configs/study_area.yaml's eligibility-mask exclusions.",
    "antar_vitality_composites": "Yearly growing-season (Jul-Aug) median kNDVI (Camps-Valls et al. 2021) from merged HLS v2.0 / Landsat C2 SR, plus per-pixel valid-observation count. Input to antar.hazard.observation.dieback_event. Tiled.",
    "antar_landtrendr_segmentation": "LandTrendr (Kennedy et al. 2010) temporal segmentation of the annual kNDVI series above. Tiled.",
    "antar_structure": "GEDI L2A RH98 (quality-filtered) canopy-top height plus Lang et al. (2023) global 10m canopy height -- training data for MERISTEM's attainable-height model.",
    "antar_disturbance_ancillary": "Per-year 'no_disturbance' boolean (NOT Hansen GFC loss OR MODIS MCD64A1 burn that year) -- the no_disturbance input antar.hazard.observation.dieback_event expects.",
}


def get_credentials():
    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    creds = Credentials(
        None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
        client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"],
    )
    creds.refresh(gareq.Request())
    return creds


def main():
    creds = get_credentials()
    ee.Initialize(creds, project=GCP_PROJECT)
    tasks = ee.data.getTaskList()
    completed_descs = {t["description"] for t in tasks if t.get("state") == "COMPLETED"}
    print(f"Completed EE task descriptions: {sorted(completed_descs)}")

    drive = build("drive", "v3", credentials=creds)
    res = drive.files().list(
        q=f"'{DRIVE_FOLDER_ID}' in parents and trashed=false",
        pageSize=200, fields="files(id,name,size,md5Checksum,modifiedTime)",
    ).execute()
    files = res.get("files", [])
    print(f"{len(files)} files in Drive folder")

    entries = []
    for f in sorted(files, key=lambda x: x["name"]):
        name = f["name"]
        if name.startswith("S2-VHM"):
            continue  # already registered separately in configs/manifests/s2_vhm.yaml
        # Match the file to its export description by stripping any tile suffix / extension.
        base = name.split("-0000")[0].removesuffix(".tif")
        if base not in DESCRIPTIONS:
            print(f"WARNING: no description mapped for {name}, skipping -- register manually")
            continue
        if base not in completed_descs:
            print(f"WARNING: {base} not in EE's own COMPLETED task list, skipping {name}")
            continue
        size_mb = round(int(f.get("size", 0)) / 1e6, 1)
        entries.append(ManifestEntry(
            variable=name.removesuffix(".tif"),
            source="Google Earth Engine batch export (antar.io.gee_export." + base.replace("antar_", "export_") + ")",
            version="n/a",
            url="https://drive.google.com/#folders/" + DRIVE_FOLDER_ID,
            citation="See per-source citations in docs/concept_v2/sec3_data.tex (Table tab:data) and antar/io/gee_export.py docstrings.",
            license="Per underlying source (SRTM/SoilGrids/ERA5-Land/MODIS/WDPA/HLS-Landsat/GEDI -- see gee_export.py docstrings); not uniform, not restated here.",
            access_date=ACCESS_DATE,
            local_path=None,
            checksum_sha256=None,
            spatial_extent="Armenia (clipped)",
            notes=(
                f"{DESCRIPTIONS[base]} Not stored locally (cloud-first by design) -- lives in "
                f"Google Drive folder antar_gee_exports, file id={f['id']}, {size_mb} MB, "
                f"Drive MD5={f.get('md5Checksum', 'n/a')} (recorded instead of a SHA-256: hashing "
                f"would mean downloading a multi-GB file for no integrity benefit Drive doesn't "
                f"already provide)."
            ),
        ))

    write_manifest(entries, "configs/manifests/gee_exports.yaml")
    print(f"Wrote {len(entries)} entries to configs/manifests/gee_exports.yaml")


if __name__ == "__main__":
    main()
