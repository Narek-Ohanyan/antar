"""Pull GHCN-Daily station records for Armenia (53 stations, country code AM) --
a real, freely downloadable observational cross-check/downscaling reference for
CHELSA-daily, per the concept note's data plan (Table tab:data: "daily variability
and QDM reference"). No registration needed; see ncei.noaa.gov/pub/data/ghcn/daily/.

Cloud-first, storage-minimal by design: station files are parsed in place (never
all held in memory at once), written out as compact gzip'd CSV, uploaded to Google
Drive, and the local copy is deleted -- mirroring the S2-VHM pull's pattern
(see configs/manifests/s2_vhm.yaml, IMPLEMENTATION_LOG.md). Nothing beyond the small
station inventory and one station's .dly at a time touches local disk for long.

Units, per GHCN-Daily's own documented convention (not assumed):
  PRCP: tenths of mm -> divided by 10 here to get mm.
  TMAX/TMIN/TAVG: tenths of degC -> divided by 10 here to get degC.
  SNOW, SNWD: already whole mm, not tenths -- left as-is.
Only unflagged values (QFLAG blank, i.e. passed NOAA's own quality control) are kept;
flagged values become NaN rather than being silently included.
"""
import csv
import gzip
import io
import json
import sys
from pathlib import Path

import requests

STATIONS_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
DLY_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/{station_id}.dly"
ELEMENTS = {"PRCP", "TMAX", "TMIN", "TAVG", "SNOW", "SNWD"}
TENTHS_ELEMENTS = {"PRCP", "TMAX", "TMIN", "TAVG"}  # divide by 10; SNOW/SNWD are already mm

STAGING_DIR = Path(__file__).resolve().parent.parent / "data" / "_staging_ghcnd"
DRIVE_FOLDER_NAME = "antar_gee_exports"


def fetch_armenia_stations() -> list[dict]:
    text = requests.get(STATIONS_URL, timeout=30).text
    stations = []
    for line in text.splitlines():
        if not line.startswith("AM"):
            continue
        stations.append({
            "station_id": line[0:11].strip(),
            "lat": float(line[12:20]),
            "lon": float(line[21:30]),
            "elevation_m": float(line[31:37]),
            "name": line[41:71].strip(),
        })
    return stations


def parse_dly(text: str, station_id: str):
    """Yields (date_iso, element, value) tuples, unit-converted, QC-filtered."""
    for line in text.splitlines():
        if len(line) < 269 or not line.startswith(station_id):
            continue
        element = line[17:21]
        if element not in ELEMENTS:
            continue
        year, month = line[11:15], line[15:17]
        for day in range(1, 32):
            start = 21 + (day - 1) * 8
            raw = line[start:start + 5]
            qflag = line[start + 6]
            try:
                value = int(raw)
            except ValueError:
                continue
            if value == -9999 or qflag != " ":
                continue
            try:
                date_iso = f"{year}-{month}-{day:02d}"
                import datetime
                datetime.date.fromisoformat(date_iso)  # validates real calendar day
            except ValueError:
                continue
            if element in TENTHS_ELEMENTS:
                value = value / 10.0
            yield date_iso, element, value


def pull_station(station_id: str) -> list[tuple]:
    resp = requests.get(DLY_URL.format(station_id=station_id), timeout=60)
    resp.raise_for_status()
    return list(parse_dly(resp.text, station_id))


def get_drive_client():
    import ee
    from google.oauth2.credentials import Credentials
    import google.auth.transport.requests as gareq
    from googleapiclient.discovery import build

    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    creds = Credentials(
        None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
        client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"],
    )
    creds.refresh(gareq.Request())
    return build("drive", "v3", credentials=creds)


def find_or_create_drive_folder(drive, name: str) -> str:
    q = f"name='{name}' and mimeType='application/vnd.google-apps.folder' and trashed=false"
    res = drive.files().list(q=q, fields="files(id,name)").execute()
    files = res.get("files", [])
    if files:
        return files[0]["id"]
    meta = {"name": name, "mimeType": "application/vnd.google-apps.folder"}
    return drive.files().create(body=meta, fields="id").execute()["id"]


def upload_to_drive(drive, folder_id: str, local_path: Path) -> str:
    from googleapiclient.http import MediaFileUpload
    media = MediaFileUpload(str(local_path), mimetype="application/gzip", resumable=False)
    meta = {"name": local_path.name, "parents": [folder_id]}
    f = drive.files().create(body=meta, media_body=media, fields="id").execute()
    return f["id"]


def main():
    import hashlib

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    stations = fetch_armenia_stations()
    print(f"=== {len(stations)} Armenia stations found ===", flush=True)

    combined_path = STAGING_DIR / "ghcnd_armenia.csv.gz"
    n_rows = 0
    with gzip.open(combined_path, "wt", newline="") as gz:
        writer = csv.writer(gz)
        writer.writerow(["station_id", "date", "element", "value"])
        for i, st in enumerate(stations, 1):
            try:
                rows = pull_station(st["station_id"])
            except requests.HTTPError as e:
                print(f"{st['station_id']} ({st['name']}): FAILED ({e})", flush=True)
                continue
            for date_iso, element, value in rows:
                writer.writerow([st["station_id"], date_iso, element, value])
            n_rows += len(rows)
            print(f"[{i}/{len(stations)}] {st['station_id']} ({st['name']}): {len(rows)} obs", flush=True)

    checksum = hashlib.sha256(combined_path.read_bytes()).hexdigest()
    size_mb = combined_path.stat().st_size / 1e6
    print(f"=== {n_rows} total observations, {size_mb:.1f} MB compressed, sha256={checksum} ===", flush=True)

    # Station inventory saved alongside (tiny, kept for the manifest's own reference).
    inventory_path = STAGING_DIR / "ghcnd_armenia_stations.json"
    inventory_path.write_text(json.dumps(stations, indent=2))

    print("Uploading to Google Drive...", flush=True)
    drive = get_drive_client()
    folder_id = find_or_create_drive_folder(drive, DRIVE_FOLDER_NAME)
    data_file_id = upload_to_drive(drive, folder_id, combined_path)
    inv_file_id = upload_to_drive(drive, folder_id, inventory_path)
    print(f"Uploaded: data file_id={data_file_id}, inventory file_id={inv_file_id}", flush=True)

    result = {
        "n_stations": len(stations), "n_observations": n_rows, "size_mb": round(size_mb, 2),
        "sha256": checksum, "drive_folder": DRIVE_FOLDER_NAME,
        "data_file_id": data_file_id, "inventory_file_id": inv_file_id,
    }
    (STAGING_DIR / "pull_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))

    # Storage-minimal by design: nothing stays on local disk once it's safely on Drive.
    combined_path.unlink()
    inventory_path.unlink()
    print("Local staging files deleted (data now lives on Drive only).", flush=True)


if __name__ == "__main__":
    sys.exit(main())
