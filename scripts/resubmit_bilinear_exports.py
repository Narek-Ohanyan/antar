"""Re-submit export_soils and export_era5land_forcing with the real .resample('bilinear') fix
just applied to antar.io.gee_export -- configs/resampling_policy.yaml's two identified real
fixes, decided earlier this session, applied to the code, now actually re-run.
"""
import json
import sys
from pathlib import Path

import ee
import yaml
from google.oauth2.credentials import Credentials
import google.auth.transport.requests as gareq

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.gee_export import export_soils, export_era5land_forcing  # noqa: E402

GCP_PROJECT = "antar-armenia-2"


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
    study_area = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "study_area.yaml"))
    bbox = tuple(study_area["bbox_wgs84"])
    crs = study_area["crs"]
    scale_m = study_area["resolution_m"]

    creds = get_credentials()
    ee.Initialize(creds, project=GCP_PROJECT)

    t1 = export_soils(bbox, crs, scale_m)
    print(f"soils (bilinear): {t1.id}, status: {t1.status()['state']}", flush=True)

    t2 = export_era5land_forcing(bbox, 2000, 2024, crs, scale_m)
    print(f"era5land_forcing (bilinear): {t2.id}, status: {t2.status()['state']}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
