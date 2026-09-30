"""Submit the real TerraClimate export (antar.io.gee_export.export_terraclimate) -- the
independent water-balance benchmark named in the concept note but never pulled this session,
ROADMAP.md S1.2. Follows the same auth/submission pattern as scripts/register_gee_exports.py.
"""
import json
import sys
from pathlib import Path

import ee
import yaml
from google.oauth2.credentials import Credentials
import google.auth.transport.requests as gareq

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.gee_export import export_terraclimate  # noqa: E402

GCP_PROJECT = "antar-armenia-2"  # per ROADMAP.md S2's own guidance: use this one for new EE work


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

    creds = get_credentials()
    ee.Initialize(creds, project=GCP_PROJECT)

    task = export_terraclimate(bbox, crs)
    print(f"Submitted: {task.id}, status: {task.status()}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
