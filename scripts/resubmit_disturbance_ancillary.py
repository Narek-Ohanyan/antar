"""Re-export the harvest/fire layer with the Hansen no-loss masking fixed (see antar.io.gee_export.export_disturbance_ancillary).

The first export read 0 ("disturbed") at every pixel without a recorded loss, so MNEME's dieback rule could never fire. This submits the corrected export to the same Drive folder under a new
name (``antar_disturbance_ancillary_v2``), on the same grid as the terrain layer; the old file is left in place and is no longer used. Prints the task id; the export takes a few minutes.
"""
import json
import sys
from pathlib import Path

import ee
import yaml
from google.oauth2.credentials import Credentials
import google.auth.transport.requests as gareq

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.gee_export import export_disturbance_ancillary  # noqa: E402

GCP_PROJECT = "pure-highlander-495708-a9"
DESCRIPTION = "antar_disturbance_ancillary_v2"


def main():
    study_area = yaml.safe_load(open(Path(__file__).resolve().parent.parent / "configs" / "study_area.yaml"))
    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    creds = Credentials(None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
                        client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"])
    creds.refresh(gareq.Request())
    ee.Initialize(creds, project=GCP_PROJECT)
    task = export_disturbance_ancillary(tuple(study_area["bbox_wgs84"]), 2000, 2024, study_area["crs"], study_area["resolution_m"], description=DESCRIPTION)
    print(f"{DESCRIPTION}: task {task.id}, state {task.status()['state']}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
