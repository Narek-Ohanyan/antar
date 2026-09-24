"""Master-grid contract and a v1-style units audit helper."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MASTER_CRS = "EPSG:32638"      # WGS 84 / UTM 38N (Armenia lies in 38N)
MASTER_RES_M = 30.0
ARMENIA_BBOX_WGS84 = (43.4, 38.8, 46.7, 41.4)   # lon_min, lat_min, lon_max, lat_max (buffered)


@dataclass(frozen=True)
class UnitAudit:
    name: str
    ok: bool
    message: str


def audit_dvpd_dtmax_ratio(delta_vpd_kpa: float, delta_tmax_k: float, lo: float = 0.03, hi: float = 0.25) -> UnitAudit:
    """Thermodynamic sanity check on paired deltas.

    At fixed relative humidity, dVPD/dTmax = (1 - RH) * de_s/dT ~ 0.05-0.20 kPa K-1 for
    growing-season conditions in the South Caucasus; RH declines can raise it modestly.
    The v1 training table implies ~0.83 kPa K-1, i.e. at least one variable is mis-scaled.
    """
    r = delta_vpd_kpa / delta_tmax_k if delta_tmax_k != 0 else np.inf
    ok = lo <= r <= hi
    return UnitAudit("dVPD/dTmax", ok, f"ratio = {r:.3f} kPa/K (expected {lo}-{hi})")


def _main() -> None:  # pragma: no cover - thin CLI used by the workflow
    import argparse
    import json
    from pathlib import Path

    ap = argparse.ArgumentParser(description="Write the master-grid contract as JSON.")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({"crs": MASTER_CRS, "resolution_m": MASTER_RES_M, "bbox_wgs84": ARMENIA_BBOX_WGS84}, indent=2))


if __name__ == "__main__":  # pragma: no cover
    _main()
