"""Placeholder stand-structure input, pending MERISTEM (concept note Sec. 5.4).

Forest ecology enters TOPOHYDRO through leaf area index (LAI), which should
follow attainable height H(a) from MERISTEM's growth model. MERISTEM is built
after TOPOHYDRO in the engine build order, so until it exists this module
loads a clearly labelled, config-driven placeholder LAI per functional group
(``configs/stand_defaults.yaml``) instead of hardcoding a number in code.
"""
from __future__ import annotations

from pathlib import Path

import yaml

PLACEHOLDER_STATUS = "ILLUSTRATIVE_PLACEHOLDER"


def load_placeholder_lai(config_path: str | Path) -> dict[str, float]:
    """Return {functional_group: LAI} from a stand_defaults.yaml-shaped file.

    Raises if the file is not explicitly labelled as a placeholder, so a real
    MERISTEM-derived LAI source can never be silently mistaken for this stub.
    """
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    if cfg.get("status") != PLACEHOLDER_STATUS:
        raise ValueError(
            f"stand-coupling config at {config_path} must be labelled status: {PLACEHOLDER_STATUS} "
            "until MERISTEM supplies real LAI"
        )
    return {str(k): float(v) for k, v in cfg["lai_by_functional_group"].items()}
