"""Data manifests: every dataset enters through one of these, never a hardcoded path.

A manifest entry records what a variable's source *is* (dataset, version, URL,
licence, the citation its provider requires) and what was actually fetched
(access date, checksum of the local file) -- source and provenance in one
place, so no raster is ever used without knowing where it came from. Nothing
here downloads anything; that is a deliberate boundary, matching the "no
silent downloads" convention.
"""
from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field
from pathlib import Path

import yaml


@dataclass(frozen=True)
class ManifestEntry:
    variable: str                  # e.g. "pr", "tas", "tasmax", "vpd_reference"
    source: str                    # dataset name, e.g. "CHELSA V2.1"
    version: str                   # e.g. "2.1"
    url: str                       # provider page or direct URL the file was retrieved from
    citation: str                  # citation string the provider requires
    license: str
    access_date: str               # ISO 8601, e.g. "2026-09-24"
    local_path: str | None = None  # path under the gitignored data/ directory, if fetched
    checksum_sha256: str | None = None
    spatial_extent: str = "Armenia (clipped)"
    notes: str = ""


def sha256_checksum(path: str | Path, chunk_size: int = 1 << 20) -> str:
    """SHA-256 of a local file, streamed so large rasters do not need to fit in memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(entries: list[ManifestEntry], path: str | Path) -> None:
    payload = [asdict(e) for e in entries]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(yaml.safe_dump(payload, sort_keys=False))


def load_manifest(path: str | Path) -> list[ManifestEntry]:
    raw = yaml.safe_load(Path(path).read_text()) or []
    return [ManifestEntry(**entry) for entry in raw]


def verify_entry(entry: ManifestEntry, base_dir: str | Path = ".") -> bool:
    """True iff ``entry.local_path`` exists under ``base_dir`` and its checksum matches.

    Records provenance only -- it never fetches a missing or mismatched file.
    """
    if entry.local_path is None or entry.checksum_sha256 is None:
        return False
    full = Path(base_dir) / entry.local_path
    if not full.is_file():
        return False
    return sha256_checksum(full) == entry.checksum_sha256
