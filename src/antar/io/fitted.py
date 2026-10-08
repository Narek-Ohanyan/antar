"""Reading and writing the fitted-result files in ``configs/fitted/``, some of which are stored gzip-compressed.

The 45-member scenario result is 57 MB as text and about 6 MB gzipped, so ``future_projections_dense.yaml`` is kept as ``future_projections_dense.yaml.gz``. Code keeps using the plain name:
:func:`resolve` returns the path that exists (the plain file if there is one, else the ``.gz`` one), and :func:`read_yaml` and :func:`write_yaml` open and close either kind. A name listed in
``COMPRESSED`` is always written compressed, so a rerun cannot bring the 57 MB text file back into the repository.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import yaml

COMPRESSED = {"future_projections_dense.yaml"}

try:
    _Loader, _Dumper = yaml.CSafeLoader, yaml.CSafeDumper
except AttributeError:  # pragma: no cover
    _Loader, _Dumper = yaml.SafeLoader, yaml.SafeDumper


def resolve(path) -> Path:
    """The file to read for ``path``: the plain file if it exists, else ``path`` + ".gz"; ``path`` itself if neither does."""
    p = Path(path)
    if p.exists():
        return p
    gz = p.with_name(p.name + ".gz")
    return gz if gz.exists() else p


def write_target(path) -> Path:
    """The file to write for ``path``: the ``.gz`` name when the file is one of ``COMPRESSED`` (or ``path`` already ends in .gz), else ``path``."""
    p = Path(path)
    if p.suffix == ".gz":
        return p
    return p.with_name(p.name + ".gz") if p.name in COMPRESSED else p


def read_yaml(path):
    p = resolve(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8") as f:
        return yaml.load(f, Loader=_Loader)


def write_yaml(path, doc, **dump_kwargs) -> Path:
    """Write ``doc`` (compressed when the name says so) and return the path written. Any plain copy of a compressed file is removed, so there is only ever one."""
    p = write_target(path)
    text = yaml.dump(doc, Dumper=dump_kwargs.pop("Dumper", _Dumper), sort_keys=False, default_flow_style=False, **dump_kwargs)
    if p.suffix == ".gz":
        with gzip.open(p, "wt", encoding="utf-8", compresslevel=9) as f:
            f.write(text)
        plain = p.with_name(p.name[:-3])
        if plain.exists():
            plain.unlink()
    else:
        p.write_text(text)
    return p
