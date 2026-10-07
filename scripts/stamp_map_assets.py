"""Writes ui/assets/map/checksums.json: an Adler-32 checksum of the decoded pixels of every data raster in ui/assets/map/*.bin.

The browser (ui/raster.js) decodes each raster and compares; a mismatch means a host or CDN changed the image (for example converted the PNG to a lossy WebP)
and the map would silently show wrong values. The checksum is over the single grey plane for id and cover rasters and over the R and G bytes (metres = 256*R + G)
for elevation, in row-major order, exactly the bytes ui/raster.js hashes. Run after scripts/build_map_assets.py or scripts/build_dem_asset.py.

    python3 scripts/stamp_map_assets.py
"""
import json
import zlib
from pathlib import Path

import numpy as np
from PIL import Image

MAP = Path(__file__).resolve().parent.parent / "ui" / "assets" / "map"
CHANNELS = {"elevation": 2}


def checksum(path: Path) -> int:
    im = Image.open(path)
    if im.format != "PNG":
        raise SystemExit(f"{path.name} is {im.format}, not a lossless PNG")
    name = path.stem
    if CHANNELS.get(name, 1) == 2:
        payload = np.asarray(im.convert("RGB"))[..., :2].tobytes()
    else:
        payload = np.asarray(im.convert("L")).tobytes()
    return zlib.adler32(payload) & 0xFFFFFFFF


def main():
    out = {p.stem: checksum(p) for p in sorted(MAP.glob("*.bin"))}
    (MAP / "checksums.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {MAP / 'checksums.json'}: {len(out)} rasters")


if __name__ == "__main__":
    main()
