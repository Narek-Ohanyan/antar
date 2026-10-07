"""ui/raster.js: the checksum a browser computes over a decoded raster is the one the build wrote, and every shipped raster matches its recorded checksum."""
import json
import shutil
import subprocess
import zlib
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
MAP = ROOT / "ui" / "assets" / "map"
JS = ROOT / "ui" / "raster.js"
needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def js_adler(rgba, channels):
    code = f"const R=require({str(JS)!r});const d=JSON.parse(require('fs').readFileSync(0,'utf8'));console.log(R.adler32(Uint8ClampedArray.from(d),{channels}))"
    return int(subprocess.run(["node", "-e", code], input=json.dumps(list(map(int, rgba.ravel()))), capture_output=True, text=True, check=True).stdout)


@needs_node
@pytest.mark.parametrize("channels", [1, 2])
def test_the_browser_checksum_equals_zlibs_over_the_same_bytes(channels):
    rng = np.random.default_rng(7)
    rgba = rng.integers(0, 256, size=(40, 50, 4), dtype=np.uint8)
    payload = rgba[..., :channels].tobytes()
    assert js_adler(rgba, channels) == zlib.adler32(payload) & 0xFFFFFFFF


@needs_node
def test_the_checksum_sees_a_single_changed_value():
    rgba = np.zeros((30, 30, 4), dtype=np.uint8)
    a = js_adler(rgba, 1)
    rgba[17, 3, 0] = 1
    assert js_adler(rgba, 1) != a
    rgba[17, 3, 0] = 0
    rgba[17, 3, 1] = 5                                   # a change in a channel that is not hashed does not matter for grey planes
    assert js_adler(rgba, 1) == a and js_adler(rgba, 2) != a


def test_every_shipped_raster_is_a_lossless_png_with_its_recorded_checksum():
    sums = json.loads((MAP / "checksums.json").read_text())
    files = sorted(MAP.glob("*.bin"))
    assert {p.stem for p in files} == set(sums) and len(files) == 9
    for p in files:
        im = Image.open(p)
        assert im.format == "PNG" and im.size == (566, 594), p.name
        payload = np.asarray(im.convert("RGB"))[..., :2].tobytes() if p.stem == "elevation" else np.asarray(im.convert("L")).tobytes()
        assert zlib.adler32(payload) & 0xFFFFFFFF == sums[p.stem], p.name


def test_no_data_raster_is_named_like_a_picture():
    """A host's image optimiser converts *.png, *.jpg and the like; the data must not carry those names."""
    assert not list(MAP.glob("*.png")) and not list(MAP.glob("*.jpg")) and not list(MAP.glob("*.webp"))


def test_the_page_loads_the_guard_before_the_map_code():
    html = (ROOT / "ui" / "index.html").read_text()
    assert html.index("raster.js") < html.index("map.js") < html.index("home.js")
    assert "new Image()" not in (ROOT / "ui" / "map.js").read_text() and "new Image()" not in (ROOT / "ui" / "home.js").read_text()   # nothing decodes a data raster behind the guard's back
    assert "Raster.load" in (ROOT / "ui" / "map.js").read_text() and "Raster.load" in (ROOT / "ui" / "home.js").read_text()
