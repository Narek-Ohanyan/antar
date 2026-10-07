"""Every photograph on the site has an author, a licence we can redistribute under, and a link to the source; none is shown without them."""
import json
import re
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
PHOTOS = ROOT / "ui" / "assets" / "photos"
CREDITS = json.loads((PHOTOS / "credits.json").read_text())
ALLOWED = {"CC BY 4.0", "CC BY 3.0", "CC BY-SA 4.0", "CC BY-SA 3.0", "CC0", "Public domain"}
FIELDS = ("file", "caption", "place", "alt", "date", "author", "licence", "licence_url", "source_page", "source", "original_px", "changes")


def test_every_image_file_has_exactly_one_credit_and_every_credit_has_its_files():
    files = {re.sub(r"-(800|1600)\.jpg$", "", f.name) for f in PHOTOS.glob("*.jpg")}
    assert files == {c["file"] for c in CREDITS} and len(CREDITS) == len({c["file"] for c in CREDITS})
    for c in CREDITS:
        for w in (800, 1600):
            with Image.open(PHOTOS / f"{c['file']}-{w}.jpg") as im:
                assert im.width == w, (c["file"], w)
                assert abs(im.height / im.width - c["original_px"][1] / c["original_px"][0]) < 0.003      # resized, not cropped, on disk


def test_every_credit_is_complete_and_the_licence_allows_reuse_with_attribution():
    for c in CREDITS:
        for f in FIELDS:
            assert c.get(f) not in (None, "", []), (c["file"], f)
        assert c["licence"] in ALLOWED, c
        assert c["licence_url"].startswith("https://creativecommons.org/") or c["licence"] == "Public domain"
        assert c["source_page"].startswith("https://commons.wikimedia.org/wiki/File:")
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", c["date"])
        assert len(c["alt"]) > 30                                        # a description, not a label


def test_images_are_light_enough_for_a_web_page():
    assert sum(f.stat().st_size for f in PHOTOS.glob("*.jpg")) < 3_000_000
    assert all(f.stat().st_size < 700_000 for f in PHOTOS.glob("*-1600.jpg"))


def test_the_build_puts_the_credits_in_the_manifest_and_the_acknowledgments_page_lists_them():
    manifest = json.loads((ROOT / "ui" / "data" / "manifest.json").read_text())
    assert [p["file"] for p in manifest["photos"]] == [c["file"] for c in CREDITS]
    app = (ROOT / "ui" / "app.js").read_text()
    ack = app[app.index("function renderAck"):app.index("/* ---------- ", app.index("function renderAck"))]
    assert "state.M.photos" in ack and "ph.licence_url" in ack and "ph.source_page" in ack
    index = (ROOT / "ui" / "index.html").read_text()
    assert "Photographs belong to their authors" in index and "Creative Commons" in index


def test_engine_states_are_written_once_and_valid():
    meth = json.loads((ROOT / "ui" / "data" / "methodology.json").read_text())
    assert all(e["state"] in {"run", "reduced", "built"} for e in meth["engines"])
    assert {e["id"]: e["state"] for e in meth["engines"]}["mneme"] == "built"           # no fitted model: it must never read "run"


def test_the_relief_artwork_exists_and_is_the_terrain_not_a_picture():
    svg = (ROOT / "ui" / "assets" / "relief.svg").read_text()
    assert 'class="minor"' in svg and 'class="major"' in svg and 'class="outline"' in svg and len(svg) < 200_000
    assert "<image" not in svg and "base64" not in svg
