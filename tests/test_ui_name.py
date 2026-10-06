"""The project is named exactly 'ANTAR — Assessment of Niche, Treeline & Analogue Refugia' -- in the browser title and in the
hero heading -- and no other name is invented for it (user's standing rule)."""
import html
import re
from pathlib import Path

UI = Path(__file__).resolve().parent.parent / "ui"
FULL = "ANTAR — Assessment of Niche, Treeline & Analogue Refugia"


def test_browser_title_is_the_full_name():
    m = re.search(r"<title>(.*?)</title>", (UI / "index.html").read_text(), re.S)
    assert html.unescape(m.group(1)) == FULL


def test_hero_heading_is_the_full_name():
    app = (UI / "app.js").read_text()
    m = re.search(r'<h1 id="hero-title">(.*?)</h1>', app, re.S)
    assert html.unescape(m.group(1)) == FULL


def test_no_other_name_for_the_project_in_the_ui():
    banned = ["refugia explorer", "Armenia reforestation refugia", "Reforestation Predictive Framework"]
    for f in list(UI.glob("*.html")) + list(UI.glob("*.js")) + list(UI.glob("data/*.json")):
        text = f.read_text()
        for b in banned:
            assert b.lower() not in text.lower(), f"{f.name} contains '{b}'"


def test_foracca_appears_only_in_the_acknowledgments_and_with_its_official_name():
    """FORACCA is a different project (Forest Restoration and Climate Change in Armenia); the author's acknowledgments name it, nothing else may."""
    for f in list(UI.glob("*.html")) + list(UI.glob("*.json")) + list(UI.glob("data/*.json")) + [UI / "map.js", UI / "place.js", UI / "charts.js", UI / "interp.js"]:
        assert "FORACCA" not in f.read_text(), f.name
    app = (UI / "app.js").read_text()
    start = app.index("function renderAck()")
    end = app.index("/* ---------- ", start)
    assert "FORACCA" not in app[:start] and "FORACCA" not in app[end:]
    ack = app[start:end]
    assert "Forest Restoration and Climate Change in Armenia (FORACCA)" in ack
    assert "Recreation" not in ack                                    # the project's official name is Restoration
    assert "Swiss Agency for Development and Cooperation" in ack
