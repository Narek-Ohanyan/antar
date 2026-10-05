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
    m = re.search(r'<div class="hero">.*?<h1>(.*?)</h1>', app, re.S)
    assert html.unescape(m.group(1)) == FULL


def test_no_other_name_for_the_project_in_the_ui():
    banned = ["refugia explorer", "Armenia reforestation refugia", "Reforestation Predictive Framework", "FORACCA"]
    for f in list(UI.glob("*.html")) + list(UI.glob("*.js")) + list(UI.glob("data/*.json")):
        text = f.read_text()
        for b in banned:
            assert b.lower() not in text.lower(), f"{f.name} contains '{b}'"
