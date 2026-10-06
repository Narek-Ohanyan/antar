"""The Method page is one card per engine (ui/data/methodology.json). Its formulas are written as LaTeX between dollar signs and typeset in
the browser by KaTeX (self-hosted); these tests keep the text honest and the formulas valid."""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
UI = ROOT / "ui"
M = json.loads((UI / "data" / "methodology.json").read_text())


def test_there_is_one_card_per_engine_with_status_and_both_views():
    ids = [e["id"] for e in M["engines"]]
    assert ids == ["topohydro", "xylem", "mneme", "meristem", "treeline", "refugium", "aegis"]
    for e in M["engines"]:
        assert e["status"] and e["technical"] and e["plain"], e["id"]


def test_formulas_are_latex_not_unicode_pseudo_maths():
    technical = " ".join(t for e in M["engines"] for t in e["technical"])
    assert technical.count("$") >= 40 and technical.count("$") % 2 == 0
    prose = re.sub(r"\$[^$]+\$", "", technical)
    assert not re.search("[\u0391-\u03a9\u03b1-\u03c9]", prose), "a Greek letter outside a formula: write the formula in LaTeX"


def test_every_formula_typesets_in_katex():
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")
    r = subprocess.run([node, str(Path(__file__).parent / "katex_check.mjs"), str(ROOT)], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-500:]
    assert re.search(r"\d+ formulas, 0 errors", r.stdout)


def test_katex_is_self_hosted_and_wired():
    index = (UI / "index.html").read_text()
    app = (UI / "app.js").read_text()
    assert re.search(r'<script src="math\.js\?v=[0-9a-f]{10}"></script>', index)
    assert "Tex.text(" in app and "Tex.typeset(" in app and "Spec." not in app
    for f in ("katex.min.js", "katex.min.css", "LICENSE"):
        assert (UI / "assets" / "katex" / f).is_file(), f
    assert len(list((UI / "assets" / "katex" / "fonts").glob("*.woff2"))) >= 15
    css = (UI / "assets" / "katex" / "katex.min.css").read_text()
    assert "woff2" in css and ".ttf" not in css and "http" not in css       # one self-hosted font format, no third-party requests


def test_the_flow_chart_matches_its_generator_and_names_every_engine():
    import importlib.util
    spec = importlib.util.spec_from_file_location("arch_svg", ROOT / "scripts" / "build_architecture_svg.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    svg = (UI / "assets" / "architecture.svg").read_text()
    assert svg == mod.build(), "run scripts/build_architecture_svg.py"
    assert "<title" in svg and "<desc" in svg and 'role="img"' in svg
    for engine in ("TOPOHYDRO", "XYLEM", "MNEME", "MERISTEM", "REFUGIUM", "AEGIS"):
        assert engine in svg
    assert "var(--svg-run" in svg and "dash" in svg                       # themed, and inactive flows are dashed
    assert "flowsvg" in (UI / "app.js").read_text()


def test_the_text_matches_what_was_built():
    text = json.dumps(M)
    assert "ISIMIP3a" in text and "seasonal" in text                       # the seasonal, scenario-dependent atmosphere
    assert "nearly flat across scenarios" not in text                      # the artefact finding was retracted
    assert "quantile delta mapping" not in text.lower() or "not" in text.lower()


EARLIER_VERSION = re.compile(r"\bv1\b(?!\.\d)|version[\s~]*(1|one)\b|Output[\s~]*1\.2|\bFVS\b|prior[- ]project|earlier (vulnerability )?model", re.I)


def test_the_methodology_never_mentions_the_earlier_version():
    """Standing instruction from the author: no mention of version 1 in the methodology. Dataset versions that happen to be called v1
    (HydroSHEDS, the GBIF API, a height-model release) are third-party and live in the data manifest, which is not scanned."""
    for f in (UI / "index.html", UI / "data" / "methodology.json", ROOT / "README.md", ROOT / "CITATION.cff", *UI.glob("*.js")):
        text = f.read_text()
        m = EARLIER_VERSION.search(text)
        assert not m, f"{f.name}: {text[max(0, m.start() - 60):m.end() + 60]!r}"
