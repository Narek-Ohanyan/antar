"""The projection charts of the Place explorer enlarge when clicked (or activated with Enter / Space): a native modal dialog, sized to the screen
so that text is drawn at a readable size on a phone, with the numbers under the chart."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

UI = Path(__file__).resolve().parent.parent / "ui"


def test_cards_are_keyboard_reachable_buttons_that_open_the_dialog():
    place = (UI / "place.js").read_text()
    assert 'class="card zoomable" role="button" tabindex="0"' in place
    assert 'aria-label="Enlarge the chart:' in place
    assert "Charts.zoom(" in place and 'e.key === "Enter" || e.key === " "' in place
    assert "click a chart to enlarge it" in place


def test_dialog_is_native_closes_on_escape_button_and_backdrop_and_restores_focus():
    charts = (UI / "charts.js").read_text()
    assert 'document.createElement("dialog")' in charts and "showModal()" in charts
    assert "e.target === dlg" in charts and 'aria-label="Close the enlarged view"' in charts
    assert "o.focus(" in charts                                            # focus returns to the chart that was opened
    css = (UI / "styles.css").read_text()
    assert "dialog.zoom" in css and "dialog.zoom::backdrop" in css and "html:has(dialog.zoom[open])" in css
    assert re.search(r"dialog\.zoom \{[^}]*max-width: calc\(100vw - 16px\)", css)       # the UA default would shrink it on a phone


def test_the_enlarged_chart_is_drawn_larger_and_sized_to_the_screen():
    place = (UI / "place.js").read_text()
    assert "window.innerWidth - 16" in place and "big: true" in place
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")
    script = f"""
const fs = require('fs'); eval(fs.readFileSync({str(UI / 'charts.js')!r}, 'utf8') + '; globalThis.Charts = Charts;');
const s = [{{name: 'SSP1-2.6', color: '#0a0', points: [[2019, 1, null, null], [2050, 2, 1, 3], [2100, 3, 2, 4]]}}];
const small = Charts.line({{series: s, xticks: [2019, 2050, 2100], width: 440, height: 270, legendOn: false}});
const big = Charts.line({{series: s, xticks: [2019, 2050, 2100], width: 316, height: 300, legendOn: true, big: true}});
if (!/class="chart"/.test(small) || /r="5.2"/.test(small)) throw new Error('small chart changed');
if (!/class="chart big"/.test(big) || !/r="5.2"/.test(big) || !/stroke-width="3.4"/.test(big)) throw new Error('big chart not drawn larger');
console.log('ok');
"""
    r = subprocess.run([node, "-e", script], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0 and "ok" in r.stdout, r.stdout + r.stderr
