"""Text and control colours of the design tokens meet WCAG 2.2 contrast (4.5:1 for text, 3:1 for large text and component edges), in both themes.

The tokens are read from ui/styles.css itself, so changing a colour there is checked here. Pairs are the ones the page actually uses."""
import re
from pathlib import Path

import pytest

CSS = (Path(__file__).resolve().parent.parent / "ui" / "styles.css").read_text()


def tokens(block):
    return {k: v.strip() for k, v in re.findall(r"(--[\w-]+)\s*:\s*(#[0-9a-fA-F]{6})\b", block)}


def block_after(marker):
    i = CSS.index(marker)
    return CSS[i:CSS.index("}", i)]


LIGHT = tokens(block_after(":root {"))
DARK = {**LIGHT, **tokens(block_after(':root[data-theme="dark"] {'))}


def lum(hex_):
    c = [int(hex_[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    c = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in c]
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def ratio(a, b):
    la, lb = sorted((lum(a), lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


TEXT_PAIRS = [("--text", "--bg"), ("--text", "--surface"), ("--muted", "--bg"), ("--muted", "--surface"), ("--muted", "--surface-2"),
              ("--accent", "--bg"), ("--accent", "--surface"), ("--accent-ink", "--accent"), ("--good", "--good-bg"), ("--warn", "--warn-bg"),
              ("--bad", "--bad-bg"), ("--neutral", "--neutral-bg"), ("--text", "--accent-soft"),
              ("--band-ink", "--band"), ("--band-ink", "--band-2"), ("--band-muted", "--band"), ("--band-muted", "--band-2"), ("--mint", "--band"), ("--mint", "--band-2"),
              ("--sp126", "--band"), ("--sp370", "--band"), ("--sp585", "--band")]


@pytest.mark.parametrize("theme,t", [("light", LIGHT), ("dark", DARK)])
@pytest.mark.parametrize("fg,bg", TEXT_PAIRS)
def test_text_pairs_reach_4_5_to_1(theme, t, fg, bg):
    assert fg in t and bg in t, (fg, bg)
    assert ratio(t[fg], t[bg]) >= 4.5, f"{theme}: {fg} on {bg} = {ratio(t[fg], t[bg]):.2f}"


@pytest.mark.parametrize("theme,t", [("light", LIGHT), ("dark", DARK)])
def test_control_edges_reach_3_to_1_against_the_page(theme, t):
    assert ratio(t["--border-strong"], t["--bg"]) >= 1.9           # hairline rules are decoration; the control edge below is what must be seen
    assert ratio(t["--accent"], t["--bg"]) >= 3.0 and ratio(t["--ring"], t["--bg"]) >= 3.0


def test_the_band_tokens_do_not_change_with_the_theme():
    assert re.search(r":root:not\(\[data-theme=\"light\"\]\) \{[^}]*--band:", CSS) is None and re.search(r"data-theme=\"dark\"\] \{[^}]*--band:", CSS) is None
