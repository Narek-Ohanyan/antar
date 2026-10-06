"""Draws the flow chart of the Method page: the six engines, what passes between them, and how far each is built.

Each box is filled by its state in the current release (run / run with reduced scope / built and unit-tested but not applied / not
built) and a dashed arrow means "not active in this release", so the chart cannot claim more than the code does. Colours are CSS
variables with literal fallbacks: when the SVG is inlined in the page it follows the light and dark themes, and it still renders
on its own.

    python3 scripts/build_architecture_svg.py        # writes ui/assets/architecture.svg
"""
from __future__ import annotations

import html
import re
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "ui" / "assets" / "architecture.svg"

S = 52.0                  # px per cm
X0, Y0 = 444.0, 44.0      # pixel position of the origin (x = 0 cm, y = 0 cm)
W, H = 888, 780

STATE_CLASS = {"run": "n-run", "reduced": "n-red", "built": "n-built", "none": "n-none", "driver": "n-drv", "cross": "n-cross"}
TOKEN = re.compile(r"([_^])\{([^}]*)\}")
DESC = ("Six engines and their state. TOPOHYDRO feeds XYLEM, MNEME and MERISTEM. XYLEM's mechanistic hazard goes to REFUGIUM; MNEME and the "
        "area-of-applicability gate are built but not active; MERISTEM's treeline is run and its growth terms are not. REFUGIUM feeds AEGIS, "
        "the decision layer. Uncertainty and validation apply to every arrow.")


def px(x, y):
    return X0 + S * x, Y0 - S * y


def rich(text, size):
    """'h^{mech}' / 'T_{max}' -> tspans with raised or lowered, smaller text."""
    out, pos = [], 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            out.append(html.escape(text[pos:m.start()]))
        up = m.group(1) == "^"
        out.append(f'<tspan dy="{-0.32 * size if up else 0.22 * size:.1f}" font-size="{size * 0.72:.1f}">{html.escape(m.group(2))}</tspan>'
                   f'<tspan dy="{0.32 * size if up else -0.22 * size:.1f}">​</tspan>')
        pos = m.end()
    if pos < len(text):
        out.append(html.escape(text[pos:]))
    return "".join(out)


class Fig:
    def __init__(self):
        self.parts = []

    def box(self, cx, cy, w, h, lines, state, dashed=False, title_size=12.6, body_size=11.2):
        x, y = px(cx - w / 2, cy + h / 2)
        cls = STATE_CLASS[state] + (" dash" if dashed else "")
        self.parts.append(f'<rect class="node {cls}" x="{x:.1f}" y="{y:.1f}" width="{w * S:.1f}" height="{h * S:.1f}" rx="5"/>')
        lh = 14.6
        top = y + h * S / 2 - (len(lines) - 1) * lh / 2 + 4
        for i, (kind, text) in enumerate(lines):
            size = title_size if kind == "t" else body_size
            klass = "ttl" if kind == "t" else ("sub" if kind == "s" else "txt")
            self.parts.append(f'<text class="{klass}" x="{x + w * S / 2:.1f}" y="{top + i * lh:.1f}" text-anchor="middle" font-size="{size}">{rich(text, size)}</text>')

    def arrow(self, pts, dashed=False, label=None, label_at=None):
        d = "M" + " L".join(f"{px(x, y)[0]:.1f},{px(x, y)[1]:.1f}" for x, y in pts)
        self.parts.append(f'<path class="arr{" dash" if dashed else ""}" d="{d}" marker-end="url(#ah)"/>')
        if label:
            lx, ly = px(*label_at)
            self.parts.append(f'<text class="lab" x="{lx:.1f}" y="{ly:.1f}" font-size="10.2">{rich(label, 10.2)}</text>')

    def legend(self, y_cm):
        items = [("run", "run, in the results"), ("reduced", "run, reduced scope"), ("built", "built and tested, not applied"), ("none", "not built")]
        x = 18.0
        _, y = px(0, y_cm)
        for state, label in items:
            self.parts.append(f'<rect class="node {STATE_CLASS[state]}" x="{x:.1f}" y="{y - 11:.1f}" width="18" height="14" rx="3"/>')
            self.parts.append(f'<text class="txt" x="{x + 25:.1f}" y="{y:.1f}" font-size="11">{html.escape(label)}</text>')
            x += 25 + 6.3 * len(label) + 22
        self.parts.append(f'<path class="arr dash" d="M{x:.1f},{y - 4:.1f} L{x + 26:.1f},{y - 4:.1f}" marker-end="url(#ah)"/>')
        self.parts.append(f'<text class="txt" x="{x + 33:.1f}" y="{y:.1f}" font-size="11">dashed: not active</text>')

    def svg(self):
        style = """
.node{stroke:var(--svg-line,#16365C);stroke-width:1.1}
.node.dash{stroke-dasharray:5 3}
.n-run{fill:var(--svg-run,#D9EFE3)}
.n-red{fill:var(--svg-reduced,#FBF0CF)}
.n-built{fill:var(--svg-built,#EAF0F6)}
.n-none{fill:var(--svg-none,#F6DAD3)}
.n-drv{fill:var(--svg-drv,#F4EFE6)}
.n-cross{fill:var(--svg-paper,#FFFFFF);stroke:var(--svg-mute,#4B5563)}
.ttl{font-weight:700;fill:var(--svg-ink,#16365C)}
.txt,.sub{fill:var(--svg-ink,#1E293B)}
.sub{font-style:italic;fill:var(--svg-mute,#4B5563)}
.lab{font-style:italic;fill:var(--svg-mute,#4B5563)}
.arr{fill:none;stroke:var(--svg-mute,#4B5563);stroke-width:1.3}
.arr.dash{stroke-dasharray:5 4}
#ah path{fill:var(--svg-mute,#4B5563)}
text{font-family:'Source Sans 3',system-ui,-apple-system,'Segoe UI',sans-serif}
"""
        defs = '<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0.5 L8,4 L0,7.5 z"/></marker></defs>'
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-labelledby="flow-t flow-d" class="fig-arch">'
                f'<title id="flow-t">How the six engines connect, and how far each is built</title><desc id="flow-d">{html.escape(DESC)}</desc>'
                f"<style>{style}</style>{defs}" + "".join(self.parts) + "</svg>\n")


def build():
    f = Fig()
    f.box(0, 0, 13.2, 0.95, [("t", "Drivers"), ("x", "CHELSA-daily 2019 · ISIMIP3a/3b · ERA5-Land · SRTM 30 m · SoilGrids · stations")], "driver")
    f.box(0, -1.85, 9.4, 1.25, [("t", "TOPOHYDRO — topoclimate and water balance"),
                                 ("x", "T, VPD, R_{n}, PET × 3, snow, soil bucket per species group"),
                                 ("x", "→ CWD, WSI, ψ_{s}, GDD, late frost, GST")], "run")
    f.arrow([(0, -0.48), (0, -1.22)])
    yb = -4.75
    f.box(-5.8, yb, 5.0, 2.15, [("t", "XYLEM — hydraulic failure"),
                                 ("x", "vulnerability curve, stomatal closure,"),
                                 ("x", "g_{min}(T) heat transition, HFI"),
                                 ("x", "50 × 200 trait draws → h^{mech}"),
                                 ("s", "traits: XFT P_{50}, S; others placeholder")], "run")
    f.box(0, yb, 5.0, 2.15, [("t", "MNEME — mortality hazard"),
                              ("x", "cloglog panel model, within/between,"),
                              ("x", "lags; GLM + monotone GBM stack → h^{stat}"),
                              ("s", "dieback labels built; 0 events in the"),
                              ("s", "sampled panel, so not fitted")], "built")
    f.box(5.8, yb, 5.0, 2.15, [("t", "MERISTEM — niche, treeline, growth"),
                                ("x", "niche per species and treeline z_{tl}: run"),
                                ("x", "height H*, growth, establishment:"),
                                ("x", "built, not fitted"),
                                ("s", "niche inputs: CHELSA-BIOCLIM+, GBIF, soils")], "reduced")
    f.arrow([(-3.4, -2.48), (-3.4, -2.8), (-5.8, -2.8), (-5.8, -3.67)], label="ψ_{s}, VPD, T", label_at=(-5.7, -2.7))
    f.arrow([(0, -2.48), (0, -3.67)], label="annual CWD, VPD, GDD, T", label_at=(0.1, -3.05))
    f.arrow([(3.4, -2.48), (3.4, -2.8), (5.8, -2.8), (5.8, -3.67)], label="GST (treeline only)", label_at=(5.9, -2.7))
    f.box(0, -7.1, 5.4, 1.15, [("t", "Area-of-applicability gate"), ("x", "h = w h^{stat} + (1 − w) h^{mech}  (link scale)"), ("s", "built; inactive until h^{stat} exists")], "built", dashed=True)
    f.arrow([(0, -5.83), (0, -6.5)], dashed=True)
    f.arrow([(-3.4, -5.83), (-3.4, -7.1), (-2.7, -7.1)], dashed=True)
    f.box(0, -9.1, 9.4, 1.25, [("t", "REFUGIUM — viability, refugia, analogues"),
                                ("x", "V^{(1)} = 1 − h^{mech}; robust refugium, criterion (a); risk-averse score"),
                                ("s", "other hazards, height factor, criteria (b), (c), analogues: not computed")], "reduced")
    f.arrow([(-5.8, -5.83), (-5.8, -8.1), (-4.0, -8.1), (-4.0, -8.47)], label="h^{mech} (used)", label_at=(-5.7, -7.35))
    f.arrow([(0, -7.68), (0, -8.47)], dashed=True, label="h^{hyd}", label_at=(0.1, -8.1))
    f.arrow([(5.8, -5.83), (5.8, -8.1), (4.0, -8.1), (4.0, -8.47)], dashed=True, label="H*, growth, h^{est}", label_at=(5.9, -7.35))
    f.box(0, -10.85, 9.4, 1.0, [("t", "AEGIS — decision layer"), ("x", "CVaR-robust species × method portfolio (MILP), 45 members")], "reduced")
    f.arrow([(0, -9.73), (0, -10.35)])
    f.box(-4.45, -12.5, 8.1, 1.15, [("t", "Uncertainty"), ("x", "SSP × GCM × horizon × traits: carried (run)"), ("x", "ANOVA and Sobol' partition: built, not applied")], "cross")
    f.box(4.45, -12.5, 8.1, 1.15, [("t", "Validation"), ("x", "niche block CV, synthetic extrapolation test: run"), ("x", "hazard CV, conformal, design-based: built")], "cross")
    f.legend(-13.6)
    return f.svg()


if __name__ == "__main__":
    OUT.write_text(build())
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
