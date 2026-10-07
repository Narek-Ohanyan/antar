"""ui/home.js and ui/palette.js: the pure functions behind the home page, run through node (skipped, not passed, if node is missing)."""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HOME, PALETTE = ROOT / "ui" / "home.js", ROOT / "ui" / "palette.js"
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def node(code, payload=None):
    r = subprocess.run(["node", "-e", f"const H=require({str(HOME)!r}),P=require({str(PALETTE)!r});const d=JSON.parse(require('fs').readFileSync(0,'utf8')||'null');{code}"],
                       input=json.dumps(payload), capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def cone(w=41, h=41, cx=20.0, cy=20.0):
    return [((x - cx) ** 2 + (y - cy) ** 2) ** 0.5 for y in range(h) for x in range(w)]


def test_a_circular_hill_gives_one_closed_loop_at_the_right_radius():
    out = node("const L=H.contours(d.z,d.w,d.h,d.level);console.log(JSON.stringify(L))", {"z": cone(), "w": 41, "h": 41, "level": 10})
    assert len(out) == 1
    loop = out[0]
    assert loop[0] == loop[-1] and len(loop) > 40
    radii = [((x - 20) ** 2 + (y - 20) ** 2) ** 0.5 for x, y in loop]
    assert max(abs(r - 10) for r in radii) < 0.05                      # linear interpolation along a cone is exact to well under a pixel


def test_every_contour_vertex_lies_on_a_cell_edge_and_levels_do_not_mix():
    z = cone()
    out = node("console.log(JSON.stringify([5,10,15].map(l=>H.contours(d.z,d.w,d.h,l))))", {"z": z, "w": 41, "h": 41})
    assert [len(x) for x in out] == [1, 1, 1]
    for level, lines in zip((5, 10, 15), out):
        for x, y in lines[0]:
            assert abs(((x - 20) ** 2 + (y - 20) ** 2) ** 0.5 - level) < 0.1


def test_missing_cells_end_a_line_instead_of_inventing_one():
    z = cone()
    for i in range(41 * 41):
        if i % 41 > 25:                                                  # no data on the east side
            z[i] = None
    out = node("const z=d.z.map(v=>v===null?NaN:v);console.log(JSON.stringify(H.contours(z,41,41,10)))", {"z": z})
    assert len(out) == 1 and out[0][0] != out[0][-1]                    # an open line, not a closed loop
    assert all(x <= 26 for x, y in out[0])                              # nothing is drawn where there is no data


def test_flat_and_empty_fields_give_no_lines():
    assert node("console.log(JSON.stringify(H.contours(d.z,5,5,3)))", {"z": [1.0] * 25}) == []
    assert node("console.log(JSON.stringify(H.contours(d.z,5,5,3)))", {"z": [9.0] * 25}) == []


def test_saddle_cells_are_resolved_without_crossing_lines():
    # 2x2 checkerboard: high, low / low, high: two separate segments whichever way the centre falls
    z = [10, 0, 0, 10]
    low = node("console.log(JSON.stringify(H.contours(d.z,2,2,6.0)))", {"z": z})
    assert len(low) == 2 and all(len(s) == 2 for s in low)


def test_downsample_keeps_means_and_blanks_mostly_empty_blocks():
    z = [1, 3, 5, 7, 2, 4, 6, 8, 10, 10, None, None, 10, 10, None, None]
    out = node("const z=d.z.map(v=>v===null?NaN:v);const r=H.downsample(z,4,4,2);console.log(JSON.stringify({W:r.W,H:r.H,z:Array.from(r.z).map(v=>Number.isNaN(v)?null:v)}))", {"z": z})
    assert out["W"] == 2 and out["H"] == 2
    assert out["z"][0] == pytest.approx((1 + 3 + 2 + 4) / 4) and out["z"][1] == pytest.approx((5 + 7 + 6 + 8) / 4)
    assert out["z"][2] == 10 and out["z"][3] is None                     # 4 of 4 finite -> mean; 0 of 4 -> blank


MANIFEST = {
    "grids": {"validation": {"n_cells": 25}, "dense": {"n_cells": 854}},
    "groups": {"a": {"short": "alpha", "label": "Alpha"}, "b": {"short": "beta", "label": "Beta"}},
    "baseline_viability_2019": {"a": 0.99, "b": 0.90},
    "scenario_summary": {g: {s: {h: {"mean": base - dx * (i + 1), "min": base - dx * (i + 1) - .01, "max": base - dx * (i + 1) + .01, "n_gcms": 5}
                                for i, h in enumerate(("2050", "2080", "2100"))} for s, dx in (("ssp126", .001), ("ssp370", .004), ("ssp585", .01))}
                         for g, base in (("a", .99), ("b", .90))},
    "treeline": {"grid": "validation", "summary": {f"{s}__{h}": {"ensemble_mean_shift_m": m * (i + 1), "ensemble_min_shift_m": m * (i + 1) - 20, "ensemble_max_shift_m": m * (i + 1) + 30, "n_gcms": 5}
                                                   for s, m in (("ssp126", 30), ("ssp370", 100), ("ssp585", 150)) for i, h in enumerate(("2050", "2080", "2100"))}},
}


def model(ssp="ssp585", hz="2100", tl=2670):
    return node("console.log(JSON.stringify(H.headlineModel(d.M,d.ssp,d.hz,d.tl)))", {"M": MANIFEST, "ssp": ssp, "hz": hz, "tl": tl})


def test_headline_numbers_are_the_manifest_numbers():
    m = model()
    a, b = m["groups"]
    assert a["base"] == .99 and a["cur"]["mean"] == pytest.approx(.99 - .03)
    assert a["deltaPp"] == pytest.approx(-3.0) and b["deltaPp"] == pytest.approx(-3.0)
    assert m["treeline"]["mean"] == 450 and m["treeline"]["min"] == 430 and m["treeline"]["max"] == 480
    assert a["series"]["ssp585"][0] == .99 and len(a["series"]["ssp126"]) == 4          # 2019 then three horizons
    assert m["treelineNowM"] == 2670 and m["treelineLaterM"] == 2670 + 450


def test_the_lead_sentence_states_the_selection_and_the_changes():
    lead = model("ssp370", "2050")["lead"]
    assert "SSP3-7.0" in lead and "2050" in lead and "by 100 m" in lead and "alpha −0.4 pp" in lead and "against 2019" in lead
    assert model()["lead"].count("−3.0 pp") == 2


def test_a_missing_scenario_summary_falls_back_to_2019_without_inventing_a_change():
    m = node("const M=JSON.parse(JSON.stringify(d.M));M.scenario_summary=null;M.treeline=null;console.log(JSON.stringify(H.headlineModel(M,'ssp585','2100',null)))", {"M": MANIFEST})
    assert m["treeline"] is None and all(g["cur"] is None and g["deltaPp"] is None for g in m["groups"])
    assert "not available" in m["lead"] and m["treelineLaterM"] is None


def test_trend_chart_draws_every_path_and_is_flat_when_the_change_is_below_the_minimum_span():
    spec = {"series": [{"ssp": s, "values": v, "selected": s == "ssp585"} for s, v in (("ssp126", [1, 1, 1, 1]), ("ssp370", [1, 1, 1, 1]), ("ssp585", [1, 1, 1, .99999]))],
            "band": {"lo": [1, 1, 1, .99999], "hi": [1, 1, 1, 1]}, "marker": 3, "minSpan": 0.04, "label": "Pine survival"}
    svg = node("console.log(JSON.stringify(H.sparkSvg(d)))", spec)
    assert svg.count('class="sp-line ctx') == 2 and svg.count('class="sp-line sel') == 1 and 'class="sp-band"' in svg
    assert 'aria-label="Pine survival"' in svg and svg.count("<circle") == 4
    ys = {round(float(y), 1) for y in re.findall(r'cy="([\d.]+)"', svg)}
    assert len(ys) == 1                                                  # 0.00001 is drawn as no change at all


def test_trend_chart_shows_a_real_change():
    spec = {"series": [{"ssp": "ssp585", "values": [.99, .98, .97, .94], "selected": True}], "marker": 3, "minSpan": 0.04, "label": "x"}
    svg = node("console.log(JSON.stringify(H.sparkSvg(d)))", spec)
    ys = [float(y) for y in re.findall(r'cy="([\d.]+)"', svg)]
    assert ys == sorted(ys) and ys[-1] - ys[0] > 30                      # falling values run down the chart


def test_cards_carry_the_numbers_a_link_to_the_map_and_no_unlabelled_colour_meaning():
    html = node("const m=H.headlineModel(d.M,'ssp585','2100',2670);console.log(JSON.stringify(m.groups.map(g=>H.viabilityCard(g,m,2)).join('')+H.treelineCard(m.treeline,m,2)))", {"M": MANIFEST})
    assert "96.0" in html and "−3.0 pp" in html and "against 2019 (99.0%)" in html
    assert "#/map?q=viab_a&amp;mode" not in html and "#/map?q=viab_a&mode=scen&ssp=ssp585&hz=2100&gcm=ens" in html
    assert "+450" in html and "not a forecast" in html
    assert "▼" in html                                              # the direction is also shown by a glyph, not by colour alone


def test_photo_markup_has_alt_text_credit_licence_and_source():
    p = {"file": "x", "caption": "Cap", "place": "Place", "alt": "A thing", "date": "2020-01-02", "author": "An Author", "licence": "CC BY-SA 4.0",
         "licence_url": "https://creativecommons.org/licenses/by-sa/4.0", "source_page": "https://commons.wikimedia.org/wiki/File:X.jpg", "source": "Wikimedia Commons", "original_px": [4000, 3000]}
    html = node("console.log(JSON.stringify(H.photoFigure(d,0)+H.photoZoomBody(d)))", p)
    for s in ("alt=\"A thing\"", "An Author", "CC BY-SA 4.0", "creativecommons.org/licenses/by-sa/4.0", "commons.wikimedia.org/wiki/File:X.jpg", "resized", "noopener", "x-1600.jpg", "x-800.jpg"):
        assert s in html, s


ITEMS = [{"g": "Pages", "t": "Map", "s": "page"}, {"g": "Pages", "t": "Treeline", "s": "page"}, {"g": "Map layers", "t": "Climatic water deficit — broadleaf", "s": "XYLEM", "k": "cwd"},
         {"g": "Map layers", "t": "Treeline shift vs 2019", "s": "Treeline"}, {"g": "Places", "t": "Syunik marz", "s": "Place explorer"}, {"g": "Places", "t": "Yerevan (city)", "s": "Place explorer"}]


def rank(q, limit=None):
    return [i["t"] for i in node("console.log(JSON.stringify(P.rank(d.items,d.q,d.limit)))", {"items": ITEMS, "q": q, "limit": limit})]


def test_search_needs_every_word_and_prefers_title_starts():
    assert rank("treeline")[:2] == ["Treeline", "Treeline shift vs 2019"]
    assert rank("water broadleaf") == ["Climatic water deficit — broadleaf"]
    assert rank("zzz") == []
    assert rank("syunik") == ["Syunik marz"]


def test_search_ignores_case_and_accents_and_lists_everything_for_an_empty_query():
    assert rank("YEREVAN") == ["Yerevan (city)"] and rank("climatic wäter") == ["Climatic water deficit — broadleaf"]
    assert len(rank("")) == len(ITEMS) and len(rank("", 3)) == 3
    assert rank("cwd") == ["Climatic water deficit — broadleaf"]        # keywords count too


def test_negative_shifts_use_a_proper_minus_in_the_sentence():
    m = node("const M=JSON.parse(JSON.stringify(d.M));M.treeline.summary['ssp126__2050'].ensemble_min_shift_m=-8.4;console.log(JSON.stringify(H.headlineModel(M,'ssp126','2050',2670).lead))", {"M": MANIFEST})
    assert "models: \u22128 to" in m and "-8" not in m


def test_the_csv_holds_every_headline_number_and_parses_back():
    import csv
    import io
    text = node("console.log(JSON.stringify(H.headlineCsv(d.M)))", {"M": MANIFEST})
    assert text.endswith("\r\n") and text.count("\r\n") == 1 + 2 * 9 + 9          # header, 2 groups x 3 paths x 3 horizons, 9 treeline rows
    rows = list(csv.DictReader(io.StringIO(text)))
    assert list(rows[0].keys()) == ["quantity", "group", "emissions_path", "horizon", "baseline_2019", "ensemble_mean", "model_min", "model_max", "n_climate_models", "unit", "grid"]
    a = next(r for r in rows if r["group"] == "alpha" and r["emissions_path"] == "SSP5-8.5" and r["horizon"] == "2100")
    assert float(a["ensemble_mean"]) == pytest.approx(0.96) and float(a["baseline_2019"]) == pytest.approx(0.99) and a["n_climate_models"] == "5" and a["unit"] == "fraction"
    t = next(r for r in rows if r["quantity"] == "climatic treeline shift" and r["emissions_path"] == "SSP3-7.0" and r["horizon"] == "2080")
    assert float(t["ensemble_mean"]) == 200 and float(t["model_min"]) == 180 and float(t["model_max"]) == 230 and t["unit"] == "m"


def test_the_csv_quotes_awkward_text():
    text = node("const M=JSON.parse(JSON.stringify(d.M));M.groups.a.short='odd, \"name\"';console.log(JSON.stringify(H.headlineCsv(M)))", {"M": MANIFEST})
    assert '"odd, ""name"""' in text


def test_each_card_names_the_grid_behind_it_even_when_they_differ():
    html = node("const m=H.headlineModel(d.M,'ssp585','2100',2670,{viability:'dense',treeline:'validation'});console.log(JSON.stringify(m.groups.map(g=>H.viabilityCard(g,m,2)).join('')+'|'+H.treelineCard(m.treeline,m,2)))", {"M": MANIFEST})
    viab, tl = html.split("|")
    assert viab.count("854 nodes") == 2 and "25 nodes" not in viab and tl.count("25 nodes") == 1 and "854 nodes" not in tl
    none = node("const m=H.headlineModel(d.M,'ssp585','2100',2670);console.log(JSON.stringify(H.viabilityCard(m.groups[0],m,2)))", {"M": MANIFEST})
    assert "cgrid" not in none
