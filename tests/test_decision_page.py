"""ui/decision.js: cell lookup, the CSV of treated cells and the marz table (node), and the exported portfolio is consistent with the optimiser's own totals."""
import csv as pycsv
import io
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JS = ROOT / "ui" / "decision.js"
DATA = ROOT / "ui" / "data"
needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")

UNITS = {"lat": [40.0, 40.1, 40.2], "lon": [44.0, 44.1, 44.2], "area_ha": [3200.0, 3190.0, 3180.0], "eligible": [True, True, False], "group_ids": ["a", "b"],
         "cell_deg": {"dlat": 0.05, "dlon": 0.05}, "viability_mean": {"a": [0.98, 0.97, 0.9], "b": [0.9, 0.95, 0.9]}, "viability_worst20": {"a": [0.96, 0.95, 0.8], "b": [0.85, 0.9, 0.8]}}


def node(code, payload=None):
    r = subprocess.run(["node", "-e", f"const D=require({str(JS)!r});const d=JSON.parse(require('fs').readFileSync(0,'utf8')||'null');{code}"],
                       input=json.dumps(payload), capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


@needs_node
def test_a_point_is_found_in_its_cell_and_nowhere_in_the_gaps():
    f = lambda lat, lon: node("console.log(D.cellAt(d.u,d.lat,d.lon))", {"u": UNITS, "lat": lat, "lon": lon})
    assert f(40.0, 44.0) == 0 and f(40.12, 44.09) == 1 and f(40.2, 44.2) == 2
    assert f(40.0, 44.02) == 0 and f(40.0, 44.0251) == -1 and f(40.0, 44.03) == -1 and f(41.0, 44.0) == -1         # half a cell is 0.025 degrees


CTX = "{groupLabel:g=>g.toUpperCase(),methodCost:m=>({cheap:186,dear:900})[m]??null,methods:['cheap','dear','unknown'],marzOf:i=>['Syunik','Syunik','Tavush'][i]}"


@needs_node
def test_rows_carry_place_area_group_method_cost_and_survival():
    rows = node(f"console.log(JSON.stringify(D.portfolioRows(d.u,d.sel,{CTX})))", {"u": UNITS, "sel": [[0, 0, 0], [1, 1, 1]]})
    assert len(rows) == 2
    a, b = rows
    assert (a["latitude"], a["longitude"], a["marz"], a["area_ha"], a["species_group"], a["method"]) == (40.0, 44.0, "Syunik", 3200.0, "A", "cheap")
    assert a["cost_per_ha_usd"] == 186 and a["cost_usd"] == 186 * 3200 and a["viability_mean"] == 0.98 and a["viability_worst20pct"] == 0.96
    assert b["species_group"] == "B" and b["cost_usd"] == round(900 * 3190) and b["viability_mean"] == 0.95
    unknown = node(f"console.log(JSON.stringify(D.portfolioRows(d.u,d.sel,{CTX})))", {"u": UNITS, "sel": [[2, 0, 2]]})[0]
    assert unknown["cost_per_ha_usd"] is None and unknown["cost_usd"] is None                                  # an unknown cost stays unknown, never zero


@needs_node
def test_the_csv_parses_back_quotes_awkward_text_and_leaves_unknowns_empty():
    rows = node(f"const r=D.portfolioRows(d.u,d.sel,{CTX});r[0].marz='Odd, \"marz\"';console.log(JSON.stringify(D.csv(r)))", {"u": UNITS, "sel": [[0, 0, 0], [2, 0, 2]]})
    parsed = list(pycsv.DictReader(io.StringIO(rows)))
    assert rows.endswith("\r\n") and len(parsed) == 2 and "group_id" not in parsed[0]
    assert parsed[0]["marz"] == 'Odd, "marz"' and parsed[0]["area_ha"] == "3200" and parsed[1]["cost_usd"] == "" and parsed[1]["cost_per_ha_usd"] == ""
    assert node("console.log(JSON.stringify(D.csv([])))") == ""


@needs_node
def test_the_marz_table_sums_area_and_cells_and_orders_by_area():
    t = node(f"console.log(JSON.stringify(D.marzTable(D.portfolioRows(d.u,d.sel,{CTX}))))", {"u": UNITS, "sel": [[0, 0, 0], [1, 1, 1], [2, 0, 0]]})
    assert [m["marz"] for m in t] == ["Syunik", "Tavush"] and t[0]["cells"] == 2 and t[0]["area_ha"] == 3200 + 3190 and t[1]["area_ha"] == 3180
    assert t[0]["byGroup"] == {"a": 3200, "b": 3190}


# ---- the exported portfolio against the optimiser's own totals (skipped until a portfolio has been built)
manifest_path, units_path = DATA / "manifest.json", DATA / "aegis_units.json"
built = manifest_path.exists() and units_path.exists() and (json.loads(manifest_path.read_text()).get("aegis") or {}).get("units_file")


@pytest.mark.skipif(not built, reason="no portfolio built")
def test_every_treated_cell_is_eligible_valid_and_counted_once():
    U = json.loads(units_path.read_text())
    u, n = U["units"], len(U["units"]["lat"])
    assert all(len(u[k]) == n for k in ("lon", "area_ha", "eligible")) and all(len(v) == n for v in u["viability_mean"].values())
    for b in U["budgets"]:
        cells = [s[0] for s in b["selected"]]
        assert len(cells) == len(set(cells))                                                        # at most one option per cell
        assert all(u["eligible"][i] for i in cells)
        assert all(0 <= s[1] < len(u["group_ids"]) and 0 <= s[2] < len(U["intervention_ids"]) for s in b["selected"])


@pytest.mark.skipif(not built, reason="no portfolio built")
def test_the_exported_cells_add_up_to_the_totals_the_page_shows():
    M, U = json.loads(manifest_path.read_text())["aegis"], json.loads(units_path.read_text())
    u, cost = U["units"], {c["name"]: c["cost_per_ha"] for c in M["cost_table"]}
    assert [b["budget_usd"] for b in U["budgets"]] == [r["budget_usd"] for r in M["budget_sweep"]]
    for b, row in zip(U["budgets"], M["budget_sweep"]):
        assert len(b["selected"]) == row["n_units_planted"]
        area = sum(u["area_ha"][i] for i, _, _ in b["selected"])
        usd = sum(u["area_ha"][i] * cost[U["intervention_ids"][m]] for i, _, m in b["selected"])
        assert area == pytest.approx(row["area_ha"], rel=1e-3) and usd == pytest.approx(row["cost_usd"], rel=1e-3)
        assert usd <= row["budget_usd"] * (1 + 1e-9)                                                  # the plan never spends more than the budget
