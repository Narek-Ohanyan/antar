/* Decision page: the robust portfolio on a map ("where to plant what"), with a budget selector.

   Pure functions (no DOM), loaded by the browser as window.Decision and by node (tests/test_decision_page.py):
     cellAt(units, lat, lon)                   index of the grid cell that contains a point, or -1
     portfolioRows(units, selected, ctx)       one row per treated cell: place, area, group, method, cost, viability
     csv(rows)                                 RFC 4180 text of those rows
     marzTable(rows)                           treated area and cells per marz, with the area per species group
   The DOM part (mount) draws the cells on the same base map as the Place explorer. Nothing is computed here that the optimiser did not decide:
   which cells are treated, and with what, is read from data/aegis_units.json. */
(function (root) {
  "use strict";
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const GROUP_COLOURS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"];      // ColorBrewer Dark2: distinguishable with the common colour-vision differences

  function cellAt(u, lat, lon) {
    const hl = u.cell_deg.dlat / 2, ho = u.cell_deg.dlon / 2;
    for (let i = 0; i < u.lat.length; i++) if (Math.abs(lat - u.lat[i]) <= hl && Math.abs(lon - u.lon[i]) <= ho) return i;
    return -1;
  }

  /* selected: [[unit, group index, method index], ...]; ctx: { groupLabel(id), methodCost(name), methods: [names], marzOf(unit) } */
  function portfolioRows(u, selected, ctx) {
    return selected.map(([i, gi, mi]) => {
      const g = u.group_ids[gi], method = ctx.methods[mi], cost = ctx.methodCost(method);
      const treated = u.open_ha ? u.open_ha[i] : u.area_ha[i];                       // restoration acts on the cell's open land only
      return {
        latitude: u.lat[i], longitude: u.lon[i], marz: ctx.marzOf(i), area_ha: treated, cell_area_ha: u.area_ha[i], forest_share: u.forest_share ? u.forest_share[i] : null,
        species_group: ctx.groupLabel(g), group_id: g, method, cost_per_ha_usd: cost, cost_usd: cost == null ? null : Math.round(cost * treated),
        viability_mean: u.viability_mean[g][i], viability_worst20pct: u.viability_worst20[g][i],
      };
    });
  }

  function csv(rows) {
    if (!rows.length) return "";
    const cols = Object.keys(rows[0]).filter((k) => k !== "group_id");
    const q = (v) => (v == null ? "" : /[",\r\n]/.test(String(v)) ? '"' + String(v).replace(/"/g, '""') + '"' : String(v));
    return [cols.join(",")].concat(rows.map((r) => cols.map((c) => q(r[c])).join(","))).join("\r\n") + "\r\n";
  }

  function marzTable(rows) {
    const by = new Map();
    for (const r of rows) {
      const m = by.get(r.marz) || { marz: r.marz, cells: 0, area_ha: 0, byGroup: {} };
      m.cells += 1; m.area_ha += r.area_ha; m.byGroup[r.group_id] = (m.byGroup[r.group_id] || 0) + r.area_ha;
      by.set(r.marz, m);
    }
    return [...by.values()].sort((a, b) => b.area_ha - a.area_ha || a.marz.localeCompare(b.marz));
  }

  /* ---------------- the page (browser only) ---------------- */
  const ha = (v) => (v >= 1e5 ? Math.round(v / 1e3).toLocaleString("en") + " thousand ha" : Math.round(v).toLocaleString("en") + " ha");
  const money = (v) => "$" + (v >= 1e6 ? (v / 1e6).toFixed(v % 1e6 ? 1 : 0) + "M" : Math.round(v).toLocaleString("en"));

  function mount(host, { A, ae, units, sw, groupLabel, onCleanup }) {
    const G = A.grid, u = units.units, methods = units.intervention_ids, nGroups = u.group_ids.length;
    const costOf = Object.fromEntries((ae.cost_table || []).map((c) => [c.name, c.cost_per_ha]));
    const regionAt = (lat, lon) => {
      const col = Math.floor((Interp.mercX(lon) - G.x0) / G.px_m), row = Math.floor((G.y_top - Interp.mercY(lat)) / G.px_m);
      const id = col >= 0 && row >= 0 && col < A.w && row < A.h ? A.region[row * A.w + col] : 0;
      const r = id ? G.regions.find((x) => x.id === id) : null;
      return r ? (r.name === "Yerevan" ? "Yerevan (city)" : r.name) : "outside Armenia";
    };
    const marzOfCell = u.lat.map((la, i) => regionAt(la, u.lon[i]));
    const ctx = { groupLabel, methodCost: (m) => (costOf[m] != null ? costOf[m] : null), methods, marzOf: (i) => marzOfCell[i] };
    const bySel = new Map(units.budgets.map((b) => [b.budget_usd, b.selected]));
    const $ = (q) => host.querySelector(q);
    const short = (g) => groupLabel(g).replace(/\s*\(.*\)/, "");
    const cssv = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

    const bounds = G.bounds_latlon;
    const map = L.map($("#dmap"), { zoomControl: true, minZoom: 7, maxZoom: 11, zoomSnap: 0.25, attributionControl: false, scrollWheelZoom: false,
      maxBounds: [[bounds[0][0] - 0.25, bounds[0][1] - 0.4], [bounds[1][0] + 0.25, bounds[1][1] + 0.4]], maxBoundsViscosity: 1 });
    if (typeof state !== "undefined") state.mapObj = map;                 // the router removes it when the page changes
    map.fitBounds(bounds);
    map.createPane("base").style.zIndex = 210; map.getPane("base").style.pointerEvents = "none";
    map.createPane("cells").style.zIndex = 330; map.getPane("cells").style.pointerEvents = "none";
    L.imageOverlay(drawBase(A, "all").toDataURL("image/png"), bounds, { pane: "base", interactive: false, className: "map-img" }).addTo(map);
    L.geoJSON(A.borders, { interactive: false, style: (f) => (f.properties.kind === "country" ? { color: cssv("--map-outline"), weight: 1.6 } : { color: cssv("--map-line"), weight: 0.8, opacity: 0.8 }) }).addTo(map);

    const S = 4;                                                            // canvas pixels per map pixel: crisp cell edges
    const toPx = (lon, lat) => [((Interp.mercX(lon) - G.x0) / G.px_m) * S, ((G.y_top - Interp.mercY(lat)) / G.px_m) * S];
    let overlay = null, current = null, pinned = -1;
    function paint(budget) {
      const selected = bySel.get(budget) || [], chosen = new Map(selected.map(([i, gi, mi]) => [i, [gi, mi]]));
      const c = document.createElement("canvas"); c.width = A.w * S; c.height = A.h * S;
      const cx = c.getContext("2d"), dla = u.cell_deg.dlat / 2, dlo = u.cell_deg.dlon / 2;
      for (let i = 0; i < u.lat.length; i++) {
        const [x0, y0] = toPx(u.lon[i] - dlo, u.lat[i] + dla), [x1, y1] = toPx(u.lon[i] + dlo, u.lat[i] - dla), w = x1 - x0, h = y1 - y0;
        if (chosen.has(i)) { cx.fillStyle = GROUP_COLOURS[chosen.get(i)[0] % GROUP_COLOURS.length]; cx.globalAlpha = 0.86; cx.fillRect(x0, y0, w, h); cx.globalAlpha = 1; cx.strokeStyle = "rgba(255,255,255,.75)"; cx.lineWidth = 1.5; cx.strokeRect(x0 + 0.75, y0 + 0.75, w - 1.5, h - 1.5); }
        else if (u.eligible[i]) { cx.fillStyle = "rgba(40,40,40,.26)"; cx.fillRect(x0, y0, w, h); cx.strokeStyle = "rgba(255,255,255,.35)"; cx.lineWidth = 1; cx.strokeRect(x0 + 0.5, y0 + 0.5, w - 1, h - 1); }
        else { cx.strokeStyle = "rgba(0,0,0,.14)"; cx.lineWidth = 1; cx.strokeRect(x0 + 0.5, y0 + 0.5, w - 1, h - 1); }
      }
      if (overlay) map.removeLayer(overlay);
      overlay = L.imageOverlay(c.toDataURL("image/png"), bounds, { pane: "cells", interactive: false, className: "map-img" }).addTo(map);
      current = { budget, selected, chosen, rows: portfolioRows(u, selected, ctx) };
    }

    function side(row) {
      const rows = current.rows, tot = rows.reduce((a, r) => a + r.area_ha, 0);
      $("#dgroups").innerHTML = rows.length ? u.group_ids.map((g, gi) => {
        const a = rows.filter((r) => r.group_id === g).reduce((s, r) => s + r.area_ha, 0);
        return `<div class="dg"><span class="sw" style="background:${GROUP_COLOURS[gi % GROUP_COLOURS.length]}"></span><span class="dgn">${esc(short(g))}</span><span class="dgb"><i style="width:${tot ? (100 * a / tot).toFixed(1) : 0}%;background:${GROUP_COLOURS[gi % GROUP_COLOURS.length]}"></i></span><span class="dgv">${ha(a)}</span></div>`;
      }).join("") : '<p class="small muted">No cell is treated at this budget.</p>';
      const used = {}; rows.forEach((r) => { used[r.method] = (used[r.method] || 0) + r.area_ha; });
      $("#dmethods").innerHTML = Object.keys(used).length ? "<ul>" + Object.entries(used).sort((a, b) => b[1] - a[1]).map(([m, a]) => `<li>${esc(m.replace(/_/g, " "))}: ${ha(a)}${costOf[m] != null ? ` at $${costOf[m].toLocaleString("en")}/ha` : ""}</li>`).join("") + "</ul>" : "";
      const mt = marzTable(rows);
      $("#dmarz").innerHTML = mt.length ? `<table><thead><tr><th>Marz</th><th class="num">Cells</th><th class="num">Treated area</th></tr></thead><tbody>${mt.map((m) => `<tr><td>${esc(m.marz)}</td><td class="num">${m.cells}</td><td class="num">${ha(m.area_ha)}</td></tr>`).join("")}</tbody></table>` : "";
    }

    function render(idx) {
      const r = sw[idx];
      paint(r.budget_usd);
      $("#dbudget").textContent = money(r.budget_usd);
      $("#dslider").setAttribute("aria-valuetext", `${money(r.budget_usd)}: ${r.n_units_planted} cells, ${ha(r.area_ha)}`);
      $("#dstats").innerHTML = `<div><dt>Cells treated</dt><dd>${r.n_units_planted}<small>of ${ae.n_eligible_units} candidates</small></dd></div><div><dt>Open land treated</dt><dd>${ha(r.area_ha)}</dd></div>
        <div><dt>Cost</dt><dd>${money(r.cost_usd)}<small>of the ${money(r.budget_usd)} budget</small></dd></div><div><dt>Expected benefit</dt><dd>${Math.round(r.expected).toLocaleString("en")}<small>$/yr</small></dd></div>
        <div><dt>Worst 20% of scenarios</dt><dd>${Math.round(r.cvar).toLocaleString("en")}<small>$/yr</small></dd></div>`;
      side();
      $("#dcell").innerHTML = '<span class="muted">Point at a cell on the map to see what it holds.</span>';
    }

    function describe(i) {
      if (i < 0) return '<span class="muted">Point at a cell on the map to see what it holds.</span>';
      const ch = current.chosen.get(i), land = u.open_ha ? `, ${ha(u.open_ha[i])} of it open land (${Math.round(100 * u.forest_share[i])}% forest, ${Math.round(100 * u.woodland_share[i])}% woodland)` : "";
      const head = `<strong>${esc(marzOfCell[i])}</strong> · ${u.lat[i].toFixed(3)}°N, ${u.lon[i].toFixed(3)}°E · cell of ${ha(u.area_ha[i])}${land}`;
      const v = u.group_ids.map((g, gi) => `<li><span class="sw" style="background:${GROUP_COLOURS[gi % GROUP_COLOURS.length]}"></span>${esc(short(g))}: survival ${(100 * u.viability_mean[g][i]).toFixed(1)}% on average, ${(100 * u.viability_worst20[g][i]).toFixed(1)}% in the worst 20% of scenarios</li>`).join("");
      const verdict = !u.eligible[i] ? "Not a candidate: protected, mostly settlements, cropland or quarries, or too little open land (mostly forest, woodland or water)." : ch ? `<b>Treated:</b> ${esc(short(u.group_ids[ch[0]]))}, ${esc(methods[ch[1]].replace(/_/g, " "))}.` : "Eligible, not treated at this budget.";
      return `${head}<br>${verdict}<ul class="dv">${v}</ul>`;
    }
    const onMove = (e) => { if (pinned < 0) $("#dcell").innerHTML = describe(cellAt(u, e.latlng.lat, e.latlng.lng)); };
    map.on("mousemove", onMove);
    map.on("click", (e) => { const i = cellAt(u, e.latlng.lat, e.latlng.lng); pinned = i === pinned ? -1 : i; $("#dcell").innerHTML = describe(i); });
    map.on("mouseout", () => { if (pinned < 0) $("#dcell").innerHTML = describe(-1); });

    const slider = $("#dslider");
    slider.addEventListener("input", () => { pinned = -1; render(+slider.value); });
    $("#dcsv").addEventListener("click", () => {
      const url = URL.createObjectURL(new Blob([csv(current.rows)], { type: "text/csv;charset=utf-8" }));
      const a = document.createElement("a"); a.href = url; a.download = `antar_portfolio_${Math.round(current.budget / 1e6)}M.csv`; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    });
    render(+slider.value);
    $("#dlegend").innerHTML = u.group_ids.map((g, gi) => `<span><i style="background:${GROUP_COLOURS[gi % GROUP_COLOURS.length]}"></i>${esc(short(g))}</span>`).join("") +
      '<span><i class="elig"></i>eligible, not treated</span><span><i class="inel"></i>not eligible</span>';
    if (onCleanup) onCleanup(() => map.off());
  }

  const api = { cellAt, portfolioRows, csv, marzTable, GROUP_COLOURS, mount };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Decision = api;
})(typeof window !== "undefined" ? window : globalThis);
