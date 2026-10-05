/* Place explorer: pick any place in Armenia (click the map) or a whole marz, and see every quantity for it,
   today and under every emissions path. No sample points are drawn anywhere: a place's value is the same
   interpolation the map uses (Interp), and a marz's value is its area-weighted mean over the same surface.
   Uses map.js helpers (loadMapAssets, interpFor, drawBase, scenarioSource, mapQuantities, fa) at call time. */

function xhairIcon() { return L.divIcon({ className: "xhair", html: "<i></i>", iconSize: [26, 26], iconAnchor: [13, 13] }); }

async function renderSite(p) {
  const M = state.M;
  view().innerHTML = `<div class="wrap"><div class="loading">Loading…</div></div>`;
  let A;
  try { A = await loadMapAssets(); } catch (e) { view().innerHTML = `<div class="wrap"><div class="callout bad"><strong>Map layers could not be loaded.</strong> ${esc(e.message)}</div></div>`; return; }
  if (!location.hash.startsWith("#/site")) return;

  const G = A.grid, groups = Object.keys(M.groups), Q = mapQuantities(), plan = M.map_plan || [];
  const sel = { group: M.groups[p.group] ? p.group : groups[0], marz: +p.marz || 0, lat: isFinite(+p.lat) && p.lat ? +p.lat : 40.74, lon: isFinite(+p.lon) && p.lon ? +p.lon : 44.86 };
  const regionAt = (lat, lon) => {
    const col = Math.floor((Interp.mercX(lon) - G.x0) / G.px_m), row = Math.floor((G.y_top - Interp.mercY(lat)) / G.px_m);
    if (col < 0 || row < 0 || col >= A.w || row >= A.h) return null;
    const id = A.region[row * A.w + col];
    return id ? G.regions.find((r) => r.id === id) : null;
  };
  if (!sel.marz && !regionAt(sel.lat, sel.lon)) { sel.lat = 40.74; sel.lon = 44.86; }
  const rname = (r) => (r ? (r.name === "Yerevan" ? "Yerevan (city)" : r.name + " marz") : "outside Armenia");

  view().innerHTML = `<div class="wrap"><h1>Place explorer</h1>${banner()}
    <p class="lead muted" style="max-width:820px">Click anywhere in Armenia, or choose a whole marz, to see every quantity for that place today and under each emissions path (lines are the mean of the five climate models; shaded bands span the lowest to the highest model). Values come from the same interpolated surface as the map.</p>
    <div class="grid" style="grid-template-columns:minmax(280px,420px) 1fr;align-items:start">
      <div class="card"><div id="pmap" style="height:380px;border-radius:8px;border:1px solid var(--border);background:var(--surface-2)"></div>
        <label for="marz">Show</label><select id="marz"><option value="0">A specific place (click the map)</option>${G.regions.map((r) => `<option value="${r.id}">${esc(r.name === "Yerevan" ? "Yerevan (city)" : r.name + " marz")} — average</option>`).join("")}</select>
        <label for="grp">Species group</label><select id="grp">${groups.map((k) => `<option value="${k}">${esc(M.groups[k].label)}</option>`).join("")}</select>
        <div id="where" class="small" style="margin-top:12px"></div></div>
      <div class="card tablewrap"><h3 style="margin-top:0">Everything computed here, 2019</h3><div id="ptable"></div></div></div>
    <h2>Projections by emissions path</h2><div id="pleg" class="small" style="margin-bottom:8px"></div><div id="pcharts" class="grid cols-3"></div><div id="ppending"></div></div>`;

  const bounds = G.bounds_latlon;
  const map = L.map("pmap", { zoomControl: true, minZoom: 6, maxZoom: 11, zoomSnap: 0.25, attributionControl: false, scrollWheelZoom: false,
    maxBounds: [[bounds[0][0] - 0.25, bounds[0][1] - 0.4], [bounds[1][0] + 0.25, bounds[1][1] + 0.4]], maxBoundsViscosity: 1 });
  state.mapObj = map;
  map.fitBounds(bounds);
  map.createPane("base").style.zIndex = 210; map.getPane("base").style.pointerEvents = "none";
  L.imageOverlay(drawBase(A, "all").toDataURL("image/png"), bounds, { pane: "base", interactive: false, className: "map-img" }).addTo(map);
  L.geoJSON(A.borders, { interactive: false, style: (f) => (f.properties.kind === "country" ? { color: cssVar("--map-outline"), weight: 1.6 } : { color: cssVar("--map-line"), weight: 0.8, opacity: 0.8 }) }).addTo(map);
  let hi = null, pin = null;

  const cache = {};
  function evaluator(gid) {
    const key = gid + "|" + sel.marz + "|" + sel.lat + "|" + sel.lon;
    if (cache[key]) return cache[key];
    const ip = interpFor(gid, A), g = state.G[gid];
    let ev;
    if (!sel.marz) {
      const it = Interp.build(ip.nx, ip.ny, [Interp.mercX(sel.lon)], [Interp.mercY(sel.lat)], IDW_K, IDW_POWER);
      ev = (row) => Interp.apply(it, row)[0];
    } else {
      const { it } = ip, W = new Float64Array(g.n_cells);
      let aSum = 0;
      for (let k = 0; k < A.n; k++) {
        if (A.region[A.inIdx[k]] !== sel.marz) continue;
        let sw = 0; for (let j = 0; j < it.k; j++) sw += it.w[k * it.k + j];
        for (let j = 0; j < it.k; j++) W[it.idx[k * it.k + j]] += (A.area[k] * it.w[k * it.k + j]) / sw;
        aSum += A.area[k];
      }
      for (let j = 0; j < W.length; j++) W[j] /= aSum;
      ev = (row) => {
        if (row.every(ok)) { let v = 0; for (let j = 0; j < W.length; j++) v += W[j] * row[j]; return v; }
        const vals = Interp.apply(it, row); let s = 0, a = 0;                       // missing node values: exact pixel mean instead
        for (let k = 0; k < A.n; k++) if (A.region[A.inIdx[k]] === sel.marz && isFinite(vals[k])) { s += A.area[k] * vals[k]; a += A.area[k]; }
        return a ? s / a : null;
      };
    }
    return (cache[key] = ev);
  }

  const relevant = (id) => { const pl = plan.find((x) => x.id === id); return !pl || !pl.group || pl.group === sel.group; };
  function update() {
    $("#marz").value = sel.marz; $("#grp").value = sel.group;
    setParams(sel.marz ? { marz: sel.marz, group: sel.group } : { lat: sel.lat.toFixed(4), lon: sel.lon.toFixed(4), group: sel.group });
    const r = sel.marz ? G.regions.find((x) => x.id === sel.marz) : regionAt(sel.lat, sel.lon);
    $("#where").innerHTML = sel.marz ? `<strong>${esc(rname(r))}</strong> — area-weighted average over ${fmt(r.area_km2, 0)} km²` : `<strong>${esc(rname(r))}</strong><br>${sel.lat.toFixed(3)}°N, ${sel.lon.toFixed(3)}°E`;
    if (hi) map.removeLayer(hi); if (pin) map.removeLayer(pin);
    hi = sel.marz ? L.geoJSON(A.borders.features.filter((f) => f.properties.kind === "marz" && f.properties.id === sel.marz), { interactive: false, style: { color: cssVar("--accent"), weight: 2.4, fillColor: cssVar("--accent"), fillOpacity: 0.22 } }).addTo(map) : null;
    pin = sel.marz ? null : L.marker([sel.lat, sel.lon], { icon: xhairIcon(), interactive: false }).addTo(map);

    // 2019 table
    const byEngine = {};
    for (const [id, meta] of Object.entries(M.layers)) {
      if (meta.engine === "Ecosystem map" || !relevant(id)) continue;
      const row = state.G[meta.grid].layers[id], v = row ? evaluator(meta.grid)(row) : null;
      if (ok(v)) (byEngine[meta.engine] = byEngine[meta.engine] || []).push([meta, v]);
    }
    $("#ptable").innerHTML = `<table><tbody>${ENGINE_ORDER.filter((e) => byEngine[e]).map((e) => `<tr><th colspan="2">${esc(ENGINE_TITLE[e] || e)}</th></tr>` + byEngine[e].map(([m, v]) => `<tr><td>${esc(m.label)}</td><td class="num">${fa(v)} ${esc(m.binary ? "(share)" : m.unit || "")}</td></tr>`).join("")).join("")}</tbody></table>`;

    // scenario charts, every quantity that has scenario results
    const cards = [], pending = [];
    for (const pl of plan) {
      if (!relevant(pl.id)) continue;
      const src = scenarioSource(pl.id);
      if (!src) { pending.push(pl); continue; }
      const ev = evaluator(src.gid), vals = src.mat.map((row) => ev(row));
      const baseRow = pl.id === "treeline_shift" ? null : state.G[src.gid].layers[pl.id];
      const base = pl.id === "treeline_shift" ? 0 : baseRow ? ev(baseRow) : null;
      const series = Object.keys(SSP).map((ssp) => ({ name: SSP[ssp], color: SSP_COLOR[ssp],
        points: [[2019, base, null, null]].concat(HORIZONS.map((h) => {
          const vs = src.members.map((m, i) => (m.ssp === ssp && m.horizon === h ? vals[i] : null)).filter(ok);
          return vs.length ? [h, mean(vs), Math.min(...vs), Math.max(...vs)] : [h, null, null, null];
        })) }));
      const meta = M.layers[pl.id] || {}, unit = pl.id === "treeline_shift" ? "m" : meta.binary ? "share of models" : meta.unit || "";
      cards.push(`<div class="card"><h3 style="font-size:.98rem">${esc(pl.label)}</h3>${Charts.line({ series, xticks: [2019, 2050, 2080, 2100], ylabel: unit, yfmt: fa, xfmt: String, width: 440, height: 270, legendOn: false, ...(meta.binary || /^(pviab|robust)_/.test(pl.id) ? { ymin: 0, ymax: 1 } : {}) })}</div>`);
    }
    $("#pcharts").innerHTML = cards.join("") || `<p class="muted">No scenario results exist for this selection.</p>`;
    $("#pleg").innerHTML = Object.keys(SSP).map((s) => `<span class="swatch" style="margin-right:14px"><i style="background:${cssVar("--c" + (s === "ssp126" ? 1 : s === "ssp370" ? 2 : 3))}"></i>${SSP[s]}</span>`).join("") + `<span class="muted">the 2019 point is the baseline run</span>`;
    $("#ppending").innerHTML = pending.length ? `<div class="callout" style="margin-top:14px"><strong>Not computed yet (${pending.length}):</strong> ${pending.map((x) => esc(x.label)).join("; ")}. ${[...new Set(pending.map((x) => x.scenario_from))].map((f) => `<code>${esc(f)}</code>`).join(" and ")} produce them in the next dense run.</div>` : "";
  }

  map.on("click", (e) => { if (!regionAt(e.latlng.lat, e.latlng.lng)) return; sel.marz = 0; sel.lat = e.latlng.lat; sel.lon = e.latlng.lng; update(); });
  $("#marz").addEventListener("change", (e) => { sel.marz = +e.target.value; update(); });
  $("#grp").addEventListener("change", (e) => { sel.group = e.target.value; update(); });
  update();
}
