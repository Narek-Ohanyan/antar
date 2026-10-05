/* ANTAR map page: Armenia only, no street basemap.

   Layers, bottom to top:
     base      the country in one colour, with mapped forest / woodland / water from the national Ecosystem Map
     surface   the chosen model quantity interpolated from the model nodes over the WHOLE country (Interp, IDW)
     hatch     optional hatching where forest stands today
     lines     marz borders and the national outline
     nodes     optional: the model nodes the surface is interpolated from

   The surface is an interpolation, not a model run at every pixel; the panel states the interpolation
   error (leave-one-out at the nodes) for whatever is on screen. Helpers (state, $, fmt, makeScale, ...) come
   from app.js and are used only at call time. */

const MAP_FILES = ["region", "forest", "woodland", "water", "forest_broadleaf", "forest_oak", "forest_pine", "forest_juniper"];
const COVER_GROUPS = { all: "Forest and woodland (all)", broadleaf: "Broadleaf forest (Fagus, other deciduous)", oak: "Oak forest (Quercus)", pine: "Pine forest (Pinus kochiana)", juniper: "Juniper woodland" };
const IDW_K = 8, IDW_POWER = 2;
/* decimals follow magnitude: 640.988 m is noise, 0.0123 probability is not */
const fa = (v) => (ok(v) ? (Math.abs(v) >= 100 ? Math.round(v).toLocaleString("en") : Math.abs(v) >= 10 ? v.toFixed(1) : Math.abs(v) >= 1 ? v.toFixed(2) : v.toFixed(4)) : "—");

function loadPng(url) {
  return new Promise((resolve, reject) => {
    const im = new Image();
    im.onload = () => {
      const c = document.createElement("canvas"); c.width = im.width; c.height = im.height;
      const cx = c.getContext("2d", { willReadFrequently: true }); cx.drawImage(im, 0, 0);
      const d = cx.getImageData(0, 0, im.width, im.height).data, out = new Uint8Array(im.width * im.height);
      for (let i = 0; i < out.length; i++) out[i] = d[i * 4];
      resolve(out);
    };
    im.onerror = () => reject(new Error("could not load " + url));
    im.src = url;
  });
}

async function loadMapAssets() {
  if (state.mapAssets) return state.mapAssets;
  const [grid, borders, ...rasters] = await Promise.all([
    fetch("assets/map/grid.json").then((r) => r.json()), fetch("assets/map/borders.geojson").then((r) => r.json()),
    ...MAP_FILES.map((f) => loadPng(`assets/map/${f}.png`)),
  ]);
  const A = { grid, borders, w: grid.width, h: grid.height };
  MAP_FILES.forEach((f, i) => (A[f] = rasters[i]));
  const idx = []; for (let i = 0; i < A.region.length; i++) if (A.region[i] > 0) idx.push(i);
  A.inIdx = Int32Array.from(idx); A.n = idx.length;
  A.rank = new Int32Array(A.region.length).fill(-1);
  A.tx = new Float64Array(A.n); A.ty = new Float64Array(A.n); A.area = new Float32Array(A.n); A.lat = new Float32Array(A.n);
  for (let k = 0; k < A.n; k++) {
    const f = A.inIdx[k], row = Math.floor(f / A.w), col = f % A.w;
    A.rank[f] = k;
    A.tx[k] = grid.x0 + (col + 0.5) * grid.px_m; A.ty[k] = grid.y_top - (row + 0.5) * grid.px_m;
    A.lat[k] = (Math.atan(Math.sinh(A.ty[k] / grid.earth_radius_m)) * 180) / Math.PI;
    A.area[k] = ((grid.px_m * Math.cos((A.lat[k] * Math.PI) / 180)) ** 2) / 1e6;       // ground km2 of this pixel
  }
  state.mapAssets = A;
  return A;
}

function interpFor(gid, A) {
  state.interp = state.interp || {};
  if (!state.interp[gid]) {
    const g = state.G[gid];
    const nx = Float64Array.from(g.lon, Interp.mercX), ny = Float64Array.from(g.lat, Interp.mercY);
    state.interp[gid] = { nx, ny, it: Interp.build(nx, ny, A.tx, A.ty, IDW_K, IDW_POWER) };
  }
  return state.interp[gid];
}

const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
function cssRgb(n) { const h = cssVar(n); return [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)); }
function lutFromScale(sc) {
  const lut = new Uint8Array(256 * 3);
  for (let i = 0; i < 256; i++) {
    const m = sc.fn(sc.lo + ((sc.hi - sc.lo) * i) / 255).match(/\d+/g);
    lut[i * 3] = +m[0]; lut[i * 3 + 1] = +m[1]; lut[i * 3 + 2] = +m[2];
  }
  return lut;
}
const toUrl = (canvas) => canvas.toDataURL("image/png");

/* ---- canvases ---- */
function drawBase(A, cover) {
  const c = document.createElement("canvas"); c.width = A.w; c.height = A.h;
  const cx = c.getContext("2d"), img = cx.createImageData(A.w, A.h), d = img.data;
  const land = cssRgb("--map-land"), forest = cssRgb("--map-forest"), wood = cssRgb("--map-woodland"), water = cssRgb("--map-water");
  const fr = cover === "all" ? A.forest : A["forest_" + cover];
  for (let k = 0; k < A.n; k++) {
    const f = A.inIdx[k];
    let rgb = land;
    const mix = (col, a) => (rgb = rgb.map((v, j) => v + (col[j] - v) * a));
    if (cover === "all") { mix(wood, A.woodland[f] / 255); mix(forest, A.forest[f] / 255); }
    else mix(forest, fr[f] / 255);
    mix(water, A.water[f] / 255);
    d[f * 4] = rgb[0]; d[f * 4 + 1] = rgb[1]; d[f * 4 + 2] = rgb[2]; d[f * 4 + 3] = 255;
  }
  cx.putImageData(img, 0, 0);
  return c;
}

function drawSurface(A, vals, sc, forestOnly) {
  const c = document.createElement("canvas"); c.width = A.w; c.height = A.h;
  const cx = c.getContext("2d"), img = cx.createImageData(A.w, A.h), d = img.data, lut = lutFromScale(sc);
  for (let k = 0; k < A.n; k++) {
    const v = vals[k], f = A.inIdx[k];
    if (!isFinite(v)) continue;
    if (forestOnly && (A.forest[f] + A.woodland[f]) / 255 < 0.25) continue;
    const t = Math.min(255, Math.max(0, Math.round(((v - sc.lo) / (sc.hi - sc.lo || 1)) * 255)));
    d[f * 4] = lut[t * 3]; d[f * 4 + 1] = lut[t * 3 + 1]; d[f * 4 + 2] = lut[t * 3 + 2]; d[f * 4 + 3] = 255;
  }
  cx.putImageData(img, 0, 0);
  return c;
}

function drawHatch(A) {
  const S = 3, c = document.createElement("canvas"); c.width = A.w * S; c.height = A.h * S;
  const cx = c.getContext("2d"), img = cx.createImageData(c.width, c.height), d = img.data;
  for (let y = 0; y < c.height; y++) {
    for (let x = 0; x < c.width; x++) {
      if (((x + y) % 7) > 1) continue;
      const f = Math.floor(y / S) * A.w + Math.floor(x / S);
      if (A.forest[f] / 255 < 0.5) continue;
      const o = (y * c.width + x) * 4; d[o] = 8; d[o + 1] = 50; d[o + 2] = 28; d[o + 3] = 150;
    }
  }
  cx.putImageData(img, 0, 0);
  return c;
}

/* ---- ensemble statistic across climate models, per node ---- */
function ensembleStat(members, mat, n, sel) {
  const idx = [];
  members.forEach((m, i) => { if (m.ssp === sel.ssp && m.horizon === +sel.hz && (sel.gcm === "ens" || m.gcm === sel.gcm)) idx.push(i); });
  return Array.from({ length: n }, (_, c) => {
    const vs = idx.map((i) => mat[i][c]).filter(ok);
    if (!vs.length) return null;
    if (sel.gcm !== "ens" || sel.stat === "mean") return mean(vs);
    const lo = Math.min(...vs), hi = Math.max(...vs);
    return sel.stat === "min" ? lo : sel.stat === "max" ? hi : hi - lo;
  });
}

/* Where does a quantity have scenario results, and on which grid? */
function scenarioSource(id) {
  for (const gid of preferDense(Object.keys(state.G))) {
    const sc = state.G[gid].scenarios || {};
    if (sc.fp_series && sc.fp_series[id]) return { gid, members: sc.members, mat: sc.fp_series[id] };
    if (sc.tl_series && sc.tl_series[id]) return { gid, members: sc.treeline_members, mat: sc.tl_series[id] };
  }
  return null;
}
function preferDense(ids) { return ids.slice().sort((a, b) => (b === "dense") - (a === "dense")); }

function mapQuantities() {
  const M = state.M, Q = {};
  for (const [id, meta] of Object.entries(M.layers)) if (meta.map !== false) Q[id] = { id, ...meta, src: scenarioSource(id) };
  const tl = scenarioSource("treeline_shift");
  if (tl) Q.treeline_shift = { id: "treeline_shift", label: "Treeline shift vs 2019 (scenario)", unit: "m", engine: "Treeline", grid: tl.gid, grid_label: M.grids[tl.gid].label, good: "high", scenarioOnly: true, src: tl };
  return Q;
}

function seriesFor(q, sel) {
  const M = state.M;
  if (sel.mode === "base" && !q.scenarioOnly) return { gid: q.grid, values: state.G[q.grid].layers[q.id], diverging: false, binary: q.binary };
  if (!q.src) return null;
  const g = state.G[q.src.gid], vals = ensembleStat(q.src.members, q.src.mat, g.n_cells, sel);
  const stops = q.good === "low" ? DIV_BAD_HIGH : q.good === "high" ? DIV_GOOD_HIGH : DIV_SHIFT;
  if (sel.stat === "spread" && sel.gcm === "ens") return { gid: q.src.gid, values: vals, diverging: false, spread: true, unit: q.unit };
  if (q.scenarioOnly) return { gid: q.src.gid, values: vals, diverging: true, stops, unit: "m" };
  if (sel.diff) {
    const base = g.layers[q.id];
    if (!base) return { gid: q.src.gid, values: vals, diverging: false, noDiff: true };
    return { gid: q.src.gid, values: vals.map((v, i) => (ok(v) && ok(base[i]) ? v - base[i] : null)), diverging: true, stops, change: true };
  }
  return { gid: q.src.gid, values: vals, diverging: false, binary: q.binary };
}

/* ---- the page ---- */
async function renderMap(p) {
  const M = state.M;
  view().innerHTML = `<div class="wrap"><div class="loading">Loading the map…</div></div>`;
  let A;
  try { A = await loadMapAssets(); } catch (e) { view().innerHTML = `<div class="wrap"><div class="callout bad"><strong>Map layers could not be loaded.</strong> ${esc(e.message)}. Run <code>python3 scripts/build_map_assets.py</code> and serve the <code>ui/</code> folder over http.</div></div>`; return; }
  if (!location.hash.startsWith("#/map")) return;                       // user navigated away while loading

  const Q = mapQuantities(), gcms = [...new Set((Object.values(state.G).find((g) => g.scenarios && g.scenarios.members) || { scenarios: { members: [] } }).scenarios.members.map((m) => m.gcm))].sort();
  const sel = {
    view: p.view === "cover" ? "cover" : "result", cover: COVER_GROUPS[p.cover] ? p.cover : "all",
    q: Q[p.q] ? p.q : "treeline_shift", mode: p.mode || "scen", ssp: p.ssp || "ssp585", hz: +(p.hz || 2100), gcm: p.gcm || "ens",
    stat: p.stat || "mean", diff: p.diff == null ? true : p.diff === "1",
    borders: p.borders !== "0", nodes: p.nodes === "1", hatch: p.hatch === "1", forestOnly: p.forestOnly === "1", opacity: +(p.op || 0.92),
  };
  if (!Q[sel.q]) sel.q = Object.keys(Q)[0];
  if (Q[sel.q].scenarioOnly || !Q[sel.q].src) sel.mode = Q[sel.q].scenarioOnly ? "scen" : sel.mode;
  const byEngine = {};
  Object.values(Q).forEach((q) => (byEngine[q.engine] = byEngine[q.engine] || []).push(q));
  const opts = ENGINE_ORDER.filter((e) => byEngine[e]).map((e) => `<optgroup label="${esc(ENGINE_TITLE[e] || e)}">${byEngine[e].map((q) => `<option value="${esc(q.id)}">${esc(q.label)}${q.src ? "" : " (2019 only)"}</option>`).join("")}</optgroup>`).join("");

  view().innerHTML = `<div class="wrap">${banner()}<div class="maplayout">
    <aside class="card side">
      <div class="seg" role="tablist"><button data-view="result" class="${sel.view === "result" ? "on" : ""}">Model result</button><button data-view="cover" class="${sel.view === "cover" ? "on" : ""}">Current forest cover</button></div>
      <div id="ctl-result"><label for="q">Quantity</label><select id="q">${opts}</select><div id="scenctl"></div></div>
      <div id="ctl-cover"><label for="cg">Cover shown</label><select id="cg">${Object.entries(COVER_GROUPS).map(([k, v]) => `<option value="${k}">${esc(v)}</option>`).join("")}</select></div>
      <div class="opts"><label class="chk"><input type="checkbox" id="o-borders"> Marz borders</label>
        <label class="chk res"><input type="checkbox" id="o-hatch"> Hatch where forest stands today</label>
        <label class="chk res"><input type="checkbox" id="o-forestonly"> Colour only mapped forest &amp; woodland</label>
        <label class="chk res"><input type="checkbox" id="o-nodes"> Show the model nodes</label>
        <label class="res" for="o-op">Surface opacity</label><input class="res" type="range" id="o-op" min="0.3" max="1" step="0.02"></div>
      <div id="qinfo" class="small" style="margin-top:12px"></div>
      <div class="legend" id="legend"></div>
      <div id="cellinfo" class="small" style="margin-top:12px"></div>
    </aside>
    <div><div id="map"></div><div id="regions" style="margin-top:12px"></div>
      <p class="small muted" style="margin-top:8px" id="mapnote"></p></div></div></div>`;

  const bounds = A.grid.bounds_latlon;
  const map = L.map("map", { zoomControl: true, scrollWheelZoom: true, minZoom: 7, maxZoom: 12, zoomSnap: 0.25, attributionControl: true,
    maxBounds: [[bounds[0][0] - 0.25, bounds[0][1] - 0.4], [bounds[1][0] + 0.25, bounds[1][1] + 0.4]], maxBoundsViscosity: 1 });
  state.mapObj = map;
  map.attributionControl.setPrefix(false);
  map.attributionControl.addAttribution("National outline and cover: Ecosystem Map of Armenia (CC BY 4.0) · marz borders: geoBoundaries (CC BY 2.5, approximate)");
  L.control.scale({ imperial: false }).addTo(map);
  map.fitBounds(bounds);
  ["base", "surface", "hatch"].forEach((n, i) => { map.createPane(n).style.zIndex = 210 + i * 10; map.getPane(n).style.pointerEvents = "none"; });
  let baseLayer, surfaceLayer, hatchLayer, bordersLayer, nodesLayer = L.layerGroup().addTo(map);
  const overlay = (canvas, pane, opacity) => L.imageOverlay(toUrl(canvas), bounds, { pane, opacity: opacity == null ? 1 : opacity, interactive: false, className: "map-img" }).addTo(map);
  bordersLayer = L.geoJSON(A.borders, { pane: "overlayPane", interactive: false,
    style: (f) => (f.properties.kind === "country" ? { color: cssVar("--map-outline"), weight: 1.8, opacity: 0.95 } : { color: cssVar("--map-line"), weight: 0.9, opacity: 0.8 }) });
  const hatchCanvas = drawHatch(A);
  let current = { vals: null, sc: null, q: null, s: null };

  function scenControls(q) {
    if (!q.src && !q.scenarioOnly) { $("#scenctl").innerHTML = `<p class="small muted" style="margin-top:8px">${q.engine === "TOPOHYDRO" || q.engine === "REFUGIUM" || q.engine === "XYLEM" ? "Scenario results for this quantity are not stored yet; they appear once the dense run completes." : "This quantity is computed for 2019 only."}</p>`; return; }
    $("#scenctl").innerHTML = `
      ${q.scenarioOnly ? "" : `<label for="mode">Period</label><select id="mode"><option value="base">2019 baseline</option><option value="scen">Future scenario</option></select>`}
      <div id="scen">
        <label for="ssp">Emissions path</label><select id="ssp">${Object.entries(SSP).map(([k, v]) => `<option value="${k}">${v}</option>`).join("")}</select>
        <label for="hz">Horizon</label><select id="hz">${HORIZONS.map((h) => `<option value="${h}">${h}</option>`).join("")}</select>
        <label for="gcm">Climate model</label><select id="gcm"><option value="ens">All 5 models</option>${gcms.map((g) => `<option value="${g}">${g}</option>`).join("")}</select>
        <div id="statrow"><label for="stat">Across the 5 models show</label><select id="stat"><option value="mean">Mean</option><option value="min">Lowest model</option><option value="max">Highest model</option><option value="spread">Disagreement (highest − lowest)</option></select></div>
        ${q.scenarioOnly ? "" : `<label class="chk" style="margin-top:12px"><input type="checkbox" id="diff"> Show change vs 2019</label>`}
      </div>`;
    if ($("#mode")) $("#mode").value = sel.mode;
    $("#ssp").value = sel.ssp; $("#hz").value = sel.hz; $("#gcm").value = sel.gcm; $("#stat").value = sel.stat;
    if ($("#diff")) $("#diff").checked = sel.diff;
    const sync = () => { $("#scen").style.display = q.scenarioOnly || sel.mode === "scen" ? "" : "none"; $("#statrow").style.display = sel.gcm === "ens" ? "" : "none"; };
    sync();
    ["mode", "ssp", "hz", "gcm", "stat"].forEach((id) => $("#" + id) && $("#" + id).addEventListener("change", (e) => { sel[id] = e.target.value; sync(); refresh(); }));
    if ($("#diff")) $("#diff").addEventListener("change", (e) => { sel.diff = e.target.checked; refresh(); });
  }

  function regionTable(vals, unit, isCover) {
    const NR = A.grid.regions.length + 1;
    const aS = new Float64Array(NR), vS = new Float64Array(NR), fS = new Float64Array(NR), wS = new Float64Array(NR), cvS = new Float64Array(NR), cS = new Float64Array(NR);
    const covF = isCover && sel.cover !== "all" ? A["forest_" + sel.cover] : null;
    for (let k = 0; k < A.n; k++) {
      const f = A.inIdx[k], r = A.region[f], a = A.area[k];
      const fo = (covF ? covF[f] : A.forest[f]) / 255, wo = covF ? 0 : A.woodland[f] / 255, cov = Math.min(1, fo + wo);
      const ok_ = vals && isFinite(vals[k]);
      for (const R of [r, 0]) {
        aS[R] += a; fS[R] += a * fo; wS[R] += a * wo;
        if (ok_) { vS[R] += a * vals[k]; cvS[R] += a * cov * vals[k]; cS[R] += a * cov; }   // cover-weighted mean = sum(a*cov*v) / sum(a*cov)
      }
    }
    const rows = [[0, "Armenia"]].concat(A.grid.regions.map((r) => [r.id, r.name === "Yerevan" ? "Yerevan (city)" : r.name]));
    if (isCover) {
      return `<div class="card tablewrap"><h3 style="margin-top:0">Cover by marz</h3><table><thead><tr><th>Marz</th><th class="num">Land (km²)</th><th class="num">${covF ? "Cover" : "Forest"} (km²)</th><th class="num">${covF ? "" : "Woodland (km²)"}</th><th class="num">Share of land</th></tr></thead><tbody>${rows.map(([id, nm]) => `<tr${id === 0 ? ' style="font-weight:600"' : ""}><td>${esc(nm)}</td><td class="num">${fmt(aS[id], 0)}</td><td class="num">${fmt(fS[id], 0)}</td><td class="num">${covF ? "" : fmt(wS[id], 0)}</td><td class="num">${pct((fS[id] + wS[id]) / aS[id], 1)}</td></tr>`).join("")}</tbody></table>
        <p class="small muted" style="margin:8px 0 0">Forest = classes 31–37 (closed forests and plantations). Woodland = classes 39, 41, 43, 44 (open plantation, subalpine, mixed and juniper woodland). Areas are from the 10 m Ecosystem Map averaged to 500 m cells.</p></div>`;
    }
    const u = esc(unit || "");
    return `<div class="card tablewrap"><h3 style="margin-top:0">By marz</h3><table><thead><tr><th>Marz</th><th class="num">Mean over all land</th><th class="num">Mean over forest &amp; woodland</th><th class="num">Forest &amp; woodland share</th></tr></thead><tbody>${rows.map(([id, nm]) => `<tr${id === 0 ? ' style="font-weight:600"' : ""}><td>${esc(nm)}</td><td class="num">${fa(vS[id] / aS[id])} ${u}</td><td class="num">${cS[id] > 0 ? fa(cvS[id] / cS[id]) + " " + u : "—"}</td><td class="num">${pct((fS[id] + wS[id]) / aS[id], 1)}</td></tr>`).join("")}</tbody></table>
      <p class="small muted" style="margin:8px 0 0">Area-weighted means of the interpolated surface. "Forest &amp; woodland" weights each pixel by the cover mapped there today.</p></div>`;
  }

  function setLayer(kind, canvas, opacity) {
    if (kind === "base") { if (baseLayer) map.removeLayer(baseLayer); baseLayer = overlay(canvas, "base"); }
    if (kind === "surface") { if (surfaceLayer) map.removeLayer(surfaceLayer); surfaceLayer = canvas ? overlay(canvas, "surface", opacity) : null; }
    if (kind === "hatch") { if (hatchLayer) map.removeLayer(hatchLayer); hatchLayer = canvas ? overlay(canvas, "hatch") : null; }
  }

  function refresh() {
    $("#ctl-result").style.display = sel.view === "result" ? "" : "none";
    $("#ctl-cover").style.display = sel.view === "cover" ? "" : "none";
    view().querySelectorAll(".res").forEach((e) => (e.style.display = sel.view === "result" ? "" : "none"));
    view().querySelectorAll(".seg button").forEach((b) => b.classList.toggle("on", b.dataset.view === sel.view));
    $("#o-borders").checked = sel.borders; $("#o-hatch").checked = sel.hatch; $("#o-forestonly").checked = sel.forestOnly; $("#o-nodes").checked = sel.nodes; $("#o-op").value = sel.opacity;
    if (bordersLayer) { if (sel.borders && !map.hasLayer(bordersLayer)) bordersLayer.addTo(map); if (!sel.borders && map.hasLayer(bordersLayer)) map.removeLayer(bordersLayer); }
    setLayer("base", drawBase(A, sel.view === "cover" ? sel.cover : "all"));
    nodesLayer.clearLayers();
    const params = { view: sel.view, cover: sel.cover, q: sel.q, mode: sel.mode, ssp: sel.ssp, hz: sel.hz, gcm: sel.gcm, stat: sel.stat, diff: sel.diff ? 1 : 0, borders: sel.borders ? 1 : 0, nodes: sel.nodes ? 1 : 0, hatch: sel.hatch ? 1 : 0, forestOnly: sel.forestOnly ? 1 : 0, op: sel.opacity };
    setParams(params);
    current = { vals: null, sc: null, q: null, s: null };

    if (sel.view === "cover") {
      setLayer("surface", null); setLayer("hatch", null);
      const cg = $("#cg"); cg.value = sel.cover;
      const grp = A.grid.cover_area_km2;
      $("#qinfo").innerHTML = `<div><strong>${esc(COVER_GROUPS[sel.cover])}</strong></div><div class="muted">Ecosystem Map of Armenia, 10 m, shown at 500 m</div>
        <div class="kv" style="margin-top:8px"><span>Forest (classes 31–37)</span><span>${fmt(grp.forest, 0)} km² · ${pct(grp.forest / A.grid.country_area_km2, 1)}</span><span>Woodland (39, 41, 43, 44)</span><span>${fmt(grp.woodland, 0)} km² · ${pct(grp.woodland / A.grid.country_area_km2, 1)}</span><span>Water bodies</span><span>${fmt(grp.water, 0)} km²</span></div>`;
      $("#legend").innerHTML = `<div><span class="swatch"><i style="background:${cssVar("--map-forest")}"></i>forest</span><span class="swatch"><i style="background:${cssVar("--map-woodland")}"></i>woodland</span><span class="swatch"><i style="background:${cssVar("--map-land")}"></i>everything else</span><span class="swatch"><i style="background:${cssVar("--map-water")}"></i>water</span></div>`;
      $("#regions").innerHTML = regionTable(null, "", true);
      $("#mapnote").textContent = "Everything that is not forest, woodland or water is drawn in a single colour. Hover the map for the cover under the pointer.";
      return;
    }

    const q = Q[sel.q], s = seriesFor(q, sel);
    const meta = [`<div><strong>${esc(q.label)}</strong></div>`, `<div class="muted">${esc(q.engine)} · ${gridChip(s ? s.gid : q.grid)}</div>`];
    if (q.placeholder) meta.push(`<div class="callout" style="margin:8px 0;padding:7px 10px">Placeholder input: ${esc(q.placeholder)}</div>`);
    if (!s || !s.values || !s.values.some(ok)) {
      setLayer("surface", null); setLayer("hatch", null);
      $("#qinfo").innerHTML = meta.join("") + `<p class="muted">No data for this selection.</p>`; $("#legend").innerHTML = ""; $("#regions").innerHTML = ""; return;
    }
    if (s.noDiff) meta.push(`<div class="callout" style="margin:8px 0;padding:7px 10px">No 2019 baseline exists on this grid, so change cannot be shown.</div>`);
    const g = state.G[s.gid], ip = interpFor(s.gid, A);
    const vals = Interp.apply(ip.it, s.values);
    const unit = s.unit != null ? s.unit : q.unit || "";
    let sc;
    if (s.binary) { const st = ["#c0583a", "#e8d9a0", "#2f8f5b"]; sc = { fn: (x) => ramp(st, x), lo: 0, hi: 1, stops: st }; }
    else { const sample = []; for (let k = 0; k < vals.length; k += 7) if (isFinite(vals[k])) sample.push(vals[k]); sc = makeScale(sample, { diverging: s.diverging, stops: s.stops, unit: s.spread ? "" : unit }); }
    setLayer("surface", drawSurface(A, vals, sc, sel.forestOnly), sel.opacity);
    setLayer("hatch", sel.hatch ? hatchCanvas : null);
    current = { vals, sc, q, s, unit };
    const loo = Interp.leaveOneOut(ip.nx, ip.ny, s.values, IDW_K, IDW_POWER);
    const present = s.values.filter(ok), sd = Math.sqrt(present.reduce((a, v) => a + (v - mean(present)) ** 2, 0) / present.length);
    const skill = loo.r2 >= 0.7 ? ["good", "good"] : loo.r2 >= 0.3 ? ["warn", "moderate"] : ["bad", "weak"];
    const label = s.change ? "change vs 2019" : s.spread ? "disagreement between the 5 models" : "";
    meta.push(`<div class="kv" style="margin-top:8px"><span>Model nodes</span><span>${present.length} of ${g.n_cells}</span><span>Node range</span><span>${fmt(Math.min(...present), 4)} – ${fmt(Math.max(...present), 4)} ${esc(s.binary ? "" : unit)}</span><span>Node mean</span><span>${fmt(mean(present), 3)}</span></div>
      <div class="callout ${skill[0] === "good" ? "info" : skill[0] === "bad" ? "bad" : ""}" style="margin:10px 0 0;padding:8px 10px"><strong>Interpolation check:</strong> predicting each node from the others gives R² = ${fmt(loo.r2, 2)}, RMSE ${fa(loo.rmse)} ${esc(s.binary ? "" : unit)} (spread of the nodes: sd ${fa(sd)}). ${chip(skill[0], skill[1])} ${skill[1] === "weak" ? "Between nodes the colours are poorly constrained; read the pattern as indicative only." : skill[1] === "moderate" ? "The surface follows the nodes only partly." : ""}</div>`);
    $("#qinfo").innerHTML = meta.join("");
    $("#legend").innerHTML = (label ? `<div class="small" style="margin-bottom:4px"><strong>${esc(label)}</strong></div>` : "") + (s.binary
      ? `<div class="bar" style="background:linear-gradient(90deg,${sc.stops.join(",")})"></div><div class="ends"><span>none of nearby nodes</span><span>share robust</span><span>all</span></div>`
      : `<div class="bar" style="background:linear-gradient(90deg,${sc.stops.join(",")})"></div><div class="ends"><span>${fa(sc.lo)}</span><span>${esc(unit)}${s.diverging ? " (centre = 0)" : ""}</span><span>${fa(sc.hi)}</span></div>${sc.widened ? `<p class="small muted" style="margin:4px 0 0"><strong>Scale widened</strong>: the values vary by less than ${fmt(MIN_SPAN[unit] || 0, 3)} ${esc(unit)}, so the colours span at least that much instead of exaggerating a negligible difference.</p>` : ""}`);
    if (sel.nodes) {
      const rad = g.n_cells < 200 ? 5 : 2.5;
      for (let i = 0; i < g.n_cells; i++) L.circleMarker([g.lat[i], g.lon[i]], { radius: rad, weight: 1, color: "#101a14", fillColor: ok(s.values[i]) ? sc.fn(s.values[i]) : "#bbb", fillOpacity: 1, pane: "markerPane" }).bindTooltip(`${g.lat[i].toFixed(3)}°N ${g.lon[i].toFixed(3)}°E · ${fmt(g.elev[i], 0)} m: ${fmt(s.values[i], 4)} ${unit}`).addTo(nodesLayer);
    }
    $("#regions").innerHTML = regionTable(vals, s.binary ? "" : unit, false);
    $("#mapnote").textContent = `Surface = inverse-distance interpolation (${IDW_K} nearest of ${g.n_cells} model nodes, power ${IDW_POWER}) clipped to Armenia; it is not a model run at every pixel. Forest cover and borders are the real national map.`;
  }

  /* pointer read-out */
  function readout(ll) {
    const pt = L.CRS.EPSG3857.project(ll), col = Math.floor((pt.x - A.grid.x0) / A.grid.px_m), row = Math.floor((A.grid.y_top - pt.y) / A.grid.px_m);
    if (col < 0 || row < 0 || col >= A.w || row >= A.h) return null;
    const f = row * A.w + col, k = A.rank[f];
    if (k < 0) return null;
    return { f, k, region: A.grid.regions.find((r) => r.id === A.region[f]), forest: A.forest[f] / 255, wood: A.woodland[f] / 255, water: A.water[f] / 255 };
  }
  function describe(ll, info) {
    const nm = info.region ? (info.region.name === "Yerevan" ? "Yerevan (city)" : info.region.name + " marz") : "";
    let h = `<strong>${esc(nm)}</strong> · ${ll.lat.toFixed(3)}°N, ${ll.lng.toFixed(3)}°E`;
    if (sel.view === "result" && current.vals) {
      const v = current.vals[info.k];
      h += `<br>${current.s.change ? "Change in " : current.s.spread ? "Model disagreement, " : ""}${esc(current.q.label)}: <strong>${fa(v)}${current.s.binary ? "" : " " + esc(current.unit || "")}</strong>`;
      const g = state.G[current.s.gid]; let best = -1, bd = Infinity;
      for (let i = 0; i < g.n_cells; i++) { const d = Math.hypot(Interp.mercX(g.lon[i]) - Interp.mercX(ll.lng), Interp.mercY(g.lat[i]) - Interp.mercY(ll.lat)); if (d < bd) { bd = d; best = i; } }
      h += `<br>Nearest model node: ${fmt(bd * Math.cos((ll.lat * Math.PI) / 180) / 1000, 1)} km away`;
      info.node = { gid: current.s.gid, i: best };
    }
    h += `<br>Mapped today: forest ${pct(info.forest, 0)} · woodland ${pct(info.wood, 0)}${info.water > 0.05 ? " · water " + pct(info.water, 0) : ""}`;
    return h;
  }
  let lastMove = 0;
  map.on("mousemove", (e) => {
    const now = Date.now(); if (now - lastMove < 40) return; lastMove = now;
    const info = readout(e.latlng); $("#cellinfo").innerHTML = info ? `<div class="card" style="padding:9px 11px">${describe(e.latlng, info)}</div>` : "";
  });
  map.on("click", (e) => {
    const info = readout(e.latlng); if (!info) return;
    const h = describe(e.latlng, info) + (info.node ? `<br><a href="#/site?grid=${info.node.gid}&i=${info.node.i}">Open nearest node in site explorer →</a>` : "");
    L.popup().setLatLng(e.latlng).setContent(h).openOn(map);
  });

  /* controls */
  view().querySelectorAll(".seg button").forEach((b) => b.addEventListener("click", () => { sel.view = b.dataset.view; refresh(); }));
  $("#q").value = sel.q;
  $("#q").addEventListener("change", (e) => { sel.q = e.target.value; if (Q[sel.q].scenarioOnly) sel.mode = "scen"; scenControls(Q[sel.q]); refresh(); });
  $("#cg").addEventListener("change", (e) => { sel.cover = e.target.value; refresh(); });
  $("#o-borders").addEventListener("change", (e) => { sel.borders = e.target.checked; refresh(); });
  $("#o-hatch").addEventListener("change", (e) => { sel.hatch = e.target.checked; refresh(); });
  $("#o-forestonly").addEventListener("change", (e) => { sel.forestOnly = e.target.checked; refresh(); });
  $("#o-nodes").addEventListener("change", (e) => { sel.nodes = e.target.checked; refresh(); });
  $("#o-op").addEventListener("input", (e) => { sel.opacity = +e.target.value; if (surfaceLayer) surfaceLayer.setOpacity(sel.opacity); });
  const mo = new MutationObserver(() => { if (sel.view) refresh(); });          // theme switch: redraw the base in the new palette
  mo.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  state.mapCleanup = () => mo.disconnect();
  scenControls(Q[sel.q]); refresh();
}
