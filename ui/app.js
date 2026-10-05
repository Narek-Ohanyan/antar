/* ANTAR explorer -- static single-page app. Reads only ui/data/*.json (built by ui/build_data.py
   from the real fitted outputs). Nothing is computed here that is not already a fitted result,
   apart from simple means / ranges across ensemble members for display. */
const { esc } = Charts;
const state = { M: null, G: {}, meth: null, map: null, mapObj: null, markers: [] };
const $ = (s, el = document) => el.querySelector(s);
const view = () => $("#view");

/* ---------- helpers ---------- */
const SSP = { ssp126: "SSP1-2.6", ssp370: "SSP3-7.0", ssp585: "SSP5-8.5" };
const SSP_COLOR = { ssp126: "var(--c1)", ssp370: "var(--c2)", ssp585: "var(--c3)" };
const HORIZONS = [2050, 2080, 2100];
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const ok = (v) => v != null && isFinite(v);
function fmt(v, d = 3) { return ok(v) ? (Math.abs(v) >= 1000 ? Math.round(v).toLocaleString("en") : +Number(v).toFixed(d) + "") : "—"; }
const pct = (v, d = 1) => (ok(v) ? (100 * v).toFixed(d) + "%" : "—");
const chip = (kind, text) => `<span class="chip ${kind}">${esc(text)}</span>`;
const gridChip = (gid) => chip(gid === "dense" ? "good" : "neutral", gid === "dense" ? "1044-pt dense grid" : "80-pt validation grid");
const groupLabel = (g) => (state.M.groups[g] ? state.M.groups[g].label : g);
function pageParams() { const q = location.hash.split("?")[1] || ""; return Object.fromEntries(new URLSearchParams(q)); }
function setParams(p) { const base = location.hash.split("?")[0]; history.replaceState(null, "", base + "?" + new URLSearchParams(p).toString()); }
function scenarioGridId() { return Object.keys(state.G).filter((g) => state.G[g].scenarios && state.G[g].scenarios.members).sort().reverse()[0]; }
function treelineGridId() { return Object.keys(state.G).filter((g) => state.G[g].scenarios && state.G[g].scenarios.treeline_members).sort().reverse()[0]; }

/* ---------- colour scales ---------- */
const SEQ = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"];
const DIV_GOOD_HIGH = ["#b2182b", "#ef8a62", "#fddbc7", "#f2f2f2", "#d1e5f0", "#67a9cf", "#2166ac"];
const DIV_BAD_HIGH = DIV_GOOD_HIGH.slice().reverse();
const DIV_SHIFT = ["#8c510a", "#d8b365", "#f6e8c3", "#f2f2f2", "#c7eae5", "#5ab4ac", "#01665e"];
function hex2rgb(h) { return [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)); }
function ramp(stops, t) {
  t = Math.min(1, Math.max(0, t)); const x = t * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(x)), f = x - i;
  const a = hex2rgb(stops[i]), b = hex2rgb(stops[i + 1]);
  return `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * f)).join(",")})`;
}
function quantile(sorted, q) { const i = (sorted.length - 1) * q, lo = Math.floor(i), hi = Math.ceil(i); return sorted[lo] + (sorted[hi] - sorted[lo]) * (i - lo); }
/* Smallest colour span that is allowed to look like "a difference". Without it, a percentile-stretched
   scale paints a 0.9999-vs-1.0 gap as dark purple against bright yellow -- exaggerating variation
   that is physically negligible. */
const MIN_SPAN = { probability: 0.1, "probability/yr": 0.02, fraction: 0.1, m: 60, "°C": 1.0, mm: 25, days: 6, "°C·d": 100, MPa: 0.2, "": 0.1 };
function makeScale(values, { diverging, stops, binary, unit }) {
  const v = values.filter(ok).sort((a, b) => a - b);
  if (!v.length) return { fn: () => "#999", lo: 0, hi: 1, stops: SEQ };
  if (binary) return { fn: (x) => (x >= 0.5 ? "#2f8f5b" : "#c0583a"), binary: true, lo: 0, hi: 1 };
  let lo = quantile(v, 0.02), hi = quantile(v, 0.98);
  if (lo === hi) { lo = v[0]; hi = v[v.length - 1]; }
  const minSpan = MIN_SPAN[unit] != null ? MIN_SPAN[unit] : 0;
  let widened = false;
  if (diverging) {
    let m = Math.max(Math.abs(lo), Math.abs(hi)) || 1e-9;
    if (m < minSpan / 2) { m = minSpan / 2; widened = true; }
    lo = -m; hi = m;
  } else if (hi - lo < minSpan) {
    const c = (hi + lo) / 2; lo = c - minSpan / 2; hi = c + minSpan / 2; widened = true;
  }
  if (lo === hi) { lo -= 1; hi += 1; }
  const st = stops || (diverging ? DIV_SHIFT : SEQ);
  return { fn: (x) => ramp(st, (x - lo) / (hi - lo)), lo, hi, stops: st, min: v[0], max: v[v.length - 1], widened };
}

/* ---------- quantities for the map ---------- */
function scenarioKind(id) {
  if (id === "treeline_shift") return { kind: "treeline" };
  if (id.startsWith("viab_")) return { kind: "viab", group: id.slice(5) };
  if (id.startsWith("hmech_") && !id.startsWith("hmech_sd_")) return { kind: "hmech", group: id.slice(6) };
  return null;
}
function buildQuantities() {
  const M = state.M, Q = {};
  for (const [id, meta] of Object.entries(M.layers)) Q[id] = { id, ...meta, scen: scenarioKind(id) };
  const tg = treelineGridId();
  if (tg) Q.treeline_shift = { id: "treeline_shift", label: "Treeline shift vs 2019 (scenario)", unit: "m", engine: "Treeline", grid: tg, grid_label: M.grids[tg].label, scen: { kind: "treeline" }, scenarioOnly: true };
  return Q;
}
const ENGINE_ORDER = ["TOPOHYDRO", "XYLEM", "REFUGIUM", "Treeline", "Ecosystem map"];
const ENGINE_TITLE = { TOPOHYDRO: "Climate & water balance (TOPOHYDRO)", XYLEM: "Hydraulic-failure hazard (XYLEM)", REFUGIUM: "Viability & refugia (REFUGIUM)", Treeline: "Treeline", "Ecosystem map": "Observed land cover (Ecosystem Map)" };

function scenarioMembers(sg, kind) { return kind === "treeline" ? sg.scenarios.treeline_members : sg.scenarios.members; }
function scenarioMatrix(sg, scen) {
  if (scen.kind === "treeline") return sg.scenarios.treeline_shift;
  return scen.kind === "viab" ? sg.scenarios.viab[scen.group] : sg.scenarios.hmech[scen.group];
}
function ensembleAt(sg, scen, sel) {
  const members = scenarioMembers(sg, scen.kind), mat = scenarioMatrix(sg, scen);
  if (!mat) return null;
  const idx = [];
  members.forEach((m, i) => { if (m.ssp === sel.ssp && m.horizon === +sel.hz && (sel.gcm === "ens" || m.gcm === sel.gcm)) idx.push(i); });
  return Array.from({ length: sg.n_cells }, (_, c) => { const vs = idx.map((i) => mat[i][c]).filter(ok); return vs.length ? mean(vs) : null; });
}
function getSeries(q, sel) {
  const M = state.M;
  if (sel.mode === "base" && !q.scenarioOnly) {
    const g = state.G[q.grid];
    return { grid: g, gid: q.grid, values: g.layers[q.id], diverging: false, binary: M.layers[q.id] && M.layers[q.id].binary };
  }
  const sgid = q.scen.kind === "treeline" ? treelineGridId() : scenarioGridId();
  if (!sgid) return null;
  const sg = state.G[sgid];
  const vals = ensembleAt(sg, q.scen, sel);
  if (!vals) return null;
  if (q.scen.kind === "treeline") return { grid: sg, gid: sgid, values: vals, diverging: true, stops: DIV_SHIFT };
  if (sel.diff) {
    const base = sg.layers[q.id];
    if (!base) return { grid: sg, gid: sgid, values: vals, diverging: false, noDiff: true };
    const d = vals.map((v, i) => (ok(v) && ok(base[i]) ? v - base[i] : null));
    return { grid: sg, gid: sgid, values: d, diverging: true, stops: q.scen.kind === "hmech" ? DIV_BAD_HIGH : DIV_GOOD_HIGH };
  }
  return { grid: sg, gid: sgid, values: vals, diverging: false };
}

/* ---------- pages ---------- */
const PAGES = {
  home: renderHome, map: renderMap, treeline: renderTreeline, site: renderSite, decision: renderDecision,
  models: renderModels, method: renderMethod, status: renderStatus, refs: renderRefs, ack: renderAck,
};
const NAV = [["home", "Overview"], ["map", "Map"], ["treeline", "Treeline"], ["site", "Site explorer"], ["decision", "Decision"], ["models", "Models"], ["method", "Method"], ["status", "Status & limits"], ["refs", "References"], ["ack", "Acknowledgments"]];

function route() {
  const page = (location.hash.replace(/^#\//, "").split("?")[0]) || "home";
  const fn = PAGES[page] || renderHome;
  $("nav").innerHTML = NAV.map(([k, t]) => `<a href="#/${k}" class="${k === page ? "active" : ""}">${t}</a>`).join("");
  if (state.mapObj) { state.mapObj.remove(); state.mapObj = null; }
  window.scrollTo(0, 0);
  fn(pageParams());
}

function banner() {
  const M = state.M, rej = Object.values(M.provenance).flatMap((p) => p.rejected || []);
  const grids = Object.keys(M.grids);
  let h = "";
  if (!grids.includes("dense")) h += `<div class="callout info"><strong>Grid in use:</strong> every layer currently comes from the 80-point validation grid. A denser 1044-point grid is being computed; a dense result is adopted only once it covers the grid almost completely.</div>`;
  if (rej.length) h += `<div class="callout"><strong>Dense result rejected:</strong> ${rej.map(esc).join("; ")}</div>`;
  return h;
}

/* ---- Home ---- */
function renderHome() {
  const M = state.M, ss = M.scenario_summary, tl = M.treeline, ae = M.aegis;
  const base = M.baseline_viability_2019 || {};
  const worst = ss ? Object.fromEntries(Object.keys(ss).map((g) => [g, ss[g].ssp585["2100"]])) : {};
  const t585 = tl && tl.summary["ssp585__2100"], t126 = tl && tl.summary["ssp126__2100"];
  const eng = state.meth.engines;
  const stat = (num, cap, sub) => `<div class="card stat"><div class="num">${num}</div><div class="cap">${cap}</div><div class="sub">${sub || ""}</div></div>`;
  view().innerHTML = `<div class="wrap">
    <div class="hero"><img class="hero-logo" src="assets/antar_logo.jpeg" alt="ANTAR"><h1>Armenia reforestation refugia</h1>
    <p class="lead">A hybrid process-statistical framework for finding climate-resilient places to restore forest in Armenia: it models water stress, hydraulic failure, species niches, treeline and scenario-robust planting decisions, and reports each result with the caveats that came with it.</p></div>
    ${banner()}
    <h2>Headline results</h2>
    <div class="grid cols-4">
      ${Object.keys(M.groups).map((g) => stat(pct(base[g], 1), `${esc(M.groups[g].short)} — mean 2019 viability`, worst[g] ? `SSP5-8.5, 2100 ensemble mean: ${pct(worst[g].mean, 1)} (GCM range ${pct(worst[g].min, 1)}–${pct(worst[g].max, 1)})` : "")).join("")}
      ${t585 ? stat(`+${fmt(t585.ensemble_mean_shift_m, 0)} m`, "Climatic treeline shift, SSP5-8.5 by 2100", `GCM range +${fmt(t585.ensemble_min_shift_m, 0)} to +${fmt(t585.ensemble_max_shift_m, 0)} m` + (t126 ? ` · SSP1-2.6: +${fmt(t126.ensemble_mean_shift_m, 0)} m` : "")) : ""}
      ${ae ? stat(String(ae.n_units), "planting units evaluated", `${ae.n_scenarios} scenarios · price of robustness ${fmt(ae.frontier.price_of_robustness, 2)}`) : ""}
    </div>
    <div class="callout"><strong>Read these with care.</strong> Viability stays near its 2019 level across scenarios because the future projections change temperature and precipitation but hold wind, radiation and humidity at 2019 values; that is a limitation of the method, not evidence that climate change is harmless. See <a href="#/status">Status &amp; limits</a>.</div>
    <h2>The six engines</h2>
    <div class="grid cols-3">${eng.map((e) => `<div class="card"><h3>${esc(e.name)}</h3><p class="small muted">${esc(e.status)}</p><a href="#/method" data-jump="${e.id}">How it works →</a></div>`).join("")}</div>
    <h2>Where to look</h2>
    <div class="grid cols-3">
      <div class="card"><h3>Map</h3><p class="small">Pick any quantity, a climate model, an emissions path and a horizon, and see it across the sampled cells.</p><a class="btn" href="#/map">Open the map</a></div>
      <div class="card"><h3>Treeline</h3><p class="small">How far uphill the climatic treeline moves in each of 45 climate-model × scenario × horizon members.</p><a class="btn secondary" href="#/treeline">See treeline change</a></div>
      <div class="card"><h3>Decision</h3><p class="small">A budget-constrained, scenario-robust planting portfolio and its efficient frontier.</p><a class="btn secondary" href="#/decision">See the portfolio</a></div>
    </div></div>`;
  view().querySelectorAll("[data-jump]").forEach((a) => a.addEventListener("click", () => { sessionStorageSafe("jump", a.dataset.jump); }));
}
function sessionStorageSafe(k, v) { try { sessionStorage.setItem(k, v); } catch (e) { /* ignore */ } }

/* ---- Map ---- */
function renderMap(p) {
  const M = state.M, Q = buildQuantities();
  const sg = scenarioGridId(), gcms = sg ? [...new Set(state.G[sg].scenarios.members.map((m) => m.gcm))].sort() : [];
  const sel = { q: Q[p.q] ? p.q : (Q.viab_pine ? "viab_pine" : Object.keys(Q)[0]), mode: p.mode || "base", ssp: p.ssp || "ssp585", hz: +(p.hz || 2100), gcm: p.gcm || "ens", diff: p.diff === "1" };
  if (Q[sel.q] && Q[sel.q].scenarioOnly) sel.mode = "scen";
  const byEngine = {};
  Object.values(Q).forEach((q) => (byEngine[q.engine] = byEngine[q.engine] || []).push(q));
  const opts = ENGINE_ORDER.filter((e) => byEngine[e]).map((e) => `<optgroup label="${esc(ENGINE_TITLE[e] || e)}">${byEngine[e].map((q) => `<option value="${esc(q.id)}">${esc(q.label)}</option>`).join("")}</optgroup>`).join("");
  view().innerHTML = `<div class="wrap">${banner()}<div class="maplayout">
    <aside class="card side"><h3>What to show</h3>
      <label for="q">Quantity</label><select id="q">${opts}</select>
      <div id="scenctl"></div>
      <div id="qinfo" class="small" style="margin-top:12px"></div>
      <div class="legend" id="legend"></div>
      <div id="cellinfo" style="margin-top:14px"></div>
    </aside>
    <div><div id="map"></div><p class="small muted" style="margin-top:8px">Each circle is one sampled grid cell, not a continuous surface. Basemap © OpenStreetMap contributors.</p></div></div></div>`;
  $("#q").value = sel.q;
  const grid0 = state.G[Object.keys(state.G)[0]];
  const map = L.map("map", { zoomControl: true, scrollWheelZoom: true });
  state.mapObj = map;
  L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", { maxZoom: 12, attribution: "© OpenStreetMap" }).addTo(map);
  const lats = grid0.lat, lons = grid0.lon;
  map.fitBounds([[Math.min(...lats) - 0.1, Math.min(...lons) - 0.1], [Math.max(...lats) + 0.1, Math.max(...lons) + 0.1]]);
  let layerGroup = L.layerGroup().addTo(map);

  function scenControls(q) {
    if (!q.scen) { $("#scenctl").innerHTML = ""; return; }
    const sgid = q.scen.kind === "treeline" ? treelineGridId() : scenarioGridId();
    if (!sgid) { $("#scenctl").innerHTML = `<p class="small muted">No scenario results are available for this quantity.</p>`; return; }
    $("#scenctl").innerHTML = `
      ${q.scenarioOnly ? "" : `<label for="mode">Period</label><select id="mode"><option value="base">2019 baseline</option><option value="scen">Future scenario</option></select>`}
      <div id="scen" style="${sel.mode === "scen" || q.scenarioOnly ? "" : "display:none"}">
        <label for="ssp">Emissions path</label><select id="ssp">${Object.entries(SSP).map(([k, v]) => `<option value="${k}">${v}</option>`).join("")}</select>
        <label for="hz">Horizon</label><select id="hz">${HORIZONS.map((h) => `<option value="${h}">${h}</option>`).join("")}</select>
        <label for="gcm">Climate model</label><select id="gcm"><option value="ens">Ensemble mean (5 models)</option>${gcms.map((g) => `<option value="${g}">${g}</option>`).join("")}</select>
        ${q.scen.kind === "treeline" ? "" : `<label class="toggle" style="margin-top:12px"><input type="checkbox" id="diff"> Show change vs 2019</label>`}
      </div>`;
    if ($("#mode")) $("#mode").value = sel.mode;
    $("#ssp").value = sel.ssp; $("#hz").value = sel.hz; $("#gcm").value = sel.gcm;
    if ($("#diff")) $("#diff").checked = sel.diff;
    ["mode", "ssp", "hz", "gcm"].forEach((id) => $("#" + id) && $("#" + id).addEventListener("change", (e) => { sel[id] = e.target.value; if (id === "mode") $("#scen").style.display = sel.mode === "scen" ? "" : "none"; refresh(); }));
    if ($("#diff")) $("#diff").addEventListener("change", (e) => { sel.diff = e.target.checked; refresh(); });
  }

  function refresh() {
    const q = Q[sel.q];
    const s = getSeries(q, sel);
    layerGroup.clearLayers();
    const meta = [`<div><strong>${esc(q.label)}</strong></div>`, `<div class="muted">${esc(q.engine)} · ${gridChip(s ? s.gid : q.grid)}</div>`];
    if (q.placeholder) meta.push(`<div class="callout" style="margin:8px 0;padding:7px 10px">Placeholder input: ${esc(q.placeholder)}</div>`);
    if (!s || !s.values) { $("#qinfo").innerHTML = meta.join("") + `<p class="muted">No data for this selection.</p>`; $("#legend").innerHTML = ""; return; }
    if (s.noDiff) meta.push(`<div class="callout" style="margin:8px 0;padding:7px 10px">No 2019 baseline exists on this grid, so change cannot be shown.</div>`);
    const unit = q.unit || "";
    const sc = makeScale(s.values, { diverging: s.diverging, stops: s.stops, binary: s.binary, unit: s.noUnit ? '' : (q.unit || '') });
    const g = s.grid, n = g.n_cells, radius = n < 200 ? 9 : n < 600 ? 6 : 4;
    const present = s.values.filter(ok);
    for (let i = 0; i < n; i++) {
      const v = s.values[i];
      const mk = L.circleMarker([g.lat[i], g.lon[i]], { radius, weight: 1, color: ok(v) ? "#0b1a12" : "#888", fillColor: ok(v) ? sc.fn(v) : "#bbb", fillOpacity: ok(v) ? 0.92 : 0.35 });
      mk.bindPopup(() => `<strong>${g.lat[i].toFixed(3)}°N, ${g.lon[i].toFixed(3)}°E</strong><br>Elevation ${fmt(g.elev[i], 0)} m<br>${esc(q.label)}: <strong>${fmt(v, 4)}${s.binary ? "" : " " + esc(unit)}</strong><br><a href="#/site?grid=${s.gid}&i=${i}">Open in site explorer →</a>`);
      mk.addTo(layerGroup);
    }
    $("#qinfo").innerHTML = meta.join("") + `<div class="kv" style="margin-top:8px"><span>Cells shown</span><span>${present.length} / ${n}</span><span>Data range</span><span>${fmt(Math.min(...present), 4)} – ${fmt(Math.max(...present), 4)} ${esc(s.binary ? "" : unit)}</span><span>Mean</span><span>${fmt(mean(present), 3)}</span></div>`;
    $("#legend").innerHTML = sc.binary
      ? `<div><span class="swatch"><i style="background:#2f8f5b"></i>yes</span><span class="swatch"><i style="background:#c0583a"></i>no</span></div>`
      : `<div class="bar" style="background:linear-gradient(90deg,${sc.stops.join(",")})"></div><div class="ends"><span>${fmt(sc.lo, 3)}</span><span>${esc(unit)}${s.diverging ? " (centre = 0)" : ""}</span><span>${fmt(sc.hi, 3)}</span></div><p class="small muted" style="margin:4px 0 0">${sc.widened ? '<strong>Scale widened</strong>: the data vary by less than ' + fmt(MIN_SPAN[unit] || 0, 3) + ' ' + esc(unit) + ', so colours are not stretched to make tiny differences look large.' : 'Colour range clipped to the 2nd–98th percentile.'}</p>`;
    setParams({ q: sel.q, mode: sel.mode, ssp: sel.ssp, hz: sel.hz, gcm: sel.gcm, diff: sel.diff ? 1 : 0 });
  }
  $("#q").addEventListener("change", (e) => { sel.q = e.target.value; if (Q[sel.q].scenarioOnly) sel.mode = "scen"; scenControls(Q[sel.q]); refresh(); });
  scenControls(Q[sel.q]); refresh();
}

/* ---- Treeline ---- */
function renderTreeline() {
  const tl = state.M.treeline;
  if (!tl) { view().innerHTML = `<div class="wrap"><p>Treeline change has not been computed.</p></div>`; return; }
  const S = tl.summary, series = Object.keys(SSP).map((ssp) => ({
    name: SSP[ssp], color: SSP_COLOR[ssp],
    points: [[2019, 0, 0, 0]].concat(HORIZONS.map((h) => { const s = S[`${ssp}__${h}`]; return [h, s.ensemble_mean_shift_m, s.ensemble_min_shift_m, s.ensemble_max_shift_m]; })),
  }));
  const rows = Object.keys(SSP).map((ssp) => `<tr><td>${SSP[ssp]}</td>${HORIZONS.map((h) => { const s = S[`${ssp}__${h}`]; return `<td class="num">${s.ensemble_mean_shift_m >= 0 ? "+" : ""}${fmt(s.ensemble_mean_shift_m, 0)} <span class="muted">(${fmt(s.ensemble_min_shift_m, 0)} to ${fmt(s.ensemble_max_shift_m, 0)})</span></td>`; }).join("")}</tr>`).join("");
  view().innerHTML = `<div class="wrap"><h1>Treeline change</h1>
    <p class="lead muted" style="max-width:780px">The climatic treeline is the elevation above which growing-season temperature falls below ${tl.threshold_c} °C. Warming moves it uphill. Below, its shift relative to the 2019 climate across five climate models, three emissions paths and three horizons.</p>
    ${banner()}
    <div class="card"><h3>Mean shift of the climatic treeline across sampled cells (metres, + = uphill)</h3>
      ${Charts.line({ series, xticks: [2019, 2050, 2080, 2100], xlabel: "Horizon", ylabel: "Shift vs 2019 (m)", yfmt: (v) => fmt(v, 0), xfmt: (v) => (v === 2019 ? "2019" : String(v)) })}
      <p class="small muted">Lines are the mean across the five climate models; shaded bands span the lowest and highest model. ${gridChip(tl.grid)}</p></div>
    <h2>Numbers</h2><div class="card tablewrap"><table><thead><tr><th>Emissions path</th>${HORIZONS.map((h) => `<th class="num">${h}</th>`).join("")}</tr></thead><tbody>${rows}</tbody></table>
      <p class="small muted" style="margin-top:8px">Ensemble mean shift in metres, with the range across the five models in brackets.</p></div>
    <div class="callout"><strong>A climatic ceiling, not a forecast.</strong> This is where temperature would permit trees, not where forest will stand: realised treelines lag climate by decades. The ${tl.threshold_c} °C threshold is the global Körner–Paulsen value, not an Armenia calibration, and the lapse rate used is the fitted April–September value (${fmt(tl.gamma_k_per_km, 2)} K/km).</div>
    <div class="callout"><strong>Why a few cells show a drop under warming.</strong> Growing-season temperature is the mean over days that clear 0.9 °C. Warming adds cold early-spring and late-autumn days to that set, which can pull the mean down even though every day warmed. ${pct(tl.frac_pairs_negative, 1)} of cell × model pairs show a negative shift, mostly under low warming. Treat small negative values as roughly "no change".</div>
    <p class="small muted">Method check: the reconstruction reproduces the full pipeline's 2019 growing-season temperature to within ${tl.baseline_check_max_abs_diff_c == null ? "n/a (no overlapping cells)" : fmt(tl.baseline_check_max_abs_diff_c, 4) + " °C"}.</p>
    <p><a class="btn" href="#/map?q=treeline_shift&mode=scen&ssp=ssp585&hz=2100&gcm=ens">See the 2100 SSP5-8.5 shift on the map</a></p></div>`;
}

/* ---- Site explorer ---- */
function renderSite(p) {
  const gid = state.G[p.grid] && state.G[p.grid].scenarios && state.G[p.grid].scenarios.members ? p.grid : scenarioGridId();
  if (!gid) { view().innerHTML = `<div class="wrap"><p>No scenario results are available.</p></div>`; return; }
  const g = state.G[gid], sc = g.scenarios, M = state.M;
  let i = +p.i; if (!(i >= 0 && i < g.n_cells)) i = Math.floor(g.n_cells / 2);
  const grp = p.group && M.groups[p.group] ? p.group : "mesic_diffuse_porous_broadleaf";
  const lonMin = Math.min(...g.lon), lonMax = Math.max(...g.lon), latMin = Math.min(...g.lat), latMax = Math.max(...g.lat);
  const loc = g.lat.map((la, k) => { const x = 10 + ((g.lon[k] - lonMin) / (lonMax - lonMin || 1)) * 280, y = 150 - ((la - latMin) / (latMax - latMin || 1)) * 140; return `<circle data-i="${k}" cx="${x}" cy="${y}" r="${g.n_cells < 200 ? 5.5 : 3.2}" fill="${k === i ? "var(--c3)" : "var(--accent)"}" stroke="var(--bg)" opacity="${k === i ? 1 : 0.65}"><title>${g.lat[k].toFixed(3)}°N ${g.lon[k].toFixed(3)}°E · ${fmt(g.elev[k], 0)} m</title></circle>`; }).join("");
  const byEngine = {};
  for (const [id, arr] of Object.entries(g.layers)) { const m = M.layers[id]; if (!m || !ok(arr[i])) continue; (byEngine[m.engine] = byEngine[m.engine] || []).push([m, arr[i]]); }
  const tbl = ENGINE_ORDER.filter((e) => byEngine[e]).map((e) => `<tr><th colspan="2">${esc(ENGINE_TITLE[e] || e)}</th></tr>` + byEngine[e].map(([m, v]) => `<tr><td>${esc(m.label)}</td><td class="num">${fmt(v, 4)} ${esc(m.binary ? "" : m.unit || "")}</td></tr>`).join("")).join("");

  function perSsp(kind, group, baselineVal) {
    const members = kind === "treeline" ? sc.treeline_members : sc.members;
    const mat = kind === "treeline" ? sc.treeline_shift : kind === "viab" ? sc.viab[group] : sc.hmech[group];
    if (!mat) return [];
    return Object.keys(SSP).map((ssp) => ({
      name: SSP[ssp], color: SSP_COLOR[ssp],
      points: [[2019, baselineVal, null, null]].concat(HORIZONS.map((h) => {
        const vs = members.map((m, k) => (m.ssp === ssp && m.horizon === h ? mat[k][i] : null)).filter(ok);
        return vs.length ? [h, mean(vs), Math.min(...vs), Math.max(...vs)] : [h, null, null, null];
      })),
    }));
  }
  const xf = (v) => String(v);
  const viabS = perSsp("viab", grp, g.layers["viab_" + grp] ? g.layers["viab_" + grp][i] : null);
  const hazS = perSsp("hmech", grp, g.layers["hmech_" + grp] ? g.layers["hmech_" + grp][i] : null);
  const tlS = sc.treeline_members ? perSsp("treeline", null, 0) : [];
  view().innerHTML = `<div class="wrap"><h1>Site explorer</h1>${banner()}
    <div class="grid cols-2" style="grid-template-columns:minmax(260px,330px) 1fr">
      <div class="card"><h3>Choose a cell</h3><svg class="locator" viewBox="0 0 300 160" style="width:100%;background:var(--surface-2);border-radius:8px">${loc}</svg>
        <p class="small muted">Click a dot, or pick below. ${gridChip(gid)}</p>
        <label for="cellsel">Cell</label><select id="cellsel">${g.lat.map((la, k) => `<option value="${k}">${la.toFixed(3)}°N, ${g.lon[k].toFixed(3)}°E · ${fmt(g.elev[k], 0)} m</option>`).join("")}</select>
        <label for="grp">Functional group</label><select id="grp">${Object.keys(M.groups).map((k) => `<option value="${k}">${esc(M.groups[k].label)}</option>`).join("")}</select>
        <div class="kv" style="margin-top:12px"><span>Latitude</span><span>${g.lat[i].toFixed(4)}</span><span>Longitude</span><span>${g.lon[i].toFixed(4)}</span><span>Elevation</span><span>${fmt(g.elev[i], 0)} m</span></div></div>
      <div class="card tablewrap"><h3>Everything computed for this cell</h3><table><tbody>${tbl || "<tr><td>No layers for this cell.</td></tr>"}</tbody></table></div></div>
    <h2>Trajectories for ${esc(groupLabel(grp))}</h2>
    <div class="grid cols-2">
      <div class="card"><h3>Viability</h3>${Charts.line({ series: viabS, xticks: [2019, 2050, 2080, 2100], xlabel: "Horizon", ylabel: "Viability", yfmt: (v) => fmt(v, 3), xfmt: xf, height: 300 })}</div>
      <div class="card"><h3>Hydraulic-failure hazard</h3>${Charts.line({ series: hazS, xticks: [2019, 2050, 2080, 2100], xlabel: "Horizon", ylabel: "Hazard (probability)", yfmt: (v) => fmt(v, 3), xfmt: xf, height: 300 })}</div>
    </div>
    ${tlS.length ? `<div class="card" style="margin-top:14px"><h3>Treeline shift at this cell (+ = uphill)</h3>${Charts.line({ series: tlS, xticks: [2019, 2050, 2080, 2100], xlabel: "Horizon", ylabel: "Shift (m)", yfmt: (v) => fmt(v, 0), xfmt: xf, height: 300 })}</div>` : ""}
    <p class="small muted">Lines are means across the five climate models; shaded bands span the lowest to the highest model. The 2019 point is the single-year baseline run.</p></div>`;
  $("#cellsel").value = i; $("#grp").value = grp;
  const go = (extra) => { location.hash = "#/site?" + new URLSearchParams({ grid: gid, i, group: grp, ...extra }).toString(); };
  $("#cellsel").addEventListener("change", (e) => go({ i: e.target.value }));
  $("#grp").addEventListener("change", (e) => go({ group: e.target.value }));
  view().querySelectorAll("svg.locator circle").forEach((c) => c.addEventListener("click", () => go({ i: c.dataset.i })));
}

/* ---- Decision ---- */
function renderDecision() {
  const ae = state.M.aegis;
  if (!ae) { view().innerHTML = `<div class="wrap"><p>No portfolio has been computed.</p></div>`; return; }
  const sw = ae.budget_sweep, fr = ae.frontier;
  const money = (v) => "$" + (v >= 1e6 ? (v / 1e6).toFixed(v % 1e6 ? 1 : 0) + "M" : v.toLocaleString("en"));
  const identical = sw.every((r) => r.n_units_planted === sw[0].n_units_planted && Math.abs(r.expected - sw[0].expected) < 1e-6);
  const front = fr.lambdas ? Charts.line({ series: [{ name: "Expected benefit", color: "var(--c1)", points: fr.lambdas.map((l, k) => [l, fr.expected[k]]) }, { name: "CVaR (worst 20%)", color: "var(--c3)", dash: "6 4", points: fr.lambdas.map((l, k) => [l, fr.cvar[k]]) }], xlabel: "Risk-aversion λ (0 = mean only, 1 = worst-case only)", ylabel: "Benefit ($/yr)", yfmt: (v) => fmt(v, 0), xticks: fr.lambdas, xfmt: (v) => v.toFixed(1) }) : "";
  const byGroup = {};
  ae.options.forEach((o) => (byGroup[o.group_label] = byGroup[o.group_label] || []).push(o));
  const methods = [...new Map(ae.options.map((o) => [o.intervention, o])).values()];
  view().innerHTML = `<div class="wrap"><h1>Decision: scenario-robust planting portfolio</h1>
    <p class="lead muted" style="max-width:780px">${ae.n_units} candidate units, ${ae.options.length} options (functional group × intervention method), evaluated against ${ae.n_scenarios} climate scenarios. Choose at most one option per unit to maximise a blend of expected benefit and worst-case (CVaR) benefit within a budget.</p>
    ${banner()}
    <div class="grid cols-3">
      <div class="card stat"><div class="num">${fmt(sw[0].expected, 0)}</div><div class="cap">expected benefit, $/yr (at every budget)</div><div class="sub">CVaR ${fmt(sw[0].cvar, 0)}</div></div>
      <div class="card stat"><div class="num">${fmt(fr.price_of_robustness, 2)}</div><div class="cap">price of robustness</div><div class="sub">expected benefit given up to hedge the worst 20% of scenarios</div></div>
      <div class="card stat"><div class="num">${sw[0].n_units_planted}</div><div class="cap">units planted at every budget</div><div class="sub">${gridChip(ae.grid)}</div></div>
    </div>
    ${identical ? `<div class="callout info"><strong>Why every budget gives the same answer.</strong> Even the most expensive option, across all eligible units, costs less than the smallest budget tested, so money is not the binding constraint at this sample size. The choice would start to depend on budget only with many more candidate units.</div>` : ""}
    <div class="callout"><strong>What is and isn't modelled.</strong> ${esc(ae.scope_note)}</div>
    <h2>Mean–CVaR frontier at a $45M budget</h2><div class="card">${front}<p class="small muted">The gap between the two lines is what hedging against bad scenarios costs. They are almost identical here because viability barely differs across the 45 scenarios, so there is little downside to hedge.</p></div>
    <h2>Budget sweep</h2><div class="card tablewrap"><table><thead><tr><th>Budget</th><th class="num">Expected benefit</th><th class="num">CVaR</th><th class="num">Units planted</th></tr></thead><tbody>${sw.map((r) => `<tr><td>${money(r.budget_usd)}</td><td class="num">${fmt(r.expected, 1)}</td><td class="num">${fmt(r.cvar, 1)}</td><td class="num">${r.n_units_planted}</td></tr>`).join("")}</tbody></table></div>
    <h2>Intervention methods and costs</h2><div class="card tablewrap"><table><thead><tr><th>Method</th><th class="num">Cost per hectare</th><th>Cost basis</th></tr></thead><tbody>${methods.map((o) => `<tr><td>${esc(o.intervention.replace(/_/g, " "))}</td><td class="num">${o.cost_per_ha_usd != null ? "$" + o.cost_per_ha_usd.toLocaleString("en") : "—"}</td><td>${chip(o.cost_basis === "blended" ? "warn" : "good", o.cost_basis === "blended" ? "blended average (no individual figure)" : o.cost_basis.replace(/_/g, " "))}</td></tr>`).join("")}</tbody></table>
      <p class="small muted" style="margin-top:8px">Benefit = viability × $${ae.value_per_ha_year_usd}/ha/yr (national ecosystem-services value). Eligible units exclude protected areas and heavily human-modified land.</p></div></div>`;
}

/* ---- Models & validation ---- */
function renderModels() {
  const M = state.M, me = M.meristem, mn = M.mneme, ev = M.ecosystem_validation, vg = M.variogram;
  const speciesRows = me ? Object.entries(me.species).map(([sp, v]) => {
    const folds = (v.boyce_folds || []).map((f) => `<span class="chip ${f >= 0.3 ? "good" : f >= 0 ? "neutral" : "bad"}">${f.toFixed(2)}</span>`).join(" ");
    return `<tr><td><em>${esc(sp)}</em></td><td>${v.status === "fitted" ? chip("good", "fitted") : chip("warn", v.status.replace(/_/g, " "))}</td><td class="num">${v.n_presence}</td><td class="num">${v.boyce_mean != null ? v.boyce_mean.toFixed(2) : "—"}</td><td>${folds || "—"}</td></tr>`;
  }).join("") : "";
  const evRows = ev ? Object.entries(ev.validation).map(([g, v]) => `<tr><td>${esc(groupLabel(g))}</td><td class="num">${ok(v.spearman_rho) ? v.spearman_rho.toFixed(2) : "undefined"}</td><td class="num">${ok(v.p_value) ? v.p_value.toFixed(3) : "—"}</td><td class="num">${v.n}</td><td class="num">${pct(v.mean_observed_forest_cover_fraction, 1)}</td></tr>`).join("") : "";
  let vgHtml = "";
  if (vg && vg.empirical_semivariogram) {
    const pts = vg.empirical_semivariogram.map((b) => [b.lag_km, b.semivariance]), ex = vg.exponential_fit, ln = vg.linear_fit;
    const curves = [];
    const xmax = Math.max(...pts.map((p) => p[0]));
    if (ex) curves.push({ name: `exponential (R² ${ex.r2.toFixed(2)})`, color: "var(--c2)", x0: 0, x1: xmax, fn: (h) => ex.nugget + ex.sill * (1 - Math.exp(-h / ex.range_param_km)) });
    if (ln) curves.push({ name: `straight line (R² ${ln.r2.toFixed(2)})`, color: "var(--c4)", dash: "6 4", x0: 0, x1: xmax, fn: (h) => ln.slope * h + ln.intercept });
    vgHtml = `<div class="card"><h3>Spatial correlation of climatic water deficit</h3>${Charts.scatter({ points: pts, curves, xlabel: "Distance between cells (km)", ylabel: "Semivariance of detrended CWD", pointLabel: (p) => `${p[0]} km: ${p[1]}` })}
      <p class="small">${gridChip(vg.grid)} ${vg.n_points} cells, ${vg.n_pairs.toLocaleString("en")} pairs. <strong>${esc(vg.verdict || "")}</strong></p></div>`;
  }
  view().innerHTML = `<div class="wrap"><h1>Models &amp; validation</h1>${banner()}
    <h2>MERISTEM — species niche</h2>
    ${me ? `<div class="card tablewrap"><table><thead><tr><th>Species</th><th>Status</th><th class="num">Presences</th><th class="num">Mean Boyce</th><th>Per-fold Boyce</th></tr></thead><tbody>${speciesRows}</tbody></table>
      <p class="small muted" style="margin-top:8px">Boyce index under spatial block cross-validation: +1 = predictions track presences, 0 = no better than random, negative = worse. Folds are coloured accordingly.</p></div>
      <div class="callout"><strong>Caveat.</strong> ${esc(me.cwd_caveat)}</div>` : "<p>MERISTEM has not been fitted.</p>"}
    <h2>MNEME — observed dieback</h2>
    ${mn ? `<div class="card"><div class="grid cols-4"><div class="stat"><div class="num">${mn.n_points_valid_kndvi ?? "—"}/${mn.n_points ?? "—"}</div><div class="cap">points with a usable satellite series</div></div><div class="stat"><div class="num">${(mn.n_person_years ?? 0).toLocaleString("en")}</div><div class="cap">person-years</div></div><div class="stat"><div class="num">${mn.n_events ?? "—"}</div><div class="cap">observed dieback events</div></div><div class="stat"><div class="num">${mn.status === "fitted" ? chip("good", "fitted") : chip("bad", "no fit")}</div><div class="cap">${esc((mn.status || "").replace(/_/g, " "))}</div></div></div>
      <p class="small muted" style="margin-top:10px">${esc(mn.scope_note || "")}</p></div>` : "<p>MNEME has not been run.</p>"}
    <h2>Does predicted viability track mapped forest?</h2>
    ${ev ? `<div class="card tablewrap"><table><thead><tr><th>Group</th><th class="num">Spearman ρ</th><th class="num">p</th><th class="num">cells</th><th class="num">mean mapped cover</th></tr></thead><tbody>${evRows}</tbody></table>
      <p class="small muted" style="margin-top:8px">REFUGIUM viability against the real Ecosystem Map of Armenia's forest classes in a ${ev.window_radius_m} m window. ${ev.n_outside} of ${ev.n_cells_total} cells fall outside Armenia's border and are excluded. ρ is undefined for pine because no mapped pine cover falls in any sampled window.</p></div>` : ""}
    <h2>Sampling design</h2>${vgHtml || "<p>No variogram has been estimated.</p>"}</div>`;
}

/* ---- Methodology ---- */
function renderMethod() {
  let plain = false; try { plain = localStorage.getItem("antar_plain") === "1"; } catch (e) { /* ignore */ }
  const m = state.meth;
  const body = () => `<div class="callout info"><strong>${plain ? "Plain-language view." : "Technical view."}</strong> ${esc(plain ? m.intro.plain : m.intro.technical)}</div>` +
    m.engines.map((e) => `<section class="card engine" id="${e.id}" style="margin-bottom:14px"><h2 style="margin-top:0">${esc(e.name)}</h2><p>${chip("neutral", "status")} <span class="small">${esc(e.status)}</span></p><ul class="${plain ? "plainbox" : ""}">${(plain ? e.plain : e.technical).map((t) => `<li>${esc(t)}</li>`).join("")}</ul></section>`).join("");
  view().innerHTML = `<div class="wrap"><h1>Methodology</h1><p class="muted" style="max-width:760px">What is actually implemented, engine by engine. Switch to the plain-language view for an everyday-analogy explanation of each mechanism; the technical view keeps the equations and parameters.</p>
    <label class="toggle"><input type="checkbox" id="plain"> Explain it simply</label><div id="mbody" style="margin-top:14px"></div></div>`;
  $("#plain").checked = plain;
  const draw = () => { $("#mbody").innerHTML = body(); };
  draw();
  $("#plain").addEventListener("change", (e) => { plain = e.target.checked; try { localStorage.setItem("antar_plain", plain ? "1" : "0"); } catch (x) { /* ignore */ } draw(); });
  let j = null; try { j = sessionStorage.getItem("jump"); sessionStorage.removeItem("jump"); } catch (e) { /* ignore */ }
  if (j && document.getElementById(j)) document.getElementById(j).scrollIntoView();
}

/* ---- Status & limits ---- */
function renderStatus() {
  const m = state.meth, M = state.M;
  const prov = Object.entries(M.provenance).map(([k, p]) => `<tr><td>${esc(p.dataset)}</td><td><code>${esc(p.file || "missing")}</code></td><td>${p.grid ? gridChip(p.grid) : "—"}</td><td class="num">${p.n_cells ?? "—"}</td><td>${(p.rejected || []).map((r) => `<span class="small">${esc(r)}</span>`).join("<br>")}</td></tr>`).join("");
  view().innerHTML = `<div class="wrap"><h1>Status &amp; limits</h1><p class="muted" style="max-width:760px">What is fitted, what is reduced in scope, and what is still a placeholder. Nothing is hidden because it is incomplete.</p>
    <h2>Engines</h2><div class="card tablewrap"><table><thead><tr><th>Engine</th><th>State</th></tr></thead><tbody>${m.engines.map((e) => `<tr><td><a href="#/method" data-jump="${e.id}">${esc(e.name)}</a></td><td>${esc(e.status)}</td></tr>`).join("")}</tbody></table></div>
    <h2>Placeholders and assumptions</h2><div class="card tablewrap"><table><thead><tr><th>Quantity</th><th>Value used</th><th>Affects</th><th>State</th></tr></thead><tbody>${m.placeholders.map((p) => `<tr><td>${esc(p.name)}</td><td>${esc(p.value)}</td><td>${esc(p.effect)}</td><td>${esc(p.state)}</td></tr>`).join("")}</tbody></table></div>
    <h2>Known gaps</h2><div class="card"><ul>${m.known_gaps.map((g) => `<li>${esc(g)}</li>`).join("")}</ul></div>
    <h2>Data lineage — which file each result comes from</h2><div class="card tablewrap"><table><thead><tr><th>Result</th><th>File used</th><th>Grid</th><th class="num">Cells</th><th>Notes</th></tr></thead><tbody>${prov}</tbody></table>
      <p class="small muted" style="margin-top:8px">Generated ${esc(M.generated)}. A dense-grid file is used only if it covers at least ${Math.round(100 * M.dense_min_coverage)}% of the dense grid; otherwise it is listed here as rejected.</p></div></div>`;
  view().querySelectorAll("[data-jump]").forEach((a) => a.addEventListener("click", () => sessionStorageSafe("jump", a.dataset.jump)));
}

/* ---- References ---- */
function renderRefs() {
  const refs = state.M.references;
  view().innerHTML = `<div class="wrap"><h1>References</h1><p class="muted">Every dataset behind the results, with the citation and licence recorded when it was obtained (${refs.length} entries).</p>
    <input type="search" id="rq" placeholder="Filter by dataset, source or citation…" style="width:min(520px,100%)"><div id="rl" class="card" style="margin-top:12px"></div></div>`;
  const draw = (q) => {
    const f = refs.filter((r) => !q || [r.dataset, r.source, r.citation].join(" ").toLowerCase().includes(q.toLowerCase()));
    $("#rl").innerHTML = f.length ? f.map((r) => `<div class="refitem"><div class="t">${esc(r.dataset || "")}</div><div class="small">${esc(r.citation)}</div><div class="small muted">${r.url ? `<a href="${esc(r.url)}" target="_blank" rel="noopener">${esc(r.url)}</a> · ` : ""}${r.license ? "Licence: " + esc(r.license) + " · " : ""}${r.accessed ? "accessed " + esc(r.accessed) : ""}</div></div>`).join("") : `<p class="muted">No matches.</p>`;
  };
  draw(""); $("#rq").addEventListener("input", (e) => draw(e.target.value));
}

/* ---- Acknowledgments ---- */
function renderAck() {
  view().innerHTML = `<div class="wrap"><h1>Acknowledgments</h1>
    <!-- prior project name/attribution to be filled in by hand -->
    <div class="card"><p>This section is intentionally left to be completed by the author.</p>
    <p class="muted small">Data providers and their licences are listed on the <a href="#/refs">References</a> page.</p></div></div>`;
}

/* ---------- boot ---------- */
async function init() {
  view().innerHTML = `<div class="loading">Loading results…</div>`;
  try {
    const get = async (u) => { const r = await fetch(u); if (!r.ok) throw new Error(u + " → HTTP " + r.status); return r.json(); };
    state.M = await get("data/manifest.json");
    state.meth = await get("data/methodology.json");
    for (const [gid, info] of Object.entries(state.M.grids)) state.G[gid] = await get("data/" + info.file);
  } catch (e) {
    view().innerHTML = `<div class="wrap"><div class="callout bad"><strong>Could not load the data.</strong> ${esc(e.message)}<br>Build it with <code>python3 ui/build_data.py</code> and serve the folder over HTTP (browsers block <code>fetch</code> on <code>file://</code>): <code>python3 -m http.server 8000 --directory ui</code>, then open <code>http://localhost:8000</code>.</div></div>`;
    return;
  }
  window.addEventListener("hashchange", route);
  route();
}
function toggleTheme() {
  const cur = document.documentElement.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  const next = cur === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem("antar_theme", next); } catch (e) { /* ignore */ }
}
try { const t = localStorage.getItem("antar_theme"); if (t) document.documentElement.dataset.theme = t; } catch (e) { /* ignore */ }
document.addEventListener("DOMContentLoaded", () => { $("#theme").addEventListener("click", toggleTheme); init(); });
