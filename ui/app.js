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
/* Which grid a result comes from. The node count is not in the label: it differs by quantity (terrain alone reaches nodes the water balance cannot), so every card and panel states its own. */
const gridChip = (gid) => chip(gid === "dense" ? "good" : "neutral", gid === "dense" ? "dense grid" : "validation grid");
const groupLabel = (g) => (state.M.groups[g] ? state.M.groups[g].label : g);
function pageParams() { const q = location.hash.split("?")[1] || ""; return Object.fromEntries(new URLSearchParams(q)); }
function setParams(p) { const base = location.hash.split("?")[0]; history.replaceState(null, "", base + "?" + new URLSearchParams(p).toString()); }
const denseFirst = (ids) => ids.slice().sort((a, b) => (b === "dense") - (a === "dense"));
function scenarioGridId() { return denseFirst(Object.keys(state.G).filter((g) => state.G[g].scenarios && state.G[g].scenarios.members))[0]; }
function treelineGridId() { return denseFirst(Object.keys(state.G).filter((g) => state.G[g].scenarios && state.G[g].scenarios.treeline_members))[0]; }

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

const ENGINE_ORDER = ["TOPOHYDRO", "XYLEM", "REFUGIUM", "Treeline", "Ecosystem map"];
const ENGINE_TITLE = { TOPOHYDRO: "Climate & water balance (TOPOHYDRO)", XYLEM: "Hydraulic-failure hazard (XYLEM)", REFUGIUM: "Viability & refugia (REFUGIUM)", Treeline: "Treeline", "Ecosystem map": "Observed land cover (Ecosystem Map)" };

/* ---------- pages ---------- */
const PAGES = {
  home: renderHome, map: renderMap, treeline: renderTreeline, site: renderSite, decision: renderDecision,
  models: renderModels, method: renderMethod, status: renderStatus, refs: renderRefs, ack: renderAck, author: renderAuthor,
};
const NAV = [["home", "Overview"], ["map", "Map"], ["treeline", "Treeline"], ["site", "Place explorer"], ["decision", "Decision"], ["models", "Models"], ["method", "Method"], ["status", "Status & limits"], ["refs", "References"], ["ack", "Acknowledgments"]];

/* Phones: the ten pages live in a menu opened by the button in the bar (large tap targets); it closes on a choice, on Escape and on a tap outside. */
function setMenu(open) {
  const btn = $("#menu"), nav = $("#mainnav");
  if (!btn || !nav) return;
  nav.classList.toggle("open", open);
  btn.setAttribute("aria-expanded", String(open));
  btn.setAttribute("aria-label", open ? "Close the menu" : "Open the menu");
  document.body.classList.toggle("menu-open", open);
}
function initMenu() {
  const btn = $("#menu");
  btn.addEventListener("click", () => setMenu(btn.getAttribute("aria-expanded") !== "true"));
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && btn.getAttribute("aria-expanded") === "true") { setMenu(false); btn.focus(); } });
  document.addEventListener("click", (e) => { if (btn.getAttribute("aria-expanded") === "true" && !e.target.closest("#mainnav, #menu")) setMenu(false); });
  window.matchMedia("(min-width: 761px)").addEventListener("change", () => setMenu(false));
}

/* Four inner pages open with a credited photograph behind their title. The page renderers are untouched: once a page has put its title in the
   view, that title is moved into the banner. Pages whose view is a tool (map, place explorer) and the record pages keep the plain header. */
const PAGE_PHOTO = { treeline: "aragats", decision: "dilijan-beech", models: "khosrov", method: "dilijan-ridge" };
function bannerize() {
  const file = state.pagePhoto, v = view();
  if (!file || !v || v.querySelector(".pageban")) return;
  const wrap = v.querySelector(":scope > .wrap"), h1 = wrap && wrap.querySelector(":scope > h1");
  const ph = (state.M.photos || []).find((x) => x.file === file);
  if (!h1 || !ph) return;
  const ban = document.createElement("figure");
  ban.className = "pageban";
  ban.innerHTML = `<img src="assets/photos/${ph.file}-800.jpg" srcset="assets/photos/${ph.file}-800.jpg 800w, assets/photos/${ph.file}-1600.jpg 1600w" sizes="100vw" width="${ph.original_px[0]}" height="${ph.original_px[1]}" alt="" decoding="async">
    <div class="wrap pageban-in"></div>
    <figcaption class="wrap"><span class="credit">${esc(ph.caption)}, ${esc(ph.place)} \u00B7 Photo: ${esc(ph.author)} \u00B7 <a href="${esc(ph.licence_url)}" target="_blank" rel="noopener">${esc(ph.licence)}</a> \u00B7 <a href="${esc(ph.source_page)}" target="_blank" rel="noopener">${esc(ph.source)}<span class="sr-only"> (opens in a new tab)</span></a> \u00B7 resized</span></figcaption>`;
  ban.querySelector(".pageban-in").appendChild(h1);
  v.insertBefore(ban, wrap);
  document.body.classList.add("has-ban");
}

function route() {
  setMenu(false);
  const page = (location.hash.replace(/^#\//, "").split("?")[0]) || "home";
  const fn = PAGES[page] || renderHome;
  $("nav").innerHTML = NAV.map(([k, t]) => `<a href="#/${k}" class="${k === page ? "active" : ""}">${t}</a>`).join("");
  const activeLink = $("nav a.active"); if (activeLink) activeLink.scrollIntoView({ inline: "center", block: "nearest" });   // phones: keep the current page visible in the scrolling nav
  if (state.banCleanup) { state.banCleanup(); state.banCleanup = null; }
  document.body.classList.remove("has-ban");
  state.pagePhoto = PAGE_PHOTO[page] || null;
  if (state.revealCleanup) { state.revealCleanup(); state.revealCleanup = null; }
  if (state.heroCleanup) { state.heroCleanup(); state.heroCleanup = null; }
  document.body.classList.toggle("is-home", fn === renderHome);
  if (state.mapCleanup) { state.mapCleanup(); state.mapCleanup = null; }
  if (state.mapObj) { state.mapObj.remove(); state.mapObj = null; }
  window.scrollTo(0, 0);
  syncTopbar();
  fn(pageParams());
  if (state.pagePhoto) { bannerize(); const mo = new MutationObserver(bannerize); mo.observe(view(), { childList: true }); state.banCleanup = () => mo.disconnect(); }
  view().focus({ preventScroll: true });          // screen readers land on the new page's content
}
function syncTopbar() { $("#topbar").classList.toggle("solid", window.scrollY > 24); }

function banner() {
  const M = state.M, rej = Object.values(M.provenance).flatMap((p) => p.rejected || []);
  const grids = Object.keys(M.grids);
  let h = "";
  if (!grids.includes("dense")) h += `<div class="callout info"><strong>Grid in use:</strong> every layer currently comes from the coarse validation grid (a sample of nodes about 37 km apart). A denser Armenia-only grid is being computed; a dense result is adopted only once it covers the grid almost completely.</div>`;
  if (rej.length) h += `<div class="callout"><strong>Dense result rejected:</strong> ${rej.map(esc).join("; ")}</div>`;
  return h;
}

/* ---- Home ---- */
/* The state of an engine in this release, as the flow chart on the Method page draws it (both read it from data/methodology.json). */
const ENGINE_STATE = { run: ["good", "Run"], reduced: ["warn", "Reduced scope"], built: ["neutral", "Built, not applied"] };
const engineState = (k) => (ENGINE_STATE[k] ? chip(ENGINE_STATE[k][0], ENGINE_STATE[k][1]) : "");
const ICON = (() => {
  const svg = (d) => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">${d}</svg>`;
  return {
    map: svg('<path d="M9 4 3 6.5v13L9 17l6 2.5 6-2.5v-13L15 7 9 4z"/><path d="M9 4v13M15 7v12.5"/>'),
    pin: svg('<path d="M12 21s7-6.2 7-11.5a7 7 0 0 0-14 0C5 14.8 12 21 12 21z"/><circle cx="12" cy="9.5" r="2.4"/>'),
    peak: svg('<path d="m3 19 6.5-11 3.5 5.5 2-3L21 19H3z"/><path d="M9.5 8 8.2 10.2l1.3-.7 1.3.9L9.5 8z"/>'),
    scale: svg('<path d="M12 4v16M7 20h10"/><path d="M5 8h14"/><path d="m5 8-3 7a3.4 3 0 0 0 6 0L5 8zM19 8l-3 7a3.4 3 0 0 0 6 0l-3-7z"/>'),
  };
})();

function renderHome() {
  const M = state.M;
  const engines = state.meth.engines.filter((e) => e.id !== "treeline");        // treeline change is part of MERISTEM, not a seventh engine
  const photos = M.photos || [];
  const sgid = scenarioGridId(), tgid = treelineGridId();
  const tlv = tgid && state.G[tgid].layers.treeline_2019 ? state.G[tgid].layers.treeline_2019.filter(ok) : [];
  const treelineMean = tlv.length ? mean(tlv) : null;
  const nMembers = M.scenario_summary ? Object.values(M.scenario_summary).reduce((n, g) => Math.max(n, Object.values(g).reduce((m, s) => m + Object.values(s).reduce((k, h) => k + (h.n_gcms || 0), 0), 0)), 0) : 0;
  const ledger = `<dl class="ledger">
    <div><dt>Model nodes</dt><dd>${(M.node_counts && M.node_counts.viability) || (sgid ? M.grids[sgid].n_cells : "\u2014")}<small>${sgid === "dense" ? "dense grid" : "validation grid"}</small></dd></div>
    <div><dt>Climate members</dt><dd>${nMembers || "\u2014"}<small>5 models \u00D7 3 paths \u00D7 3 horizons</small></dd></div>
    <div><dt>Species groups</dt><dd>${Object.keys(M.groups).length}<small>hydraulic traits</small></dd></div>
    <div><dt>Engines</dt><dd>${engines.length}<small>chained, each tested</small></dd></div>
    <div><dt>Datasets cited</dt><dd>${(M.references || []).length}<small>with licence and date</small></dd></div>
  </dl>`;
  view().innerHTML = `
  <section class="hero-full" aria-labelledby="hero-title">
    <div class="hero-media"><div class="hero-frame">
      <video id="hero-video" class="hero-video" poster="assets/hero-poster.jpg" muted loop playsinline preload="auto" disablepictureinpicture disableremoteplayback aria-hidden="true" tabindex="-1">
        <source src="assets/hero.mp4" type="video/mp4"><img src="assets/hero-poster.jpg" alt="ANTAR">
      </video>
    </div></div>
    <div class="hero-copy">
      <div class="hero-text">
        <span class="rule" aria-hidden="true"></span>
        <h1 id="hero-title">ANTAR — Assessment of Niche, Treeline &amp; Analogue Refugia</h1>
        <p class="hero-sub">Climate-resilient places to restore forest in Armenia: water stress, hydraulic failure, species niches, treeline change and robust planting decisions.</p>
      </div>
      <div class="hero-cta">
        <a class="btn" href="#/map">Open the map</a><a class="btn secondary" href="#/method">How it works</a>
        <button type="button" id="hero-toggle" class="hero-toggle" aria-label="Pause the logo animation" title="Pause the logo animation"></button>
      </div>
    </div>
    <button type="button" id="scroll-cue" class="scroll-cue" aria-label="Scroll to the overview"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M6 9l6 6 6-6"/></svg></button>
  </section>
  <div class="wrap">
    <section class="sec reveal" id="about">
      <div class="sec-head"><span class="sec-no">01</span><h2>Overview</h2></div>
      <div class="about">
        <p class="lead">A hybrid process-statistical framework for finding climate-resilient places to restore forest in Armenia: it models water stress, hydraulic failure, species niches, treeline and scenario-robust planting decisions, and reports each result with the caveats that came with it.</p>
        ${ledger}
        <div>${banner()}</div>
      </div>
    </section>
  </div>
  <section class="band reveal" id="results" aria-labelledby="res-h"></section>
  <div class="wrap">
    ${photos.length ? `<section class="sec reveal" id="landscapes" aria-labelledby="land-h">
      <div class="sec-head"><span class="sec-no">03</span><h2 id="land-h">Landscapes</h2></div>
      <p class="lead">Four kinds of Armenian landscape that the results speak about, from closed broadleaf forest to the open mountain above it. The photographs illustrate these landscape types: they are not model outputs, and none was used as model input.</p>
      <div class="landscapes">${photos.map((ph, i) => Home.photoFigure(ph, i)).join("")}</div>
      <p class="small muted">Third-party photographs under Creative Commons licences, resized for the page. The author, licence and source are under each image and listed on the <a href="#/ack">Acknowledgments</a> page.</p>
    </section>` : ""}
    <section class="sec reveal">
      <div class="sec-head"><span class="sec-no">04</span><h2>The ${["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"][engines.length] || engines.length} engines</h2></div>
      <div class="grid cols-3 engines-grid">${engines.map((e, i) => { const [ttl, sub] = e.name.split(" \u2014 "); return `<div class="card eng"><div class="eng-top"><span class="eng-no">E${i + 1}</span>${engineState(e.state)}</div><h3>${esc(ttl)}</h3>${sub ? `<p class="eng-sub">${esc(sub)}</p>` : ""}<p class="small">${esc(e.status)}</p><a href="#/method" data-jump="${e.id}">How it works \u2192</a></div>`; }).join("")}</div>
    </section>
    <section class="sec reveal">
      <div class="sec-head"><span class="sec-no">05</span><h2>Where to look</h2></div>
      <div class="grid cols-4 look-grid">
        <a class="card look" href="#/map"><span class="look-ic">${ICON.map}</span><h3>Map</h3><p class="small">Pick any quantity, a climate model, an emissions path and a horizon, and see it coloured across the whole country, with today's forest cover and marz borders.</p><span class="look-go">Open the map \u2192</span></a>
        <a class="card look" href="#/site"><span class="look-ic">${ICON.pin}</span><h3>Place explorer</h3><p class="small">Click anywhere in Armenia, or choose a marz, to see every quantity for that place today and under each emissions path.</p><span class="look-go">Explore a place \u2192</span></a>
        <a class="card look" href="#/treeline"><span class="look-ic">${ICON.peak}</span><h3>Treeline</h3><p class="small">How far uphill the climatic treeline moves in each of 45 climate-model \u00D7 scenario \u00D7 horizon members.</p><span class="look-go">See treeline change \u2192</span></a>
        <a class="card look" href="#/decision"><span class="look-ic">${ICON.scale}</span><h3>Decision</h3><p class="small">A budget-constrained, scenario-robust planting portfolio and its efficient frontier.</p><span class="look-go">See the portfolio \u2192</span></a>
      </div>
    </section>
    <section class="sec reveal">
      <div class="sec-head"><span class="sec-no">06</span><h2>About the author</h2></div>
      <div class="author-card">
        <img src="assets/author.jpg" width="720" height="720" alt="Portrait of Narek Ohanyan" loading="lazy">
        <div><p class="author-name">Narek Ohanyan</p>
          <p>Narek Ohanyan is a young climate leader from Armenia and a climate &amp; environmental researcher at the AUA Acopian Center for the Environment. His current research focuses on modeling forest climate resilience.</p>
          <p><a href="#/author">Read the full biography \u2192</a></p></div>
      </div>
    </section>
  </div>`;
  view().querySelectorAll("[data-jump]").forEach((a) => a.addEventListener("click", () => { sessionStorageSafe("jump", a.dataset.jump); }));
  const q0 = pageParams();
  Home.mountBand($("#results"), { M, treelineMean, gridHtml: sgid ? gridChip(sgid) : "", start: { ssp: q0.ssp, hz: q0.hz }, onChange: (c) => setParams({ ssp: c.ssp, hz: c.hz }) });
  Home.mountPhotos(view(), photos);
  initHeroVideo();
  initReveal();
}
function sessionStorageSafe(k, v) { try { sessionStorage.setItem(k, v); } catch (e) { /* ignore */ } }

/* Sections fade up as they enter the viewport. If the observer never reports (hidden tab, unsupported), everything is shown. */
function initReveal() {
  const els = [...view().querySelectorAll(".reveal")];
  if (!els.length) return;
  if (matchMedia("(prefers-reduced-motion: reduce)").matches || !("IntersectionObserver" in window)) { els.forEach((e) => e.classList.add("in")); return; }
  let fired = false;
  const io = new IntersectionObserver((es) => { fired = true; es.forEach((e) => { if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); } }); }, { rootMargin: "0px 0px -8% 0px", threshold: 0.05 });
  els.forEach((e) => io.observe(e));
  const t = setTimeout(() => { if (!fired) els.forEach((e) => e.classList.add("in")); }, 1500);
  state.revealCleanup = () => { io.disconnect(); clearTimeout(t); };
}

/* Hero video: muted, looping, never the only way to see the name (the heading says it). Autoplay is skipped when the
   visitor asks for reduced motion or Save-Data; the visitor can always pause or play; it also pauses off-screen. */
function initHeroVideo() {
  const v = $("#hero-video"), btn = $("#hero-toggle");
  if (!v || !btn) return;
  const mq = matchMedia("(prefers-reduced-motion: reduce)"), conn = navigator.connection || {};
  const PLAY = '<svg viewBox="0 0 24 24" width="20" height="20" aria-hidden="true" focusable="false"><path d="M8 5.5v13a1 1 0 0 0 1.5.86l10-6.5a1 1 0 0 0 0-1.72l-10-6.5A1 1 0 0 0 8 5.5z" fill="currentColor"/></svg>';
  const PAUSE = '<svg viewBox="0 0 24 24" width="20" height="20" aria-hidden="true" focusable="false"><rect x="6" y="5" width="4" height="14" rx="1.2" fill="currentColor"/><rect x="14" y="5" width="4" height="14" rx="1.2" fill="currentColor"/></svg>';
  let wantPlay = !(mq.matches || conn.saveData), visible = true;
  const label = (playing) => { const t = playing ? "Pause the logo animation" : "Play the logo animation"; btn.innerHTML = playing ? PAUSE : PLAY; btn.setAttribute("aria-label", t); btn.title = t; };
  // a rejected play() (hidden tab, iOS low-power mode) only shows the play button; the wish to play is kept and retried
  const sync = () => { if (wantPlay && visible && document.visibilityState !== "hidden") v.play().catch(() => label(false)); else v.pause(); };
  v.muted = v.defaultMuted = true;
  v.addEventListener("play", () => label(true));
  v.addEventListener("pause", () => label(false));
  const cue = $("#scroll-cue");
  if (cue) cue.addEventListener("click", () => $("#about").scrollIntoView({ behavior: mq.matches ? "auto" : "smooth", block: "start" }));
  const src = v.querySelector("source");
  if (src) src.addEventListener("error", () => { btn.hidden = true; });          // file missing: the poster image stays
  btn.addEventListener("click", () => { wantPlay = v.paused; sync(); });
  const onMotion = () => { if (mq.matches) { wantPlay = false; sync(); } };
  mq.addEventListener("change", onMotion);
  const io = new IntersectionObserver((es) => { visible = es[0].isIntersecting; sync(); }, { threshold: 0.25 });
  io.observe(v);
  const onVis = () => sync();
  document.addEventListener("visibilitychange", onVis);
  label(false);
  sync();
  state.heroCleanup = () => { io.disconnect(); document.removeEventListener("visibilitychange", onVis); mq.removeEventListener("change", onMotion); v.pause(); };
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
    <p class="lead muted">The climatic treeline is the elevation above which growing-season temperature falls below ${tl.threshold_c} °C. Warming moves it uphill. Below, its shift relative to the 2019 climate across five climate models, three emissions paths and three horizons.</p>
    ${banner()}
    <div class="card"><h3>Mean shift of the climatic treeline across sampled cells (metres, + = uphill)</h3>
      ${Charts.line({ series, xticks: [2019, 2050, 2080, 2100], xlabel: "Horizon", ylabel: "Shift vs 2019 (m)", yfmt: (v) => fmt(v, 0), xfmt: (v) => (v === 2019 ? "2019" : String(v)) })}
      <p class="small muted">Lines are the mean across the five climate models; shaded bands span the lowest and highest model. ${gridChip(tl.grid)}</p></div>
    <h2>Numbers</h2><div class="card tablewrap"><table><thead><tr><th>Emissions path</th>${HORIZONS.map((h) => `<th class="num">${h}</th>`).join("")}</tr></thead><tbody>${rows}</tbody></table>
      <p class="small muted" style="margin-top:8px">Ensemble mean shift in metres, with the range across the five models in brackets.</p></div>
    <div class="callout"><strong>A climatic ceiling, not a forecast.</strong> This is where temperature would permit trees, not where forest will stand: realised treelines lag climate by decades. The ${tl.threshold_c} °C threshold is the global Körner–Paulsen value, not an Armenia calibration, and the lapse rate used is the fitted April–September value (${fmt(tl.gamma_k_per_km, 2)} K/km).</div>
    <div class="callout"><strong>Why a few cells show a drop under warming.</strong> Growing-season temperature is the mean over days that clear 0.9 °C. Warming adds cold early-spring and late-autumn days to that set, which can pull the mean down even though every day warmed. ${pct(tl.frac_pairs_negative, 1)} of cell × model pairs show a negative shift, mostly under low warming. Treat small negative values as roughly "no change".</div>
    <p><a class="btn" href="#/map?q=treeline_shift&mode=scen&ssp=ssp585&hz=2100&gcm=ens">See the 2100 SSP5-8.5 shift on the map</a></p></div>`;
}

/* ---- Decision ---- */
function renderDecision() {
  const ae = state.M.aegis;
  if (!ae) {
    const why = (state.M.provenance.aegis && state.M.provenance.aegis.withheld) || [];
    view().innerHTML = `<div class="wrap"><h1>Decision</h1>
      <div class="callout"><strong>The planting portfolio is being recomputed.</strong> The earlier portfolio was built on scenario results that have since been corrected, so it is withheld rather than shown from superseded inputs.${why.length ? ` <span class="small">(${why.map((w) => esc(w.replace(/^[\w.-]+\.yaml: /, ""))).join("; ")})</span>` : ""}</div>
      <h2>What it will show</h2>
      <p class="lead">For each planting unit, a choice of species group and planting method that maximises a mix of the expected benefit and the average benefit in the worst 20% of the 45 climate-model \u00D7 path \u00D7 horizon scenarios (conditional value at risk, \u03B1 = 0.8, \u03BB = 0.5), under a budget. It reports how much expected benefit the robust choice gives up (the price of robustness) and the efficient frontier between the two.</p>
      <p>The assumptions behind benefit and cost, and what is still a placeholder, are on the <a href="#/method" data-jump="aegis">Method</a> and <a href="#/status">Status &amp; limits</a> pages.</p></div>`;
    view().querySelectorAll("[data-jump]").forEach((a) => a.addEventListener("click", () => sessionStorageSafe("jump", a.dataset.jump)));
    return;
  }
  const sw = ae.budget_sweep, fr = ae.frontier;
  const money = (v) => "$" + (v >= 1e6 ? (v / 1e6).toFixed(v % 1e6 ? 1 : 0) + "M" : v.toLocaleString("en"));
  const identical = sw.every((r) => r.n_units_planted === sw[0].n_units_planted && Math.abs(r.expected - sw[0].expected) < 1e-6);
  const costs = (ae.cost_table || []).map((c) => c.cost_per_ha).filter(ok), maxCost = costs.length ? Math.max(...costs) : null;
  const eligible = ae.n_eligible_units != null ? ae.n_eligible_units : sw[0].n_units_planted, unitHa = ae.unit_area_ha != null ? ae.unit_area_ha : 1;
  const minBudget = Math.min(...sw.map((r) => r.budget_usd)), maxTotal = maxCost != null ? eligible * unitHa * maxCost : null;
  const nonBinding = maxTotal != null && maxTotal < minBudget;
  const costBasis = { sourced_armenian_cost_data: ["good", "Armenian cost data, range shown"], sourced_cheapest_mix: ["neutral", "World Bank cheapest option mix, shared"], sourced_most_expensive: ["neutral", "World Bank most expensive option"], blended: ["warn", "blended average, no individual figure"] };
  const ct = Object.fromEntries((ae.cost_table || []).map((c) => [c.name, c]));
  const front = fr.lambdas ? Charts.line({ series: [{ name: "Expected benefit", color: "var(--c1)", points: fr.lambdas.map((l, k) => [l, fr.expected[k]]) }, { name: "CVaR (worst 20%)", color: "var(--c3)", dash: "6 4", points: fr.lambdas.map((l, k) => [l, fr.cvar[k]]) }], xlabel: "Risk-aversion λ (0 = mean only, 1 = worst-case only)", ylabel: "Benefit ($/yr)", yfmt: (v) => fmt(v, 0), xticks: fr.lambdas, xfmt: (v) => v.toFixed(1) }) : "";
  const byGroup = {};
  ae.options.forEach((o) => (byGroup[o.group_label] = byGroup[o.group_label] || []).push(o));
  const methods = [...new Map(ae.options.map((o) => [o.intervention, o])).values()];
  view().innerHTML = `<div class="wrap"><h1>Decision: scenario-robust planting portfolio</h1>
    <p class="lead muted">${ae.n_units} candidate units, ${ae.options.length} options (functional group × intervention method), evaluated against ${ae.n_scenarios} climate scenarios. Choose at most one option per unit to maximise a blend of expected benefit and worst-case (CVaR) benefit within a budget.</p>
    ${banner()}
    <div class="grid cols-3">
      <div class="card stat"><div class="num">${fmt(sw[0].expected, 0)}</div><div class="cap">expected benefit, $/yr (at every budget)</div><div class="sub">CVaR ${fmt(sw[0].cvar, 0)}</div></div>
      <div class="card stat"><div class="num">${fmt(fr.price_of_robustness, 2)}</div><div class="cap">price of robustness</div><div class="sub">expected benefit given up to hedge the worst 20% of scenarios</div></div>
      <div class="card stat"><div class="num">${sw[0].n_units_planted}</div><div class="cap">units given an option${identical ? ", at every budget" : ""}</div><div class="sub">${eligible} of ${ae.n_units} candidate units are eligible · ${gridChip(ae.grid)}</div></div>
    </div>
    ${identical ? `<div class="callout info"><strong>Why every budget gives the same answer.</strong> Each unit is ${fmt(unitHa, 0)} ha (one model node) and ${eligible} are eligible. ${nonBinding ? `Even the most expensive method on all of them costs ${money(maxTotal)}, less than the smallest budget tested (${money(minBudget)}); those budgets are the scale of a 50,000 ha programme (World Bank 2023), so money never binds here.` : "The budget does not change the allocation here."} Benefit depends only on the species group, never on the method, so the methods are interchangeable in this result: which method is reported for a unit is a tie, not a finding. Budgets would start to bind only if each unit stood for a larger area, and methods could be told apart only with data that link a method to survival.</div>` : ""}
    <div class="callout"><strong>What is and isn't modelled.</strong> ${esc(ae.scope_note)}</div>
    <h2>Mean–CVaR frontier at a $45M budget</h2><div class="card">${front}<p class="small muted">The gap between the two lines is what hedging against bad scenarios costs. ${fr.price_of_robustness != null && fr.price_of_robustness < 0.005 ? `Here it is ${fmt(fr.price_of_robustness, 2)}: the same allocation is best on average and in the worst 20% of the ${ae.n_scenarios} scenarios, because viability differs little between them, so there is little downside to hedge.` : `Here it is ${fmt(fr.price_of_robustness, 2)}.`}</p></div>
    <h2>Budget sweep</h2><div class="card tablewrap"><table><thead><tr><th>Budget</th><th class="num">Expected benefit</th><th class="num">CVaR</th><th class="num">Units planted</th></tr></thead><tbody>${sw.map((r) => `<tr><td>${money(r.budget_usd)}</td><td class="num">${fmt(r.expected, 1)}</td><td class="num">${fmt(r.cvar, 1)}</td><td class="num">${r.n_units_planted}</td></tr>`).join("")}</tbody></table></div>
    <h2>Intervention methods and costs</h2><div class="card tablewrap"><table><thead><tr><th>Method</th><th class="num">Cost per hectare</th><th>Cost basis</th></tr></thead><tbody>${methods.map((o) => { const c = ct[o.intervention] || {}, b = costBasis[o.cost_basis] || ["warn", String(o.cost_basis || "").replace(/_/g, " ")]; return `<tr><td>${esc(o.intervention.replace(/_/g, " "))}</td><td class="num">${o.cost_per_ha_usd != null ? "$" + o.cost_per_ha_usd.toLocaleString("en") : "\u2014"}${c.cost_low != null ? `<br><span class="small muted">$${c.cost_low.toLocaleString("en")}\u2013$${c.cost_high.toLocaleString("en")}</span>` : ""}</td><td>${chip(b[0], b[1])}${c.source ? `<br><span class="small muted">${esc(c.source)}</span>` : ""}</td></tr>`; }).join("")}</tbody></table>
      <p class="small muted" style="margin-top:8px">Benefit = viability × $${ae.value_per_ha_year_usd}/ha/yr (national ecosystem-services value)${ae.incremental_share != null ? ` × ${Math.round(100 * ae.incremental_share)}% (the share a restored degraded hectare adds) = $${fmt(ae.net_value_per_ha_year_usd, 2)}/ha/yr at full viability` : ""}. Eligible units exclude protected areas and heavily human-modified land.</p></div></div>`;
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
      <p class="small muted" style="margin-top:8px">Boyce index under spatial block cross-validation: +1 = predictions track presences, 0 = no better than random, negative = worse. Folds are coloured accordingly.</p></div>` : "<p>MERISTEM has not been fitted.</p>"}
    <h2>MNEME — observed dieback</h2>
    ${mn ? `<div class="card"><div class="grid cols-4"><div class="stat"><div class="num">${mn.n_points_valid_kndvi ?? "—"}/${mn.n_points ?? "—"}</div><div class="cap">points with a usable satellite series</div></div><div class="stat"><div class="num">${(mn.n_person_years ?? 0).toLocaleString("en")}</div><div class="cap">person-years</div></div><div class="stat"><div class="num">${mn.n_events ?? "—"}</div><div class="cap">observed dieback events</div></div><div class="stat"><div class="num">${mn.status === "fitted" ? chip("good", "fitted") : chip("bad", "no fit")}</div><div class="cap">${esc((mn.status || "").replace(/_/g, " "))}</div></div></div>
      <p class="small muted" style="margin-top:10px">${esc(mn.scope_note || "")}</p></div>` : "<p>MNEME has not been run.</p>"}
    <h2>Does predicted viability track mapped forest?</h2>
    ${ev ? `<div class="card tablewrap"><table><thead><tr><th>Group</th><th class="num">Spearman ρ</th><th class="num">p</th><th class="num">cells</th><th class="num">mean mapped cover</th></tr></thead><tbody>${evRows}</tbody></table>
      <p class="small muted" style="margin-top:8px">REFUGIUM viability against the real Ecosystem Map of Armenia's forest classes in a ${ev.window_radius_m} m window. ${ev.n_outside} of ${ev.n_cells_total} cells fall outside Armenia's border and are excluded. ρ is undefined for pine because no mapped pine cover falls in any sampled window.</p></div>` : ""}
    <h2>Sampling design</h2>${vgHtml || "<p>No variogram has been estimated.</p>"}</div>`;
}

/* ---- Methodology ---- */
/* One card per engine, with an everyday-analogy toggle. Formulas in the technical text are written as $...$ and typeset (ui/math.js). */
function renderMethod() {
  let plain = false; try { plain = localStorage.getItem("antar_plain") === "1"; } catch (e) { /* ignore */ }
  const m = state.meth;
  const flow = `<figure class="flow" aria-label="Flow chart of the six engines"><div class="flow-svg" id="flowsvg"></div></figure>`;
  const body = () => `<div class="callout info"><strong>${plain ? "Plain-language view." : "Technical view."}</strong> ${esc(plain ? m.intro.plain : m.intro.technical)}</div>` + flow +
    m.engines.map((e) => `<section class="card engine" id="${e.id}" style="margin-bottom:14px"><h2 style="margin-top:0">${esc(e.name)}</h2><p>${chip("neutral", "status")} <span class="small">${esc(e.status)}</span></p><ul class="${plain ? "plainbox" : ""}">${(plain ? e.plain : e.technical).map((t) => `<li>${plain ? esc(t) : Tex.text(t)}</li>`).join("")}</ul></section>`).join("");
  view().innerHTML = `<div class="wrap"><h1>Methodology</h1><p class="muted">What is actually implemented, engine by engine. Switch to the plain-language view for an everyday-analogy explanation of each mechanism; the technical view keeps the equations and parameters.</p>
    <label class="toggle"><input type="checkbox" id="plain"> Explain it simply</label><div id="mbody" style="margin-top:14px"></div></div>`;
  $("#plain").checked = plain;
  const draw = () => {
    $("#mbody").innerHTML = body();
    if (!plain) Tex.load().then(() => Tex.typeset($("#mbody"))).catch(() => {});
    if (!state.flowSvg) state.flowSvg = fetch("assets/architecture.svg?v=" + window.ANTAR_BUILD).then((r) => (r.ok ? r.text() : "")).catch(() => "");
    state.flowSvg.then((svg) => { const box = $("#flowsvg"); if (box) box.innerHTML = svg || '<p class="muted small">The flow chart could not be loaded.</p>'; });
  };
  draw();
  $("#plain").addEventListener("change", (e) => { plain = e.target.checked; try { localStorage.setItem("antar_plain", plain ? "1" : "0"); } catch (x) { /* ignore */ } draw(); });
  let j = null; try { j = sessionStorage.getItem("jump"); sessionStorage.removeItem("jump"); } catch (e) { /* ignore */ }
  if (j && document.getElementById(j)) document.getElementById(j).scrollIntoView();
}

/* ---- Status & limits ---- */
function renderStatus() {
  const m = state.meth;
  view().innerHTML = `<div class="wrap"><h1>Status &amp; limits</h1><h2>Engines</h2><div class="card tablewrap"><table><thead><tr><th>Engine</th><th>State</th></tr></thead><tbody>${m.engines.map((e) => `<tr><td><a href="#/method" data-jump="${e.id}">${esc(e.name)}</a></td><td>${esc(e.status)}</td></tr>`).join("")}</tbody></table></div>
    <h2>Placeholders and assumptions</h2><div class="card tablewrap"><table><thead><tr><th>Quantity</th><th>Value used</th><th>Affects</th><th>State</th></tr></thead><tbody>${m.placeholders.map((p) => `<tr><td>${esc(p.name)}</td><td>${esc(p.value)}</td><td>${esc(p.effect)}</td><td>${esc(p.state)}</td></tr>`).join("")}</tbody></table></div>
    <h2>Known gaps</h2><div class="card"><ul>${m.known_gaps.map((g) => `<li>${esc(g)}</li>`).join("")}</ul></div></div>`;
  view().querySelectorAll("[data-jump]").forEach((a) => a.addEventListener("click", () => sessionStorageSafe("jump", a.dataset.jump)));
}

/* ---- References ---- */
function renderRefs() {
  const refs = state.M.references;
  view().innerHTML = `<div class="wrap"><h1>References</h1><input type="search" id="rq" placeholder="Filter by dataset, source or citation…" style="width:min(520px,100%)"><div id="rl" class="card" style="margin-top:12px"></div></div>`;
  const draw = (q) => {
    const f = refs.filter((r) => !q || [r.dataset, r.source, r.citation].join(" ").toLowerCase().includes(q.toLowerCase()));
    $("#rl").innerHTML = f.length ? f.map((r) => `<div class="refitem"><div class="t">${esc(r.dataset || "")}</div><div class="small">${esc(r.citation)}</div><div class="small muted">${r.url ? `<a href="${esc(r.url)}" target="_blank" rel="noopener">${esc(r.url)}</a> · ` : ""}${r.license ? "Licence: " + esc(r.license) + " · " : ""}${r.accessed ? "accessed " + esc(r.accessed) : ""}</div></div>`).join("") : `<p class="muted">No matches.</p>`;
  };
  draw(""); $("#rq").addEventListener("input", (e) => draw(e.target.value));
}

/* ---- About the author ---- */
function renderAuthor() {
  view().innerHTML = `<div class="wrap"><h1>About the author</h1>
    <div class="author">
      <figure class="author-photo"><img src="assets/author.jpg" width="720" height="720" alt="Portrait of Narek Ohanyan"></figure>
      <div class="author-bio">
        <h2>Narek Ohanyan</h2>
        <p class="author-role">Climate &amp; environmental researcher, AUA Acopian Center for the Environment</p>
        <p>Narek Ohanyan is a young climate leader from Armenia. He is a climate &amp; environmental researcher at the AUA Acopian Center for the Environment. He has also served as a Guest Scientist at the Swiss Federal Institute for Forest, Snow and Landscape Research WSL. His current research focuses on modeling forest climate resilience.</p>
        <p>Focusing on climate analytics and resilience, he authored the book <em>The Overshoot: Life After the 1.5°C Limit</em>. Across 10 years of dedication and activism, he has managed 50+ projects and secured tens of thousands of dollars in grants, empowering 10k+ Youth for Climate Action.</p>
        <p>Nationally, Narek serves as the Lead Organizer for LCOY Armenia (2025-2026). He advises the Ministry of Environment and UNICEF through the Youth Climate Council, and previously served as a UNFCCC COP27 Party Delegate and UN Youth and Children High-Level Climate Champion. Narek holds a BS in Computer Science (’26) from the American University of Armenia.</p>
        <p><a class="btn" href="https://www.narekohanyan.com" target="_blank" rel="noopener">Learn more at www.narekohanyan.com<span class="sr-only"> (opens in a new tab)</span></a></p>
      </div>
    </div></div>`;
}

/* ---- Acknowledgments ---- */
function renderAck(p) {
  const ext = (href, text) => `<a href="${href}" target="_blank" rel="noopener">${text}<span class="sr-only"> (opens in a new tab)</span></a>`;
  view().innerHTML = `<div class="wrap"><h1>Acknowledgments</h1>
    <section class="ack-rows" aria-label="Acknowledgments">
      <article class="ack-row">
        <header class="ack-label"><span class="ack-no">01</span><h2>Data, methodology and the research visit</h2></header>
        <div class="ack-text"><p>The author would like to express sincere gratitude to the Swiss Federal Institute for Forest, Snow and Landscape Research WSL, and in particular to Franziska Zilker, Tobias Kühnhanss and Dr. Michael James McCarthy of the Dynamic Macroecology group, for providing datasets and bias-corrected environmental data and for their technical feedback on the methodology. The author also thanks his supervisor, PD Dr. Marco Pütz, and the coordinator, Dr. Dominik Braunschweiger, for making his guest scientist visit at WSL possible.</p></div>
      </article>
      <article class="ack-row">
        <header class="ack-label"><span class="ack-no">02</span><h2>Conceptualization and logistics</h2></header>
        <div class="ack-text"><p>The author further acknowledges Alen Amirkhanian, Director of the AUA Acopian Center for the Environment, and the wider team of the Forest Restoration and Climate Change in Armenia (FORACCA) project for their logistical support and collaborative insights during the conceptualization of ANTAR. <a href="#/ack?s=foracca" data-scroll="foracca">About the project ↓</a></p></div>
      </article>
    </section>

    <section class="foracca" id="foracca" aria-labelledby="foracca-h">
      <div class="sec-head"><span class="sec-no">03</span><h2 id="foracca-h">About the FORACCA project</h2></div>
      <div class="foracca-grid">
        <div class="foracca-text">
          <p class="foracca-lead">Forest Restoration and Climate Change in Armenia (FORACCA) is a Swiss-funded programme that supports reforestation on community lands, climate-smart forest management and climate-resilient development in Armenia.</p>
          <h3>What it sets out to do</h3>
          <ol class="foracca-aims">
            <li><span>1</span>Advance scientific understanding of Armenia’s capacity to address climate change and sustainably manage its forests.</li>
            <li><span>2</span>Promote climate-smart practices in rural areas.</li>
            <li><span>3</span>Ensure evidence-based policymaking for climate adaptation and efficient forest management.</li>
          </ol>
          <h3>Climate services</h3>
          <p>The project provides new climate services for Armenia, including high-resolution climate scenarios and local climate impact profiles for every municipality.</p>
          <p class="small muted">Facts as published by ${ext("https://www.wsl.ch/en/projects/foracca/", "WSL")} and ${ext("https://armenpress.am/en/article/1126549", "Armenpress")}; checked 2026-10-06.</p>
        </div>
        <dl class="factsheet">
          <div><dt>Full name</dt><dd>Forest Restoration and Climate Change in Armenia (FORACCA)</dd></div>
          <div><dt>Funder</dt><dd>Swiss Agency for Development and Cooperation (SDC)</dd></div>
          <div><dt>Programme</dt><dd>10 years, 2023–2033, CHF 10 million</dd></div>
          <div><dt>Main phase</dt><dd>2025–2028</dd></div>
          <div><dt>Implemented by</dt><dd>The Forest Alliance, a consortium of Armenian NGOs led by Shen NGO; the Swiss Federal Research Institute WSL; and the Food and Agriculture Organization of the United Nations (FAO)</dd></div>
          <div><dt>More</dt><dd>${ext("https://www.wsl.ch/en/projects/foracca/", "WSL project page")}</dd></div>
        </dl>
      </div>
    </section>
    ${(state.M.photos || []).length ? `<section class="ack-rows photo-credits" aria-label="Photographs">
      <article class="ack-row">
        <header class="ack-label"><span class="ack-no">04</span><h2>Photographs</h2></header>
        <div class="ack-text">
          <p>The landscape photographs on the Overview page are third-party works, shown under Creative Commons licences. Each copy is resized for the page and is otherwise unchanged apart from cropping by the page layout. They illustrate landscape types; they are not model outputs.</p>
          <ul>${state.M.photos.map((ph) => `<li><strong>${esc(ph.caption)}</strong>, ${esc(ph.place)} (${esc(ph.date.slice(0, 4))}). Photo: ${esc(ph.author)}; ${ext(ph.licence_url, esc(ph.licence))}; ${ext(ph.source_page, esc(ph.source))}.</li>`).join("")}</ul>
        </div>
      </article>
    </section>` : ""}
    <p class="muted small ack-foot">Data providers and their licences are listed on the <a href="#/refs">References</a> page.</p></div>`;
  const go = () => $("#foracca").scrollIntoView({ behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth", block: "start" });
  view().querySelectorAll("[data-scroll]").forEach((l) => l.addEventListener("click", (e) => { e.preventDefault(); go(); }));
  if (p && p.s === "foracca") go();
}

/* ---------- footer: copy the reference ---------- */
function initCite() {
  const status = document.querySelector(".cite-status");
  document.querySelectorAll("[data-copy]").forEach((b) => b.addEventListener("click", async () => {
    const src = document.getElementById(b.dataset.copy);
    if (!src) return;
    const text = src.innerText.trim();
    let done = false;
    try { await navigator.clipboard.writeText(text); done = true; } catch (e) { /* clipboard blocked: select the text instead */ }
    if (!done) {
      const sel = getSelection(), range = document.createRange();
      range.selectNodeContents(src); sel.removeAllRanges(); sel.addRange(range);
      try { done = document.execCommand("copy"); } catch (e) { done = false; }
    }
    status.textContent = done ? (b.dataset.copy === "cite-bib" ? "BibTeX copied" : "Reference copied") : "Press Ctrl/Cmd+C to copy the selected text";
    setTimeout(() => { status.textContent = ""; }, 3000);
  }));
}

/* ---------- boot ---------- */
async function init() {
  view().innerHTML = `<div class="loading">Loading results…</div>`;
  try {
    const get = async (u) => { const r = await fetch(u + "?v=" + (window.ANTAR_BUILD || "dev")); if (!r.ok) throw new Error(u + " → HTTP " + r.status); return r.json(); };
    state.M = await get("data/manifest.json");
    state.meth = await get("data/methodology.json");
    for (const [gid, info] of Object.entries(state.M.grids)) state.G[gid] = await get("data/" + info.file);
  } catch (e) {
    view().innerHTML = `<div class="wrap"><div class="callout bad"><strong>Could not load the data.</strong> ${esc(e.message)}<br>Build it with <code>python3 ui/build_data.py</code> and serve the folder over HTTP (browsers block <code>fetch</code> on <code>file://</code>): <code>python3 -m http.server 8000 --directory ui</code>, then open <code>http://localhost:8000</code>.</div></div>`;
    return;
  }
  window.addEventListener("hashchange", route);
  window.addEventListener("scroll", syncTopbar, { passive: true });
  route();
}
function toggleTheme() {
  const cur = document.documentElement.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  const next = cur === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem("antar_theme", next); } catch (e) { /* ignore */ }
}
try { const t = localStorage.getItem("antar_theme"); if (t) document.documentElement.dataset.theme = t; } catch (e) { /* ignore */ }
document.addEventListener("DOMContentLoaded", () => { $("#theme").addEventListener("click", toggleTheme); initMenu(); initCite(); Palette.init(); init(); });
