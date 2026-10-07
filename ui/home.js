/* ANTAR home page: the headline-results band (scenario selector, trend lines, relief with the treeline contours) and the landscape photographs.

   Pure functions (no DOM) are loaded by the browser as window.Home and by node (tests/test_home.py):
     contours(z, W, H, level)         marching squares on a W x H grid; non-finite cells are skipped; returns polylines [[x, y], ...]
     downsample(z, W, H, k)           k x k block mean over the finite cells, NaN where fewer than 3/4 of the block is finite
     headlineModel(M, ssp, hz, tl)    the numbers on the cards for one emissions path and horizon, straight from the manifest
     sparkSvg(opts)                   a small trend chart as an SVG string
   Nothing here computes a model result: means and ranges across ensemble members are already in the manifest. */
(function (root) {
  "use strict";
  const fin = (v) => v != null && isFinite(v);
  const SSPS = ["ssp126", "ssp370", "ssp585"];
  const SSP_NAME = { ssp126: "SSP1-2.6", ssp370: "SSP3-7.0", ssp585: "SSP5-8.5" };
  const SSP_KEY = { ssp126: "s126", ssp370: "s370", ssp585: "s585" };
  const HORIZONS = ["2050", "2080", "2100"];
  const YEARS = [2019, 2050, 2080, 2100];
  const MINUS = "−";
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const thousands = (v) => (v < 0 ? MINUS : "") + Math.abs(Math.round(v)).toLocaleString("en");

  /* ---------------- marching squares ---------------- */
  function contours(z, W, H, level) {
    const adj = new Map(), pt = new Map();
    const eid = (vertical, x, y) => (y * W + x) * 2 + vertical;      // vertical = 0: edge (x,y)-(x+1,y); 1: edge (x,y)-(x,y+1)
    const link = (e1, e2) => {
      (adj.get(e1) || adj.set(e1, []).get(e1)).push(e2);
      (adj.get(e2) || adj.set(e2, []).get(e2)).push(e1);
    };
    const cross = (from, to) => (level - from) / (to - from);
    for (let y = 0; y < H - 1; y++) {
      for (let x = 0; x < W - 1; x++) {
        const a = z[y * W + x], b = z[y * W + x + 1], c = z[(y + 1) * W + x + 1], d = z[(y + 1) * W + x];
        if (!(fin(a) && fin(b) && fin(c) && fin(d))) continue;
        const k = (a >= level ? 8 : 0) | (b >= level ? 4 : 0) | (c >= level ? 2 : 0) | (d >= level ? 1 : 0);
        if (k === 0 || k === 15) continue;
        const T = eid(0, x, y), B = eid(0, x, y + 1), L = eid(1, x, y), R = eid(1, x + 1, y);
        if (!pt.has(T)) pt.set(T, [x + cross(a, b), y]);
        if (!pt.has(B)) pt.set(B, [x + cross(d, c), y + 1]);
        if (!pt.has(L)) pt.set(L, [x, y + cross(a, d)]);
        if (!pt.has(R)) pt.set(R, [x + 1, y + cross(b, c)]);
        const centre = (a + b + c + d) / 4 >= level;
        switch (k) {
          case 1: link(L, B); break;
          case 2: link(B, R); break;
          case 3: link(L, R); break;
          case 4: link(T, R); break;
          case 5: if (centre) { link(L, T); link(B, R); } else { link(T, R); link(L, B); } break;
          case 6: link(T, B); break;
          case 7: link(L, T); break;
          case 8: link(L, T); break;
          case 9: link(T, B); break;
          case 10: if (centre) { link(T, R); link(L, B); } else { link(L, T); link(B, R); } break;
          case 11: link(T, R); break;
          case 12: link(L, R); break;
          case 13: link(B, R); break;
          case 14: link(L, B); break;
        }
      }
    }
    const seen = new Set(), lines = [];
    const walk = (start) => {
      const line = [pt.get(start)];
      seen.add(start);
      let prev = -1, cur = start;
      for (;;) {
        const next = adj.get(cur).find((n) => n !== prev && !(seen.has(n) && n !== start));
        if (next === undefined) break;
        line.push(pt.get(next));
        if (next === start) break;                                  // the loop is closed
        seen.add(next); prev = cur; cur = next;
      }
      return line;
    };
    for (const [e, nb] of adj) if (nb.length === 1 && !seen.has(e)) lines.push(walk(e));     // open lines first (they end at missing data)
    for (const e of adj.keys()) if (!seen.has(e)) lines.push(walk(e));                         // then closed loops
    return lines;
  }

  function pathD(lines, scale, offset, nd) {
    const f = (v) => (v * scale + offset).toFixed(nd == null ? 1 : nd);
    return lines.map((l) => "M" + l.map((p) => f(p[0]) + " " + f(p[1])).join("L")).join("");
  }

  function downsample(z, W, H, k) {
    const w = Math.floor(W / k), h = Math.floor(H / k), out = new Float32Array(w * h);
    for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
      let s = 0, n = 0;
      for (let j = 0; j < k; j++) for (let i = 0; i < k; i++) { const v = z[(y * k + j) * W + x * k + i]; if (fin(v)) { s += v; n++; } }
      out[y * w + x] = n >= 0.75 * k * k ? s / n : NaN;
    }
    return { z: out, W: w, H: h };
  }

  /* ---------------- headline numbers ---------------- */
  const signed = (v, d, unit) => { const r = +v.toFixed(d); return (r > 0 ? "+" : r < 0 ? MINUS : "") + Math.abs(r).toFixed(d) + (unit || ""); };

  /* One group's numbers for the chosen path and horizon, and its trend lines for all three paths (2019, then the three horizons). */
  function groupModel(M, g, ssp, hz) {
    const base = M.baseline_viability_2019 && M.baseline_viability_2019[g];
    const ss = M.scenario_summary && M.scenario_summary[g];
    const cur = ss && ss[ssp] && ss[ssp][hz];
    const series = {}, lo = {}, hi = {};
    if (ss && fin(base)) for (const s of SSPS) {
      if (!ss[s]) continue;
      series[s] = [base].concat(HORIZONS.map((h) => (ss[s][h] ? ss[s][h].mean : null)));
      lo[s] = [base].concat(HORIZONS.map((h) => (ss[s][h] ? ss[s][h].min : null)));
      hi[s] = [base].concat(HORIZONS.map((h) => (ss[s][h] ? ss[s][h].max : null)));
    }
    return {
      id: g, short: M.groups[g].short, label: M.groups[g].label, base: fin(base) ? base : null, cur: cur || null,
      deltaPp: cur && fin(base) ? 100 * (cur.mean - base) : null, series, lo, hi,
    };
  }

  function treelineModel(M, ssp, hz) {
    const sm = M.treeline && M.treeline.summary;
    if (!sm || !sm[ssp + "__" + hz]) return null;
    const cur = sm[ssp + "__" + hz], series = {}, lo = {}, hi = {};
    for (const s of SSPS) {
      if (!HORIZONS.every((h) => sm[s + "__" + h])) continue;
      series[s] = [0].concat(HORIZONS.map((h) => sm[s + "__" + h].ensemble_mean_shift_m));
      lo[s] = [0].concat(HORIZONS.map((h) => sm[s + "__" + h].ensemble_min_shift_m));
      hi[s] = [0].concat(HORIZONS.map((h) => sm[s + "__" + h].ensemble_max_shift_m));
    }
    return { mean: cur.ensemble_mean_shift_m, min: cur.ensemble_min_shift_m, max: cur.ensemble_max_shift_m, n: cur.n_gcms, series, lo, hi, grid: M.treeline.grid };
  }

  function headlineModel(M, ssp, hz, treelineMean) {
    const groups = Object.keys(M.groups).map((g) => groupModel(M, g, ssp, hz)).filter((g) => g.base != null);
    const tl = treelineModel(M, ssp, hz);
    const withNow = groups.filter((g) => g.deltaPp != null);
    const parts = withNow.map((g) => `${g.short} ${signed(g.deltaPp, 1, " pp")}`).join(", ");
    let lead = "";
    if (tl && withNow.length) lead = `Under ${SSP_NAME[ssp]} by ${hz}, the climatic treeline rises by ${thousands(tl.mean)} m on average (models: ${thousands(tl.min)} to ${thousands(tl.max)} m), and one-year hydraulic survival changes by ${parts} against 2019.`;
    else if (tl) lead = `Under ${SSP_NAME[ssp]} by ${hz}, the climatic treeline rises by ${thousands(tl.mean)} m on average (models: ${thousands(tl.min)} to ${thousands(tl.max)} m).`;
    else if (withNow.length) lead = `Under ${SSP_NAME[ssp]} by ${hz}, one-year hydraulic survival changes by ${parts} against 2019.`;
    else lead = "Scenario results are not available yet; the 2019 values are shown.";
    const treelineNow = fin(treelineMean) ? treelineMean : null;
    const nc = M.node_counts || {}, nodes = (n) => (n ? `${n} nodes` : "");
    return { ssp, hz, groups, treeline: tl, lead, gridLabel: { viability: nodes(nc.viability), treeline: nodes(nc.treeline), aegis: nodes(nc.aegis) }, treelineNowM: treelineNow, treelineLaterM: treelineNow != null && tl ? treelineNow + tl.mean : null };
  }

  /* ---------------- small trend chart ---------------- */
  /* series: [{ssp, values: [4], selected}], band: {lo: [4], hi: [4]} for the selected one, marker: index of the chosen horizon. */
  function sparkSvg({ series, band, marker, fmtEnd, w = 240, h = 92, label = "", minSpan = 0 }) {
    const padL = 6, padR = 52, padT = 10, padB = 22;
    const all = [];
    series.forEach((s) => s.values.forEach((v) => { if (fin(v)) all.push(v); }));
    if (band) band.lo.concat(band.hi).forEach((v) => { if (fin(v)) all.push(v); });
    if (!all.length) return "";
    let lo = Math.min(...all), hi = Math.max(...all);
    if (hi - lo < Math.max(minSpan, 1e-9)) { const c = (hi + lo) / 2, m = (minSpan > 0 ? minSpan : 1) / 2; lo = c - m; hi = c + m; }     // a change smaller than minSpan is drawn as flat, as on the map
    const span = hi - lo; lo -= span * 0.12; hi += span * 0.12;
    const X = (i) => padL + ((YEARS[i] - YEARS[0]) / (YEARS[3] - YEARS[0])) * (w - padL - padR);
    const Y = (v) => padT + (1 - (v - lo) / (hi - lo)) * (h - padT - padB);
    const line = (vals) => vals.map((v, i) => (fin(v) ? [X(i), Y(v)] : null)).filter(Boolean).map((p, i) => (i ? "L" : "M") + p[0].toFixed(1) + " " + p[1].toFixed(1)).join("");
    let out = `<svg class="spark" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(label)}" focusable="false">`;
    out += `<line class="sp-axis" x1="${padL}" y1="${h - padB}" x2="${w - padR}" y2="${h - padB}"/>`;
    YEARS.forEach((yr, i) => { out += `<line class="sp-tick" x1="${X(i).toFixed(1)}" y1="${h - padB}" x2="${X(i).toFixed(1)}" y2="${h - padB + 3}"/><text class="sp-x" x="${X(i).toFixed(1)}" y="${h - 4}" text-anchor="${i === 0 ? "start" : i === 3 ? "end" : "middle"}">${yr}</text>`; });
    if (band) {
      const up = band.hi.map((v, i) => (fin(v) ? [X(i), Y(v)] : null)).filter(Boolean), dn = band.lo.map((v, i) => (fin(v) ? [X(i), Y(v)] : null)).filter(Boolean).reverse();
      if (up.length && dn.length) out += `<path class="sp-band" d="M${up.map((p) => p[0].toFixed(1) + " " + p[1].toFixed(1)).join("L")}L${dn.map((p) => p[0].toFixed(1) + " " + p[1].toFixed(1)).join("L")}Z"/>`;
    }
    series.filter((s) => !s.selected).forEach((s) => { out += `<path class="sp-line ctx ${SSP_KEY[s.ssp]}" d="${line(s.values)}"/>`; });
    series.filter((s) => s.selected).forEach((s) => {
      out += `<path class="sp-line sel ${SSP_KEY[s.ssp]}" d="${line(s.values)}"/>`;
      s.values.forEach((v, i) => { if (fin(v)) out += `<circle class="sp-dot ${i === marker ? "on" : ""} ${SSP_KEY[s.ssp]}" cx="${X(i).toFixed(1)}" cy="${Y(v).toFixed(1)}" r="${i === marker ? 4.2 : 2.4}"/>`; });
      const last = s.values[marker];
      if (fin(last) && fmtEnd) out += `<text class="sp-end ${SSP_KEY[s.ssp]}" x="${(X(marker) + 8).toFixed(1)}" y="${(Y(last) + 4).toFixed(1)}">${esc(fmtEnd(last))}</text>`;
    });
    return out + "</svg>";
  }

  function trendLabel(name, ssp, series, fmtVal) {
    const v = series[ssp];
    return v ? `${name}, ${SSP_NAME[ssp]}: ${YEARS.map((y, i) => (fin(v[i]) ? `${y} ${fmtVal(v[i])}` : null)).filter(Boolean).join(", ")}.` : name;
  }

  /* Every number behind the headline cards, for all emissions paths and horizons, as CSV (RFC 4180: comma, double quotes, CRLF). */
  function headlineCsv(M) {
    const q = (v) => (/[",\r\n]/.test(String(v)) ? '"' + String(v).replace(/"/g, '""') + '"' : String(v));
    const rows = [["quantity", "group", "emissions_path", "horizon", "baseline_2019", "ensemble_mean", "model_min", "model_max", "n_climate_models", "unit", "grid"]];
    const grid = M.treeline && M.treeline.grid ? M.treeline.grid : "";
    for (const g of Object.keys(M.groups)) {
      const base = M.baseline_viability_2019 && M.baseline_viability_2019[g], ss = M.scenario_summary && M.scenario_summary[g];
      if (!ss || !fin(base)) continue;
      for (const s of SSPS) for (const h of HORIZONS) {
        const c = ss[s] && ss[s][h];
        if (c) rows.push(["one-year hydraulic survival", M.groups[g].short, SSP_NAME[s], h, base.toFixed(4), c.mean.toFixed(4), c.min.toFixed(4), c.max.toFixed(4), c.n_gcms, "fraction", grid]);
      }
    }
    if (M.treeline && M.treeline.summary) for (const s of SSPS) for (const h of HORIZONS) {
      const c = M.treeline.summary[s + "__" + h];
      if (c) rows.push(["climatic treeline shift", "all", SSP_NAME[s], h, "0", c.ensemble_mean_shift_m.toFixed(1), c.ensemble_min_shift_m.toFixed(1), c.ensemble_max_shift_m.toFixed(1), c.n_gcms, "m", grid]);
    }
    return rows.map((r) => r.map(q).join(",")).join("\r\n") + "\r\n";
  }

  /* ---------------- cards ---------------- */
  function viabilityCard(g, m, hzIndex) {
    const pct = (v) => (100 * v).toFixed(1);
    const sel = g.series[m.ssp];
    const spark = sel ? sparkSvg({
      series: SSPS.filter((s) => g.series[s]).map((s) => ({ ssp: s, values: g.series[s], selected: s === m.ssp })),
      band: { lo: g.lo[m.ssp], hi: g.hi[m.ssp] }, marker: hzIndex + 1, fmtEnd: (v) => pct(v) + "%", minSpan: 0.04,
      label: trendLabel(`${g.short} one-year hydraulic survival`, m.ssp, g.series, (v) => pct(v) + "%"),
    }) : "";
    const dir = g.deltaPp == null ? "" : g.deltaPp < -0.05 ? "down" : g.deltaPp > 0.05 ? "up" : "flat";
    const arrow = dir === "down" ? "▼" : dir === "up" ? "▲" : "▬";
    const mapLink = `#/map?q=viab_${g.id}&mode=scen&ssp=${m.ssp}&hz=${m.hz}&gcm=ens`;
    return `<article class="fcard">
      <header><span class="eyebrow">XYLEM \u00B7 REFUGIUM</span>${m.gridLabel && m.gridLabel.viability ? `<span class="cgrid">${esc(m.gridLabel.viability)}</span>` : ""}</header>
      <h3>${esc(g.short)}</h3>
      <p class="fnum" aria-label="${g.cur ? pct(g.cur.mean) + " percent" : pct(g.base) + " percent"}">${g.cur ? pct(g.cur.mean) : pct(g.base)}<small>%</small></p>
      <p class="fwhat">${g.cur ? `mean one-year hydraulic survival, ${SSP_NAME[m.ssp]} in ${m.hz}` : "mean one-year hydraulic survival, 2019"}</p>
      ${g.cur ? `<p class="fdelta ${dir}"><span aria-hidden="true">${arrow}</span> <b>${signed(g.deltaPp, 1, " pp")}</b> against 2019 (${pct(g.base)}%)</p>` : ""}
      ${spark}
      ${g.cur ? `<p class="frange">Range across the ${g.cur.n_gcms} climate models: <b>${pct(g.cur.min)}–${pct(g.cur.max)}%</b></p>` : ""}
      <a class="flink" href="${mapLink}">See it on the map <span aria-hidden="true">→</span></a>
    </article>`;
  }

  function treelineCard(tl, m, hzIndex) {
    const sel = tl.series[m.ssp];
    const spark = sel ? sparkSvg({
      series: SSPS.filter((s) => tl.series[s]).map((s) => ({ ssp: s, values: tl.series[s], selected: s === m.ssp })),
      band: { lo: tl.lo[m.ssp], hi: tl.hi[m.ssp] }, marker: hzIndex + 1, fmtEnd: (v) => signed(v, 0, " m"), minSpan: 120,
      label: trendLabel("Climatic treeline shift", m.ssp, tl.series, (v) => signed(v, 0, " m")),
    }) : "";
    return `<article class="fcard tl">
      <header><span class="eyebrow">MERISTEM \u00B7 TREELINE</span>${m.gridLabel && m.gridLabel.treeline ? `<span class="cgrid">${esc(m.gridLabel.treeline)}</span>` : ""}</header>
      <h3>Climatic treeline</h3>
      <p class="fnum" aria-label="${signed(tl.mean, 0, " metres")}">${signed(tl.mean, 0)}<small> m</small></p>
      <p class="fwhat">uphill shift of the potential treeline, ${SSP_NAME[m.ssp]} in ${m.hz}</p>
      <p class="fdelta flat"><span aria-hidden="true">▲</span> a climatic ceiling, not a forecast of where forest will stand</p>
      ${spark}
      <p class="frange">Range across the ${tl.n} climate models: <b>${signed(tl.min, 0)} to ${signed(tl.max, 0, " m")}</b></p>
      <a class="flink" href="#/treeline">See treeline change <span aria-hidden="true">→</span></a>
    </article>`;
  }

  /* The price of robustness is shown as a number only when it is larger than the solver can resolve (the Decision page uses the same rule). */
  function priceLine(f) {
    if (!f || !fin(f.price_of_robustness)) return "";
    const resolved = fin(f.resolution_usd) ? Math.abs(f.price_of_robustness) > f.resolution_usd : Math.abs(f.price_of_robustness) >= 0.005;
    return `<p class="fdelta flat"><span aria-hidden="true">\u25AC</span> price of robustness ${resolved ? `<b>${thousands(f.price_of_robustness)} $/yr</b>` : "<b>\u2248 0</b> (below what the solver can resolve)"}</p>`;
  }

  function aegisCard(ae, m) {
    if (!ae) return "";
    return `<article class="fcard">
      <header><span class="eyebrow">AEGIS · DECISION</span>${m && m.gridLabel && m.gridLabel.aegis ? `<span class="cgrid">${esc(m.gridLabel.aegis)}</span>` : ""}</header>
      <h3>Planting portfolio</h3>
      <p class="fnum">${esc(String(ae.n_units))}</p>
      <p class="fwhat">candidate planting units${ae.n_eligible_units != null ? `, ${esc(String(ae.n_eligible_units))} eligible,` : ""} evaluated against ${esc(String(ae.n_scenarios))} climate scenarios</p>
      ${priceLine(ae.frontier)}
      <a class="flink" href="#/decision">See the portfolio <span aria-hidden="true">→</span></a>
    </article>`;
  }

  /* ---------------- landscape photographs ---------------- */
  function photoFigure(p, i) {
    const src = (w) => `assets/photos/${p.file}-${w}.jpg`;
    return `<figure class="photo p${i + 1}">
      <button type="button" class="photo-open" data-photo="${esc(p.file)}" aria-label="Enlarge the photograph: ${esc(p.caption)}, ${esc(p.place)}">
        <img src="${src(800)}" srcset="${src(800)} 800w, ${src(1600)} 1600w" sizes="(min-width: 900px) 40vw, 92vw" width="${p.original_px[0]}" height="${p.original_px[1]}" alt="${esc(p.alt)}" loading="lazy" decoding="async">
      </button>
      <figcaption><strong>${esc(p.caption)}</strong><span>${esc(p.place)}</span>${creditLine(p)}</figcaption>
    </figure>`;
  }
  function creditLine(p) {
    return `<span class="credit">Photo: ${esc(p.author)} · <a href="${esc(p.licence_url)}" target="_blank" rel="noopener">${esc(p.licence)}</a> · <a href="${esc(p.source_page)}" target="_blank" rel="noopener">${esc(p.source)}<span class="sr-only"> (opens in a new tab)</span></a> · resized</span>`;
  }
  function photoZoomBody(p) {
    return `<figure class="photo-big"><img src="assets/photos/${p.file}-1600.jpg" width="${p.original_px[0]}" height="${p.original_px[1]}" alt="${esc(p.alt)}">
      <figcaption>Taken ${esc(p.date.slice(0, 4))}. ${creditLine(p)} The copy shown is resized; the original is on the source page.</figcaption></figure>`;
  }

  /* ---------------- the band: markup and behaviour (browser only) ---------------- */
  const sel = (r, q) => r.querySelector(q);
  let reliefPromise = null;
  function loadRelief() {
    if (reliefPromise) return reliefPromise;
    const build = window.ANTAR_BUILD || "";
    reliefPromise = Promise.all([fetch("assets/relief.svg?v=" + build).then((r) => (r.ok ? r.text() : null)).catch(() => null), Raster.load("elevation", 2), Raster.load("region", 1)]).then(([svg, e, r]) => {
      if (!svg) return null;
      const W = e.w, H = e.h, ed = e.data, rd = r.data, z = new Float32Array(W * H);
      for (let i = 0; i < W * H; i++) z[i] = rd[i * 4] > 0 ? ed[i * 4] * 256 + ed[i * 4 + 1] : NaN;      // metres inside the country, NaN outside
      return { svg, ds: downsample(z, W, H, 2), k: 2 };
    }).catch(() => null);                                                        // a missing or altered raster hides the figure; it never draws wrong contours
    return reliefPromise;
  }

  function bandMarkup(M, gridHtml, start) {
    const radio = (name, val, text, checked) => `<label class="opt"><input type="radio" name="${name}" value="${val}"${checked ? " checked" : ""}><span>${text}</span></label>`;
    return `<div class="wrap band-grid">
      <div class="band-intro">
        <div class="sec-head on-band"><span class="sec-no">02</span><h2 id="res-h">Headline results</h2></div>
        <p class="band-lead" id="res-lead" aria-live="polite"></p>
        <div class="scen" role="group" aria-label="Choose a scenario">
          <div class="scen-row"><span class="scen-label" id="l-ssp">Emissions path</span><div class="seg-on-band" role="radiogroup" aria-labelledby="l-ssp">${SSPS.map((s) => radio("ssp", s, SSP_NAME[s], s === start.ssp)).join("")}</div></div>
          <div class="scen-row"><span class="scen-label" id="l-hz">Horizon</span><div class="seg-on-band" role="radiogroup" aria-labelledby="l-hz">${HORIZONS.map((h) => radio("hz", h, h, h === start.hz)).join("")}</div></div>
        </div>
        <figure class="relief" id="relief">
          <div class="relief-art" id="relief-art"></div>
          <figcaption id="relief-cap"></figcaption>
        </figure>
      </div>
      <div class="band-main">
        <div class="fcards" id="fcards"></div>
        <div class="band-foot">
          <p class="pathkey" aria-hidden="true"><span><i class="s126"></i>SSP1-2.6</span><span><i class="s370"></i>SSP3-7.0</span><span><i class="s585"></i>SSP5-8.5</span><span>solid line and shading: the chosen path and its model range; dashed: the other two</span></p>
          <p class="small">${gridHtml} Each card names the number of model nodes behind it. Means and ranges are taken across 5 climate models for the chosen emissions path and horizon; the 2019 value is the model's own baseline year.</p>
          <p class="band-dl"><button type="button" class="dl-btn" id="dl-csv"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M12 4v11m0 0-4-4m4 4 4-4M5 19h14"/></svg>Download these numbers (CSV)</button><span class="small"> all three paths and all horizons, with the model range</span></p>
          <details class="howto"><summary>How to read these numbers</summary>
            <dl>
              <dt>One-year hydraulic survival</dt><dd>One minus the modelled probability that a tree suffers hydraulic failure (xylem embolism) in a single year, averaged over the model nodes and over 50 draws of the uncertain trait values. It is the first, reduced version of viability: growth, frost and other mortality causes are not in it, and it is not a forecast of forest cover.</dd>
              <dt>SSP1-2.6, SSP3-7.0, SSP5-8.5</dt><dd>Low, intermediate and high emissions paths (Shared Socioeconomic Pathways) from the ISIMIP3b climate-model ensemble.</dd>
              <dt>Range across climate models</dt><dd>The lowest and highest of the five climate models (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL); it shows model disagreement, not a confidence interval.</dd>
              <dt>Climatic treeline</dt><dd>The elevation where the growing-season mean temperature reaches the global 6.45 °C threshold of Körner and Paulsen. It marks what temperature would allow, not where forest will stand.</dd>
              <dt>The relief</dt><dd>Contours every 250 m from the terrain model the results are computed on. The two highlighted lines mark the mean potential treeline elevation across the model nodes today and after the chosen shift; the real treeline differs from place to place.</dd>
            </dl>
          </details>
        </div>
      </div>
    </div>`;
  }

  function mountBand(host, { M, treelineMean, gridHtml, start, onChange }) {
    const first = { ssp: SSPS.includes(start && start.ssp) ? start.ssp : "ssp585", hz: HORIZONS.includes(start && start.hz) ? start.hz : "2100" };
    host.innerHTML = bandMarkup(M, gridHtml, first);
    const cards = sel(host, "#fcards"), lead = sel(host, "#res-lead"), fig = sel(host, "#relief"), art = sel(host, "#relief-art"), cap = sel(host, "#relief-cap");
    const hlCache = new Map();
    let relief = null, hl0 = null, hl1 = null, token = 0;
    const current = () => ({ ssp: host.querySelector("input[name=ssp]:checked").value, hz: host.querySelector("input[name=hz]:checked").value });
    const hlPath = (level) => {
      if (!hlCache.has(level)) hlCache.set(level, pathD(contours(relief.ds.z, relief.ds.W, relief.ds.H, level), relief.k, relief.k / 2, 1));
      return hlCache.get(level);
    };
    function draw() {
      const { ssp, hz } = current();
      const m = headlineModel(M, ssp, hz, treelineMean);
      const idx = HORIZONS.indexOf(hz);
      cards.innerHTML = m.groups.map((g) => viabilityCard(g, m, idx)).join("") + (m.treeline ? treelineCard(m.treeline, m, idx) : "") + aegisCard(M.aegis, m);
      cards.classList.remove("tick"); void cards.offsetWidth; cards.classList.add("tick");
      lead.textContent = m.lead;
      if (relief && m.treelineNowM != null) {
        hl0.setAttribute("d", hlPath(m.treelineNowM));
        if (m.treelineLaterM != null) { hl1.setAttribute("d", hlPath(m.treelineLaterM)); hl1.setAttribute("class", "hl1 " + SSP_KEY[ssp]); hl1.removeAttribute("hidden"); } else hl1.setAttribute("hidden", "");
        cap.innerHTML = `<span class="key k0"></span> Mean potential treeline, 2019 climate: <b>${thousands(m.treelineNowM)} m</b>` +
          (m.treelineLaterM != null ? `<br><span class="key k1 ${SSP_KEY[ssp]}"></span> After the ${SSP_NAME[ssp]} shift to ${hz}: <b>${thousands(m.treelineLaterM)} m</b>` : "") +
          `<br><span class="small muted">Contours every 250 m; the figure is the relief the models run on, not a model result.</span>`;
      }
    }
    host.addEventListener("change", (e) => { if (e.target.name === "ssp" || e.target.name === "hz") { draw(); if (onChange) onChange(current()); } });
    const dl = sel(host, "#dl-csv");
    if (dl) dl.addEventListener("click", () => {
      const url = URL.createObjectURL(new Blob([headlineCsv(M)], { type: "text/csv;charset=utf-8" }));
      const a = document.createElement("a"); a.href = url; a.download = "antar_headline_results.csv"; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 2000);
    });
    draw();
    const mine = ++token;
    loadRelief().then((r) => {
      if (mine !== token || !host.isConnected) return;
      if (!r || treelineMean == null) { fig.hidden = true; return; }
      relief = r;
      art.innerHTML = r.svg;
      const svg = art.querySelector("svg");
      if (!svg) return;
      svg.setAttribute("role", "img"); svg.setAttribute("aria-label", "Relief of Armenia with contour lines every 250 metres and the mean potential treeline elevation highlighted");
      const ns = "http://www.w3.org/2000/svg";
      hl0 = document.createElementNS(ns, "path"); hl0.setAttribute("class", "hl0");
      hl1 = document.createElementNS(ns, "path"); hl1.setAttribute("class", "hl1");
      svg.append(hl0, hl1);
      draw();
    });
  }

  function mountPhotos(host, photos) {
    host.querySelectorAll(".photo-open").forEach((btn) => btn.addEventListener("click", () => {
      const p = photos.find((x) => x.file === btn.dataset.photo);
      if (p && typeof Charts !== "undefined") Charts.zoom({ title: `${p.caption}, ${p.place}`, body: photoZoomBody(p), opener: btn });
    }));
  }

  const api = { headlineCsv, mountBand, mountPhotos, contours, pathD, downsample, groupModel, treelineModel, headlineModel, sparkSvg, viabilityCard, treelineCard, aegisCard, photoFigure, photoZoomBody, signed, SSPS, SSP_NAME, HORIZONS, YEARS };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Home = api;
})(typeof window !== "undefined" ? window : globalThis);
