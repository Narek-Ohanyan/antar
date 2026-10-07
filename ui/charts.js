/* Minimal dependency-free SVG charts. Every function returns an SVG string sized by viewBox, so it
   scales to its container. Colours come from CSS custom properties (see styles.css) so charts follow
   the light/dark theme. */
const Charts = (() => {
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  function niceTicks(lo, hi, n = 5) {
    if (!isFinite(lo) || !isFinite(hi)) return [0, 1];
    if (lo === hi) { lo -= 1; hi += 1; }
    const span = hi - lo, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || raw;
    const start = Math.ceil(lo / step) * step, out = [];
    for (let v = start; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  const fmtTick = (v) => (Math.abs(v) >= 1000 ? Math.round(v).toLocaleString("en") : +v.toPrecision(4) + "");

  function frame(w, h, m) { return { w, h, m, iw: w - m.l - m.r, ih: h - m.t - m.b }; }

  function axes(f, xs, ys, xlabel, ylabel, xfmt, yfmt) {
    let g = "";
    for (const t of ys) {
      const y = f.m.t + f.ih - ((t - ys._lo) / (ys._hi - ys._lo)) * f.ih;
      g += `<line class="grid" x1="${f.m.l}" x2="${f.m.l + f.iw}" y1="${y}" y2="${y}"/>` +
           `<text class="tick" x="${f.m.l - 8}" y="${y + 4}" text-anchor="end">${esc((yfmt || fmtTick)(t))}</text>`;
    }
    for (const t of xs) {
      const x = f.m.l + ((t - xs._lo) / (xs._hi - xs._lo)) * f.iw;
      g += `<line class="axis" x1="${x}" x2="${x}" y1="${f.m.t + f.ih}" y2="${f.m.t + f.ih + 5}"/>` +
           `<text class="tick" x="${x}" y="${f.m.t + f.ih + 20}" text-anchor="middle">${esc((xfmt || fmtTick)(t))}</text>`;
    }
    g += `<line class="axis" x1="${f.m.l}" x2="${f.m.l + f.iw}" y1="${f.m.t + f.ih}" y2="${f.m.t + f.ih}"/>`;
    if (xlabel) g += `<text class="label" x="${f.m.l + f.iw / 2}" y="${f.h - 6}" text-anchor="middle">${esc(xlabel)}</text>`;
    if (ylabel) g += `<text class="label" transform="translate(14 ${f.m.t + f.ih / 2}) rotate(-90)" text-anchor="middle">${esc(ylabel)}</text>`;
    return g;
  }

  function legend(items, x, y, k = 1) {
    let g = "", dx = 0;
    for (const it of items) {
      g += `<g transform="translate(${x + dx} ${y})"><line x1="0" x2="18" y1="0" y2="0" stroke="${it.color}" stroke-width="3" ${it.dash ? `stroke-dasharray="${it.dash}"` : ""}/>` +
           `<text class="tick" x="24" y="4">${esc(it.name)}</text></g>`;
      dx += (34 + it.name.length * 6.4) * k;
    }
    return g;
  }

  /* series: [{name,color,dash,points:[[x,y,lo?,hi?]], marker?}]  */
  function line({ series, xlabel, ylabel, xticks, yfmt, xfmt, ymin, ymax, width = 640, height = 340, legendOn = true, big = false }) {
    const f = frame(width, height, { t: legendOn ? 34 : 14, r: 16, b: 44, l: big ? (width < 560 ? 66 : 78) : 64 });
    const all = series.flatMap((s) => s.points);
    const xv = all.map((p) => p[0]);
    const yv = all.flatMap((p) => [p[1], p[2] ?? p[1], p[3] ?? p[1]]).filter((v) => v != null && isFinite(v));
    if (!yv.length) return `<p class="muted">No data to chart.</p>`;
    let lo = ymin ?? Math.min(...yv), hi = ymax ?? Math.max(...yv);
    if (ymin == null || ymax == null) { const pad = (hi - lo || Math.abs(hi) || 1) * 0.08; if (ymin == null) lo -= pad; if (ymax == null) hi += pad; }
    const ys = niceTicks(lo, hi, 5); ys._lo = Math.min(lo, ys[0]); ys._hi = Math.max(hi, ys[ys.length - 1]);
    const xs = xticks ? Object.assign([...xticks], { _lo: Math.min(...xticks), _hi: Math.max(...xticks) })
                      : Object.assign(niceTicks(Math.min(...xv), Math.max(...xv)), { _lo: Math.min(...xv), _hi: Math.max(...xv) });
    const X = (v) => f.m.l + ((v - xs._lo) / (xs._hi - xs._lo)) * f.iw;
    const Y = (v) => f.m.t + f.ih - ((v - ys._lo) / (ys._hi - ys._lo)) * f.ih;
    let body = axes(f, xs, ys, xlabel, ylabel, xfmt, yfmt);
    for (const s of series) {
      const band = s.points.filter((p) => p[2] != null && p[3] != null);
      if (band.length > 1) {
        const up = band.map((p) => `${X(p[0])},${Y(p[3])}`), dn = band.slice().reverse().map((p) => `${X(p[0])},${Y(p[2])}`);
        body += `<polygon points="${up.concat(dn).join(" ")}" fill="${s.color}" opacity="0.16"/>`;
      }
      const pts = s.points.filter((p) => p[1] != null);
      body += `<polyline fill="none" stroke="${s.color}" stroke-width="${big ? 3.4 : 2.4}" stroke-linejoin="round" ${s.dash ? `stroke-dasharray="${s.dash}"` : ""} points="${pts.map((p) => `${X(p[0])},${Y(p[1])}`).join(" ")}"/>`;
      for (const p of pts) body += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="${big ? 5.2 : 3.4}" fill="${s.color}"><title>${esc(s.name)}: ${esc((yfmt || fmtTick)(p[1]))}${p[2] != null ? ` (range ${esc((yfmt || fmtTick)(p[2]))} – ${esc((yfmt || fmtTick)(p[3]))})` : ""}</title></circle>`;
    }
    if (legendOn) body += legend(series.map((s) => ({ name: s.name, color: s.color, dash: s.dash })), f.m.l, big ? 16 : 14, big ? (width < 560 ? 1.05 : 1.32) : 1);
    return `<svg class="chart${big ? " big" : ""}" viewBox="0 0 ${width} ${height}" role="img">${body}</svg>`;
  }

  /* items: [{label, value, lo?, hi?, color?, note?}]  horizontal bars from a zero line */
  function bars({ items, xlabel, xfmt, width = 640, rowH = 34, xmin, xmax }) {
    const m = { t: 10, r: 24, b: 40, l: 180 };
    const height = m.t + m.b + items.length * rowH, f = frame(width, height, m);
    const vals = items.flatMap((i) => [i.value, i.lo ?? i.value, i.hi ?? i.value, 0]).filter((v) => v != null);
    let lo = xmin ?? Math.min(...vals), hi = xmax ?? Math.max(...vals);
    const ticks = niceTicks(lo, hi, 5); const xs = Object.assign(ticks, { _lo: Math.min(lo, ticks[0]), _hi: Math.max(hi, ticks[ticks.length - 1]) });
    const X = (v) => f.m.l + ((v - xs._lo) / (xs._hi - xs._lo)) * f.iw;
    let g = "";
    for (const t of xs) g += `<line class="grid" x1="${X(t)}" x2="${X(t)}" y1="${m.t}" y2="${m.t + items.length * rowH}"/><text class="tick" x="${X(t)}" y="${height - 18}" text-anchor="middle">${esc((xfmt || fmtTick)(t))}</text>`;
    g += `<line class="axis" x1="${X(0)}" x2="${X(0)}" y1="${m.t}" y2="${m.t + items.length * rowH}"/>`;
    items.forEach((it, i) => {
      const y = m.t + i * rowH + 6, h = rowH - 14;
      g += `<text class="tick" x="${m.l - 8}" y="${y + h / 2 + 4}" text-anchor="end">${esc(it.label)}</text>`;
      if (it.value != null) {
        const x0 = X(0), x1 = X(it.value);
        g += `<rect x="${Math.min(x0, x1)}" y="${y}" width="${Math.abs(x1 - x0)}" height="${h}" rx="3" fill="${it.color || "var(--accent)"}" opacity="0.85"><title>${esc(it.label)}: ${esc((xfmt || fmtTick)(it.value))}${it.note ? " — " + esc(it.note) : ""}</title></rect>`;
        if (it.lo != null && it.hi != null) g += `<line class="whisker" x1="${X(it.lo)}" x2="${X(it.hi)}" y1="${y + h / 2}" y2="${y + h / 2}"/>`;
      } else g += `<text class="tick" x="${X(0) + 8}" y="${y + h / 2 + 4}">${esc(it.note || "no value")}</text>`;
    });
    if (xlabel) g += `<text class="label" x="${f.m.l + f.iw / 2}" y="${height - 2}" text-anchor="middle">${esc(xlabel)}</text>`;
    return `<svg class="chart" viewBox="0 0 ${width} ${height}" role="img">${g}</svg>`;
  }

  /* points: [[x,y]] drawn as dots; curves: [{name,color,dash,fn:(x)=>y,x0,x1}] */
  function scatter({ points, curves = [], xlabel, ylabel, width = 640, height = 340, pointLabel }) {
    const f = frame(width, height, { t: 34, r: 16, b: 44, l: 70 });
    const xmax = Math.max(...points.map((p) => p[0]), ...curves.map((c) => c.x1 || 0));
    const sample = curves.flatMap((c) => Array.from({ length: 41 }, (_, i) => c.fn((c.x0 || 0) + ((c.x1 - (c.x0 || 0)) * i) / 40)));
    const yv = points.map((p) => p[1]).concat(sample).filter(isFinite);
    const ylo = Math.min(0, ...yv), yhi = Math.max(...yv) * 1.06;
    const ys = niceTicks(ylo, yhi, 5); ys._lo = Math.min(ylo, ys[0]); ys._hi = Math.max(yhi, ys[ys.length - 1]);
    const xs = niceTicks(0, xmax, 6); xs._lo = 0; xs._hi = Math.max(xmax, xs[xs.length - 1]);
    const X = (v) => f.m.l + ((v - xs._lo) / (xs._hi - xs._lo)) * f.iw, Y = (v) => f.m.t + f.ih - ((v - ys._lo) / (ys._hi - ys._lo)) * f.ih;
    let g = axes(f, xs, ys, xlabel, ylabel);
    for (const c of curves) {
      const pts = Array.from({ length: 81 }, (_, i) => { const x = (c.x0 || 0) + ((c.x1 - (c.x0 || 0)) * i) / 80; return `${X(x)},${Y(c.fn(x))}`; });
      g += `<polyline fill="none" stroke="${c.color}" stroke-width="2.2" ${c.dash ? `stroke-dasharray="${c.dash}"` : ""} points="${pts.join(" ")}"/>`;
    }
    for (const p of points) g += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="4.2" fill="var(--text)" opacity="0.85"><title>${esc(pointLabel ? pointLabel(p) : p.join(", "))}</title></circle>`;
    g += legend([{ name: "empirical", color: "var(--text)" }].concat(curves.map((c) => ({ name: c.name, color: c.color, dash: c.dash }))), f.m.l, 14);
    return `<svg class="chart" viewBox="0 0 ${width} ${height}" role="img">${g}</svg>`;
  }

  /* Enlarged view of a chart in a modal dialog (native <dialog>: Escape closes it, focus is trapped and returns to the opener,
     a click outside or on the close button closes it). `body` is trusted HTML built by the caller from its own data. */
  function closeZoom(dlg) {
    const o = dlg._opener;
    if (dlg.open) dlg.close();
    if (o && o.isConnected) o.focus({ preventScroll: true });      // Escape gets this from the browser; the close button and the backdrop get it here
  }

  function zoom({ title, body, opener }) {
    let dlg = document.getElementById("zoomdlg");
    if (!dlg) {
      dlg = document.createElement("dialog");
      dlg.id = "zoomdlg"; dlg.className = "zoom"; dlg.setAttribute("aria-labelledby", "zoomttl");
      document.body.appendChild(dlg);
      dlg.addEventListener("click", (e) => { if (e.target === dlg) closeZoom(dlg); });              // a click on the backdrop
    }
    dlg._opener = opener || document.activeElement;
    dlg.innerHTML = `<div class="zoom-head"><h3 id="zoomttl">${esc(title)}</h3><button type="button" class="zoom-x" aria-label="Close the enlarged chart">&times;</button></div><div class="zoom-body">${body}</div>`;
    dlg.querySelector(".zoom-x").addEventListener("click", () => closeZoom(dlg));
    if (typeof dlg.showModal === "function") dlg.showModal(); else dlg.setAttribute("open", "");
    dlg.scrollTop = 0;
  }

  return { line, bars, scatter, niceTicks, esc, zoom };
})();
