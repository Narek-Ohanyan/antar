/* Search palette: Ctrl K / ⌘K (or the Search button, or "/") opens a box that jumps to any page, map layer or marz.
   rank() is a pure function (tested in node, tests/test_home.py); the rest is the dialog. Nothing here reads a result: it only builds links
   from the page list, data/manifest.json's map plan and the marz names in assets/map/grid.json. */
(function (root) {
  "use strict";
  const norm = (s) => String(s).toLowerCase().normalize("NFD").replace(/[̀-ͯ]/g, "");
  const escRe = (s) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  /* Every query word must occur in the item's title, subtitle or keywords; a word at the start of the title counts most, then at the start of a title word. */
  function score(item, tokens) {
    const title = norm(item.t), hay = norm(item.t + " " + (item.s || "") + " " + (item.k || ""));
    let total = 0;
    for (const tok of tokens) {
      if (!hay.includes(tok)) return -1;
      if (title.startsWith(tok)) total += 6;
      else if (new RegExp("(^|[^a-z0-9])" + escRe(tok)).test(title)) total += 4;
      else if (title.includes(tok)) total += 2;
      else total += 1;
    }
    return total;
  }
  function rank(items, query, limit) {
    const max = limit || 14, tokens = norm(query).split(/\s+/).filter(Boolean);
    if (!tokens.length) return items.slice(0, max);
    return items.map((it, i) => ({ it, i, sc: score(it, tokens) })).filter((x) => x.sc >= 0).sort((a, b) => b.sc - a.sc || a.i - b.i).slice(0, max).map((x) => x.it);
  }

  /* ---------------- dialog (browser only) ---------------- */
  let dlg = null, items = [], shown = [], active = 0, regions = null;

  async function loadItems() {
    const out = [];
    for (const [k, t] of (typeof NAV !== "undefined" ? NAV : [])) out.push({ g: "Pages", t, s: "page", href: "/" + k, k: k });
    const M = typeof state !== "undefined" ? state.M : null, plan = (M && M.map_plan) || [];
    for (const p of plan) {
      const grp = p.group && M.groups[p.group] ? M.groups[p.group].short : "";
      out.push({ g: "Map layers", t: grp && !norm(p.label).includes(norm(grp)) ? `${p.label} \u2014 ${grp}` : p.label, s: p.engine || "", href: `/map?q=${encodeURIComponent(p.id)}`, k: p.id });
    }
    if (!regions) {
      try { regions = (await (await fetch("assets/map/grid.json?v=" + (window.ANTAR_BUILD || ""))).json()).regions || []; } catch (e) { regions = []; }
    }
    for (const r of regions) out.push({ g: "Places", t: r.name === "Yerevan" ? "Yerevan (city)" : `${r.name} marz`, s: "Place explorer", href: `/site?marz=${r.id}`, k: "marz province region" });
    out.push({ g: "Actions", t: "Switch light / dark theme", s: "action", act: () => document.getElementById("theme").click(), k: "dark light mode appearance" });
    out.push({ g: "Actions", t: "How to cite ANTAR", s: "action", act: () => { const c = document.getElementById("cite"); if (c) c.scrollIntoView({ behavior: "smooth" }); }, k: "reference bibtex citation" });
    return out;
  }

  function paint() {
    const list = dlg.querySelector(".pal-list"), input = dlg.querySelector("input");
    if (!shown.length) { list.innerHTML = `<li class="pal-empty" role="presentation">No page, layer or place matches “${esc(input.value)}”.</li>`; input.removeAttribute("aria-activedescendant"); return; }
    let html = "", last = "";
    shown.forEach((it, i) => {
      if (it.g !== last) { html += `<li class="pal-group" role="presentation">${esc(it.g)}</li>`; last = it.g; }
      html += `<li class="pal-item" role="option" id="pal-o${i}" data-i="${i}" aria-selected="${i === active}"><span class="pal-t">${esc(it.t)}</span><span class="pal-s">${esc(it.s)}</span></li>`;
    });
    list.innerHTML = html;
    input.setAttribute("aria-activedescendant", "pal-o" + active);
    const cur = list.querySelector(`#pal-o${active}`); if (cur) cur.scrollIntoView({ block: "nearest" });
  }
  function refresh() { shown = rank(items, dlg.querySelector("input").value); active = 0; paint(); }
  function go(i) {
    const it = shown[i]; if (!it) return;
    dlg.close();
    if (it.act) it.act(); else if (it.href) location.hash = it.href;
  }
  function build() {
    dlg = document.createElement("dialog");
    dlg.className = "palette"; dlg.setAttribute("aria-label", "Search the site");
    dlg.innerHTML = `<div class="pal-head"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.9" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><circle cx="11" cy="11" r="6.5"/><path d="m16 16 4.5 4.5"/></svg>
      <input type="text" role="combobox" aria-expanded="true" aria-controls="pal-list" aria-autocomplete="list" aria-label="Search pages, map layers and places" placeholder="Search pages, map layers and places…" autocomplete="off" spellcheck="false">
      <button type="button" class="pal-x" aria-label="Close the search">Esc</button></div>
      <ul class="pal-list" id="pal-list" role="listbox" aria-label="Results"></ul>
      <div class="pal-foot"><span><kbd>↑</kbd> <kbd>↓</kbd> move</span><span><kbd>Enter</kbd> open</span><span><kbd>Esc</kbd> close</span></div>`;
    document.body.appendChild(dlg);
    const input = dlg.querySelector("input");
    input.addEventListener("input", refresh);
    input.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown") { e.preventDefault(); active = Math.min(shown.length - 1, active + 1); paint(); }
      else if (e.key === "ArrowUp") { e.preventDefault(); active = Math.max(0, active - 1); paint(); }
      else if (e.key === "Home") { e.preventDefault(); active = 0; paint(); }
      else if (e.key === "End") { e.preventDefault(); active = Math.max(0, shown.length - 1); paint(); }
      else if (e.key === "Enter") { e.preventDefault(); go(active); }
    });
    dlg.querySelector(".pal-list").addEventListener("click", (e) => { const li = e.target.closest(".pal-item"); if (li) go(+li.dataset.i); });
    dlg.querySelector(".pal-list").addEventListener("mousemove", (e) => { const li = e.target.closest(".pal-item"); if (li && +li.dataset.i !== active) { active = +li.dataset.i; paint(); } });
    dlg.querySelector(".pal-x").addEventListener("click", () => dlg.close());
    dlg.addEventListener("click", (e) => { if (e.target === dlg) dlg.close(); });
  }
  async function open() {
    if (!dlg) build();
    if (dlg.open) return;
    dlg.querySelector("input").value = "";
    items = await loadItems();
    refresh();
    dlg.showModal();
    dlg.querySelector("input").focus();
  }
  function init() {
    const btn = document.getElementById("search");
    if (!btn) return;
    const mac = /Mac|iPhone|iPad/i.test(navigator.platform || navigator.userAgent || "");
    const kbd = document.getElementById("search-kbd"); if (kbd) kbd.textContent = mac ? "⌘ K" : "Ctrl K";
    btn.addEventListener("click", open);
    document.addEventListener("keydown", (e) => {
      const typing = /^(INPUT|TEXTAREA|SELECT)$/.test((e.target && e.target.tagName) || "") || (e.target && e.target.isContentEditable);
      if ((e.key === "k" || e.key === "K") && (e.ctrlKey || e.metaKey)) { e.preventDefault(); open(); }
      else if (e.key === "/" && !typing && !e.ctrlKey && !e.metaKey && !e.altKey) { e.preventDefault(); open(); }
    });
  }

  const api = { rank, norm, score, init };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Palette = api;
})(typeof window !== "undefined" ? window : globalThis);
