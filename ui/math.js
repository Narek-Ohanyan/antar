/* Typesets the formulas written as $...$ in the engine descriptions (data/methodology.json) with KaTeX (MIT, self-hosted
   in assets/katex, loaded only when a page that has formulas is opened). Nothing here computes a result. */
(function (root) {
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  let ready = null;
  const build = () => window.ANTAR_BUILD || "";

  function load() {
    if (ready) return ready;
    ready = new Promise((resolve, reject) => {
      const css = document.createElement("link");
      css.rel = "stylesheet"; css.href = "assets/katex/katex.min.css?v=" + build();
      document.head.appendChild(css);
      const s = document.createElement("script");
      s.src = "assets/katex/katex.min.js?v=" + build();
      s.onload = () => resolve();
      s.onerror = () => reject(new Error("KaTeX could not be loaded"));
      document.head.appendChild(s);
    });
    return ready;
  }

  /* Escapes a plain string and wraps each "$...$" segment as a formula element for typeset(). */
  function text(s) {
    return String(s).split(/(\$[^$]+\$)/).map((p) => (p.length > 2 && p[0] === "$" && p.endsWith("$") ? `<span class="tex">${esc(p.slice(1, -1))}</span>` : esc(p))).join("");
  }

  /* Renders every unrendered formula element under rootEl; a formula KaTeX cannot parse stays visible as its LaTeX source. */
  function typeset(rootEl) {
    rootEl.querySelectorAll(".tex:not(.done)").forEach((el) => {
      const tex = el.textContent;
      try { window.katex.render(tex, el, { displayMode: el.dataset.display === "1", throwOnError: false, strict: "ignore", trust: false }); }
      catch (e) { el.textContent = tex; el.classList.add("tex-fail"); }
      el.classList.add("done");
    });
  }

  root.Tex = { load, text, typeset };
})(window);
