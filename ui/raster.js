/* Loads the map's data rasters (region ids, forest cover, elevation) exactly as they were built, and refuses them if anything altered them on the way.

   Why: these files are data, not pictures. A host or CDN that "optimises" images can swap a PNG for a lossy WebP when the browser says it supports it, and
   the values (a region id, 256*R+G metres) are then wrong everywhere. So the files are lossless PNGs stored under a .bin name, fetched as bytes, decoded without
   colour conversion, and compared with the Adler-32 checksum written at build time (assets/map/checksums.json, scripts/stamp_map_assets.py). A mismatch throws.
   adler32() is pure and is tested against zlib (tests/test_raster.py). */
(function (root) {
  "use strict";
  const build = () => (typeof window !== "undefined" && window.ANTAR_BUILD) || "dev";

  /* Adler-32 over the first `channels` bytes of every RGBA pixel of `d`, in order (the same bytes scripts/stamp_map_assets.py hashes). */
  function adler32(d, channels) {
    let a = 1, b = 0;
    const n = d.length / 4;
    for (let i = 0; i < n; i++) {
      for (let c = 0; c < channels; c++) {
        a += d[i * 4 + c]; if (a >= 65521) a -= 65521;
        b += a; if (b >= 65521) b -= 65521;
      }
    }
    return ((b << 16) | a) >>> 0;
  }

  let sums = null;
  function checksums() {
    if (!sums) sums = fetch("assets/map/checksums.json?v=" + build()).then((r) => { if (!r.ok) throw new Error("could not load the raster checksums (HTTP " + r.status + ")"); return r.json(); });
    return sums;
  }

  async function decode(url) {
    const r = await fetch(url);
    if (!r.ok) { const e = new Error("could not load " + url + " (HTTP " + r.status + ")"); e.missing = r.status === 404; throw e; }
    const blob = new Blob([await r.arrayBuffer()], { type: "image/png" });
    let source = null, w = 0, h = 0;
    if (typeof createImageBitmap === "function") {
      try { source = await createImageBitmap(blob, { colorSpaceConversion: "none", premultiplyAlpha: "none" }); w = source.width; h = source.height; } catch (e) { source = null; }
    }
    if (!source) {
      const u = URL.createObjectURL(blob);
      try {
        source = await new Promise((res, rej) => { const im = new Image(); im.onload = () => res(im); im.onerror = () => rej(new Error("could not decode " + url)); im.src = u; });
        w = source.naturalWidth; h = source.naturalHeight;
      } finally { URL.revokeObjectURL(u); }
    }
    const c = document.createElement("canvas"); c.width = w; c.height = h;
    const cx = c.getContext("2d", { willReadFrequently: true });
    cx.drawImage(source, 0, 0);
    return { w, h, data: cx.getImageData(0, 0, w, h).data };
  }

  /* name: the file under assets/map without extension; channels: 1 for an id/percentage plane, 2 for elevation (R and G). Returns {w, h, data (RGBA)}. */
  async function load(name, channels) {
    const [r, want] = await Promise.all([decode(`assets/map/${name}.bin?v=${build()}`), checksums().then((s) => s[name])]);
    if (want == null) throw new Error(`no checksum is recorded for the raster ${name}`);
    const got = adler32(r.data, channels);
    if (got !== want) throw new Error(`the raster ${name} does not match its checksum (read ${got}, expected ${want}): something between the server and the browser altered the image`);
    return r;
  }

  const api = { adler32, load };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Raster = api;
})(typeof window !== "undefined" ? window : globalThis);
