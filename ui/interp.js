/* Inverse-distance-weighted (Shepard) interpolation of node values onto map pixels.
   Pure functions, no DOM: loaded by the browser as window.Interp and by node (tests/test_ui_interp.py,
   which checks this file against an independent numpy implementation).

   Definition (the reference the test implements):
     for each target t: take the K nearest nodes by Euclidean distance in projected metres
     (ties broken by lower node index); w_j = 1 / d_j^p, and the single node if d_j < 1 m;
     value = sum(w_j v_j) / sum(w_j) over the neighbours whose value is finite
     (NaN if none are). */
(function (root) {
  const R = 6378137;
  const fin = (v) => v != null && isFinite(v);      // isFinite(null) is true in JS, so null must be excluded explicitly
  const mercX = (lon) => (R * lon * Math.PI) / 180;
  const mercY = (lat) => R * Math.log(Math.tan(Math.PI / 4 + (lat * Math.PI) / 360));

  /* K nearest nodes per target via a uniform bucket grid. exclude[t] (optional) is a node index to skip. */
  function knn(nx, ny, tx, ty, k, exclude) {
    const n = nx.length, nt = tx.length;
    k = Math.min(k, exclude ? n - 1 : n);
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (let i = 0; i < n; i++) { minX = Math.min(minX, nx[i]); maxX = Math.max(maxX, nx[i]); minY = Math.min(minY, ny[i]); maxY = Math.max(maxY, ny[i]); }
    const span = Math.max(maxX - minX, maxY - minY, 1);
    const cell = Math.max(span / Math.max(1, Math.sqrt(n / 2)), 1);
    const gw = Math.floor((maxX - minX) / cell) + 1, gh = Math.floor((maxY - minY) / cell) + 1;
    const buckets = Array.from({ length: gw * gh }, () => []);
    for (let i = 0; i < n; i++) buckets[Math.floor((ny[i] - minY) / cell) * gw + Math.floor((nx[i] - minX) / cell)].push(i);
    const idx = new Int32Array(nt * k), dist = new Float64Array(nt * k);
    const cand = [];
    for (let t = 0; t < nt; t++) {
      const cx = Math.floor((tx[t] - minX) / cell), cy = Math.floor((ty[t] - minY) / cell);
      const maxRing = Math.max(gw, gh, Math.abs(cx), Math.abs(cy), Math.abs(cx - gw), Math.abs(cy - gh)) + 2;
      let order = null;
      cand.length = 0;
      for (let ring = 0; ring <= maxRing; ring++) {
        for (let gy = cy - ring; gy <= cy + ring; gy++) {
          for (let gx = cx - ring; gx <= cx + ring; gx++) {
            if (Math.max(Math.abs(gx - cx), Math.abs(gy - cy)) !== ring) continue;
            if (gx < 0 || gy < 0 || gx >= gw || gy >= gh) continue;
            for (const i of buckets[gy * gw + gx]) if (!exclude || exclude[t] !== i) cand.push(i);
          }
        }
        if (cand.length >= k) {
          /* every node not yet visited is at least `ring * cell` from the target, so once the k-th nearest
             candidate is within that, the k nearest are final */
          order = cand.map((i) => [Math.hypot(nx[i] - tx[t], ny[i] - ty[t]), i]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
          if (order[k - 1][0] <= ring * cell) break;
        }
      }
      for (let j = 0; j < k; j++) { idx[t * k + j] = order[j][1]; dist[t * k + j] = order[j][0]; }
    }
    return { idx, dist, k };
  }

  function weights(nb, power) {
    const { dist, k } = nb, w = new Float64Array(dist.length);
    for (let t = 0; t < dist.length / k; t++) {
      let exact = -1;
      for (let j = 0; j < k; j++) if (dist[t * k + j] < 1) { exact = j; break; }
      for (let j = 0; j < k; j++) w[t * k + j] = exact >= 0 ? (j === exact ? 1 : 0) : 1 / Math.pow(dist[t * k + j], power);
    }
    return w;
  }

  /* Build once per (node set, target set); reuse for every layer. */
  function build(nodeX, nodeY, targetX, targetY, k = 8, power = 2) {
    const nb = knn(nodeX, nodeY, targetX, targetY, k);
    return { idx: nb.idx, dist: nb.dist, w: weights(nb, power), k: nb.k, nt: targetX.length };
  }

  function apply(interp, values, out) {
    const { idx, w, k, nt } = interp;
    out = out || new Float32Array(nt);
    for (let t = 0; t < nt; t++) {
      let s = 0, sw = 0;
      for (let j = 0; j < k; j++) {
        const v = values[idx[t * k + j]];
        if (fin(v)) { s += w[t * k + j] * v; sw += w[t * k + j]; }
      }
      out[t] = sw > 0 ? s / sw : NaN;
    }
    return out;
  }

  /* Leave-one-out check at the nodes themselves: predict each node from the others. */
  function leaveOneOut(nodeX, nodeY, values, k = 8, power = 2) {
    const n = nodeX.length, ex = new Int32Array(n);
    for (let i = 0; i < n; i++) ex[i] = i;
    const nb = knn(nodeX, nodeY, nodeX, nodeY, k, ex);
    const it = { idx: nb.idx, w: weights(nb, power), k: nb.k, nt: n };
    const pred = apply(it, values);
    let sse = 0, sst = 0, m = 0, mean = 0;
    for (let i = 0; i < n; i++) if (fin(values[i]) && fin(pred[i])) { mean += values[i]; m++; }
    mean /= m || 1;
    for (let i = 0; i < n; i++) if (fin(values[i]) && fin(pred[i])) { sse += (values[i] - pred[i]) ** 2; sst += (values[i] - mean) ** 2; }
    return { n: m, rmse: m ? Math.sqrt(sse / m) : NaN, r2: sst > 0 ? 1 - sse / sst : NaN, pred };
  }

  const api = { R, mercX, mercY, knn, weights, build, apply, leaveOneOut };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Interp = api;
})(typeof window !== "undefined" ? window : globalThis);
