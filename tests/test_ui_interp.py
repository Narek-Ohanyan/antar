"""ui/interp.js (the map's interpolation) against an independent numpy implementation of the same
definition, run through node. If node is missing the test is skipped, not passed."""
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

JS = Path(__file__).resolve().parent.parent / "ui" / "interp.js"
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def run_js(payload):
    code = ("const I=require(%r);const d=JSON.parse(require('fs').readFileSync(0,'utf8'));"
            "const it=I.build(d.nx,d.ny,d.tx,d.ty,d.k,d.p);const out=Array.from(I.apply(it,d.v));"
            "const loo=I.leaveOneOut(d.nx,d.ny,d.v,d.k,d.p);"
            "console.log(JSON.stringify({out,rmse:loo.rmse,r2:loo.r2,pred:Array.from(loo.pred)}));") % str(JS)
    r = subprocess.run(["node", "-e", code], input=json.dumps(payload), capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def idw_ref(nx, ny, v, tx, ty, k, p, exclude_self=False):
    out = np.full(len(tx), np.nan)
    for t in range(len(tx)):
        d = np.hypot(nx - tx[t], ny - ty[t])
        if exclude_self:
            d[t] = np.inf
        order = np.lexsort((np.arange(len(d)), d))[:k]
        dn, vn = d[order], v[order]
        w = np.zeros(len(order))
        if (dn < 1).any():
            w[int(np.argmax(dn < 1))] = 1.0
        else:
            w = 1.0 / dn ** p
        good = np.isfinite(vn)
        if w[good].sum() > 0:
            out[t] = (w[good] * vn[good]).sum() / w[good].sum()
    return out


@pytest.mark.parametrize("n,k", [(25, 8), (60, 8), (200, 6), (9, 8)])
def test_matches_numpy_reference(n, k):
    rng = np.random.default_rng(n)
    nx, ny = rng.uniform(0, 3e5, n), rng.uniform(0, 3e5, n)
    v = np.sin(nx / 5e4) + 0.1 * rng.normal(size=n)
    v[rng.integers(0, n)] = np.nan                       # one missing node value
    tx, ty = rng.uniform(-2e4, 3.2e5, 300), rng.uniform(-2e4, 3.2e5, 300)
    got = run_js({"nx": nx.tolist(), "ny": ny.tolist(), "tx": tx.tolist(), "ty": ty.tolist(),
                  "v": [None if np.isnan(a) else a for a in v], "k": k, "p": 2})
    ref = idw_ref(nx, ny, v, tx, ty, min(k, n), 2)
    np.testing.assert_allclose(np.array(got["out"], dtype=float), ref, rtol=2e-5, atol=2e-6, equal_nan=True)


def test_leave_one_out_matches_reference_and_exact_at_nodes():
    rng = np.random.default_rng(7)
    n = 40
    nx, ny = rng.uniform(0, 1e5, n), rng.uniform(0, 1e5, n)
    v = nx / 1e5 + 0.05 * rng.normal(size=n)
    got = run_js({"nx": nx.tolist(), "ny": ny.tolist(), "tx": nx.tolist(), "ty": ny.tolist(), "v": v.tolist(), "k": 8, "p": 2})
    np.testing.assert_allclose(got["out"], v, rtol=1e-5)          # interpolation reproduces the nodes
    pred = idw_ref(nx, ny, v, nx, ny, 8, 2, exclude_self=True)
    np.testing.assert_allclose(got["pred"], pred, rtol=2e-5)
    r2 = 1 - ((v - pred) ** 2).sum() / ((v - v.mean()) ** 2).sum()
    assert got["r2"] == pytest.approx(r2, abs=1e-5)
    assert got["rmse"] == pytest.approx(np.sqrt(((v - pred) ** 2).mean()), rel=1e-5)


def test_leave_one_out_ignores_missing_nodes():
    """null node values must be skipped, not read as 0 (isFinite(null) is true in JS)."""
    rng = np.random.default_rng(3)
    n = 30
    nx, ny = rng.uniform(0, 1e5, n), rng.uniform(0, 1e5, n)
    v = 5.0 + nx / 1e5
    vm = v.copy()
    vm[[4, 11]] = np.nan
    got = run_js({"nx": nx.tolist(), "ny": ny.tolist(), "tx": nx.tolist(), "ty": ny.tolist(),
                  "v": [None if np.isnan(a) else a for a in vm], "k": 8, "p": 2})
    pred = idw_ref(nx, ny, vm, nx, ny, 8, 2, exclude_self=True)
    good = np.isfinite(vm)
    assert got["rmse"] == pytest.approx(np.sqrt(((vm[good] - pred[good]) ** 2).mean()), rel=1e-5)
    assert got["rmse"] < 0.5          # would be ~5 if the missing nodes were read as 0

def run_js_elev(payload):
    code = ("const I=require(%r);const d=JSON.parse(require('fs').readFileSync(0,'utf8'));"
            "const it=I.build(d.nx,d.ny,d.tx,d.ty,d.k,d.p);const b=I.olsSlope(d.nz,d.v);"
            "const out=Array.from(I.applyElev(it,d.v,d.nz,d.tz,b));"
            "const loo=I.leaveOneOutElev(d.nx,d.ny,d.nz,d.v,d.k,d.p);"
            "console.log(JSON.stringify({b,out,r2:loo.r2,rmse:loo.rmse,pred:Array.from(loo.pred)}));") % str(JS)
    r = subprocess.run(["node", "-e", code], input=json.dumps(payload), capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def ols_ref(z, v, skip=None):
    m = np.isfinite(z) & np.isfinite(v)
    if skip is not None:
        m[skip] = False
    if m.sum() < 3:
        return 0.0
    zc, vc = z[m] - z[m].mean(), v[m] - v[m].mean()
    return float((zc * vc).sum() / (zc ** 2).sum()) if (zc ** 2).sum() > 0 else 0.0


def idw_elev_ref(nx, ny, nz, v, tx, ty, tz, k, b, skip_self=False):
    out = np.full(len(tx), np.nan)
    for t in range(len(tx)):
        d = np.hypot(nx - tx[t], ny - ty[t])
        if skip_self:
            d[t] = np.inf
        order = np.lexsort((np.arange(len(d)), d))[:k]
        dn = d[order]
        w = np.zeros(len(order))
        if (dn < 1).any():
            w[int(np.argmax(dn < 1))] = 1.0
        else:
            w = 1.0 / dn ** 2
        adj = v[order] + b * (tz[t] - nz[order])
        good = np.isfinite(adj)
        out[t] = (w[good] * adj[good]).sum() / w[good].sum() if w[good].sum() > 0 else np.nan
    return out


def test_elevation_adjusted_matches_numpy_reference():
    rng = np.random.default_rng(11)
    n = 50
    nx, ny, nz = rng.uniform(0, 2e5, n), rng.uniform(0, 2e5, n), rng.uniform(500, 3500, n)
    v = 15 - 0.0065 * nz + 0.3 * rng.normal(size=n)
    v[7] = np.nan
    tx, ty, tz = rng.uniform(0, 2e5, 200), rng.uniform(0, 2e5, 200), rng.uniform(500, 4000, 200)
    got = run_js_elev({"nx": nx.tolist(), "ny": ny.tolist(), "nz": nz.tolist(), "v": [None if np.isnan(a) else a for a in v],
                       "tx": tx.tolist(), "ty": ty.tolist(), "tz": tz.tolist(), "k": 8, "p": 2})
    b = ols_ref(nz, v)
    assert got["b"] == pytest.approx(b, rel=1e-9)
    np.testing.assert_allclose(got["out"], idw_elev_ref(nx, ny, nz, v, tx, ty, tz, 8, b), rtol=3e-5, atol=3e-5)
    pred = np.array([idw_elev_ref(nx, ny, nz, v, nx, ny, nz, 8, ols_ref(nz, v, i), skip_self=True)[i] for i in range(n)])
    np.testing.assert_allclose(got["pred"], pred, rtol=3e-5, atol=3e-5, equal_nan=True)


def test_elevation_adjustment_recovers_a_lapse_rate_that_plain_idw_cannot():
    """Temperature falling 6.5 K per km, nodes scattered over terrain of very different heights."""
    rng = np.random.default_rng(5)
    n = 40
    nx, ny, nz = rng.uniform(0, 2e5, n), rng.uniform(0, 2e5, n), rng.uniform(500, 3500, n)
    v = 20 - 0.0065 * nz
    got = run_js_elev({"nx": nx.tolist(), "ny": ny.tolist(), "nz": nz.tolist(), "v": v.tolist(),
                       "tx": nx.tolist(), "ty": ny.tolist(), "tz": nz.tolist(), "k": 8, "p": 2})
    plain = run_js({"nx": nx.tolist(), "ny": ny.tolist(), "tx": nx.tolist(), "ty": ny.tolist(), "v": v.tolist(), "k": 8, "p": 2})
    assert got["b"] == pytest.approx(-0.0065, rel=1e-6)
    assert got["r2"] > 0.999 and plain["r2"] < 0.5
