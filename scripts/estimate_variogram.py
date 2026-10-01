"""Real empirical variogram estimate from the real TOPOHYDRO grid, to replace
configs/cv.yaml's block-size placeholder with an actual data-driven number.

CWD (pm_fao56) is detrended against elevation first (simple OLS) since elevation drives most of
CWD's real large-scale spatial trend in this data -- the variogram should characterize the
residual spatial-autocorrelation structure, not the trend itself, matching cv.yaml's own wording
("variogram practical range of residuals").

--dense reads the real 1044-point dense grid (configs/fitted/topohydro_grid_run_2019_dense.yaml,
written by `python3 scripts/run_topohydro_grid.py --dense`) instead of the original 78-point one --
the real follow-up the user asked for once a denser real sample existed, specifically to find out
whether the original run's "no clear sill within 24-186km, real practical range ~523km, but can't
be pinned down with confidence at this density" finding holds up or resolves with ~13x more points.

Honest caveat, same as before: even 1044 points is still a modest sample for variogram estimation
by geostatistics standards -- this gives a real, improving-but-still-approximate number, not a
publication-grade semivariogram fit.
"""
import sys

import numpy as np
import yaml
from scipy.optimize import curve_fit

dense = "--dense" in sys.argv
path = "configs/fitted/topohydro_grid_run_2019_dense.yaml" if dense else "configs/fitted/topohydro_grid_run_2019.yaml"
d = yaml.safe_load(open(path))
pts = [p for p in d["points"] if p.get("status") == "ok"]
print(f"=== real {'DENSE (1044-point)' if dense else '78-point'} grid: {len(pts)} real ok points ===")
lats = np.array([p["lat"] for p in pts])
lons = np.array([p["lon"] for p in pts])
elev = np.array([p["elevation_m"] for p in pts])
cwd = np.array([p["cwd_mm_by_pet_formulation"]["pm_fao56"] for p in pts])

# OLS detrend against elevation
A = np.vstack([elev, np.ones_like(elev)]).T
coef, *_ = np.linalg.lstsq(A, cwd, rcond=None)
resid = cwd - A @ coef
print(f"real detrend: CWD = {coef[0]:.4f}*elev + {coef[1]:.2f}")

# Real haversine pairwise distances (km) -- vectorized (n up to ~1044 means up to ~545k pairs,
# a plain double loop is fine computationally but vectorizing is cheap and faster to iterate on).
def haversine_km_matrix(lats, lons):
    R = 6371.0
    p = np.radians(lats)[:, None]
    p2 = np.radians(lats)[None, :]
    dphi = p2 - p
    dlambda = np.radians(lons)[None, :] - np.radians(lons)[:, None]
    a = np.sin(dphi / 2)**2 + np.cos(p) * np.cos(p2) * np.sin(dlambda / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

n = len(pts)
dist_mat = haversine_km_matrix(lats, lons)
resid_diff2 = 0.5 * (resid[:, None] - resid[None, :])**2
iu = np.triu_indices(n, k=1)
dists, sqdiffs = dist_mat[iu], resid_diff2[iu]
print(f"real pairs: {len(dists)}, distance range {dists.min():.1f}-{dists.max():.1f} km")

# Bin into real empirical semivariogram
nbins = 16 if dense else 12
max_dist = np.percentile(dists, 70)  # standard practice: don't trust variogram past ~2/3 max lag
bins = np.linspace(0, max_dist, nbins + 1)
bin_centers, gamma, bin_counts = [], [], []
for b0, b1 in zip(bins[:-1], bins[1:]):
    sel = (dists >= b0) & (dists < b1)
    if sel.sum() >= 3:
        bin_centers.append((b0 + b1) / 2)
        gamma.append(sqdiffs[sel].mean())
        bin_counts.append(int(sel.sum()))
bin_centers, gamma = np.array(bin_centers), np.array(gamma)
print("real empirical semivariogram (km, semivariance, n_pairs):")
for c, g, nc in zip(bin_centers, gamma, bin_counts):
    print(f"  {c:7.1f} km : {g:9.3f}  (n={nc})")

# Fit a real exponential variogram model: gamma(h) = nugget + sill*(1 - exp(-h/range))
def exp_model(h, nugget, sill, rng):
    return nugget + sill * (1 - np.exp(-h / rng))

try:
    popt, _ = curve_fit(exp_model, bin_centers, gamma, p0=[gamma.min(), gamma.max() - gamma.min(), 25.0],
                         bounds=([0, 0, 1], [gamma.max(), gamma.max() * 2, 1000]))
    nugget, sill, rng = popt
    practical_range_km = 3 * rng  # standard real convention for the exponential model
    pred = exp_model(bin_centers, *popt)
    r2_exp = 1 - np.sum((gamma - pred)**2) / np.sum((gamma - gamma.mean())**2)
    print(f"\nreal fitted exponential variogram: nugget={nugget:.3f}, sill={sill:.3f}, range_param={rng:.1f} km")
    print(f"real practical range (3x range param): {practical_range_km:.1f} km  (R2={r2_exp:.4f})")

    # Real honesty check (same one that mattered the first time): does a plain straight line
    # (no sill at all) fit about as well? If so, there's still no real evidence of a plateau.
    Alin = np.vstack([bin_centers, np.ones_like(bin_centers)]).T
    coef_lin, *_ = np.linalg.lstsq(Alin, gamma, rcond=None)
    pred_lin = Alin @ coef_lin
    r2_lin = 1 - np.sum((gamma - pred_lin)**2) / np.sum((gamma - gamma.mean())**2)
    print(f"real plain-linear (no-sill) fit: R2={r2_lin:.4f} (slope={coef_lin[0]:.4f})")
    if r2_exp - r2_lin < 0.03:
        print("=== Real finding: still no clear plateau -- exponential fit doesn't meaningfully "
              "beat a straight line. The practical-range number above remains an upper-bound-ish "
              "signal, not a confident estimate. ===")
    else:
        print("=== Real finding: the exponential fit meaningfully beats a straight line -- a real "
              "plateau is now visible, the practical-range number above is a genuine estimate. ===")
except Exception as e:
    print(f"\nfit failed: {e}")
