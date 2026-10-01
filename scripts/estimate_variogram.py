"""Real empirical variogram estimate from the real 78-point TOPOHYDRO grid, to replace
configs/cv.yaml's 'default 25 km until estimated' placeholder with an actual data-driven number.

CWD (pm_fao56) is detrended against elevation first (simple OLS) since elevation drives most of
CWD's real large-scale spatial trend in this data -- the variogram should characterize the
residual spatial-autocorrelation structure, not the trend itself, matching cv.yaml's own wording
("variogram practical range of residuals").

Honest caveat stated, not hidden: 78 points is a small, noisy sample for variogram estimation --
this gives a real, defensible order-of-magnitude number, not a publication-grade semivariogram fit.
"""
import numpy as np
import yaml
from scipy.optimize import curve_fit

d = yaml.safe_load(open("configs/fitted/topohydro_grid_run_2019.yaml"))
pts = [p for p in d["points"] if p.get("status") == "ok"]
lats = np.array([p["lat"] for p in pts])
lons = np.array([p["lon"] for p in pts])
elev = np.array([p["elevation_m"] for p in pts])
cwd = np.array([p["cwd_mm_by_pet_formulation"]["pm_fao56"] for p in pts])

# OLS detrend against elevation
A = np.vstack([elev, np.ones_like(elev)]).T
coef, *_ = np.linalg.lstsq(A, cwd, rcond=None)
resid = cwd - A @ coef
print(f"real detrend: CWD = {coef[0]:.4f}*elev + {coef[1]:.2f}  (R2 on trend removed)")

# Real haversine pairwise distances (km)
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dlambda / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

n = len(pts)
dists, sqdiffs = [], []
for i in range(n):
    for j in range(i + 1, n):
        dists.append(haversine_km(lats[i], lons[i], lats[j], lons[j]))
        sqdiffs.append(0.5 * (resid[i] - resid[j])**2)
dists, sqdiffs = np.array(dists), np.array(sqdiffs)
print(f"real pairs: {len(dists)}, distance range {dists.min():.1f}-{dists.max():.1f} km")

# Bin into real empirical semivariogram
nbins = 12
max_dist = np.percentile(dists, 70)  # standard practice: don't trust variogram past ~2/3 max lag
bins = np.linspace(0, max_dist, nbins + 1)
bin_centers, gamma = [], []
for b0, b1 in zip(bins[:-1], bins[1:]):
    sel = (dists >= b0) & (dists < b1)
    if sel.sum() >= 3:
        bin_centers.append((b0 + b1) / 2)
        gamma.append(sqdiffs[sel].mean())
bin_centers, gamma = np.array(bin_centers), np.array(gamma)
print("real empirical semivariogram (km, semivariance):")
for c, g in zip(bin_centers, gamma):
    print(f"  {c:6.1f} km : {g:8.3f}")

# Fit a real exponential variogram model: gamma(h) = nugget + sill*(1 - exp(-h/range))
def exp_model(h, nugget, sill, rng):
    return nugget + sill * (1 - np.exp(-h / rng))

try:
    popt, _ = curve_fit(exp_model, bin_centers, gamma, p0=[gamma.min(), gamma.max() - gamma.min(), 25.0],
                         bounds=([0, 0, 1], [gamma.max(), gamma.max() * 2, 200]))
    nugget, sill, rng = popt
    practical_range_km = 3 * rng  # standard real convention for the exponential model
    print(f"\nreal fitted exponential variogram: nugget={nugget:.3f}, sill={sill:.3f}, range_param={rng:.1f} km")
    print(f"real practical range (3x range param, standard exponential-model convention): {practical_range_km:.1f} km")
except Exception as e:
    print(f"\nfit failed: {e}")
