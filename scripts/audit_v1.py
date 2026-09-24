"""Reproducible audit of the prior-project v1 training table and Random-Forest model.

Everything the v2 concept note says about v1 is computed here from
``Armenia_ML_Training_Data.parquet``; nothing is quoted from memory.

    python scripts/audit_v1.py --parquet path/to/Armenia_ML_Training_Data.parquet --out docs/audit_v1.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold, train_test_split

from antar.io.grid import audit_dvpd_dtmax_ratio
from antar.validation.splits import block_kfold

FEATS = ["Delta_VPD_GS", "Delta_Tmax_GS", "Delta_P_GS"]
TARGET = "vhm_std"


def rf(n=500, depth=15, seed=42):
    return RandomForestRegressor(n_estimators=n, max_depth=depth, random_state=seed, n_jobs=-1)


def to_km(df):
    lat0 = df["y"].mean()
    kx = 111.32 * math.cos(math.radians(lat0))
    ky = 110.57
    return (df["x"] - df["x"].min()).values * kx, (df["y"] - df["y"].min()).values * ky


def cv_scores(df, feats, gx, gy, block_km, k=5, seed=0, n=200, buffer_km=0.0):
    y = df[TARGET].values
    pred = np.full(len(df), np.nan)
    for tr, te in block_kfold(gx, gy, block_km, k=k, buffer=buffer_km, seed=seed):
        m = rf(n).fit(df.iloc[tr][feats], y[tr])
        pred[te] = m.predict(df.iloc[te][feats])
    ok = ~np.isnan(pred)
    return r2_score(y[ok], pred[ok]), mean_absolute_error(y[ok], pred[ok])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--out", default="docs/audit_v1.json")
    a = ap.parse_args()

    df = pd.read_parquet(a.parquet)
    y = df[TARGET]
    X = df[FEATS]
    gx, gy = to_km(df)
    R: dict = {"n_pixels": int(len(df))}

    # 1. replication of the v1 evaluation (random 80/20, RF 500 trees, depth 15, seed 42)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42)
    m = rf().fit(Xtr, ytr)
    p = m.predict(Xte)
    R["v1_random_split"] = {"r2": r2_score(yte, p), "mae": mean_absolute_error(yte, p)}
    R["v1_importances_pct"] = dict(zip(FEATS, (m.feature_importances_ * 100).round(1).tolist()))
    R["null_train_mean"] = {"mae": mean_absolute_error(yte, np.full(len(yte), ytr.mean())),
                             "r2": r2_score(yte, np.full(len(yte), ytr.mean()))}

    # 2. coordinates-only baseline, random split
    XYtr, XYte, ytr2, yte2 = train_test_split(df[["x", "y"]], y, test_size=0.2, random_state=42)
    mxy = rf().fit(XYtr, ytr2)
    R["coords_only_random_split_r2"] = r2_score(yte2, mxy.predict(XYte))

    # 3. blocked CV versus block size, for v1 features and for coordinates only
    R["blocked_cv"] = []
    for bs in (5, 10, 25, 50, 75):
        rows = {"block_km": bs}
        for label, feats in (("v1_features", FEATS), ("coords_only", ["x", "y"])):
            rs = [cv_scores(df, feats, gx, gy, bs, seed=s) for s in range(3)]
            rows[f"{label}_r2_mean"] = float(np.mean([r[0] for r in rs]))
            rows[f"{label}_r2_sd"] = float(np.std([r[0] for r in rs]))
            rows[f"{label}_mae_mean"] = float(np.mean([r[1] for r in rs]))
        R["blocked_cv"].append(rows)
    pred = np.zeros(len(df))
    for tr, te in KFold(5, shuffle=True, random_state=1).split(df):
        pred[te] = rf(200).fit(X.iloc[tr], y.iloc[tr]).predict(X.iloc[te])
    R["random_5fold"] = {"r2": r2_score(y, pred), "mae": mean_absolute_error(y, pred)}

    # 4. are the delta features location proxies?
    R["features_predict_location_r2"] = {}
    for tgt in ("x", "y"):
        a_, b_, c_, d_ = train_test_split(X, df[tgt], test_size=0.2, random_state=42)
        R["features_predict_location_r2"][tgt] = r2_score(d_, rf(200).fit(a_, c_).predict(b_))

    # 5. sigma(H) versus mean height
    s, mu = df["vhm_std"].values, df["vhm_mean"].values
    Xm = df[["vhm_mean"]]
    a_, b_, c_, d_ = train_test_split(Xm, y, test_size=0.2, random_state=42)
    R["sigma_vs_meanH"] = {"corr": float(np.corrcoef(mu, s)[0, 1]),
                            "rf_r2_from_meanH_only": r2_score(d_, rf(200, 8).fit(a_, c_).predict(b_)),
                            "median_cv": float(np.median(s / mu))}

    # 6. extrapolation behaviour of the deployed model (clipping rules copied from app.py)
    def app_predict(dt, dp, dv):
        Xin = pd.DataFrame({"Delta_VPD_GS": (df["Delta_VPD_GS"] + dv).clip(upper=1.5),
                             "Delta_Tmax_GS": (df["Delta_Tmax_GS"] + dt).clip(upper=10.0),
                             "Delta_P_GS": (df["Delta_P_GS"] + dp).clip(lower=-150.0)})
        return m.predict(Xin)

    scen = {"SSP1-2.6 2041-60": (1.02, -1.94, 0.168), "SSP2-4.5 2081-2100": (2.41, -4.44, 0.417),
            "SSP3-7.0 2081-2100": (5.44, -13.50, 1.073), "SSP5-8.5 2081-2100": (5.72, -9.83, 1.045)}
    R["scenarios"] = {}
    for name, (dt, dp, dv) in scen.items():
        pp = app_predict(dt, dp, dv)
        fvs = ((y - pp) / y * 100).clip(0, 100)
        R["scenarios"][name] = {"mean_fvs": float(fvs.mean()), "mean_pred_m": float(pp.mean()),
                                 "frac_fvs_gt50": float((fvs > 50).mean())}
    p0 = m.predict(X)
    fvs0 = ((y - p0) / y * 100).clip(0, 100)
    R["fvs_at_zero_forcing"] = {"mean": float(fvs0.mean()), "frac_gt50": float((fvs0 > 50).mean())}
    R["flat_response"] = {
        "max_abs_diff_dT_1.0_vs_5.72": float(np.abs(app_predict(1.0, -5, 0.3) - app_predict(5.72, -5, 0.3)).max()),
        "max_abs_diff_dT_5.72_vs_9.0": float(np.abs(app_predict(5.72, -5, 0.3) - app_predict(9.0, -5, 0.3)).max()),
        "train_range_dTmax": [float(df.Delta_Tmax_GS.min()), float(df.Delta_Tmax_GS.max())],
        "train_range_dVPD": [float(df.Delta_VPD_GS.min()), float(df.Delta_VPD_GS.max())],
        "train_range_dP": [float(df.Delta_P_GS.min()), float(df.Delta_P_GS.max())],
    }
    # sweep for the figure: prediction versus additional warming, the other two deltas held at zero
    sweep = np.linspace(0, 8, 33)
    R["sweep_dTmax"] = {"dT": sweep.tolist(), "mean_pred": [float(app_predict(t, 0, 0).mean()) for t in sweep]}

    # 7. units audit and spatial dependence
    R["dvpd_dtmax_ratio"] = float(df.Delta_VPD_GS.mean() / df.Delta_Tmax_GS.mean())
    R["dvpd_dtmax_audit"] = audit_dvpd_dtmax_ratio(df.Delta_VPD_GS.mean(), df.Delta_Tmax_GS.mean()).message
    xy = np.c_[gx, gy]
    d, i = cKDTree(xy).query(xy, k=2)
    R["median_nn_distance_km"] = float(np.median(d[:, 1]))
    R["nn_corr"] = {c: float(np.corrcoef(df[c].values, df[c].values[i[:, 1]])[0, 1]) for c in FEATS + [TARGET]}
    R["cell_area_km2"] = float(111.32 * math.cos(math.radians(df.y.mean())) / 120 * 110.57 / 120)

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(R, indent=2))
    print(json.dumps(R, indent=2))


if __name__ == "__main__":
    main()
