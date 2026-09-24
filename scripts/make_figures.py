"""Figures for the v2 concept note.

Fig. 1  (real data)     v1 audit: blocked-CV skill vs. block size; response of the v1 model beyond its training support.
Fig. 3  (illustrative)  plant-hydraulic engine: vulnerability curves, g_min phase transition, two-phase failure index.
Fig. 4  (illustrative)  controlled experiment with known truth: why a hybrid, AOA-gated model beats a pure learner under extrapolation.

Figures 3 and 4 use placeholder trait values and synthetic data *by design*; they demonstrate mechanisms and failure modes,
not skill.  Fig. 1 is computed from the v1 training table (run scripts/audit_v1.py first).

    python scripts/make_figures.py --audit docs/audit_v1.json --out docs/figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from antar.climate.vapour import saturation_vapour_pressure as es
from antar.hazard.gating import aoa_weight, blend_hazards
from antar.hazard.models import CloglogGLM
from antar.hydraulics import twophase
from antar.hydraulics.gmin import gmin_temperature
from antar.hydraulics.vulnerability import plc, psi_at_plc
from antar.validation.aoa import AreaOfApplicability

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "verm": "#D55E00", "sky": "#56B4E9",
      "pink": "#CC79A7", "yellow": "#F0E442", "grey": "#666666", "black": "#111111"}

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7, "legend.frameon": False, "legend.fontsize": 7.2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.bbox": "tight", "mathtext.fontset": "dejavusans",
})


def panel_label(ax, s, dx=-0.13, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.pdf")
    fig.savefig(out / f"{name}.png", dpi=220)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ Fig. 1
def fig_v1_audit(audit: dict, out: Path):
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.75), gridspec_kw={"wspace": 0.38})
    bc = audit["blocked_cv"]
    bs = [r["block_km"] for r in bc]
    a = ax[0]
    a.axhline(0, color=OI["black"], lw=0.6)
    a.axhline(audit["v1_random_split"]["r2"], color=OI["grey"], ls="--", lw=0.9)
    a.text(76, audit["v1_random_split"]["r2"] + 0.02, "v1 random 80/20 split", ha="right", va="bottom", fontsize=7, color=OI["grey"])
    a.axhline(audit["coords_only_random_split_r2"], color=OI["orange"], ls="--", lw=0.9)
    a.text(76, audit["coords_only_random_split_r2"] + 0.02, "coordinates only, random split", ha="right", va="bottom", fontsize=7, color=OI["orange"])
    for key, col, lab in (("v1_features", OI["blue"], "v1 features, blocked CV"), ("coords_only", OI["verm"], "coordinates only, blocked CV")):
        m = np.array([r[f"{key}_r2_mean"] for r in bc])
        s = np.array([r[f"{key}_r2_sd"] for r in bc])
        a.errorbar(bs, m, yerr=s, color=col, marker="o", ms=3.5, lw=1.4, capsize=2, label=lab)
    a.set_xscale("log")
    a.set_xticks(bs)
    a.set_xticklabels([str(b) for b in bs])
    a.set_xlabel("spatial block size (km)")
    a.set_ylabel(r"out-of-block $R^2$ for $\sigma(H)$")
    a.set_ylim(-0.78, 0.62)
    a.legend(loc="lower left", handlelength=1.6)
    panel_label(a, "a")

    b = ax[1]
    sw = audit["sweep_dTmax"]
    lo, hi = audit["flat_response"]["train_range_dTmax"]
    b.axvspan(lo, hi, color=OI["blue"], alpha=0.55, lw=0)
    b.annotate("training support of $\\Delta T_{max}$\n(0.05-0.10 K)", xy=(0.075, 5.9), xytext=(0.75, 6.05), fontsize=7, color=OI["blue"], va="center",
               arrowprops=dict(arrowstyle="-", lw=0.6, color=OI["blue"]))
    b.plot(sw["dT"], sw["mean_pred"], color=OI["verm"], lw=1.6)
    b.plot(sw["dT"][:2], sw["mean_pred"][:2], color=OI["verm"], lw=1.6, marker="o", ms=3.2)
    for x0, lab in ((1.02, "SSP1-2.6\nmid-century"), (5.72, "SSP5-8.5\nend-century")):
        b.axvline(x0, color=OI["grey"], lw=0.7, ls=":")
        b.text(x0 + 0.1, 3.72, lab, fontsize=6.6, color=OI["grey"], va="bottom")
    b.set_xlabel(r"warming added to $\Delta T_{max}$ (K)")
    b.set_ylabel(r"mean predicted $\sigma(H)$ (m)")
    b.set_ylim(3.6, 6.4)
    b.annotate("one step outside the range moves the\nmean by 1.5 m; nothing changes after that", xy=(3.4, 4.47), xytext=(2.0, 5.15),
               fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6, color=OI["grey"]), color=OI["black"])
    panel_label(b, "b")
    save(fig, out, "fig1_v1_audit")


# ------------------------------------------------------------------------------------------------ Fig. 3
def fig_hydraulics(out: Path):
    groups = {
        "mesic broadleaf": dict(p50=-3.0, slope=40, psi_close=-2.2, g25=3.0, tp=38, lethal=88, cap=30000, col=OI["blue"]),
        "ring-porous oak": dict(p50=-2.8, slope=35, psi_close=-2.4, g25=3.5, tp=38, lethal=88, cap=30000, col=OI["green"]),
        "pine": dict(p50=-3.6, slope=30, psi_close=-2.6, g25=1.5, tp=40, lethal=50, cap=20000, col=OI["orange"]),
        "juniper": dict(p50=-9.0, slope=25, psi_close=-5.0, g25=1.0, tp=42, lethal=50, cap=15000, col=OI["verm"]),
    }
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"wspace": 0.42})
    psi = np.linspace(0, -12, 400)
    a = ax[0]
    for name, g in groups.items():
        a.plot(-psi, plc(psi, g["p50"], g["slope"]), color=g["col"], lw=1.4, label=name)
        a.plot(-g["psi_close"], plc(g["psi_close"], g["p50"], g["slope"]), marker="v", ms=4.5, color=g["col"], mec="white", mew=0.5)
    a.set_xlabel(r"water potential $-\psi$ (MPa)")
    a.set_ylabel("percent loss of conductivity (%)")
    a.set_xlim(0, 12)
    a.text(11.8, 14, r"$\blacktriangledown$ = $\psi_{close}$", fontsize=7, color=OI["grey"], ha="right")
    handles, labels = a.get_legend_handles_labels()
    panel_label(a, "a")

    b = ax[1]
    T = np.linspace(20, 48, 300)
    vpd = 0.35 * es(T)                                            # 24-h VPD at fixed RH = 30 % (illustrative)
    for name, g in groups.items():
        t = twophase.Traits(p50=g["p50"], slope=g["slope"], psi_close=g["psi_close"], capacitance=g["cap"], g25=g["g25"], tp=g["tp"], lethal_plc=g["lethal"])
        e_min = gmin_temperature(T, t.g25, t.tp) * vpd / 82.0
        b.semilogy(T, t.buffer / (e_min * 86400.0), color=g["col"], lw=1.4)
        b.plot(g["tp"], t.buffer / (gmin_temperature(g["tp"], t.g25, t.tp) * 0.35 * es(g["tp"]) / 82.0 * 86400.0), marker="|", ms=8, color=g["col"], mew=1.4)
    b.set_xlabel(r"hot-spell $T_{max}$ ($^\circ$C)")
    b.set_ylabel(r"phase-2 time to failure $t_{crit}$ (days)")
    b.text(0.97, 0.97, r"$|$ marks $T_p$" + "\n" + r"$g_{min}$ rises $\times$4.8" + "\n" + r"per 10 K above $T_p$", transform=b.transAxes, fontsize=6.6, color=OI["grey"], ha="right", va="top", linespacing=1.25)
    panel_label(b, "b")

    c = ax[2]
    n = 45
    d = np.arange(n)
    psi_soil = np.maximum(-0.4 - 0.11 * d, -6.0)
    t = twophase.Traits(**{k: groups["mesic broadleaf"][k2] for k, k2 in (("p50", "p50"), ("slope", "slope"), ("psi_close", "psi_close"), ("g25", "g25"), ("tp", "tp"), ("capacitance", "cap"))})
    warm_T = np.where((d >= 24) & (d <= 36), 41.0, 34.0)
    warm_v = 0.35 * es(warm_T)
    cool_T = np.full(n, 24.0)
    cool_v = 0.35 * es(cool_T)
    for TT, VV, col, lab in ((warm_T, warm_v, OI["pink"], "warm, low elevation\n(41 $^\\circ$C spell, d24-36)"), (cool_T, cool_v, OI["black"], "cool, high elevation\n(24 $^\\circ$C)")):
        sim = twophase.simulate_two_phase(psi_soil, TT, VV, 82.0, t)
        h = sim["hfi"].copy()
        over = np.where(h > 1.0)[0]
        if over.size:
            h[over[0] + 2:] = np.nan                      # the index has no meaning after failure
        c.plot(d, h, color=col, lw=1.5, label=lab)
    c.axhline(1.0, color=OI["grey"], lw=0.8, ls="--")
    c.text(0.0, 1.03, "failure at HFI = 1", fontsize=6.2, va="bottom", ha="left", color=OI["grey"])
    c.set_xlim(-1, 45)
    c.set_ylim(0, 1.85)
    c.set_xlabel("day of identical soil dry-down")
    c.set_ylabel("hydraulic-failure index, HFI")
    c.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.4, fontsize=6.4, labelspacing=0.5, frameon=False)
    panel_label(c, "c")
    fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.36, 1.07), handlelength=1.6, columnspacing=1.4)
    fig.text(0.995, -0.11, "illustrative: placeholder trait values, fixed RH = 30 %", ha="right", va="bottom", fontsize=6.2, color=OI["grey"], style="italic")
    save(fig, out, "fig3_hydraulics")


# ------------------------------------------------------------------------------------------------ Fig. 4
def _sample_traits(rng, n, tp_mean, bias=False):
    """Tree-level trait draws.  ``bias=True`` gives the (deliberately imperfect) prior used by the mechanistic sub-model:
    a +2 C systematic error in the lab-measured T_p, wider g_min spread, +15 % capacitance."""
    p50 = rng.normal(-3.0, 0.45 if not bias else 0.6, n)
    slope = np.full(n, 40.0)
    psi_close = rng.normal(-2.2, 0.25, n)
    g25 = 3.0 * np.exp(rng.normal(0.0, 0.30 if not bias else 0.4, n))
    cap = 30000.0 * np.exp(rng.normal(0.0, 0.25, n)) * (1.0 if not bias else 1.15)
    tp = rng.normal(tp_mean + (2.0 if bias else 0.0), 2.0, n)
    return p50, slope, psi_close, g25, cap, tp


def hazard_grid(T, W, tp_mean=38.0, n_draw=400, seed=0, bias=False):
    """Annual hazard h(T, W) = p_spell * E_traits[ exp(-t_crit / mu_D(W)) ]  (closed-stomata spell length ~ Exponential)."""
    rng = np.random.default_rng(seed)
    p50, slope, psi_close, g25, cap, tp = _sample_traits(rng, n_draw, tp_mean, bias)
    psi_crit = psi_at_plc(88.0, p50, slope)
    buf = np.maximum(cap * (psi_close - psi_crit), 1.0)
    T = np.atleast_1d(T).astype(float)
    W = np.atleast_1d(W).astype(float)
    vpd = 0.35 * es(T)
    e_min = gmin_temperature(T[:, None], g25[None, :], tp[None, :]) * vpd[:, None] / 82.0
    tcrit = buf[None, :] / (e_min * 86400.0)                                # (nT, ndraw)
    mu = (1.2 + 1.8 * (1.0 - W)) if not bias else (1.5 + 1.2 * (1.0 - W))    # mean closed-stomata spell, days
    return 0.6 * np.mean(np.exp(-tcrit / mu[:, None]), axis=1)


MODELS = ["random_forest", "cloglog_glm", "mechanistic_only", "hybrid_gated"]


def _world(tp_mean):
    from scipy.interpolate import RegularGridInterpolator
    Tg = np.linspace(23.0, 45.0, 89)
    Wg = np.linspace(0, 1, 21)
    truth = RegularGridInterpolator((Wg, Tg), np.array([hazard_grid(Tg, np.full_like(Tg, w), tp_mean) for w in Wg]))
    mech = RegularGridInterpolator((Wg, Tg), np.array([hazard_grid(Tg, np.full_like(Tg, w), tp_mean, seed=5, bias=True) for w in Wg]))
    return truth, mech


def _fit_and_score(world, seed, T_lo=24.0, T_hi=34.0, T_max=44.0, n=40000, keep_models=False):
    truth, mech = world
    rng = np.random.default_rng(seed)
    T = rng.uniform(T_lo, T_hi, n)
    W = rng.uniform(0, 1, n)
    y = rng.random(n) < truth(np.c_[W, T])
    X = np.c_[T, W]
    rf = RandomForestClassifier(200, min_samples_leaf=400, random_state=seed, n_jobs=-1).fit(X, y)
    glm = CloglogGLM(l2=1e-6).fit(X, y)
    folds = np.digitize(T, np.linspace(T_lo, T_hi, 6)[1:-1])           # leave-temperature-band-out folds
    aoa = AreaOfApplicability().fit(X, folds)

    def predict_all(Tq, Wq):
        Xq = np.c_[Tq, Wq]
        h_rf = np.clip(rf.predict_proba(Xq)[:, 1], 1e-5, 0.999)
        h_gl = np.clip(glm.predict_hazard(Xq), 1e-5, 0.999)
        h_me = np.clip(mech(np.c_[Wq, Tq]), 1e-5, 0.999)
        di = aoa.dissimilarity_index(Xq)
        w = aoa_weight(di, aoa.threshold_, kappa=2.0)
        h_hy = blend_hazards(h_gl, h_me, w)          # statistical component = cloglog GLM; physics outside the AOA
        return dict(zip(MODELS, (h_rf, h_gl, h_me, h_hy))), di, w

    def mae(lo, hi, m=6000):
        Tt, Wt = rng.uniform(lo, hi, m), rng.uniform(0, 1, m)
        ht = truth(np.c_[Wt, Tt])
        preds, _, _ = predict_all(Tt, Wt)
        return {k: float(np.mean(np.abs(v - ht))) for k, v in preds.items()}, float(ht.mean())

    e_in, m_in = mae(T_lo, T_hi)
    e_out, m_out = mae(T_hi + 0.5, T_max)
    res = dict(mae_inside=e_in, mae_beyond=e_out, true_inside=m_in, true_beyond=m_out, aoa_threshold=float(aoa.threshold_), n_events=int(y.sum()))
    if keep_models:
        res["predict_all"] = predict_all
    return res


def fig_extrapolation(out: Path, report: dict):
    T_lo, T_hi, T_max = 24.0, 34.0, 44.0
    worlds = {"A": (38.0, "world A: $T_p$ = 38 $^\\circ$C\n(g$_{min}$ phase transition in range)"),
              "B": (47.0, "world B: $T_p$ = 47 $^\\circ$C\n(no phase transition in range)")}
    W_ = {k: _world(v[0]) for k, v in worlds.items()}
    n_rep = 10
    stats = {}
    for k in worlds:
        runs = [_fit_and_score(W_[k], seed=100 * (k == "B") + s) for s in range(n_rep)]
        stats[k] = {
            "mae_beyond_mean": {m: float(np.mean([r["mae_beyond"][m] for r in runs])) for m in MODELS},
            "mae_beyond_sd": {m: float(np.std([r["mae_beyond"][m] for r in runs])) for m in MODELS},
            "mae_inside_mean": {m: float(np.mean([r["mae_inside"][m] for r in runs])) for m in MODELS},
            "mean_true_hazard_inside": float(np.mean([r["true_inside"] for r in runs])),
            "mean_true_hazard_beyond": float(np.mean([r["true_beyond"] for r in runs])),
            "aoa_threshold_mean": float(np.mean([r["aoa_threshold"] for r in runs])),
            "n_events_mean": float(np.mean([r["n_events"] for r in runs])),
        }
    report["extrapolation_experiment"] = {"n_replicates": n_rep, "worlds": stats}

    fig = plt.figure(figsize=(7.1, 5.3))
    gs = fig.add_gridspec(2, 2, wspace=0.62, hspace=0.62)
    axs = {"A": fig.add_subplot(gs[0, 0]), "B": fig.add_subplot(gs[0, 1])}
    Tq = np.linspace(T_lo, T_max, 160)
    Wq = np.full_like(Tq, 0.5)
    di_w = None
    for lab, (k, (tp, title)) in zip("ab", worlds.items()):
        run = _fit_and_score(W_[k], seed=1 + (k == "B"), keep_models=True)
        preds, di, w = run["predict_all"](Tq, Wq)
        if di_w is None:
            di_w = (di, w, run["aoa_threshold"])
        a = axs[k]
        a.axvspan(T_lo, T_hi, color=OI["sky"], alpha=0.32, lw=0)
        a.plot(Tq, W_[k][0](np.c_[Wq, Tq]), color=OI["black"], lw=2.0, label="truth")
        a.plot(Tq, preds["random_forest"], color=OI["verm"], lw=1.3, label="random forest")
        a.plot(Tq, preds["cloglog_glm"], color=OI["orange"], lw=1.3, label="cloglog GLM")
        a.plot(Tq, preds["hybrid_gated"], color=OI["green"], lw=1.6, label="hybrid, AOA-gated")
        a.set_title(title, fontsize=7.2, loc="left", pad=4)
        a.set_xlabel(r"hot-spell $T_{max}$ ($^\circ$C), $W=0.5$")
        a.set_ylabel("annual mortality hazard")
        a.set_ylim(-0.005, 0.27)
        a.text(29, 0.255, "training\nsupport", ha="center", va="top", fontsize=6.8, color=OI["blue"])
        panel_label(a, lab, dx=-0.16, dy=1.1)
        if k == "A":
            a.legend(loc="upper left", bbox_to_anchor=(0.02, 0.70), handlelength=1.4, fontsize=6.6)

    c = fig.add_subplot(gs[1, 0])
    di, w, thr = di_w
    c.plot(Tq, di, color=OI["black"], lw=1.3)
    c.axhline(thr, color=OI["grey"], ls="--", lw=0.8)
    c.text(T_lo + 0.3, thr * 1.08, "AOA threshold", ha="left", va="bottom", fontsize=6.6, color=OI["grey"])
    c.set_xlabel(r"hot-spell $T_{max}$ ($^\circ$C)")
    c.set_ylabel("dissimilarity index DI")
    c2 = c.twinx()
    c2.spines["right"].set_visible(True)
    c2.plot(Tq, w, color=OI["green"], lw=1.4)
    c2.set_ylabel(r"weight on learner $w$", color=OI["green"])
    c2.set_ylim(-0.02, 1.05)
    panel_label(c, "c", dx=-0.16, dy=1.1)

    d = fig.add_subplot(gs[1, 1])
    names = ["RF", "GLM", "mech.", "hybrid"]
    cols = [OI["verm"], OI["orange"], OI["sky"], OI["green"]]
    xpos = np.arange(2)
    wd = 0.2
    for i, (m, nm, col) in enumerate(zip(MODELS, names, cols)):
        vals = [stats[k]["mae_beyond_mean"][m] for k in worlds]
        sds = [stats[k]["mae_beyond_sd"][m] for k in worlds]
        d.bar(xpos + (i - 1.5) * wd, vals, wd * 0.92, yerr=sds, color=col, label=nm, capsize=1.5, error_kw=dict(lw=0.7))
    d.set_xticks(xpos)
    d.set_xticklabels(["world A", "world B"])
    d.set_ylabel("hazard MAE beyond support")
    d.legend(loc="upper right", ncol=2, handlelength=0.9, columnspacing=0.8, fontsize=6.6)
    panel_label(d, "d", dx=-0.16, dy=1.1)
    fig.text(0.99, 0.0, "synthetic data with known truth (10 replicates); mechanistic prior carries a deliberate +2 $^\\circ$C $T_p$ error", ha="right", va="bottom", fontsize=6.2, color=OI["grey"], style="italic")
    save(fig, out, "fig4_extrapolation")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", default="docs/audit_v1.json")
    ap.add_argument("--out", default="docs/figures")
    ap.add_argument("--report", default="docs/figure_report.json")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    report: dict = {}
    if Path(a.audit).exists():
        fig_v1_audit(json.loads(Path(a.audit).read_text()), out)
    else:
        print("audit json not found - skipping Fig. 1 (run scripts/audit_v1.py)")
    fig_hydraulics(out)
    fig_extrapolation(out, report)
    Path(a.report).write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
