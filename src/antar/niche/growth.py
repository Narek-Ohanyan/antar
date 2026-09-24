"""Height growth with climate-modified *physiological age*.

Chapman-Richards:  H(a) = H* (1 - exp(-k a))^p

Time-varying climate is handled by letting physiological age advance by
phi_y in (0, 1] each calendar year (phi = 1 in an unstressed year):

    a_hat_{y+1} = a_hat_y + phi_y,     phi_y = f_W f_T f_F f_tl   (Eq. 8.2)

so the closed-form curve is retained and stress years simply delay development.
f_tl (:func:`thermal_treeline_modifier`) is a fourth, physically distinct
modifier from f_T (:func:`gdd_modifier`): above a species' realised range,
growth *stops* (a meristematic/carbon-sink-limitation ceiling), it does not
merely slow for lack of accumulated heat. The concept note is explicit that
these two must stay separate, independently testable functions, never merged
or aliased -- so they are.
"""
from __future__ import annotations

import numpy as np
from scipy.special import expit


def chapman_richards(age, h_star, k, p: float = 1.5):
    age = np.maximum(np.asarray(age, dtype=float), 0.0)
    return h_star * (1.0 - np.exp(-k * age)) ** p


def physiological_age(phi_series, a0: float = 0.0):
    """Cumulative physiological age after each year, given annual growth-effectiveness phi in (0, 1]."""
    return a0 + np.cumsum(np.clip(np.asarray(phi_series, dtype=float), 0.0, 1.0))


def water_modifier(cwd_mm, cwd50_mm, m: float = 3.0):
    """f_W = 1 / (1 + (CWD/CWD50)^m): equals 0.5 at CWD = CWD50."""
    return 1.0 / (1.0 + (np.asarray(cwd_mm, dtype=float) / cwd50_mm) ** m)


def gdd_modifier(gdd_y, gdd_ref: float):
    """f_T: the growing-degree-day (Liebig-type) growth modifier (Sec. 8.1; "as in 3-PG",
    Landsberg & Waring 1997).

    No closed form is given beyond naming it a growing-degree-day modifier;
    implemented as a Michaelis-Menten/Monod saturating ramp -- the increasing-
    saturating counterpart of this module's own :func:`water_modifier` (0 at
    GDD=0, 0.5 at GDD=gdd_ref, -> 1 as GDD grows). A documented choice, not a
    quoted formula; ``gdd_ref`` (species/site-specific) has no default.
    """
    gdd = np.asarray(gdd_y, dtype=float)
    return gdd / (gdd + gdd_ref)


def late_frost_modifier(n_frost_days, n_frost50: float, m: float = 3.0):
    """f_F: the late-frost (Liebig-type) growth modifier (Sec. 8.1).

    No closed form is given beyond naming it a late-frost modifier;
    implemented with the same saturating form as :func:`water_modifier`
    (1 at zero frost days, 0.5 at n_frost50, -> 0 as frost days grow) applied
    to :func:`antar.climate.indices.late_frost_days` instead of CWD. A
    documented choice, not a quoted formula.
    """
    return 1.0 / (1.0 + (np.asarray(n_frost_days, dtype=float) / n_frost50) ** m)


# Koerner & Paulsen (2004): pan-biome tree-form growing-season thresholds, cited in the
# concept note as T_tl ~ 6.4-6.5 degC and LGS_min ~ 90-94 days; midpoints used as defaults.
_KORNER_PAULSEN_T_TL_C = 6.45
_KORNER_PAULSEN_LGS_MIN_DAYS = 92.0


def thermal_treeline_modifier(gst, lgs, s_t: float, s_l: float,
                               t_tl_c: float = _KORNER_PAULSEN_T_TL_C,
                               lgs_min_days: float = _KORNER_PAULSEN_LGS_MIN_DAYS):
    """f_tl (Eq. 8.3): a smooth version of the Koerner-Paulsen treeline rule.

    ``gst``/``lgs`` are the growing-season mean temperature and length
    (:func:`antar.climate.indices.growing_season_mean_temperature` /
    :func:`.growing_season_length`). ``t_tl_c``/``lgs_min_days`` default to the
    concept note's own cited values; ``s_t``/``s_l`` (the transition widths)
    are explicitly "fitted, not assumed, from growth-collapse tree-ring
    chronologies" (Sec. 8.1) and have no default. As f_tl -> 0, physiological
    age stops advancing regardless of water status -- this is why it
    multiplies phi_y directly rather than folding into :func:`gdd_modifier`.
    """
    gst = np.asarray(gst, dtype=float)
    lgs = np.asarray(lgs, dtype=float)
    return expit((gst - t_tl_c) / s_t) * expit((lgs - lgs_min_days) / s_l)


def potential_treeline_elevation(gst_ref_c, z_ref_m, gamma_k_per_m,
                                  t_tl_c: float = _KORNER_PAULSEN_T_TL_C):
    """Diagnostic potential treeline elevation z_tl(x, t) (Sec. 8.1): invert the
    Sec. 5.1 lapse relation T(z) = T_ref + Gamma(z - z_ref) for the elevation at
    which growing-season mean temperature equals T_tl.

    A statement of the model's climatic ceiling, to be recomputed per ensemble
    member and reported directly -- explicitly *not* a forecast of where forest
    will actually stand: "realised treelines lag their climatic potential by
    decades" (Sec. 8.1).
    """
    gamma = np.asarray(gamma_k_per_m, dtype=float)
    return z_ref_m + (t_tl_c - np.asarray(gst_ref_c, dtype=float)) / gamma


def height_trajectory(h_star, k, p, phi_series):
    return chapman_richards(physiological_age(phi_series), h_star, k, p)


def years_to_height(h_target, h_star, k, p: float = 1.5):
    """Physiological age at which H reaches h_target (inf if h_target >= h_star)."""
    if h_target >= h_star:
        return np.inf
    return -np.log(1.0 - (h_target / h_star) ** (1.0 / p)) / k


def probability_reaches_height(phi_series_ensemble, h_star_ensemble, k_ensemble, p_ensemble,
                                h_min: float, horizon_years: int):
    """P[H_T >= H_min] = P[sum_{y<=T} phi_y >= a*], evaluated over an ensemble of
    scenario/parameter draws (Sec. 8.1), a* = :func:`years_to_height`.

    ``phi_series_ensemble`` (M, >=horizon_years): each member's annual
    growth-effectiveness series phi_y = f_W f_T f_F f_tl; ``h_star_ensemble``/
    ``k_ensemble``/``p_ensemble`` (M,): that member's Chapman-Richards
    parameters. Returns the fraction of the M members reaching H_min by
    ``horizon_years``.
    """
    phi = np.clip(np.asarray(phi_series_ensemble, dtype=float)[:, :horizon_years], 0.0, 1.0)
    a_t = phi.sum(axis=1)
    h_star_ensemble = np.asarray(h_star_ensemble, dtype=float)
    k_ensemble = np.asarray(k_ensemble, dtype=float)
    p_ensemble = np.asarray(p_ensemble, dtype=float)
    a_star = np.array([
        years_to_height(h_min, hs, k, p) for hs, k, p in zip(h_star_ensemble, k_ensemble, p_ensemble)
    ])
    return float(np.mean(a_t >= a_star))
