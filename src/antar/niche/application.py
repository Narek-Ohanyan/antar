"""Putting the fitted species niche models to work at grid nodes, under the present climate and under each scenario member.

A presence-background model gives a relative suitability, not a probability of occurrence, so it is turned into a yes/no by the usual presence threshold: a site is within a species' niche
when its score is at least the ``omission`` quantile of the scores at the species' own occurrence records (10 % omission: one record in ten falls below it). A group is niche-supported at a
site when at least one of its species models that passed the quality gate says so. A species with no usable model gives no support anywhere, which is different from "suitable everywhere":
the options of a group without niche support are not ranked.

Scenario members need the niche predictors under the member's climate. Winter minimum temperature and growing degree days come from CHELSA-BIOCLIM+ for the period and climate model of
the member; the water deficit is the member's own TOPOHYDRO result; soils do not change. CHELSA-BIOCLIM+ has no future vapour-pressure deficit, so it is scaled from its present value by the
rise in saturation vapour pressure for the member's warming at constant relative humidity (Magnus formula), which is a lower bound on the rise in a drying continental climate.
"""
from __future__ import annotations

import numpy as np

PERIOD_OF_HORIZON = {2050: "2041-2070", 2080: "2071-2100", 2100: "2071-2100"}   # CHELSA-BIOCLIM+ has three 30-year periods; 2080 and 2100 share the last one
GCM_NAME = {"gfdl-esm4": "GFDL-ESM4", "ipsl-cm6a-lr": "IPSL-CM6A-LR", "mpi-esm1-2-hr": "MPI-ESM1-2-HR", "mri-esm2-0": "MRI-ESM2-0", "ukesm1-0-ll": "UKESM1-0-LL"}


def chelsa_key(horizon: int, gcm: str, ssp: str) -> str:
    """The CHELSA-BIOCLIM+ layer key of a scenario member, e.g. ``2041-2070|GFDL-ESM4|ssp126``."""
    return f"{PERIOD_OF_HORIZON[int(horizon)]}|{GCM_NAME[gcm.lower()]}|{ssp.lower()}"


def saturation_vapour_pressure_kpa(t_c):
    """Magnus formula: 0.6108 exp(17.27 T / (T + 237.3)) kPa, T in degrees Celsius."""
    t = np.asarray(t_c, dtype=float)
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def vpd_after_warming(vpd_present, t_present_c, warming_c):
    """Vapour-pressure deficit after a warming of ``warming_c`` degrees, relative humidity held constant: VPD scales with the saturation vapour pressure, e_s(T + dT) / e_s(T)."""
    return np.asarray(vpd_present, dtype=float) * saturation_vapour_pressure_kpa(np.asarray(t_present_c, dtype=float) + warming_c) / saturation_vapour_pressure_kpa(t_present_c)


def logistic_score(coef, intercept, X):
    """Fitted suitability: 1 / (1 + exp(-(intercept + X . coef))). Rows with a missing predictor are NaN."""
    X = np.atleast_2d(np.asarray(X, dtype=float))
    z = float(intercept) + X @ np.asarray(coef, dtype=float)
    return 1.0 / (1.0 + np.exp(-z))


def presence_threshold(scores_at_presence, omission: float = 0.10) -> float:
    """The score below which a share ``omission`` of the species' own records falls (the minimum-training-presence rule relaxed to ``omission``)."""
    s = np.asarray(scores_at_presence, dtype=float)
    s = s[np.isfinite(s)]
    if s.size == 0:
        raise ValueError("no occurrence scores to set a threshold from")
    return float(np.quantile(s, omission))


def passes_gate(status: str, boyce_mean, min_boyce: float) -> bool:
    """A species model is used only if it was fitted and its mean spatial-block Boyce index is at least ``min_boyce``."""
    return status == "fitted" and boyce_mean is not None and np.isfinite(boyce_mean) and float(boyce_mean) >= min_boyce


def group_support(suitable_by_species: dict, passing: list):
    """Boolean support of a group from its species: True where any passing species is within its niche. ``suitable_by_species`` maps species to a boolean array; returns None when no
    species of the group passed (no niche model, so no claim either way)."""
    arrays = [np.asarray(suitable_by_species[s], dtype=bool) for s in passing if s in suitable_by_species]
    if not arrays:
        return None
    return np.any(np.vstack(arrays), axis=0)
