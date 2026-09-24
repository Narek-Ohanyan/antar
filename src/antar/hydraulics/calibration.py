"""Pattern-oriented recalibration (concept note Sec. 6.4).

After the trait-hyperparameter priors are narrowed by approximate Bayesian
computation against observed patterns (Boyce index of predicted viability,
rank correlation with tree-ring growth collapses, dieback frequency by
elevation band -- Sec. 6.4), whatever mismatch remains is absorbed by a
two-parameter recalibration on the link scale, fitted per functional group
with hierarchical shrinkage:

    cloglog(h~_mech) = a_g + b_g * cloglog(h_mech),   b_g > 0

:func:`recalibrate_hazard` is that transform. b_g > 0 keeps the map monotone,
so it can only rescale the mechanistic response, never change its shape --
which is exactly why extrapolation outside the training climate keeps its
mechanistic form (Sec. 6.4's own stated reason for choosing this recipe).

The ABC step itself -- ``run a draw on 2000-2024 forcing at reference sites,
accept if summary statistics fall within tolerance of the observed ones`` --
is not implemented here: two of its three acceptance criteria need data this
repository does not yet have (tree-ring chronologies and dieback histories,
blocked on the hyperion data pull) and the third (the Boyce index of
predicted viability) needs MERISTEM's niche model, built later in the engine
order. See IMPLEMENTATION_LOG.md.
"""
from __future__ import annotations

import numpy as np


def recalibrate_hazard(h_mech, a_g: float, b_g: float):
    """cloglog(h~) = a_g + b_g * cloglog(h_mech). Identity when a_g=0, b_g=1."""
    h = np.clip(np.asarray(h_mech, dtype=float), 1e-12, 1 - 1e-12)
    cloglog_h = np.log(-np.log1p(-h))
    eta = a_g + b_g * cloglog_h
    return -np.expm1(-np.exp(np.clip(eta, -30, 30)))
