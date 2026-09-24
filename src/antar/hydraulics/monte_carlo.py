"""Mechanistic hazard = probability that hydraulic failure occurs under parameter uncertainty.

    h_mech(i, s, y) = (1/M) * sum_m 1[ HFI^{(m)}_{i,s,y} >= 1 ],   theta^{(m)} ~ p(theta | trait database)

This converts a deterministic index into a probability that already carries
trait and soil-parameter uncertainty (no arbitrary fragility curve).

Sec. 6.3 specifies *two* nested levels of sampling, not one: "an outer loop of
50 draws takes hyper-parameters (group means and spreads of P50, S,
psi_close, C, g25, Tp) from posteriors built on XFT, TRY and the g_min
compilation ... an inner loop of 200 draws samples individuals." The
mechanistic hazard :func:`failure_probability` is the Eq. 6.5 average over one
inner loop; :func:`two_level_failure_probability` repeats that for each of the
outer draws and returns the resulting array, whose *spread* is the engine's
knowledge uncertainty -- also stated explicitly in Sec. 6.3.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from .twophase import Traits


def sample_traits(rng: np.random.Generator, base: Traits, sd: dict[str, float], n: int) -> list[Traits]:
    """Draw ``n`` trait sets: normal perturbations (sd per field) with sensible positivity constraints."""
    out = []
    for _ in range(n):
        upd = {}
        for k, s in sd.items():
            v = getattr(base, k) + rng.normal(0.0, s)
            if k in ("capacitance", "g25", "slope"):
                v = max(v, 0.05 * getattr(base, k))
            upd[k] = v
        t = base.with_updates(**upd)
        # keep closure potential above (less negative than) P50 so that the buffer is well defined
        if t.psi_close < t.p50:
            t = t.with_updates(psi_close=t.p50)
        out.append(t)
    return out


def failure_probability(simulate: Callable[[Traits], float], traits_draws: list[Traits]) -> float:
    """``simulate`` maps a Traits draw to the annual max HFI."""
    hfi = np.array([simulate(t) for t in traits_draws])
    return float(np.mean(hfi >= 1.0))


def sample_trait_hyperparameters(rng: np.random.Generator, base: Traits, hyper_sd: dict[str, float]) -> Traits:
    """One outer-loop draw of the trait-group *hyperparameters* (Sec. 6.3).

    Same perturbation and positivity-guard logic as :func:`sample_traits`, but
    representing knowledge uncertainty in the group-level posterior (built on
    XFT/TRY/g_min-compilation data, later updated with local measurements),
    not individual-to-individual variation.
    """
    upd = {}
    for k, s in hyper_sd.items():
        v = getattr(base, k) + rng.normal(0.0, s)
        if k in ("capacitance", "g25", "slope"):
            v = max(v, 0.05 * getattr(base, k))
        upd[k] = v
    hyper = base.with_updates(**upd)
    if hyper.psi_close < hyper.p50:
        hyper = hyper.with_updates(psi_close=hyper.p50)
    return hyper


def two_level_failure_probability(
    simulate: Callable[[Traits], float],
    base: Traits,
    hyper_sd: dict[str, float],
    individual_sd: dict[str, float],
    outer_draws: int = 50,
    inner_draws: int = 200,
    seed: int = 0,
) -> np.ndarray:
    """The full two-level Monte Carlo of Sec. 6.3 / Eq. 6.5.

    Returns an array of length ``outer_draws``: h_mech at each hyperparameter
    draw (itself the Eq. 6.5 average over ``inner_draws`` individuals). The
    spread of the returned array across the outer loop is the framework's
    trait-knowledge uncertainty, to be reported separately from individual
    variation, not averaged away.
    """
    rng = np.random.default_rng(seed)
    h = np.empty(outer_draws)
    for m in range(outer_draws):
        hyper = sample_trait_hyperparameters(rng, base, hyper_sd)
        individuals = sample_traits(rng, hyper, individual_sd, inner_draws)
        h[m] = failure_probability(simulate, individuals)
    return h
