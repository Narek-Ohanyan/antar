"""Two-phase time to hydraulic failure.

Phase 1 - *time to stomatal closure*: the soil dries until predawn water
potential reaches the species' closure potential (psi_close ~ Psi_gs90).  This
phase carries most of the survival time (Gauthey 2026 summarising Waite et al.
2026; Martin-StPaul et al. 2017).

Phase 2 - *time to critical failure*: with stomata shut, the plant drains its
own water store through the cuticle.  Mass balance per unit leaf area:

    C  d(psi)/dt = - g_min(T) * VPD / P_atm

so the buffer B = C (psi_close - psi_crit)   [mmol m-2]   is consumed at
E_min = g_min(T) VPD/P_atm   [mmol m-2 s-1].  The hydraulic-failure index

    HFI = sum_{closed days in spell} 86400 * E_min,d / B

reaches 1 when the plant has reached psi_crit (P88 for angiosperms, P50 for
conifers).  ``HFI >= 1`` is the model's definition of hydraulic failure.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ..constants import SECONDS_PER_DAY
from .gmin import gmin_temperature
from .vulnerability import psi_at_plc


@dataclass(frozen=True)
class Traits:
    p50: float = -3.0           # MPa
    slope: float = 40.0         # % MPa-1 at P50
    psi_close: float = -2.5     # MPa, predawn potential at ~stomatal closure (Psi_gs90)
    capacitance: float = 30000.0  # mmol m-2 (leaf area) MPa-1, whole-tree effective (~0.5 kg m-2 MPa-1); seedlings are far lower
    g25: float = 3.0            # mmol m-2 s-1, residual conductance at 25 C
    tp: float = 38.0            # C, phase-transition temperature
    lethal_plc: float = 88.0    # 88 angiosperm, 50 conifer

    @property
    def psi_crit(self) -> float:
        return float(psi_at_plc(self.lethal_plc, self.p50, self.slope))

    @property
    def buffer(self) -> float:
        """B = C (psi_close - psi_crit), mmol m-2 (0 if closure is already beyond the lethal limit)."""
        return max(self.capacitance * (self.psi_close - self.psi_crit), 0.0)

    def with_updates(self, **kw) -> "Traits":
        return replace(self, **kw)


def t_crit_days_constant(traits: Traits, vpd_kpa: float, p_atm_kpa: float, t_c: float) -> float:
    """Closed-form phase-2 duration (days) for constant VPD and temperature."""
    e_min = float(gmin_temperature(t_c, traits.g25, traits.tp)) * vpd_kpa / p_atm_kpa   # mmol m-2 s-1
    return traits.buffer / (e_min * SECONDS_PER_DAY)


def simulate_two_phase(psi_soil, t_c_leaf, vpd24, p_atm_kpa, traits: Traits, refill_fraction: float = 1.0):
    """Run the two-phase index over a daily series (arrays of equal length).

    Parameters
    ----------
    psi_soil : predawn ~ soil water potential (MPa)
    t_c_leaf : leaf/air temperature driving g_min (C), typically T_max
    vpd24    : 24 h VPD (kPa)
    refill_fraction : fraction of the accumulated deficit recovered per open day (1 = full refill)

    Returns dict(hfi, closed, day_first_closed, day_first_failure, hfi_max).
    """
    psi_soil = np.asarray(psi_soil, dtype=float)
    t_c_leaf = np.asarray(t_c_leaf, dtype=float)
    vpd24 = np.asarray(vpd24, dtype=float)
    n = psi_soil.size
    buf = traits.buffer
    e_min_day = SECONDS_PER_DAY * gmin_temperature(t_c_leaf, traits.g25, traits.tp) * vpd24 / p_atm_kpa
    hfi = np.zeros(n)
    closed = psi_soil <= traits.psi_close
    acc = 0.0
    for d in range(n):
        if closed[d]:
            acc += e_min_day[d]
        else:
            acc *= (1.0 - refill_fraction)
        hfi[d] = acc / buf if buf > 0 else np.inf if closed[d] else 0.0
    first_closed = int(np.argmax(closed)) if closed.any() else -1
    fail = np.where(hfi >= 1.0)[0]
    return {
        "hfi": hfi,
        "closed": closed,
        "day_first_closed": first_closed,
        "day_first_failure": int(fail[0]) if fail.size else -1,
        "hfi_max": float(hfi.max()) if n else 0.0,
    }
