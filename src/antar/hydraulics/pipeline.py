"""XYLEM's entry point from TOPOHYDRO output (Fig. 2: Module A -> Module B).

:func:`mechanistic_hazard_for_cell` reads directly off a
:class:`antar.climate.forcing.CellTopoclimate` -- psi_soil_mpa (Phase 1's
closure trigger, Sec. 6.2), T_max as the leaf-temperature proxy (Sec. 6.4's
own documented simplification: "leaf temperature equal to T_max, tested
against a leaf-air offset"), VPD24 (Phase 2's mass balance, Eq. 6.2) and
pressure_kpa -- and runs the Sec. 6.3 / Eq. 6.5 two-level Monte Carlo on it.
This is the only new logic here; everything it calls is already tested
elsewhere in ``antar.hydraulics`` and ``antar.climate``.
"""
from __future__ import annotations

import numpy as np

from ..climate.forcing import CellTopoclimate
from .monte_carlo import two_level_failure_probability
from .twophase import Traits, simulate_two_phase


def mechanistic_hazard_for_cell(
    cell: CellTopoclimate,
    base: Traits,
    hyper_sd: dict,
    individual_sd: dict,
    pet_formulation: str = "pm_fao56",
    refill_fraction: float = 1.0,
    outer_draws: int = 50,
    inner_draws: int = 200,
    seed: int = 0,
) -> np.ndarray:
    """Run XYLEM for one cell-year using one member of TOPOHYDRO's PET ensemble.

    Returns the outer-loop h_mech array (length ``outer_draws``); its spread
    is the framework's trait-knowledge uncertainty (Sec. 6.3).
    """
    if pet_formulation not in cell.psi_soil_mpa:
        raise KeyError(f"'{pet_formulation}' not in this cell's PET ensemble: {sorted(cell.psi_soil_mpa)}")
    psi_soil = cell.psi_soil_mpa[pet_formulation]

    def simulate(traits: Traits) -> float:
        sim = simulate_two_phase(psi_soil, cell.t_max_c, cell.vpd_24h_kpa, cell.pressure_kpa, traits, refill_fraction)
        return sim["hfi_max"]

    return two_level_failure_probability(
        simulate, base, hyper_sd, individual_sd, outer_draws=outer_draws, inner_draws=inner_draws, seed=seed
    )
