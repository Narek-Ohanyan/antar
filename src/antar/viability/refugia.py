"""Climate-change refugia, defined operationally and robustly.

A cell is a *robust refugium for species s at horizon T* when

 (a) P_c[ V_s(c) >= v* ] >= rho           (viable in at least a fraction rho of the ensemble),
 (b) B_i > 0                              (topographically decoupled from the regional exposure),
 (c) it lies inside the model's area of applicability (else flagged, not asserted),
 (d) it is connected to / large enough for a functional stand (checked downstream on the raster).

Buffer index  B_i = (E_reg,i - E_i) / sd(E)  where E is a standardised composite exposure
(CWD, growing-season VPD, heat) and E_reg,i its regional (neighbourhood) mean.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter


def buffer_index(exposure_grid, radius_cells: int = 25):
    """Regional-minus-local exposure in units of the national SD (2-D grid; NaN-safe)."""
    e = np.asarray(exposure_grid, dtype=float)
    valid = np.isfinite(e)
    filled = np.where(valid, e, 0.0)
    size = 2 * radius_cells + 1
    num = uniform_filter(filled, size=size, mode="nearest")
    den = uniform_filter(valid.astype(float), size=size, mode="nearest")
    reg = np.divide(num, den, out=np.full_like(e, np.nan), where=den > 0)
    sd = np.nanstd(e)
    return (reg - e) / sd


def robust_refugium(v_ensemble, v_star: float = 0.6, rho: float = 0.8, buffer_index_grid=None, inside_aoa=None):
    """R_ik (Eq. 8.6): criterion (a) alone, or the full conjunction with (b)/(c) when supplied.

    ``v_ensemble``: (n_members, ...) viability array -> criterion (a),
    ``P_omega[V >= v_star] >= rho``. ``buffer_index_grid`` (this module's own
    :func:`buffer_index` output) adds criterion (b), ``B_i > 0``.
    ``inside_aoa`` (e.g. ``antar.validation.aoa.AreaOfApplicability.inside``)
    adds criterion (c). Both default to ``None`` (criterion skipped) so a
    caller with only viability draws still gets (a) alone; the note's full
    Eq. 8.6 definition of a *robust refugium* requires supplying both --
    criterion (d), stand connectivity/size, is checked downstream on the
    raster and is out of scope for this function.
    """
    v = np.asarray(v_ensemble, dtype=float)
    mask = np.mean(v >= v_star, axis=0) >= rho
    if buffer_index_grid is not None:
        mask = mask & (np.asarray(buffer_index_grid, dtype=float) > 0.0)
    if inside_aoa is not None:
        mask = mask & np.asarray(inside_aoa, dtype=bool)
    return mask


def refugium_score(v_ensemble, lam: float = 0.5):
    """Risk-averse score: median viability minus lam * interquartile range across the ensemble."""
    v = np.asarray(v_ensemble, dtype=float)
    q25, q50, q75 = np.percentile(v, [25, 50, 75], axis=0)
    return q50 - lam * (q75 - q25)
