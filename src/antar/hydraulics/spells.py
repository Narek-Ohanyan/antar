"""Stress-spell summaries: the emulator's climate-side input features (Sec. 6.3, "Emulation").

"A Latin hypercube over stress-spell summaries (longest closed spell, mean
VPD24 and T during it, number of closed days) and trait-group descriptors is
passed through the full engine" -- :func:`stress_spell_summary` computes those
four summaries from a closed/open sequence (as returned by
:func:`antar.hydraulics.twophase.simulate_two_phase`'s ``closed`` array) plus
the daily VPD24/temperature series that drove it.
"""
from __future__ import annotations

import numpy as np


def longest_true_run(mask) -> tuple[int, int, int]:
    """(length, start, end) of the longest contiguous run of True in a 1-D boolean array.

    ``end`` is exclusive, so ``mask[start:end]`` is the run itself. (0, -1, -1)
    if ``mask`` is never True.
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return 0, -1, -1
    padded = np.concatenate(([0], mask.astype(int), [0]))
    edges = np.diff(padded)
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    lengths = ends - starts
    i = int(np.argmax(lengths))
    return int(lengths[i]), int(starts[i]), int(ends[i])


def stress_spell_summary(closed, vpd24, t_c_leaf) -> dict:
    """Concept-note Sec. 6.3 emulator features for one cell-year.

    Returns ``{longest_closed_spell_days, mean_vpd24_during_spell,
    mean_t_during_spell, n_closed_days}``. The mean VPD/T are over the days
    within the *longest* closed spell specifically (the note's "during it"
    refers to the longest spell), while ``n_closed_days`` counts closed days
    over the whole record.
    """
    closed = np.asarray(closed, dtype=bool)
    vpd24 = np.asarray(vpd24, dtype=float)
    t_c_leaf = np.asarray(t_c_leaf, dtype=float)
    length, start, end = longest_true_run(closed)
    if length == 0:
        return {
            "longest_closed_spell_days": 0,
            "mean_vpd24_during_spell": 0.0,
            "mean_t_during_spell": 0.0,
            "n_closed_days": 0,
        }
    return {
        "longest_closed_spell_days": length,
        "mean_vpd24_during_spell": float(np.mean(vpd24[start:end])),
        "mean_t_during_spell": float(np.mean(t_c_leaf[start:end])),
        "n_closed_days": int(np.sum(closed)),
    }
