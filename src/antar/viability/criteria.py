"""Criteria (b) and (c) of the robust refugium on the dense node lattice, and the full conjunction.

The dense grid takes every ``stride``-th CHELSA pixel inside Armenia, so its nodes sit on a regular lattice with gaps outside the country. The buffer index of ``refugia.buffer_index``
(regional minus local exposure, in national standard deviations) is evaluated on that lattice with a radius counted in lattice cells. The area of applicability (Meyer & Pebesma 2021,
``validation.aoa``) is built, group by group, on the climate and water-balance conditions of the nodes where the group's forest is actually mapped; a node, or a scenario member, outside
it is not supported by any observed stand and is flagged, never asserted.
"""
from __future__ import annotations

import numpy as np

from antar.validation.aoa import AreaOfApplicability
from antar.validation.splits import spatial_block_ids
from antar.viability.refugia import buffer_index


def lattice_index(lat, lon, bbox, grid_shape, stride: int):
    """Lattice (row, column) of nodes that sit on the CHELSA pixel centres of a ``grid_shape`` raster covering ``bbox`` = (lon_min, lat_min, lon_max, lat_max), every ``stride``-th pixel."""
    lon_min, lat_min, lon_max, lat_max = bbox
    rows, cols = grid_shape
    chelsa_row = np.rint((lat_max - np.asarray(lat, dtype=float)) / ((lat_max - lat_min) / rows) - 0.5).astype(int)
    chelsa_col = np.rint((np.asarray(lon, dtype=float) - lon_min) / ((lon_max - lon_min) / cols) - 0.5).astype(int)
    if np.any(chelsa_row % stride) or np.any(chelsa_col % stride):
        raise ValueError("nodes do not sit on every stride-th CHELSA pixel of the given grid")
    return chelsa_row // stride, chelsa_col // stride


def exposure_composite(columns, reference):
    """Mean of the standardised exposure variables (CWD, water stress, heat): each column minus its reference mean over the reference sd. ``reference`` is a list of (mean, sd)."""
    z = [(np.asarray(c, dtype=float) - m) / s for c, (m, s) in zip(columns, reference)]
    return np.mean(z, axis=0)


def node_buffer_index(exposure, ii, jj, lattice_shape, radius_cells: int):
    """Buffer index of every node: the node's exposure is placed on the lattice (NaN where there is no node), ``buffer_index`` is evaluated, and each node's value is read back."""
    grid = np.full(lattice_shape, np.nan)
    grid[ii, jj] = exposure
    return buffer_index(grid, radius_cells=radius_cells)[ii, jj]


def class_fractions(classes, classes_by_group: dict):
    """Share of the pixels of a window that fall in each group's ecosystem classes; ``classes`` is an array of class codes."""
    arr = np.asarray(classes)
    return {g: float(np.isin(arr, list(cl)).mean()) for g, cl in classes_by_group.items()}


def group_aoa(train_X, train_lat, train_lon, X_new, block_deg: float = 1.0, min_train: int = 20):
    """Area of applicability of one group's climate/water-balance model, from the nodes where the group's forest is mapped.

    Returns ``{"status", "n_train", "threshold", "inside" (bool per row of X_new), "di" (dissimilarity index)}``. With fewer than ``min_train`` training nodes (or fewer than two
    spatial blocks) no area of applicability can be established: status is "insufficient_training_sample" and nothing is inside.
    """
    train_X = np.asarray(train_X, dtype=float)
    X_new = np.asarray(X_new, dtype=float)
    ok = np.all(np.isfinite(train_X), axis=1)
    blocks = spatial_block_ids(np.asarray(train_lon)[ok], np.asarray(train_lat)[ok], block_deg) if ok.any() else np.array([])
    if ok.sum() < min_train or np.unique(blocks).size < 2:
        return {"status": "insufficient_training_sample", "n_train": int(ok.sum()), "threshold": None,
                "inside": np.zeros(len(X_new), dtype=bool), "di": np.full(len(X_new), np.nan)}
    aoa = AreaOfApplicability().fit(train_X[ok], folds=blocks)
    finite = np.all(np.isfinite(X_new), axis=1)
    di = np.full(len(X_new), np.nan)
    di[finite] = aoa.dissimilarity_index(X_new[finite])
    return {"status": "ok", "n_train": int(ok.sum()), "threshold": float(aoa.threshold_), "inside": finite & (di <= aoa.threshold_), "di": di}


def robust_conjunction(criterion_a, buffer, inside_aoa):
    """Criteria (a), (b) B > 0 and (c) inside the area of applicability, all required; a NaN buffer index fails (b)."""
    b = np.asarray(buffer, dtype=float)
    return np.asarray(criterion_a, dtype=bool) & np.isfinite(b) & (b > 0.0) & np.asarray(inside_aoa, dtype=bool)
