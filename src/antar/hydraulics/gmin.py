"""Minimum (residual, cuticular) leaf conductance with a temperature phase transition.

Bi-phasic response (Cochard 2019; Schuster et al. 2016; Duursma et al. 2019):
    g_min(T) = g25 * Q10a ** ((T - 25)/10)                          T <= Tp
    g_min(T) = g25 * Q10a ** ((Tp - 25)/10) * Q10b ** ((T - Tp)/10)  T >  Tp
Cochard (2019) reports Q10a ~ 1.2 and Q10b ~ 4.8 with Tp species-specific
(order 30-45 C).  Units of g_min are inherited from g25 (mmol m-2 s-1).
"""
from __future__ import annotations

import numpy as np


def gmin_temperature(t_c, g25, tp_c, q10a: float = 1.2, q10b: float = 4.8):
    t = np.asarray(t_c, dtype=float)
    below = g25 * q10a ** ((t - 25.0) / 10.0)
    at_tp = g25 * q10a ** ((tp_c - 25.0) / 10.0)
    above = at_tp * q10b ** ((t - tp_c) / 10.0)
    return np.where(t <= tp_c, below, above)
