"""Discrete-time hazard learners.

* :class:`CloglogGLM` - the interpretable backbone.  h = 1 - exp(-exp(eta)) so that
  cumulative hazards add over years and survival is S(T) = exp(-sum_y exp(eta_y)).
* :func:`monotone_hgb`  - gradient boosting with monotone (physically signed) constraints.
* :func:`stack_nnls`    - non-negative stacking of out-of-fold predictions on the link scale.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import nnls
from sklearn.ensemble import HistGradientBoostingClassifier


def cloglog_inverse(eta):
    return -np.expm1(-np.exp(np.clip(eta, -30, 30)))


def cloglog(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-12, 1 - 1e-12)
    return np.log(-np.log1p(-p))


def _soft_threshold(x, lam):
    return np.sign(x) * np.maximum(np.abs(x) - lam, 0.0)


def _weighted_penalized_wls(A, z, w, l1: float, l2: float, fit_intercept: bool, beta_init,
                             cd_max_iter: int = 200, cd_tol: float = 1e-8):
    """Solve  argmin_beta  sum_i w_i (z_i - A_i.beta)^2 + l2||beta_s||_2^2 + l1||beta_s||_1
    (beta_s = beta excluding the unpenalised intercept), the weighted-least-squares
    subproblem each IRLS step needs to solve for a fixed working response ``z``.

    l1 == 0 uses the closed-form ridge normal equations (exact, matches the
    pre-elastic-net code path). l1 > 0 has no closed form; cyclic coordinate
    descent with soft-thresholding is the standard algorithm for this exact
    problem (Friedman, Hastie & Tibshirani 2010), applied here to the current
    IRLS working response rather than to y directly (the IRLS-outer / CD-inner
    scheme this generalises to for GLMs, same reference).
    """
    p = A.shape[1]
    if l1 == 0.0:
        pen = l2 * np.eye(p)
        if fit_intercept:
            pen[0, 0] = 0.0
        H = A.T @ (A * w[:, None]) + pen
        return np.linalg.solve(H, A.T @ (w * z))

    beta = beta_init.copy()
    resid = z - A @ beta
    col_ss = (A**2 * w[:, None]).sum(axis=0)
    start = 1 if fit_intercept else 0
    for _ in range(cd_max_iter):
        beta_old = beta.copy()
        if fit_intercept:
            r0 = resid + beta[0]           # A[:, 0] == 1
            new0 = np.sum(w * r0) / np.sum(w)
            resid += beta[0] - new0
            beta[0] = new0
        for j in range(start, p):
            aj = A[:, j]
            rj = resid + aj * beta[j]
            rho = np.sum(w * aj * rj)
            denom = col_ss[j] + l2
            new_j = _soft_threshold(rho, l1 / 2.0) / denom if denom > 0 else 0.0
            resid += aj * (beta[j] - new_j)
            beta[j] = new_j
        if np.max(np.abs(beta - beta_old)) < cd_tol:
            break
    return beta


class CloglogGLM:
    """Elastic-net-penalised binomial GLM with complementary log-log link, fitted by IRLS
    (l1 == 0, the default, reduces to the exact closed-form ridge solve)."""

    def __init__(self, l1: float = 0.0, l2: float = 1e-4, max_iter: int = 100, tol: float = 1e-8,
                 fit_intercept: bool = True, cd_max_iter: int = 200, cd_tol: float = 1e-8):
        self.l1 = l1
        self.l2 = l2
        self.max_iter = max_iter
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.cd_max_iter = cd_max_iter
        self.cd_tol = cd_tol
        self.coef_ = None
        self.intercept_ = 0.0

    def _design(self, X):
        X = np.asarray(X, dtype=float)
        return np.column_stack([np.ones(len(X)), X]) if self.fit_intercept else X

    def fit(self, X, y, offset=None, sample_weight=None):
        A = self._design(X)
        y = np.asarray(y, dtype=float)
        off = np.zeros(len(y)) if offset is None else np.asarray(offset, dtype=float)
        sw = np.ones(len(y)) if sample_weight is None else np.asarray(sample_weight, dtype=float)
        p = A.shape[1]
        beta = np.zeros(p)
        if self.fit_intercept:
            beta[0] = float(cloglog(np.clip(y.mean(), 1e-4, 1 - 1e-4)))
        for _ in range(self.max_iter):
            eta = A @ beta + off
            e = np.exp(np.clip(eta, -30, 30))
            mu = -np.expm1(-e)
            mu = np.clip(mu, 1e-10, 1 - 1e-10)
            dmu = np.maximum((1 - mu) * e, 1e-10)
            w = sw * dmu**2 / (mu * (1 - mu))
            z = eta - off + (y - mu) / dmu
            b_new = _weighted_penalized_wls(
                A, z, w, self.l1, self.l2, self.fit_intercept, beta, self.cd_max_iter, self.cd_tol
            )
            if np.max(np.abs(b_new - beta)) < self.tol:
                beta = b_new
                break
            beta = b_new
        if self.fit_intercept:
            self.intercept_, self.coef_ = float(beta[0]), beta[1:]
        else:
            self.intercept_, self.coef_ = 0.0, beta
        return self

    def linear_predictor(self, X, offset=None):
        off = 0.0 if offset is None else np.asarray(offset, dtype=float)
        return self.intercept_ + np.asarray(X, dtype=float) @ self.coef_ + off

    def predict_hazard(self, X, offset=None):
        return cloglog_inverse(self.linear_predictor(X, offset))


def monotone_hgb(n_features: int, increasing=(), decreasing=(), **kw) -> HistGradientBoostingClassifier:
    """Gradient boosting classifier with monotone constraints (feature index lists)."""
    cst = np.zeros(n_features, dtype=int)
    for j in increasing:
        cst[j] = 1
    for j in decreasing:
        cst[j] = -1
    params = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=1.0, random_state=0)
    params.update(kw)
    return HistGradientBoostingClassifier(monotonic_cst=cst.tolist(), **params)


def stack_nnls(oof_link_preds, y_link_target):
    """Non-negative weights (sum to 1) combining base learners on the link scale.

    ``oof_link_preds`` (n, m): out-of-fold linear predictors from blocked CV.
    ``y_link_target``  (n,)  : a smoothed empirical link target (e.g. cloglog of a
                                kernel-smoothed event rate) or pseudo-response.
    """
    w, _ = nnls(np.asarray(oof_link_preds, dtype=float), np.asarray(y_link_target, dtype=float))
    s = w.sum()
    return w / s if s > 0 else np.full(w.shape, 1.0 / len(w))
