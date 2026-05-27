from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import linalg, stats


def as_1d(x, name: str) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.ndim != 1:
        arr = np.ravel(arr)
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty.")
    return arr


def as_2d(x, name: str, ncol: int | None = None) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional array.")
    if ncol is not None and arr.shape[1] != ncol:
        raise ValueError(f"{name} must have exactly {ncol} columns.")
    return arr


def check_lengths(y: np.ndarray, *arrays: np.ndarray) -> None:
    n = len(y)
    for arr in arrays:
        if len(arr) != n:
            raise ValueError("Input vectors and rows of matrix inputs must have the same length.")


def complete_cases(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        arr = np.asarray(arr)
        try:
            finite = np.isfinite(arr.astype(float))
        except (TypeError, ValueError):
            finite = ~pd.isna(arr)
        if arr.ndim == 1:
            mask &= finite
        else:
            mask &= np.all(finite, axis=1)
    return mask


def prepare_covs_eff(covs_eff, n: int) -> tuple[np.ndarray | None, list[str]]:
    if covs_eff is None:
        return None, []
    if isinstance(covs_eff, pd.DataFrame):
        names = [str(col) for col in covs_eff.columns]
        arr = covs_eff.to_numpy(dtype=float)
    else:
        arr = np.asarray(covs_eff, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        if arr.ndim != 2:
            raise ValueError("covs_eff must be None, a numeric vector, matrix, or data frame.")
        names = [f"z.{j + 1}" for j in range(arr.shape[1])]
    if arr.shape[0] != n:
        raise ValueError("covs_eff must have the same number of observations as Y.")
    if arr.shape[1] < 1:
        raise ValueError("covs_eff must contain at least one covariate column.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("covs_eff must contain finite numeric covariates.")
    return arr, names


def validate_fitmethod(fitmethod: str) -> str:
    fitmethod = str(fitmethod).lower()
    if fitmethod not in {"joint", "separate"}:
        raise ValueError("fitmethod must be 'joint' or 'separate'.")
    return fitmethod


def validate_covs_rank_options(covs_drop: bool, covs_tol: float) -> tuple[bool, float]:
    if not isinstance(covs_drop, (bool, np.bool_)):
        raise ValueError("covs_drop must be True or False.")
    covs_tol = float(covs_tol)
    if not np.isfinite(covs_tol) or covs_tol <= 0:
        raise ValueError("covs_tol must be a positive finite numeric value.")
    return bool(covs_drop), covs_tol


@dataclass
class CovariateRankInfo:
    supplied: int
    used: int
    redundant: list[str]
    dropped: list[str]
    rank_deficient: bool
    drop: bool
    tol: float


class CovariateRankTracker:
    def __init__(self, names: list[str], covs_drop: bool, covs_tol: float) -> None:
        self.names = list(names)
        self.supplied = len(names)
        self.used = len(names)
        self.redundant: list[str] = []
        self.dropped: list[str] = []
        self.rank_deficient = False
        self.covs_drop = bool(covs_drop)
        self.covs_tol = float(covs_tol)
        self._warned = False

    def update(self, keep: np.ndarray) -> None:
        keep_set = set(int(i) for i in keep)
        dropped = [name for i, name in enumerate(self.names) if i not in keep_set]
        self.rank_deficient = True
        self.used = min(self.used, len(keep_set))
        for name in dropped:
            if name not in self.redundant:
                self.redundant.append(name)
            if self.covs_drop and name not in self.dropped:
                self.dropped.append(name)
        if not self._warned:
            action = "dropping redundant covariate columns" if self.covs_drop else "using a generalized inverse"
            warnings.warn(
                "covs_eff is rank deficient after residualizing on the local polynomial basis; "
                f"{action} for numerical stability.",
                RuntimeWarning,
                stacklevel=3,
            )
            self._warned = True

    def summary(self) -> CovariateRankInfo:
        return CovariateRankInfo(
            supplied=self.supplied,
            used=self.used,
            redundant=list(self.redundant),
            dropped=list(self.dropped),
            rank_deficient=self.rank_deficient,
            drop=self.covs_drop,
            tol=self.covs_tol,
        )


def solve_covariate_gamma(
    zwz: np.ndarray,
    zwy: np.ndarray,
    names: list[str],
    tracker: CovariateRankTracker | None,
    covs_drop: bool,
    covs_tol: float,
) -> tuple[np.ndarray, int]:
    zwz = np.asarray(zwz, dtype=float)
    zwy = np.asarray(zwy, dtype=float)
    p = zwz.shape[0]
    if p == 0:
        return np.empty((0, zwy.shape[1])), 0
    _, r, piv = linalg.qr(zwz, mode="economic", pivoting=True)
    diag = np.abs(np.diag(r))
    threshold = covs_tol * max(1.0, diag[0] if diag.size else 0.0)
    rank = int(np.sum(diag > threshold))
    keep = np.asarray(piv[:rank], dtype=int)

    if rank < p:
        if tracker is not None:
            tracker.update(keep)
        gamma = np.zeros((p, zwy.shape[1]), dtype=float)
        if covs_drop and rank > 0:
            gamma_keep = np.linalg.pinv(zwz[np.ix_(keep, keep)], rcond=1e-20) @ zwy[keep, :]
            gamma[keep, :] = gamma_keep
        elif not covs_drop:
            gamma = np.linalg.pinv(zwz, rcond=1e-20) @ zwy
    else:
        try:
            gamma = np.linalg.solve(zwz, zwy)
        except np.linalg.LinAlgError:
            gamma = np.linalg.pinv(zwz, rcond=1e-20) @ zwy
    return np.asarray(gamma, dtype=float), rank


def normalize_binary(x, name: str) -> np.ndarray:
    arr = np.asarray(x)
    if arr.dtype == bool:
        return arr.astype(bool)
    vals = np.unique(arr[np.isfinite(arr.astype(float))].astype(float))
    if not set(vals.tolist()).issubset({0.0, 1.0}):
        raise ValueError(f"{name} must be logical or contain only 0 and 1.")
    return arr.astype(float).astype(bool)


def validate_order(value, name: str) -> int:
    if value is None:
        raise ValueError(f"{name} must not be None.")
    value = float(value)
    if not np.isfinite(value) or value < 0 or abs(value - round(value)) > np.sqrt(np.finfo(float).eps):
        raise ValueError(f"{name} must be a nonnegative integer.")
    return int(round(value))


def validate_deriv(deriv: Iterable[float], p: int) -> tuple[int, int]:
    arr = np.asarray(tuple(deriv), dtype=float)
    if arr.shape != (2,) or np.any(~np.isfinite(arr)) or np.any(arr < 0):
        raise ValueError("deriv must be a nonnegative integer vector of length 2.")
    if np.any(np.abs(arr - np.round(arr)) > np.sqrt(np.finfo(float).eps)):
        raise ValueError("deriv must be a nonnegative integer vector of length 2.")
    out = tuple(int(v) for v in np.round(arr))
    if sum(out) > p:
        raise ValueError("sum(deriv) must be less than or equal to p.")
    return out


def kernel_weights(u: np.ndarray, kernel: str) -> np.ndarray:
    kernel = kernel.lower()
    if kernel in {"tri", "triangular"}:
        return np.maximum(1.0 - np.abs(u), 0.0) * (np.abs(u) <= 1.0)
    if kernel in {"epa", "epanechnikov"}:
        return 0.75 * (1.0 - u**2) * (np.abs(u) <= 1.0)
    if kernel in {"uni", "uniform"}:
        return 0.5 * (np.abs(u) <= 1.0)
    if kernel in {"gau", "gaussian"}:
        return stats.norm.pdf(u)
    raise ValueError("kernel must be one of tri, epa, uni, or gau.")


def multi_indices_2d(p: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for deg in range(p + 1):
        for ypow in range(deg + 1):
            xpow = deg - ypow
            out.append((xpow, ypow))
    return out


def basis_2d(centered: np.ndarray, p: int) -> np.ndarray:
    idx = multi_indices_2d(p)
    x1 = centered[:, 0]
    x2 = centered[:, 1]
    return np.column_stack([(x1**a) * (x2**b) for a, b in idx])


def basis_1d(distance: np.ndarray, p: int) -> np.ndarray:
    return np.column_stack([distance**j for j in range(p + 1)])


def target_2d(p: int, deriv: tuple[int, int], tangvec: np.ndarray | None, row: int) -> np.ndarray:
    target = np.zeros(len(multi_indices_2d(p)))
    if tangvec is not None:
        if p < 1:
            raise ValueError("tangvec requires p >= 1.")
        idx = multi_indices_2d(p)
        target[idx.index((1, 0))] = tangvec[row, 0]
        target[idx.index((0, 1))] = tangvec[row, 1]
        return target
    if deriv in multi_indices_2d(p):
        pos = multi_indices_2d(p).index(deriv)
        target[pos] = float(math.factorial(deriv[0]) * math.factorial(deriv[1]))
    return target


def target_1d(p: int, deriv: int = 0) -> np.ndarray:
    if deriv > p:
        raise ValueError("deriv must be less than or equal to p.")
    target = np.zeros(p + 1)
    target[deriv] = float(math.factorial(deriv))
    return target


def weighted_pinv_design(design: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    keep = weights > 0
    X = design[keep, :]
    w = weights[keep]
    if X.shape[0] == 0:
        raise ValueError("No observations inside the bandwidth.")
    WX = X * w[:, None]
    gram = X.T @ WX
    inv_gram = np.linalg.pinv(gram, rcond=1e-12)
    return keep, X, inv_gram


def cluster_indices(cluster_values: np.ndarray, groups=None) -> tuple[np.ndarray, np.ndarray]:
    cluster_values = np.asarray(cluster_values)
    if groups is None:
        try:
            groups, codes = np.unique(cluster_values, return_inverse=True)
            return np.asarray(groups), np.asarray(codes, dtype=int)
        except TypeError:
            groups = pd.unique(cluster_values)
    else:
        groups = np.asarray(groups)

    try:
        order = np.argsort(groups, kind="mergesort")
        sorted_groups = groups[order]
        pos = np.searchsorted(sorted_groups, cluster_values)
        valid = (pos >= 0) & (pos < len(sorted_groups))
        if np.any(valid):
            valid_idx = np.nonzero(valid)[0]
            valid[valid_idx] = sorted_groups[pos[valid_idx]] == cluster_values[valid_idx]
        codes = np.full(cluster_values.shape[0], -1, dtype=int)
        codes[valid] = order[pos[valid]]
        return groups, codes
    except (TypeError, ValueError):
        codes = pd.Categorical(cluster_values, categories=groups).codes
        return groups, np.asarray(codes, dtype=int)


def cluster_sums(values: np.ndarray, cluster_values: np.ndarray, groups=None, codes: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    if codes is None:
        groups, codes = cluster_indices(cluster_values, groups)
    else:
        groups = np.asarray(groups)
        codes = np.asarray(codes, dtype=int)
    out = np.zeros((len(groups), values.shape[1]))
    if values.shape[0] == 0:
        return out
    valid = codes >= 0
    np.add.at(out, codes[valid], values[valid, :])
    return out


@dataclass
class LocalFit:
    estimate: np.ndarray
    se: np.ndarray
    influence: np.ndarray
    n_eff: int


def local_fit_targets(
    design_full: np.ndarray,
    weights_full: np.ndarray,
    outcomes_full: np.ndarray,
    target: np.ndarray,
    *,
    vce: str = "hc1",
    cluster: np.ndarray | None = None,
    cluster_groups: np.ndarray | None = None,
    scale_override: float | None = None,
) -> LocalFit:
    outcomes_full = np.asarray(outcomes_full, dtype=float)
    if outcomes_full.ndim == 1:
        outcomes_full = outcomes_full.reshape(-1, 1)

    keep, X, inv_gram = weighted_pinv_design(design_full, weights_full)
    w = weights_full[keep]
    Y = outcomes_full[keep, :]
    beta = inv_gram @ (X.T @ (w[:, None] * Y))
    fitted = X @ beta
    resid = Y - fitted
    n_eff = X.shape[0]
    k = X.shape[1]

    leverage = np.sum((X @ inv_gram) * X, axis=1) * w
    adj = np.ones(n_eff)
    vce = vce.lower()
    if vce == "hc2":
        adj = 1.0 / np.maximum(1.0 - leverage, 1e-8)
    elif vce == "hc3":
        adj = 1.0 / np.maximum(1.0 - leverage, 1e-8) ** 2

    # Influence contribution for the requested linear functional.
    row = target @ inv_gram
    score_base = (X * w[:, None]) @ row
    infl_kept = score_base[:, None] * resid * np.sqrt(adj)[:, None]

    if cluster is not None:
        cluster_kept = np.asarray(cluster)[keep]
        if cluster_groups is None:
            groups, codes = cluster_indices(cluster_kept)
            active_groups = len(groups)
        else:
            groups = np.asarray(cluster_groups)
            _, codes = cluster_indices(cluster_kept, groups)
            active_groups = int(np.unique(codes[codes >= 0]).size)
        summed = cluster_sums(infl_kept, cluster_kept, groups, codes)
        scale = 1.0
        if scale_override is not None:
            scale = float(scale_override)
        elif n_eff > k:
            if vce == "hc1":
                scale *= n_eff / (n_eff - k)
            if active_groups > 1:
                scale *= ((n_eff - 1.0) / (n_eff - k)) * (active_groups / (active_groups - 1.0))
        cov = scale * (summed.T @ summed)
        infl_source = np.sqrt(scale) * summed
    else:
        scale = 1.0
        if scale_override is not None:
            scale = float(scale_override)
        elif vce == "hc1" and n_eff > k:
            scale = n_eff / (n_eff - k)
        cov = scale * (infl_kept.T @ infl_kept)
        infl_source = np.sqrt(scale) * infl_kept

    estimate = target @ beta
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))

    infl_full = np.zeros((design_full.shape[0], outcomes_full.shape[1]))
    if cluster is None:
        infl_full[keep, :] = infl_source
    elif cluster_groups is not None:
        infl_full = infl_source
    else:
        # For cross-evaluation covariance, keep cluster-level sums in rows
        # matching the first occurrence of each cluster. This preserves sums
        # without needing a second representation.
        cluster_all = np.asarray(cluster)
        all_groups, _ = cluster_indices(cluster_all)
        infl_full = np.zeros((len(all_groups), outcomes_full.shape[1]))
        _, group_rows = cluster_indices(groups, all_groups)
        valid = group_rows >= 0
        infl_full[group_rows[valid], :] = infl_source[valid, :]

    return LocalFit(estimate=np.asarray(estimate), se=se, influence=infl_full, n_eff=int(n_eff))


def ci_columns(est: np.ndarray, se: np.ndarray, level: float, side: str) -> tuple[np.ndarray, np.ndarray]:
    if side == "two":
        cval = stats.norm.ppf((level + 100.0) / 200.0)
        return est - cval * se, est + cval * se
    cval = stats.norm.ppf(level / 100.0)
    if side == "left":
        return np.repeat(-np.inf, len(est)), est + cval * se
    if side == "right":
        return est - cval * se, np.repeat(np.inf, len(est))
    raise ValueError("side must be two, left, or right.")


def p_values(tvalues: np.ndarray) -> np.ndarray:
    return 2.0 * stats.norm.sf(np.abs(tvalues))


def bandwidth_floor(values: np.ndarray, bwcheck: int | None) -> float:
    values = np.sort(np.asarray(values, dtype=float)[np.isfinite(values)])
    if values.size == 0:
        return np.nan
    if bwcheck is None:
        return 0.0
    k = min(max(int(bwcheck), 1), values.size)
    return float(values[k - 1])


def cer_factor(n: int, p: int) -> float:
    n = max(int(n), 1)
    return n ** (1.0 / (2.0 * p + 4.0) - 1.0 / (p + 4.0))


def ensure_dataframe(x: np.ndarray, columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(x, columns=columns)


def covariance_from_influence(influence: np.ndarray) -> np.ndarray:
    influence = np.asarray(influence, dtype=float)
    return influence @ influence.T
