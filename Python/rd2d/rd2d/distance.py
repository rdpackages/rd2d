from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ._utils import (
    as_1d,
    as_2d,
    bandwidth_floor,
    basis_1d,
    cer_factor,
    check_lengths,
    ci_columns,
    complete_cases,
    covariance_from_influence,
    kernel_weights,
    local_fit_targets,
    p_values,
    target_1d,
    validate_order,
)
from .location import _bwselect_base, _is_cer, _is_common, _se_from_influence
from .results import RD2DResult


DISTANCE_MAIN_COLUMNS = [
    "b1",
    "b2",
    "estimate.p",
    "std.err.p",
    "estimate.q",
    "std.err.q",
    "t.value",
    "p.value",
    "ci.lower",
    "ci.upper",
    "h0",
    "h1",
    "h0.rbc",
    "h1.rbc",
    "N.Co",
    "N.Tr",
]


def _clean_distance_inputs(y, distance, b=None, fuzzy=None, cluster=None):
    y = as_1d(y, "Y")
    distance = as_2d(distance, "distance")
    check_lengths(y, distance)
    if b is None:
        b_arr = np.full((distance.shape[1], 2), np.nan)
    else:
        b_arr = as_2d(b, "b", ncol=2)
        if b_arr.shape[0] != distance.shape[1]:
            raise ValueError("b must have one row for each distance column.")
    arrays = [y, distance]
    if fuzzy is not None:
        fuzzy = as_1d(fuzzy, "fuzzy")
        check_lengths(y, fuzzy)
        arrays.append(fuzzy)
    if cluster is not None:
        cluster = np.asarray(cluster)
        check_lengths(y, cluster)
        arrays.append(cluster.astype(str))
    keep = complete_cases(*arrays)
    y = y[keep]
    distance = distance[keep, :]
    if fuzzy is not None:
        fuzzy = fuzzy[keep]
    if cluster is not None:
        cluster = np.asarray(cluster)[keep]
    return y, distance, b_arr, fuzzy, cluster


def _validate_kink_unknown(kink_unknown) -> tuple[bool, bool]:
    arr = np.asarray(kink_unknown, dtype=bool)
    if arr.ndim == 0 or arr.size == 1:
        arr = np.repeat(bool(arr.ravel()[0]), 2)
    if arr.size != 2:
        raise ValueError("kink_unknown must be a scalar or length-2 logical vector.")
    if bool(arr[1]) and not bool(arr[0]):
        raise ValueError("kink_unknown[1] can be True only when kink_unknown[0] is True.")
    return bool(arr[0]), bool(arr[1])


def _validate_kink_position(kink_position, neval: int, b) -> np.ndarray:
    out = np.zeros(neval, dtype=bool)
    if kink_position is None:
        return out
    arr = np.asarray(kink_position)
    if arr.dtype == bool:
        if arr.size != neval:
            raise ValueError("kink_position must have one value for each boundary point.")
        return arr.astype(bool)
    idx = np.asarray(arr, dtype=int)
    if np.any(idx < 0) or np.any(idx >= neval):
        raise ValueError("Python kink_position uses zero-based indices between 0 and neval - 1.")
    out[idx] = True
    if np.any(out) and b is None:
        raise ValueError("kink_position requires b.")
    return out


def _distance_to_known_kink(b: np.ndarray, kink_position: np.ndarray) -> np.ndarray:
    if not np.any(kink_position) or np.any(~np.isfinite(b)):
        return np.repeat(np.inf, b.shape[0])
    kink = b[kink_position, :]
    dist = np.sqrt(((b[:, None, :] - kink[None, :, :]) ** 2).sum(axis=2))
    return np.min(dist, axis=1)


def _qr_crossprod_inv(x: np.ndarray) -> np.ndarray:
    gram = x.T @ x
    try:
        chol = np.linalg.cholesky(gram)
        inv_chol = np.linalg.inv(chol)
        return inv_chol.T @ inv_chol
    except np.linalg.LinAlgError:
        return np.linalg.pinv(gram, rcond=1e-20)


def _kernel_label(kernel: str) -> str:
    kernel = kernel.lower()
    if kernel in {"tri", "triangular"}:
        return "Triangular"
    if kernel in {"uni", "uniform"}:
        return "Uniform"
    if kernel in {"gau", "gaussian"}:
        return "Gaussian"
    return "Epanechnikov"


def _distance_rot(distance_abs: np.ndarray, kernel: str) -> float:
    constants = {
        "Epanechnikov": (1.0 / 6.0, 4.0 / (3.0 * np.pi)),
        "Triangular": (3.0 / 20.0, 3.0 / (2.0 * np.pi)),
        "Uniform": (1.0 / 4.0, 1.0 / np.pi),
        "Gaussian": (1.0, 1.0 / (4.0 * np.pi)),
    }
    mu2K_squared, l2K_squared = constants[_kernel_label(kernel)]
    var_hat = 0.5 * float(np.mean(distance_abs ** 2))
    trace_const = 1.0 / (2.0 * np.pi * var_hat ** 3)
    return float(((2.0 * l2K_squared) / (len(distance_abs) * mu2K_squared * trace_const)) ** (1.0 / 6.0))


def _poly_x(distance: np.ndarray, degree: int) -> np.ndarray:
    distance = np.asarray(distance, dtype=float)
    return np.column_stack([distance ** j for j in range(degree + 1)])


def _poly_lm(y: np.ndarray, distance: np.ndarray, degree: int) -> dict[str, Any]:
    X = _poly_x(distance, degree)
    invXX = _qr_crossprod_inv(X)
    coef = invXX @ (X.T @ y)
    res = y - X @ coef
    sigma2 = float(np.sum(res ** 2) / max(X.shape[0] - X.shape[1], 1))
    se = np.sqrt(np.maximum(np.diag(invXX) * sigma2, 0.0))
    return {"coef": coef, "se": se, "predict": lambda newdistance: _poly_x(np.asarray(newdistance, dtype=float), degree) @ coef}


def _basis_1d_scaled(distance_abs: np.ndarray, h: float, p: int) -> np.ndarray:
    u = np.asarray(distance_abs, dtype=float) / h
    return np.column_stack([u ** j for j in range(p + 1)])


def _vce_const_1d(w_R: np.ndarray, resd: np.ndarray, h: float) -> np.ndarray:
    scores = np.asarray(w_R, dtype=float) * np.asarray(resd, dtype=float)[:, None]
    return scores.T @ scores * h ** 2


def _local_intercepts_multi(y: np.ndarray, fuzzy: np.ndarray, distance_abs: np.ndarray, h: float, p: int, kernel: str) -> tuple[float, float]:
    w = kernel_weights(distance_abs / h, kernel) / h ** 2
    keep = w > 0
    if np.sum(keep) <= p + 1:
        return np.nan, np.nan
    R = _basis_1d_scaled(distance_abs[keep], h, p)
    sqrtw_R = np.sqrt(w[keep])[:, None] * R
    invG = _qr_crossprod_inv(sqrtw_R)
    beta = invG @ (R.T @ (w[keep, None] * np.column_stack([y[keep], fuzzy[keep]])))
    return float(beta[0, 0]), float(beta[0, 1])


def _distance_bw_constants(
    y: np.ndarray,
    dist: np.ndarray,
    p: int,
    kernel: str,
    vce: str,
    bwcheck: int | None,
    scaleregul: float,
    cqt: float,
    fuzzy: np.ndarray | None,
    bwparam: str,
) -> pd.DataFrame:
    rows = []
    for j in range(dist.shape[1]):
        signed = dist[:, j]
        side1 = signed >= 0
        u = np.abs(signed)
        u0 = u[~side1]
        u1 = u[side1]
        y0 = y[~side1].copy()
        y1 = y[side1].copy()
        dn = _distance_rot(u, kernel)
        dn0 = dn1 = dn
        bw_min0 = bw_min1 = 0.0
        bw_max0 = float(np.max(u0))
        bw_max1 = float(np.max(u1))
        if bwcheck is not None:
            su0 = np.sort(u0)
            su1 = np.sort(u1)
            bw_min0 = float(su0[min(bwcheck, len(su0)) - 1])
            bw_min1 = float(su1[min(bwcheck, len(su1)) - 1])
            bw_max0 = float(su0[-1])
            bw_max1 = float(su1[-1])
            dn0 = min(max(dn, bw_min0), bw_max0)
            dn1 = min(max(dn, bw_min1), bw_max1)

        degree = p + 1
        thr0 = float(np.quantile(u0, cqt))
        thr1 = float(np.quantile(u1, cqt))
        fit0 = _poly_lm(y0[u0 <= thr0], u0[u0 <= thr0], degree)
        fit1 = _poly_lm(y1[u1 <= thr1], u1[u1 <= thr1], degree)

        if fuzzy is not None and bwparam == "main":
            f0 = fuzzy[~side1]
            f1 = fuzzy[side1]
            mu_y0, mu_f0 = _local_intercepts_multi(y0, f0, u0, dn0, p, kernel)
            mu_y1, mu_f1 = _local_intercepts_multi(y1, f1, u1, dn1, p, kernel)
            tau_itt = mu_y1 - mu_y0
            tau_fs = mu_f1 - mu_f0
            if np.isfinite(tau_fs) and abs(tau_fs) > np.sqrt(np.finfo(float).eps):
                grad_itt = 1.0 / tau_fs
                grad_fs = -tau_itt / tau_fs ** 2
                y0 = grad_itt * y0 + grad_fs * f0
                y1 = grad_itt * y1 + grad_fs * f1
                fit0 = _poly_lm(y0[u0 <= thr0], u0[u0 <= thr0], degree)
                fit1 = _poly_lm(y1[u1 <= thr1], u1[u1 <= thr1], degree)

        vec = np.zeros(p + 1)
        vec[0] = 1.0

        def side_constants(us: np.ndarray, ys: np.ndarray, h: float, model: dict[str, Any]):
            w = kernel_weights(us / h, kernel) / h ** 2
            keep = w > 0
            eX = us[keep]
            eY = ys[keep]
            ew = w[keep]
            R = _basis_1d_scaled(eX, h, p)
            sqrtw_R = np.sqrt(ew)[:, None] * R
            invG = _qr_crossprod_inv(sqrtw_R)
            sigma_half = ew[:, None] * R
            resd = eY - model["predict"](eX)
            eN = int(np.sum(keep))
            k = p + 1
            if vce == "hc1" and eN > k:
                resd = resd * np.sqrt(eN / (eN - k))
            elif vce in {"hc2", "hc3"}:
                hii = np.sum((sqrtw_R @ invG) * sqrtw_R, axis=1)
                if vce == "hc2":
                    resd = resd * np.sqrt(1.0 / np.maximum(1.0 - hii, 1e-12))
                else:
                    resd = resd * (1.0 / np.maximum(1.0 - hii, 1e-12))
            sigma = _vce_const_1d(sigma_half, resd, h)
            pmatrix = R.T @ (((eX / h) ** (p + 1) * ew).reshape(-1, 1))
            coeff = np.zeros(p + 2)
            coeff[-1] = float(np.asarray(vec.reshape(1, -1) @ invG @ pmatrix).item())
            B = float(coeff @ model["coef"])
            Reg = float(coeff[-1] ** 2 * model["se"][-1] ** 2)
            V = float(np.asarray(vec.reshape(1, -1) @ invG @ sigma @ invG @ vec.reshape(-1, 1)).item())
            return B, V, Reg, eN

        B0, V0, R0, eN0 = side_constants(u0, y0, dn0, fit0)
        B1, V1, R1, eN1 = side_constants(u1, y1, dn1, fit1)
        hn = (2 * (V0 + V1) / ((2 * p + 2) * ((B0 - B1) ** 2 + scaleregul * R0 + scaleregul * R1))) ** (1.0 / (2 * p + 4))
        rows.append([hn, hn, B0, B1, V0, V1, R0, R1, eN0, eN1, bw_min0, bw_min1, bw_max0, bw_max1])
    return pd.DataFrame(rows, columns=["h.0", "h.1", "b.0", "b.1", "v.0", "v.1", "r.0", "r.1", "N.Co", "N.Tr", "bw.min.0", "bw.min.1", "bw.max.0", "bw.max.1"])


def _distance_final_bwcheck(h0: np.ndarray, h1: np.ndarray, dist: np.ndarray, bwcheck: int | None) -> tuple[np.ndarray, np.ndarray]:
    if bwcheck is None:
        return h0, h1
    h0 = h0.copy()
    h1 = h1.copy()
    for j in range(dist.shape[1]):
        signed = dist[:, j]
        side1 = signed >= 0
        u = np.abs(signed)
        u0 = np.sort(u[~side1])
        u1 = np.sort(u[side1])
        h0[j] = min(max(h0[j], u0[min(bwcheck, len(u0)) - 1]), u0[-1])
        h1[j] = min(max(h1[j], u1[min(bwcheck, len(u1)) - 1]), u1[-1])
    return h0, h1


def rdbw2d_distance(
    Y,
    distance,
    *,
    b=None,
    p: int = 1,
    kink_unknown=(False, False),
    kink_position=None,
    kernel: str = "tri",
    bwselect: str = "mserd",
    bwparam: str = "main",
    vce: str = "hc1",
    bwcheck: int | None = None,
    masspoints: str = "check",
    cluster=None,
    scaleregul: float = 1.0,
    cqt: float = 0.5,
    fuzzy=None,
) -> RD2DResult:
    """Bandwidth selection for distance-based boundary discontinuity designs.

    Parameters
    ----------
    Y : array-like
        Outcome variable.
    distance : array-like, shape (n_obs, n_eval)
        Signed distance scores. Nonnegative values identify observations on the
        treated side and negative values identify observations on the control
        side.
    b : array-like, optional
        Boundary evaluation points with two columns.
    p : int, default 1
        Local polynomial order.
    kink_unknown : bool or pair of bool, default ``(False, False)``
        Controls unknown-kink bandwidth adjustments for point estimation and
        inference.
    kink_position : array-like, optional
        Known kink locations, supplied as zero-based indices or a logical vector
        with one entry per evaluation point.
    fuzzy : array-like, optional
        Treatment receipt/status variable for fuzzy designs.

    Returns
    -------
    RD2DResult
        Result object with a ``bws`` table containing ``b1``, ``b2``, ``h0``,
        ``h1``, ``N.Co``, and ``N.Tr``.
    """

    y, dist, b_arr, fuzzy, cluster = _clean_distance_inputs(Y, distance, b, fuzzy, cluster)
    p = validate_order(p, "p")
    kink_unknown = _validate_kink_unknown(kink_unknown)
    kink_position_arr = _validate_kink_position(kink_position, dist.shape[1], b)
    neval = dist.shape[1]
    bwselect = bwselect.lower()
    base_select = _bwselect_base(bwselect)
    if base_select not in {"mserd", "msetwo", "imserd", "imsetwo"}:
        raise ValueError("unsupported bwselect.")
    mseconsts = _distance_bw_constants(y, dist, p, kernel, vce, bwcheck, scaleregul, cqt, fuzzy, bwparam if fuzzy is not None else "main")
    if base_select == "mserd":
        hn = (2 * (mseconsts["v.0"] + mseconsts["v.1"]) / ((2 * p + 2) * ((mseconsts["b.0"] - mseconsts["b.1"]) ** 2 + scaleregul * mseconsts["r.0"] + scaleregul * mseconsts["r.1"]))) ** (1.0 / (2 * p + 4))
        if bwcheck is not None:
            hn = np.maximum(hn, np.maximum(mseconsts["bw.min.0"], mseconsts["bw.min.1"]))
            hn = np.minimum(hn, np.maximum(mseconsts["bw.max.0"], mseconsts["bw.max.1"]))
        mseconsts["h.0"] = hn
        mseconsts["h.1"] = hn
    elif base_select == "msetwo":
        h0 = (2 * mseconsts["v.0"] / ((2 * p + 2) * (mseconsts["b.0"] ** 2 + scaleregul * mseconsts["r.0"]))) ** (1.0 / (2 * p + 4))
        h1 = (2 * mseconsts["v.1"] / ((2 * p + 2) * (mseconsts["b.1"] ** 2 + scaleregul * mseconsts["r.1"]))) ** (1.0 / (2 * p + 4))
        if bwcheck is not None:
            h0 = np.minimum(np.maximum(h0, mseconsts["bw.min.0"]), mseconsts["bw.max.0"])
            h1 = np.minimum(np.maximum(h1, mseconsts["bw.min.1"]), mseconsts["bw.max.1"])
        mseconsts["h.0"] = h0
        mseconsts["h.1"] = h1
    elif base_select == "imserd":
        V = float(np.mean(mseconsts["v.0"]) + np.mean(mseconsts["v.1"]))
        B = float(np.mean((mseconsts["b.0"] - mseconsts["b.1"]) ** 2 + scaleregul * mseconsts["r.0"] + scaleregul * mseconsts["r.1"]))
        h = (2 * V / ((2 * p + 2) * B)) ** (1.0 / (2 * p + 4))
        hn = np.repeat(h, neval)
        if bwcheck is not None:
            hn = np.maximum(hn, np.maximum(mseconsts["bw.min.0"], mseconsts["bw.min.1"]))
            hn = np.minimum(hn, np.maximum(mseconsts["bw.max.0"], mseconsts["bw.max.1"]))
        mseconsts["h.0"] = hn
        mseconsts["h.1"] = hn
    else:
        V0 = float(np.mean(mseconsts["v.0"]))
        V1 = float(np.mean(mseconsts["v.1"]))
        B0 = float(np.mean(mseconsts["b.0"] ** 2 + scaleregul * mseconsts["r.0"]))
        B1 = float(np.mean(mseconsts["b.1"] ** 2 + scaleregul * mseconsts["r.1"]))
        h0 = np.repeat((2 * V0 / ((2 * p + 2) * B0)) ** (1.0 / (2 * p + 4)), neval)
        h1 = np.repeat((2 * V1 / ((2 * p + 2) * B1)) ** (1.0 / (2 * p + 4)), neval)
        if bwcheck is not None:
            h0 = np.minimum(np.maximum(h0, mseconsts["bw.min.0"]), mseconsts["bw.max.0"])
            h1 = np.minimum(np.maximum(h1, mseconsts["bw.min.1"]), mseconsts["bw.max.1"])
        mseconsts["h.0"] = h0
        mseconsts["h.1"] = h1

    smooth_exp = 1.0 / (p + 4.0) if _is_cer(bwselect) else 1.0 / (2.0 * p + 4.0)
    M = np.empty(neval, dtype=int)
    M0 = np.empty(neval, dtype=int)
    M1 = np.empty(neval, dtype=int)
    for j in range(neval):
        unique = np.unique(dist[:, j])
        M0[j] = int(np.sum(unique < 0))
        M1[j] = int(np.sum(unique >= 0))
        M[j] = M0[j] + M1[j]
    if _is_cer(bwselect):
        if _is_common(bwselect):
            factor = np.array([cer_factor(int(m), p) for m in M])
            mseconsts["h.0"] *= factor
            mseconsts["h.1"] *= factor
        else:
            mseconsts["h.0"] *= np.array([cer_factor(int(m), p) for m in M0])
            mseconsts["h.1"] *= np.array([cer_factor(int(m), p) for m in M1])

    if kink_unknown[0]:
        if _is_common(bwselect):
            factor = M.astype(float) ** (smooth_exp - 1.0 / 4.0)
            mseconsts["h.0"] *= factor
            mseconsts["h.1"] *= factor
        else:
            mseconsts["h.0"] *= M0.astype(float) ** (smooth_exp - 1.0 / 4.0)
            mseconsts["h.1"] *= M1.astype(float) ** (smooth_exp - 1.0 / 4.0)

    bws = pd.DataFrame(
        {
            "b1": b_arr[:, 0],
            "b2": b_arr[:, 1],
            "h0": mseconsts["h.0"].to_numpy(dtype=float),
            "h1": mseconsts["h.1"].to_numpy(dtype=float),
            "N.Co": mseconsts["N.Co"].to_numpy(dtype=float),
            "N.Tr": mseconsts["N.Tr"].to_numpy(dtype=float),
        }
    )
    known_dist = _distance_to_known_kink(b_arr, kink_position_arr)
    finite = np.isfinite(known_dist)
    if np.any(finite):
        smooth_exp = 1.0 / (p + 4.0) if _is_cer(bwselect) else 1.0 / (2.0 * p + 4.0)
        if _is_common(bwselect):
            factor = M[finite].astype(float) ** (smooth_exp - 1.0 / 4.0)
            h_rate = bws.loc[finite, "h0"].to_numpy(dtype=float) * factor
            h_bound = np.maximum(h_rate, known_dist[finite])
            bws.loc[finite, "h0"] = np.minimum(bws.loc[finite, "h0"], h_bound)
            bws.loc[finite, "h1"] = np.minimum(bws.loc[finite, "h1"], h_bound)
        else:
            h_rate0 = bws.loc[finite, "h0"].to_numpy(dtype=float) * M0[finite].astype(float) ** (smooth_exp - 1.0 / 4.0)
            h_rate1 = bws.loc[finite, "h1"].to_numpy(dtype=float) * M1[finite].astype(float) ** (smooth_exp - 1.0 / 4.0)
            bws.loc[finite, "h0"] = np.minimum(bws.loc[finite, "h0"], np.maximum(h_rate0, known_dist[finite]))
            bws.loc[finite, "h1"] = np.minimum(bws.loc[finite, "h1"], np.maximum(h_rate1, known_dist[finite]))
    return RD2DResult(
        bws=bws,
        mseconsts=mseconsts.copy(),
        opt={
            "N": len(y),
            "N.0": int(np.sum(dist[:, 0] < 0)),
            "N.1": int(np.sum(dist[:, 0] >= 0)),
            "neval": neval,
            "p": p,
            "kernel": kernel,
            "bwselect": bwselect,
            "bwparam": "main" if fuzzy is None else bwparam,
            "vce": vce,
            "bwcheck": bwcheck,
            "masspoints": masspoints,
            "cluster": cluster,
            "clustered": cluster is not None,
            "scaleregul": scaleregul,
            "cqt": cqt,
            "fuzzy": fuzzy is not None,
            "kink_unknown": kink_unknown,
            "kink_position": kink_position_arr,
        },
        rdmodel="rdbw2d_distance",
    )


def rdbw2d_dist(*args, **kwargs) -> RD2DResult:
    return rdbw2d_distance(*args, **kwargs)


def _parse_distance_h(h, neval: int) -> tuple[np.ndarray, np.ndarray]:
    h = np.asarray(h, dtype=float)
    if h.ndim == 0 or h.size == 1:
        value = float(np.ravel(h)[0])
        if value <= 0:
            raise ValueError("scalar h must be positive.")
        return np.repeat(value, neval), np.repeat(value, neval)
    if h.shape != (neval, 2):
        raise ValueError("h must be a positive scalar or a J x 2 matrix.")
    if np.any(h <= 0):
        raise ValueError("all bandwidths in h must be positive.")
    return h[:, 0], h[:, 1]


def _fit_distance_order(
    outcomes: np.ndarray,
    dist: np.ndarray,
    h0: np.ndarray,
    h1: np.ndarray,
    p: int,
    kernel: str,
    vce: str,
    cluster,
) -> dict[str, Any]:
    nout = outcomes.shape[1]
    neval = dist.shape[1]
    mu0 = np.full((neval, nout), np.nan)
    mu1 = np.full((neval, nout), np.nan)
    se0 = np.full((neval, nout), np.nan)
    se1 = np.full((neval, nout), np.nan)
    n0 = np.zeros(neval, dtype=int)
    n1 = np.zeros(neval, dtype=int)
    infl0 = [[None for _ in range(nout)] for _ in range(neval)]
    infl1 = [[None for _ in range(nout)] for _ in range(neval)]
    target = target_1d(p, 0)
    for j in range(neval):
        signed = dist[:, j]
        side1 = signed >= 0
        u = np.abs(signed)
        design = basis_1d(u, p)
        w0 = kernel_weights(u / h0[j], kernel) / max(h0[j] ** 2, np.finfo(float).eps)
        w1 = kernel_weights(u / h1[j], kernel) / max(h1[j] ** 2, np.finfo(float).eps)
        w0 = np.where(~side1, w0, 0.0)
        w1 = np.where(side1, w1, 0.0)
        try:
            fit0 = local_fit_targets(design, w0, outcomes, target, vce=vce, cluster=cluster)
            mu0[j, :] = fit0.estimate
            se0[j, :] = fit0.se
            n0[j] = fit0.n_eff
            for k in range(nout):
                infl0[j][k] = fit0.influence[:, k]
        except ValueError:
            pass
        try:
            fit1 = local_fit_targets(design, w1, outcomes, target, vce=vce, cluster=cluster)
            mu1[j, :] = fit1.estimate
            se1[j, :] = fit1.se
            n1[j] = fit1.n_eff
            for k in range(nout):
                infl1[j][k] = fit1.influence[:, k]
        except ValueError:
            pass
    return {"mu0": mu0, "mu1": mu1, "se0": se0, "se1": se1, "N0": n0, "N1": n1, "infl0": infl0, "infl1": infl1}


def _stack_influence(infl: list[list[np.ndarray | None]], outcome: int) -> np.ndarray | None:
    rows = []
    for row in infl:
        value = row[outcome]
        if value is None:
            return None
        rows.append(np.asarray(value, dtype=float))
    return np.vstack(rows)


def _contrast_influence(fit: dict[str, Any], outcome: int) -> np.ndarray | None:
    a = _stack_influence(fit["infl1"], outcome)
    b = _stack_influence(fit["infl0"], outcome)
    if a is None or b is None:
        return None
    return a - b


def _side_influence(fit: dict[str, Any], outcome: int, side: int) -> np.ndarray | None:
    return _stack_influence(fit["infl1" if side == 1 else "infl0"], outcome)


def _make_distance_table(
    b: np.ndarray,
    est_p: np.ndarray,
    se_p: np.ndarray,
    est_q: np.ndarray,
    se_q: np.ndarray,
    h0: np.ndarray,
    h1: np.ndarray,
    h0_rbc: np.ndarray,
    h1_rbc: np.ndarray,
    n0: np.ndarray,
    n1: np.ndarray,
    level: float,
    side: str,
) -> pd.DataFrame:
    tval = est_q / se_q
    ci_l, ci_u = ci_columns(est_q, se_q, level, side)
    data = np.column_stack(
        [
            b[:, 0],
            b[:, 1],
            est_p,
            se_p,
            est_q,
            se_q,
            tval,
            p_values(tval),
            ci_l,
            ci_u,
            h0,
            h1,
            h0_rbc,
            h1_rbc,
            n0,
            n1,
        ]
    )
    return pd.DataFrame(data, columns=DISTANCE_MAIN_COLUMNS)


def rd2d_distance(
    Y,
    distance,
    *,
    h=None,
    b=None,
    p: int = 1,
    q: int | None = None,
    kink_unknown=(False, False),
    kink_position=None,
    kernel: str = "tri",
    level: float = 95,
    cbands: bool = True,
    side: str = "two",
    bwselect: str = "mserd",
    bwparam: str = "main",
    params_other: list[str] | tuple[str, ...] | str | None = None,
    params_cov: list[str] | tuple[str, ...] | str | None = None,
    vce: str = "hc1",
    bwcheck: int | None = None,
    masspoints: str = "check",
    cluster=None,
    scaleregul: float = 1.0,
    cqt: float = 0.5,
    fuzzy=None,
) -> RD2DResult:
    """Distance-based local polynomial boundary discontinuity estimator.

    Parameters
    ----------
    Y : array-like
        Outcome variable.
    distance : array-like, shape (n_obs, n_eval)
        Signed distance scores. Each column corresponds to one boundary
        evaluation point.
    h : float or array-like, optional
        User-supplied bandwidths. A scalar uses the same bandwidth for all
        sides and evaluation points. A matrix must have columns ``h0`` and
        ``h1``.
    b : array-like, optional
        Boundary evaluation points with two columns.
    p, q : int
        Polynomial orders for point estimation and robust bias-corrected
        inference.
    cbands : bool, default True
        If true, stores the main covariance matrix needed for uniform
        confidence bands and aggregate inference.
    fuzzy : array-like, optional
        Treatment receipt/status variable for fuzzy designs.

    Returns
    -------
    RD2DResult
        Result object with ``main`` and ``bw`` tables. Sharp fits include
        side-specific ``main.0`` and ``main.1`` tables. Fuzzy fits include
        ``itt`` and ``fs`` tables.
    """

    y, dist, b_arr, fuzzy, cluster = _clean_distance_inputs(Y, distance, b, fuzzy, cluster)
    p = validate_order(p, "p")
    kink_unknown = _validate_kink_unknown(kink_unknown)
    if q is None:
        q = p if kink_unknown[0] else p + 1
    q = validate_order(q, "q")
    if q < p:
        raise ValueError("q must be greater than or equal to p.")
    side = side.lower()
    if side not in {"two", "left", "right"}:
        raise ValueError("side must be two, left, or right.")
    if params_other is None:
        params_other_set: set[str] = set()
    elif isinstance(params_other, str):
        params_other_set = {params_other}
    else:
        params_other_set = set(params_other)
    if params_cov is None:
        params_cov_set: set[str] = set()
    elif isinstance(params_cov, str):
        params_cov_set = {params_cov}
    else:
        params_cov_set = set(params_cov)
    if cbands:
        params_cov_set.add("main")
    kink_position_arr = _validate_kink_position(kink_position, dist.shape[1], b)
    neval = dist.shape[1]

    fit_bwcheck = 50 + p + 1 if bwcheck is None else bwcheck
    if h is None:
        bws = rdbw2d_distance(
            y,
            dist,
            b=b_arr,
            p=p,
            kink_unknown=kink_unknown,
            kink_position=kink_position_arr,
            kernel=kernel,
            bwselect=bwselect,
            bwparam=bwparam,
            vce=vce,
            bwcheck=fit_bwcheck,
            masspoints=masspoints,
            cluster=cluster,
            scaleregul=scaleregul,
            cqt=cqt,
            fuzzy=fuzzy,
        ).bws
        h0 = bws["h0"].to_numpy(dtype=float)
        h1 = bws["h1"].to_numpy(dtype=float)
    else:
        bwselect = "user provided"
        h0, h1 = _parse_distance_h(h, neval)

    h0, h1 = _distance_final_bwcheck(h0, h1, dist, fit_bwcheck)
    h0_rbc = h0.copy()
    h1_rbc = h1.copy()
    if kink_unknown[1]:
        smooth_exp = 1.0 / (p + 4.0) if _is_cer(bwselect) else 1.0 / (2.0 * p + 4.0)
        for j in range(neval):
            signed = dist[:, j]
            n0 = max(int(np.sum(signed < 0)), 1)
            n1 = max(int(np.sum(signed >= 0)), 1)
            if _is_common(bwselect):
                factor = max(len(y), 1) ** (smooth_exp - 1.0 / 4.0)
                h0_rbc[j] *= factor
                h1_rbc[j] *= factor
            else:
                h0_rbc[j] *= n0 ** (smooth_exp - 1.0 / 4.0)
                h1_rbc[j] *= n1 ** (smooth_exp - 1.0 / 4.0)
        h0_rbc, h1_rbc = _distance_final_bwcheck(h0_rbc, h1_rbc, dist, fit_bwcheck)

    outcomes = y.reshape(-1, 1) if fuzzy is None else np.column_stack([y, fuzzy])
    fit_p = _fit_distance_order(outcomes, dist, h0, h1, p, kernel, vce, cluster)
    fit_q = fit_p if q == p and np.array_equal(h0, h0_rbc) and np.array_equal(h1, h1_rbc) else _fit_distance_order(outcomes, dist, h0_rbc, h1_rbc, q, kernel, vce, cluster)

    infl_tables: dict[str, np.ndarray] = {}
    is_fuzzy = fuzzy is not None
    if not is_fuzzy:
        tau_p = fit_p["mu1"][:, 0] - fit_p["mu0"][:, 0]
        tau_q = fit_q["mu1"][:, 0] - fit_q["mu0"][:, 0]
        infl_main_p = _contrast_influence(fit_p, 0)
        infl_main_q = _contrast_influence(fit_q, 0)
        se_p = _se_from_influence(infl_main_p)
        se_q = _se_from_influence(infl_main_q)
        main = _make_distance_table(b_arr, tau_p, se_p, tau_q, se_q, h0, h1, h0_rbc, h1_rbc, fit_p["N0"], fit_p["N1"], level, side)
        main0 = _make_distance_table(
            b_arr,
            fit_p["mu0"][:, 0],
            fit_p["se0"][:, 0],
            fit_q["mu0"][:, 0],
            fit_q["se0"][:, 0],
            h0,
            np.repeat(np.nan, neval),
            h0_rbc,
            np.repeat(np.nan, neval),
            fit_p["N0"],
            np.repeat(np.nan, neval),
            level,
            side,
        )
        main1 = _make_distance_table(
            b_arr,
            fit_p["mu1"][:, 0],
            fit_p["se1"][:, 0],
            fit_q["mu1"][:, 0],
            fit_q["se1"][:, 0],
            np.repeat(np.nan, neval),
            h1,
            np.repeat(np.nan, neval),
            h1_rbc,
            np.repeat(np.nan, neval),
            fit_p["N1"],
            level,
            side,
        )
        result_tables = {"main": main, "bw": main[["b1", "b2", "h0", "h1", "N.Co", "N.Tr"]].copy(), "main.0": main0, "main.1": main1}
        if infl_main_q is not None:
            infl_tables["main"] = infl_main_q
        for name, side_id in {"main.0": 0, "main.1": 1}.items():
            infl = _side_influence(fit_q, 0, side_id)
            if infl is not None:
                infl_tables[name] = infl
        tau_itt = tau_itt_q = tau_fs = tau_fs_q = None
    else:
        itt_p = fit_p["mu1"][:, 0] - fit_p["mu0"][:, 0]
        fs_p = fit_p["mu1"][:, 1] - fit_p["mu0"][:, 1]
        itt_q = fit_q["mu1"][:, 0] - fit_q["mu0"][:, 0]
        fs_q = fit_q["mu1"][:, 1] - fit_q["mu0"][:, 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            tau_p = itt_p / fs_p
            tau_q = itt_q / fs_q
        infl_itt_p = _contrast_influence(fit_p, 0)
        infl_fs_p = _contrast_influence(fit_p, 1)
        infl_itt_q = _contrast_influence(fit_q, 0)
        infl_fs_q = _contrast_influence(fit_q, 1)
        infl_main_p = infl_itt_p / fs_p[:, None] - (itt_p / (fs_p**2))[:, None] * infl_fs_p
        infl_main_q = infl_itt_q / fs_q[:, None] - (itt_q / (fs_q**2))[:, None] * infl_fs_q
        se_p = _se_from_influence(infl_main_p)
        se_q = _se_from_influence(infl_main_q)
        se_itt_p = _se_from_influence(infl_itt_p)
        se_itt_q = _se_from_influence(infl_itt_q)
        se_fs_p = _se_from_influence(infl_fs_p)
        se_fs_q = _se_from_influence(infl_fs_q)
        main = _make_distance_table(b_arr, tau_p, se_p, tau_q, se_q, h0, h1, h0_rbc, h1_rbc, fit_p["N0"], fit_p["N1"], level, side)
        itt = _make_distance_table(b_arr, itt_p, se_itt_p, itt_q, se_itt_q, h0, h1, h0_rbc, h1_rbc, fit_p["N0"], fit_p["N1"], level, side)
        fs = _make_distance_table(b_arr, fs_p, se_fs_p, fs_q, se_fs_q, h0, h1, h0_rbc, h1_rbc, fit_p["N0"], fit_p["N1"], level, side)
        result_tables = {"main": main, "bw": main[["b1", "b2", "h0", "h1", "N.Co", "N.Tr"]].copy(), "itt": itt, "itt.0": None, "itt.1": None, "fs": fs, "fs.0": None, "fs.1": None}
        infl_tables.update({"main": infl_main_q, "itt": infl_itt_q, "fs": infl_fs_q})
        side_map = {"itt.0": (0, 0), "itt.1": (0, 1), "fs.0": (1, 0), "fs.1": (1, 1)}
        for name, (outcome, side_id) in side_map.items():
            if name in params_other_set:
                mu_p = fit_p["mu1" if side_id == 1 else "mu0"][:, outcome]
                mu_q = fit_q["mu1" if side_id == 1 else "mu0"][:, outcome]
                se_side_p = fit_p["se1" if side_id == 1 else "se0"][:, outcome]
                se_side_q = fit_q["se1" if side_id == 1 else "se0"][:, outcome]
                h0_side = h0 if side_id == 0 else np.repeat(np.nan, neval)
                h1_side = h1 if side_id == 1 else np.repeat(np.nan, neval)
                h0r_side = h0_rbc if side_id == 0 else np.repeat(np.nan, neval)
                h1r_side = h1_rbc if side_id == 1 else np.repeat(np.nan, neval)
                n0_side = fit_p["N0"] if side_id == 0 else np.repeat(np.nan, neval)
                n1_side = fit_p["N1"] if side_id == 1 else np.repeat(np.nan, neval)
                result_tables[name] = _make_distance_table(b_arr, mu_p, se_side_p, mu_q, se_side_q, h0_side, h1_side, h0r_side, h1r_side, n0_side, n1_side, level, side)
            infl = _side_influence(fit_q, outcome, side_id)
            if infl is not None:
                infl_tables[name] = infl
        tau_itt, tau_itt_q, tau_fs, tau_fs_q = itt_p, itt_q, fs_p, fs_q

    cov_tables = {
        name: covariance_from_influence(infl)
        for name, infl in infl_tables.items()
        if name in params_cov_set and infl is not None
    }

    return RD2DResult(
        **result_tables,
        opt={
            "b": b_arr,
            "p": p,
            "q": q,
            "kernel": kernel,
            "kink_unknown": kink_unknown,
            "kink_position": kink_position_arr,
            "N": len(y),
            "N.0": int(np.sum(dist[:, 0] < 0)),
            "N.1": int(np.sum(dist[:, 0] >= 0)),
            "M": len(y),
            "M.0": int(np.sum(dist[:, 0] < 0)),
            "M.1": int(np.sum(dist[:, 0] >= 0)),
            "neval": neval,
            "bwselect": bwselect,
            "bwparam": "main" if not is_fuzzy else bwparam,
            "vce": vce,
            "bwcheck": bwcheck,
            "masspoints": masspoints,
            "cluster": cluster,
            "clustered": cluster is not None,
            "scaleregul": scaleregul,
            "cqt": cqt,
            "level": level,
            "side": side,
            "cbands": cbands,
            "params.other": sorted(params_other_set),
            "params.cov": sorted(params_cov_set),
            "fuzzy": is_fuzzy,
        },
        tau_hat=tau_p,
        tau_hat_q=tau_q,
        se_hat=se_p,
        se_hat_q=se_q,
        params_cov=cov_tables,
        ci={"CI.l": main["ci.lower"].to_numpy(), "CI.r": main["ci.upper"].to_numpy()},
        pvalues=main["p.value"].to_numpy(),
        tvalues=main["t.value"].to_numpy(),
        tau_itt=tau_itt,
        tau_itt_q=tau_itt_q,
        tau_fs=tau_fs,
        tau_fs_q=tau_fs_q,
        rdmodel="fuzzy rd2d_distance" if is_fuzzy else "rd2d_distance",
    )


def rd2d_dist(*args, **kwargs) -> RD2DResult:
    return rd2d_distance(*args, **kwargs)
