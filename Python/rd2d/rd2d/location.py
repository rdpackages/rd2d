from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from ._utils import (
    as_1d,
    as_2d,
    bandwidth_floor,
    basis_2d,
    cer_factor,
    check_lengths,
    ci_columns,
    complete_cases,
    covariance_from_influence,
    kernel_weights,
    local_fit_targets,
    normalize_binary,
    p_values,
    target_2d,
    validate_deriv,
    validate_order,
)
from .results import RD2DResult


LOCATION_MAIN_COLUMNS = [
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
    "h01",
    "h02",
    "h11",
    "h12",
    "N.Co",
    "N.Tr",
]


def _bwselect_base(bwselect: str) -> str:
    return {
        "cerrd": "mserd",
        "certwo": "msetwo",
        "icerrd": "imserd",
        "icertwo": "imsetwo",
    }.get(bwselect, bwselect)


def _is_cer(bwselect: str) -> bool:
    return bwselect in {"cerrd", "certwo", "icerrd", "icertwo"}


def _is_common(bwselect: str) -> bool:
    return _bwselect_base(bwselect) in {"mserd", "imserd"}


def _clean_location_inputs(y, x, assignment, b, fuzzy=None, cluster=None):
    y = as_1d(y, "Y")
    x = as_2d(x, "X", ncol=2)
    assignment = normalize_binary(assignment, "assignment")
    b = as_2d(b, "b", ncol=2)
    check_lengths(y, x, assignment)

    arrays = [y, x, assignment.astype(float)]
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
    x = x[keep, :]
    assignment = assignment[keep]
    if fuzzy is not None:
        fuzzy = fuzzy[keep]
    if cluster is not None:
        cluster = np.asarray(cluster)[keep]
    return y, x, assignment, b, fuzzy, cluster


def _parse_location_h(h, neval: int) -> tuple[np.ndarray, np.ndarray]:
    if h is None:
        raise ValueError("h must not be None here.")
    h = np.asarray(h, dtype=float)
    if h.ndim == 0 or h.size == 1:
        value = float(np.ravel(h)[0])
        if value <= 0:
            raise ValueError("scalar h must be positive.")
        return np.full((neval, 2), value), np.full((neval, 2), value)
    if h.shape != (neval, 4):
        raise ValueError("h must be a positive scalar or a J x 4 matrix.")
    if np.any(h <= 0):
        raise ValueError("all bandwidths in h must be positive.")
    return h[:, 0:2], h[:, 2:4]


def _qr_crossprod_inv(x: np.ndarray) -> np.ndarray:
    gram = x.T @ x
    try:
        chol = np.linalg.cholesky(gram)
        inv_chol = np.linalg.inv(chol)
        return inv_chol.T @ inv_chol
    except np.linalg.LinAlgError:
        return np.linalg.pinv(gram, rcond=1e-20)


def _basis_count_2d(p: int) -> int:
    return (p + 1) * (p + 2) // 2


def _hxy(h) -> np.ndarray:
    h = np.asarray(h, dtype=float).ravel()
    if h.size == 1:
        return np.array([h[0], h[0]], dtype=float)
    return h[:2].astype(float)


def _h_normalize(h, kernel_type: str) -> np.ndarray:
    h = np.asarray(h, dtype=float).ravel()
    if kernel_type == "prod":
        return _hxy(h)
    if h.size == 2:
        return np.array([float(np.sqrt(h[0] ** 2 + h[1] ** 2))])
    return np.array([float(h[0])])


def _exact_kernel_weights(x: np.ndarray, distance: np.ndarray, h, kernel: str, kernel_type: str) -> np.ndarray:
    h_norm = _h_normalize(h, kernel_type)
    if kernel_type == "prod":
        return kernel_weights(x[:, 0] / h_norm[0], kernel) * kernel_weights(x[:, 1] / h_norm[1], kernel) / (h_norm[0] * h_norm[1])
    return kernel_weights(distance / h_norm[0], kernel) / (h_norm[0] ** 2)


def _get_H(h, p: int) -> np.ndarray:
    hxy = _hxy(h)
    diag = []
    for deg in range(p + 1):
        for ypow in range(deg + 1):
            xpow = deg - ypow
            diag.append((hxy[0] ** xpow) * (hxy[1] ** ypow))
    return np.diag(diag)


def _get_invH(h, p: int) -> np.ndarray:
    diag = np.diag(_get_H(h, p)).astype(float)
    return np.diag(1.0 / diag)


def _vce_const(w_R: np.ndarray, resd: np.ndarray, h, cluster=None) -> np.ndarray:
    hxy = _hxy(h)
    scores = np.asarray(w_R, dtype=float) * np.asarray(resd, dtype=float)[:, None]
    if cluster is None:
        return scores.T @ scores * hxy[0] * hxy[1]

    cluster = np.asarray(cluster)
    groups = pd.unique(cluster)
    summed = np.zeros((len(groups), scores.shape[1]))
    for i, group in enumerate(groups):
        summed[i, :] = scores[cluster == group, :].sum(axis=0)
    n = len(cluster)
    k = scores.shape[1]
    g = len(groups)
    weight = ((n - 1) / (n - k)) * (g / (g - 1)) if g > 1 and n > k else 1.0
    return summed.T @ summed * hxy[0] * hxy[1] * weight


def _local_design_exact(
    x: np.ndarray,
    y: np.ndarray,
    distance: np.ndarray,
    h,
    p: int,
    kernel: str,
    kernel_type: str,
    outcomes: np.ndarray | None = None,
):
    h_norm = _h_normalize(h, kernel_type)
    hxy = _hxy(h_norm)
    w = _exact_kernel_weights(x, distance, h_norm, kernel, kernel_type)
    keep = w > 0
    ew = w[keep]
    eu = np.column_stack([x[keep, 0] / hxy[0], x[keep, 1] / hxy[1]])
    eR = basis_2d(eu, p)
    sqrt_ew = np.sqrt(ew)
    sqrtw_R = sqrt_ew[:, None] * eR
    invG = _qr_crossprod_inv(sqrtw_R)
    if outcomes is None:
        eY = y[keep].reshape(-1, 1)
    else:
        eY = np.asarray(outcomes, dtype=float)[keep, :]
    return {
        "h": h_norm,
        "hxy": hxy,
        "keep": keep,
        "eN": int(np.sum(keep)),
        "ew": ew,
        "sqrt_ew": sqrt_ew,
        "eY": eY,
        "eR": eR,
        "sqrtw_R": sqrtw_R,
        "w_R": ew[:, None] * eR,
        "k": eR.shape[1],
        "invG": invG,
        "H": _get_H(hxy, p),
        "invH": _get_invH(hxy, p),
    }


def _lm_exact(
    x: np.ndarray,
    y: np.ndarray,
    distance: np.ndarray,
    h,
    p: int,
    vce: str,
    kernel: str,
    kernel_type: str,
    cluster=None,
    varr: bool = False,
):
    loc = _local_design_exact(x, y, distance, h, p, kernel, kernel_type)
    sqrtw_Y = loc["sqrt_ew"][:, None] * loc["eY"]
    beta = loc["invH"] @ loc["invG"] @ loc["sqrtw_R"].T @ sqrtw_Y
    cov_const = None
    if varr:
        resd = loc["eY"][:, 0] - loc["eR"] @ (loc["H"] @ beta)[:, 0]
        if vce == "hc1" and loc["eN"] > loc["k"]:
            resd = resd * np.sqrt(loc["eN"] / (loc["eN"] - loc["k"]))
        elif vce in {"hc2", "hc3"}:
            hii = np.sum((loc["sqrtw_R"] @ loc["invG"]) * loc["sqrtw_R"], axis=1)
            if vce == "hc2":
                resd = resd * np.sqrt(1.0 / np.maximum(1.0 - hii, 1e-12))
            else:
                resd = resd * (1.0 / np.maximum(1.0 - hii, 1e-12))
        cov_const = loc["invG"].T @ _vce_const(loc["w_R"], resd, loc["hxy"], cluster=cluster[loc["keep"]] if cluster is not None else None) @ loc["invG"]
    return {"beta": beta[:, 0], "cov_const": cov_const, "eN": loc["eN"]}


def _lm_multi_exact(
    x: np.ndarray,
    y: np.ndarray,
    distance: np.ndarray,
    outcomes: np.ndarray,
    h,
    p: int,
    vce: str,
    kernel: str,
    kernel_type: str,
):
    loc = _local_design_exact(x, y, distance, h, p, kernel, kernel_type, outcomes=outcomes)
    sqrtw_Y = loc["sqrt_ew"][:, None] * loc["eY"]
    beta = loc["invH"] @ loc["invG"] @ loc["sqrtw_R"].T @ sqrtw_Y
    return {"beta": beta, "eN": loc["eN"]}


def _get_coeff_exact(x: np.ndarray, distance: np.ndarray, vec: np.ndarray, p: int, dn: float, kernel: str, kernel_type: str) -> np.ndarray:
    w = _exact_kernel_weights(x, distance, dn, kernel, kernel_type)
    keep = w > 0
    ew = w[keep]
    hxy = _hxy(_h_normalize(dn, kernel_type))
    eu = np.column_stack([x[keep, 0] / hxy[0], x[keep, 1] / hxy[1]])
    eR_aug = basis_2d(eu, p + 1)
    p_count = _basis_count_2d(p)
    p1_count = _basis_count_2d(p + 1)
    eR = eR_aug[:, :p_count]
    eS = eR_aug[:, p_count:p1_count]
    sqrt_ew = np.sqrt(ew)
    sqrtw_R = sqrt_ew[:, None] * eR
    sqrtw_eS = sqrt_ew[:, None] * eS
    invG = _qr_crossprod_inv(sqrtw_R)
    tail = np.asarray(vec, dtype=float).reshape(1, -1) @ invG @ sqrtw_R.T @ sqrtw_eS
    return np.concatenate([np.zeros(p_count), np.asarray(tail).ravel()])


def _bw_v2_exact(
    x: np.ndarray,
    y: np.ndarray,
    distance: np.ndarray,
    p: int,
    vec: np.ndarray,
    dn: float,
    bn1: float,
    bn2: float | None,
    vce: str,
    kernel: str,
    kernel_type: str,
    cluster=None,
) -> dict[str, float]:
    w = _exact_kernel_weights(x, distance, dn, kernel, kernel_type)
    keep = w > 0
    ew = w[keep]
    eY = y[keep]
    hxy = _hxy(_h_normalize(dn, kernel_type))
    eu = np.column_stack([x[keep, 0] / hxy[0], x[keep, 1] / hxy[1]])

    if bn2 is None:
        eR_aug = basis_2d(eu, p + 1)
        p_count = _basis_count_2d(p)
        p1_count = _basis_count_2d(p + 1)
        eR = eR_aug[:, :p_count]
        eS = eR_aug[:, p_count:p1_count]
        eT = None
    else:
        eR_aug = basis_2d(eu, p + 2)
        p_count = _basis_count_2d(p)
        p1_count = _basis_count_2d(p + 1)
        p2_count = _basis_count_2d(p + 2)
        eR = eR_aug[:, :p_count]
        eS = eR_aug[:, p_count:p1_count]
        eT = eR_aug[:, p1_count:p2_count]

    sqrt_ew = np.sqrt(ew)
    sqrtw_R = sqrt_ew[:, None] * eR
    sqrtw_eS = sqrt_ew[:, None] * eS
    sqrtw_Y = sqrt_ew * eY
    w_R = ew[:, None] * eR
    invG = _qr_crossprod_inv(sqrtw_R)
    vec = np.asarray(vec, dtype=float)
    vec_q_tail = vec.reshape(1, -1) @ invG @ sqrtw_R.T @ sqrtw_eS
    vec_q = np.concatenate([np.zeros(p_count), np.asarray(vec_q_tail).ravel()])
    vec_t = None
    if eT is not None:
        sqrtw_eT = sqrt_ew[:, None] * eT
        vec_t_tail = vec.reshape(1, -1) @ invG @ sqrtw_R.T @ sqrtw_eT
        vec_t = np.concatenate([np.zeros(p1_count), np.asarray(vec_t_tail).ravel()])

    invH = _get_invH(hxy, p)
    H = _get_H(hxy, p)
    beta = invH @ invG @ sqrtw_R.T @ sqrtw_Y.reshape(-1, 1)
    resd = eY - (eR @ (H @ beta))[:, 0]

    eN = int(np.sum(keep))
    k = eR.shape[1]
    if vce == "hc1" and eN > k:
        resd = resd * np.sqrt(eN / (eN - k))
    elif vce in {"hc2", "hc3"}:
        hii = np.sum((sqrtw_R @ invG) * sqrtw_R, axis=1)
        if vce == "hc2":
            resd = resd * np.sqrt(1.0 / np.maximum(1.0 - hii, 1e-12))
        else:
            resd = resd * (1.0 / np.maximum(1.0 - hii, 1e-12))

    sigma = _vce_const(w_R, resd, hxy, cluster=cluster[keep] if cluster is not None else None)
    V = float(np.asarray(vec.reshape(1, -1) @ invG.T @ sigma @ invG @ vec.reshape(-1, 1)).item())

    fit_p1 = _lm_exact(x, y, distance, bn1, p + 1, vce, kernel, kernel_type, cluster=cluster, varr=True)
    B = float(np.asarray(vec_q.reshape(1, -1) @ fit_p1["beta"].reshape(-1, 1)).item())
    Reg1 = float(np.asarray(vec_q.reshape(1, -1) @ fit_p1["cov_const"] @ vec_q.reshape(-1, 1) / (bn1 ** (2 + 2 * (p + 1)))).item())

    Reg2 = np.nan
    if bn2 is not None and vec_t is not None:
        fit_p2 = _lm_exact(x, y, distance, bn2, p + 2, vce, kernel, kernel_type, cluster=cluster, varr=False)
        Reg2 = float(np.asarray((dn * vec_t.reshape(1, -1)) @ fit_p2["beta"].reshape(-1, 1)).item())

    return {"B": B, "V": V, "Reg.2": Reg2, "Reg.1": Reg1}


def _rot_location(x: np.ndarray, kernel_type_label: str, M: int | None) -> float:
    constants = {
        "Epanechnikov": (1.0 / 6.0, 4.0 / (3.0 * np.pi)),
        "Triangular": (3.0 / 20.0, 3.0 / (2.0 * np.pi)),
        "Uniform": (1.0 / 4.0, 1.0 / np.pi),
        "Gaussian": (1.0, 1.0 / (4.0 * np.pi)),
    }
    mu2K_squared, l2K_squared = constants[kernel_type_label]
    cov = np.cov(x.T, ddof=1)
    inv_cov = np.linalg.pinv(cov, rcond=1e-20)
    det_cov = np.linalg.det(cov)
    sqrt_det = np.sqrt(det_cov) if np.isfinite(det_cov) and det_cov > 0 else np.sqrt(np.maximum(np.prod(np.linalg.eigvalsh(cov)), 0.0))
    trace_const = 1.0 / (2.0 ** 4 * np.pi * sqrt_det) * (
        2.0 * np.trace(inv_cov @ inv_cov) + (np.trace(inv_cov) ** 2)
    )
    n_eff = len(x) if M is None else M
    return float(((2.0 * l2K_squared) / (n_eff * mu2K_squared * trace_const)) ** (1.0 / 6.0))


def _unique_location_count(x: np.ndarray, d: np.ndarray) -> tuple[int, int, int]:
    frame = pd.DataFrame({"x1": x[:, 0], "x2": x[:, 1], "d": d.astype(int)})
    unique = frame.drop_duplicates()
    m0 = int(np.sum(unique["d"].to_numpy() == 0))
    m1 = int(np.sum(unique["d"].to_numpy() == 1))
    return m0 + m1, m0, m1


def _bwcheck_bounds_rad(x: np.ndarray, d: np.ndarray, b: np.ndarray, bwcheck: int | None) -> np.ndarray | None:
    if bwcheck is None:
        return None
    out = np.empty((len(b), 6), dtype=float)
    side0 = ~d
    side1 = d
    for j, point in enumerate(b):
        dist = np.sqrt(np.sum((x - point) ** 2, axis=1))
        vals0 = np.sort(dist[side0])
        vals1 = np.sort(dist[side1])
        if len(vals0) < bwcheck or len(vals1) < bwcheck:
            raise ValueError("bwcheck is larger than the available observations at an evaluation point.")
        out[j, :] = [vals0[bwcheck - 1], vals1[bwcheck - 1], vals0[-1], vals1[-1], len(vals0), len(vals1)]
    return out


def _fit_bwcheck_bounds(
    x: np.ndarray,
    d: np.ndarray,
    b: np.ndarray,
    bwcheck: int | None,
    kernel_type: str,
    masspoints: str,
) -> tuple[np.ndarray, np.ndarray] | None:
    if bwcheck is None:
        return None
    source_x = x
    source_d = d
    if masspoints == "adjust":
        frame = pd.DataFrame({"x1": x[:, 0], "x2": x[:, 1], "d": d.astype(int)})
        frame = frame.drop_duplicates()
        source_x = frame[["x1", "x2"]].to_numpy(dtype=float)
        source_d = frame["d"].to_numpy(dtype=int).astype(bool)

    scale = np.nanstd(x, axis=0, ddof=1)
    scale[~np.isfinite(scale) | (scale <= 0)] = 1.0
    side0 = ~source_d
    side1 = source_d
    mins = np.empty((len(b), 2, 2), dtype=float)
    maxs = np.empty((len(b), 2, 2), dtype=float)
    for j, point in enumerate(b):
        dx = source_x - point
        if kernel_type == "prod":
            dist = np.maximum(np.abs(dx[:, 0] / scale[0]), np.abs(dx[:, 1] / scale[1]))
            multiplier = scale
        else:
            dist = np.sqrt(np.sum(dx ** 2, axis=1))
            multiplier = np.ones(2)
        vals0 = np.sort(dist[side0])
        vals1 = np.sort(dist[side1])
        if len(vals0) < bwcheck or len(vals1) < bwcheck:
            raise ValueError("bwcheck is larger than the available observations at an evaluation point.")
        min0 = vals0[bwcheck - 1] * multiplier
        min1 = vals1[bwcheck - 1] * multiplier
        max0 = vals0[-1] * multiplier
        max1 = vals1[-1] * multiplier
        mins[j, 0, :] = min0
        mins[j, 1, :] = min1
        maxs[j, 0, :] = max0
        maxs[j, 1, :] = max1
    return mins, maxs


def _apply_fit_bwcheck(
    h0: np.ndarray,
    h1: np.ndarray,
    x: np.ndarray,
    d: np.ndarray,
    b: np.ndarray,
    bwcheck: int | None,
    kernel_type: str,
    masspoints: str,
) -> tuple[np.ndarray, np.ndarray]:
    bounds = _fit_bwcheck_bounds(x, d, b, bwcheck, kernel_type, masspoints)
    if bounds is None:
        return h0, h1
    mins, maxs = bounds
    h0_out = np.minimum(np.maximum(h0, mins[:, 0, :]), maxs[:, 0, :])
    h1_out = np.minimum(np.maximum(h1, mins[:, 1, :]), maxs[:, 1, :])
    return h0_out, h1_out


def rdbw2d(
    Y,
    X,
    assignment,
    b,
    *,
    p: int = 1,
    deriv=(0, 0),
    tangvec=None,
    kernel: str = "tri",
    kernel_type: str = "prod",
    bwselect: str = "mserd",
    bwparam: str = "main",
    method: str = "dpi",
    vce: str = "hc1",
    bwcheck: int | None = 20,
    masspoints: str = "check",
    cluster=None,
    scaleregul: float = 1.0,
    scalebiascrct: float = 1.0,
    stdvars: bool = True,
    fuzzy=None,
) -> RD2DResult:
    """Bandwidth selection for location-based boundary discontinuity designs.

    Parameters are accepted as NumPy arrays, pandas objects, or array-like
    values. The returned object contains a ``bws`` table with one row per
    boundary evaluation point and columns ``b1``, ``b2``, ``h01``, ``h02``,
    ``h11``, ``h12``, ``N.Co``, and ``N.Tr``.
    """

    y, x, d, b, fuzzy, cluster = _clean_location_inputs(Y, X, assignment, b, fuzzy, cluster)
    p = validate_order(p, "p")
    deriv = validate_deriv(deriv, p)
    tangvec_arr = None if tangvec is None else as_2d(tangvec, "tangvec", ncol=2)
    if tangvec_arr is not None and tangvec_arr.shape[0] != b.shape[0]:
        raise ValueError("tangvec must have the same number of rows as b.")
    kernel_type = kernel_type.lower()
    if kernel_type not in {"prod", "rad"}:
        raise ValueError("kernel_type must be 'prod' or 'rad'.")
    bwselect = bwselect.lower()
    base_select = _bwselect_base(bwselect)
    if base_select not in {"mserd", "msetwo", "imserd", "imsetwo"}:
        raise ValueError("unsupported bwselect.")

    n = len(y)
    n0 = int(np.sum(~d))
    n1 = int(np.sum(d))
    neval = b.shape[0]
    kernel_label = "Epanechnikov"
    if kernel.lower() in {"tri", "triangular"}:
        kernel_label = "Triangular"
    elif kernel.lower() in {"uni", "uniform"}:
        kernel_label = "Uniform"
    elif kernel.lower() in {"gau", "gaussian"}:
        kernel_label = "Gaussian"

    x_work = x.copy()
    b_work = b.copy()
    if stdvars:
        sd = np.nanstd(x_work, axis=0, ddof=1)
        sd[~np.isfinite(sd) | (sd <= 0)] = 1.0
        x_work = x_work / sd
        b_work = b_work / sd
    else:
        sd = np.ones(2)

    M, M0, M1 = (n, n0, n1)
    if masspoints in {"check", "adjust"}:
        M, M0, M1 = _unique_location_count(x_work, d)

    dn = _rot_location(x_work, kernel_label, M)
    bwcheck_bounds = _bwcheck_bounds_rad(x_work, d, b_work, bwcheck)

    deriv_sum = int(sum(deriv))
    deriv_denom = 2 * p + 2 - 2 * deriv_sum
    e_deriv = np.zeros((_basis_count_2d(p),), dtype=float)
    if tangvec_arr is not None:
        e_deriv[1] = tangvec_arr[0, 0]
        e_deriv[2] = tangvec_arr[0, 1]
        deriv_sum = 1
        deriv_denom = 2 * p
    elif deriv_sum >= 1:
        indices = [(a, bpow) for deg in range(p + 1) for bpow in range(deg + 1) for a in [deg - bpow]]
        pos = indices.index(tuple(deriv))
        e_deriv[pos] = float(math.factorial(deriv[0]) * math.factorial(deriv[1]))
    else:
        e_deriv[0] = 1.0

    rows: list[list[float]] = []
    mse_rows: list[list[float]] = []
    for j, point in enumerate(b_work):
        vec = e_deriv.copy()
        if tangvec_arr is not None:
            vec[1] = tangvec_arr[j, 0]
            vec[2] = tangvec_arr[j, 1]

        centered = x_work - point
        radial = np.sqrt(np.sum(centered ** 2, axis=1))
        dn0 = dn1 = dn
        if bwcheck_bounds is not None:
            bw_min0, bw_min1, bw_max0, bw_max1 = bwcheck_bounds[j, 0], bwcheck_bounds[j, 1], bwcheck_bounds[j, 2], bwcheck_bounds[j, 3]
            dn0 = min(max(dn, bw_min0), bw_max0)
            dn1 = min(max(dn, bw_min1), bw_max1)
        else:
            bw_min0 = bw_min1 = 0.0
            bw_max0 = bw_max1 = np.inf

        y_bw = y.copy()
        if fuzzy is not None and bwparam == "main":
            outcomes = np.column_stack([y, fuzzy])
            side0 = ~d
            side1 = d
            fit_grad0 = _lm_multi_exact(centered[side0], y[side0], radial[side0], outcomes[side0], dn0, p, vce, kernel, kernel_type)
            fit_grad1 = _lm_multi_exact(centered[side1], y[side1], radial[side1], outcomes[side1], dn1, p, vce, kernel, kernel_type)
            tau_itt_grad = float(vec @ fit_grad1["beta"][:, 0] - vec @ fit_grad0["beta"][:, 0])
            tau_fs_grad = float(vec @ fit_grad1["beta"][:, 1] - vec @ fit_grad0["beta"][:, 1])
            if np.isfinite(tau_fs_grad) and abs(tau_fs_grad) > np.sqrt(np.finfo(float).eps):
                grad_itt = 1.0 / tau_fs_grad
                grad_fs = -tau_itt_grad / (tau_fs_grad ** 2)
                y_bw = grad_itt * y + grad_fs * fuzzy

        side0 = ~d
        side1 = d
        w0 = _exact_kernel_weights(centered[side0], radial[side0], dn0, kernel, kernel_type)
        w1 = _exact_kernel_weights(centered[side1], radial[side1], dn1, kernel, kernel_type)
        eN0 = int(np.sum(w0 > 0))
        eN1 = int(np.sum(w1 > 0))

        vec_q0 = _get_coeff_exact(centered[side0], radial[side0], vec, p, dn0, kernel, kernel_type)
        vec_q1 = _get_coeff_exact(centered[side1], radial[side1], vec, p, dn1, kernel, kernel_type)
        thr0 = float(np.median(radial[side0]))
        thr1 = float(np.median(radial[side1]))
        bn0 = thr0
        bn1v = thr1

        if method == "dpi":
            bn_const0 = _bw_v2_exact(centered[side0], y_bw[side0], radial[side0], p + 1, vec_q0, dn0, thr0, None, vce, kernel, kernel_type)
            bn_const1 = _bw_v2_exact(centered[side1], y_bw[side1], radial[side1], p + 1, vec_q1, dn1, thr1, None, vce, kernel, kernel_type)
            bn0 = ((2 + 2 * (p + 1)) * bn_const0["V"] / (2 * (bn_const0["B"] ** 2 + scaleregul * bn_const0["Reg.1"]))) ** (1.0 / (2 * p + 6))
            bn1v = ((2 + 2 * (p + 1)) * bn_const1["V"] / (2 * (bn_const1["B"] ** 2 + scaleregul * bn_const1["Reg.1"]))) ** (1.0 / (2 * p + 6))
            if bwcheck_bounds is not None:
                bn0 = min(max(bn0, bw_min0), bw_max0)
                bn1v = min(max(bn1v, bw_min1), bw_max1)

        hn_const0 = _bw_v2_exact(centered[side0], y_bw[side0], radial[side0], p, vec, dn0, bn0, thr0, vce, kernel, kernel_type)
        hn_const1 = _bw_v2_exact(centered[side1], y_bw[side1], radial[side1], p, vec, dn1, bn1v, thr1, vce, kernel, kernel_type)

        if base_select in {"mserd", "imserd"}:
            denom = deriv_denom * (
                (hn_const0["B"] + scalebiascrct * hn_const0["Reg.2"] - hn_const1["B"] - scalebiascrct * hn_const1["Reg.2"]) ** 2
                + scaleregul * hn_const0["Reg.1"]
                + scaleregul * hn_const1["Reg.1"]
            )
            hn = ((2 + 2 * deriv_sum) * (hn_const0["V"] + hn_const1["V"]) / denom) ** (1.0 / (2 * p + 4))
            if bwcheck_bounds is not None:
                hn = min(max(hn, bw_min0, bw_min1), max(bw_max0, bw_max1))
            hn0 = hn1 = hn
        else:
            hn0 = ((2 + 2 * deriv_sum) * hn_const0["V"] / (
                deriv_denom * ((hn_const0["B"] + scalebiascrct * hn_const0["Reg.2"]) ** 2 + scaleregul * hn_const0["Reg.1"])
            )) ** (1.0 / (2 * p + 4))
            hn1 = ((2 + 2 * deriv_sum) * hn_const1["V"] / (
                deriv_denom * ((hn_const1["B"] + scalebiascrct * hn_const1["Reg.2"]) ** 2 + scaleregul * hn_const1["Reg.1"])
            )) ** (1.0 / (2 * p + 4))
            if bwcheck_bounds is not None:
                hn0 = min(max(hn0, bw_min0), bw_max0)
                hn1 = min(max(hn1, bw_min1), bw_max1)

        rows.append([b[j, 0], b[j, 1], hn0 * sd[0], hn0 * sd[1], hn1 * sd[0], hn1 * sd[1], eN0, eN1])
        mse_rows.append([
            b_work[j, 0],
            b_work[j, 1],
            hn0,
            hn0,
            hn1,
            hn1,
            eN0,
            eN1,
            hn_const0["B"],
            hn_const1["B"],
            hn_const0["V"],
            hn_const1["V"],
            hn_const0["Reg.2"],
            hn_const1["Reg.2"],
            hn_const0["Reg.1"],
            hn_const1["Reg.1"],
        ])

    mseconsts = pd.DataFrame(
        mse_rows,
        columns=[
            "b1",
            "b2",
            "h01",
            "h02",
            "h11",
            "h12",
            "N.Co",
            "N.Tr",
            "bias.0",
            "bias.1",
            "var.0",
            "var.1",
            "reg.bias.0",
            "reg.bias.1",
            "reg.var.0",
            "reg.var.1",
        ],
    )

    if base_select == "imserd":
        V = float(np.nanmean(mseconsts["var.0"]) + np.nanmean(mseconsts["var.1"]))
        B = float(np.nanmean(
            (mseconsts["bias.0"] + scalebiascrct * mseconsts["reg.bias.0"] - mseconsts["bias.1"] - scalebiascrct * mseconsts["reg.bias.1"]) ** 2
            + scaleregul * mseconsts["reg.var.0"]
            + scaleregul * mseconsts["reg.var.1"]
        ))
        h_imse = ((2 + 2 * deriv_sum) * V / (deriv_denom * B)) ** (1.0 / (2 * p + 4))
        for row in rows:
            row[2:6] = [h_imse * sd[0], h_imse * sd[1], h_imse * sd[0], h_imse * sd[1]]
    elif base_select == "imsetwo":
        V0 = float(np.nanmean(mseconsts["var.0"]))
        V1 = float(np.nanmean(mseconsts["var.1"]))
        B0 = float(np.nanmean((mseconsts["bias.0"] + scalebiascrct * mseconsts["reg.bias.0"]) ** 2 + scaleregul * mseconsts["reg.var.0"]))
        B1 = float(np.nanmean((mseconsts["bias.1"] + scalebiascrct * mseconsts["reg.bias.1"]) ** 2 + scaleregul * mseconsts["reg.var.1"]))
        h0_imse = ((2 + 2 * deriv_sum) * V0 / (deriv_denom * B0)) ** (1.0 / (2 * p + 4))
        h1_imse = ((2 + 2 * deriv_sum) * V1 / (deriv_denom * B1)) ** (1.0 / (2 * p + 4))
        for row in rows:
            row[2:6] = [h0_imse * sd[0], h0_imse * sd[1], h1_imse * sd[0], h1_imse * sd[1]]

    if _is_cer(bwselect):
        if _is_common(bwselect):
            factor = cer_factor(M, p)
            for row in rows:
                row[2] *= factor
                row[3] *= factor
                row[4] *= factor
                row[5] *= factor
        else:
            factor0 = cer_factor(M0, p)
            factor1 = cer_factor(M1, p)
            for row in rows:
                row[2] *= factor0
                row[3] *= factor0
                row[4] *= factor1
                row[5] *= factor1

    bws = pd.DataFrame(rows, columns=["b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr"])
    out = RD2DResult(
        bws=bws,
        mseconsts=mseconsts,
        opt={
            "N": n,
            "N.0": n0,
            "N.1": n1,
            "M.0": n0,
            "M.1": n1,
            "neval": neval,
            "p": p,
            "deriv": deriv,
            "tangvec": tangvec_arr,
            "kernel": kernel,
            "kernel_type": kernel_type,
            "bwselect": bwselect,
            "bwparam": "main" if fuzzy is None else bwparam,
            "method": method,
            "bwcheck": bwcheck,
            "stdvars": stdvars,
            "fuzzy": fuzzy is not None,
            "cluster": cluster,
            "clustered": cluster is not None,
            "vce": vce,
            "masspoints": masspoints,
            "scaleregul": scaleregul,
            "scalebiascrct": scalebiascrct,
        },
        rdmodel="rdbw2d",
    )
    return out


def _location_weights(centered: np.ndarray, h: np.ndarray, kernel: str, kernel_type: str) -> np.ndarray:
    if kernel_type == "prod":
        w = kernel_weights(centered[:, 0] / h[0], kernel) * kernel_weights(centered[:, 1] / h[1], kernel)
    else:
        u = np.sqrt((centered[:, 0] / h[0]) ** 2 + (centered[:, 1] / h[1]) ** 2)
        w = kernel_weights(u, kernel)
    return w / max(float(h[0] * h[1]), np.finfo(float).eps)


def _fit_location_order(
    y: np.ndarray,
    x: np.ndarray,
    d: np.ndarray,
    b: np.ndarray,
    outcomes: np.ndarray,
    h0: np.ndarray,
    h1: np.ndarray,
    p: int,
    deriv: tuple[int, int],
    tangvec: np.ndarray | None,
    kernel: str,
    kernel_type: str,
    vce: str,
    cluster,
) -> dict[str, Any]:
    nout = outcomes.shape[1]
    neval = b.shape[0]
    mu0 = np.full((neval, nout), np.nan)
    mu1 = np.full((neval, nout), np.nan)
    se0 = np.full((neval, nout), np.nan)
    se1 = np.full((neval, nout), np.nan)
    n0 = np.zeros(neval, dtype=int)
    n1 = np.zeros(neval, dtype=int)
    infl0 = [[None for _ in range(nout)] for _ in range(neval)]
    infl1 = [[None for _ in range(nout)] for _ in range(neval)]

    for j, point in enumerate(b):
        centered = x - point
        design = basis_2d(centered, p)
        target = target_2d(p, deriv, tangvec, j)
        w0 = _location_weights(centered, h0[j, :], kernel, kernel_type)
        w1 = _location_weights(centered, h1[j, :], kernel, kernel_type)
        w0 = np.where(~d, w0, 0.0)
        w1 = np.where(d, w1, 0.0)
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

    return {
        "mu0": mu0,
        "mu1": mu1,
        "se0": se0,
        "se1": se1,
        "N0": n0,
        "N1": n1,
        "infl0": infl0,
        "infl1": infl1,
    }


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


def _se_from_influence(infl: np.ndarray | None) -> np.ndarray:
    if infl is None:
        return np.repeat(np.nan, 0)
    return np.sqrt(np.maximum(np.diag(covariance_from_influence(infl)), 0.0))


def _make_table(
    b: np.ndarray,
    est_p: np.ndarray,
    se_p: np.ndarray,
    est_q: np.ndarray,
    se_q: np.ndarray,
    h0: np.ndarray,
    h1: np.ndarray,
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
            h0[:, 0],
            h0[:, 1],
            h1[:, 0],
            h1[:, 1],
            n0,
            n1,
        ]
    )
    return pd.DataFrame(data, columns=LOCATION_MAIN_COLUMNS)


def rd2d(
    Y,
    X,
    assignment,
    b,
    *,
    h=None,
    deriv=(0, 0),
    tangvec=None,
    p: int = 1,
    q: int | None = None,
    kernel: str = "tri",
    kernel_type: str = "prod",
    vce: str = "hc1",
    masspoints: str = "check",
    cluster=None,
    level: float = 95,
    params_other: list[str] | tuple[str, ...] | str | None = None,
    params_cov: list[str] | tuple[str, ...] | str | None = None,
    side: str = "two",
    bwselect: str = "mserd",
    bwparam: str = "main",
    method: str = "dpi",
    bwcheck: int | None = None,
    scaleregul: float = 3.0,
    scalebiascrct: float = 1.0,
    stdvars: bool = True,
    fuzzy=None,
) -> RD2DResult:
    """Location-based local polynomial boundary discontinuity estimator.

    Parameters
    ----------
    Y : array-like
        Outcome variable.
    X : array-like, shape (n_obs, 2)
        Bivariate running variable.
    assignment : array-like
        Treatment assignment indicator. Values must be logical or binary.
    b : array-like, shape (n_eval, 2)
        Boundary evaluation points.
    h : float or array-like, optional
        User-supplied bandwidths. A scalar uses the same bandwidth for every
        coordinate, side, and evaluation point. A matrix must have columns
        ``h01``, ``h02``, ``h11``, and ``h12``.
    fuzzy : array-like, optional
        Treatment receipt/status variable for fuzzy designs.
    params_cov : str or sequence of str, optional
        Result tables for which covariance matrices should be stored. Stored
        covariance enables uniform confidence bands, WBATE, and LBATE in
        ``summary()``.

    Returns
    -------
    RD2DResult
        Result object with ``main`` and ``bw`` tables. Fuzzy fits also include
        ``itt`` and ``fs`` tables; requested side-specific tables are included
        when available.
    """

    y, x, d, b, fuzzy, cluster = _clean_location_inputs(Y, X, assignment, b, fuzzy, cluster)
    p = validate_order(p, "p")
    if q is None:
        q = p + 1
    q = validate_order(q, "q")
    if q < p:
        raise ValueError("q must be greater than or equal to p.")
    deriv = validate_deriv(deriv, p)
    if sum(deriv) > q:
        raise ValueError("sum(deriv) must be less than or equal to q.")
    tangvec_arr = None if tangvec is None else as_2d(tangvec, "tangvec", ncol=2)
    if tangvec_arr is not None and tangvec_arr.shape[0] != b.shape[0]:
        raise ValueError("tangvec must have the same number of rows as b.")
    kernel_type = kernel_type.lower()
    if kernel_type not in {"prod", "rad"}:
        raise ValueError("kernel_type must be 'prod' or 'rad'.")
    side = side.lower()
    if side not in {"two", "left", "right"}:
        raise ValueError("side must be 'two', 'left', or 'right'.")
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

    is_fuzzy = fuzzy is not None
    fit_bwcheck = 50 + p + 1 if bwcheck is None else bwcheck
    if h is None:
        bws = rdbw2d(
            y,
            x,
            d.astype(float),
            b,
            p=p,
            deriv=deriv,
            tangvec=tangvec_arr,
            kernel=kernel,
            kernel_type=kernel_type,
            bwselect=bwselect,
            bwparam=bwparam,
            method=method,
            vce=vce,
            bwcheck=fit_bwcheck,
            masspoints=masspoints,
            cluster=cluster,
            scaleregul=scaleregul,
            scalebiascrct=scalebiascrct,
            stdvars=stdvars,
            fuzzy=fuzzy,
        ).bws
        h0 = bws[["h01", "h02"]].to_numpy(dtype=float)
        h1 = bws[["h11", "h12"]].to_numpy(dtype=float)
    else:
        bwselect = "user provided"
        h0, h1 = _parse_location_h(h, b.shape[0])

    h0, h1 = _apply_fit_bwcheck(h0, h1, x, d, b, fit_bwcheck, kernel_type, masspoints)

    outcomes = y.reshape(-1, 1) if not is_fuzzy else np.column_stack([y, fuzzy])
    fit_p = _fit_location_order(y, x, d, b, outcomes, h0, h1, p, deriv, tangvec_arr, kernel, kernel_type, vce, cluster)
    fit_q = fit_p if q == p else _fit_location_order(y, x, d, b, outcomes, h0, h1, q, deriv, tangvec_arr, kernel, kernel_type, vce, cluster)

    infl_tables: dict[str, np.ndarray] = {}
    if not is_fuzzy:
        tau_p = fit_p["mu1"][:, 0] - fit_p["mu0"][:, 0]
        tau_q = fit_q["mu1"][:, 0] - fit_q["mu0"][:, 0]
        infl_main_p = _contrast_influence(fit_p, 0)
        infl_main_q = _contrast_influence(fit_q, 0)
        se_p = _se_from_influence(infl_main_p)
        se_q = _se_from_influence(infl_main_q)
        main = _make_table(b, tau_p, se_p, tau_q, se_q, h0, h1, fit_p["N0"], fit_p["N1"], level, side)
        bw = main[["b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr"]].copy()
        main0 = main1 = None
        if "main.0" in params_other_set:
            main0 = _make_table(
                b,
                fit_p["mu0"][:, 0],
                fit_p["se0"][:, 0],
                fit_q["mu0"][:, 0],
                fit_q["se0"][:, 0],
                h0,
                np.full_like(h1, np.nan),
                fit_p["N0"],
                np.repeat(np.nan, len(b)),
                level,
                side,
            )
        if "main.1" in params_other_set:
            main1 = _make_table(
                b,
                fit_p["mu1"][:, 0],
                fit_p["se1"][:, 0],
                fit_q["mu1"][:, 0],
                fit_q["se1"][:, 0],
                np.full_like(h0, np.nan),
                h1,
                np.repeat(np.nan, len(b)),
                fit_p["N1"],
                level,
                side,
            )
        if infl_main_q is not None:
            infl_tables["main"] = infl_main_q
        infl0 = _side_influence(fit_q, 0, 0)
        infl1 = _side_influence(fit_q, 0, 1)
        if infl0 is not None:
            infl_tables["main.0"] = infl0
        if infl1 is not None:
            infl_tables["main.1"] = infl1
        result_tables = {"main": main, "bw": bw, "main.0": main0, "main.1": main1}
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
        main = _make_table(b, tau_p, se_p, tau_q, se_q, h0, h1, fit_p["N0"], fit_p["N1"], level, side)
        itt = _make_table(b, itt_p, se_itt_p, itt_q, se_itt_q, h0, h1, fit_p["N0"], fit_p["N1"], level, side)
        fs = _make_table(b, fs_p, se_fs_p, fs_q, se_fs_q, h0, h1, fit_p["N0"], fit_p["N1"], level, side)
        bw = main[["b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr"]].copy()
        result_tables = {"main": main, "bw": bw, "itt": itt, "itt.0": None, "itt.1": None, "fs": fs, "fs.0": None, "fs.1": None}
        infl_tables.update({"main": infl_main_q, "itt": infl_itt_q, "fs": infl_fs_q})
        side_map = {
            "itt.0": (0, 0),
            "itt.1": (0, 1),
            "fs.0": (1, 0),
            "fs.1": (1, 1),
        }
        for name, (outcome, side_id) in side_map.items():
            if name in params_other_set:
                mu_p = fit_p["mu1" if side_id == 1 else "mu0"][:, outcome]
                mu_q = fit_q["mu1" if side_id == 1 else "mu0"][:, outcome]
                se_side_p = fit_p["se1" if side_id == 1 else "se0"][:, outcome]
                se_side_q = fit_q["se1" if side_id == 1 else "se0"][:, outcome]
                h0_side = h0 if side_id == 0 else np.full_like(h0, np.nan)
                h1_side = h1 if side_id == 1 else np.full_like(h1, np.nan)
                n0_side = fit_p["N0"] if side_id == 0 else np.repeat(np.nan, len(b))
                n1_side = fit_p["N1"] if side_id == 1 else np.repeat(np.nan, len(b))
                result_tables[name] = _make_table(b, mu_p, se_side_p, mu_q, se_side_q, h0_side, h1_side, n0_side, n1_side, level, side)
            infl = _side_influence(fit_q, outcome, side_id)
            if infl is not None:
                infl_tables[name] = infl
        tau_itt, tau_itt_q, tau_fs, tau_fs_q = itt_p, itt_q, fs_p, fs_q

    cov_tables = {
        name: covariance_from_influence(infl)
        for name, infl in infl_tables.items()
        if name in params_cov_set and infl is not None
    }

    out = RD2DResult(
        **result_tables,
        opt={
            "b": b,
            "deriv": deriv,
            "tangvec": tangvec_arr,
            "p": p,
            "q": q,
            "kernel": kernel,
            "kernel_type": kernel_type,
            "N": len(y),
            "N.0": int(np.sum(~d)),
            "N.1": int(np.sum(d)),
            "M": len(y),
            "M.0": int(np.sum(~d)),
            "M.1": int(np.sum(d)),
            "neval": len(b),
            "bwselect": bwselect,
            "bwparam": "main" if not is_fuzzy else bwparam,
            "method": method,
            "vce": vce,
            "bwcheck": bwcheck,
            "masspoints": masspoints,
            "cluster": cluster,
            "clustered": cluster is not None,
            "scaleregul": scaleregul,
            "scalebiascrct": scalebiascrct,
            "stdvars": stdvars,
            "fuzzy": is_fuzzy,
            "level": level,
            "side": side,
            "params.other": sorted(params_other_set),
            "params.cov": sorted(params_cov_set),
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
        rdmodel="fuzzy rd2d" if is_fuzzy else "rd2d",
    )
    return out
