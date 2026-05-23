from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


class AttrDict(dict):
    """Dictionary with attribute-style reads for result objects."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


class SummaryResult(AttrDict):
    """Returned by :func:`summary` and result ``summary()`` methods."""


class RD2DResult(AttrDict):
    """Result object for Python rd2d fits."""

    def summary(self, **kwargs: Any) -> SummaryResult:
        return summary(self, **kwargs)

    def __repr__(self) -> str:
        model = self.get("rdmodel", "rd2d")
        opt = self.get("opt", {})
        neval = opt.get("neval", "?")
        nobs = opt.get("N", "?")
        return f"<{model}: n={nobs}, evaluation_points={neval}>"


def _as_cov(cov: Any, n: int) -> np.ndarray | None:
    if cov is None:
        return None
    cov = np.asarray(cov, dtype=float)
    if cov.shape != (n, n):
        return None
    if not np.all(np.isfinite(cov)):
        return None
    if np.any(np.diag(cov) <= 0):
        return None
    return cov


def _normal_cval(level: float, side: str) -> float:
    if side == "two":
        return float(stats.norm.ppf((level + 100.0) / 200.0))
    return float(stats.norm.ppf(level / 100.0))


def _interval(center: float, se: float, level: float, side: str) -> tuple[float, float]:
    cval = _normal_cval(level, side)
    if not np.isfinite(center) or not np.isfinite(se):
        return (np.nan, np.nan)
    if side == "left":
        return (-np.inf, center + cval * se)
    if side == "right":
        return (center - cval * se, np.inf)
    return (center - cval * se, center + cval * se)


def _cov_to_corr(cov: np.ndarray) -> np.ndarray:
    sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    denom = np.outer(sd, sd)
    corr = np.divide(cov, denom, out=np.zeros_like(cov), where=denom > 0)
    corr = (corr + corr.T) / 2.0
    np.fill_diagonal(corr, 1.0)
    eigval, eigvec = np.linalg.eigh(corr)
    eigval = np.maximum(eigval, 1e-10)
    corr = (eigvec * eigval) @ eigvec.T
    d = np.sqrt(np.diag(corr))
    corr = corr / np.outer(d, d)
    np.fill_diagonal(corr, 1.0)
    return corr


def _simulated_cval(cov: np.ndarray, level: float, side: str, reps: int) -> float:
    corr = _cov_to_corr(cov)
    rng = np.random.default_rng(20260519)
    draws = rng.multivariate_normal(np.zeros(corr.shape[0]), corr, size=int(reps))
    if side == "two":
        stat = np.max(np.abs(draws), axis=1)
    elif side == "left":
        stat = np.max(draws, axis=1)
    else:
        stat = np.max(-draws, axis=1)
    return float(np.quantile(stat, level / 100.0, method="averaged_inverted_cdf"))


def _validate_repp(repp: int) -> int:
    try:
        value = float(repp)
    except (TypeError, ValueError):
        raise ValueError("repp must be a positive integer.") from None
    reps = int(value)
    if not np.isfinite(value) or reps < 1 or reps != value:
        raise ValueError("repp must be a positive integer.")
    return reps


def _uniform_band(est: np.ndarray, cov: np.ndarray, level: float, side: str, reps: int) -> tuple[np.ndarray, np.ndarray]:
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    cval = _simulated_cval(cov, level, side, reps)
    if side == "left":
        return np.repeat(-np.inf, len(est)), est + cval * se
    if side == "right":
        return est - cval * se, np.repeat(np.inf, len(est))
    return est - cval * se, est + cval * se


def _append_aggregate_rows(
    table: pd.DataFrame,
    cov: np.ndarray | None,
    weights: np.ndarray | None,
    lbate: bool,
    level: float,
    side: str,
    reps: int,
) -> pd.DataFrame:
    rows: list[pd.Series] = []
    base_cols = table.columns

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if len(weights) != len(table):
            raise ValueError("WBATE weights must match the number of evaluation points.")
        if not np.all(np.isfinite(weights)):
            raise ValueError("WBATE weights must be finite.")
        if np.sum(weights) == 0:
            raise ValueError("WBATE weights must have a nonzero sum.")
        weights = weights / np.sum(weights)
        row = pd.Series(np.nan, index=base_cols, name="WBATE", dtype=float)
        row["estimate.p"] = float(np.sum(weights * table["estimate.p"].to_numpy()))
        row["estimate.q"] = float(np.sum(weights * table["estimate.q"].to_numpy()))
        if cov is not None:
            se = float(np.sqrt(np.maximum(weights @ cov @ weights, 0.0)))
            row["std.err.q"] = se
            if se > 0 and np.isfinite(row["estimate.q"]):
                row["t.value"] = row["estimate.q"] / se
                row["p.value"] = 2.0 * stats.norm.sf(abs(row["t.value"]))
                row["ci.lower"], row["ci.upper"] = _interval(row["estimate.q"], se, level, side)
        rows.append(row)

    if lbate:
        row = pd.Series(np.nan, index=base_cols, name="LBATE", dtype=float)
        row["estimate.p"] = float(np.nanmax(table["estimate.p"].to_numpy()))
        row["estimate.q"] = float(np.nanmax(table["estimate.q"].to_numpy()))
        if cov is not None:
            se = np.sqrt(np.maximum(np.diag(cov), 0.0))
            cval = _simulated_cval(cov, level, side, reps)
            est_q = table["estimate.q"].to_numpy(dtype=float)
            if side == "left":
                row["ci.lower"] = -np.inf
                row["ci.upper"] = float(np.nanmax(est_q + cval * se))
            elif side == "right":
                row["ci.lower"] = float(np.nanmax(est_q - cval * se))
                row["ci.upper"] = np.inf
            else:
                row["ci.lower"] = float(np.nanmax(est_q - cval * se))
                row["ci.upper"] = float(np.nanmax(est_q + cval * se))
        rows.append(row)

    if rows:
        table = pd.concat([table, pd.DataFrame(rows)], axis=0)
    return table


def summary(
    result: RD2DResult,
    *,
    output: str | list[str] | tuple[str, ...] | None = None,
    cbands: str | list[str] | tuple[str, ...] | None = None,
    repp: int = 1000,
    WBATE: np.ndarray | list[float] | None = None,
    LBATE: bool = False,
    subset: list[int] | np.ndarray | None = None,
) -> SummaryResult:
    """Build summary tables with optional uniform bands and aggregate rows.

    Parameters
    ----------
    result : RD2DResult
        Result returned by ``rd2d()``, ``rd2d_distance()``, or ``rd2d_dist()``.
    output : str or sequence of str, optional
        Result table names to include. Defaults to ``"main"``.
    cbands : str or sequence of str, optional
        Result table names for which uniform confidence bands should be added.
        The original fit must store the matching covariance matrix.
    repp : int, default 1000
        Number of Gaussian simulation repetitions for uniform confidence bands
        and LBATE critical values.
    WBATE : array-like, optional
        Weights for a weighted boundary average treatment effect row.
    LBATE : bool, default False
        If true, appends a largest boundary average treatment effect row.
    subset : array-like of int, optional
        Zero-based row indices to display from pointwise tables. Aggregate rows
        continue to use the full boundary table.

    Returns
    -------
    SummaryResult
        Dictionary-like object with ``tables`` and ``cbands`` entries.
    """

    outputs = output
    if outputs is None:
        outputs = ["main"]
    elif isinstance(outputs, str):
        outputs = [outputs]
    else:
        outputs = list(outputs)

    cband_outputs: set[str]
    if cbands is None:
        cband_outputs = set()
    elif isinstance(cbands, str):
        cband_outputs = {cbands}
    else:
        cband_outputs = set(cbands)

    opt = result.get("opt", {})
    level = float(opt.get("level", 95.0))
    side = str(opt.get("side", "two"))
    reps = _validate_repp(repp)

    tables: dict[str, pd.DataFrame] = {}
    bands: dict[str, pd.DataFrame] = {}
    for name in outputs:
        if name not in result or not isinstance(result[name], pd.DataFrame):
            raise ValueError(f"output '{name}' is not available in this rd2d result.")
        table_all = result[name].copy()
        is_bw = name == "bw"
        cov_store = result.get("params_cov", result.get("params.cov", {}))
        cov = None if is_bw else _as_cov(cov_store.get(name), len(table_all))

        if name in cband_outputs:
            if cov is None:
                raise ValueError(f"Uniform confidence bands require params_cov='{name}' in the original fit.")
            lower, upper = _uniform_band(
                table_all["estimate.q"].to_numpy(dtype=float),
                cov,
                level,
                side,
                reps,
            )
            table_all["cb.lower"] = lower
            table_all["cb.upper"] = upper
            bands[name] = table_all[["cb.lower", "cb.upper"]].copy()

        aggregate_rows = pd.DataFrame()
        if not is_bw and (WBATE is not None or LBATE):
            aggregate_table = _append_aggregate_rows(
                table_all.copy(),
                cov,
                None if WBATE is None else np.asarray(WBATE, dtype=float),
                bool(LBATE),
                level,
                side,
                reps,
            )
            aggregate_rows = aggregate_table.iloc[len(table_all) :, :].copy()

        if subset is not None:
            idx = np.asarray(subset, dtype=int)
            if np.any(idx < 0) or np.any(idx >= len(table_all)):
                raise ValueError("subset uses zero-based row indices and must be within the result table.")
            table = table_all.iloc[idx].copy()
        else:
            table = table_all.copy()

        if not aggregate_rows.empty:
            table = pd.concat([table, aggregate_rows], axis=0)
        tables[name] = table

    return SummaryResult(model=result.get("rdmodel"), tables=tables, cbands=bands, outputs=outputs)
