from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=True)


def compare_file(reference: Path, candidate: Path) -> dict[str, object]:
    ref = read_csv(reference)
    cand = read_csv(candidate)
    out: dict[str, object] = {
        "file": reference.name,
        "status": "ok",
        "rows.reference": len(ref),
        "rows.candidate": len(cand),
        "max.abs.diff": np.nan,
        "column": "",
    }
    if not candidate.exists():
        out["status"] = "missing"
        return out
    if list(ref.columns) != list(cand.columns):
        out["status"] = "columns"
        out["column"] = "schema mismatch"
        return out
    if ref.shape != cand.shape:
        out["status"] = "shape"
        return out

    numeric_cols = [col for col in ref.columns if pd.api.types.is_numeric_dtype(ref[col]) and pd.api.types.is_numeric_dtype(cand[col])]
    max_diff = -np.inf
    max_col = ""
    for col in numeric_cols:
        diff = np.nanmax(np.abs(ref[col].to_numpy(dtype=float) - cand[col].to_numpy(dtype=float)))
        if np.isfinite(diff) and diff > max_diff:
            max_diff = float(diff)
            max_col = col
    if max_diff == -np.inf:
        max_diff = np.nan

    text_cols = [col for col in ref.columns if col not in numeric_cols]
    for col in text_cols:
        if not ref[col].fillna("").astype(str).equals(cand[col].fillna("").astype(str)):
            out["status"] = "text"
            out["column"] = col
            out["max.abs.diff"] = max_diff
            return out

    out["max.abs.diff"] = max_diff
    out["column"] = max_col
    if np.isfinite(max_diff) and max_diff > 1e-8:
        out["status"] = "diff"
    return out


def compare_dirs(reference: Path, candidate: Path) -> pd.DataFrame:
    rows = []
    for ref_file in sorted(reference.glob("empapp_*.csv")):
        rows.append(compare_file(ref_file, candidate / ref_file.name))
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    report = compare_dirs(args.reference, args.candidate)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(args.report, index=False)
    print(report.to_string(index=False))
    bad = report.loc[report["status"] != "ok"]
    if len(bad):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
