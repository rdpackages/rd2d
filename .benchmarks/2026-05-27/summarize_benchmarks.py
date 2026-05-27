from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


KEYS = ["case", "call", "output", "row_norm", "column_norm"]
RESULT_COLUMNS = [
    "comparison",
    "lhs",
    "rhs",
    "n_common",
    "n_lhs_only",
    "n_rhs_only",
    "max_abs_diff",
    "max_rel_diff",
    "p99_abs_diff",
    "p99_rel_diff",
    "median_abs_diff",
    "n_gt_1e_12",
    "n_gt_1e_10",
    "n_gt_1e_8",
    "worst_case",
    "worst_output",
    "worst_row",
    "worst_column",
    "lhs_value",
    "rhs_value",
]


def norm_label(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    match = re.fullmatch(r"[rc](\d+)", text)
    if match:
        return str(int(match.group(1)))
    match = re.fullmatch(r"(\d+)\.0", text)
    if match:
        return str(int(match.group(1)))
    return text


def shift_zero_based(labels: pd.Series) -> pd.Series:
    normalized = labels.map(norm_label)
    nums = pd.to_numeric(normalized, errors="coerce")
    mask = nums.notna()
    if mask.any() and mask.all() and nums.min() == 0:
        normalized = (nums.astype(int) + 1).astype(str)
    return normalized


def normalize_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"row": str, "column": str})
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["column_norm"] = df["column"].map(norm_label).str.replace("_", ".", regex=False)
    df["row_norm"] = (
        df.groupby(["platform", "source", "case", "call", "output"], group_keys=False)["row"]
        .apply(shift_zero_based)
        .astype(str)
    )
    df["column_norm"] = (
        df.groupby(["platform", "source", "case", "call", "output"], group_keys=False)["column_norm"]
        .apply(shift_zero_based)
        .astype(str)
    )
    return df


def load_results(result_dir: Path) -> dict[tuple[str, str], pd.DataFrame]:
    results: dict[tuple[str, str], pd.DataFrame] = {}
    for path in sorted(result_dir.glob("*_results.csv")):
        parts = path.stem.split("_")
        platform = parts[0]
        source = "_".join(parts[1:-1])
        df = normalize_results(path)
        df["dataset"] = f"{platform}:{source}"
        results[(platform, source)] = df
    return results


def compare_frames(
    frames: dict[tuple[str, str], pd.DataFrame],
    lhs_key: tuple[str, str],
    rhs_key: tuple[str, str],
    comparison: str,
) -> tuple[dict[str, object], pd.DataFrame]:
    lhs = frames[lhs_key].copy()
    rhs = frames[rhs_key].copy()
    lhs_name = f"{lhs_key[0]}:{lhs_key[1]}"
    rhs_name = f"{rhs_key[0]}:{rhs_key[1]}"

    lhs = lhs[KEYS + ["value"]].rename(columns={"value": "lhs_value"})
    rhs = rhs[KEYS + ["value"]].rename(columns={"value": "rhs_value"})
    merged = lhs.merge(rhs, on=KEYS, how="outer", indicator=True)

    common = merged[merged["_merge"] == "both"].copy()
    common["abs_diff"] = (common["lhs_value"] - common["rhs_value"]).abs()
    denom = np.maximum(1.0, np.maximum(common["lhs_value"].abs(), common["rhs_value"].abs()))
    common["rel_diff"] = common["abs_diff"] / denom
    finite = common[np.isfinite(common["abs_diff"])].copy()

    row: dict[str, object] = {
        "comparison": comparison,
        "lhs": lhs_name,
        "rhs": rhs_name,
        "n_common": int(len(common)),
        "n_lhs_only": int((merged["_merge"] == "left_only").sum()),
        "n_rhs_only": int((merged["_merge"] == "right_only").sum()),
        "max_abs_diff": np.nan,
        "max_rel_diff": np.nan,
        "p99_abs_diff": np.nan,
        "p99_rel_diff": np.nan,
        "median_abs_diff": np.nan,
        "n_gt_1e_12": 0,
        "n_gt_1e_10": 0,
        "n_gt_1e_8": 0,
        "worst_case": "",
        "worst_output": "",
        "worst_row": "",
        "worst_column": "",
        "lhs_value": np.nan,
        "rhs_value": np.nan,
    }

    if not finite.empty:
        idx = finite["abs_diff"].idxmax()
        worst = finite.loc[idx]
        row.update(
            {
                "max_abs_diff": float(finite["abs_diff"].max()),
                "max_rel_diff": float(finite["rel_diff"].max()),
                "p99_abs_diff": float(finite["abs_diff"].quantile(0.99)),
                "p99_rel_diff": float(finite["rel_diff"].quantile(0.99)),
                "median_abs_diff": float(finite["abs_diff"].median()),
                "n_gt_1e_12": int((finite["abs_diff"] > 1e-12).sum()),
                "n_gt_1e_10": int((finite["abs_diff"] > 1e-10).sum()),
                "n_gt_1e_8": int((finite["abs_diff"] > 1e-8).sum()),
                "worst_case": worst["case"],
                "worst_output": worst["output"],
                "worst_row": worst["row_norm"],
                "worst_column": worst["column_norm"],
                "lhs_value": float(worst["lhs_value"]),
                "rhs_value": float(worst["rhs_value"]),
            }
        )

    finite.insert(0, "rhs", rhs_name)
    finite.insert(0, "lhs", lhs_name)
    finite.insert(0, "comparison", comparison)
    return row, finite


def build_comparisons(frames: dict[tuple[str, str], pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs: list[tuple[str, tuple[str, str], tuple[str, str]]] = []

    for platform in ["r", "python", "stata"]:
        pairs.append(("version_public_vs_local_separate", (platform, "public"), (platform, "local_separate")))

    for source in ["public", "local_separate"]:
        for lhs, rhs in [("r", "python"), ("r", "stata"), ("python", "stata")]:
            pairs.append((f"cross_platform_{source}", (lhs, source), (rhs, source)))

    for platform in ["r", "python", "stata"]:
        pairs.append(("local_joint_vs_separate", (platform, "local_separate"), (platform, "local_joint")))

    rows: list[dict[str, object]] = []
    diffs: list[pd.DataFrame] = []
    for comparison, lhs_key, rhs_key in pairs:
        if lhs_key not in frames or rhs_key not in frames:
            continue
        row, diff = compare_frames(frames, lhs_key, rhs_key, comparison)
        rows.append(row)
        diffs.append(diff)

    summary = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    all_diffs = pd.concat(diffs, ignore_index=True) if diffs else pd.DataFrame()
    if not all_diffs.empty:
        all_diffs = all_diffs.sort_values("abs_diff", ascending=False)
    return summary, all_diffs


def load_timings(result_dir: Path) -> pd.DataFrame:
    frames = []
    for path in sorted(result_dir.glob("*_timings.csv")):
        df = pd.read_csv(path)
        if df.empty:
            continue
        df["elapsed_seconds"] = pd.to_numeric(df["elapsed_seconds"], errors="coerce")
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def timing_summary(timings: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if timings.empty:
        return pd.DataFrame(), pd.DataFrame()

    grouped = (
        timings.groupby(["platform", "source", "version", "fitmethod", "case"], dropna=False)
        .agg(
            n=("elapsed_seconds", "count"),
            median_seconds=("elapsed_seconds", "median"),
            min_seconds=("elapsed_seconds", "min"),
            max_seconds=("elapsed_seconds", "max"),
        )
        .reset_index()
    )

    ratio_rows = []
    for platform, dfp in grouped.groupby("platform"):
        pub = dfp[dfp["source"] == "public"][["case", "median_seconds"]].rename(
            columns={"median_seconds": "public_median_seconds"}
        )
        for source in ["local_separate", "local_joint"]:
            loc = dfp[dfp["source"] == source][["case", "median_seconds"]].rename(
                columns={"median_seconds": "local_median_seconds"}
            )
            merged = loc.merge(pub, on="case", how="inner")
            if merged.empty:
                continue
            merged.insert(0, "platform", platform)
            merged.insert(1, "source", source)
            merged["local_over_public"] = merged["local_median_seconds"] / merged["public_median_seconds"]
            merged["speedup_public_over_local"] = merged["public_median_seconds"] / merged["local_median_seconds"]
            ratio_rows.append(merged)
    ratios = pd.concat(ratio_rows, ignore_index=True) if ratio_rows else pd.DataFrame()
    return grouped, ratios


def main() -> int:
    repo_root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    bench_dir = repo_root / ".benchmarks" / "2026-05-27"
    result_dir = bench_dir / "results"

    frames = load_results(result_dir)
    summary, diffs = build_comparisons(frames)
    noncluster_frames = {
        key: df[~df["case"].str.contains("cluster", regex=False)].copy()
        for key, df in frames.items()
    }
    summary_noncluster, diffs_noncluster = build_comparisons(noncluster_frames)
    timings = load_timings(result_dir)
    time_summary, time_ratios = timing_summary(timings)

    summary.to_csv(result_dir / "comparison_summary.csv", index=False)
    summary_noncluster.to_csv(result_dir / "comparison_summary_noncluster.csv", index=False)
    diffs.head(250).to_csv(result_dir / "worst_diffs_top250.csv", index=False)
    diffs_noncluster.head(250).to_csv(result_dir / "worst_diffs_noncluster_top250.csv", index=False)
    time_summary.to_csv(result_dir / "timing_summary.csv", index=False)
    time_ratios.to_csv(result_dir / "timing_speed_ratios.csv", index=False)

    pd.set_option("display.max_columns", 30)
    print("Comparison summary")
    print(summary.to_string(index=False))
    print()
    print("Non-cluster comparison summary")
    print(summary_noncluster.to_string(index=False))
    print()
    if not time_ratios.empty:
        speed = (
            time_ratios.groupby(["platform", "source"])
            .agg(
                cases=("case", "count"),
                median_local_over_public=("local_over_public", "median"),
                median_public_over_local=("speedup_public_over_local", "median"),
            )
            .reset_index()
        )
        print("Timing ratio summary")
        print(speed.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
