from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


CHECKS = [
    {
        "file": "comparison_summary.csv",
        "comparison": "cross_platform_local_separate",
        "max_rel_diff": 1e-8,
        "max_abs_diff": None,
        "label": "local separate cross-platform agreement",
    },
    {
        "file": "comparison_summary_noncluster.csv",
        "comparison": "version_public_vs_local_separate",
        "max_rel_diff": 1e-8,
        "max_abs_diff": 1e-8,
        "label": "non-cluster public-vs-local separate agreement",
    },
]


def run_summarizer(repo_root: Path) -> None:
    script = repo_root / ".benchmarks" / "2026-05-27" / "summarize_benchmarks.py"
    result = subprocess.run([sys.executable, str(script), str(repo_root)], capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        result.check_returncode()


def check_row(row: pd.Series, spec: dict) -> list[str]:
    failures: list[str] = []
    if int(row["n_common"]) == 0:
        failures.append("no common cells")
    max_rel = float(row["max_rel_diff"])
    if max_rel > spec["max_rel_diff"]:
        failures.append(f"max_rel_diff={max_rel:.3g} > {spec['max_rel_diff']:.3g}")
    if spec["max_abs_diff"] is not None:
        max_abs = float(row["max_abs_diff"])
        if max_abs > spec["max_abs_diff"]:
            failures.append(f"max_abs_diff={max_abs:.3g} > {spec['max_abs_diff']:.3g}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Check rd2d benchmark numerical tolerances.")
    parser.add_argument("repo_root", nargs="?", default=".", help="repository root")
    parser.add_argument("--no-refresh", action="store_true", help="read existing summary CSVs without rerunning the summarizer")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    result_dir = repo_root / ".benchmarks" / "2026-05-27" / "results"
    if not args.no_refresh:
        run_summarizer(repo_root)

    all_failures: list[str] = []
    checked = 0
    for spec in CHECKS:
        path = result_dir / spec["file"]
        if not path.exists():
            all_failures.append(f"{spec['label']}: missing {path}")
            continue
        summary = pd.read_csv(path)
        rows = summary[summary["comparison"] == spec["comparison"]]
        if rows.empty:
            all_failures.append(f"{spec['label']}: no rows for {spec['comparison']}")
            continue
        for _, row in rows.iterrows():
            checked += 1
            failures = check_row(row, spec)
            if failures:
                pair = f"{row['lhs']} vs {row['rhs']}"
                where = f"{row['worst_case']} / {row['worst_output']} / {row['worst_row']} / {row['worst_column']}"
                all_failures.append(f"{spec['label']}: {pair}: {'; '.join(failures)}; worst={where}")

    if all_failures:
        print("Benchmark tolerance check FAILED")
        for failure in all_failures:
            print(f"- {failure}")
        return 1

    print(f"Benchmark tolerance check passed ({checked} comparison rows).")
    for spec in CHECKS:
        print(f"- {spec['label']}: rel <= {spec['max_rel_diff']:.0e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
