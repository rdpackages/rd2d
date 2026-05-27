from __future__ import annotations

import argparse
from pathlib import Path


SCRIPT_FILES = [
    "README.md",
    "make_inputs.R",
    "bench_r.R",
    "bench_python.py",
    "bench_stata.do",
    "summarize_benchmarks.py",
    "check_tolerances.py",
    "check_numerics.py",
    "report_status.py",
    "profile_r.R",
    "profile_python.py",
    "profile_stata.do",
]

RESULT_FILES = [
    "comparison_summary.csv",
    "comparison_summary_noncluster.csv",
    "timing_summary.csv",
    "timing_speed_ratios.csv",
]

GENERATED_DIRS = [
    "inputs",
    "results",
    "profiles",
    "logs",
]


def missing_paths(base: Path, names: list[str]) -> list[Path]:
    return [base / name for name in names if not (base / name).exists()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the rd2d benchmark workspace layout.")
    parser.add_argument("bench_dir", nargs="?", default=Path(__file__).resolve().parent, type=Path)
    parser.add_argument("--require-results", action="store_true", help="also require generated summary CSVs")
    args = parser.parse_args()

    bench_dir = args.bench_dir.resolve()
    failures: list[str] = []

    if not bench_dir.exists():
        failures.append(f"benchmark directory not found: {bench_dir}")
    else:
        for path in missing_paths(bench_dir, SCRIPT_FILES):
            failures.append(f"missing benchmark source file: {path.name}")

        if args.require_results:
            result_dir = bench_dir / "results"
            if not result_dir.exists():
                failures.append("missing generated results directory")
            else:
                for path in missing_paths(result_dir, RESULT_FILES):
                    failures.append(f"missing generated result file: results/{path.name}")

        for dirname in GENERATED_DIRS:
            path = bench_dir / dirname
            if path.exists() and not path.is_dir():
                failures.append(f"generated path is not a directory: {dirname}")

    if failures:
        print("Benchmark workspace validation FAILED")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"Benchmark workspace validation passed: {bench_dir}")
    if args.require_results:
        print("Generated result summaries are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
