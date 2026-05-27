from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


DEFAULT_STATA_EXE = Path(r"C:\Program Files\Stata16\StataMP-64.exe")
MODES = ("public", "local_separate", "local_joint")
PROFILE_FITMETHODS = ("separate", "joint")


def find_repo_root(path: Path) -> Path:
    path = path.resolve()
    for candidate in (path, *path.parents):
        if (candidate / "R" / "rd2d" / "DESCRIPTION").exists():
            return candidate
    raise SystemExit("Could not find rd2d repository root.")


def command_text(command: list[str | Path]) -> str:
    return " ".join(f'"{str(part)}"' if " " in str(part) else str(part) for part in command)


def run(command: list[str | Path], cwd: Path) -> None:
    print(f"\n==> {command_text(command)}", flush=True)
    subprocess.run([str(part) for part in command], cwd=str(cwd), check=True)


def run_quiet(command: list[str | Path], cwd: Path) -> None:
    print(f"\n==> {command_text(command)}", flush=True)
    result = subprocess.run([str(part) for part in command], cwd=str(cwd), capture_output=True, text=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        result.check_returncode()


def resolve_stata_exe(value: str | None) -> Path:
    if value:
        return Path(value)
    if os.environ.get("STATA_EXE"):
        return Path(os.environ["STATA_EXE"])
    return DEFAULT_STATA_EXE


def run_numerical_benchmarks(
    repo_root: Path,
    bench_dir: Path,
    platforms: set[str],
    modes: tuple[str, ...],
    stata_exe: Path,
    skip_inputs: bool,
) -> None:
    if not skip_inputs:
        run(["Rscript", bench_dir / "make_inputs.R", repo_root], repo_root)

    if "r" in platforms:
        for mode in modes:
            run(["Rscript", bench_dir / "bench_r.R", repo_root, mode], repo_root)

    if "python" in platforms:
        for mode in modes:
            run([sys.executable, bench_dir / "bench_python.py", repo_root, mode], repo_root)

    if "stata" in platforms:
        if not stata_exe.exists():
            raise SystemExit(f"Stata executable not found: {stata_exe}")
        for mode in modes:
            run([stata_exe, "/e", "do", bench_dir / "bench_stata.do", repo_root, mode], repo_root)


def run_profiles(repo_root: Path, bench_dir: Path, platforms: set[str], stata_exe: Path) -> None:
    if "r" in platforms:
        run(["Rscript", bench_dir / "profile_r.R", repo_root], repo_root)

    if "python" in platforms:
        for fitmethod in PROFILE_FITMETHODS:
            run([sys.executable, bench_dir / "profile_python.py", repo_root, fitmethod], repo_root)

    if "stata" in platforms:
        if not stata_exe.exists():
            raise SystemExit(f"Stata executable not found: {stata_exe}")
        for fitmethod in PROFILE_FITMETHODS:
            run([stata_exe, "/e", "do", bench_dir / "profile_stata.do", repo_root, fitmethod], repo_root)


def validate_workspace(bench_dir: Path, require_results: bool) -> None:
    command: list[str | Path] = [sys.executable, bench_dir / "validate_workspace.py", bench_dir]
    if require_results:
        command.append("--require-results")
    run_quiet(command, bench_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run rd2d cross-language numerical benchmarks.")
    parser.add_argument("--repo-root", default=".", help="repository root; auto-detected from this path")
    parser.add_argument("--bench-dir", default=".benchmarks/2026-05-27", help="benchmark workspace relative to repo root")
    parser.add_argument(
        "--platforms",
        nargs="+",
        choices=("r", "python", "stata"),
        default=("r", "python", "stata"),
        help="platform benchmarks to run",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=MODES,
        default=MODES,
        help="benchmark modes to run",
    )
    parser.add_argument("--stata-exe", default=None, help="path to Stata executable; defaults to STATA_EXE or Stata 16 MP")
    parser.add_argument("--skip-inputs", action="store_true", help="reuse existing input CSVs")
    parser.add_argument("--check-only", action="store_true", help="only refresh summaries and run tolerance checks")
    parser.add_argument("--profiles", action="store_true", help="also run larger speed profile workloads")
    parser.add_argument("--no-tolerances", action="store_true", help="skip tolerance checks")
    parser.add_argument("--verbose-summary", action="store_true", help="print full summary tables from summarize_benchmarks.py")
    args = parser.parse_args()

    repo_root = find_repo_root(Path(args.repo_root))
    bench_dir = (repo_root / args.bench_dir).resolve()
    if not bench_dir.exists():
        raise SystemExit(f"Benchmark directory not found: {bench_dir}")

    platforms = set(args.platforms)
    modes = tuple(args.modes)
    stata_exe = resolve_stata_exe(args.stata_exe)

    print(f"Repository: {repo_root}")
    print(f"Benchmark:  {bench_dir}")
    print(f"Platforms:  {', '.join(sorted(platforms))}")
    print(f"Modes:      {', '.join(modes)}")

    validate_workspace(bench_dir, require_results=args.check_only)

    if not args.check_only:
        run_numerical_benchmarks(repo_root, bench_dir, platforms, modes, stata_exe, args.skip_inputs)

    summary_cmd = [sys.executable, bench_dir / "summarize_benchmarks.py", repo_root]
    if args.verbose_summary:
        run(summary_cmd, repo_root)
    else:
        run_quiet(summary_cmd, repo_root)
    if not args.no_tolerances:
        run([sys.executable, bench_dir / "check_tolerances.py", repo_root, "--no-refresh"], repo_root)

    if args.profiles:
        run_profiles(repo_root, bench_dir, platforms, stata_exe)

    run_quiet([sys.executable, bench_dir / "report_status.py"], repo_root)

    print("\nNumerical benchmark workflow completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
