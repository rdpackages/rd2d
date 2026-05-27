from __future__ import annotations

import cProfile
import pstats
import sys
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd


def make_location(n: int = 4500, neval: int = 15, seed: int = 2026052701) -> dict:
    rng = np.random.default_rng(seed)
    x1 = np.round(rng.normal(scale=1.4, size=n), 2)
    x2 = np.round(0.35 * x1 + rng.normal(scale=1.2, size=n), 2)
    assignment = (x1 + 0.25 * x2 >= 0).astype(float)
    fuzzy = (rng.uniform(size=n) < np.where(assignment == 1, 0.80, 0.20)).astype(float)
    cluster = rng.integers(1, 451, size=n)
    y = 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * assignment + rng.normal(scale=0.6, size=n)
    y_fuzzy = 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rng.normal(scale=0.6, size=n)
    theta = np.linspace(-0.9, 0.9, neval)
    b = pd.DataFrame({"x1": 0.25 * theta, "x2": theta})
    x = pd.DataFrame({"x1": x1, "x2": x2})
    return {
        "y": y,
        "y_fuzzy": y_fuzzy,
        "x": x,
        "assignment": assignment,
        "fuzzy": fuzzy,
        "cluster": cluster,
        "b": b,
    }


def make_distance(n: int = 7000, neval: int = 15, seed: int = 2026052702) -> dict:
    rng = np.random.default_rng(seed)
    d0 = rng.uniform(-1.25, 1.25, size=n)
    assignment = (d0 >= 0).astype(float)
    fuzzy = (rng.uniform(size=n) < np.where(assignment == 1, 0.80 + 0.04 * d0, 0.20 + 0.04 * d0)).astype(float)
    cluster = rng.integers(1, 551, size=n)
    y = 1 + 2 * d0 + 2.4 * assignment + rng.normal(scale=0.3, size=n)
    y_fuzzy = 1 + 2 * d0 + 1.5 * fuzzy + rng.normal(scale=0.3, size=n)
    shifts = np.linspace(-0.3, 0.3, neval)
    distance = pd.DataFrame({f"d{j + 1}": d0 - shift for j, shift in enumerate(shifts)})
    b = pd.DataFrame({"x1": np.zeros(neval), "x2": np.zeros(neval)})
    return {
        "y": y,
        "y_fuzzy": y_fuzzy,
        "distance": distance,
        "fuzzy": fuzzy,
        "cluster": cluster,
        "b": b,
    }


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: python profile_python.py <repo_root> [separate|joint]")
    repo_root = Path(sys.argv[1]).resolve()
    fitmethod = sys.argv[2] if len(sys.argv) >= 3 else "separate"
    if fitmethod not in {"separate", "joint"}:
        raise SystemExit("fitmethod must be separate or joint")

    sys.path.insert(0, str(repo_root / "Python" / "rd2d"))
    import rd2d  # noqa: WPS433

    profile_dir = repo_root / ".benchmarks" / "2026-05-27" / "profiles"
    profile_dir.mkdir(parents=True, exist_ok=True)

    location = make_location()
    distance = make_distance()

    cases: dict[str, Callable[[], object]] = {
        "loc_sharp_cov": lambda: rd2d.rd2d(
            location["y"],
            location["x"],
            location["assignment"],
            location["b"],
            h=0.95,
            vce="hc0",
            fitmethod=fitmethod,
            masspoints="off",
            bwcheck=None,
            params_other=["main.0", "main.1"],
            params_cov=["main", "main.0", "main.1"],
        ),
        "loc_fuzzy_cluster_cov": lambda: rd2d.rd2d(
            location["y_fuzzy"],
            location["x"],
            location["assignment"],
            location["b"],
            h=0.95,
            fuzzy=location["fuzzy"],
            cluster=location["cluster"],
            vce="hc1",
            fitmethod=fitmethod,
            masspoints="off",
            bwcheck=None,
            params_cov=["main", "itt", "fs"],
        ),
        "loc_bw": lambda: rd2d.rdbw2d(
            location["y"],
            location["x"],
            location["assignment"],
            location["b"],
            vce="hc1",
            fitmethod=fitmethod,
            masspoints="off",
            bwcheck=None,
            stdvars=False,
        ),
        "dist_sharp_cov": lambda: rd2d.rd2d_distance(
            distance["y"],
            distance["distance"],
            b=distance["b"],
            h=0.5,
            p=1,
            q=2,
            vce="hc0",
            fitmethod=fitmethod,
            kernel="tri",
            masspoints="off",
            bwcheck=None,
            cbands=False,
            params_other=["main.0", "main.1"],
            params_cov=["main", "main.0", "main.1"],
        ),
        "dist_fuzzy_cluster_cov": lambda: rd2d.rd2d_distance(
            distance["y_fuzzy"],
            distance["distance"],
            b=distance["b"],
            h=0.5,
            p=1,
            q=2,
            fuzzy=distance["fuzzy"],
            cluster=distance["cluster"],
            vce="hc1",
            fitmethod=fitmethod,
            kernel="tri",
            masspoints="off",
            bwcheck=None,
            cbands=False,
            params_cov=["main", "itt", "fs"],
        ),
        "dist_bw": lambda: rd2d.rdbw2d_distance(
            distance["y"],
            distance["distance"],
            b=distance["b"],
            p=1,
            vce="hc1",
            fitmethod=fitmethod,
            kernel="tri",
            masspoints="off",
            bwcheck=None,
        ),
    }

    timing_rows: list[dict] = []
    for case_name, fn in cases.items():
        print(f"profiling Python {fitmethod}: {case_name}", flush=True)
        fn()
        for rep in range(1, 4):
            start = time.perf_counter()
            fn()
            elapsed = time.perf_counter() - start
            timing_rows.append({"case": case_name, "fitmethod": fitmethod, "rep": rep, "elapsed_seconds": elapsed})

        profile_path = profile_dir / f"python_{fitmethod}_{case_name}.prof"
        summary_path = profile_dir / f"python_{fitmethod}_{case_name}_summary.txt"
        profiler = cProfile.Profile()
        profiler.enable()
        fn()
        profiler.disable()
        profiler.dump_stats(profile_path)
        with summary_path.open("w", encoding="utf-8") as fh:
            stats = pstats.Stats(profiler, stream=fh).strip_dirs().sort_stats("cumtime")
            stats.print_stats(40)

    timing_df = pd.DataFrame(timing_rows)
    timing_df.to_csv(profile_dir / f"python_{fitmethod}_profile_timings.csv", index=False)
    print(timing_df.groupby("case", as_index=False)["elapsed_seconds"].median())


if __name__ == "__main__":
    main()
