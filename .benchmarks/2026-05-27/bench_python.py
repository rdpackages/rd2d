from __future__ import annotations

import importlib.metadata as metadata
import inspect
import sys
import time
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit("usage: python bench_python.py <repo_root> <mode: public|local_separate|local_joint>")
    repo_root = Path(sys.argv[1]).resolve()
    mode = sys.argv[2]
    if mode not in {"public", "local_separate", "local_joint"}:
        raise SystemExit(f"unknown mode: {mode}")

    if mode != "public":
        sys.path.insert(0, str(repo_root / "Python" / "rd2d"))

    import rd2d  # noqa: WPS433

    bench_dir = repo_root / ".benchmarks" / "2026-05-27"
    input_dir = bench_dir / "inputs"
    result_dir = bench_dir / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    if mode == "public":
        version = metadata.version("rd2d")
        fitmethod = None
        source = "public"
    else:
        with (repo_root / "Python" / "rd2d" / "pyproject.toml").open("rb") as fh:
            version = tomllib.load(fh)["project"]["version"]
        fitmethod = "joint" if mode == "local_joint" else "separate"
        source = mode

    def read_input(name: str) -> pd.DataFrame:
        return pd.read_csv(input_dir / name)

    def with_fitmethod(fun, kwargs: dict) -> dict:
        out = dict(kwargs)
        if fitmethod is not None and "fitmethod" in inspect.signature(fun).parameters:
            out["fitmethod"] = fitmethod
        return out

    def with_covs(fun, kwargs: dict, covs) -> dict:
        out = dict(kwargs)
        if covs is not None and "covs_eff" in inspect.signature(fun).parameters:
            out["covs_eff"] = covs
        return out

    def rows_for_table(case: str, call: str, output: str, table) -> list[dict]:
        if table is None:
            return []
        if isinstance(table, pd.DataFrame):
            frame = table
        else:
            arr = np.asarray(table)
            if arr.size == 0:
                return []
            frame = pd.DataFrame(arr)
        rows = []
        for i, (_, row) in enumerate(frame.iterrows(), start=1):
            row_name = str(frame.index[i - 1])
            for col in frame.columns:
                try:
                    val = float(row[col])
                except (TypeError, ValueError):
                    continue
                rows.append(
                    {
                        "platform": "Python",
                        "source": source,
                        "version": version,
                        "fitmethod": fitmethod or "implicit_separate",
                        "case": case,
                        "call": call,
                        "output": output,
                        "row": row_name,
                        "column": str(col),
                        "value": val,
                    }
                )
        return rows

    def get_attr(obj, name: str):
        if hasattr(obj, name):
            return getattr(obj, name)
        if isinstance(obj, dict):
            return obj.get(name)
        return None

    def append_result(rows: list[dict], case: str, call: str, obj) -> None:
        for name in ["main", "bw", "bws", "mseconsts", "itt", "fs", "main_0", "main_1", "itt_0", "itt_1", "fs_0", "fs_1"]:
            value = get_attr(obj, name)
            if value is not None:
                rows.extend(rows_for_table(case, call, name.replace("_", "."), value))
        params_cov = get_attr(obj, "params_cov")
        if isinstance(params_cov, dict):
            for name, value in params_cov.items():
                rows.extend(rows_for_table(case, call, f"cov.{name}", value))

    def cases() -> dict[str, callable]:
        loc = read_input("synthetic_location.csv")
        loc_b = read_input("synthetic_location_eval.csv")
        dist = read_input("synthetic_distance.csv")
        dist_b = read_input("synthetic_distance_eval.csv")
        jasa = read_input("jasa_location.csv")
        jasa_b = read_input("jasa_location_eval.csv")
        joe = read_input("joe_distance.csv")
        joe_b = read_input("joe_distance_eval.csv")
        jss_l = read_input("jss_location.csv")
        jss_l_b = read_input("jss_location_eval.csv")
        jss_d = read_input("jss_distance.csv")
        jss_d_b = read_input("jss_distance_eval.csv")

        out = {
            "synth_loc_sharp_fixed": lambda: rd2d.rd2d(**with_fitmethod(rd2d.rd2d, dict(
                Y=loc["y"], X=loc[["x1", "x2"]], assignment=loc["assignment"], b=loc_b,
                h=0.95, vce="hc0", masspoints="off", bwcheck=1,
                params_other=["main.0", "main.1"], params_cov=["main", "main.0", "main.1"]
            ))),
            "synth_loc_fuzzy_fixed": lambda: rd2d.rd2d(**with_fitmethod(rd2d.rd2d, dict(
                Y=loc["y_fuzzy"], X=loc[["x1", "x2"]], assignment=loc["assignment"], b=loc_b,
                h=0.95, fuzzy=loc["fuzzy"], vce="hc0", masspoints="off", bwcheck=1,
                params_other=["itt.0", "itt.1", "fs.0", "fs.1"],
                params_cov=["main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"]
            ))),
            "synth_loc_cluster_fixed": lambda: rd2d.rd2d(**with_fitmethod(rd2d.rd2d, dict(
                Y=loc["y"], X=loc[["x1", "x2"]], assignment=loc["assignment"], b=loc_b,
                h=0.95, cluster=loc["cluster"], vce="hc1", masspoints="off", bwcheck=1,
                params_cov="main"
            ))),
            "synth_loc_bw": lambda: rd2d.rdbw2d(**with_fitmethod(rd2d.rdbw2d, dict(
                Y=loc["y"], X=loc[["x1", "x2"]], assignment=loc["assignment"], b=loc_b,
                vce="hc1", masspoints="off", bwcheck=1, stdvars=False
            ))),
            "jasa_loc_fuzzy_fixed": lambda: rd2d.rd2d(**with_fitmethod(rd2d.rd2d, dict(
                Y=jasa["y"], X=jasa[["x1", "x2"]], assignment=jasa["assignment"], b=jasa_b,
                h=20, fuzzy=jasa["fuzzy"], vce="hc0", masspoints="off", bwcheck=1,
                params_cov=["main", "itt", "fs"]
            ))),
            "jss_loc_fuzzy_fixed": lambda: rd2d.rd2d(**with_fitmethod(rd2d.rd2d, dict(
                Y=jss_l["y"], X=jss_l[["x1", "x2"]], assignment=jss_l["assignment"], b=jss_l_b,
                h=20, fuzzy=jss_l["fuzzy"], vce="hc0", masspoints="off", bwcheck=1,
                params_cov=["main", "itt", "fs"]
            ))),
            "synth_dist_sharp_fixed": lambda: rd2d.rd2d_distance(**with_fitmethod(rd2d.rd2d_distance, dict(
                Y=dist["y"], distance=dist[["d1"]], b=dist_b, h=0.5, p=1, q=2,
                cbands=True, vce="hc0", kernel="tri", masspoints="off", bwcheck=1,
                params_other=["main.0", "main.1"], params_cov=["main", "main.0", "main.1"]
            ))),
            "synth_dist_fuzzy_fixed": lambda: rd2d.rd2d_distance(**with_fitmethod(rd2d.rd2d_distance, dict(
                Y=dist["y_fuzzy"], distance=dist[["d1"]], b=dist_b, h=0.5, p=1, q=2,
                fuzzy=dist["fuzzy"], cbands=False, vce="hc0", kernel="tri", masspoints="off",
                bwcheck=1, params_cov=["main", "itt", "fs"]
            ))),
            "synth_dist_cluster_fixed": lambda: rd2d.rd2d_distance(**with_fitmethod(rd2d.rd2d_distance, dict(
                Y=dist["y"], distance=dist[["d1"]], b=dist_b, h=0.5, p=1, q=2,
                cluster=dist["cluster"], cbands=True, vce="hc1", kernel="tri",
                masspoints="off", bwcheck=1, params_cov="main"
            ))),
            "synth_dist_bw": lambda: rd2d.rdbw2d_distance(**with_fitmethod(rd2d.rdbw2d_distance, dict(
                Y=dist["y"], distance=dist[["d1"]], b=dist_b, p=1, vce="hc1",
                kernel="tri", masspoints="off", bwcheck=1
            ))),
            "joe_dist_fuzzy_fixed": lambda: rd2d.rd2d_distance(**with_fitmethod(rd2d.rd2d_distance, dict(
                Y=joe["y"], distance=joe[["d1", "d2", "d3"]], b=joe_b, h=20,
                p=1, q=2, fuzzy=joe["fuzzy"], cbands=False, vce="hc0",
                kernel="tri", masspoints="off", bwcheck=1, params_cov=["main", "itt", "fs"]
            ))),
            "jss_dist_fuzzy_fixed": lambda: rd2d.rd2d_distance(**with_fitmethod(rd2d.rd2d_distance, dict(
                Y=jss_d["y"], distance=jss_d[["d1", "d2", "d3"]], b=jss_d_b, h=20,
                p=1, q=2, fuzzy=jss_d["fuzzy"], cbands=False, vce="hc0",
                kernel="tri", masspoints="off", bwcheck=1, params_cov=["main", "itt", "fs"]
            ))),
        }
        if mode != "public":
            out["synth_loc_fuzzy_covs_fixed"] = lambda: rd2d.rd2d(**with_covs(rd2d.rd2d, with_fitmethod(rd2d.rd2d, dict(
                Y=loc["y_fuzzy"], X=loc[["x1", "x2"]], assignment=loc["assignment"], b=loc_b,
                h=0.95, fuzzy=loc["fuzzy"], vce="hc0", masspoints="off", bwcheck=1,
                params_cov=["main", "itt", "fs"]
            )), loc[["z1", "z2"]]))
            out["synth_dist_fuzzy_covs_fixed"] = lambda: rd2d.rd2d_distance(**with_covs(rd2d.rd2d_distance, with_fitmethod(rd2d.rd2d_distance, dict(
                Y=dist["y_fuzzy"], distance=dist[["d1"]], b=dist_b, h=0.5,
                p=1, q=2, fuzzy=dist["fuzzy"], cbands=False, vce="hc0",
                kernel="tri", masspoints="off", bwcheck=1, params_cov=["main", "itt", "fs"]
            )), dist[["z1", "z2"]]))
        return out

    result_rows: list[dict] = []
    timing_rows: list[dict] = []
    for case_name, fn in cases().items():
        print(f"Python {source}: {case_name}", flush=True)
        obj = fn()
        append_result(result_rows, case_name, "main_call", obj)
        for rep in range(1, 5):
            start = time.perf_counter()
            fn()
            elapsed = time.perf_counter() - start
            timing_rows.append(
                {
                    "platform": "Python",
                    "source": source,
                    "version": version,
                    "fitmethod": fitmethod or "implicit_separate",
                    "case": case_name,
                    "rep": rep,
                    "elapsed_seconds": elapsed,
                }
            )

    pd.DataFrame(result_rows).to_csv(result_dir / f"python_{mode}_results.csv", index=False)
    pd.DataFrame(timing_rows).to_csv(result_dir / f"python_{mode}_timings.csv", index=False)
    pd.DataFrame([{
        "platform": "Python",
        "source": source,
        "mode": mode,
        "version": version,
        "fitmethod": fitmethod or "implicit_separate",
    }]).to_csv(result_dir / f"python_{mode}_metadata.csv", index=False)


if __name__ == "__main__":
    main()
