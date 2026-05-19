################################################################################
# rd2d Python Package
# Illustration: simulation and estimation
################################################################################

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
PKG_DIR = SCRIPT_DIR / "rd2d"
OUTPUT_DIR = SCRIPT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

if str(PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PKG_DIR))

from rd2d import rdbw2d, rdbw2d_dist, rd2d, rd2d_dist  # noqa: E402


def design_matrix(dat: pd.DataFrame) -> np.ndarray:
    x1 = dat["x.1"].to_numpy()
    x2 = dat["x.2"].to_numpy()
    return np.column_stack(
        [
            np.ones(len(dat)),
            x1,
            x2,
            x1**2,
            x1 * x2,
            x2**2,
            x1**3,
            x1**2 * x2,
            x1 * x2**2,
            x2**3,
        ]
    )


DGP = {
    "beta_y_0": np.array(
        [
            0.369916579111109,
            0.00430768720228995,
            -0.00245733885625568,
            1.62105590793036e-05,
            7.94581926007163e-06,
            4.24074450908172e-05,
            2.29593705661776e-08,
            1.59000624539961e-07,
            3.45504841426239e-07,
            2.81256828567388e-07,
        ]
    ),
    "beta_y_1": np.array(
        [
            0.736166509744787,
            0.000756347351213138,
            -0.00154115603887117,
            3.49990029700921e-05,
            8.61468650013817e-05,
            -0.000155449166992341,
            -2.82846014355806e-07,
            1.41889707426739e-07,
            -1.6593205900769e-06,
            3.94207017318509e-06,
        ]
    ),
    "beta_fuzzy_0": np.array(
        [
            -26.5660685226902,
            3.23372932760228e-14,
            1.40086045914661e-13,
            -6.50288763217046e-16,
            7.6156958476755e-16,
            -1.30170689826067e-14,
            -2.72635230089794e-18,
            -3.95577524930275e-18,
            -8.48892391384158e-17,
            8.50472091483012e-17,
        ]
    ),
    "beta_fuzzy_1": np.array(
        [
            0.328585902510212,
            0.00259026946365757,
            -0.00265595841237584,
            0.000215463378801299,
            -6.62666277106809e-06,
            -0.000563004965776261,
            -1.56069812328922e-06,
            1.21170156753277e-07,
            2.88676468236169e-06,
            1.28517906890237e-05,
        ]
    ),
    "lambda_0": 0.0442625515378338,
    "lambda_1": 0.581534010672373,
    "sigma_y_0": 0.330524283558143,
    "sigma_y_1": 0.329017540357162,
}


def logistic(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def make_eval_grid(neval: int = 40) -> pd.DataFrame:
    half = int(np.ceil(neval / 2))
    first = pd.DataFrame(
        {
            "x.1": np.zeros(half),
            "x.2": 40 - np.arange(half) * 40 / half,
        }
    )
    second = pd.DataFrame(
        {
            "x.1": np.arange(neval - half) * 56 / half,
            "x.2": np.zeros(neval - half),
        }
    )
    return pd.concat([first, second], ignore_index=True)


def make_signed_distances(X: pd.DataFrame, eval_points: pd.DataFrame, assignment: np.ndarray) -> np.ndarray:
    distance = np.column_stack(
        [
            np.sqrt((X["x.1"] - row["x.1"]) ** 2 + (X["x.2"] - row["x.2"]) ** 2)
            for _, row in eval_points.iterrows()
        ]
    )
    return distance * (2 * assignment[:, None] - 1)


def simulate_spp_cubic(n: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(
        {
            "x.1": 100 * rng.beta(3, 4, size=n) - 25,
            "x.2": 100 * rng.beta(3, 4, size=n) - 25,
        }
    )
    assignment = (X["x.1"].to_numpy() >= 0) & (X["x.2"].to_numpy() >= 0)
    assignment = assignment.astype(float)
    design = design_matrix(X)
    mu_y_0 = design @ DGP["beta_y_0"]
    mu_y_1 = design @ DGP["beta_y_1"]
    mu_fuzzy_0 = logistic(design @ DGP["beta_fuzzy_0"])
    mu_fuzzy_1 = logistic(design @ DGP["beta_fuzzy_1"])
    fuzzy0 = (rng.uniform(size=n) <= mu_fuzzy_0).astype(float)
    fuzzy1 = (rng.uniform(size=n) <= mu_fuzzy_1).astype(float)
    y0 = mu_y_0 + DGP["lambda_0"] * (fuzzy0 - mu_fuzzy_0) + rng.normal(scale=DGP["sigma_y_0"], size=n)
    y1 = mu_y_1 + DGP["lambda_1"] * (fuzzy1 - mu_fuzzy_1) + rng.normal(scale=DGP["sigma_y_1"], size=n)
    return pd.DataFrame(
        {
            "x.1": X["x.1"],
            "x.2": X["x.2"],
            "assignment": assignment,
            "fuzzy": np.where(assignment == 1, fuzzy1, fuzzy0),
            "Y": np.where(assignment == 1, y1, y0),
        }
    )


def main() -> None:
    n = 6000
    seed = 20260508
    repp = 499
    selected = [0, 4, 9, 14, 20, 24, 29, 34, 39]

    dat = simulate_spp_cubic(n=n, seed=seed)
    eval_points = make_eval_grid()
    X = dat[["x.1", "x.2"]]
    distance = make_signed_distances(X, eval_points, dat["assignment"].to_numpy())
    wbate_weights = np.ones(len(eval_points))

    dat.to_csv(OUTPUT_DIR / "rd2d_illustration_data.csv", index=False)
    eval_points.to_csv(OUTPUT_DIR / "rd2d_illustration_eval.csv", index=False)
    pd.DataFrame(distance).to_csv(OUTPUT_DIR / "rd2d_illustration_distances.csv", index=False)

    print("\nLocation-based bandwidth selection with rdbw2d().")
    bw_location = rdbw2d(dat["Y"], X, dat["assignment"], eval_points, masspoints="off")
    print(bw_location.bws.iloc[selected])

    print("\nLocation-based fuzzy estimation with rd2d().")
    fit_location = rd2d(
        dat["Y"],
        X,
        dat["assignment"],
        eval_points,
        fuzzy=dat["fuzzy"],
        masspoints="off",
        params_cov=["main", "itt", "fs"],
        repp=repp,
    )
    summ_location = fit_location.summary(
        output=["main", "itt", "fs"],
        cbands=["main", "itt", "fs"],
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(summ_location.tables["main"].iloc[selected])

    print("\nDistance-based bandwidth selection with rdbw2d_dist().")
    bw_distance = rdbw2d_dist(dat["Y"], distance, b=eval_points, masspoints="off")
    print(bw_distance.bws.iloc[selected])

    print("\nDistance-based fuzzy estimation with rd2d_dist().")
    fit_distance = rd2d_dist(
        dat["Y"],
        distance,
        b=eval_points,
        fuzzy=dat["fuzzy"],
        masspoints="off",
        params_cov=["main", "itt", "fs"],
        repp=repp,
    )
    summ_distance = fit_distance.summary(
        output=["main", "itt", "fs"],
        cbands=["main", "itt", "fs"],
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(summ_distance.tables["main"].iloc[selected])

    with open(OUTPUT_DIR / "rd2d_illustration_results.pkl", "wb") as handle:
        pickle.dump(
            {
                "data": dat,
                "eval": eval_points,
                "distance": distance,
                "bw_location": bw_location,
                "fit_location": fit_location,
                "summary_location": summ_location,
                "bw_distance": bw_distance,
                "fit_distance": fit_distance,
                "summary_distance": summ_distance,
            },
            handle,
        )

    print(f"\nSaved illustration results to {OUTPUT_DIR / 'rd2d_illustration_results.pkl'}")


if __name__ == "__main__":
    main()
