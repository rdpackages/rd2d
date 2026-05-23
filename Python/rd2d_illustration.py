################################################################################
# RD2D Python Package
# Numerical Illustration
################################################################################

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


if Path("rd2d").exists() and str(Path("rd2d").resolve()) not in sys.path:
    sys.path.insert(0, str(Path("rd2d").resolve()))

from rd2d import rdbw2d, rdbw2d_dist, rd2d, rd2d_dist  # noqa: E402


# Generate boundary evaluation points.
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


# Generate signed distances to each boundary evaluation point.
def make_signed_distances(
    X: pd.DataFrame,
    eval_points: pd.DataFrame,
    assignment: np.ndarray,
) -> np.ndarray:
    distance = np.column_stack(
        [
            np.sqrt(
                (X["x.1"] - row["x.1"]) ** 2
                + (X["x.2"] - row["x.2"]) ** 2
            )
            for _, row in eval_points.iterrows()
        ]
    )
    return distance * (2 * assignment[:, None] - 1)


# Select displayed point rows and aggregate rows.
def display_table(table: pd.DataFrame, selected: list[int]) -> pd.DataFrame:
    point_rows = table.iloc[selected]
    aggregate_rows = table.loc[table.index.isin(["WBATE", "LBATE"])]
    return pd.concat([point_rows, aggregate_rows])


# Run the illustration.
def main() -> None:
    dat = pd.read_csv("rd2d_data.csv")
    X = dat[["x.1", "x.2"]]
    Y = dat["Y"]
    A = dat["assignment"]
    D = dat["fuzzy"]

    neval = int(os.getenv("RD2D_ILLUSTRATION_NEVAL", "40"))
    repp = int(os.getenv("RD2D_ILLUSTRATION_REPP", "499"))
    eval_points = make_eval_grid(neval)
    distance = make_signed_distances(X, eval_points, A.to_numpy())
    wbate_weights = np.ones(len(eval_points))
    selected = [i - 1 for i in [1, 5, 10, 15, 21, 25, 30, 35, 40] if i <= neval]

    # Location-based bandwidth selection.
    bw_location = rdbw2d(Y, X, A, eval_points, masspoints="off")
    print(bw_location.bws.iloc[selected])

    # Location-based fuzzy estimation.
    fit_location = rd2d(
        Y,
        X,
        A,
        eval_points,
        fuzzy=D,
        params_other="itt.0",
        params_cov=["main", "itt", "fs", "itt.0"],
        masspoints="off",
    )

    # Location-based fuzzy main effect.
    summary_location_main = fit_location.summary(
        output="main",
        cbands="main",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_location_main.tables["main"], selected))

    # Location-based reduced-form effect.
    summary_location_itt = fit_location.summary(
        output="itt",
        cbands="itt",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_location_itt.tables["itt"], selected))

    # Location-based first-stage effect.
    summary_location_fs = fit_location.summary(
        output="fs",
        cbands="fs",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_location_fs.tables["fs"], selected))

    # Location-based control-side reduced-form effect.
    summary_location_itt0 = fit_location.summary(
        output="itt.0",
        cbands="itt.0",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_location_itt0.tables["itt.0"], selected))

    # Distance-based bandwidth selection.
    bw_distance = rdbw2d_dist(Y, distance, b=eval_points, masspoints="off")
    print(bw_distance.bws.iloc[selected])

    # Distance-based sharp estimation.
    fit_distance = rd2d_dist(
        Y,
        distance,
        b=eval_points,
        masspoints="off",
        cbands=True,
    )
    summary_distance = fit_distance.summary(
        output="main",
        cbands="main",
        repp=repp,
    )
    print(display_table(summary_distance.tables["main"], selected))

    # Distance-based fuzzy bandwidth selection.
    bw_distance_fuzzy = rdbw2d_dist(
        Y,
        distance,
        b=eval_points,
        fuzzy=D,
        bwparam="itt",
        masspoints="off",
    )
    print(bw_distance_fuzzy.bws.iloc[selected])

    # Distance-based fuzzy estimation.
    fit_distance_fuzzy = rd2d_dist(
        Y,
        distance,
        b=eval_points,
        fuzzy=D,
        bwparam="itt",
        params_cov=["main", "itt", "fs"],
        masspoints="off",
    )

    # Distance-based fuzzy main effect.
    summary_distance_main = fit_distance_fuzzy.summary(
        output="main",
        cbands="main",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_distance_main.tables["main"], selected))

    # Distance-based reduced-form effect.
    summary_distance_itt = fit_distance_fuzzy.summary(
        output="itt",
        cbands="itt",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_distance_itt.tables["itt"], selected))

    # Distance-based first-stage effect.
    summary_distance_fs = fit_distance_fuzzy.summary(
        output="fs",
        cbands="fs",
        repp=repp,
        WBATE=wbate_weights,
        LBATE=True,
    )
    print(display_table(summary_distance_fs.tables["fs"], selected))


if __name__ == "__main__":
    main()
