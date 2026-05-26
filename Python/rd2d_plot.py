################################################################################
# RD2D Python Package
# Plot Illustration
################################################################################

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

if Path("rd2d").exists() and str(Path("rd2d").resolve()) not in sys.path:
    sys.path.insert(0, str(Path("rd2d").resolve()))

from rd2d import rd2d, rd2d_dist  # noqa: E402


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


# Label each estimand for plots.
def label_for_output(output: str) -> str:
    return {
        "main": "Fuzzy",
        "itt": "ITT",
        "fs": "First stage",
        "itt.0": "ITT, control side",
    }.get(output, output)


# Keep boundary rows separate from aggregate rows.
def point_rows(table: pd.DataFrame) -> pd.DataFrame:
    return table.loc[~table.index.isin(["WBATE", "LBATE"])].copy()


# Keep WBATE and LBATE rows for reference lines.
def aggregate_rows(table: pd.DataFrame) -> pd.DataFrame:
    return table.loc[table.index.isin(["WBATE", "LBATE"])].copy()


# Plot point estimates, CIs, and CBs.
def make_inference_plot(summary_object, output: str):
    table = summary_object.tables[output]
    points = point_rows(table)
    aggregates = aggregate_rows(table)
    idx = np.arange(1, len(points) + 1)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    if {"cb.lower", "cb.upper"}.issubset(points.columns):
        ax.fill_between(
            idx,
            points["cb.lower"],
            points["cb.upper"],
            color="#1b6ca8",
            alpha=0.16,
            label="95% CB",
        )
    ax.vlines(
        idx,
        points["ci.lower"],
        points["ci.upper"],
        color="#1b6ca8",
        lw=0.7,
        label="95% CI",
    )
    ax.scatter(
        idx,
        points["estimate.p"],
        s=14,
        color="#1b6ca8",
        label=label_for_output(output),
        zorder=3,
    )
    ax.axhline(0, color="0.45", lw=0.6)
    ax.axvline(21, color="0.85", lw=0.7)
    for name, row in aggregates.iterrows():
        ax.axhline(
            row["estimate.p"],
            color="0.25",
            lw=0.8,
            ls="--" if name == "WBATE" else "-.",
            label=name,
        )
    ax.set_xlabel("Boundary evaluation point")
    ax.set_ylabel(label_for_output(output))
    ax.set_xticks([1, 5, 10, 15, 21, 25, 30, 35, 40])
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    fig.tight_layout()
    return fig


# Plot point estimates over the boundary.
def make_effect_heatmap(summary_object, output: str):
    table = point_rows(summary_object.tables[output])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    vals = table["estimate.p"].to_numpy()
    lim = np.nanmax(np.abs(vals))
    sc = ax.scatter(
        table["b1"],
        table["b2"],
        c=vals,
        cmap="RdBu_r",
        vmin=-lim,
        vmax=lim,
        s=220,
        marker="s",
        edgecolor="white",
        linewidth=0.7,
    )
    for i, (_, row) in enumerate(table.iterrows(), start=1):
        ax.text(row["b1"], row["b2"], str(i), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Score 1")
    ax.set_ylabel("Score 2")
    fig.colorbar(sc, ax=ax, label=label_for_output(output))
    fig.tight_layout()
    return fig


# Plot p-values over the boundary.
def make_pvalue_heatmap(summary_object, output: str):
    table = point_rows(summary_object.tables[output])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    vals = table["p.value"].to_numpy()
    sc = ax.scatter(
        table["b1"],
        table["b2"],
        c=vals,
        cmap="Blues_r",
        vmin=0,
        vmax=0.1,
        s=220,
        marker="s",
        edgecolor="white",
        linewidth=0.7,
    )
    for i, (_, row) in enumerate(table.iterrows(), start=1):
        ax.text(row["b1"], row["b2"], str(i), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Score 1")
    ax.set_ylabel("Score 2")
    fig.colorbar(sc, ax=ax, label="p-value")
    fig.tight_layout()
    return fig


# Plot distance-based point estimates, CIs, and CBs.
def make_distance_inference_plot(
    summary_object,
    output: str = "main",
    y_label: str = "distance-based effect",
    color: str = "#7b3294",
):
    table = point_rows(summary_object.tables[output])
    idx = np.arange(1, len(table) + 1)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    if {"cb.lower", "cb.upper"}.issubset(table.columns):
        ax.fill_between(
            idx,
            table["cb.lower"],
            table["cb.upper"],
            color=color,
            alpha=0.14,
        )
    ax.vlines(idx, table["ci.lower"], table["ci.upper"], color=color, lw=0.7)
    ax.scatter(idx, table["estimate.p"], color=color, s=14, zorder=3)
    ax.axhline(0, color="0.45", lw=0.6)
    ax.axvline(21, color="0.85", lw=0.7)
    ax.set_xlabel("Boundary evaluation point")
    ax.set_ylabel(y_label)
    fig.tight_layout()
    return fig


# Run the plot illustration.
def main(show: bool = True) -> dict[str, object]:
    dat = pd.read_csv("rd2d_data.csv")
    X = dat[["x.1", "x.2"]]
    Y = dat["Y"]
    A = dat["assignment"]
    D = dat["fuzzy"]
    Z = dat[["Z.1", "Z.2"]]

    neval = int(os.getenv("RD2D_ILLUSTRATION_NEVAL", "40"))
    repp = int(os.getenv("RD2D_ILLUSTRATION_REPP", "499"))
    eval_points = make_eval_grid(neval)
    distance = make_signed_distances(X, eval_points, A.to_numpy())
    wbate_weights = np.ones(len(eval_points))

    # Location-based fuzzy estimation.
    fit_location = rd2d(
        Y,
        X,
        A,
        eval_points,
        fuzzy=D,
        params_other="itt.0",
        params_cov=["main", "itt", "fs", "itt.0"],
        covs_eff=Z,
        fitmethod="joint",
        masspoints="off",
    )

    # Distance-based sharp estimation.
    fit_distance = rd2d_dist(
        Y,
        distance,
        b=eval_points,
        covs_eff=Z,
        fitmethod="joint",
        masspoints="off",
        cbands=True,
    )

    # Distance-based fuzzy estimation.
    fit_distance_fuzzy = rd2d_dist(
        Y,
        distance,
        b=eval_points,
        fuzzy=D,
        covs_eff=Z,
        fitmethod="joint",
        bwparam="itt",
        params_cov=["main", "itt", "fs"],
        masspoints="off",
    )

    # Compute summaries used by the plots.
    summaries = {
        "main": fit_location.summary(
            output="main",
            cbands="main",
            repp=repp,
            WBATE=wbate_weights,
            LBATE=True,
        ),
        "itt": fit_location.summary(
            output="itt",
            cbands="itt",
            repp=repp,
            WBATE=wbate_weights,
            LBATE=True,
        ),
        "fs": fit_location.summary(
            output="fs",
            cbands="fs",
            repp=repp,
            WBATE=wbate_weights,
            LBATE=True,
        ),
        "itt.0": fit_location.summary(
            output="itt.0",
            cbands="itt.0",
            repp=repp,
            WBATE=wbate_weights,
            LBATE=True,
        ),
        "distance_sharp": fit_distance.summary(
            output="main",
            cbands="main",
            repp=repp,
        ),
        "distance_fuzzy": fit_distance_fuzzy.summary(
            output="main",
            cbands="main",
            repp=repp,
            WBATE=wbate_weights,
            LBATE=True,
        ),
    }

    # Build plots.
    plots = {}
    for output in ["main", "itt", "fs", "itt.0"]:
        prefix = "fuzzy" if output == "main" else output.replace(".", "")
        plots[f"{prefix}_inference"] = make_inference_plot(summaries[output], output)
        plots[f"{prefix}_heatmap"] = make_effect_heatmap(summaries[output], output)
        plots[f"{prefix}_pvalue_heatmap"] = make_pvalue_heatmap(summaries[output], output)

    plots["distance_sharp_inference"] = make_distance_inference_plot(
        summaries["distance_sharp"],
        output="main",
        y_label="distance-based sharp effect",
        color="#7b3294",
    )
    plots["distance_fuzzy_inference"] = make_distance_inference_plot(
        summaries["distance_fuzzy"],
        output="main",
        y_label="distance-based fuzzy effect",
        color="#a6611a",
    )

    if show:
        plt.show()
    return plots


if __name__ == "__main__":
    main()
