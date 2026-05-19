################################################################################
# rd2d Python Package
# Illustration: plots
################################################################################

from __future__ import annotations

import pickle
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "output"
PLOT_DIR = OUTPUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)
PKG_DIR = SCRIPT_DIR / "rd2d"

os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".matplotlib-cache"))
if str(PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PKG_DIR))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def label_for_output(output: str) -> str:
    return {
        "main": "Fuzzy",
        "itt": "ITT",
        "fs": "First stage",
        "itt.0": "ITT, control side",
    }.get(output, output)


def point_rows(table: pd.DataFrame) -> pd.DataFrame:
    return table.loc[~table.index.isin(["WBATE", "LBATE"])].copy()


def aggregate_rows(table: pd.DataFrame) -> pd.DataFrame:
    return table.loc[table.index.isin(["WBATE", "LBATE"])].copy()


def make_inference_plot(summary_object, output: str, file: str):
    table = summary_object.tables[output]
    points = point_rows(table)
    aggregates = aggregate_rows(table)
    idx = np.arange(1, len(points) + 1)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    if {"cb.lower", "cb.upper"}.issubset(points.columns):
        ax.fill_between(idx, points["cb.lower"], points["cb.upper"], color="#1b6ca8", alpha=0.16, label="95% CB")
    ax.vlines(idx, points["ci.lower"], points["ci.upper"], color="#1b6ca8", lw=0.7, label="95% CI")
    ax.scatter(idx, points["estimate.p"], s=14, color="#1b6ca8", label=label_for_output(output), zorder=3)
    ax.axhline(0, color="0.45", lw=0.6)
    ax.axvline(21, color="0.85", lw=0.7)
    for name, row in aggregates.iterrows():
        ax.axhline(row["estimate.p"], color="0.25", lw=0.8, ls="--" if name == "WBATE" else "-.", label=name)
    ax.set_xlabel("Boundary evaluation point")
    ax.set_ylabel(label_for_output(output))
    ax.set_xticks([1, 5, 10, 15, 21, 25, 30, 35, 40])
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    fig.tight_layout()
    fig.savefig(PLOT_DIR / file, dpi=200)
    return fig


def make_effect_heatmap(summary_object, output: str, file: str):
    table = point_rows(summary_object.tables[output])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    vals = table["estimate.p"].to_numpy()
    lim = np.nanmax(np.abs(vals))
    sc = ax.scatter(table["b1"], table["b2"], c=vals, cmap="RdBu_r", vmin=-lim, vmax=lim, s=220, marker="s", edgecolor="white", linewidth=0.7)
    for i, (_, row) in enumerate(table.iterrows(), start=1):
        ax.text(row["b1"], row["b2"], str(i), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Score 1")
    ax.set_ylabel("Score 2")
    fig.colorbar(sc, ax=ax, label=label_for_output(output))
    fig.tight_layout()
    fig.savefig(PLOT_DIR / file, dpi=200)
    return fig


def make_pvalue_heatmap(summary_object, output: str, file: str):
    table = point_rows(summary_object.tables[output])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    vals = table["p.value"].to_numpy()
    sc = ax.scatter(table["b1"], table["b2"], c=vals, cmap="Blues_r", vmin=0, vmax=0.1, s=220, marker="s", edgecolor="white", linewidth=0.7)
    for i, (_, row) in enumerate(table.iterrows(), start=1):
        ax.text(row["b1"], row["b2"], str(i), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Score 1")
    ax.set_ylabel("Score 2")
    fig.colorbar(sc, ax=ax, label="p-value")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / file, dpi=200)
    return fig


def main() -> None:
    results_file = OUTPUT_DIR / "rd2d_illustration_results.pkl"
    if not results_file.exists():
        raise FileNotFoundError("Run Python/rd2d_illustration.py before Python/rd2d_plot.py.")

    with open(results_file, "rb") as handle:
        obj = pickle.load(handle)

    make_inference_plot(obj["summary_location"], "main", "location_main_inference.png")
    make_inference_plot(obj["summary_location"], "itt", "location_itt_inference.png")
    make_inference_plot(obj["summary_location"], "fs", "location_fs_inference.png")
    make_effect_heatmap(obj["summary_location"], "main", "location_main_effect_heatmap.png")
    make_pvalue_heatmap(obj["summary_location"], "main", "location_main_pvalue_heatmap.png")

    make_inference_plot(obj["summary_distance"], "main", "distance_main_inference.png")
    make_effect_heatmap(obj["summary_distance"], "main", "distance_main_effect_heatmap.png")
    make_pvalue_heatmap(obj["summary_distance"], "main", "distance_main_pvalue_heatmap.png")

    print(f"Saved plots to {PLOT_DIR}")


if __name__ == "__main__":
    main()
