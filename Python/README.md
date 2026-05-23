# Python Implementation

This directory contains the Python implementation and related illustration
material.

- `rd2d/` is the Python package source root. It contains `pyproject.toml`, the
  importable `rd2d/` package, and tests.
- `rd2d_data.csv` is the simulated SPP-calibrated cubic fuzzy BD dataset used
  by the illustration scripts.
- `rd2d_illustration.py` reads `rd2d_data.csv` and runs location- and
  distance-based examples.
- `rd2d_plot.py` reads `rd2d_data.csv`, recomputes the needed summaries, and
  builds inference plots, effect heatmaps, and p-value heatmaps.

The illustration workflow is:

```sh
python rd2d_illustration.py
python rd2d_plot.py
```

Run the scripts from this directory. They use the current working directory and
do not create an `output/` folder or save fitted objects, tables, or plot files.
