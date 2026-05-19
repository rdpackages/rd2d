# Python Implementation

This directory contains the Python implementation and related illustration
material.

- `rd2d/` is the Python package source root. It contains `pyproject.toml`, the
  importable `rd2d/` package, and tests.
- `rd2d_illustration.py` simulates one fuzzy boundary-discontinuity dataset and
  runs location- and distance-based examples.
- `rd2d_plot.py` reads the saved illustration objects and builds inference
  plots, effect heatmaps, and p-value heatmaps.
- `output/` contains generated data, fitted illustration objects, and figures
  from the illustration scripts. This directory is ignored by git.

The illustration workflow is:

```sh
python Python/rd2d_illustration.py
python Python/rd2d_plot.py
```

The scripts can be run from either the repository root or this `Python/`
directory.
