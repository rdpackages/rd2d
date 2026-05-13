# R Files

This directory contains the R implementation and related illustration material.

- `rd2d/` is the R package source root. It contains `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, and tests.
- `rd2d_illustration.R` simulates one SPP-calibrated cubic fuzzy BD dataset and runs location- and distance-based examples.
- `rd2d_plot.R` reads the saved illustration objects and builds inference plots, effect heatmaps, and p-value heatmaps.
- `output/` contains generated data, fitted illustration objects, and figures from the illustration scripts. This directory is ignored by git.

The illustration workflow is:

```r
source("R/rd2d_illustration.R")
source("R/rd2d_plot.R")
```

The scripts can be run from either the repository root or this `R/` directory.
