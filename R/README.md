# R Files

This directory contains the R implementation and related illustration material.

- `rd2d/` is the R package source root. It contains `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, and tests.
- `rd2d_data.csv` is the simulated SPP-calibrated cubic fuzzy BD dataset used
  by the illustration scripts, including simulated pre-intervention covariates
  `Z.1` and `Z.2`.
- `rd2d_illustration.R` reads `rd2d_data.csv` and runs covariate-adjusted
  location- and distance-based examples using the joint fit method.
- `rd2d_plot.R` reads `rd2d_data.csv`, recomputes the needed covariate-adjusted
  summaries, and builds inference plots, effect heatmaps, and p-value heatmaps.
- `rd2d_simulate.R` regenerates `rd2d_data.csv` and is included for completeness.

Run the illustration workflow from this directory:

```r
source("rd2d_illustration.R")
source("rd2d_plot.R")
```

The scripts use the current working directory and do not create an `output/`
folder or save fitted objects, tables, or plot files.
