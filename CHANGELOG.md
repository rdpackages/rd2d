# Changelog

## 2026-05-23

### Stata Package Organization

- Bumped R, Python, and Stata package metadata to version `0.2.0` for
  submission preparation.
- Consolidated compiled Mata helpers into a single `lrd2d.mlib` library for
  cleaner Stata package layout.
- Updated Stata commands to load the Mata library when available and fall back
  to `rd2d_functions.do` during development.

## 2026-05-22

### Cross-Platform Illustration Cleanup

- Replaced generated illustration-output workflows with shared
  `rd2d_data.csv` datasets for R, Python, and Stata.
- Simplified the R, Python, and Stata numerical illustration files to read the
  dataset from the current working directory and print command output only.
- Kept the R simulator as a separate completeness script for regenerating
  `R/rd2d_data.csv`.

### Package Behavior and Documentation

- Moved Gaussian simulation repetitions (`repp`) to the R and Python summary
  calls that use it for confidence bands and aggregate boundary inference.
- Updated summary output labels to distinguish pointwise confidence intervals
  from uniform confidence bands.
- Updated README files, R Rd help files, Stata SMCL help files, and Stata
  package metadata for the dataset-based examples and current options.

### Stata Parity

- Ported the newer bandwidth-selection functionality to Stata, including
  full distance-based automatic bandwidths, fuzzy bandwidth targets, derivative
  options for `rdbw2d`, and known-kink distance bandwidth adjustment.
- Rebuilt the compiled Mata helper files and updated the Stata package manifest.
- Switched Stata examples and illustrations to explicit double-precision data
  generation.

## 2026-05-19

### Stata Package Creation

- Added the initial Stata package under `stata`, including location- and
  distance-based commands, help files and generated PDFs, package metadata,
  compiled Mata helpers, illustration do-files, and scripts for rebuilding
  Stata help PDFs.

### Python Package Creation

- Added the initial Python package under `Python/rd2d`, including package
  metadata, location- and distance-based estimation entry points, bandwidth
  helpers, covariance-backed summary tables, pytest coverage, Python
  illustration and plotting scripts, repository documentation, CI checks, and a
  GitHub Actions workflow for publishing the Python package to PyPI.

## 2026-05-11

### R Package Modernization

- Modernized the R package for the `0.1.0` development line, including fuzzy
  location- and distance-based boundary discontinuity support, covariance-backed
  summary inference, uniform confidence bands, WBATE/LBATE summaries, CER
  bandwidth selectors, distance-kink bandwidth controls, updated return-table
  names, focused numerical tests, and streamlined illustration files.
