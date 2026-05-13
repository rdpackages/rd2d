# Contributing

The R package source lives in `R/rd2d`. Run package-development commands from that directory, or pass `R/rd2d` explicitly from the repository root.

## Local Checks

For routine development checks from the repository root:

```sh
Rscript scripts/check-local.R --dev
```

For the release-style local check matching the CRAN gate:

```sh
Rscript scripts/check-local.R --release
```

The release check regenerates roxygen documentation, runs the testthat suite,
builds the source tarball, and runs `R CMD check --no-manual --as-cran`.
Generated check directories and tarballs are removed after a successful run.

To enable the optional pre-push hook for future local development:

```sh
git config core.hooksPath .githooks
```

The hook runs a non-mutating `R CMD check --no-manual R/rd2d` before each push.
Set `RD2D_SKIP_PRE_PUSH=1` only for emergency pushes where the local check has
already been run separately.

## Repository Layout

- `R/rd2d/`: R package source submitted to CRAN.
- `R/rd2d_illustration.R`: generates an illustration dataset and fitted examples.
- `R/rd2d_plot.R`: builds the illustration figures from saved objects.
- `R/output/`: generated illustration data, fitted objects, and figures. This directory is ignored by git.

Generated package bundles, check directories, and manuals should not be committed. Use GitHub releases for archived release artifacts when needed.
