# Contributing

The R package source lives in `R/rd2d`. Run package-development commands from that directory, or pass `R/rd2d` explicitly from the repository root.

## Local Checks

From the repository root:

```sh
R CMD check --no-manual R/rd2d
```

From the package root:

```sh
cd R/rd2d
R CMD check --no-manual .
```

## Repository Layout

- `R/rd2d/`: R package source submitted to CRAN.
- `R/rd2d_illustration.R`: generates an illustration dataset and fitted examples.
- `R/rd2d_plot.R`: builds the illustration figures from saved objects.
- `R/output/`: generated illustration data, fitted objects, and figures. This directory is ignored by git.

Generated package bundles, check directories, and manuals should not be committed. Use GitHub releases for archived release artifacts when needed.
