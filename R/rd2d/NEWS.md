# rd2d 1.0.0

- Adds `fitmethod = c("joint", "separate")` to `rd2d()`, `rdbw2d()`,
  `rd2d.distance()`, and `rdbw2d.distance()`. The new default
  `fitmethod = "joint"` preserves legacy point estimates while using joint
  variance corrections, including joint HC1 degrees-of-freedom factors and
  cluster-robust covariance calculations that account for clusters observed on
  both treatment sides. `fitmethod = "separate"` retains the previous
  side-specific variance calculation. Bandwidth selection keeps the existing
  plug-in bias estimator and now applies the same joint/separate convention to
  the selector's variance constants.
- Adds `covs.eff` to `rd2d()`, `rdbw2d()`, `rd2d.distance()`, and
  `rdbw2d.distance()` for efficiency adjustment using pre-intervention
  covariates with a common covariate coefficient across treatment sides,
  including sharp, fuzzy, joint, separate, fixed-bandwidth, and
  automatic-bandwidth paths.
- Adds `covs.drop` and `covs.tol` safeguards for rank-deficient
  `covs.eff`. Redundant covariate columns are dropped by default after
  residualizing on the local polynomial basis, with diagnostics stored in
  `opt`.
- Adds S3 methods requested for the JSS resubmission: `plot()` for estimation
  and bandwidth-selection objects, plus `coef()`, `vcov()`, and `confint()` for
  `rd2d()` and `rd2d.distance()` estimation objects.
- Improves exported R argument validation so invalid inputs raise informative
  errors directly rather than printing a message followed by an empty error.
- Adds fuzzy boundary RD support for both location-based and distance-based methods, including fuzzy main effects, ITT, first-stage, and optional one-sided outputs.
- Adds `params.other` and `params.cov` controls to `rd2d()` for optional companion tables and covariance storage.
- Moves the Gaussian simulation count `repp` from `rd2d()` and `rd2d.distance()` to the summary methods that use it for uniform confidence bands and LBATE critical values.
- Reports pointwise CI columns and optional uniform CB columns separately in `summary.rd2d()` and `summary.rd2d.distance()`.
- Speeds up location-based and distance-based fits by reducing repeated low-order basis, kernel, covariance projection, and masspoint-counting work.
- Moves location-based uniform confidence band, WBATE, and LBATE construction to `summary.rd2d()`.
- Adds distance-based WBATE and LBATE construction to `summary.rd2d.distance()`, including fuzzy main, ITT, and first-stage outputs.
- Updates location-based return tables to use lowercase column names such as `estimate.p`, `std.err.q`, `t.value`, and `p.value`.
- Fixes location-based covariance construction for signed cross-evaluation covariance, clustered finite-sample scaling, and right-sided simulated critical values.
- Fixes distance-based covariance construction to use signed local-polynomial residuals, corrected cluster covariance halves, and side-specific bandwidth plug-in matrices.
- Adds CER-optimal bandwidth selectors `cerrd`, `certwo`, `icerrd`, and
  `icertwo` for location-based and distance-based methods.
- Replaces the distance-based `kink`/`rbc` options with `kink.unknown` and
  `kink.position`, including adaptive known-kink bandwidths based on boundary
  point locations.
- Uses smooth-boundary robust bias-corrected inference as the distance-based default, and uses the same stabilized Gram-matrix inversion helper as the location-based methods.
- Uses `q = p` by default for distance-based unknown-kink specifications to
  avoid compounding the unknown-kink bandwidth shrinkage with an additional
  polynomial-order change.
- Aligns distance-based return tables with location-based naming conventions, including `main`, `main.0`, `main.1`, `bw`, `estimate.p`, `std.err.q`, `N.Co`, and `N.Tr`.
- Updates public notation to use `fuzzy` for treatment receipt/status,
  `tau.itt`/`tau.itt.q` for reduced-form outcome estimates,
  `tau.fs`/`tau.fs.q` for first-stage estimates, `assignment` for assignment,
  `cluster` for cluster identifiers, and `distance` for signed distance scores.
- Removes user-facing derivative notation from distance-based help files and printed output.
- Adds focused tests for location-based returns, lazy confidence bands, aggregate inference, covariance regularization, clustered covariance consistency, and numerical preservation.
- Adds focused tests for distance-based return names, fuzzy outputs, heteroskedastic and clustered standard errors, covariance diagonals, aggregate inference, and polynomial-order validation.
- Simplifies repository illustration files to `R/rd2d_illustration.R` and `R/rd2d_plot.R`.

# rd2d 0.0.3

Version published on CRAN on 2025-10-24.
