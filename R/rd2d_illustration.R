################################################################################
# rd2d R Package
# Illustration: simulation and estimation
################################################################################

rm(list = ls(all = TRUE))

# This script simulates one fuzzy boundary-discontinuity dataset from the
# SPP-calibrated cubic design, runs location- and distance-based examples, and
# saves the fitted objects used by rd2d_plot.R.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/"))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(dirname(frame_files[length(frame_files)]), winslash = "/"))
  }

  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
repo_dir <- if (dir.exists(file.path(script_dir, "rd2d"))) {
  dirname(script_dir)
} else {
  script_dir
}
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_rd2d <- function() {
  local_pkg <- file.path(repo_dir, "R", "rd2d")
  if (dir.exists(local_pkg)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop(
        paste(
          "Package pkgload is required to load the local rd2d source.",
          "Install pkgload, or install rd2d and replace this block with",
          "library(rd2d)."
        ),
        call. = FALSE
      )
    }
    suppressPackageStartupMessages(
      pkgload::load_all(local_pkg, export_all = FALSE, helpers = FALSE, quiet = TRUE)
    )
    return(invisible(TRUE))
  }

  library(rd2d)
}

load_rd2d()

################################## DGP #########################################

design_matrix <- function(dat) {
  x1 <- dat$x.1
  x2 <- dat$x.2
  cbind(
    "(Intercept)" = 1,
    "x.1" = x1,
    "x.2" = x2,
    "I(x.1^2)" = x1^2,
    "I(x.1 * x.2)" = x1 * x2,
    "I(x.2^2)" = x2^2,
    "I(x.1^3)" = x1^3,
    "I(x.1^2 * x.2)" = x1^2 * x2,
    "I(x.1 * x.2^2)" = x1 * x2^2,
    "I(x.2^3)" = x2^3
  )
}

# Cubic SPP calibration from CTY_2026_JASA--simuls.R, design s = 3.
dgp <- list(
  beta_y_0 = c(
    0.369916579111109,
    0.00430768720228995,
    -0.00245733885625568,
    1.62105590793036e-05,
    7.94581926007163e-06,
    4.24074450908172e-05,
    2.29593705661776e-08,
    1.59000624539961e-07,
    3.45504841426239e-07,
    2.81256828567388e-07
  ),
  beta_y_1 = c(
    0.736166509744787,
    0.000756347351213138,
    -0.00154115603887117,
    3.49990029700921e-05,
    8.61468650013817e-05,
    -0.000155449166992341,
    -2.82846014355806e-07,
    1.41889707426739e-07,
    -1.6593205900769e-06,
    3.94207017318509e-06
  ),
  beta_fuzzy_0 = c(
    -26.5660685226902,
    3.23372932760228e-14,
    1.40086045914661e-13,
    -6.50288763217046e-16,
    7.6156958476755e-16,
    -1.30170689826067e-14,
    -2.72635230089794e-18,
    -3.95577524930275e-18,
    -8.48892391384158e-17,
    8.50472091483012e-17
  ),
  beta_fuzzy_1 = c(
    0.328585902510212,
    0.00259026946365757,
    -0.00265595841237584,
    0.000215463378801299,
    -6.62666277106809e-06,
    -0.000563004965776261,
    -1.56069812328922e-06,
    1.21170156753277e-07,
    2.88676468236169e-06,
    1.28517906890237e-05
  ),
  lambda_0 = 0.0442625515378338,
  lambda_1 = 0.581534010672373,
  sigma_y_0 = 0.330524283558143,
  sigma_y_1 = 0.329017540357162
)

make_eval_grid <- function(neval = 40) {
  half <- ceiling(neval / 2)
  rbind(
    data.frame(
      x.1 = rep(0, half),
      x.2 = 40 - (seq_len(half) - 1) * 40 / half
    ),
    data.frame(
      x.1 = (seq_len(neval - half) - 1) * 56 / half,
      x.2 = rep(0, neval - half)
    )
  )
}

make_signed_distances <- function(X, eval, assignment) {
  distance <- sapply(seq_len(nrow(eval)), function(j) {
    sqrt((X$x.1 - eval$x.1[j])^2 + (X$x.2 - eval$x.2[j])^2)
  })
  distance * matrix(rep(2 * assignment - 1, times = ncol(distance)),
                nrow = nrow(distance), ncol = ncol(distance))
}

simulate_spp_cubic <- function(n, seed) {
  set.seed(seed)

  X <- data.frame(
    x.1 = 100 * rbeta(n, 3, 4) - 25,
    x.2 = 100 * rbeta(n, 3, 4) - 25
  )
  assignment <- as.numeric(X$x.1 >= 0 & X$x.2 >= 0)
  R <- design_matrix(X)

  mu_y_0 <- as.numeric(R %*% dgp$beta_y_0)
  mu_y_1 <- as.numeric(R %*% dgp$beta_y_1)
  mu_fuzzy_0 <- plogis(as.numeric(R %*% dgp$beta_fuzzy_0))
  mu_fuzzy_1 <- plogis(as.numeric(R %*% dgp$beta_fuzzy_1))

  fuzzy0 <- as.numeric(runif(n) <= mu_fuzzy_0)
  fuzzy1 <- as.numeric(runif(n) <= mu_fuzzy_1)
  Y0 <- mu_y_0 + dgp$lambda_0 * (fuzzy0 - mu_fuzzy_0) + rnorm(n, sd = dgp$sigma_y_0)
  Y1 <- mu_y_1 + dgp$lambda_1 * (fuzzy1 - mu_fuzzy_1) + rnorm(n, sd = dgp$sigma_y_1)

  data.frame(
    x.1 = X$x.1,
    x.2 = X$x.2,
    assignment = assignment,
    fuzzy = ifelse(assignment == 1, fuzzy1, fuzzy0),
    Y = ifelse(assignment == 1, Y1, Y0)
  )
}

################################ Estimation ####################################

n <- 6000
seed <- 20260508
repp <- 499
selected <- c(1, 5, 10, 15, 21, 25, 30, 35, 40)

dat <- simulate_spp_cubic(n = n, seed = seed)
eval <- make_eval_grid()
X <- dat[, c("x.1", "x.2")]
distance <- make_signed_distances(X, eval, dat$assignment)
wbate_weights <- rep(1, nrow(eval))

write.csv(dat, file.path(output_dir, "rd2d_illustration_data.csv"), row.names = FALSE)
write.csv(eval, file.path(output_dir, "rd2d_illustration_eval.csv"), row.names = FALSE)
write.csv(distance, file.path(output_dir, "rd2d_illustration_distances.csv"), row.names = FALSE)

cat("\nLocation-based bandwidth selection with rdbw2d().\n")
cat("This selects MSE-optimal bandwidths for the reduced-form outcome.\n")
bw_location <- rdbw2d(
  Y = dat$Y,
  X = X,
  assignment = dat$assignment,
  b = eval,
  masspoints = "off",
  stdvars = FALSE
)
summary(bw_location, subset = selected)

cat("\nLocation-based fuzzy fit with rd2d().\n")
cat("The main table is the fuzzy Wald effect; ITT and FS are also returned.\n")
fit_location <- rd2d(
  Y = dat$Y,
  X = X,
  assignment = dat$assignment,
  b = eval,
  fuzzy = dat$fuzzy,
  params.other = "itt.0",
  params.cov = c("main", "itt", "fs", "itt.0"),
  masspoints = "off",
  stdvars = FALSE,
  repp = repp
)

set.seed(seed + 1)
summary(fit_location, output = "main", subset = selected,
        cbands = "main", WBATE = wbate_weights, LBATE = TRUE)

cat("\nReduced-form/ITT summary from the same fuzzy rd2d() fit.\n")
set.seed(seed + 2)
summary(fit_location, output = "itt", subset = selected,
        cbands = "itt", WBATE = wbate_weights, LBATE = TRUE)

cat("\nFirst-stage summary from the same fuzzy rd2d() fit.\n")
set.seed(seed + 3)
summary(fit_location, output = "fs", subset = selected,
        cbands = "fs", WBATE = wbate_weights, LBATE = TRUE)

cat("\nControl-side reduced-form summary requested with params.other = 'itt.0'.\n")
set.seed(seed + 4)
summary(fit_location, output = "itt.0", subset = selected,
        cbands = "itt.0", WBATE = wbate_weights, LBATE = TRUE)

cat("\ndistance-based bandwidth selection with rdbw2d.distance().\n")
cat("Signed distances encode the side of the boundary for each evaluation point.\n")
bw_distance <- rdbw2d.distance(
  Y = dat$Y,
  distance = distance,
  b = eval,
  masspoints = "off"
)
summary(bw_distance, subset = selected)

cat("\nSharp distance-based estimation with rd2d.distance().\n")
cat("This illustrates the distance-score interface for the same simulated data.\n")
set.seed(seed + 5)
fit_distance <- rd2d.distance(
  Y = dat$Y,
  distance = distance,
  b = eval,
  masspoints = "off",
  cbands = TRUE,
  repp = repp
)
set.seed(seed + 6)
summary(fit_distance, subset = selected, cbands = "main")

cat("\nFuzzy distance-based bandwidth selection with rdbw2d.distance().\n")
cat("The bandwidth target is the reduced-form outcome, matching the sharp fit.\n")
bw_distance_fuzzy <- rdbw2d.distance(
  Y = dat$Y,
  distance = distance,
  b = eval,
  fuzzy = dat$fuzzy,
  bwparam = "itt",
  masspoints = "off"
)
summary(bw_distance_fuzzy, subset = selected)

cat("\nFuzzy distance-based estimation with rd2d.distance().\n")
cat("The main table is the fuzzy Wald effect; ITT and FS are also returned.\n")
set.seed(seed + 7)
fit_distance_fuzzy <- rd2d.distance(
  Y = dat$Y,
  distance = distance,
  b = eval,
  fuzzy = dat$fuzzy,
  bwparam = "itt",
  params.cov = c("main", "itt", "fs"),
  masspoints = "off",
  repp = repp
)
set.seed(seed + 8)
summary(fit_distance_fuzzy, output = "main", subset = selected,
        cbands = "main", WBATE = wbate_weights, LBATE = TRUE)

cat("\ndistance-based reduced-form/ITT summary from the same fuzzy rd2d.distance() fit.\n")
set.seed(seed + 9)
summary(fit_distance_fuzzy, output = "itt", subset = selected,
        cbands = "itt", WBATE = wbate_weights, LBATE = TRUE)

cat("\ndistance-based first-stage summary from the same fuzzy rd2d.distance() fit.\n")
set.seed(seed + 10)
summary(fit_distance_fuzzy, output = "fs", subset = selected,
        cbands = "fs", WBATE = wbate_weights, LBATE = TRUE)

quiet_summary <- function(...) {
  set.seed(seed + 20)
  invisible(utils::capture.output(out <- summary(...)))
  out
}

summaries <- list(
  main = quiet_summary(
    fit_location, output = "main", cbands = "main",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  itt = quiet_summary(
    fit_location, output = "itt", cbands = "itt",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  fs = quiet_summary(
    fit_location, output = "fs", cbands = "fs",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  itt.0 = quiet_summary(
    fit_location, output = "itt.0", cbands = "itt.0",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  distance = quiet_summary(
    fit_distance, output = "main", cbands = "main"
  ),
  distance_sharp = quiet_summary(
    fit_distance, output = "main", cbands = "main"
  ),
  distance_fuzzy = quiet_summary(
    fit_distance_fuzzy, output = "main", cbands = "main",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  distance_fuzzy_itt = quiet_summary(
    fit_distance_fuzzy, output = "itt", cbands = "itt",
    WBATE = wbate_weights, LBATE = TRUE
  ),
  distance_fuzzy_fs = quiet_summary(
    fit_distance_fuzzy, output = "fs", cbands = "fs",
    WBATE = wbate_weights, LBATE = TRUE
  )
)

saveRDS(
  list(
    data = dat,
    eval = eval,
    distance = distance,
    selected = selected,
    dgp = dgp,
    meta = data.frame(
      name = c("n", "seed", "dgp.s", "repp"),
      value = c(n, seed, 3, repp)
    ),
    bw_location = bw_location,
    fit_location = fit_location,
    summaries = summaries,
    bw_distance = bw_distance,
    bw_distance_fuzzy = bw_distance_fuzzy,
    fit_distance = fit_distance,
    fit_distance_fuzzy = fit_distance_fuzzy
  ),
  file.path(output_dir, "rd2d_illustration_results.rds")
)

cat(sprintf("\nIllustration objects saved to: %s\n", output_dir))
