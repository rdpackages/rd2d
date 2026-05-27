args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript profile_r.R <repo_root>", call. = FALSE)
}

repo_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
bench_dir <- file.path(repo_root, ".benchmarks", "2026-05-27")
profile_dir <- file.path(bench_dir, "profiles")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)

lib <- tempfile("rd2d-profile-r-lib-")
dir.create(lib)
utils::install.packages(
  file.path(repo_root, "R", "rd2d"),
  lib = lib, repos = NULL, type = "source", quiet = TRUE,
  INSTALL_opts = "--no-multiarch"
)
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(rd2d))

make_location <- function(n = 4500, neval = 15, seed = 2026052701) {
  set.seed(seed)
  x1 <- round(rnorm(n, sd = 1.4), 2)
  x2 <- round(0.35 * x1 + rnorm(n, sd = 1.2), 2)
  assignment <- as.integer(x1 + 0.25 * x2 >= 0)
  fuzzy <- as.integer(runif(n) < ifelse(assignment == 1, 0.80, 0.20))
  cluster <- sample(seq_len(450), n, replace = TRUE)
  y <- 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * assignment + rnorm(n, sd = 0.6)
  y_fuzzy <- 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rnorm(n, sd = 0.6)
  theta <- seq(-0.9, 0.9, length.out = neval)
  b <- data.frame(x1 = 0.25 * theta, x2 = theta)
  list(
    y = y, y_fuzzy = y_fuzzy, x = data.frame(x1 = x1, x2 = x2),
    assignment = assignment, fuzzy = fuzzy, cluster = cluster, b = b
  )
}

make_distance <- function(n = 7000, neval = 15, seed = 2026052702) {
  set.seed(seed)
  d0 <- runif(n, -1.25, 1.25)
  assignment <- as.integer(d0 >= 0)
  fuzzy <- as.integer(runif(n) < ifelse(assignment == 1, 0.80 + 0.04 * d0, 0.20 + 0.04 * d0))
  cluster <- sample(seq_len(550), n, replace = TRUE)
  y <- 1 + 2 * d0 + 2.4 * assignment + rnorm(n, sd = 0.3)
  y_fuzzy <- 1 + 2 * d0 + 1.5 * fuzzy + rnorm(n, sd = 0.3)
  shifts <- seq(-0.3, 0.3, length.out = neval)
  distance <- as.data.frame(vapply(shifts, function(s) d0 - s, numeric(n)))
  names(distance) <- paste0("d", seq_len(neval))
  b <- data.frame(x1 = rep(0, neval), x2 = rep(0, neval))
  list(y = y, y_fuzzy = y_fuzzy, distance = distance, fuzzy = fuzzy, cluster = cluster, b = b)
}

location <- make_location()
distance <- make_distance()

cases <- list(
  loc_sharp_cov = function() rd2d(
    location$y, location$x, location$assignment, location$b,
    h = 0.95, vce = "hc0", fitmethod = "separate",
    masspoints = "off", bwcheck = NULL,
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ),
  loc_fuzzy_cluster_cov = function() rd2d(
    location$y_fuzzy, location$x, location$assignment, location$b,
    h = 0.95, fuzzy = location$fuzzy, cluster = location$cluster,
    vce = "hc1", fitmethod = "separate", masspoints = "off",
    bwcheck = NULL, params.cov = c("main", "itt", "fs")
  ),
  loc_bw = function() rdbw2d(
    location$y, location$x, location$assignment, location$b,
    vce = "hc1", fitmethod = "separate", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE
  ),
  dist_sharp_cov = function() rd2d.distance(
    distance$y, distance$distance, b = distance$b, h = 0.5,
    p = 1, q = 2, vce = "hc0", fitmethod = "separate",
    kernel = "tri", masspoints = "off", bwcheck = NULL,
    cbands = FALSE, params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ),
  dist_fuzzy_cluster_cov = function() rd2d.distance(
    distance$y_fuzzy, distance$distance, b = distance$b, h = 0.5,
    p = 1, q = 2, fuzzy = distance$fuzzy, cluster = distance$cluster,
    vce = "hc1", fitmethod = "separate", kernel = "tri",
    masspoints = "off", bwcheck = NULL, cbands = FALSE,
    params.cov = c("main", "itt", "fs")
  ),
  dist_bw = function() rdbw2d.distance(
    distance$y, distance$distance, b = distance$b, p = 1,
    vce = "hc1", fitmethod = "separate", kernel = "tri",
    masspoints = "off", bwcheck = NULL
  )
)

timings <- list()
for (case_name in names(cases)) {
  message("profiling R: ", case_name)
  invisible(suppressWarnings(cases[[case_name]]()))
  for (rep in seq_len(3)) {
    gc()
    elapsed <- system.time(suppressWarnings(cases[[case_name]]()))[["elapsed"]]
    timings[[length(timings) + 1]] <- data.frame(
      case = case_name, rep = rep, elapsed_seconds = elapsed
    )
  }

  profile_file <- file.path(profile_dir, paste0("r_", case_name, ".out"))
  summary_file <- file.path(profile_dir, paste0("r_", case_name, "_summary.txt"))
  Rprof(profile_file, interval = 0.002, memory.profiling = FALSE)
  invisible(suppressWarnings(cases[[case_name]]()))
  Rprof(NULL)
  summary <- utils::summaryRprof(profile_file)
  capture.output(summary, file = summary_file)
}

timing_df <- do.call(rbind, timings)
utils::write.csv(
  timing_df,
  file.path(profile_dir, "r_profile_timings.csv"),
  row.names = FALSE
)
print(stats::aggregate(elapsed_seconds ~ case, timing_df, median))
