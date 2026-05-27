args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript bench_r.R <repo_root> <mode: public|local_separate|local_joint>", call. = FALSE)
}

repo_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
mode <- args[[2]]
bench_dir <- file.path(repo_root, ".benchmarks", "2026-05-27")
input_dir <- file.path(bench_dir, "inputs")
result_dir <- file.path(bench_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!mode %in% c("public", "local_separate", "local_joint")) {
  stop("unknown mode: ", mode, call. = FALSE)
}

if (mode == "public") {
  user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", paste(R.version$major, sub("\\..*", "", R.version$minor), sep = "."))
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
  suppressPackageStartupMessages(library(rd2d))
  version <- as.character(utils::packageVersion("rd2d"))
  fitmethod <- NA_character_
  source_label <- "public"
} else {
  lib <- tempfile("rd2d-r-lib-")
  dir.create(lib)
  status <- utils::install.packages(
    file.path(repo_root, "R", "rd2d"),
    lib = lib, repos = NULL, type = "source", quiet = TRUE,
    INSTALL_opts = "--no-multiarch"
  )
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(rd2d))
  version <- as.character(utils::packageVersion("rd2d"))
  fitmethod <- if (mode == "local_joint") "joint" else "separate"
  source_label <- mode
}

has_arg <- function(fun, arg) arg %in% names(formals(fun))
add_fitmethod <- function(fun, args) {
  if (!is.na(fitmethod) && has_arg(fun, "fitmethod")) args$fitmethod <- fitmethod
  args
}
add_covs <- function(fun, args, covs) {
  if (!is.null(covs) && has_arg(fun, "covs.eff")) args$covs.eff <- covs
  args
}

as_matrix <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.data.frame(x)) return(as.matrix(x))
  if (is.matrix(x)) return(x)
  NULL
}

append_table <- function(rows, case, call, output, tab) {
  M <- as_matrix(tab)
  if (is.null(M)) return(rows)
  cn <- colnames(M)
  if (is.null(cn)) cn <- paste0("V", seq_len(ncol(M)))
  rn <- rownames(M)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(M)))
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
      val <- suppressWarnings(as.numeric(M[i, j]))
      if (is.finite(val) || is.na(val)) {
        rows[[length(rows) + 1]] <- data.frame(
          platform = "R", source = source_label, version = version,
          fitmethod = ifelse(is.na(fitmethod), "implicit_separate", fitmethod),
          case = case, call = call, output = output, row = rn[[i]],
          column = cn[[j]], value = val, stringsAsFactors = FALSE
        )
      }
    }
  }
  rows
}

append_result <- function(rows, case, call, obj) {
  for (nm in c("main", "bw", "bws", "mseconsts", "itt", "fs", "main.0", "main.1", "itt.0", "itt.1", "fs.0", "fs.1")) {
    if (!is.null(obj[[nm]]) && !all(is.na(obj[[nm]]))) {
      rows <- append_table(rows, case, call, nm, obj[[nm]])
    }
  }
  if (!is.null(obj$params.cov)) {
    for (nm in names(obj$params.cov)) {
      rows <- append_table(rows, case, call, paste0("cov.", nm), obj$params.cov[[nm]])
    }
  }
  rows
}

read_input <- function(name) utils::read.csv(file.path(input_dir, name), check.names = FALSE)

location_cases <- function() {
  dat <- read_input("synthetic_location.csv")
  b <- read_input("synthetic_location_eval.csv")
  jasa <- read_input("jasa_location.csv")
  b_jasa <- read_input("jasa_location_eval.csv")
  jss <- read_input("jss_location.csv")
  b_jss <- read_input("jss_location_eval.csv")
  list(
    synth_loc_sharp_fixed = function() do.call(rd2d, add_fitmethod(rd2d, list(
      Y = dat$y, X = dat[, c("x1", "x2")], assignment = dat$assignment, b = b,
      h = 0.95, vce = "hc0", masspoints = "off", bwcheck = NULL,
      params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1")
    ))),
    synth_loc_fuzzy_fixed = function() do.call(rd2d, add_fitmethod(rd2d, list(
      Y = dat$y_fuzzy, X = dat[, c("x1", "x2")], assignment = dat$assignment,
      b = b, h = 0.95, fuzzy = dat$fuzzy, vce = "hc0", masspoints = "off",
      bwcheck = NULL, params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
      params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
    ))),
    synth_loc_cluster_fixed = function() do.call(rd2d, add_fitmethod(rd2d, list(
      Y = dat$y, X = dat[, c("x1", "x2")], assignment = dat$assignment, b = b,
      h = 0.95, cluster = dat$cluster, vce = "hc1", masspoints = "off",
      bwcheck = NULL, params.cov = "main"
    ))),
    synth_loc_bw = function() do.call(rdbw2d, add_fitmethod(rdbw2d, list(
      Y = dat$y, X = dat[, c("x1", "x2")], assignment = dat$assignment, b = b,
      vce = "hc1", masspoints = "off", bwcheck = NULL, stdvars = FALSE
    ))),
    jasa_loc_fuzzy_fixed = function() do.call(rd2d, add_fitmethod(rd2d, list(
      Y = jasa$y, X = jasa[, c("x1", "x2")], assignment = jasa$assignment,
      b = b_jasa, h = 20, fuzzy = jasa$fuzzy, vce = "hc0",
      masspoints = "off", bwcheck = NULL, params.cov = c("main", "itt", "fs")
    ))),
    jss_loc_fuzzy_fixed = function() do.call(rd2d, add_fitmethod(rd2d, list(
      Y = jss$y, X = jss[, c("x1", "x2")], assignment = jss$assignment,
      b = b_jss, h = 20, fuzzy = jss$fuzzy, vce = "hc0",
      masspoints = "off", bwcheck = NULL, params.cov = c("main", "itt", "fs")
    )))
  )
}

distance_cases <- function() {
  dat <- read_input("synthetic_distance.csv")
  b <- read_input("synthetic_distance_eval.csv")
  joe <- read_input("joe_distance.csv")
  b_joe <- read_input("joe_distance_eval.csv")
  jss <- read_input("jss_distance.csv")
  b_jss <- read_input("jss_distance_eval.csv")
  list(
    synth_dist_sharp_fixed = function() do.call(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = dat$y, distance = dat[, "d1", drop = FALSE], b = b, h = 0.5,
      p = 1, q = 2, cbands = TRUE, vce = "hc0", kernel = "tri",
      masspoints = "off", bwcheck = NULL, params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1")
    ))),
    synth_dist_fuzzy_fixed = function() do.call(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = dat$y_fuzzy, distance = dat[, "d1", drop = FALSE], b = b, h = 0.5,
      p = 1, q = 2, fuzzy = dat$fuzzy, cbands = FALSE, vce = "hc0",
      kernel = "tri", masspoints = "off", bwcheck = NULL,
      params.cov = c("main", "itt", "fs")
    ))),
    synth_dist_cluster_fixed = function() do.call(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = dat$y, distance = dat[, "d1", drop = FALSE], b = b, h = 0.5,
      p = 1, q = 2, cluster = dat$cluster, cbands = TRUE, vce = "hc1",
      kernel = "tri", masspoints = "off", bwcheck = NULL, params.cov = "main"
    ))),
    synth_dist_bw = function() do.call(rdbw2d.distance, add_fitmethod(rdbw2d.distance, list(
      Y = dat$y, distance = dat[, "d1", drop = FALSE], b = b, p = 1,
      vce = "hc1", kernel = "tri", masspoints = "off", bwcheck = NULL
    ))),
    joe_dist_fuzzy_fixed = function() do.call(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = joe$y, distance = joe[, c("d1", "d2", "d3")], b = b_joe, h = 20,
      p = 1, q = 2, fuzzy = joe$fuzzy, cbands = FALSE, vce = "hc0",
      kernel = "tri", masspoints = "off", bwcheck = NULL,
      params.cov = c("main", "itt", "fs")
    ))),
    jss_dist_fuzzy_fixed = function() do.call(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = jss$y, distance = jss[, c("d1", "d2", "d3")], b = b_jss, h = 20,
      p = 1, q = 2, fuzzy = jss$fuzzy, cbands = FALSE, vce = "hc0",
      kernel = "tri", masspoints = "off", bwcheck = NULL,
      params.cov = c("main", "itt", "fs")
    )))
  )
}

local_only_cases <- function() {
  if (mode == "public") return(list())
  dat_l <- read_input("synthetic_location.csv")
  b_l <- read_input("synthetic_location_eval.csv")
  dat_d <- read_input("synthetic_distance.csv")
  b_d <- read_input("synthetic_distance_eval.csv")
  list(
    synth_loc_fuzzy_covs_fixed = function() do.call(rd2d, add_covs(rd2d, add_fitmethod(rd2d, list(
      Y = dat_l$y_fuzzy, X = dat_l[, c("x1", "x2")], assignment = dat_l$assignment,
      b = b_l, h = 0.95, fuzzy = dat_l$fuzzy, vce = "hc0", masspoints = "off",
      bwcheck = NULL, params.cov = c("main", "itt", "fs")
    )), dat_l[, c("z1", "z2")])),
    synth_dist_fuzzy_covs_fixed = function() do.call(rd2d.distance, add_covs(rd2d.distance, add_fitmethod(rd2d.distance, list(
      Y = dat_d$y_fuzzy, distance = dat_d[, "d1", drop = FALSE], b = b_d,
      h = 0.5, p = 1, q = 2, fuzzy = dat_d$fuzzy, cbands = FALSE,
      vce = "hc0", kernel = "tri", masspoints = "off", bwcheck = NULL,
      params.cov = c("main", "itt", "fs")
    )), dat_d[, c("z1", "z2")]))
  )
}

cases <- c(location_cases(), distance_cases(), local_only_cases())
result_rows <- list()
timing_rows <- list()

for (case_name in names(cases)) {
  message("R ", source_label, ": ", case_name)
  suppressWarnings(obj <- cases[[case_name]]())
  result_rows <- append_result(result_rows, case_name, "main_call", obj)
  for (rep in seq_len(4)) {
    gc()
    elapsed <- system.time(suppressWarnings(cases[[case_name]]()))[["elapsed"]]
    timing_rows[[length(timing_rows) + 1]] <- data.frame(
      platform = "R", source = source_label, version = version,
      fitmethod = ifelse(is.na(fitmethod), "implicit_separate", fitmethod),
      case = case_name, rep = rep, elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  }
}

results <- if (length(result_rows)) do.call(rbind, result_rows) else data.frame()
timings <- if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame()
utils::write.csv(results, file.path(result_dir, paste0("r_", mode, "_results.csv")), row.names = FALSE)
utils::write.csv(timings, file.path(result_dir, paste0("r_", mode, "_timings.csv")), row.names = FALSE)
utils::write.csv(data.frame(platform = "R", source = source_label, mode = mode, version = version,
                            fitmethod = ifelse(is.na(fitmethod), "implicit_separate", fitmethod),
                            stringsAsFactors = FALSE),
                 file.path(result_dir, paste0("r_", mode, "_metadata.csv")), row.names = FALSE)
