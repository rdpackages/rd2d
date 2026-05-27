args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
bench_dir <- file.path(repo_root, ".benchmarks", "2026-05-27")
input_dir <- file.path(bench_dir, "inputs")
if (!dir.exists(input_dir) && !dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create input directory: ", input_dir, call. = FALSE)
}

jasa_dir <- "C:/Users/cattaneo/Dropbox/software/rdpackages-replication/CTY_2026_JASA"
jss_dir <- "C:/Users/cattaneo/Dropbox/research/Cattaneo-Titiunik-Yu_2026_rd2d/overleaf"

write_csv <- function(x, file) {
  utils::write.csv(x, file.path(input_dir, file), row.names = FALSE, na = "")
}

make_distance <- function(x1, x2, eval, assignment) {
  D <- vapply(
    seq_len(nrow(eval)),
    function(j) sqrt((x1 - eval$x1[j])^2 + (x2 - eval$x2[j])^2),
    numeric(length(x1))
  )
  D * matrix(rep(2 * assignment - 1, times = ncol(D)), nrow = nrow(D))
}

set.seed(20260527)
n <- 900
x1 <- round(rnorm(n, sd = 1.4), 2)
x2 <- round(0.35 * x1 + rnorm(n, sd = 1.2), 2)
assignment <- as.integer(x1 + 0.25 * x2 >= 0)
fuzzy <- as.integer(runif(n) < ifelse(assignment == 1, 0.80, 0.20))
cluster <- sample(seq_len(120), n, replace = TRUE)
y <- 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * assignment + rnorm(n, sd = 0.6)
y_fuzzy <- 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rnorm(n, sd = 0.6)
z1 <- 0.5 * x1 - 0.2 * x2 + rnorm(n, sd = 0.4)
z2 <- x1^2 - 0.3 * x2 + rnorm(n, sd = 0.5)
write_csv(data.frame(y, y_fuzzy, x1, x2, assignment, fuzzy, z1, z2, cluster), "synthetic_location.csv")
write_csv(data.frame(x1 = c(-0.2, 0, 0.3), x2 = c(-0.3, 0, 0.2)), "synthetic_location_eval.csv")

set.seed(20260528)
n <- 1400
d <- runif(n, -1.25, 1.25)
assignment <- as.integer(d >= 0)
fuzzy <- as.integer(runif(n) < ifelse(assignment == 1, 0.80 + 0.04 * d, 0.20 + 0.04 * d))
cluster <- sample(seq_len(160), n, replace = TRUE)
y <- 1 + 2 * d + 2.4 * assignment + rnorm(n, sd = 0.3)
y_fuzzy <- 1 + 2 * d + 1.5 * fuzzy + rnorm(n, sd = 0.3)
z1 <- 0.6 * d + rnorm(n, sd = 0.2)
z2 <- d^2 + rnorm(n, sd = 0.2)
write_csv(data.frame(y, y_fuzzy, d1 = d, fuzzy, z1, z2, cluster), "synthetic_distance.csv")
write_csv(data.frame(x1 = 0, x2 = 0), "synthetic_distance_eval.csv")

load_spp <- function(path, max_n = 9000) {
  raw <- utils::read.csv(path)
  out <- data.frame(
    x1 = raw$running_saber11,
    x2 = raw$running_sisben,
    assignment = raw$eligible_spp,
    fuzzy = raw$beneficiary_spp,
    y = raw$spadies_any,
    z1 = raw$icfes_educm1
  )
  out <- out[stats::complete.cases(out), , drop = FALSE]
  if (nrow(out) > max_n) {
    set.seed(20260529)
    out <- out[sort(sample(seq_len(nrow(out)), max_n)), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

spp <- load_spp(file.path(jasa_dir, "spp.csv"))
jasa_eval <- data.frame(x1 = c(0, 0, 28), x2 = c(40, 20, 0))
write_csv(spp, "jasa_location.csv")
write_csv(jasa_eval, "jasa_location_eval.csv")

joe_dat <- spp
x2_scale <- stats::sd(joe_dat$x1) / stats::sd(joe_dat$x2)
joe_dat$x2 <- joe_dat$x2 * x2_scale
joe_eval <- data.frame(x1 = c(0, 0, 28), x2 = c(40, 20, 0) * x2_scale)
joe_dist <- make_distance(joe_dat$x1, joe_dat$x2, joe_eval, joe_dat$assignment)
write_csv(data.frame(
  y = joe_dat$y, y_fuzzy = joe_dat$y, fuzzy = joe_dat$fuzzy, z1 = joe_dat$z1,
  d1 = joe_dist[, 1], d2 = joe_dist[, 2], d3 = joe_dist[, 3]
), "joe_distance.csv")
write_csv(joe_eval, "joe_distance_eval.csv")

jss_file <- file.path(jss_dir, "CTY_2026_rd2d_data.csv")
if (file.exists(jss_file)) {
  jss_raw <- utils::read.csv(jss_file)
  units_raw <- jss_raw[jss_raw$row_type == "unit", , drop = FALSE]
  boundary <- jss_raw[jss_raw$row_type != "unit" & is.finite(jss_raw$b_1) & is.finite(jss_raw$b_2), , drop = FALSE]
  units <- data.frame(
    x1 = units_raw$X_1, x2 = units_raw$X_2, assignment = units_raw$T,
    fuzzy = units_raw$W, y = units_raw$Y
  )
  units <- units[stats::complete.cases(units), , drop = FALSE]
  if (nrow(units) > 9000) {
    set.seed(20260530)
    units <- units[sort(sample(seq_len(nrow(units)), 9000)), , drop = FALSE]
  }
  jss_eval <- unique(data.frame(x1 = boundary$b_1, x2 = boundary$b_2))[c(1, 5, 10), , drop = FALSE]
  rownames(jss_eval) <- NULL
  jss_dist <- make_distance(units$x1, units$x2, jss_eval, units$assignment)
  write_csv(units, "jss_location.csv")
  write_csv(jss_eval, "jss_location_eval.csv")
  write_csv(data.frame(y = units$y, y_fuzzy = units$y, fuzzy = units$fuzzy,
                       d1 = jss_dist[, 1], d2 = jss_dist[, 2], d3 = jss_dist[, 3]),
            "jss_distance.csv")
  write_csv(jss_eval, "jss_distance_eval.csv")
}

files <- list.files(input_dir, full.names = TRUE)
manifest <- data.frame(
  input = basename(files),
  rows = vapply(files, function(f) nrow(utils::read.csv(f)), integer(1))
)
write_csv(manifest, "manifest.csv")
print(manifest)
