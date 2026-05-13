make_regression_data <- function(n = 240, seed = 20260507) {
  set.seed(seed)
  x1 <- round(rnorm(n, sd = 1.4), 1)
  x2 <- round(0.35 * x1 + rnorm(n, sd = 1.2), 1)
  z <- as.numeric(x1 + 0.25 * x2 >= 0)
  fuzzy <- as.numeric(runif(n) < ifelse(z == 1, 0.8, 0.2))
  y <- 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * z + rnorm(n, sd = 0.6)
  y.fuzzy <- 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rnorm(n, sd = 0.6)
  b <- matrix(
    c(
      -0.2, -0.3,
       0.0,  0.0,
       0.3,  0.2
    ),
    ncol = 2,
    byrow = TRUE
  )

  list(y = y, y.fuzzy = y.fuzzy, x = cbind(x1, x2), z = z, fuzzy = fuzzy, b = b)
}

expect_rd2d_tables_equal <- function(auto, manual, outputs) {
  for (output in outputs) {
    expect_equal(auto[[output]], manual[[output]], tolerance = 1e-10, ignore_attr = TRUE)
  }
  expect_equal(auto$bw, manual$bw, tolerance = 1e-10, ignore_attr = TRUE)
}

test_that("sharp estimates are unchanged when automatic bandwidths are supplied manually", {
  dat <- make_regression_data()
  scenarios <- expand.grid(
    masspoints = c("check", "adjust", "off"),
    kernel_type = c("prod", "rad"),
    bwselect = c("mserd", "msetwo"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    auto <- suppressWarnings(rd2d(
      dat$y, dat$x, dat$z, dat$b,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwselect = scenario$bwselect,
      bwcheck = 8,
      params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1"),
      repp = 19
    ))
    h <- as.matrix(auto$bw[, c("h01", "h02", "h11", "h12")])
    manual <- suppressWarnings(rd2d(
      dat$y, dat$x, dat$z, dat$b,
      h = h,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwcheck = 8,
      params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1"),
      repp = 19
    ))

    expect_rd2d_tables_equal(auto, manual, c("main", "main.0", "main.1"))
    expect_equal(auto$params.cov, manual$params.cov, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("fuzzy estimates are unchanged when automatic bandwidths are supplied manually", {
  dat <- make_regression_data(seed = 20260508)
  scenarios <- expand.grid(
    masspoints = c("check", "adjust"),
    kernel_type = c("prod", "rad"),
    bwselect = c("mserd", "msetwo"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    auto <- suppressWarnings(rd2d(
      dat$y.fuzzy, dat$x, dat$z, dat$b,
      fuzzy = dat$fuzzy,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwselect = scenario$bwselect,
      bwcheck = 8,
      params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
      params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"),
      repp = 19
    ))
    h <- as.matrix(auto$bw[, c("h01", "h02", "h11", "h12")])
    manual <- suppressWarnings(rd2d(
      dat$y.fuzzy, dat$x, dat$z, dat$b,
      h = h,
      fuzzy = dat$fuzzy,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwcheck = 8,
      params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
      params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"),
      repp = 19
    ))

    expect_rd2d_tables_equal(
      auto, manual,
      c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
    )
    expect_equal(auto$params.cov, manual$params.cov, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("direct rdbw2d calls are deterministic across masspoint and kernel modes", {
  dat <- make_regression_data(seed = 20260509)
  scenarios <- expand.grid(
    masspoints = c("check", "adjust", "off"),
    kernel_type = c("prod", "rad"),
    bwselect = c("mserd", "msetwo"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    first <- suppressWarnings(rdbw2d(
      dat$y, dat$x, dat$z, dat$b,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwselect = scenario$bwselect,
      bwcheck = 8
    ))
    second <- suppressWarnings(rdbw2d(
      dat$y, dat$x, dat$z, dat$b,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwselect = scenario$bwselect,
      bwcheck = 8
    ))

    expect_equal(first$bws, second$bws, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(first$mseconsts, second$mseconsts, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("direct rdbw2d calls remain self-contained under complete-case filtering", {
  dat <- make_regression_data(seed = 20260510)
  y <- dat$y
  x <- dat$x
  z <- dat$z

  y[c(4, 11)] <- NA
  x[c(8, 19), 1] <- NA
  x[27, 2] <- NA

  keep <- complete.cases(y, x[, 1], x[, 2], z)
  raw <- suppressWarnings(rdbw2d(
    y, x, z, dat$b,
    masspoints = "adjust",
    kernel_type = "prod",
    bwcheck = 8
  ))
  clean <- suppressWarnings(rdbw2d(
    y[keep], x[keep, , drop = FALSE], z[keep], dat$b,
    masspoints = "adjust",
    kernel_type = "prod",
    bwcheck = 8
  ))

  expect_equal(raw$bws, clean$bws, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(raw$mseconsts, clean$mseconsts, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(raw$opt$N, sum(keep))
})