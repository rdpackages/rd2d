test_that("main user-facing functions are exported", {
  expect_true(is.function(rd2d))
  expect_true(is.function(rd2d.distance))
  expect_true(is.function(rdbw2d))
  expect_true(is.function(rdbw2d.distance))

  expect_true("assignment" %in% names(formals(rd2d)))
  expect_true("assignment" %in% names(formals(rdbw2d)))
  expect_true("fuzzy" %in% names(formals(rd2d)))
  expect_true("fuzzy" %in% names(formals(rdbw2d)))
  expect_true("cluster" %in% names(formals(rd2d)))
  expect_true("cluster" %in% names(formals(rdbw2d)))
  expect_false("t" %in% names(formals(rd2d)))
  expect_false("t" %in% names(formals(rdbw2d)))
  expect_false("W" %in% names(formals(rd2d)))
  expect_false("W" %in% names(formals(rdbw2d)))
  expect_false("C" %in% names(formals(rd2d)))
  expect_false("C" %in% names(formals(rdbw2d)))
})

make_location_data <- function(n = 500, seed = 101) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z <- as.numeric(x1 >= 0)
  y <- 1 + x1 + 0.5 * x2 + z + rnorm(n, sd = 0.7)
  fuzzy <- as.numeric(runif(n) < ifelse(z == 1, 0.85, 0.15))
  y.fuzzy <- 1 + x1 + 0.5 * x2 + 1.2 * fuzzy + rnorm(n, sd = 0.7)
  b <- matrix(c(0, -0.5, 0, 0, 0, 0.5), ncol = 2, byrow = TRUE)
  list(y = y, y.fuzzy = y.fuzzy, x = cbind(x1, x2), z = z, fuzzy = fuzzy, b = b)
}

test_that("sharp rd2d returns pointwise tables without covariance by default", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off", repp = 49
  )

  expect_s3_class(fit, "rd2d")
  expect_equal(
    names(fit$main),
    c(
      "b1", "b2", "estimate.p", "std.err.p", "estimate.q", "std.err.q",
      "t.value", "p.value",
      "ci.lower", "ci.upper", "h01", "h02", "h11", "h12", "N.Co", "N.Tr"
    )
  )
  expect_equal(
    names(fit$bw),
    c("b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr")
  )
  expect_false(any(c("cb.lower", "cb.upper") %in% names(fit$main)))
  expect_equal(length(fit$params.cov), 0)
})

test_that("params.cov stores only requested sharp covariance matrices", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    params.other = "main.0", params.cov = c("main", "main.0"), repp = 49
  )

  expect_true(is.data.frame(fit$main.0))
  expect_equal(names(fit$params.cov), c("main", "main.0"))

  expect_error(
    rd2d(
      dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
      params.cov = "main.1"
    )
  )
})

test_that("fuzzy rd2d returns main, itt, fs, and requested one-sided outputs", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    fuzzy = dat$fuzzy, params.other = c("itt.0", "fs.1"),
    params.cov = c("main", "itt", "fs", "itt.0", "fs.1"), repp = 49
  )

  expect_true(is.data.frame(fit$main))
  expect_true(is.data.frame(fit$itt))
  expect_true(is.data.frame(fit$fs))
  expect_true(is.data.frame(fit$itt.0))
  expect_true(is.data.frame(fit$fs.1))
  expect_equal(
    sort(names(fit$params.cov)),
    sort(c("main", "itt", "fs", "itt.0", "fs.1"))
  )
})

test_that("fuzzy rd2d automatic bandwidths use fuzzy rdbw2d selector", {
  dat <- make_location_data()
  bw <- rdbw2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, fuzzy = dat$fuzzy,
    masspoints = "off", bwcheck = 20, scaleregul = 3
  )
  fit <- rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, fuzzy = dat$fuzzy,
    masspoints = "off", bwcheck = 20, scaleregul = 3, repp = 19
  )

  expect_true(isTRUE(bw$opt$fuzzy))
  expect_equal(bw$opt$bwparam, "main")
  expect_equal(
    fit$bw[, c("h01", "h02", "h11", "h12")],
    bw$bws[, c("h01", "h02", "h11", "h12")],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("location CER bandwidth selectors are rate-corrected MSE selectors", {
  dat <- make_location_data(n = 600, seed = 20260530)
  p <- 1

  mse <- rdbw2d(
    dat$y, dat$x, dat$z, dat$b, p = p, bwselect = "mserd",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE
  )
  cer <- rdbw2d(
    dat$y, dat$x, dat$z, dat$b, p = p, bwselect = "cerrd",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE
  )
  two <- rdbw2d(
    dat$y, dat$x, dat$z, dat$b, p = p, bwselect = "msetwo",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE
  )
  certwo <- rdbw2d(
    dat$y, dat$x, dat$z, dat$b, p = p, bwselect = "certwo",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE
  )

  expect_equal(
    cer$bws[, c("h01", "h02", "h11", "h12")],
    mse$bws[, c("h01", "h02", "h11", "h12")] * rd2d_cer_factor(length(dat$y), p)
  )
  expect_equal(
    certwo$bws[, c("h01", "h02")],
    two$bws[, c("h01", "h02")] * rd2d_cer_factor(sum(dat$z == 0), p)
  )
  expect_equal(
    certwo$bws[, c("h11", "h12")],
    two$bws[, c("h11", "h12")] * rd2d_cer_factor(sum(dat$z == 1), p)
  )
})

test_that("fuzzy bwparam can select reduced-form bandwidth target", {
  dat <- make_location_data()
  bw.itt <- rdbw2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, fuzzy = dat$fuzzy, bwparam = "itt",
    masspoints = "off", bwcheck = 20, scaleregul = 3
  )
  bw.rf <- rdbw2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b,
    masspoints = "off", bwcheck = 20, scaleregul = 3
  )
  fit.itt <- rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, fuzzy = dat$fuzzy, bwparam = "itt",
    masspoints = "off", bwcheck = 20, scaleregul = 3, repp = 19
  )
  sharp.ignored <- rdbw2d(
    dat$y, dat$x, dat$z, dat$b, bwparam = "itt",
    masspoints = "off", bwcheck = 20
  )

  expect_true(isTRUE(bw.itt$opt$fuzzy))
  expect_equal(bw.itt$opt$bwparam, "itt")
  expect_equal(sharp.ignored$opt$bwparam, "main")
  expect_equal(bw.itt$bws, bw.rf$bws, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(
    fit.itt$bw[, c("h01", "h02", "h11", "h12")],
    bw.rf$bws[, c("h01", "h02", "h11", "h12")],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("stored clustered covariance diagonals match pointwise standard errors", {
  dat <- make_location_data(n = 700, seed = 20260511)
  set.seed(20260511)
  cluster <- sample(seq_len(90), length(dat$y), replace = TRUE)

  sharp <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, cluster = cluster, vce = "hc1",
    masspoints = "off", params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1"), repp = 19
  ))
  for (output in c("main", "main.0", "main.1")) {
    expect_equal(
      diag(sharp$params.cov[[output]]),
      sharp[[output]][["std.err.q"]]^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }

  fuzzy <- suppressWarnings(rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.9, cluster = cluster, vce = "hc1",
    masspoints = "off", fuzzy = dat$fuzzy,
    params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
    params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"),
    repp = 19
  ))
  for (output in c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")) {
    expect_equal(
      diag(fuzzy$params.cov[[output]]),
      fuzzy[[output]][["std.err.q"]]^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }
})

test_that("combined clustered fuzzy covariance matches component calculations", {
  dat <- make_location_data(n = 500, seed = 20260512)
  set.seed(20260512)
  cluster <- sample(seq_len(70), length(dat$y), replace = TRUE)

  fit <- suppressWarnings(rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.9, cluster = cluster, vce = "hc1",
    masspoints = "off", fuzzy = dat$fuzzy,
    params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
    params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"),
    repp = 19
  ))

  q <- fit$opt$q
  eval <- data.frame(x.1 = dat$b[, 1], x.2 = dat$b[, 2])
  deriv.q <- matrix(
    0,
    nrow = nrow(eval),
    ncol = factorial(q + 2) / (factorial(q) * 2)
  )
  deriv.q[, 1] <- 1
  hgrid <- cbind(fit$opt$h01, fit$opt$h02)
  hgrid.1 <- cbind(fit$opt$h11, fit$opt$h12)
  dat.Y <- data.frame(
    x.1 = dat$x[, 1], x.2 = dat$x[, 2], y = dat$y.fuzzy,
    d = dat$z, fuzzy = dat$fuzzy
  )
  dat.fs <- dat.Y
  dat.fs$y <- dat.fs$fuzzy

  cov_fuzzy_tables <- getFromNamespace("rdbw2d_cov_fuzzy_tables", "rd2d")
  cov_fuzzy <- getFromNamespace("rdbw2d_cov_fuzzy", "rd2d")
  cov_main <- getFromNamespace("rdbw2d_cov", "rd2d")
  cov_side <- getFromNamespace("rdbw2d_cov_side", "rd2d")

  combined <- cov_fuzzy_tables(
    dat.Y, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
    "hc1", cluster, fit$tau.itt.q, fit$tau.fs.q,
    outputs = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
  )
  component.cov <- list(
    main = cov_fuzzy(
      dat.Y, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster, fit$tau.itt.q, fit$tau.fs.q
    ),
    itt = cov_main(
      dat.Y, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster
    ),
    fs = cov_main(
      dat.fs, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster
    ),
    itt.0 = cov_side(
      dat.Y, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster, side = 0
    ),
    itt.1 = cov_side(
      dat.Y, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster, side = 1
    ),
    fs.0 = cov_side(
      dat.fs, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster, side = 0
    ),
    fs.1 = cov_side(
      dat.fs, eval, deriv.q, c(0, 0), q, hgrid, hgrid.1, "epa", "prod",
      "hc1", cluster, side = 1
    )
  )

  for (output in names(component.cov)) {
    expect_equal(combined[[output]], component.cov[[output]], tolerance = 1e-10)
  }
})

test_that("summary computes uniform bands lazily from stored covariance", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    params.cov = "main", repp = 49
  )

  printed <- capture.output(
    summ <- summary(fit, cbands = "main", WBATE = rep(1, 3), LBATE = TRUE)
  )

  expect_s3_class(summ, "summary.rd2d")
  expect_true(any(grepl("Unif. CB", printed)))
  expect_true(all(c("cb.lower", "cb.upper") %in% names(summ$tables$main)))
  expect_equal(tail(rownames(summ$tables$main), 2), c("WBATE", "LBATE"))
  expect_true(all(is.finite(summ$tables$main[c("WBATE", "LBATE"), "ci.lower"])))
  expect_true(all(is.finite(summ$tables$main[c("WBATE", "LBATE"), "ci.upper"])))
  expect_true(all(is.na(summ$tables$main[c("WBATE", "LBATE"), "cb.lower"])))
  expect_true(all(is.na(summ$tables$main[c("WBATE", "LBATE"), "cb.upper"])))
  expect_false("aggregates" %in% names(summ))
})

test_that("summary aggregates use all evaluation points when display is subset", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    params.cov = "main", repp = 49
  )

  set.seed(20260507)
  capture.output(
    summ.all <- summary(fit, WBATE = rep(1, nrow(fit$main)), LBATE = TRUE)
  )
  set.seed(20260507)
  capture.output(
    summ.subset <- summary(
      fit, WBATE = rep(1, nrow(fit$main)), LBATE = TRUE, subset = 1
    )
  )

  aggregate.cols <- c(
    "estimate.p", "estimate.q", "std.err.q", "t.value", "p.value",
    "ci.lower", "ci.upper"
  )
  expect_equal(
    summ.subset$tables$main[c("WBATE", "LBATE"), aggregate.cols],
    summ.all$tables$main[c("WBATE", "LBATE"), aggregate.cols]
  )
  expect_equal(rownames(summ.subset$tables$main), c("1", "WBATE", "LBATE"))
})

test_that("summary WBATE uses pointwise Gaussian critical value", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    params.cov = "main", repp = 49
  )

  capture.output(
    summ.pointwise <- summary(fit, WBATE = rep(1, nrow(fit$main)))
  )
  capture.output(
    summ.with.bands <- summary(
      fit, cbands = "main", WBATE = rep(1, nrow(fit$main))
    )
  )

  expect_equal(
    summ.with.bands$tables$main["WBATE", c("ci.lower", "ci.upper")],
    summ.pointwise$tables$main["WBATE", c("ci.lower", "ci.upper")]
  )
})

test_that("summary confidence bands use all evaluation points when display is subset", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    params.cov = "main", repp = 49
  )

  set.seed(20260507)
  capture.output(summ.all <- summary(fit, cbands = "main"))
  set.seed(20260507)
  capture.output(summ.subset <- summary(fit, cbands = "main", subset = 1))

  expect_equal(
    summ.subset$tables$main["1", c("cb.lower", "cb.upper")],
    summ.all$tables$main["1", c("cb.lower", "cb.upper")]
  )
})

test_that("summary fuzzy aggregates use all evaluation points", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off",
    fuzzy = dat$fuzzy, params.cov = "itt", repp = 49
  )

  set.seed(20260507)
  capture.output(
    summ.all <- summary(
      fit, output = "itt", WBATE = rep(1, nrow(fit$itt)), LBATE = TRUE
    )
  )
  set.seed(20260507)
  capture.output(
    summ.subset <- summary(
      fit, output = "itt", WBATE = rep(1, nrow(fit$itt)), LBATE = TRUE,
      subset = 1
    )
  )

  aggregate.cols <- c(
    "estimate.p", "estimate.q", "std.err.q", "t.value", "p.value",
    "ci.lower", "ci.upper"
  )
  expect_equal(
    summ.subset$tables$itt[c("WBATE", "LBATE"), aggregate.cols],
    summ.all$tables$itt[c("WBATE", "LBATE"), aggregate.cols]
  )
})

test_that("summary errors clearly when requested inference lacks covariance", {
  dat <- make_location_data()
  fit <- rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.9, masspoints = "off", repp = 49
  )

  expect_error(capture.output(summary(fit, cbands = "main")), "params.cov")
  expect_error(capture.output(summary(fit, WBATE = rep(1, 3))), "params.cov")
  expect_error(capture.output(summary(fit, LBATE = TRUE)), "params.cov")
  expect_error(capture.output(summary(fit, CBuniform = TRUE)), "Unsupported")
})

test_that("second-order derivative estimates include multi-index factorials", {
  grid <- expand.grid(
    x1 = seq(-1, 1, length.out = 13),
    x2 = seq(-1, 1, length.out = 13)
  )
  z <- as.numeric(grid$x1 >= 0)
  y <- 1 + 0.5 * grid$x1 - 0.2 * grid$x2 + z * (1 + 2 * grid$x1^2)
  fit <- rd2d(
    y, as.matrix(grid), z, matrix(c(0, 0), ncol = 2),
    h = 2, deriv = c(2, 0), p = 2, q = 2, kernel = "uni",
    vce = "hc0", masspoints = "off", repp = 19
  )

  expect_equal(fit$main$estimate.p, 4, tolerance = 1e-8)
  expect_equal(fit$main$estimate.q, 4, tolerance = 1e-8)
})

test_that("deriv must be a nonnegative integer multi-index", {
  dat <- make_location_data()

  expect_output(
    expect_error(
      rd2d(
        dat$y, dat$x, dat$z, dat$b, h = 0.9, deriv = c(0.5, 0),
        masspoints = "off", repp = 19
      )
    ),
    "nonnegative integer"
  )
  expect_output(
    expect_error(
      rdbw2d(
        dat$y, dat$x, dat$z, dat$b, p = 2, deriv = c(1, -1),
        masspoints = "off"
      )
    ),
    "nonnegative integer"
  )
})

test_that("covariance regularization protects Gaussian simulation", {
  cval <- getFromNamespace("rd2d_cval", "rd2d")
  bad.cov <- matrix(c(1, 1.02, 1.02, 1), nrow = 2)

  expect_warning(
    crit <- cval(bad.cov, rep = 25, side = "two", alpha = 95),
    "regularized"
  )
  expect_true(is.finite(crit))
})

test_that("right-sided simulated critical values are positive", {
  cval <- getFromNamespace("rd2d_cval", "rd2d")

  set.seed(20260508)
  crit <- cval(diag(40), rep = 200, side = "right", alpha = 95)

  expect_true(is.finite(crit))
  expect_gt(crit, 0)
})

test_that("rdbw2d standardizes by default", {
  dat <- make_location_data()
  bw <- rdbw2d(dat$y, dat$x, dat$z, dat$b, masspoints = "off")

  expect_true(isTRUE(bw$opt$stdvars))
  expect_equal(
    as.matrix(bw$bws[, c("b1", "b2")]),
    dat$b,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("polynomial orders must be nonnegative integers", {
  dat <- make_location_data()

  expect_output(
    expect_error(
      rd2d(
        dat$y, dat$x, dat$z, dat$b, h = 0.9, p = 1.5,
        masspoints = "off", repp = 19
      )
    ),
    "p must be a nonnegative integer"
  )
  expect_output(
    expect_error(
      rd2d(
        dat$y, dat$x, dat$z, dat$b, h = 0.9, q = 1.5,
        masspoints = "off", repp = 19
      )
    ),
    "q must be a nonnegative integer"
  )
  expect_output(
    expect_error(
      rdbw2d(
        dat$y, dat$x, dat$z, dat$b, p = 1.5,
        masspoints = "off"
      )
    ),
    "p must be a nonnegative integer"
  )
})
