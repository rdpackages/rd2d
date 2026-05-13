make_distance_data <- function(n = 500, seed = 20260509) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  d <- as.numeric(x >= 0)
  y <- 1 + 2 * x + 3 * d + rnorm(n, sd = 0.25)
  fuzzy <- as.numeric(
    runif(n) < ifelse(d == 1, 0.8 + 0.05 * x, 0.2 + 0.05 * x)
  )
  y.fuzzy <- 1 + 2 * x + 1.5 * fuzzy + rnorm(n, sd = 0.25)
  list(y = y, y.fuzzy = y.fuzzy, distance = matrix(x, ncol = 1), x = x, fuzzy = fuzzy)
}

manual_dist_se <- function(y, x, h, p, side, kernel = "tri") {
  ind.side <- (x >= 0) == side
  u <- abs(x[ind.side])
  y <- y[ind.side]
  w <- kernel_weight(u / h, kernel) / h^2
  ind <- w > 0
  u <- u[ind]
  y <- y[ind]
  w <- w[ind]
  R <- sapply(0:p, function(j) (u / h)^j)
  G.inv <- MASS::ginv(crossprod(sqrt(w) * R), 1e-20)
  beta <- G.inv %*% crossprod(R, w * y)
  res <- as.numeric(y - R %*% beta)
  meat <- crossprod(res * (w * R)) * h^2
  cov <- G.inv %*% meat %*% G.inv
  sqrt(cov[1, 1] / h^2)
}

test_that("distance bandwidth polynomial helper matches lm", {
  set.seed(20260513)
  distance <- runif(140, 0, 1.5)
  y <- 1 + 0.5 * distance - 0.2 * distance^2 + 0.1 * distance^3 +
    rnorm(length(distance), sd = 0.2)

  fast <- rd2d_distance_poly_lm(y, distance, degree = 3)
  slow <- lm(y ~ distance + I(distance^2) + I(distance^3))

  expect_equal(fast$coef, unname(coef(slow)), tolerance = 1e-12)
  expect_equal(fast$se, unname(summary(slow)$coefficients[, 2]), tolerance = 1e-12)
  expect_equal(fast$predict(distance), unname(predict(slow)), tolerance = 1e-12)
})

test_that("distance public API uses distance argument", {
  expect_true("distance" %in% names(formals(rd2d.distance)))
  expect_true("distance" %in% names(formals(rdbw2d.distance)))
  expect_true("kink.unknown" %in% names(formals(rd2d.distance)))
  expect_true("kink.unknown" %in% names(formals(rdbw2d.distance)))
  expect_true("kink.position" %in% names(formals(rd2d.distance)))
  expect_true("kink.position" %in% names(formals(rdbw2d.distance)))
  expect_false("D" %in% names(formals(rd2d.distance)))
  expect_false("D" %in% names(formals(rdbw2d.distance)))
  expect_false("kink" %in% names(formals(rd2d.distance)))
  expect_false("kink" %in% names(formals(rdbw2d.distance)))
  expect_false("kink.known" %in% names(formals(rd2d.distance)))
  expect_false("kink.known" %in% names(formals(rdbw2d.distance)))
  expect_false("rbc" %in% names(formals(rd2d.distance)))

  dat <- make_distance_data()
  b <- data.frame(x.1 = 0, x.2 = 0)
  fit <- rd2d.distance(
    Y = dat$y, distance = dat$distance, h = 0.45, b = b,
    p = 1, q = 1, cbands = FALSE, masspoints = "off",
    bwcheck = NULL
  )
  bws <- rdbw2d.distance(
    Y = dat$y, distance = dat$distance, b = b,
    masspoints = "off", bwcheck = NULL
  )

  expect_s3_class(fit, "rd2d.distance")
  expect_s3_class(bws, "rdbw2d.distance")
})

test_that("rd2d.distance returns location-style names", {
  dat <- make_distance_data()
  fit <- rd2d.distance(
    dat$y, dat$distance, h = 0.45, b = data.frame(x.1 = 0, x.2 = 0),
    p = 1, q = 1, cbands = FALSE, masspoints = "off",
    bwcheck = NULL
  )

  expect_equal(
    names(fit$main),
    c(
      "b1", "b2", "estimate.p", "std.err.p", "estimate.q", "std.err.q",
      "t.value", "p.value", "ci.lower", "ci.upper", "h0", "h1",
      "h0.rbc", "h1.rbc", "N.Co", "N.Tr"
    )
  )
  expect_equal(names(fit$bw), c("b1", "b2", "h0", "h1", "N.Co", "N.Tr"))
  expect_equal(
    names(fit),
    c(
      "main", "bw", "main.0", "main.1", "opt", "tau.hat", "tau.hat.q",
      "se.hat", "se.hat.q", "params.cov", "cb", "pvalues", "tvalues",
      "tau.itt", "tau.itt.q", "tau.fs", "tau.fs.q", "rdmodel", "call"
    )
  )
  expect_true(is.data.frame(fit$main.0))
  expect_true(is.data.frame(fit$main.1))
  expect_equal(length(fit$params.cov), 0)
  expect_false(any(c(
    "results", "results.A0", "results.A1", "main.A0", "main.A1",
    "cov.q", "cov.us", "zvalues"
  ) %in% names(fit)))
  expect_false(any(c("cb.lower", "cb.upper") %in% names(fit$main)))
})

test_that("rdbw2d.distance automatic bandwidths reproduce manual rd2d.distance", {
  dat <- make_distance_data(n = 600, seed = 20260514)
  b <- data.frame(x.1 = 0, x.2 = 0)

  bws <- rdbw2d.distance(
    dat$y, dat$distance, b = b, p = 1, kernel = "tri", vce = "hc1",
    masspoints = "off", bwcheck = NULL
  )
  auto <- rd2d.distance(
    dat$y, dat$distance, b = b, p = 1, kernel = "tri", vce = "hc1",
    masspoints = "off", bwcheck = NULL, cbands = TRUE, repp = 49
  )
  manual <- rd2d.distance(
    dat$y, dat$distance, h = as.matrix(bws$bws[, c("h0", "h1")]),
    b = b, p = 1, kernel = "tri", vce = "hc1",
    masspoints = "off", bwcheck = NULL, cbands = TRUE, repp = 49
  )

  expect_equal(auto$bw, manual$bw, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(auto$main, manual$main, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(auto$params.cov, manual$params.cov, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("rd2d.distance smooth-boundary default uses q equal to p plus one", {
  dat <- make_distance_data(n = 500, seed = 20260512)
  fit <- rd2d.distance(
    dat$y, dat$distance, h = 0.65, b = data.frame(x.1 = 0, x.2 = 0),
    p = 2, cbands = TRUE, masspoints = "off", bwcheck = NULL, repp = 49
  )

  expect_equal(fit$opt$kink.unknown, c(FALSE, FALSE))
  expect_equal(fit$opt$q, fit$opt$p + 1)
  expect_output(summary(fit, cbands = "main"), "Unif. CB")

  capture.output(summ <- summary(fit, cbands = "main"))
  expect_s3_class(summ, "summary.rd2d.distance")
  expect_equal(summ$outputs, "main")
  expect_equal(names(summ$tables), "main")
  expect_equal(names(summ$cbands), "main")
  expect_equal(
    names(summ$tables$main),
    c(names(fit$main), "cb.lower", "cb.upper")
  )
  expect_equal(
    summ$tables$main[, c("cb.lower", "cb.upper"), drop = FALSE],
    summ$cbands$main,
    ignore_attr = TRUE
  )

  capture.output(bw.summ <- summary(fit, output = "bw"))
  expect_s3_class(bw.summ, "summary.rd2d.distance")
  expect_equal(bw.summ$outputs, "bw")
  expect_equal(bw.summ$tables$bw, fit$bw, ignore_attr = TRUE)

  expect_error(summary(fit, extra = TRUE), "Unsupported summary.rd2d.distance")
  expect_error(summary(fit, CBuniform = TRUE), "Unsupported summary.rd2d.distance")
  expect_error(
    summary(fit, cbands = "itt"),
    "cbands contains unavailable output"
  )

  no.cov <- rd2d.distance(
    dat$y, dat$distance, h = 0.65, b = data.frame(x.1 = 0, x.2 = 0),
    p = 2, cbands = FALSE, masspoints = "off", bwcheck = NULL
  )
  expect_error(
    capture.output(summary(no.cov, cbands = "main")),
    "Uniform confidence bands require a stored covariance matrix"
  )
})

test_that("rd2d.distance unknown kink defaults to q equal to p", {
  dat <- make_distance_data(n = 500, seed = 20260512)
  b <- data.frame(x.1 = 0, x.2 = 0)

  kink.default <- rd2d.distance(
    dat$y, dat$distance, b = b, p = 1, kink.unknown = c(TRUE, FALSE),
    cbands = FALSE, masspoints = "off", bwcheck = NULL
  )
  manual.pq <- rd2d.distance(
    dat$y, dat$distance, h = 0.65, b = b, p = 1, q = 1,
    cbands = FALSE, masspoints = "off", bwcheck = NULL
  )
  manual.q2 <- rd2d.distance(
    dat$y, dat$distance, b = b, p = 1, q = 2,
    kink.unknown = c(TRUE, FALSE), cbands = FALSE,
    masspoints = "off", bwcheck = NULL
  )

  expect_equal(kink.default$opt$kink.unknown, c(TRUE, FALSE))
  expect_equal(kink.default$opt$q, kink.default$opt$p)
  expect_equal(kink.default$main$estimate.q, kink.default$main$estimate.p)
  expect_equal(kink.default$main$std.err.q, kink.default$main$std.err.p)
  expect_equal(manual.pq$opt$q, manual.pq$opt$p)
  expect_equal(manual.pq$main$estimate.q, manual.pq$main$estimate.p)
  expect_equal(manual.pq$main$std.err.q, manual.pq$main$std.err.p)
  expect_equal(manual.q2$opt$q, 2)
})

test_that("distance methods default to smooth-boundary bandwidths", {
  dat <- make_distance_data(n = 500, seed = 20260516)
  b <- data.frame(x.1 = 0, x.2 = 0)

  bws <- rdbw2d.distance(
    dat$y, dat$distance, b = b, masspoints = "off", bwcheck = NULL
  )
  fit <- rd2d.distance(
    dat$y, dat$distance, b = b, masspoints = "off", bwcheck = NULL,
    cbands = FALSE
  )

  expect_equal(bws$opt$kink.unknown, c(FALSE, FALSE))
  expect_equal(fit$opt$kink.unknown, c(FALSE, FALSE))
  expect_equal(fit$opt$q, fit$opt$p + 1)
  expect_equal(fit$main$h0.rbc, fit$main$h0, tolerance = 1e-12)
  expect_equal(fit$main$h1.rbc, fit$main$h1, tolerance = 1e-12)
})

test_that("distance methods support optional unknown-kink bandwidths", {
  dat <- make_distance_data(n = 500, seed = 20260517)
  b <- data.frame(x.1 = 0, x.2 = 0)

  bws <- rdbw2d.distance(
    dat$y, dat$distance, b = b, kink.unknown = c(TRUE, FALSE),
    masspoints = "off", bwcheck = NULL
  )
  fit <- rd2d.distance(
    dat$y, dat$distance, b = b, kink.unknown = c(TRUE, TRUE),
    masspoints = "off", bwcheck = NULL, cbands = FALSE
  )

  expect_equal(bws$opt$kink.unknown, c(TRUE, FALSE))
  expect_equal(fit$opt$kink.unknown, c(TRUE, TRUE))
  expect_equal(fit$opt$q, fit$opt$p)
  expect_true(all(fit$main$h0.rbc <= fit$main$h0 + sqrt(.Machine$double.eps)))
  expect_true(all(fit$main$h1.rbc <= fit$main$h1 + sqrt(.Machine$double.eps)))

  scalar <- rdbw2d.distance(
    dat$y, dat$distance, b = b, kink.unknown = TRUE,
    masspoints = "off", bwcheck = NULL
  )
  expect_equal(scalar$opt$kink.unknown, c(TRUE, TRUE))
})

test_that("distance methods support known-kink adaptive bandwidths", {
  dat <- make_distance_data(n = 500, seed = 20260528)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)
  kink.position <- c(FALSE, TRUE, FALSE)

  smooth <- rdbw2d.distance(
    dat$y, distance, b = b, p = 1, masspoints = "off", bwcheck = NULL
  )
  adaptive <- rdbw2d.distance(
    dat$y, distance, b = b, p = 1, kink.position = kink.position,
    masspoints = "off", bwcheck = NULL
  )
  adaptive.integer <- rdbw2d.distance(
    dat$y, distance, b = b, p = 1, kink.position = 2,
    masspoints = "off", bwcheck = NULL
  )
  adaptive.multiple <- rdbw2d.distance(
    dat$y, distance, b = b, p = 1, kink.position = c(1, 3),
    masspoints = "off", bwcheck = NULL
  )

  distance.to.kink <- abs(b$x.1)
  kink.rate <- smooth$bws$h0 * nrow(distance)^(-1 / 4) /
    nrow(distance)^(-1 / (2 * smooth$opt$p + 4))
  expected <- pmin(smooth$bws$h0, pmax(kink.rate, distance.to.kink))

  expect_equal(adaptive$opt$kink.position, kink.position)
  expect_equal(adaptive$bws$h0, expected, tolerance = 1e-12)
  expect_equal(adaptive$bws$h1, expected, tolerance = 1e-12)
  expect_equal(adaptive.integer$bws, adaptive$bws, tolerance = 1e-12)
  expect_equal(adaptive.multiple$opt$kink.position, c(TRUE, FALSE, TRUE))
})

test_that("distance kink options validate incompatible inputs", {
  dat <- make_distance_data(n = 300, seed = 20260529)
  b <- data.frame(x.1 = 0, x.2 = 0)

  expect_error(
    rdbw2d.distance(
      dat$y, dat$distance, b = b, kink.unknown = c(FALSE, TRUE),
      masspoints = "off", bwcheck = NULL
    ),
    "kink.unknown\\[2\\]"
  )
  expect_error(
    rdbw2d.distance(
      dat$y, dat$distance, kink.position = TRUE,
      masspoints = "off", bwcheck = NULL
    ),
    "requires b"
  )
  expect_error(
    rdbw2d.distance(
      dat$y, dat$distance, b = b, kink.position = c(TRUE, FALSE),
      masspoints = "off", bwcheck = NULL
    ),
    "one TRUE/FALSE"
  )
  expect_error(
    rdbw2d.distance(
      dat$y, dat$distance, b = b, kink.position = 2,
      masspoints = "off", bwcheck = NULL
    ),
    "between 1"
  )
  expect_error(
    rd2d.distance(
      dat$y, dat$distance, h = 0.45, b = b, kink.unknown = c(TRUE, FALSE),
      cbands = FALSE, masspoints = "off", bwcheck = NULL
    ),
    "automatic bandwidth"
  )
})

test_that("distance CER bandwidth selectors are rate-corrected MSE selectors", {
  dat <- make_distance_data(n = 500, seed = 20260530)
  b <- data.frame(x.1 = 0, x.2 = 0)
  p <- 1

  mse <- rdbw2d.distance(
    dat$y, dat$distance, b = b, p = p, bwselect = "mserd",
    masspoints = "off", bwcheck = NULL
  )
  cer <- rdbw2d.distance(
    dat$y, dat$distance, b = b, p = p, bwselect = "cerrd",
    masspoints = "off", bwcheck = NULL
  )
  two <- rdbw2d.distance(
    dat$y, dat$distance, b = b, p = p, bwselect = "msetwo",
    masspoints = "off", bwcheck = NULL
  )
  certwo <- rdbw2d.distance(
    dat$y, dat$distance, b = b, p = p, bwselect = "certwo",
    masspoints = "off", bwcheck = NULL
  )

  expect_equal(cer$bws[, c("h0", "h1")], mse$bws[, c("h0", "h1")] * rd2d_cer_factor(length(dat$y), p))
  expect_equal(certwo$bws$h0, two$bws$h0 * rd2d_cer_factor(sum(dat$x < 0), p))
  expect_equal(certwo$bws$h1, two$bws$h1 * rd2d_cer_factor(sum(dat$x >= 0), p))
})

test_that("rd2d.distance standard errors use signed scaled-basis residuals", {
  dat <- make_distance_data(n = 400, seed = 20260510)
  h <- 0.5
  fit <- rd2d.distance(
    dat$y, dat$distance, h = h, b = data.frame(x.1 = 0, x.2 = 0),
    p = 1, q = 1, cbands = TRUE, masspoints = "off",
    bwcheck = NULL, kernel = "tri", vce = "hc0"
  )

  manual.se <- sqrt(
    manual_dist_se(dat$y, dat$x, h, 1, FALSE)^2 +
      manual_dist_se(dat$y, dat$x, h, 1, TRUE)^2
  )

  expect_equal(fit$main$std.err.p, manual.se, tolerance = 1e-10)
  expect_equal(fit$main$std.err.q, manual.se, tolerance = 1e-10)
  expect_equal(diag(fit$params.cov$main), fit$main$std.err.q^2, tolerance = 1e-10)
})

test_that("rd2d.distance clustered covariance diagonal matches standard errors", {
  dat <- make_distance_data(n = 600, seed = 20260511)
  set.seed(20260511)
  cluster <- sample(seq_len(90), length(dat$y), replace = TRUE)

  fit <- rd2d.distance(
    dat$y, dat$distance, h = 0.45, b = data.frame(x.1 = 0, x.2 = 0),
    p = 1, q = 2, cbands = TRUE,
    masspoints = "off", bwcheck = NULL, kernel = "tri", vce = "hc1",
    cluster = cluster
  )

  expect_equal(diag(fit$params.cov$main), fit$main$std.err.q^2, tolerance = 1e-10)
})

test_that("summary.rd2d.distance reports WBATE and LBATE rows", {
  dat <- make_distance_data(n = 650, seed = 20260518)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)
  weights <- c(1, 2, 3)

  fit <- rd2d.distance(
    dat$y, distance, h = 0.55, b = b, p = 1, q = 1,
    cbands = TRUE, masspoints = "off", bwcheck = NULL, kernel = "tri",
    repp = 49
  )

  set.seed(20260518)
  capture.output(summ <- summary(fit, WBATE = weights, LBATE = TRUE))
  tab <- summ$tables$main

  expect_equal(tail(rownames(tab), 2), c("WBATE", "LBATE"))

  normalized.weights <- weights / sum(weights)
  wbate.se <- sqrt(as.numeric(
    t(normalized.weights) %*% fit$params.cov$main %*% normalized.weights
  ))
  wbate.center <- sum(normalized.weights * fit$main$estimate.q)
  wbate.cval <- qnorm((fit$opt$level + 100) / 200)

  expect_equal(
    tab["WBATE", "estimate.p"],
    sum(normalized.weights * fit$main$estimate.p),
    tolerance = 1e-12
  )
  expect_equal(tab["WBATE", "estimate.q"], wbate.center, tolerance = 1e-12)
  expect_equal(tab["WBATE", "std.err.q"], wbate.se, tolerance = 1e-12)
  expect_equal(tab["WBATE", "t.value"], wbate.center / wbate.se, tolerance = 1e-12)
  expect_equal(
    tab["WBATE", "p.value"],
    2 * pnorm(abs(wbate.center / wbate.se), lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    unname(unlist(tab["WBATE", c("ci.lower", "ci.upper")])),
    c(wbate.center - wbate.cval * wbate.se,
      wbate.center + wbate.cval * wbate.se),
    tolerance = 1e-12
  )

  set.seed(20260518)
  lbate.cval <- rd2d_cval(
    fit$params.cov$main, rep = fit$opt$repp, side = fit$opt$side,
    alpha = fit$opt$level, lp = Inf
  )
  pointwise.se <- sqrt(diag(fit$params.cov$main))

  expect_equal(tab["LBATE", "estimate.p"], max(fit$main$estimate.p), tolerance = 1e-12)
  expect_equal(tab["LBATE", "estimate.q"], max(fit$main$estimate.q), tolerance = 1e-12)
  expect_true(is.na(tab["LBATE", "std.err.q"]))
  expect_true(is.na(tab["LBATE", "t.value"]))
  expect_true(is.na(tab["LBATE", "p.value"]))
  expect_equal(
    unname(unlist(tab["LBATE", c("ci.lower", "ci.upper")])),
    c(max(fit$main$estimate.q - lbate.cval * pointwise.se),
      max(fit$main$estimate.q + lbate.cval * pointwise.se)),
    tolerance = 1e-12
  )
})

test_that("summary.rd2d.distance aggregates use the full boundary under subset", {
  dat <- make_distance_data(n = 650, seed = 20260519)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)
  weights <- c(1, 2, 3)

  fit <- rd2d.distance(
    dat$y, distance, h = 0.55, b = b, p = 1, q = 1,
    cbands = TRUE, masspoints = "off", bwcheck = NULL, kernel = "tri",
    repp = 49
  )

  set.seed(20260519)
  capture.output(summ.all <- summary(fit, WBATE = weights, LBATE = TRUE))
  set.seed(20260519)
  capture.output(summ.sub <- summary(
    fit, WBATE = weights, LBATE = TRUE, subset = 1
  ))

  aggregate.cols <- c(
    "estimate.p", "estimate.q", "std.err.q", "t.value", "p.value",
    "ci.lower", "ci.upper"
  )
  expect_equal(rownames(summ.sub$tables$main), c("1", "WBATE", "LBATE"))
  expect_equal(
    summ.sub$tables$main[c("WBATE", "LBATE"), aggregate.cols],
    summ.all$tables$main[c("WBATE", "LBATE"), aggregate.cols]
  )
})

test_that("summary.rd2d.distance aggregate inference requires stored covariance", {
  dat <- make_distance_data(n = 650, seed = 20260520)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)

  fit <- rd2d.distance(
    dat$y, distance, h = 0.55, b = b, p = 1, q = 1,
    cbands = FALSE, masspoints = "off", bwcheck = NULL, kernel = "tri"
  )

  expect_error(
    capture.output(summary(fit, WBATE = rep(1, 3))),
    "WBATE inference requires.*stored covariance matrix"
  )
  expect_error(
    capture.output(summary(fit, LBATE = TRUE)),
    "LBATE inference requires.*stored covariance matrix"
  )
})

test_that("summary.rd2d.distance ignores WBATE and LBATE for bandwidth tables", {
  dat <- make_distance_data(n = 650, seed = 20260521)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)

  fit <- rd2d.distance(
    dat$y, distance, h = 0.55, b = b, p = 1, q = 1,
    cbands = FALSE, masspoints = "off", bwcheck = NULL, kernel = "tri"
  )

  expect_warning(
    capture.output(summ <- summary(
      fit, output = "bw", WBATE = rep(1, 3), LBATE = TRUE
    )),
    "WBATE and LBATE are not available for output='bw'"
  )
  expect_equal(summ$outputs, "bw")
  expect_false(any(c("WBATE", "LBATE") %in% rownames(summ$tables$bw)))
})

test_that("fuzzy rd2d.distance returns main, itt, fs, and requested side outputs", {
  dat <- make_distance_data(n = 700, seed = 20260522)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)

  fit <- rd2d.distance(
    dat$y.fuzzy, distance, h = 0.55, b = b, p = 1, q = 1,
    fuzzy = dat$fuzzy, params.other = c("itt.0", "itt.1", "fs.0", "fs.1"),
    params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1"),
    cbands = FALSE, masspoints = "off", bwcheck = NULL, kernel = "tri"
  )

  expect_true(isTRUE(fit$opt$fuzzy))
  expect_equal(
    names(fit),
    c(
      "main", "bw", "itt", "itt.0", "itt.1", "fs", "fs.0", "fs.1",
      "opt", "tau.hat", "tau.hat.q", "se.hat", "se.hat.q",
      "params.cov", "cb", "pvalues", "tvalues", "tau.itt", "tau.itt.q",
      "tau.fs", "tau.fs.q", "rdmodel", "call"
    )
  )
  for (output in c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")) {
    expect_s3_class(fit[[output]], "data.frame")
    expect_equal(names(fit[[output]]), names(fit$main))
    expect_equal(
      diag(fit$params.cov[[output]]),
      fit[[output]]$std.err.q^2,
      tolerance = 1e-10
    )
  }
  expect_equal(fit$itt$estimate.p, fit$tau.itt, tolerance = 1e-12)
  expect_equal(fit$fs$estimate.p, fit$tau.fs, tolerance = 1e-12)
  expect_equal(fit$main$estimate.p, fit$tau.itt / fit$tau.fs, tolerance = 1e-12)
  expect_equal(fit$main$t.value, fit$main$estimate.q / fit$main$std.err.q)
})

test_that("fuzzy rdbw2d.distance supports main and itt bandwidth targets", {
  dat <- make_distance_data(n = 700, seed = 20260523)
  b <- data.frame(x.1 = 0, x.2 = 0)

  bw.main <- rdbw2d.distance(
    dat$y.fuzzy, dat$distance, b = b, p = 1, fuzzy = dat$fuzzy,
    bwparam = "main", masspoints = "off", bwcheck = NULL
  )
  bw.itt <- rdbw2d.distance(
    dat$y.fuzzy, dat$distance, b = b, p = 1, fuzzy = dat$fuzzy,
    bwparam = "itt", masspoints = "off", bwcheck = NULL
  )
  sharp <- rdbw2d.distance(
    dat$y.fuzzy, dat$distance, b = b, p = 1, bwparam = "itt",
    masspoints = "off", bwcheck = NULL
  )

  expect_true(isTRUE(bw.main$opt$fuzzy))
  expect_true(isTRUE(bw.itt$opt$fuzzy))
  expect_equal(bw.main$opt$bwparam, "main")
  expect_equal(bw.itt$opt$bwparam, "itt")
  expect_equal(sharp$opt$bwparam, "main")
  expect_true(all(is.finite(as.matrix(bw.main$bws[, c("h0", "h1")]))))
  expect_true(all(is.finite(as.matrix(bw.itt$bws[, c("h0", "h1")]))))
})

test_that("summary.rd2d.distance reports fuzzy WBATE and LBATE rows", {
  dat <- make_distance_data(n = 700, seed = 20260524)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)
  weights <- c(1, 2, 3)

  fit <- rd2d.distance(
    dat$y.fuzzy, distance, h = 0.55, b = b, p = 1, q = 1,
    fuzzy = dat$fuzzy, params.cov = c("main", "itt", "fs"),
    cbands = FALSE, masspoints = "off", bwcheck = NULL, kernel = "tri",
    repp = 49
  )

  capture.output(summ <- summary(
    fit, output = c("main", "itt", "fs"), WBATE = weights, LBATE = TRUE
  ))
  expect_equal(names(summ$tables), c("main", "itt", "fs"))
  for (output in c("main", "itt", "fs")) {
    tab <- summ$tables[[output]]
    expect_equal(tail(rownames(tab), 2), c("WBATE", "LBATE"))

    normalized.weights <- weights / sum(weights)
    wbate.se <- sqrt(as.numeric(
      t(normalized.weights) %*% fit$params.cov[[output]] %*% normalized.weights
    ))
    wbate.center <- sum(normalized.weights * fit[[output]]$estimate.q)

    expect_equal(
      tab["WBATE", "estimate.p"],
      sum(normalized.weights * fit[[output]]$estimate.p),
      tolerance = 1e-12
    )
    expect_equal(tab["WBATE", "estimate.q"], wbate.center, tolerance = 1e-12)
    expect_equal(tab["WBATE", "std.err.q"], wbate.se, tolerance = 1e-12)
    expect_equal(tab["WBATE", "t.value"], wbate.center / wbate.se, tolerance = 1e-12)
    expect_equal(tab["LBATE", "estimate.p"], max(fit[[output]]$estimate.p), tolerance = 1e-12)
    expect_equal(tab["LBATE", "estimate.q"], max(fit[[output]]$estimate.q), tolerance = 1e-12)
  }
})

test_that("fuzzy rd2d.distance aggregate inference requires requested covariance", {
  dat <- make_distance_data(n = 700, seed = 20260525)
  distance <- cbind(dat$x, 1.1 * dat$x, 0.9 * dat$x)
  b <- data.frame(x.1 = c(-0.2, 0, 0.2), x.2 = 0)

  fit <- rd2d.distance(
    dat$y.fuzzy, distance, h = 0.55, b = b, p = 1, q = 1,
    fuzzy = dat$fuzzy, cbands = FALSE, masspoints = "off", bwcheck = NULL,
    kernel = "tri"
  )

  expect_error(
    capture.output(summary(fit, WBATE = rep(1, 3))),
    "WBATE inference requires.*params.cov = \"main\""
  )
  expect_error(
    capture.output(summary(fit, output = "itt", WBATE = rep(1, 3))),
    "WBATE inference requires.*params.cov = \"itt\""
  )
})

test_that("distance polynomial orders are nonnegative integers", {
  dat <- make_distance_data()

  expect_output(
    expect_error(
      rd2d.distance(dat$y, dat$distance, h = 0.45, p = 1.5, masspoints = "off")
    ),
    "p must be a nonnegative integer"
  )
  expect_output(
    expect_error(
      rd2d.distance(dat$y, dat$distance, h = 0.45, p = 2, q = 1, masspoints = "off")
    ),
    "q must be greater than or equal to p"
  )
  expect_output(
    expect_error(
      rdbw2d.distance(dat$y, dat$distance, p = 1.5, masspoints = "off")
    ),
    "p must be a nonnegative integer"
  )
})
