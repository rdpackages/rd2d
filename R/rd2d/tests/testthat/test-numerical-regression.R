make_regression_data <- function(n = 240, seed = 20260507) {
  set.seed(seed)
  x1 <- round(rnorm(n, sd = 1.4), 1)
  x2 <- round(0.35 * x1 + rnorm(n, sd = 1.2), 1)
  z <- as.numeric(x1 + 0.25 * x2 >= 0)
  fuzzy <- as.numeric(runif(n) < ifelse(z == 1, 0.8, 0.2))
  y <- 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * z + rnorm(n, sd = 0.6)
  y.fuzzy <- 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rnorm(n, sd = 0.6)
  cov <- cbind(
    cov1 = 0.6 * x1 - 0.2 * x2 + rnorm(n, sd = 0.5),
    cov2 = x1^2 - 0.3 * x2 + rnorm(n, sd = 0.5)
  )
  b <- matrix(
    c(
      -0.2, -0.3,
       0.0,  0.0,
       0.3,  0.2
    ),
    ncol = 2,
    byrow = TRUE
  )

  list(
    y = y, y.fuzzy = y.fuzzy, x = cbind(x1, x2), z = z, fuzzy = fuzzy,
    cov = cov, b = b
  )
}

expect_rd2d_tables_equal <- function(auto, manual, outputs) {
  for (output in outputs) {
    expect_equal(auto[[output]], manual[[output]], tolerance = 1e-10, ignore_attr = TRUE)
  }
  expect_equal(auto$bw, manual$bw, tolerance = 1e-10, ignore_attr = TRUE)
}

expect_shared_columns_equal <- function(joint, separate, output) {
  stable.cols <- c(
    "b1", "b2", "estimate.p", "estimate.q", "h01", "h02",
    "h11", "h12", "N.Co", "N.Tr"
  )
  expect_equal(
    joint[[output]][, stable.cols],
    separate[[output]][, stable.cols],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
}

test_that("joint and separate fitmethods agree without clusters", {
  dat <- make_regression_data(seed = 20260531)
  h <- cbind(
    c(0.85, 0.90, 0.95),
    c(0.75, 0.80, 0.85),
    c(0.95, 1.00, 1.05),
    c(0.70, 0.76, 0.82)
  )

  joint <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = h, deriv = c(1, 0),
    p = 2, q = 3, vce = "hc0", masspoints = "off", bwcheck = NULL,
    fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  separate <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = h, deriv = c(1, 0),
    p = 2, q = 3, vce = "hc0", masspoints = "off", bwcheck = NULL,
    fitmethod = "separate",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  expect_rd2d_tables_equal(joint, separate, c("main", "main.0", "main.1"))
  expect_equal(joint$params.cov, separate$params.cov, tolerance = 1e-10)
})

test_that("joint HC1 uses joint degrees-of-freedom correction without clusters", {
  dat <- make_regression_data(seed = 20260600)

  joint <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  separate <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "separate",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  for (output in c("main", "main.0", "main.1")) {
    expect_shared_columns_equal(joint, separate, output)
    expect_false(isTRUE(all.equal(
      joint[[output]]$std.err.q, separate[[output]]$std.err.q,
      tolerance = 1e-10
    )))
    expect_equal(
      diag(joint$params.cov[[output]]),
      joint[[output]]$std.err.q^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }
})

test_that("joint clustered fits use joint degrees-of-freedom with side-disjoint clusters", {
  dat <- make_regression_data(seed = 20260601)
  set.seed(20260601)
  cluster.base <- sample(seq_len(60), length(dat$y), replace = TRUE)
  cluster <- ifelse(dat$z == 1, paste0("tr", cluster.base), paste0("co", cluster.base))

  joint <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  separate <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "separate",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  for (output in c("main", "main.0", "main.1")) {
    expect_shared_columns_equal(joint, separate, output)
    expect_false(isTRUE(all.equal(
      joint$params.cov[[output]], separate$params.cov[[output]],
      tolerance = 1e-10
    )))
    expect_equal(
      diag(joint$params.cov[[output]]),
      joint[[output]]$std.err.q^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }
})

test_that("joint clustered covariance uses cross-side cluster scores", {
  dat <- make_regression_data(seed = 20260602)
  set.seed(20260602)
  cluster <- sample(seq_len(45), length(dat$y), replace = TRUE)
  cluster.effect <- rnorm(45, sd = 1.5)
  dat$y <- dat$y + cluster.effect[cluster]
  crosses.side <- tapply(dat$z, cluster, function(x) length(unique(x)) > 1)
  expect_true(any(unlist(crosses.side)))

  joint <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  separate <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "separate",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  expect_shared_columns_equal(joint, separate, "main")
  expect_false(isTRUE(all.equal(
    joint$params.cov$main, separate$params.cov$main, tolerance = 1e-10
  )))
  expect_equal(
    diag(joint$params.cov$main),
    joint$main$std.err.q^2,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  for (output in c("main.0", "main.1")) {
    expect_shared_columns_equal(joint, separate, output)
    expect_false(isTRUE(all.equal(
      joint$params.cov[[output]], separate$params.cov[[output]],
      tolerance = 1e-10
    )))
  }
})

test_that("joint fuzzy clustered covariance changes only combined outputs", {
  dat <- make_regression_data(seed = 20260603)
  set.seed(20260603)
  cluster <- sample(seq_len(50), length(dat$y.fuzzy), replace = TRUE)
  cluster.effect <- rnorm(50, sd = 1.4)
  dat$y.fuzzy <- dat$y.fuzzy + cluster.effect[cluster]

  outputs <- c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
  side.outputs <- c("itt.0", "itt.1", "fs.0", "fs.1")
  joint <- suppressWarnings(rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.95,
    fuzzy = dat$fuzzy, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = side.outputs, params.cov = outputs
  ))
  separate <- suppressWarnings(rd2d(
    dat$y.fuzzy, dat$x, dat$z, dat$b, h = 0.95,
    fuzzy = dat$fuzzy, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, fitmethod = "separate",
    params.other = side.outputs, params.cov = outputs
  ))

  for (output in c("main", "itt", "fs")) {
    expect_shared_columns_equal(joint, separate, output)
    expect_false(isTRUE(all.equal(
      joint$params.cov[[output]], separate$params.cov[[output]],
      tolerance = 1e-10
    )))
    expect_equal(
      diag(joint$params.cov[[output]]),
      joint[[output]]$std.err.q^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }
  for (output in side.outputs) {
    expect_shared_columns_equal(joint, separate, output)
    expect_false(isTRUE(all.equal(
      joint$params.cov[[output]], separate$params.cov[[output]],
      tolerance = 1e-10
    )))
    expect_equal(
      diag(joint$params.cov[[output]]),
      joint[[output]]$std.err.q^2,
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
  }
})

test_that("rdbw2d joint and separate agree for nonclustered HC0", {
  dat <- make_regression_data(seed = 20260604)

  joint <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc0", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, fitmethod = "joint"
  ))
  separate <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc0", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, fitmethod = "separate"
  ))

  expect_equal(joint$bws, separate$bws, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(
    joint$mseconsts, separate$mseconsts, tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("rdbw2d joint HC1 changes variance constants but not bias constants", {
  dat <- make_regression_data(seed = 20260605)

  joint <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc1", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, method = "rot", fitmethod = "joint"
  ))
  separate <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc1", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, method = "rot", fitmethod = "separate"
  ))

  stable.cols <- c("bias.0", "bias.1", "reg.bias.0", "reg.bias.1")
  expect_equal(
    joint$mseconsts[, stable.cols],
    separate$mseconsts[, stable.cols],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_false(isTRUE(all.equal(
    joint$mseconsts[, c("var.0", "var.1", "reg.var.0", "reg.var.1")],
    separate$mseconsts[, c("var.0", "var.1", "reg.var.0", "reg.var.1")],
    tolerance = 1e-10
  )))
  expect_equal(joint$mseconsts$var.01, rep(0, nrow(joint$mseconsts)))
  expect_equal(joint$mseconsts$reg.var.01, rep(0, nrow(joint$mseconsts)))
})

test_that("rdbw2d joint clustered selector stores cross-side score covariance", {
  dat <- make_regression_data(seed = 20260606)
  set.seed(20260606)
  cluster <- sample(seq_len(30), length(dat$y), replace = TRUE)
  cluster.effect <- rnorm(30, sd = 1.2)
  dat$y <- dat$y + cluster.effect[cluster]
  crosses.side <- tapply(dat$z, cluster, function(x) length(unique(x)) > 1)
  expect_true(any(unlist(crosses.side)))

  joint <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    fitmethod = "joint"
  ))
  separate <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, cluster = cluster, vce = "hc1",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    fitmethod = "separate"
  ))

  expect_true(any(abs(joint$mseconsts$var.01) > 1e-12))
  expect_equal(separate$mseconsts$var.01, rep(0, nrow(separate$mseconsts)))
})

test_that("rd2d automatic bandwidths pass fitmethod through to rdbw2d", {
  dat <- make_regression_data(seed = 20260607)

  bw <- suppressWarnings(rdbw2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc1", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, scaleregul = 3,
    fitmethod = "separate"
  ))
  fit <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, vce = "hc1", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, fitmethod = "separate"
  ))

  expect_equal(fit$opt$fitmethod, "separate")
  expect_equal(
    fit$bw[, c("h01", "h02", "h11", "h12")],
    bw$bws[, c("h01", "h02", "h11", "h12")],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("zero covs.eff leaves fixed-bandwidth HC0 sharp results unchanged", {
  dat <- make_regression_data(seed = 20260608)
  z0 <- matrix(0, nrow = length(dat$y), ncol = 1)

  plain <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, vce = "hc0",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  adjusted <- suppressWarnings(rd2d(
    dat$y, dat$x, dat$z, dat$b, h = 0.95, covs.eff = z0, vce = "hc0",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  expect_rd2d_tables_equal(adjusted, plain, c("main", "main.0", "main.1"))
  expect_equal(adjusted$params.cov, plain$params.cov, tolerance = 1e-10)
  expect_true(adjusted$opt$covs.eff)
  expect_equal(adjusted$opt$N.covs.eff, 1)
})

test_that("covs.eff changes covariate-predictive sharp estimates", {
  dat <- make_regression_data(seed = 20260609)
  y.cov <- dat$y + 1.2 * dat$cov[, 1] - 0.8 * dat$cov[, 2]

  plain <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, vce = "hc0",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint"
  ))
  adjusted <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, covs.eff = dat$cov, vce = "hc0",
    masspoints = "off", bwcheck = NULL, fitmethod = "joint"
  ))

  expect_false(isTRUE(all.equal(
    plain$main$estimate.p, adjusted$main$estimate.p, tolerance = 1e-8
  )))
  expect_true(adjusted$opt$covs.eff)
  expect_equal(adjusted$opt$N.covs.eff, ncol(dat$cov))
})

test_that("joint and separate covs.eff fits agree for nonclustered HC0", {
  dat <- make_regression_data(seed = 20260610)
  y.cov <- dat$y + as.numeric(dat$cov %*% c(0.7, -0.4))

  joint <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, covs.eff = dat$cov,
    vce = "hc0", masspoints = "off", bwcheck = NULL,
    fitmethod = "joint", params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  separate <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, covs.eff = dat$cov,
    vce = "hc0", masspoints = "off", bwcheck = NULL,
    fitmethod = "separate", params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))

  expect_rd2d_tables_equal(joint, separate, c("main", "main.0", "main.1"))
  expect_equal(joint$params.cov, separate$params.cov, tolerance = 1e-10)
})

test_that("rank-deficient covs.eff is dropped without changing location fits", {
  dat <- make_regression_data(seed = 20260615)
  y.cov <- dat$y + as.numeric(dat$cov %*% c(0.8, -0.6))
  cov.unique <- dat$cov
  cov.dup <- cbind(dat$cov[, 1], dat$cov[, 1], dat$cov[, 2])

  unique.fit <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, covs.eff = cov.unique,
    vce = "hc1", masspoints = "off", bwcheck = NULL, fitmethod = "joint",
    params.other = c("main.0", "main.1"),
    params.cov = c("main", "main.0", "main.1")
  ))
  expect_warning(
    dup.fit <- rd2d(
      y.cov, dat$x, dat$z, dat$b, h = 0.95, covs.eff = cov.dup,
      vce = "hc1", masspoints = "off", bwcheck = NULL, fitmethod = "joint",
      params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1")
    ),
    "covs.eff is rank deficient"
  )

  expect_true(dup.fit$opt$covs.rank.deficient)
  expect_equal(dup.fit$opt$N.covs.eff, 3)
  expect_equal(dup.fit$opt$N.covs.used, 2)
  expect_length(dup.fit$opt$covs.redundant, 1)
  expect_length(dup.fit$opt$covs.dropped, 1)
  expect_rd2d_tables_equal(dup.fit, unique.fit, c("main", "main.0", "main.1"))
  expect_equal(dup.fit$params.cov, unique.fit$params.cov, tolerance = 1e-10)
})

test_that("rank-deficient covs.eff can use generalized inverse in location fits", {
  dat <- make_regression_data(seed = 20260616)
  cov.dup <- cbind(dat$cov[, 1], dat$cov[, 1], dat$cov[, 2])

  expect_warning(
    fit <- rd2d(
      dat$y, dat$x, dat$z, dat$b, h = 0.95, covs.eff = cov.dup,
      covs.drop = FALSE, vce = "hc0", masspoints = "off",
      bwcheck = NULL, fitmethod = "joint"
    ),
    "generalized inverse"
  )

  expect_true(fit$opt$covs.rank.deficient)
  expect_equal(fit$opt$N.covs.used, 2)
  expect_length(fit$opt$covs.redundant, 1)
  expect_length(fit$opt$covs.dropped, 0)
})

test_that("automatic bandwidths pass covs.eff through to rdbw2d", {
  dat <- make_regression_data(seed = 20260611)
  y.cov <- dat$y + as.numeric(dat$cov %*% c(0.9, -0.5))

  bw <- suppressWarnings(rdbw2d(
    y.cov, dat$x, dat$z, dat$b, covs.eff = dat$cov, vce = "hc0",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    scaleregul = 3, fitmethod = "joint"
  ))
  fit <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, covs.eff = dat$cov, vce = "hc0",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    fitmethod = "joint"
  ))

  expect_equal(fit$opt$covs.eff, TRUE)
  expect_equal(
    fit$bw[, c("h01", "h02", "h11", "h12")],
    bw$bws[, c("h01", "h02", "h11", "h12")],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("rdbw2d covs.eff affects bandwidth constants when covariates matter", {
  dat <- make_regression_data(seed = 20260612)
  y.cov <- dat$y + as.numeric(dat$cov %*% c(1.1, -0.7))

  plain <- suppressWarnings(rdbw2d(
    y.cov, dat$x, dat$z, dat$b, vce = "hc0", masspoints = "off",
    bwcheck = NULL, stdvars = FALSE, fitmethod = "joint"
  ))
  adjusted <- suppressWarnings(rdbw2d(
    y.cov, dat$x, dat$z, dat$b, covs.eff = dat$cov, vce = "hc0",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    fitmethod = "joint"
  ))

  expect_false(isTRUE(all.equal(
    plain$mseconsts[, c("bias.0", "bias.1", "var.0", "var.1")],
    adjusted$mseconsts[, c("bias.0", "bias.1", "var.0", "var.1")],
    tolerance = 1e-8
  )))
  expect_true(adjusted$opt$covs.eff)
})

test_that("rank-deficient covs.eff is dropped without changing location bandwidths", {
  dat <- make_regression_data(seed = 20260617)
  y.cov <- dat$y + as.numeric(dat$cov %*% c(1.1, -0.7))
  cov.dup <- cbind(dat$cov[, 1], dat$cov[, 1], dat$cov[, 2])

  unique.bw <- suppressWarnings(rdbw2d(
    y.cov, dat$x, dat$z, dat$b, covs.eff = dat$cov, vce = "hc1",
    masspoints = "off", bwcheck = NULL, stdvars = FALSE,
    fitmethod = "joint"
  ))
  expect_warning(
    dup.bw <- rdbw2d(
      y.cov, dat$x, dat$z, dat$b, covs.eff = cov.dup, vce = "hc1",
      masspoints = "off", bwcheck = NULL, stdvars = FALSE,
      fitmethod = "joint"
    ),
    "covs.eff is rank deficient"
  )

  expect_true(dup.bw$opt$covs.rank.deficient)
  expect_equal(dup.bw$opt$N.covs.used, 2)
  expect_equal(dup.bw$bws, unique.bw$bws, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(
    dup.bw$mseconsts, unique.bw$mseconsts, tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("fuzzy covs.eff fixed-bandwidth fit returns finite component estimates", {
  dat <- make_regression_data(seed = 20260613)
  y.cov <- dat$y.fuzzy + as.numeric(dat$cov %*% c(0.8, -0.6))

  fit <- suppressWarnings(rd2d(
    y.cov, dat$x, dat$z, dat$b, h = 0.95, fuzzy = dat$fuzzy,
    covs.eff = dat$cov, vce = "hc0", masspoints = "off",
    bwcheck = NULL, fitmethod = "joint"
  ))

  expect_true(all(is.finite(fit$main$estimate.p)))
  expect_true(all(is.finite(fit$itt$estimate.p)))
  expect_true(all(is.finite(fit$fs$estimate.p)))
  expect_true(fit$opt$covs.eff)
})

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
      params.cov = c("main", "main.0", "main.1")
    ))
    h <- as.matrix(auto$bw[, c("h01", "h02", "h11", "h12")])
    manual <- suppressWarnings(rd2d(
      dat$y, dat$x, dat$z, dat$b,
      h = h,
      masspoints = scenario$masspoints,
      kernel_type = scenario$kernel_type,
      bwcheck = 8,
      params.other = c("main.0", "main.1"),
      params.cov = c("main", "main.0", "main.1")
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
      params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
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
      params.cov = c("main", "itt", "fs", "itt.0", "itt.1", "fs.0", "fs.1")
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
