
########################### Inverting Matrices #################################

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  mat <- crossprod(x)

  invMatrix <- tryCatch({
    chol2inv(chol(mat))
  },
  error = function(e) {
    if (grepl("leading minor of order", e$message)) {
      # If error is due to non-invertibility, use generalized inverse.
      warning(
        paste(
          "Calculating inverse of (t(X)%*%X), matrix is not",
          "positive-definite. Using generalized inverse."
        )
      )
      return(ginv(mat))
    } else {
      # If it's another error, just stop and propagate the error
      stop(e)
    }
  })

  return(invMatrix)
}

rd2d_warn_vce_small_en <- function(vce, eN, k, context = NULL) {
  if (!(vce %in% c("hc1", "hc2", "hc3"))) return(invisible(FALSE))
  if (is.finite(eN) && is.finite(k) && eN <= k) {
    where <- if (is.null(context)) "" else paste0(" ", context)
    warning(
      sprintf(
        paste(
          "vce='%s' may be undefined%s because the effective sample size",
          "eN=%d is not larger than the number of basis terms k=%d.",
          "Increase the bandwidth or bwcheck, lower p/q, or use vce='hc0'."
        ),
        vce, where, eN, k
      ),
      call. = FALSE
    )
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

rd2d_warn_vce_high_leverage <- function(vce, hii, context = NULL) {
  if (!(vce %in% c("hc2", "hc3"))) return(invisible(FALSE))
  bad <- !is.finite(hii) | hii >= 1
  if (any(bad, na.rm = TRUE)) {
    where <- if (is.null(context)) "" else paste0(" ", context)
    warning(
      sprintf(
        paste(
          "vce='%s' may be undefined%s because at least one leverage value",
          "is greater than or equal to 1. Increase the bandwidth or bwcheck,",
          "lower p/q, or use vce='hc0'."
        ),
        vce, where
      ),
      call. = FALSE
    )
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

rd2d_bwselect_base <- function(bwselect) {
  switch(
    bwselect,
    cerrd = "mserd",
    certwo = "msetwo",
    icerrd = "imserd",
    icertwo = "imsetwo",
    bwselect
  )
}

rd2d_bwselect_is_cer <- function(bwselect) {
  bwselect %in% c("cerrd", "certwo", "icerrd", "icertwo")
}

rd2d_bwselect_is_common <- function(bwselect) {
  rd2d_bwselect_base(bwselect) %in% c("mserd", "imserd")
}

rd2d_bw_rate_factor <- function(n, from.exp, to.exp) {
  n^(from.exp - to.exp)
}

rd2d_cer_factor <- function(n, p) {
  rd2d_bw_rate_factor(n, 1 / (2 * p + 4), 1 / (p + 4))
}

rd2d_validate_kink_unknown <- function(kink.unknown) {
  if (is.null(kink.unknown)) kink.unknown <- c(FALSE, FALSE)
  if (!is.logical(kink.unknown) && !is.numeric(kink.unknown)) {
    stop("kink.unknown must be a logical value or logical vector of length 2.", call. = FALSE)
  }
  if (!(length(kink.unknown) %in% c(1, 2)) || anyNA(kink.unknown)) {
    stop("kink.unknown must be a logical value or logical vector of length 2.", call. = FALSE)
  }
  kink.unknown <- as.logical(kink.unknown)
  if (length(kink.unknown) == 1) {
    kink.unknown <- rep(kink.unknown, 2)
  }
  if (!kink.unknown[1] && kink.unknown[2]) {
    stop("kink.unknown[2] can be TRUE only when kink.unknown[1] is TRUE.", call. = FALSE)
  }
  kink.unknown
}

rd2d_validate_kink_position <- function(kink.position, neval, b) {
  if (is.null(kink.position)) return(rep(FALSE, neval))
  if (is.logical(kink.position)) {
    if (length(kink.position) != neval || anyNA(kink.position)) {
      stop("kink.position must have one TRUE/FALSE value for each boundary point.", call. = FALSE)
    }
    out <- kink.position
  } else if (is.numeric(kink.position)) {
    if (anyNA(kink.position) || any(!is.finite(kink.position)) ||
        any(kink.position < 1) || any(kink.position > neval) ||
        any(abs(kink.position - round(kink.position)) > sqrt(.Machine$double.eps))) {
      stop("kink.position integer indices must be between 1 and the number of boundary points.", call. = FALSE)
    }
    out <- rep(FALSE, neval)
    out[as.integer(round(kink.position))] <- TRUE
  } else {
    stop("kink.position must be a logical vector or integer index vector.", call. = FALSE)
  }
  if (any(out) && is.null(b)) {
    stop("kink.position requires b so distances between boundary points can be computed.", call. = FALSE)
  }
  out
}

rd2d_distance_to_known_kink <- function(eval, kink.position) {
  if (!any(kink.position)) return(rep(Inf, nrow(eval)))
  kink.eval <- eval[kink.position, , drop = FALSE]
  dist.to.kink <- sqrt(
    outer(eval$x.1, kink.eval$x.1, "-")^2 +
      outer(eval$x.2, kink.eval$x.2, "-")^2
  )
  apply(dist.to.kink, 1, min)
}

########################### Weight Calculation #################################

kernel_weight = function(u,kernel){
  if (kernel == "triangular" || kernel == "tri") {
    abs.u <- abs(u)
    return((1 - abs.u) * (abs.u <= 1))
  }
  if (kernel == "uniform" || kernel == "uni") {
    abs.u <- abs(u)
    return(0.5 * (abs.u <= 1))
  }
  if (kernel == "gaussian" || kernel == "gau") {
    return(dnorm(u))
  }
  if (kernel == "epanechnikov" || kernel == "epa") {
    abs.u <- abs(u)
    return(0.75 * (1 - u^2) * (abs.u <= 1))
  }
  kernel <- tolower(kernel)
  if (kernel == "triangular" || kernel == "tri") {
    abs.u <- abs(u)
    return((1 - abs.u) * (abs.u <= 1))
  }
  if (kernel == "uniform" || kernel == "uni") {
    abs.u <- abs(u)
    return(0.5 * (abs.u <= 1))
  }
  if (kernel == "gaussian" || kernel == "gau") {
    return(dnorm(u))
  }
  if (kernel == "epanechnikov" || kernel == "epa") {
    abs.u <- abs(u)
    return(0.75 * (1 - u^2) * (abs.u <= 1))
  }
  abs.u <- abs(u)
  0.75 * (1 - u^2) * (abs.u <= 1)
}

##################### Rule of Thumb Bandwidth Selection ########################

rdbw2d_rot <- function(x,kernel.type, M){

  mu2K.squared <- NA
  l2K.squared <- NA

  if (kernel.type == "Epanechnikov"){
    mu2K.squared <- 1/6
    l2K.squared <- 4/(3 * pi)
  }
  if (kernel.type == "Triangular"){
    mu2K.squared <- 3/20
    l2K.squared <- 3/(2 * pi)
  }
  if (kernel.type == "Uniform"){mu2K.squared <- 1/4; l2K.squared <- 1/(pi)}
  if (kernel.type == "Gaussian"){mu2K.squared <- 1; l2K.squared <- 1/(4 * pi)}

  # Data cleaning

  x <- x[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
  x <- x[na.ok,]
  N <- dim(x)[1]

  # Estimate sample variance.

  cov.matrix <- cov(x[,c("x.1", "x.2")])

  dim.x <- 2

  inv.cov <- tryCatch(
    solve(cov.matrix),
    error = function(e) ginv(cov.matrix, 1e-20)
  )
  det.cov <- det(cov.matrix)
  sqrt.det.cov <- if (is.finite(det.cov) && det.cov > 0) {
    sqrt(det.cov)
  } else {
    det(sqrtm(cov.matrix))
  }
  trace.const <- 1 / (2^(dim.x + 2) * pi^(dim.x / 2) * sqrt.det.cov) *
    (2 * sum(diag(inv.cov %*% inv.cov)) + (sum(diag(inv.cov)))^2)

  if (is.null(M)){
    hROT <- ((dim.x * l2K.squared) /
      (N * mu2K.squared * trace.const))^(1 / (4 + dim.x))
  } else{
    hROT <- ((dim.x * l2K.squared) /
      (M * mu2K.squared * trace.const))^(1 / (4 + dim.x))
  }

  return(hROT)
}

###### Bandwidth Selection using Global Fit for Higher Order Derivatives #######

rdbw2d_effective_weights <- function(dat, h, kernel, kernel_type) {
  h <- rd2d_h_normalize(h, kernel_type)
  if (kernel_type == "prod") {
    kernel_weight(dat$x.1 / h[1], kernel) *
      kernel_weight(dat$x.2 / h[2], kernel) /
      (h[1] * h[2])
  } else {
    kernel_weight(dat$distance / h, kernel) / h^2
  }
}

rdbw2d_joint_effective_info <- function(dat.0, dat.1, h.0, h.1,
                                        kernel, kernel_type,
                                        cluster.0 = NULL,
                                        cluster.1 = NULL) {
  w.0 <- rdbw2d_effective_weights(dat.0, h.0, kernel, kernel_type)
  w.1 <- rdbw2d_effective_weights(dat.1, h.1, kernel, kernel_type)
  ind.0 <- as.logical(w.0 > 0)
  ind.1 <- as.logical(w.1 > 0)
  eC <- NULL
  if (!is.null(cluster.0) || !is.null(cluster.1)) {
    eC <- c(cluster.0[ind.0], cluster.1[ind.1])
  }

  list(eN = sum(ind.0) + sum(ind.1), eC = eC)
}

rdbw2d_bw_v2 <- function(dat.centered, p, vec, dn, bn.1,
                         bn.2 = NULL, vce, kernel, kernel_type, cluster,
                         fitmethod = c("separate", "joint"),
                         clusters = NULL, joint.info.v = NULL,
                         joint.info.b = NULL, cov.names = NULL,
                         gamma.v = NULL, gamma.b = NULL,
                         gamma.t = NULL){
  fitmethod <- match.arg(fitmethod)
  joint <- fitmethod == "joint"
  clustered <- !is.null(cluster)
  if (is.null(clusters)) clusters <- unique(cluster)
  covariate.adjusted <- !is.null(cov.names) && length(cov.names) > 0
  n.cov.v <- rd2d_covs_gamma_rank(gamma.v, cov.names)
  n.cov.b <- rd2d_covs_gamma_rank(gamma.b, cov.names)
  vce.local <- if ((joint || covariate.adjusted) && vce == "hc1") {
    "hc0"
  } else {
    vce
  }
  cluster.df <- !(joint && clustered) && !(covariate.adjusted && clustered)

  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "distance", cov.names)]
  dat.v <- rd2d_apply_covariate_adjustment(
    dat.centered, "y", cov.names, gamma.v
  )

  # Variance and coefficients for a linear combination of (p+1)-th derivatives.
  if (kernel_type == "prod"){
    w.v <- kernel_weight(dat.centered$x.1/c(dn), kernel) *
      kernel_weight(dat.centered$x.2/c(dn), kernel) / c(dn^2)
  }
  else{
    w.v <- kernel_weight(dat.centered$distance/c(dn), kernel)/c(dn^2)
  }

  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)

  ew.v <- w.v[ind.v]
  eY.v <- dat.v$y[ind.v]

  eC.v <- cluster[ind.v]

  eu.x1.v <- dat.centered$x.1[ind.v] / dn
  eu.x2.v <- dat.centered$x.2[ind.v] / dn

  if (is.null(bn.2)){
    eR.v.aug <- get_basis_xy(eu.x1.v, eu.x2.v, p + 1)
    p.count <- factorial(p+2)/(factorial(p)*2)
    p1.count <- factorial(p+1+2)/(factorial(p+1)*2)
    eR.v <- eR.v.aug[,1:p.count, drop = FALSE]
    eS.v <- eR.v.aug[,(p.count + 1):p1.count, drop = FALSE]
  } else {
    eR.v.aug <- get_basis_xy(eu.x1.v, eu.x2.v, p + 2)
    p.count <- factorial(p+2)/(factorial(p)*2)
    p1.count <- factorial(p+1+2)/(factorial(p+1)*2)
    p2.count <- factorial(p+2+2)/(factorial(p+2)*2)
    eR.v <- eR.v.aug[,1:p.count, drop = FALSE]
    eS.v <- eR.v.aug[,(p.count + 1):p1.count, drop = FALSE]
    eT.v <- eR.v.aug[,(p1.count + 1):p2.count, drop = FALSE]
  }

  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v

  w_R.v <- ew.v * eR.v
  k.v <- ncol(eR.v)

  invG.v <- qrXXinv(sqrtw_R.v)

  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)

  if (!is.null(bn.2)){
    sqrtw_eT.v <- sqrt(ew.v) * eT.v
    vec.t <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eT.v
    vec.t <- vec.t[1,]
    vec.t <- c(rep(0, factorial(p + 1 + 2)/(factorial(p + 1)*2)), vec.t)
  }

  invH.p <- get_invH(dn,p)

  H.p <- get_H(dn,p)

  beta.v <- invH.p %*% invG.v %*% t(sqrtw_R.v) %*% matrix(sqrtw_Y.v, ncol = 1)

  resd.v <- (eY.v - eR.v %*% (H.p %*% beta.v))[,1]

  if (vce.local=="hc0") {
    w.vce = 1
  } else if (vce.local=="hc1") {
    rd2d_warn_vce_small_en(vce.local, eN.v, k.v, "in the bandwidth selector")
    w.vce = sqrt(eN.v/(eN.v - k.v))
  } else if (vce.local=="hc2") {
    rd2d_warn_vce_small_en(vce.local, eN.v, k.v, "in the bandwidth selector")
    hii <- rd2d_leverage(sqrtw_R.v, invG.v)
    rd2d_warn_vce_high_leverage(vce.local, hii, "in the bandwidth selector")
    w.vce = sqrt(1/(1-hii))
  } else if (vce.local=="hc3"){
    rd2d_warn_vce_small_en(vce.local, eN.v, k.v, "in the bandwidth selector")
    hii <- rd2d_leverage(sqrtw_R.v, invG.v)
    rd2d_warn_vce_high_leverage(vce.local, hii, "in the bandwidth selector")
    w.vce = 1/(1-hii)
  }

  resd.v <- resd.v * w.vce

  if (clustered) {
    g.eff <- length(unique(eC.v))
    cluster.scale <- if (cluster.df) {
      ((eN.v - 1) / (eN.v - k.v)) * (g.eff / (g.eff - 1))
    } else {
      1
    }
    half.v <- rd2d_cluster_sums(w_R.v * resd.v, eC.v, clusters) *
      sqrt(dn^2 * cluster.scale)
    half.v <- half.v %*% invG.v
  } else {
    half.v <- sweep(as.matrix(sqrtw_R.v), 1, sqrt(ew.v) * resd.v, "*") *
      dn
    half.v <- half.v %*% invG.v
  }
  if (joint && (vce == "hc1" || clustered)) {
    if (is.null(joint.info.v)) joint.info.v <- list(eN = eN.v, eC = eC.v)
    half.v <- half.v * rd2d_joint_scale(
      vce, joint.info.v$eN, 2 * k.v + n.cov.v, joint.info.v$eC, clustered
    )
  } else if (covariate.adjusted && (vce == "hc1" || clustered)) {
    half.v <- half.v * rd2d_joint_scale(
      vce, eN.v, k.v + n.cov.v, eC.v, clustered
    )
  }
  V.half <- half.v %*% as.matrix(vec)
  V.V <- as.numeric(crossprod(V.half))

  # Bias

  dat.b <- rd2d_apply_covariate_adjustment(
    dat.centered, "y", cov.names, gamma.b
  )
  fit.ppls1 <- rd2d_lm(dat.b, bn.1, p + 1, vce.local, kernel = kernel,
                       kernel_type = kernel_type, cluster = cluster, varr = FALSE)

  deriv.ppls1 <- fit.ppls1$beta
  B.B <- vec.q %*% deriv.ppls1
  B.B <- B.B[1,1]
  bias.half.obj <- get_cov_half_v2(
    dat.b, bn.1, p + 1, vce.local, kernel, kernel_type, cluster,
    clusters, cluster.df = cluster.df
  )
  bias.half <- bias.half.obj$cov.half.const
  if (joint && (vce == "hc1" || clustered)) {
    if (is.null(joint.info.b)) {
      joint.info.b <- list(
        eN = sum(bias.half.obj$ind),
        eC = if (clustered) cluster[bias.half.obj$ind] else NULL
      )
    }
    bias.half <- bias.half * rd2d_joint_scale(
      vce, joint.info.b$eN, 2 * ncol(bias.half) + n.cov.b, joint.info.b$eC,
      clustered
    )
  } else if (covariate.adjusted && (vce == "hc1" || clustered)) {
    bias.half <- bias.half * rd2d_joint_scale(
      vce, sum(bias.half.obj$ind), ncol(bias.half) + n.cov.b,
      if (clustered) cluster[bias.half.obj$ind] else NULL, clustered
    )
  }
  Reg.1.half <- bias.half %*% matrix(vec.q, ncol = 1) /
    (bn.1^(p + 2))
  V.B <- as.numeric(crossprod(Reg.1.half))

  Reg.v <- V.B

  Reg.b <- NA

  if (!is.null(bn.2)){
    dat.t <- rd2d_apply_covariate_adjustment(
      dat.centered, "y", cov.names, gamma.t
    )
    fit.ppls2 <- rd2d_lm(dat.t, bn.2, p + 2, vce.local, kernel = kernel,
                         kernel_type = kernel_type, cluster = cluster, varr = FALSE)
    deriv.ppls2 <- fit.ppls2$beta
    Reg.b <- dn * vec.t %*% deriv.ppls2
  }

  return(list(
    "B" = B.B, "V" = V.V, "Reg.2" = Reg.b, "Reg.1" = Reg.v,
    "V.half" = V.half, "Reg.1.half" = Reg.1.half
  ))
}


##################### Fitting Local Polynomial Estimators ######################

# kernel_type = "prod" or "rad"

# o = 2 means taking out all second order derivatives.

# hgrid.0 either neval or neval by 2.

# If user does not provide hgrid.1, use one bandwidth for both groups.

rd2d_kth_max <- function(x, k) {
  n <- length(x)
  if (n < k) {
    return(c(kth = NA_real_, max = NA_real_, n = n))
  }
  kth <- if (k == n) max(x) else sort.int(x, partial = k)[k]
  c(kth = kth, max = max(x), n = n)
}

rd2d_bwcheck_bounds <- function(dat, eval, bwcheck, masspoints,
                                unique = NULL, metric = c("prod", "rad"),
                                scale = c(1, 1)) {
  if (is.null(bwcheck)) return(NULL)

  metric <- match.arg(metric)
  source <- if (masspoints == "adjust") {
    if (is.null(unique)) rd2d_unique(dat)$unique else unique
  } else {
    dat
  }
  source <- source[, c("x.1", "x.2", "d")]

  side.0 <- source$d == 0
  side.1 <- source$d == 1
  neval <- nrow(eval)
  bounds <- matrix(NA_real_, nrow = neval, ncol = 6)
  colnames(bounds) <- c("min.0", "min.1", "max.0", "max.1", "n.0", "n.1")

  for (i in seq_len(neval)) {
    dx1 <- source$x.1 - eval$x.1[i]
    dx2 <- source$x.2 - eval$x.2[i]
    distance <- if (metric == "prod") {
      pmax(abs(dx1 / scale[1]), abs(dx2 / scale[2]))
    } else {
      sqrt(dx1^2 + dx2^2)
    }

    b0 <- rd2d_kth_max(distance[side.0], bwcheck)
    b1 <- rd2d_kth_max(distance[side.1], bwcheck)
    if (is.na(b0["kth"]) || is.na(b1["kth"])) {
      unit.label <- ifelse(masspoints == "adjust", "unique mass points",
                           "observations")
      stop(
        sprintf(
          paste(
            "bwcheck=%d is larger than the available %s at evaluation",
            "point %d: %d on the control side and %d on the treatment side.",
            "Decrease bwcheck or set bwcheck=NULL."
          ),
          bwcheck, unit.label, i, b0["n"], b1["n"]
        ),
        call. = FALSE
      )
    }

    bounds[i, ] <- c(b0["kth"], b1["kth"], b0["max"], b1["max"],
                    b0["n"], b1["n"])
  }

  bounds
}

rd2d_fit <- function(dat, eval, deriv = NULL,o = 0, p = 1,
                     hgrid.0, hgrid.1 = NULL,
                     kernel = "epa", kernel_type = "prod", vce = "hc1",
                     masspoints = "adjust", cluster = NULL,
                     bwcheck = 50 + p + 1, unique = NULL,
                     bwcheck.bounds = NULL, cov.names = NULL){

  dat <- dat[,c("x.1", "x.2", "y", "d", cov.names)]
  na.ok <- complete.cases(dat)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]
  N <- dim(dat)[1]

  if (!is.null(bwcheck)) {
    if (!is.numeric(bwcheck) || length(bwcheck) != 1 ||
        !is.finite(bwcheck) || bwcheck < 1 || bwcheck != as.integer(bwcheck)) {
      stop("bwcheck must be NULL or a positive integer", call. = FALSE)
    }
    bwcheck <- as.integer(bwcheck)
  }

  sd.x1 <- sd(dat$x.1)
  sd.x2 <- sd(dat$x.2)

  neval <- dim(eval)[1]

  if (!is.null(bwcheck) && is.null(bwcheck.bounds)) {
    bwcheck.bounds <- rd2d_bwcheck_bounds(
      dat, eval, bwcheck, masspoints, unique = unique, metric = kernel_type,
      scale = c(sd.x1, sd.x2)
    )
  }

  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }

  if (ncol(hgrid.0) == 1){
    results <- data.frame(matrix(NA, ncol = 10, nrow = neval))
    colnames(results) <- c(
      'ev.x.1', 'ev.x.2', 'h.0', 'h.1', 'mu.0', 'mu.1', 'se.0',
      'se.1', 'eN.0', 'eN.1'
    )
  } else {
    results <- data.frame(matrix(NA, ncol = 12, nrow = neval))
    colnames(results) <- c(
      'ev.x.1', 'ev.x.2', 'h.0.x', 'h.0.y', 'h.1.x', 'h.1.y',
      'mu.0', 'mu.1', 'se.0', 'se.1', 'eN.0', 'eN.1'
    )
  }

  side.0 <- dat$d == 0
  side.1 <- dat$d == 1
  cluster.0 <- cluster[side.0]
  cluster.1 <- cluster[side.1]
  dat.centered <- dat
  covariate.adjusted <- !is.null(cov.names)
  if (!covariate.adjusted) {
    dat.0 <- dat[side.0, , drop = FALSE]
    dat.1 <- dat[side.1, , drop = FALSE]
    x1.0 <- dat.0$x.1
    x2.0 <- dat.0$x.2
    x1.1 <- dat.1$x.1
    x2.1 <- dat.1$x.2
  }

  for (i in 1:neval){

    ev <- eval[i,]
    vec <- deriv[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]

    if (!is.null(bwcheck)){
      bw.min.0 <- bwcheck.bounds[i, "min.0"]
      bw.min.1 <- bwcheck.bounds[i, "min.1"]
      bw.max.0 <- bwcheck.bounds[i, "max.0"]
      bw.max.1 <- bwcheck.bounds[i, "max.1"]

      # convert to the original if using product kernel
      if (kernel_type == "prod"){
        multiplier <- c(sd.x1, sd.x2)
      } else {
        multiplier <- c(1,1)
      }

      bw.min.0 <- bw.min.0 * multiplier
      bw.min.1 <- bw.min.1 * multiplier
      bw.max.0 <- bw.max.0 * multiplier
      bw.max.1 <- bw.max.1 * multiplier

      if (!is.null(hgrid.1)){
        h.0     <- pmax(h.0, bw.min.0)
        h.1     <- pmax(h.1, bw.min.1)
        h.0     <- pmin(h.0, bw.max.0)
        h.1     <- pmin(h.1, bw.max.1)
      } else{
        h.0 <- pmax(h.0, bw.min.0,bw.min.1)
        h.0 <- pmin(h.0, pmax(bw.max.0, bw.max.1))
        h.1 <- h.0
      }
    }

    if (covariate.adjusted) {
      dat.centered$x.1 <- dat$x.1 - ev$x.1
      dat.centered$x.2 <- dat$x.2 - ev$x.2
      # For product kernel, standardize covariates to have the same sd.
      dat.centered$distance <- pmax(
        abs(dat.centered$x.1/sd.x1),
        abs(dat.centered$x.2/sd.x2)
      )
      if (kernel_type == "rad") {
        dat.centered$distance <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)
      }
      gamma <- rd2d_covariate_gamma(
        dat.centered, h.0, h.1, p, kernel, kernel_type, "y", cov.names
      )
      dat.fit <- rd2d_apply_covariate_adjustment(
        dat.centered, "y", cov.names, gamma
      )
      fit.0.p <- rd2d_lm(
        dat.fit[side.0,], h.0, p, vce, kernel = kernel,
        cluster = cluster.0, varr = TRUE, kernel_type = kernel_type
      )
      fit.1.p <- rd2d_lm(
        dat.fit[side.1,], h.1, p, vce, kernel = kernel,
        cluster = cluster.1, varr = TRUE, kernel_type = kernel_type
      )
    } else {
      dat.0$x.1 <- x1.0 - ev$x.1
      dat.0$x.2 <- x2.0 - ev$x.2
      dat.1$x.1 <- x1.1 - ev$x.1
      dat.1$x.2 <- x2.1 - ev$x.2
      dat.0$distance <- pmax(abs(dat.0$x.1/sd.x1), abs(dat.0$x.2/sd.x2))
      dat.1$distance <- pmax(abs(dat.1$x.1/sd.x1), abs(dat.1$x.2/sd.x2))
      if (kernel_type == "rad") {
        dat.0$distance <- sqrt(dat.0$x.1^2 + dat.0$x.2^2)
        dat.1$distance <- sqrt(dat.1$x.1^2 + dat.1$x.2^2)
      }
      fit.0.p <- rd2d_lm(
        dat.0, h.0, p, vce, kernel = kernel,
        cluster = cluster.0, varr = TRUE, kernel_type = kernel_type
      )
      fit.1.p <- rd2d_lm(
        dat.1, h.1, p, vce, kernel = kernel,
        cluster = cluster.1, varr = TRUE, kernel_type = kernel_type
      )
    }
    mu.0 <- (vec %*% fit.0.p$beta)[1,1]
    mu.1 <- (vec %*% fit.1.p$beta)[1,1]

    # standardize
    if (length(h.0) == 1){
      h.0.x <- as.numeric(h.0)
      h.0.y <- as.numeric(h.0)
    }  else {
      h.0.x <- as.numeric(h.0[1])
      h.0.y <- as.numeric(h.0[2])
    }
    if (length(h.1) == 1){
      h.1.x <- as.numeric(h.1)
      h.1.y <- as.numeric(h.1)
    }  else {
      h.1.x <- as.numeric(h.1[1])
      h.1.y <- as.numeric(h.1[2])
    }

    # standard deviation
    invH.0 <- get_invH(c(h.0.x, h.0.y),p)
    se.0 <- matrix(vec, nrow = 1) %*%
      invH.0 %*%
      fit.0.p$cov.const %*%
      invH.0 %*%
      matrix(vec, ncol = 1) /
      (h.0.x * h.0.y)
    # se.0 <- matrix(vec, nrow = 1) %*% fit.0.p$cov.const %*%
    #   matrix(vec, ncol = 1) / (N * h.0^(2 + 2 * o))
    se.0 <- sqrt(se.0[1,1])
    invH.1 <- get_invH(c(h.1.x, h.1.y),p)
    se.1 <- matrix(vec, nrow = 1) %*%
      invH.1 %*%
      fit.1.p$cov.const %*%
      invH.1 %*%
      matrix(vec, ncol = 1) /
      (h.1.x * h.1.y)
    # se.1 <- matrix(vec, nrow = 1) %*% fit.1.p$cov.const %*%
    #   matrix(vec, ncol = 1) / (N * h.1^(2 + 2 * o))
    se.1 <- sqrt(se.1[1,1])

    # effective sample size
    eN.0 <- fit.0.p$eN
    eN.1 <- fit.1.p$eN

    if (ncol(hgrid.0) == 1){
      results[i,] <- c(
        ev[1], ev[2], h.0, h.1, mu.0, mu.1, se.0, se.1, eN.0, eN.1
      )
    } else {
      results[i,] <- c(
        ev[1], ev[2], h.0.x, h.0.y, h.1.x, h.1.y, mu.0, mu.1,
        se.0, se.1, eN.0, eN.1
      )
    }
  }

  return(results)
}

rd2d_fit_multi <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                           hgrid.0, hgrid.1 = NULL,
                           kernel = "epa", kernel_type = "prod", vce = "hc1",
                           masspoints = "adjust", cluster = NULL,
                           bwcheck = 50 + p + 1, unique = NULL,
                           bwcheck.bounds = NULL,
                           outcomes = c("y", "fuzzy"),
                           cov.names = NULL){

  dat <- dat[, c("x.1", "x.2", "d", outcomes, cov.names)]
  na.ok <- complete.cases(dat)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]

  if (!is.null(bwcheck)) {
    if (!is.numeric(bwcheck) || length(bwcheck) != 1 ||
        !is.finite(bwcheck) || bwcheck < 1 || bwcheck != as.integer(bwcheck)) {
      stop("bwcheck must be NULL or a positive integer", call. = FALSE)
    }
    bwcheck <- as.integer(bwcheck)
  }

  sd.x1 <- sd(dat$x.1)
  sd.x2 <- sd(dat$x.2)
  neval <- dim(eval)[1]

  if (!is.null(bwcheck) && is.null(bwcheck.bounds)) {
    bwcheck.bounds <- rd2d_bwcheck_bounds(
      dat, eval, bwcheck, masspoints, unique = unique, metric = kernel_type,
      scale = c(sd.x1, sd.x2)
    )
  }

  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) hgrid.1 <- as.matrix(hgrid.1)

  make_results <- function() {
    if (ncol(hgrid.0) == 1){
      results <- data.frame(matrix(NA, ncol = 10, nrow = neval))
      colnames(results) <- c(
        'ev.x.1', 'ev.x.2', 'h.0', 'h.1', 'mu.0', 'mu.1', 'se.0',
        'se.1', 'eN.0', 'eN.1'
      )
    } else {
      results <- data.frame(matrix(NA, ncol = 12, nrow = neval))
      colnames(results) <- c(
        'ev.x.1', 'ev.x.2', 'h.0.x', 'h.0.y', 'h.1.x', 'h.1.y',
        'mu.0', 'mu.1', 'se.0', 'se.1', 'eN.0', 'eN.1'
      )
    }
    results
  }
  results <- vector("list", length(outcomes))
  names(results) <- outcomes
  for (outcome in outcomes) results[[outcome]] <- make_results()

  side.0 <- dat$d == 0
  side.1 <- dat$d == 1
  cluster.0 <- cluster[side.0]
  cluster.1 <- cluster[side.1]
  dat.centered <- dat
  covariate.adjusted <- !is.null(cov.names)
  if (!covariate.adjusted) {
    dat.0 <- dat[side.0, , drop = FALSE]
    dat.1 <- dat[side.1, , drop = FALSE]
    x1.0 <- dat.0$x.1
    x2.0 <- dat.0$x.2
    x1.1 <- dat.1$x.1
    x2.1 <- dat.1$x.2
  }

  for (i in 1:neval){
    ev <- eval[i,]
    vec <- deriv[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]

    if (!is.null(bwcheck)){
      bw.min.0 <- bwcheck.bounds[i, "min.0"]
      bw.min.1 <- bwcheck.bounds[i, "min.1"]
      bw.max.0 <- bwcheck.bounds[i, "max.0"]
      bw.max.1 <- bwcheck.bounds[i, "max.1"]

      if (kernel_type == "prod"){
        multiplier <- c(sd.x1, sd.x2)
      } else {
        multiplier <- c(1,1)
      }

      bw.min.0 <- bw.min.0 * multiplier
      bw.min.1 <- bw.min.1 * multiplier
      bw.max.0 <- bw.max.0 * multiplier
      bw.max.1 <- bw.max.1 * multiplier

      if (!is.null(hgrid.1)){
        h.0 <- pmax(h.0, bw.min.0)
        h.1 <- pmax(h.1, bw.min.1)
        h.0 <- pmin(h.0, bw.max.0)
        h.1 <- pmin(h.1, bw.max.1)
      } else{
        h.0 <- pmax(h.0, bw.min.0,bw.min.1)
        h.0 <- pmin(h.0, pmax(bw.max.0, bw.max.1))
        h.1 <- h.0
      }
    }

    if (covariate.adjusted) {
      dat.centered$x.1 <- dat$x.1 - ev$x.1
      dat.centered$x.2 <- dat$x.2 - ev$x.2
      dat.centered$distance <- pmax(
        abs(dat.centered$x.1/sd.x1),
        abs(dat.centered$x.2/sd.x2)
      )
      if (kernel_type == "rad") {
        dat.centered$distance <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)
      }
      gamma <- rd2d_covariate_gamma(
        dat.centered, h.0, h.1, p, kernel, kernel_type, outcomes, cov.names
      )
      dat.fit <- rd2d_apply_covariate_adjustment(
        dat.centered, outcomes, cov.names, gamma
      )
      fit.0.p <- rd2d_lm_multi(
        dat.fit[side.0,], h.0, p, vce, kernel = kernel,
        cluster = cluster.0, varr = TRUE, kernel_type = kernel_type,
        outcomes = outcomes
      )
      fit.1.p <- rd2d_lm_multi(
        dat.fit[side.1,], h.1, p, vce, kernel = kernel,
        cluster = cluster.1, varr = TRUE, kernel_type = kernel_type,
        outcomes = outcomes
      )
    } else {
      dat.0$x.1 <- x1.0 - ev$x.1
      dat.0$x.2 <- x2.0 - ev$x.2
      dat.1$x.1 <- x1.1 - ev$x.1
      dat.1$x.2 <- x2.1 - ev$x.2
      dat.0$distance <- pmax(abs(dat.0$x.1/sd.x1), abs(dat.0$x.2/sd.x2))
      dat.1$distance <- pmax(abs(dat.1$x.1/sd.x1), abs(dat.1$x.2/sd.x2))
      if (kernel_type == "rad") {
        dat.0$distance <- sqrt(dat.0$x.1^2 + dat.0$x.2^2)
        dat.1$distance <- sqrt(dat.1$x.1^2 + dat.1$x.2^2)
      }
      fit.0.p <- rd2d_lm_multi(
        dat.0, h.0, p, vce, kernel = kernel,
        cluster = cluster.0, varr = TRUE, kernel_type = kernel_type,
        outcomes = outcomes
      )
      fit.1.p <- rd2d_lm_multi(
        dat.1, h.1, p, vce, kernel = kernel,
        cluster = cluster.1, varr = TRUE, kernel_type = kernel_type,
        outcomes = outcomes
      )
    }

    if (length(h.0) == 1){
      h.0.x <- as.numeric(h.0)
      h.0.y <- as.numeric(h.0)
    } else {
      h.0.x <- as.numeric(h.0[1])
      h.0.y <- as.numeric(h.0[2])
    }
    if (length(h.1) == 1){
      h.1.x <- as.numeric(h.1)
      h.1.y <- as.numeric(h.1)
    } else {
      h.1.x <- as.numeric(h.1[1])
      h.1.y <- as.numeric(h.1[2])
    }

    invH.0 <- get_invH(c(h.0.x, h.0.y), p)
    invH.1 <- get_invH(c(h.1.x, h.1.y), p)

    for (outcome in outcomes) {
      mu.0 <- (vec %*% fit.0.p$beta[, outcome, drop = FALSE])[1,1]
      mu.1 <- (vec %*% fit.1.p$beta[, outcome, drop = FALSE])[1,1]

      se.0 <- matrix(vec, nrow = 1) %*%
        invH.0 %*%
        fit.0.p$cov.const[[outcome]] %*%
        invH.0 %*%
        matrix(vec, ncol = 1) /
        (h.0.x * h.0.y)
      se.0 <- sqrt(se.0[1,1])
      se.1 <- matrix(vec, nrow = 1) %*%
        invH.1 %*%
        fit.1.p$cov.const[[outcome]] %*%
        invH.1 %*%
        matrix(vec, ncol = 1) /
        (h.1.x * h.1.y)
      se.1 <- sqrt(se.1[1,1])

      if (ncol(hgrid.0) == 1){
        results[[outcome]][i,] <- c(
          ev[1], ev[2], h.0, h.1, mu.0, mu.1, se.0, se.1,
          fit.0.p$eN, fit.1.p$eN
        )
      } else {
        results[[outcome]][i,] <- c(
          ev[1], ev[2], h.0.x, h.0.y, h.1.x, h.1.y, mu.0, mu.1,
          se.0, se.1, fit.0.p$eN, fit.1.p$eN
        )
      }
    }
  }

  results
}

######################### Estimating Covariance matrix  ########################

rd2d_hxy <- function(h){
  if (length(h) == 1) {
    c(as.numeric(h), as.numeric(h))
  } else {
    c(as.numeric(h[1]), as.numeric(h[2]))
  }
}

rd2d_prepare_cov_eff <- function(covs.eff, n) {
  if (is.null(covs.eff)) {
    return(list(matrix = NULL, names = NULL))
  }

  if (is.null(dim(covs.eff))) {
    cov.raw <- matrix(covs.eff, ncol = 1)
  } else if (is.matrix(covs.eff) || is.data.frame(covs.eff)) {
    cov.raw <- as.matrix(covs.eff)
  } else {
    stop("covs.eff must be NULL, a numeric vector, matrix, or data frame.",
         call. = FALSE)
  }

  if (nrow(cov.raw) != n) {
    stop("covs.eff must have the same number of observations as Y.",
         call. = FALSE)
  }
  if (ncol(cov.raw) < 1) {
    stop("covs.eff must contain at least one covariate column.", call. = FALSE)
  }

  cov.num <- suppressWarnings(matrix(
    as.numeric(cov.raw),
    nrow = nrow(cov.raw),
    ncol = ncol(cov.raw)
  ))
  bad.numeric <- is.na(cov.num) & !is.na(cov.raw)
  if (any(bad.numeric)) {
    stop("covs.eff must contain numeric covariates.", call. = FALSE)
  }

  cov.names <- paste0("z.", seq_len(ncol(cov.num)))
  colnames(cov.num) <- cov.names
  list(matrix = cov.num, names = cov.names)
}

rd2d_stop_errors <- function(errors) {
  errors <- errors[nzchar(errors)]
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }
}

.rd2d_covs_rank_env <- new.env(parent = emptyenv())
.rd2d_covs_rank_env$depth <- 0L
.rd2d_covs_rank_env$covs.drop <- TRUE
.rd2d_covs_rank_env$covs.tol <- 1e-12

rd2d_validate_covs_rank_options <- function(covs.drop, covs.tol) {
  errors <- character(0)
  if (!is.logical(covs.drop) || length(covs.drop) != 1 || is.na(covs.drop)) {
    errors <- c(errors, "covs.drop must be TRUE or FALSE.")
  }
  if (!is.numeric(covs.tol) || length(covs.tol) != 1 ||
      !is.finite(covs.tol) || covs.tol <= 0) {
    errors <- c(errors, "covs.tol must be a positive finite numeric value.")
  }
  errors
}

rd2d_covs_rank_begin <- function(cov.names, covs.drop, covs.tol) {
  depth <- .rd2d_covs_rank_env$depth
  if (is.null(depth) || depth == 0L) {
    supplied <- if (is.null(cov.names)) 0L else length(cov.names)
    .rd2d_covs_rank_env$supplied <- supplied
    .rd2d_covs_rank_env$used.min <- supplied
    .rd2d_covs_rank_env$redundant <- character(0)
    .rd2d_covs_rank_env$dropped <- character(0)
    .rd2d_covs_rank_env$rank.deficient <- FALSE
    .rd2d_covs_rank_env$warned <- FALSE
    .rd2d_covs_rank_env$covs.drop <- covs.drop
    .rd2d_covs_rank_env$covs.tol <- covs.tol
  }
  .rd2d_covs_rank_env$depth <- if (is.null(depth)) 1L else depth + 1L
  invisible(NULL)
}

rd2d_covs_rank_end <- function() {
  depth <- .rd2d_covs_rank_env$depth
  .rd2d_covs_rank_env$depth <- max(0L, if (is.null(depth)) 0L else depth - 1L)
  invisible(NULL)
}

rd2d_covs_rank_summary <- function() {
  supplied <- .rd2d_covs_rank_env$supplied
  if (is.null(supplied)) supplied <- 0L
  used <- .rd2d_covs_rank_env$used.min
  if (is.null(used)) used <- supplied
  list(
    supplied = supplied,
    used = used,
    redundant = .rd2d_covs_rank_env$redundant,
    dropped = .rd2d_covs_rank_env$dropped,
    rank.deficient = isTRUE(.rd2d_covs_rank_env$rank.deficient),
    drop = isTRUE(.rd2d_covs_rank_env$covs.drop),
    tol = .rd2d_covs_rank_env$covs.tol
  )
}

rd2d_covs_rank_option <- function(name, default) {
  value <- .rd2d_covs_rank_env[[name]]
  if (is.null(value)) default else value
}

rd2d_covs_rank_update <- function(cov.names, keep, covs.drop) {
  supplied <- length(cov.names)
  used <- length(keep)
  dropped <- cov.names[setdiff(seq_len(supplied), keep)]
  if (is.null(.rd2d_covs_rank_env$supplied)) {
    .rd2d_covs_rank_env$supplied <- supplied
  }
  if (is.null(.rd2d_covs_rank_env$used.min)) {
    .rd2d_covs_rank_env$used.min <- supplied
  }
  if (is.null(.rd2d_covs_rank_env$redundant)) {
    .rd2d_covs_rank_env$redundant <- character(0)
  }
  if (is.null(.rd2d_covs_rank_env$dropped)) {
    .rd2d_covs_rank_env$dropped <- character(0)
  }
  .rd2d_covs_rank_env$rank.deficient <- TRUE
  .rd2d_covs_rank_env$used.min <- min(.rd2d_covs_rank_env$used.min, used)
  .rd2d_covs_rank_env$redundant <- union(.rd2d_covs_rank_env$redundant, dropped)
  if (isTRUE(covs.drop)) {
    .rd2d_covs_rank_env$dropped <- union(.rd2d_covs_rank_env$dropped, dropped)
  }

  if (!isTRUE(.rd2d_covs_rank_env$warned)) {
    action <- if (isTRUE(covs.drop)) {
      "dropping redundant covariate columns"
    } else {
      "using a generalized inverse"
    }
    warning(
      paste(
        "covs.eff is rank deficient after residualizing on the local polynomial basis;",
        action,
        "for numerical stability."
      ),
      call. = FALSE
    )
    .rd2d_covs_rank_env$warned <- TRUE
  }
  invisible(NULL)
}

rd2d_covs_gamma_rank <- function(gamma, cov.names = NULL) {
  if (is.null(cov.names) || length(cov.names) == 0) return(0L)
  rank <- attr(gamma, "rank")
  if (is.null(rank) || !is.finite(rank)) length(cov.names) else as.integer(rank)
}

rd2d_print_covs_eff <- function(opt) {
  if (!isTRUE(opt$covs.eff)) return(invisible(NULL))
  cat(sprintf("Covariates             %d\n", opt$N.covs.eff))
  if (isTRUE(opt$covs.rank.deficient)) {
    cat(sprintf("Covariate rank         %d\n", opt$N.covs.used))
    if (length(opt$covs.dropped) > 0) {
      cat(sprintf("Covariates dropped     %s\n", paste(opt$covs.dropped, collapse = ", ")))
    }
  }
  invisible(NULL)
}

rd2d_covs_rank_solve <- function(ZWZ, ZWY, cov.names, covs.drop = NULL,
                                 covs.tol = NULL, ginv.tol = 1e-20) {
  covs.drop <- if (is.null(covs.drop)) {
    rd2d_covs_rank_option("covs.drop", TRUE)
  } else {
    covs.drop
  }
  covs.tol <- if (is.null(covs.tol)) {
    rd2d_covs_rank_option("covs.tol", 1e-12)
  } else {
    covs.tol
  }

  p <- ncol(ZWZ)
  qrZ <- qr(ZWZ, tol = covs.tol)
  rank <- qrZ$rank
  keep <- if (rank > 0L) qrZ$pivot[seq_len(rank)] else integer(0)

  if (rank < p) {
    rd2d_covs_rank_update(cov.names, keep, covs.drop)
    gamma <- matrix(0, nrow = p, ncol = ncol(ZWY))
    if (isTRUE(covs.drop)) {
      if (rank > 0L) {
        gamma.keep <- tryCatch(
          solve(ZWZ[keep, keep, drop = FALSE], ZWY[keep, , drop = FALSE]),
          error = function(e) {
            ginv(ZWZ[keep, keep, drop = FALSE], tol = ginv.tol) %*%
              ZWY[keep, , drop = FALSE]
          }
        )
        gamma[keep, ] <- gamma.keep
      }
    } else {
      gamma <- ginv(ZWZ, tol = ginv.tol) %*% ZWY
    }
  } else {
    gamma <- tryCatch(
      solve(ZWZ, ZWY),
      error = function(e) ginv(ZWZ, tol = ginv.tol) %*% ZWY
    )
  }

  gamma <- as.matrix(gamma)
  attr(gamma, "rank") <- rank
  attr(gamma, "kept") <- cov.names[keep]
  attr(gamma, "dropped") <- cov.names[setdiff(seq_len(p), keep)]
  gamma
}

rd2d_cluster_sums <- function(values, eC, clusters) {
  values <- as.matrix(values)
  out <- matrix(0, nrow = length(clusters), ncol = ncol(values))
  if (length(eC) == 0) return(out)

  groups <- match(eC, clusters)
  sums <- rowsum(values, group = groups, reorder = FALSE)
  out[as.integer(rownames(sums)), ] <- sums
  out
}

rd2d_covariate_gamma <- function(dat.centered, h.0, h.1, p, kernel,
                                 kernel_type, outcomes, cov.names,
                                 ginv.tol = 1e-20, covs.drop = NULL,
                                 covs.tol = NULL) {
  if (is.null(cov.names) || length(cov.names) == 0) return(NULL)

  outcomes <- as.character(outcomes)
  cov.names <- as.character(cov.names)
  side.0 <- dat.centered$d == 0
  side.1 <- dat.centered$d == 1

  side_components <- function(dat.side, h.side) {
    loc <- rd2d_local_design(
      dat.side, h.side, p, kernel, kernel_type, outcomes = c(outcomes, cov.names)
    )
    Y <- loc$eY[, outcomes, drop = FALSE]
    Z <- loc$eY[, cov.names, drop = FALSE]
    U <- crossprod(loc$w_R, loc$eY)
    U.y <- U[, outcomes, drop = FALSE]
    U.z <- U[, cov.names, drop = FALSE]

    list(
      ZWZ = crossprod(Z * loc$ew, Z) -
        crossprod(U.z, loc$invG %*% U.z),
      ZWY = crossprod(Z * loc$ew, Y) -
        crossprod(U.z, loc$invG %*% U.y)
    )
  }

  comp.0 <- side_components(dat.centered[side.0, , drop = FALSE], h.0)
  comp.1 <- side_components(dat.centered[side.1, , drop = FALSE], h.1)
  ZWZ <- comp.0$ZWZ + comp.1$ZWZ
  ZWY <- comp.0$ZWY + comp.1$ZWY

  gamma <- rd2d_covs_rank_solve(
    ZWZ, ZWY, cov.names, covs.drop = covs.drop, covs.tol = covs.tol,
    ginv.tol = ginv.tol
  )
  rownames(gamma) <- cov.names
  colnames(gamma) <- outcomes
  gamma
}

rd2d_apply_covariate_adjustment <- function(dat, outcomes, cov.names, gamma) {
  if (is.null(cov.names) || length(cov.names) == 0 || is.null(gamma)) {
    return(dat)
  }
  Z <- as.matrix(dat[, cov.names, drop = FALSE])
  for (outcome in outcomes) {
    dat[[outcome]] <- dat[[outcome]] -
      as.numeric(Z %*% gamma[, outcome, drop = FALSE])
  }
  dat
}

rd2d_project_cov_halves <- function(halves, inds, deriv, hgrid, p,
                                    clustered = FALSE, block = NULL) {
  neval <- length(halves)
  if (neval == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  n.rows <- if (clustered) nrow(halves[[1]]) else length(inds[[1]])
  projected <- matrix(0, nrow = n.rows, ncol = neval)

  for (i in seq_len(neval)) {
    half <- halves[[i]]
    if (is.null(half) || length(half) == 0) next

    h.xy <- rd2d_hxy(hgrid[i,])
    cols <- if (is.null(block)) seq_len(ncol(half)) else block
    score <- half[, cols, drop = FALSE] %*%
      get_invH(h.xy, p) %*%
      matrix(deriv[i,], ncol = 1)
    score <- as.numeric(score) / sqrt(prod(h.xy))

    if (clustered) {
      projected[, i] <- score
    } else {
      projected[inds[[i]], i] <- score
    }
  }

  projected
}

rd2d_joint_scale <- function(vce, eN, k, eC, clustered) {
  scale <- 1
  if (vce == "hc1") {
    rd2d_warn_vce_small_en(vce, eN, k, "in joint rd2d covariance")
    scale <- scale * sqrt(eN / (eN - k))
  }
  if (clustered) {
    g.eff <- length(unique(eC))
    scale <- scale * sqrt(((eN - 1) / (eN - k)) * (g.eff / (g.eff - 1)))
  }
  scale
}

rd2d_cov_project_sides <- function(dat, eval, deriv, p, hgrid.0,
                                   hgrid.1 = NULL, kernel = "epa",
                                   kernel_type = "prod", vce = "hc2",
                                   cluster = NULL,
                                   fitmethod = c("separate", "joint"),
                                   multi = FALSE, block = NULL,
                                   cov.names = NULL) {
  fitmethod <- match.arg(fitmethod)
  dat <- if (multi) {
    dat[, c("x.1", "x.2", "y", "d", "fuzzy", cov.names)]
  } else {
    dat[, c("x.1", "x.2", "y", "d", cov.names)]
  }
  na.ok <- complete.cases(dat)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]

  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) hgrid.1 <- as.matrix(hgrid.1)

  neval <- nrow(eval)
  halves.0 <- halves.1 <- vector("list", neval)
  inds.0 <- inds.1 <- vector("list", neval)
  eC.0 <- eC.1 <- vector("list", neval)
  eN.0 <- eN.1 <- integer(neval)

  side.0 <- dat$d == 0
  side.1 <- dat$d == 1
  cluster.0 <- cluster[side.0]
  cluster.1 <- cluster[side.1]
  clusters <- unique(cluster)
  clustered <- !is.null(cluster)
  joint <- fitmethod == "joint"
  covariate.adjusted <- !is.null(cov.names) && length(cov.names) > 0
  vce.half <- if ((joint || covariate.adjusted) && vce == "hc1") "hc0" else vce
  cluster.df <- !(joint && clustered) && !(covariate.adjusted && clustered)
  get_half <- if (multi) get_cov_half_multi_v2 else get_cov_half_v2
  dat.centered <- dat
  if (!covariate.adjusted) {
    dat.0 <- dat[side.0, , drop = FALSE]
    dat.1 <- dat[side.1, , drop = FALSE]
    x1.0 <- dat.0$x.1
    x2.0 <- dat.0$x.2
    x1.1 <- dat.1$x.1
    x2.1 <- dat.1$x.2
  }

  for (i in seq_len(neval)) {
    ev <- eval[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]

    n.cov <- 0L
    if (covariate.adjusted) {
      dat.centered$x.1 <- dat$x.1 - ev$x.1
      dat.centered$x.2 <- dat$x.2 - ev$x.2
      dat.centered$distance <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)

      outcomes <- if (multi) c("y", "fuzzy") else "y"
      gamma <- rd2d_covariate_gamma(
        dat.centered, h.0, h.1, p, kernel, kernel_type, outcomes, cov.names
      )
      n.cov <- rd2d_covs_gamma_rank(gamma, cov.names)
      dat.fit <- rd2d_apply_covariate_adjustment(
        dat.centered, outcomes, cov.names, gamma
      )
      half.0 <- get_half(
        dat.fit[side.0,], h.0, p, vce.half, kernel, kernel_type,
        cluster.0, clusters, cluster.df = cluster.df
      )
      half.1 <- get_half(
        dat.fit[side.1,], h.1, p, vce.half, kernel, kernel_type,
        cluster.1, clusters, cluster.df = cluster.df
      )
    } else {
      dat.0$x.1 <- x1.0 - ev$x.1
      dat.0$x.2 <- x2.0 - ev$x.2
      dat.0$distance <- sqrt(dat.0$x.1^2 + dat.0$x.2^2)
      dat.1$x.1 <- x1.1 - ev$x.1
      dat.1$x.2 <- x2.1 - ev$x.2
      dat.1$distance <- sqrt(dat.1$x.1^2 + dat.1$x.2^2)

      half.0 <- get_half(
        dat.0, h.0, p, vce.half, kernel, kernel_type,
        cluster.0, clusters, cluster.df = cluster.df
      )
      half.1 <- get_half(
        dat.1, h.1, p, vce.half, kernel, kernel_type,
        cluster.1, clusters, cluster.df = cluster.df
      )
    }

    inds.0[[i]] <- half.0$ind
    inds.1[[i]] <- half.1$ind
    halves.0[[i]] <- half.0$cov.half.const
    halves.1[[i]] <- half.1$cov.half.const
    eN.0[i] <- sum(half.0$ind)
    eN.1[i] <- sum(half.1$ind)
    if (clustered) {
      eC.0[[i]] <- cluster.0[half.0$ind]
      eC.1[[i]] <- cluster.1[half.1$ind]
    }
    if ((joint || covariate.adjusted) && (vce == "hc1" || clustered)) {
      eC.i <- if (clustered) c(eC.0[[i]], eC.1[[i]]) else NULL
      if (joint) {
        scale <- rd2d_joint_scale(
          vce, eN.0[i] + eN.1[i], 2 * ncol(deriv) + n.cov,
          eC.i, clustered
        )
        halves.0[[i]] <- halves.0[[i]] * scale
        halves.1[[i]] <- halves.1[[i]] * scale
      } else {
        scale.0 <- rd2d_joint_scale(
          vce, eN.0[i], ncol(deriv) + n.cov, eC.0[[i]], clustered
        )
        scale.1 <- rd2d_joint_scale(
          vce, eN.1[i], ncol(deriv) + n.cov, eC.1[[i]], clustered
        )
        halves.0[[i]] <- halves.0[[i]] * scale.0
        halves.1[[i]] <- halves.1[[i]] * scale.1
      }
    }
  }

  hgrid.1.use <- if (is.null(hgrid.1)) hgrid.0 else hgrid.1
  projected.0 <- projected.1 <- NULL
  if (!multi || !is.null(block)) {
    projected.0 <- rd2d_project_cov_halves(
      halves.0, inds.0, deriv, hgrid.0, p, clustered = clustered,
      block = block
    )
    projected.1 <- rd2d_project_cov_halves(
      halves.1, inds.1, deriv, hgrid.1.use, p, clustered = clustered,
      block = block
    )
  }

  list(
    projected.0 = projected.0,
    projected.1 = projected.1,
    halves.0 = halves.0,
    halves.1 = halves.1,
    inds.0 = inds.0,
    inds.1 = inds.1,
    clustered = clustered,
    joint = joint
  )
}

rdbw2d_cov <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                       hgrid.0, hgrid.1 = NULL, kernel = "epa",
                       kernel_type = "prod", vce = "hc2", cluster = NULL,
                       fitmethod = c("separate", "joint"),
                       cov.names = NULL){
  fitmethod <- match.arg(fitmethod)

  projects <- rd2d_cov_project_sides(
    dat, eval, deriv, p, hgrid.0, hgrid.1, kernel, kernel_type, vce,
    cluster, fitmethod, multi = FALSE, cov.names = cov.names
  )

  if (projects$clustered && projects$joint) {
    covs <- crossprod(projects$projected.1 - projects$projected.0)
  } else {
    covs <- crossprod(projects$projected.0) + crossprod(projects$projected.1)
  }
  covs <- (covs + t(covs)) / 2
  return(covs)
}


####################### Side-specific covariance matrix #######################

rdbw2d_cov_side <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                            hgrid.0, hgrid.1 = NULL, kernel = "epa",
                            kernel_type = "prod", vce = "hc2", cluster = NULL,
                            side = 0,
                            fitmethod = c("separate", "joint"),
                            cov.names = NULL){
  fitmethod <- match.arg(fitmethod)

  if (!(side %in% c(0, 1))) {
    stop("side must be 0 or 1.", call. = FALSE)
  }

  projects <- rd2d_cov_project_sides(
    dat, eval, deriv, p, hgrid.0, hgrid.1, kernel, kernel_type, vce,
    cluster, fitmethod, multi = FALSE, cov.names = cov.names
  )
  projected <- if (side == 0) projects$projected.0 else projects$projected.1
  covs <- crossprod(projected)
  covs <- (covs + t(covs)) / 2

  return(covs)
}


#################### Fuzzy RD covariance matrix estimation ####################

rdbw2d_cov_fuzzy <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                             hgrid.0, hgrid.1 = NULL, kernel = "epa",
                             kernel_type = "prod", vce = "hc2", cluster = NULL,
                             tau.itt, tau.fs,
                             denom.tol = sqrt(.Machine$double.eps),
                             fitmethod = c("separate", "joint"),
                             cov.names = NULL){
  fitmethod <- match.arg(fitmethod)

  dat <- dat[,c("x.1", "x.2", "y", "d", "fuzzy", cov.names)]
  na.ok <- complete.cases(dat)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]

  neval <- dim(eval)[1]
  covs <- matrix(NA_real_, nrow = neval, ncol = neval)
  valid <- is.finite(tau.itt) & is.finite(tau.fs) & abs(tau.fs) > denom.tol

  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }

  hgrid.1.use <- if (is.null(hgrid.1)) hgrid.0 else hgrid.1
  k <- ncol(deriv)
  y.idx <- seq_len(k)
  d.idx <- (k + 1):(2 * k)
  projects <- rd2d_cov_project_sides(
    dat, eval, deriv, p, hgrid.0, hgrid.1, kernel, kernel_type, vce,
    cluster, fitmethod, multi = TRUE, cov.names = cov.names
  )
  clustered <- projects$clustered

  y.0 <- rd2d_project_cov_halves(
    projects$halves.0, projects$inds.0, deriv, hgrid.0, p, clustered,
    block = y.idx
  )
  d.0 <- rd2d_project_cov_halves(
    projects$halves.0, projects$inds.0, deriv, hgrid.0, p, clustered,
    block = d.idx
  )
  y.1 <- rd2d_project_cov_halves(
    projects$halves.1, projects$inds.1, deriv, hgrid.1.use, p, clustered,
    block = y.idx
  )
  d.1 <- rd2d_project_cov_halves(
    projects$halves.1, projects$inds.1, deriv, hgrid.1.use, p, clustered,
    block = d.idx
  )

  if (clustered && fitmethod == "joint") {
    y <- y.1 - y.0
    d <- d.1 - d.0
    yy <- crossprod(y)
    yd <- crossprod(y, d)
    dy <- crossprod(d, y)
    dd <- crossprod(d)
  } else {
    yy <- crossprod(y.0) + crossprod(y.1)
    yd <- crossprod(y.0, d.0) + crossprod(y.1, d.1)
    dy <- crossprod(d.0, y.0) + crossprod(d.1, y.1)
    dd <- crossprod(d.0) + crossprod(d.1)
  }

  a <- 1 / tau.fs
  b <- -tau.itt / tau.fs^2
  covs <- outer(a, a) * yy +
    outer(a, b) * yd +
    outer(b, a) * dy +
    outer(b, b) * dd
  covs[!outer(valid, valid)] <- NA_real_
  covs <- (covs + t(covs)) / 2

  return(covs)
}

rdbw2d_cov_fuzzy_tables <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                                    hgrid.0, hgrid.1 = NULL, kernel = "epa",
                                    kernel_type = "prod", vce = "hc2",
                                    cluster = NULL, tau.itt, tau.fs, outputs,
                                    denom.tol = sqrt(.Machine$double.eps),
                                    fitmethod = c("separate", "joint"),
                                    cov.names = NULL){
  fitmethod <- match.arg(fitmethod)
  dat <- dat[,c("x.1", "x.2", "y", "d", "fuzzy", cov.names)]
  na.ok <- complete.cases(dat)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]

  neval <- dim(eval)[1]
  outputs <- unique(outputs)
  covs <- lapply(outputs, function(x) matrix(NA_real_, nrow = neval, ncol = neval))
  names(covs) <- outputs
  valid <- is.finite(tau.itt) & is.finite(tau.fs) & abs(tau.fs) > denom.tol

  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) hgrid.1 <- as.matrix(hgrid.1)

  hgrid.1.use <- if (is.null(hgrid.1)) hgrid.0 else hgrid.1
  k <- ncol(deriv)
  y.idx <- seq_len(k)
  d.idx <- (k + 1):(2 * k)
  projects <- rd2d_cov_project_sides(
    dat, eval, deriv, p, hgrid.0, hgrid.1, kernel, kernel_type, vce,
    cluster, fitmethod, multi = TRUE, cov.names = cov.names
  )
  clustered <- projects$clustered

  y.0 <- rd2d_project_cov_halves(
    projects$halves.0, projects$inds.0, deriv, hgrid.0, p, clustered,
    block = y.idx
  )
  d.0 <- rd2d_project_cov_halves(
    projects$halves.0, projects$inds.0, deriv, hgrid.0, p, clustered,
    block = d.idx
  )
  y.1 <- rd2d_project_cov_halves(
    projects$halves.1, projects$inds.1, deriv, hgrid.1.use, p, clustered,
    block = y.idx
  )
  d.1 <- rd2d_project_cov_halves(
    projects$halves.1, projects$inds.1, deriv, hgrid.1.use, p, clustered,
    block = d.idx
  )

  yy.0 <- crossprod(y.0)
  yy.1 <- crossprod(y.1)
  dd.0 <- crossprod(d.0)
  dd.1 <- crossprod(d.1)
  if (clustered && fitmethod == "joint") {
    y <- y.1 - y.0
    d <- d.1 - d.0
    yy <- crossprod(y)
    dd <- crossprod(d)
  } else {
    yy <- yy.0 + yy.1
    dd <- dd.0 + dd.1
  }

  if ("main" %in% outputs) {
    if (clustered && fitmethod == "joint") {
      yd <- crossprod(y, d)
      dy <- crossprod(d, y)
    } else {
      yd <- crossprod(y.0, d.0) + crossprod(y.1, d.1)
      dy <- crossprod(d.0, y.0) + crossprod(d.1, y.1)
    }
    a <- 1 / tau.fs
    b <- -tau.itt / tau.fs^2
    covs$main <- outer(a, a) * yy +
      outer(a, b) * yd +
      outer(b, a) * dy +
      outer(b, b) * dd
    covs$main[!outer(valid, valid)] <- NA_real_
    covs$main <- (covs$main + t(covs$main)) / 2
  }
  if ("itt" %in% outputs) covs$itt <- (yy + t(yy)) / 2
  if ("fs" %in% outputs) covs$fs <- (dd + t(dd)) / 2
  if ("itt.0" %in% outputs) covs$itt.0 <- (yy.0 + t(yy.0)) / 2
  if ("itt.1" %in% outputs) covs$itt.1 <- (yy.1 + t(yy.1)) / 2
  if ("fs.0" %in% outputs) covs$fs.0 <- (dd.0 + t(dd.0)) / 2
  if ("fs.1" %in% outputs) covs$fs.1 <- (dd.1 + t(dd.1)) / 2

  covs
}


#################### Fuzzy RD standard error estimation ########################

rdbw2d_se_fuzzy <- function(dat, eval, deriv = NULL, o = 0, p = 1,
                            hgrid.0, hgrid.1 = NULL, kernel = "epa",
                            kernel_type = "prod", vce = "hc2", cluster = NULL,
                            tau.itt, tau.fs, se.itt, se.fs,
                            denom.tol = sqrt(.Machine$double.eps),
                            fitmethod = c("separate", "joint"),
                            cov.names = NULL){
  fitmethod <- match.arg(fitmethod)

  if (fitmethod == "joint" || !is.null(cov.names)) {
    cov.main <- rdbw2d_cov_fuzzy(
      dat, eval, deriv, o, p, hgrid.0, hgrid.1, kernel, kernel_type, vce,
      cluster, tau.itt, tau.fs, denom.tol = denom.tol,
      fitmethod = fitmethod, cov.names = cov.names
    )
    return(sqrt(pmax(diag(cov.main), 0)))
  }

  dat <- dat[,c("x.1", "x.2", "y", "d", "fuzzy")]
  na.ok <- complete.cases(dat$x.1) &
    complete.cases(dat$x.2) &
    complete.cases(dat$y) &
    complete.cases(dat$d) &
    complete.cases(dat$fuzzy)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]

  neval <- dim(eval)[1]
  se <- rep(NA_real_, neval)
  cov.itt.fs <- rep(0, neval)
  valid <- is.finite(tau.itt) &
    is.finite(tau.fs) &
    is.finite(se.itt) &
    is.finite(se.fs) &
    abs(tau.fs) > denom.tol

  clusters <- unique(cluster)

  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }

  hxy <- function(h){
    if (length(h) == 1) {
      c(as.numeric(h), as.numeric(h))
    } else {
      c(as.numeric(h[1]), as.numeric(h[2]))
    }
  }

  side.0 <- dat$d == 0
  side.1 <- dat$d == 1
  dat.0 <- dat[side.0, c("x.1", "x.2", "y", "d", "fuzzy"), drop = FALSE]
  dat.1 <- dat[side.1, c("x.1", "x.2", "y", "d", "fuzzy"), drop = FALSE]
  x1.0 <- dat.0$x.1
  x2.0 <- dat.0$x.2
  x1.1 <- dat.1$x.1
  x2.1 <- dat.1$x.2
  cluster.0 <- if (is.null(cluster)) NULL else cluster[side.0]
  cluster.1 <- if (is.null(cluster)) NULL else cluster[side.1]

  for (i in 1:neval){
    if (!valid[i]) next

    ev.a <- eval[i,]
    h.a.0 <- hgrid.0[i,]
    h.a.1 <- h.a.0
    if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]

    dat.0$x.1 <- x1.0 - ev.a$x.1
    dat.0$x.2 <- x2.0 - ev.a$x.2
    dat.0$distance <- sqrt(dat.0$x.1^2 + dat.0$x.2^2)
    dat.1$x.1 <- x1.1 - ev.a$x.1
    dat.1$x.2 <- x2.1 - ev.a$x.2
    dat.1$distance <- sqrt(dat.1$x.1^2 + dat.1$x.2^2)

    cov.half.consts.0 <- get_cov_half_multi_v2(
      dat.0,
      h.a.0, p, vce, kernel, kernel_type, cluster.0, clusters
    )
    cov.half.consts.1 <- get_cov_half_multi_v2(
      dat.1,
      h.a.1, p, vce, kernel, kernel_type, cluster.1, clusters
    )

    deriv.vec <- deriv[i,]
    k <- length(deriv.vec)

    for (side in 0:1){
      if (side == 0){
        half <- cov.half.consts.0$cov.half.const
        h.xy <- hxy(h.a.0)
      } else {
        half <- cov.half.consts.1$cov.half.const
        h.xy <- hxy(h.a.1)
      }

      meat <- crossprod(half)
      invH <- get_invH(h.xy,p)
      scale <- prod(h.xy)
      y.idx <- 1:k
      w.idx <- (k + 1):(2 * k)
      cov.block <- matrix(deriv.vec, nrow = 1) %*% invH %*%
        meat[y.idx, w.idx, drop = FALSE] %*% invH %*%
        matrix(deriv.vec, ncol = 1)
      cov.itt.fs[i] <- cov.itt.fs[i] + cov.block[1,1] / scale
    }
  }

  var.ratio <- se.itt^2 / tau.fs^2 +
    tau.itt^2 * se.fs^2 / tau.fs^4 -
    2 * tau.itt * cov.itt.fs / tau.fs^3
  se[valid] <- sqrt(pmax(var.ratio[valid], 0))

  return(se)
}


################################## infl ########################################

infl <- function(x, invG){
  result <- t(matrix(x, ncol = 1)) %*% invG %*% matrix(x, ncol = 1)
  return(result[1,1])
}

rd2d_leverage <- function(sqrtw_R, invG) {
  rowSums((sqrtw_R %*% invG) * sqrtw_R)
}

rd2d_h_normalize <- function(h, kernel_type) {
  h <- if (is.data.frame(h)) {
    unlist(h, use.names = FALSE)
  } else {
    as.numeric(h)
  }
  if (kernel_type == "prod") {
    if (length(h) == 1) h <- c(h, h)
  } else {
    if (length(h) == 2) h <- sqrt(h[1]^2 + h[2]^2)
  }
  h
}

rd2d_hxy <- function(h) {
  if (length(h) == 1) {
    c(as.numeric(h), as.numeric(h))
  } else {
    c(as.numeric(h[1]), as.numeric(h[2]))
  }
}

rd2d_kernel_weights <- function(dat, h, kernel, kernel_type) {
  if (kernel_type == "prod") {
    kernel_weight(dat$x.1 / c(h[1]), kernel) *
      kernel_weight(dat$x.2 / c(h[2]), kernel) /
      c(h[1] * h[2])
  } else {
    kernel_weight(dat$distance / c(h), kernel) / c(h^2)
  }
}

rd2d_extract_outcomes <- function(dat, outcomes, ind, n = NULL) {
  outcomes <- as.character(outcomes)
  if (is.null(n)) n <- sum(ind)
  if (length(outcomes) == 1) {
    eY <- matrix(dat[[outcomes]][ind], ncol = 1)
    colnames(eY) <- outcomes
    return(eY)
  }

  eY <- matrix(NA_real_, nrow = n, ncol = length(outcomes))
  colnames(eY) <- outcomes
  for (j in seq_along(outcomes)) {
    eY[, j] <- dat[[outcomes[j]]][ind]
  }
  eY
}

rd2d_local_design <- function(dat, h, p, kernel, kernel_type, cluster = NULL,
                              outcomes = "y") {
  h <- rd2d_h_normalize(h, kernel_type)
  h.xy <- rd2d_hxy(h)
  x.1 <- dat$x.1
  x.2 <- dat$x.2
  if (kernel_type == "prod") {
    w <- kernel_weight(x.1 / h[1], kernel) *
      kernel_weight(x.2 / h[2], kernel) / (h[1] * h[2])
  } else {
    w <- kernel_weight(dat$distance / h, kernel) / h^2
  }
  ind <- as.logical(w > 0)
  ew <- w[ind]
  sqrt.ew <- sqrt(ew)
  eN <- length(ew)
  eY <- rd2d_extract_outcomes(dat, outcomes, ind, n = eN)
  eC <- cluster[ind]

  eR <- get_basis_xy(x.1[ind] / h.xy[1], x.2[ind] / h.xy[2], p)
  sqrtw_R <- sqrt.ew * eR
  w_R <- ew * eR
  invG <- qrXXinv(sqrtw_R)

  list(
    h = h,
    h.xy = h.xy,
    ind = ind,
    eN = eN,
    ew = ew,
    sqrt.ew = sqrt.ew,
    eY = eY,
    ecluster = eC,
    eR = eR,
    sqrtw_R = sqrtw_R,
    w_R = w_R,
    k = ncol(eR),
    invG = invG,
    H = get_H(h.xy, p),
    invH = get_invH(h.xy, p)
  )
}

######################### get_cov_half (efficient) #############################

get_cov_half_v2 <- function(dat, h, p, vce = "hc1", kernel = "epa",
                            kernel_type = "prod",
                            cluster = NULL, clusters = NULL,
                            cluster.df = TRUE){

  loc <- rd2d_local_design(dat, h, p, kernel, kernel_type, cluster, outcomes = "y")

  ind <- loc$ind
  eN <- loc$eN
  ew <- loc$ew
  eY <- loc$eY[, 1]
  eC <- loc$ecluster
  eR <- loc$eR
  sqrtw_R <- loc$sqrtw_R
  sqrtw_Y <- loc$sqrt.ew * eY
  w_R <- loc$w_R
  k <- loc$k
  invG <- loc$invG
  invH.p <- loc$invH
  H.p <- loc$H
  h.x <- loc$h.xy[1]
  h.y <- loc$h.xy[2]

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)

  resd <- (eY - eR %*% (H.p %*% beta))[,1]

    if (vce=="hc0") {
      w.vce <- 1
    } else if (vce=="hc1") { # TODO:Set to default
      rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_v2()")
      w.vce <- sqrt(eN/(eN - k))
    } else if (vce=="hc2") {
      rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_v2()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_v2()")
      w.vce <- sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_v2()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_v2()")
      w.vce <- 1/(1-hii)
    }

  resd <- resd * w.vce

  if (is.null(cluster)){

    sigma.half.const <- (sqrt(ew) * resd * as.matrix(sqrtw_R)) *
      sqrt(h.x * h.y)

    cov.half.const <- sigma.half.const %*% invG
  }

  if (!is.null(cluster)){

    n <- length(eC)
    g <- length(clusters)
    g.eff <- length(unique(eC))
    w.w <- if (cluster.df) {
      ((n-1)/(n-k))*(g.eff/(g.eff-1))
    } else {
      1
    }

    cov.half.const <- rd2d_cluster_sums(w_R * resd, eC, clusters) *
      sqrt(h.x * h.y * w.w)
    cov.half.const <- cov.half.const %*% invG
  }

  return(list("ind" = ind, "cov.half.const" = cov.half.const))
}

get_cov_half_multi_v2 <- function(dat, h, p, vce = "hc1",
                                  kernel = "epa", kernel_type = "prod",
                                  cluster = NULL, clusters = NULL,
                                  cluster.df = TRUE){

  loc <- rd2d_local_design(
    dat, h, p, kernel, kernel_type, cluster, outcomes = c("y", "fuzzy")
  )

  ind <- loc$ind
  eN <- loc$eN
  ew <- loc$ew
  eY <- loc$eY
  eC <- loc$ecluster
  eR <- loc$eR
  sqrtw_R <- loc$sqrtw_R
  sqrtw_Y <- sweep(eY, 1, loc$sqrt.ew, "*")
  w_R <- loc$w_R
  k <- loc$k
  invG <- loc$invG
  H.p <- loc$H
  invH.p <- loc$invH
  h.x <- loc$h.xy[1]
  h.y <- loc$h.xy[2]

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% sqrtw_Y
  resd <- eY - eR %*% (H.p %*% beta)

  if (vce=="hc0") {
    w.vce <- 1
  } else if (vce=="hc1") {
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_v2()")
    w.vce <- sqrt(eN/(eN - k))
  } else if (vce=="hc2") {
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_v2()")
    hii <- rd2d_leverage(sqrtw_R, invG)
    rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_multi_v2()")
    w.vce <- sqrt(1/(1-hii))
  } else if (vce == "hc3"){
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_v2()")
    hii <- rd2d_leverage(sqrtw_R, invG)
    rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_multi_v2()")
    w.vce <- 1/(1-hii)
  }

  resd <- sweep(resd, 1, w.vce, "*")

  n.outcomes <- 2
  cov.half.const <- NULL

  if (is.null(cluster)){
    cov.half.const <- matrix(NA_real_, nrow = eN, ncol = k * n.outcomes)
    for (j in 1:n.outcomes){
      sigma.half.const <- sweep(sqrtw_R, 1, sqrt(ew) * resd[,j], "*") *
        sqrt(h.x * h.y)
      cov.half.const[,((j - 1) * k + 1):(j * k)] <-
        sigma.half.const %*% invG
    }
  }

  if (!is.null(cluster)){
    n <- length(eC)
    g <- length(clusters)
    g.eff <- length(unique(eC))
    w.w <- if (cluster.df) {
      ((n-1)/(n-k))*(g.eff/(g.eff-1))
    } else {
      1
    }

    cov.half.const <- matrix(0, nrow = g, ncol = k * n.outcomes)

    for (j in 1:n.outcomes){
      j.idx <- ((j - 1) * k + 1):(j * k)
      sums <- rd2d_cluster_sums(w_R * resd[,j], eC, clusters) *
        sqrt(h.x * h.y * w.w)
      cov.half.const[,j.idx] <- sums %*% invG
    }
  }

  return(list("ind" = ind, "cov.half.const" = cov.half.const))
}

get_cov_half_multi_both_v2 <- function(dat, h, p, vce = "hc1",
                                       kernel = "epa", kernel_type = "prod",
                                       cluster = NULL, clusters = NULL){

  loc <- rd2d_local_design(
    dat, h, p, kernel, kernel_type, cluster, outcomes = c("y", "fuzzy")
  )

  ind <- loc$ind
  eN <- loc$eN
  ew <- loc$ew
  eY <- loc$eY
  eC <- loc$ecluster
  eR <- loc$eR
  sqrtw_R <- loc$sqrtw_R
  sqrtw_Y <- sweep(eY, 1, loc$sqrt.ew, "*")
  w_R <- loc$w_R
  k <- loc$k
  invG <- loc$invG
  H.p <- loc$H
  invH.p <- loc$invH
  h.x <- loc$h.xy[1]
  h.y <- loc$h.xy[2]

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% sqrtw_Y
  resd <- eY - eR %*% (H.p %*% beta)

  if (vce=="hc0") {
    w.vce <- 1
  } else if (vce=="hc1") {
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_both_v2()")
    w.vce <- sqrt(eN/(eN - k))
  } else if (vce=="hc2") {
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_both_v2()")
    hii <- rd2d_leverage(sqrtw_R, invG)
    rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_multi_both_v2()")
    w.vce <- sqrt(1/(1-hii))
  } else if (vce == "hc3"){
    rd2d_warn_vce_small_en(vce, eN, k, "in get_cov_half_multi_both_v2()")
    hii <- rd2d_leverage(sqrtw_R, invG)
    rd2d_warn_vce_high_leverage(vce, hii, "in get_cov_half_multi_both_v2()")
    w.vce <- 1/(1-hii)
  }

  resd.signed <- sweep(resd, 1, w.vce, "*")
  resd.abs <- abs(resd)
  resd.abs <- sweep(resd.abs, 1, w.vce, "*")

  n.outcomes <- 2
  signed.half.const <- abs.half.const <- NULL

  if (is.null(cluster)){
    signed.half.const <- matrix(NA_real_, nrow = eN, ncol = k * n.outcomes)
    abs.half.const <- matrix(NA_real_, nrow = eN, ncol = k * n.outcomes)

    for (j in 1:n.outcomes){
      j.idx <- ((j - 1) * k + 1):(j * k)
      sigma.signed <- sweep(
        sqrtw_R, 1, sqrt(ew) * resd.signed[,j], "*"
      ) * sqrt(h.x * h.y)
      sigma.abs <- sweep(
        sqrtw_R, 1, sqrt(ew) * resd.abs[,j], "*"
      ) * sqrt(h.x * h.y)

      signed.half.const[,j.idx] <- sigma.signed %*% invG
      abs.half.const[,j.idx] <- sigma.abs %*% invG
    }
  }

  if (!is.null(cluster)){
    n <- length(eC)
    g <- length(clusters)
    g.eff <- length(unique(eC))
    w.w <- ((n-1)/(n-k))*(g.eff/(g.eff-1))

    signed.half.const <- matrix(0, nrow = g, ncol = k * n.outcomes)
    abs.half.const <- matrix(0, nrow = g, ncol = k * n.outcomes)

    for (j in 1:n.outcomes){
      j.idx <- ((j - 1) * k + 1):(j * k)
      signed.half.const[,j.idx] <- rd2d_cluster_sums(
        w_R * resd.signed[,j], eC, clusters
      ) * sqrt(h.x * h.y * w.w)
      abs.half.const[,j.idx] <- rd2d_cluster_sums(
        w_R * resd.abs[,j], eC, clusters
      ) * sqrt(h.x * h.y * w.w)
      signed.half.const[,j.idx] <- signed.half.const[,j.idx] %*% invG
      abs.half.const[,j.idx] <- abs.half.const[,j.idx] %*% invG
    }
  }

  return(list(
    "ind" = ind,
    "signed.half.const" = signed.half.const,
    "abs.half.const" = abs.half.const
  ))
}

############################## Critical Value ##################################

rd2d_cov_to_corr <- function(cov, tol = 1e-10) {
  if (!is.matrix(cov) || nrow(cov) != ncol(cov)) {
    stop("cov must be a square matrix.", call. = FALSE)
  }
  if (!all(is.finite(cov)) || !all(is.finite(diag(cov))) ||
      any(diag(cov) <= 0)) {
    stop(
      "cov must have finite entries and strictly positive diagonal entries.",
      call. = FALSE
    )
  }

  se <- sqrt(diag(cov))
  corr <- cov / tcrossprod(se)
  corr <- (corr + t(corr)) / 2
  diag(corr) <- 1

  eig <- eigen(corr, symmetric = TRUE)
  if (!all(is.finite(eig$values))) {
    stop("covariance correlation matrix has non-finite eigenvalues.",
         call. = FALSE)
  }

  min.eig <- min(eig$values)
  if (min.eig < tol) {
    eig.values <- pmax(eig$values, tol)
    corr.reg <- eig$vectors %*% diag(eig.values, nrow = length(eig.values)) %*%
      t(eig$vectors)
    corr.reg <- (corr.reg + t(corr.reg)) / 2
    corr.reg <- corr.reg / sqrt(tcrossprod(diag(corr.reg)))
    corr.reg <- (corr.reg + t(corr.reg)) / 2
    diag(corr.reg) <- 1

    max.adjust <- max(abs(corr.reg - corr))
    if (min.eig < -sqrt(tol) || max.adjust > sqrt(tol)) {
      warning(
        paste(
          "Estimated covariance correlation matrix was not positive",
          "semidefinite; using an eigenvalue-regularized approximation."
        ),
        call. = FALSE
      )
    }
    corr <- corr.reg
  }

  corr
}

rd2d_row_max <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (ncol(x) == 1) return(x[, 1])
  x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]
}

rd2d_validate_repp <- function(repp) {
  if (!is.numeric(repp) || length(repp) != 1 || !is.finite(repp) ||
      repp < 1 || repp != as.integer(repp)) {
    stop("repp must be a positive integer", call. = FALSE)
  }
  as.integer(repp)
}

rd2d_cval <- function(cov, rep, side="two", alpha, lp=Inf) {
  tvec <- c()

  cval <- NA
  m <- dim(cov)[1]
  cov.t <- rd2d_cov_to_corr(cov)

  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)

  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- rd2d_row_max(abs(sim))}
      else {
        tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})
      }
    } else if (side == "left") {
      tvec <- rd2d_row_max(sim)
    } else if (side == "right") {
      tvec <- rd2d_row_max(-sim)
    }
  }

  if (!is.null(side)) {
    cval <- quantile(tvec, alpha/100, na.rm = TRUE, names = FALSE, type = 2)
  }
  return(cval)
}

############################## p Value ##################################

rd2d_pval <- function(tstat, cov, rep, side="two", lp=Inf) {
  tvec <- c()

  pval <- NA
  m <- dim(cov)[1]
  cov.t <- rd2d_cov_to_corr(cov)

  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)

  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- rd2d_row_max(abs(sim))}
      else {
        tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})
      }
    } else if (side == "left") {
      tvec <- rd2d_row_max(sim)
    } else if (side == "right") {
      tvec <- rd2d_row_max(-sim)
    }
  }

  if (!is.null(side)) {
    pval <- mean(tvec >= abs(tstat))
  }
  return(pval)
}

############################## Confidence Bands ################################

rd2d_cb <- function(mu.hat, cov.us, rep, side, alpha){

  # mu.hat: estimated (derivatives) of treatment effect
  # cov.us: estimated covariance matrix for effects at all evaluation points
  # rep: number of repetitions for Gaussian simulation
  # side: "pos", "neg", or "two"
  # alpha: confidence level
  # If side == "two", returns upper and lower bounds of CI and CB
  # If side == "pos", returns upper bounds of CI and CB
  # If side == "neg", returns lower bounds of CI and CB

  se.hat <- sqrt(diag(cov.us))
  cval <- rd2d_cval(cov.us, rep = rep, side=side, alpha = alpha, lp=Inf)

  if (side == "two"){
    zval <- qnorm((alpha + 100)/ 200)
    CI.l <- mu.hat - zval * se.hat; CI.r <- mu.hat + zval * se.hat
    CB.l <- mu.hat - cval * se.hat; CB.r <- mu.hat + cval * se.hat
  }

  if (side == "left"){
    zval <- qnorm(alpha / 100)
    CI.r <- mu.hat + zval * se.hat; CI.l <- rep(-Inf, length(CI.r))
    CB.r <- mu.hat + cval * se.hat; CB.l <- rep(-Inf, length(CB.r))
  }
  if (side == "right"){
    zval <- qnorm(alpha / 100)
    CI.l <- mu.hat - zval * se.hat; CI.r <- rep(Inf, length(CI.l))
    CB.l <- mu.hat - cval * se.hat; CB.r <- rep(Inf, length(CB.l))
  }

  return(list(CI.l = CI.l, CI.r = CI.r, CB.l = CB.l, CB.r = CB.r))
}

############################# Get Basis ########################################

get_basis_xy <- function(u.x.1, u.x.2, p) {
  n <- length(u.x.1)
  if (p == 0) return(matrix(1, nrow = n, ncol = 1))
  if (p == 1) return(cbind(1, u.x.1, u.x.2))
  if (p == 2) {
    return(cbind(1, u.x.1, u.x.2, u.x.1^2, u.x.1 * u.x.2, u.x.2^2))
  }
  if (p == 3) {
    return(cbind(
      1, u.x.1, u.x.2,
      u.x.1^2, u.x.1 * u.x.2, u.x.2^2,
      u.x.1^3, u.x.1^2 * u.x.2, u.x.1 * u.x.2^2, u.x.2^3
    ))
  }

  result <- matrix(
    NA_real_,
    nrow = n,
    ncol = factorial(p+2)/(factorial(p) * 2)
  )
  result[,1] <- 1
  count <- 2
  for (j in 1:p){
    for (k in 0:j){
      result[,count] <- u.x.1^(j-k) * u.x.2^k
      count <- count + 1
    }
  }
  result
}

get_basis <- function(u,p){
  get_basis_xy(u[,1], u[,2], p)
}

############################### Get H ##########################################

# old version with one h

# get_H <- function(h, p){
#   diags <- c()
#   for (i in 0:p){
#     diags <- c(diags, rep(h^i, i+1))
#   }
#   result <- diag(diags)
#   return(result)
# }

# new version with on h for each coordinate
get_H <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1

  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- h.x.1^(j-k) * h.x.2^k
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

########################### Get Inverse of H ###################################

# old version with one h

# get_invH <- function(h, p){
#   diags <- c()
#   for (i in 0:p){
#     diags <- c(diags, rep(h^(-i), i+1))
#   }
#   result <- diag(diags)
#   return(result)
# }

# new version with one h for each dimension
get_invH <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1

  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- 1/(h.x.1^(j-k) * h.x.2^k)
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

############################## lm fit ##########################################

# kernel_type = "prod" or "rad"

rd2d_lm <- function(dat, h, p, vce = "hc1", kernel = "epa",
                    kernel_type = "prod",
                    cluster = NULL, varr = FALSE){

  loc <- rd2d_local_design(dat, h, p, kernel, kernel_type, cluster, outcomes = "y")

  eN <- loc$eN
  eY <- loc$eY[, 1]
  eC <- loc$ecluster
  eR <- loc$eR
  sqrtw_R <- loc$sqrtw_R
  sqrtw_Y <- loc$sqrt.ew * eY
  w_R <- loc$w_R
  k <- loc$k
  invG <- loc$invG
  invH.p <- loc$invH
  H.p <- loc$H
  h.x <- loc$h.xy[1]
  h.y <- loc$h.xy[2]

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)

  cov.const <- NA

  if (varr){

    resd <- (eY - (eR %*% H.p) %*% beta)[,1]

    if (vce=="hc0") {
      w.vce = 1
    } else if (vce=="hc1") {
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm()")
      w.vce = sqrt(eN/(eN - k))
    } else if (vce=="hc2") {
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in rd2d_lm()")
      w.vce = sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in rd2d_lm()")
      w.vce = 1/(1-hii)
    }

    resd <- resd * w.vce

    sigma <- rd2d_vce(w_R, resd, eC, c(h.x, h.y))

    # sigma <- t(resd * as.matrix(sqrtw_R)) %*%
    #   (ew * resd * as.matrix(sqrtw_R)) * h^2

    cov.const <- crossprod(invG, sigma %*% invG)
  }

  return(list("beta" = beta, "cov.const" = cov.const, "eN" = eN))
}

rd2d_lm_multi <- function(dat, h, p, vce = "hc1", kernel = "epa",
                          kernel_type = "prod", cluster = NULL, varr = FALSE,
                          outcomes = c("y", "fuzzy")){

  loc <- rd2d_local_design(dat, h, p, kernel, kernel_type, cluster, outcomes)

  eN <- loc$eN
  eY <- loc$eY
  eC <- loc$ecluster
  eR <- loc$eR
  sqrtw_R <- loc$sqrtw_R
  sqrtw_Y <- sweep(eY, 1, loc$sqrt.ew, "*")
  w_R <- loc$w_R
  k <- loc$k
  invG <- loc$invG
  invH.p <- loc$invH
  H.p <- loc$H
  h.x <- loc$h.xy[1]
  h.y <- loc$h.xy[2]

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% sqrtw_Y
  colnames(beta) <- outcomes

  cov.const <- vector("list", length(outcomes))
  names(cov.const) <- outcomes

  if (varr){
    resd <- eY - eR %*% (H.p %*% beta)

    if (vce=="hc0") {
      w.vce = 1
    } else if (vce=="hc1") {
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm_multi()")
      w.vce = sqrt(eN/(eN - k))
    } else if (vce=="hc2") {
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm_multi()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in rd2d_lm_multi()")
      w.vce = sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      rd2d_warn_vce_small_en(vce, eN, k, "in rd2d_lm_multi()")
      hii <- rd2d_leverage(sqrtw_R, invG)
      rd2d_warn_vce_high_leverage(vce, hii, "in rd2d_lm_multi()")
      w.vce = 1/(1-hii)
    }

    resd <- sweep(resd, 1, w.vce, "*")

    for (outcome in outcomes) {
      cov.const[[outcome]] <- rd2d_vce(
        w_R, resd[, outcome], eC, c(h.x, h.y)
      )
      cov.const[[outcome]] <- crossprod(
        invG, cov.const[[outcome]] %*% invG
      )
    }
  }

  list("beta" = beta, "cov.const" = cov.const, "eN" = eN)
}

###### Get coefficients of a linear combination of (p+1)-th derivatives ########

get_coeff <- function(dat.centered,vec, p,dn, kernel, kernel_type){

  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "distance")]
  if (kernel_type == "prod"){
    w.v <- kernel_weight(dat.centered$x.1/c(dn), kernel) *
      kernel_weight(dat.centered$x.2/c(dn), kernel) /
      c(dn * dn)
  }
  else{
    w.v <- kernel_weight(dat.centered$distance/c(dn), kernel)/c(dn^2)
  }

  # w.v <- kernel_weight(dat.centered$distance/c(dn), kernel)/c(dn^2)

  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)

  ew.v <- w.v[ind.v]
  eY.v <- dat.centered$y[ind.v]

  eR.v.aug <- get_basis_xy(
    dat.centered$x.1[ind.v] / dn,
    dat.centered$x.2[ind.v] / dn,
    p + 1
  )
  eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2)), drop = FALSE]
  p.count <- factorial(p+2)/(factorial(p)*2)
  p1.count <- factorial(p+1+2)/(factorial(p+1)*2)
  eS.v <- eR.v.aug[,(p.count + 1):p1.count, drop = FALSE]

  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v

  invG.v <- qrXXinv(sqrtw_R.v)

  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)

  return(vec.q)
}


###################### Robust Covariance Estimation ############################

# old version with one h
# rd2d_vce <- function(w_R, resd, eC, h){
#   n <- length(eC)
#   k <- dim(w_R)[2]
#   M <- matrix(0, nrow = k, ncol = k)
#   if (is.null(eC)){
#     w.w <- 1
#     M <-  crossprod(resd * as.matrix(w_R)) * h^2
#   }
#   else{
#     clusters = unique(eC)
#     g     = length(clusters)
#     w.w =((n-1)/(n-k))*(g/(g-1))
#     for (i in 1:g) {
#       ind=ecluster==clusters[i]
#       w_R_i = w_R[ind,,drop=FALSE]
#       resd_i = resd[ind]
#       resd_i <- matrix(resd_i, ncol = 1)
#       w_R_resd_i = t(crossprod(w_R_i,resd_i))
#       M = M + crossprod(w_R_resd_i,w_R_resd_i) * h^2
#     }
#   }
#   return(M * w.w)
# }

# new version with one h for each coordinate
rd2d_vce <- function(w_R, resd, eC, h){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }

  n <- length(eC)
  k <- dim(w_R)[2]
  M <- matrix(0, nrow = k, ncol = k)
  if (is.null(eC)){
    w.w <- 1
    M <-  crossprod(resd * as.matrix(w_R)) * h.x.1 * h.x.2
  }
  else{
    clusters = unique(eC)
    g     = length(clusters)
    w.w =((n-1)/(n-k))*(g/(g-1))
    scores <- rd2d_cluster_sums(w_R * resd, eC, clusters)
    M <- crossprod(scores) * h.x.1 * h.x.2
  }
  return(M * w.w)
}

############################## rd2d unique #####################################

rd2d_unique <- function(dat){

  dat <- dat[,c("x.1", "x.2", "y", "d")]

  ord <- order(dat[,1], dat[,2])

  dat <- dat[ord,]

  N <- dim(dat)[1]

  # if x has one or no element
  if (N == 0) return(list(unique = NULL, freq = c(), index = c()))
  if (N == 1) return(list(unique = dat, freq = 1, index = 1))

  # else
  uniqueIndex <- c(
    c(dat[2:N, 1] != dat[1:(N-1),1] |
      dat[2:N, 2] != dat[1:(N-1),2]),
    TRUE
  )
  unique <- dat[uniqueIndex,]
  nUnique <- dim(unique)[1]

  # all are distinct
  if (nUnique == N) return(list(unique=unique, freq=rep(1,N), index=1:N))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=N, index=N))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, index=(1:N)[uniqueIndex]))
}

############################# output formatting ################################

# TODO: merge kernel and kernel.type

# TODO: merge cluster and vce
