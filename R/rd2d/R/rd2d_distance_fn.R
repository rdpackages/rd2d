
############################ Rule of Thumb #####################################

# Use Scott's Rule: n^{-1/(d+4)} * sample variance

rdbw2d_distance_rot <- function(distance, kernel.type){

  mu2K.squared <- NA
  l2K.squared <- NA

  if (kernel.type == "Epanechnikov"){mu2K.squared <- 1/6; l2K.squared <- 4/(3 * pi)}
  if (kernel.type == "Triangular"){mu2K.squared <- 3/20; l2K.squared <- 3/(2 * pi)}
  if (kernel.type == "Uniform"){mu2K.squared <- 1/4; l2K.squared <- 1/(pi)}
  if (kernel.type == "Gaussian"){mu2K.squared <- 1; l2K.squared <- 1/(4 * pi)}

  # Estimate sample variance.

  N <- length(distance)

  var.hat <- 1 / 2 * mean(distance^2)
  dim.x <- 2
  trace.const <- 1 / (2 * pi * var.hat^3)

  hROT <- ((dim.x * l2K.squared) /
    (N * mu2K.squared * trace.const))^(1 / (4 + dim.x))

  return(hROT)
}

rd2d_distance_poly_x <- function(distance, degree) {
  outer(distance, 0:degree, `^`)
}

rd2d_distance_poly_lm_fallback <- function(Y, distance, degree) {
  dat <- data.frame(y = Y, distance = distance)
  if (degree >= 2) {
    for (j in 2:degree) {
      dat[[paste0("V", j - 1)]] <- distance^j
    }
  }

  rhs <- if (degree == 0) "1" else "distance"
  if (degree >= 2) {
    rhs <- paste(rhs, paste0("V", seq_len(degree - 1)), sep = " + ")
  }
  fit <- lm(as.formula(paste("y ~", rhs)), data = dat)
  coef <- as.numeric(fit$coefficients)
  se <- rep(NA_real_, degree + 1)
  sum.fit <- summary(fit)$coefficients
  se[match(rownames(sum.fit), names(fit$coefficients))] <- sum.fit[, 2]

  list(
    coef = coef,
    se = se,
    predict = function(newdistance) as.vector(
      rd2d_distance_poly_x(newdistance, degree) %*% coef
    )
  )
}

rd2d_distance_poly_lm <- function(Y, distance, degree) {
  X <- rd2d_distance_poly_x(distance, degree)
  k <- ncol(X)

  if (nrow(X) <= k) {
    return(rd2d_distance_poly_lm_fallback(Y, distance, degree))
  }

  invXX <- tryCatch(
    qrXXinv(X),
    error = function(e) NULL
  )
  if (is.null(invXX) || any(!is.finite(invXX))) {
    return(rd2d_distance_poly_lm_fallback(Y, distance, degree))
  }

  coef <- as.vector(invXX %*% crossprod(X, Y))
  if (anyNA(coef) || any(!is.finite(coef))) {
    return(rd2d_distance_poly_lm_fallback(Y, distance, degree))
  }

  res <- Y - as.vector(X %*% coef)
  df.residual <- nrow(X) - k
  sigma2 <- sum(res^2) / df.residual
  if (!is.finite(sigma2)) {
    return(rd2d_distance_poly_lm_fallback(Y, distance, degree))
  }
  se <- sqrt(pmax(0, diag(invXX) * sigma2))

  list(
    coef = coef,
    se = se,
    predict = function(newdistance) as.vector(
      rd2d_distance_poly_x(newdistance, degree) %*% coef
    )
  )
}

############# One Step Bandwidth Selection for Distance Based Estimator ################

rdbw2d_distance_local_intercepts_multi <- function(Y, fuzzy, distance, h, p, kernel) {
  w <- kernel_weight(distance / h, kernel) / h^2
  ind <- w > 0
  if (sum(ind) <= p + 1) return(c(y = NA_real_, fuzzy = NA_real_))

  R <- outer(distance[ind] / h, 0:p, `^`)
  sqrtw_R <- sqrt(w[ind]) * as.matrix(R)
  inv.gamma <- tryCatch(qrXXinv(sqrtw_R), error = function(e) NULL)
  if (is.null(inv.gamma) || any(!is.finite(inv.gamma))) {
    return(c(y = NA_real_, fuzzy = NA_real_))
  }

  XtWY <- cbind(
    y = crossprod(R, w[ind] * Y[ind]),
    fuzzy = crossprod(R, w[ind] * fuzzy[ind])
  )
  beta <- inv.gamma %*% XtWY
  c(y = beta[1, 1], fuzzy = beta[1, 2])
}

rd2d_distance_covariate_gamma <- function(Y, Z, distance, h.0, h.1, p,
                                          kernel, ginv.tol = 1e-20,
                                          covs.drop = NULL, covs.tol = NULL) {
  if (is.null(Z) || ncol(as.matrix(Z)) == 0) return(NULL)

  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  cov.names <- colnames(Z)
  if (is.null(cov.names)) cov.names <- paste0("z.", seq_len(ncol(Z)))
  colnames(Z) <- cov.names
  side.0 <- distance < 0
  side.1 <- distance >= 0
  abs.distance <- abs(distance)

  side_components <- function(side, h.side) {
    u <- abs.distance[side]
    Y.side <- Y[side, , drop = FALSE]
    Z.side <- Z[side, , drop = FALSE]
    w <- kernel_weight(u / h.side, kernel) / h.side^2
    ind <- w > 0
    u <- u[ind]
    Y.side <- Y.side[ind, , drop = FALSE]
    Z.side <- Z.side[ind, , drop = FALSE]
    ew <- w[ind]
    R <- outer(u / h.side, 0:p, `^`)
    w_R <- ew * R
    sqrtw_R <- sqrt(ew) * R
    invG <- qrXXinv(sqrtw_R)
    U.y <- crossprod(w_R, Y.side)
    U.z <- crossprod(w_R, Z.side)

    list(
      ZWZ = crossprod(Z.side * ew, Z.side) -
        crossprod(U.z, invG %*% U.z),
      ZWY = crossprod(Z.side * ew, Y.side) -
        crossprod(U.z, invG %*% U.y)
    )
  }

  comp.0 <- side_components(side.0, h.0)
  comp.1 <- side_components(side.1, h.1)
  ZWZ <- comp.0$ZWZ + comp.1$ZWZ
  ZWY <- comp.0$ZWY + comp.1$ZWY

  gamma <- rd2d_covs_rank_solve(
    ZWZ, ZWY, cov.names, covs.drop = covs.drop, covs.tol = covs.tol,
    ginv.tol = ginv.tol
  )
  as.matrix(gamma)
}

rd2d_distance_apply_covariate_adjustment <- function(Y, Z, gamma) {
  if (is.null(Z) || is.null(gamma)) return(Y)
  Y <- as.matrix(Y)
  adjusted <- Y - as.matrix(Z) %*% gamma
  colnames(adjusted) <- colnames(Y)
  adjusted
}

rd2d_distance_vce <- function(w_R, resd, eC, h, k.df = NULL,
                              cluster.df = TRUE, clusters = NULL) {
  w_R <- as.matrix(w_R)
  if (is.null(k.df)) k.df <- ncol(w_R)
  if (is.null(eC)) {
    return(crossprod(resd * w_R) * h^2)
  }

  if (is.null(clusters)) clusters <- unique(eC)
  n <- length(eC)
  g <- length(unique(eC))
  scores <- rd2d_cluster_sums(w_R * resd, eC, clusters)
  scale <- if (cluster.df) ((n - 1) / (n - k.df)) * (g / (g - 1)) else 1
  crossprod(scores) * h^2 * scale
}

rd2d_distance_bwcheck_limits <- function(distance, bwcheck) {
  n <- length(distance)
  if (n == 0) return(c(min = NA_real_, max = NA_real_))
  k <- min(as.integer(bwcheck), n)
  kth <- if (k == n) max(distance) else sort.int(distance, partial = k)[k]
  c(min = kth, max = max(distance))
}

rd2d_distance_unique_counts <- function(distance) {
  n <- length(distance)
  if (n == 0) return(c(M = 0, M.0 = 0, M.1 = 0))

  if (anyDuplicated(distance) == 0) {
    m.0 <- sum(distance < 0)
    m.1 <- n - m.0
  } else {
    unique.distance <- unique(distance)
    m.0 <- sum(unique.distance < 0)
    m.1 <- length(unique.distance) - m.0
  }

  c(M = m.0 + m.1, M.0 = m.0, M.1 = m.1)
}

rdbw2d_distance_bw <- function(Y, distance, p = 1, kernel, target.vec = NULL,
                  rot = NULL, vce = "hc0", cluster = NULL,
                  bwcheck = 50 + p + 1, scaleregul = 1, cqt = 0.5,
                  verbose = FALSE, fuzzy = NULL, bwparam = "main",
                  fitmethod = c("separate", "joint"), covs.eff = NULL){

  fitmethod <- match.arg(fitmethod)

  # Data Cleaning

  distance_mat <- distance
  neval <- ncol(distance_mat)
  if (is.null(target.vec)) {
    target.vec <- matrix(0, nrow = neval, ncol = p + 1)
    target.vec[, 1] <- 1
  }
  N <- length(Y)
  if (!is.null(cluster)){
    M <- length(unique(cluster))
  } else {
    M <- N
  }
  is.fuzzy <- !is.null(fuzzy)
  warned.weak <- FALSE
  joint <- fitmethod == "joint"
  covariate.adjusted <- !is.null(covs.eff)
  covs.eff <- if (covariate.adjusted) as.matrix(covs.eff) else NULL
  n.cov <- if (covariate.adjusted) ncol(covs.eff) else 0
  clustered <- !is.null(cluster)
  clusters <- unique(cluster)

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Loop over points of evaluations

  results <- data.frame(matrix(NA, ncol = 15, nrow = neval))
  colnames(results) <- c('h.0', 'h.1', 'b.0', 'b.1', 'v.0', 'v.1',
                         'v.01', 'r.0', 'r.1', 'N.Co', 'N.Tr',
                         'bw.min.0', 'bw.min.1', 'bw.max.0', 'bw.max.1')

  for (i in 1:neval){

    distance <- distance_mat[,i]
    n.cov.i <- n.cov
    d <- (distance >= 0)
    distance <- abs(distance)
    side.0 <- !d
    side.1 <- d

    if (!is.null(rot)){
      dn <- rot[i]
    } else {
      dn <- rdbw2d_distance_rot(distance, kernel.type)
    }

    vec <- target.vec[i,]

    distance.0 <- distance[side.0]
    distance.1 <- distance[side.1]
    y.0 <- Y[side.0]
    y.1 <- Y[side.1]
    fuzzy.0 <- if (is.fuzzy) fuzzy[side.0] else NULL
    fuzzy.1 <- if (is.fuzzy) fuzzy[side.1] else NULL
    covs.0 <- if (covariate.adjusted) covs.eff[side.0, , drop = FALSE] else NULL
    covs.1 <- if (covariate.adjusted) covs.eff[side.1, , drop = FALSE] else NULL
    cluster.0 <- cluster[side.0]
    cluster.1 <- cluster[side.1]

    # Take effective sample
    dn.0 <- dn; dn.1 <- dn

    bw.min.0 <- NA
    bw.min.1 <- NA
    bw.max.0 <- NA
    bw.max.1 <- NA

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      bounds.0 <- rd2d_distance_bwcheck_limits(distance.0, bwcheck)
      bounds.1 <- rd2d_distance_bwcheck_limits(distance.1, bwcheck)
      bw.min.0   <- bounds.0["min"]
      bw.min.1   <- bounds.1["min"]
      bw.max.0   <- bounds.0["max"]
      bw.max.1   <- bounds.1["max"]
      dn.0     <- max(dn, bw.min.0)
      dn.1     <- max(dn, bw.min.1)
      dn.0     <- min(dn.0, bw.max.0)
      dn.1     <- min(dn.1, bw.max.1)
    }

    if (covariate.adjusted) {
      Y.gamma <- if (is.fuzzy) {
        cbind(y = c(y.0, y.1), fuzzy = c(fuzzy.0, fuzzy.1))
      } else {
        cbind(y = c(y.0, y.1))
      }
      gamma <- rd2d_distance_covariate_gamma(
        Y.gamma, rbind(covs.0, covs.1), c(-distance.0, distance.1),
        dn.0, dn.1, p, kernel
      )
      n.cov.i <- rd2d_covs_gamma_rank(gamma, colnames(covs.eff))
      Y.adj <- rd2d_distance_apply_covariate_adjustment(
        Y.gamma, rbind(covs.0, covs.1), gamma
      )
      n0 <- length(y.0)
      y.0 <- Y.adj[seq_len(n0), "y"]
      y.1 <- Y.adj[n0 + seq_len(length(y.1)), "y"]
      if (is.fuzzy) {
        fuzzy.0 <- Y.adj[seq_len(n0), "fuzzy"]
        fuzzy.1 <- Y.adj[n0 + seq_len(length(fuzzy.1)), "fuzzy"]
      }
    }

    if (is.fuzzy && bwparam == "main") {
      fit.grad.0 <- rdbw2d_distance_local_intercepts_multi(
        y.0, fuzzy.0, distance.0, dn.0, p, kernel
      )
      fit.grad.1 <- rdbw2d_distance_local_intercepts_multi(
        y.1, fuzzy.1, distance.1, dn.1, p, kernel
      )
      tau.itt.grad <- fit.grad.1["y"] - fit.grad.0["y"]
      tau.fs.grad <- fit.grad.1["fuzzy"] - fit.grad.0["fuzzy"]
      if (is.finite(tau.fs.grad) &&
          abs(tau.fs.grad) > sqrt(.Machine$double.eps)) {
        grad.itt <- 1 / tau.fs.grad
        grad.fs <- -tau.itt.grad / tau.fs.grad^2
        y.0 <- grad.itt * y.0 + grad.fs * fuzzy.0
        y.1 <- grad.itt * y.1 + grad.fs * fuzzy.1
      } else if (!warned.weak) {
        warning(
          "Weak or zero first-stage fuzzy RD estimate detected in bandwidth selection; using reduced-form outcome bandwidth."
        )
        warned.weak <- TRUE
      }
    }

    # Estimate (p+1)-th derivatives using cqt fraction of the data.
    degree.deriv <- p + 1
    thrshd.0 <- quantile(distance.0, cqt)
    thrshd.1 <- quantile(distance.1, cqt)
    fit.ind.0 <- distance.0 <= thrshd.0
    fit.ind.1 <- distance.1 <= thrshd.1

    model.deriv.0 <- rd2d_distance_poly_lm(y.0[fit.ind.0], distance.0[fit.ind.0], degree.deriv)
    model.deriv.1 <- rd2d_distance_poly_lm(y.1[fit.ind.1], distance.1[fit.ind.1], degree.deriv)

    deriv.ppls1.0 <- model.deriv.0$coef[length(model.deriv.0$coef)]
    deriv.ppls1.1 <- model.deriv.1$coef[length(model.deriv.1$coef)]

    if (verbose){
      print("deriv.ppls1.0"); print(deriv.ppls1.0)
      print("deriv.ppls1.1"); print(deriv.ppls1.1)
    }

    w.0   <- kernel_weight(distance.0/dn.0, kernel)/(dn.0^2)
    w.1   <- kernel_weight(distance.1/dn.1, kernel)/(dn.1^2)
    ind.0 <- (w.0 > 0); ind.1 <- (w.1 > 0)
    eN.0 <- sum(ind.0); eN.1 <- sum(ind.1)

    eX.0  <- distance.0[ind.0]
    eX.1  <- distance.1[ind.1]
    eY.0 <- y.0[ind.0]
    eY.1 <- y.1[ind.1]
    eweight.0 <- w.0[ind.0]
    eweight.1 <- w.1[ind.1]
    u.0   <- eX.0
    u.1   <- eX.1

    # -- Sd from residuals of local polynomial fit

    sd.0 <- eY.0 - model.deriv.0$predict(eX.0)
    sd.1 <- eY.1 - model.deriv.1$predict(eX.1)

    # Compute bread, meat and half-bread matrices for p-th degree model

    R.0.p <- as.matrix(get_basis_v1(u.0/dn.0,p))
    R.1.p <- as.matrix(get_basis_v1(u.1/dn.1,p))
    sqrtw_R.0 <- sqrt(eweight.0) * as.matrix(R.0.p)
    sqrtw_R.1 <- sqrt(eweight.1) * as.matrix(R.1.p)
    inv.gamma0.p <- qrXXinv(sqrtw_R.0) # Bread matrix
    inv.gamma1.p <- qrXXinv(sqrtw_R.1) # Bread matrix
    sigma0.half.p <- eweight.0 * as.matrix(R.0.p)
    sigma1.half.p <- eweight.1 * as.matrix(R.1.p) # Half-bread matrix

    # vce methods

    if (vce=="hc0") {
      w.vce.0 = 1
      w.vce.1 = 1
    } else if (vce=="hc1") {
      if (clustered && (joint || covariate.adjusted)) {
        w.vce.0 = 1
        w.vce.1 = 1
      } else if (joint) {
        w.vce.0 = sqrt((eN.0 + eN.1) / (eN.0 + eN.1 - 2 * (p + 1) - n.cov.i))
        w.vce.1 = w.vce.0
      } else {
        w.vce.0 = sqrt(eN.0/(eN.0-p-1-n.cov.i))
        w.vce.1 = sqrt(eN.1/(eN.1-p-1-n.cov.i))
      }
    } else if (vce=="hc2") {
      hii.0 <- rowSums((sqrtw_R.0 %*% inv.gamma0.p) * sqrtw_R.0)
      hii.1 <- rowSums((sqrtw_R.1 %*% inv.gamma1.p) * sqrtw_R.1)
      w.vce.0 = sqrt(1/(1-hii.0))
      w.vce.1 = sqrt(1/(1-hii.1))
    } else if (vce=="hc3"){
      hii.0 <- rowSums((sqrtw_R.0 %*% inv.gamma0.p) * sqrtw_R.0)
      hii.1 <- rowSums((sqrtw_R.1 %*% inv.gamma1.p) * sqrtw_R.1)
      w.vce.0 = 1/(1-hii.0)
      w.vce.1 = 1/(1-hii.1)
    }
    sd.0 <- sd.0 * w.vce.0
    sd.1 <- sd.1 * w.vce.1

    # (cluster-robust) meat matrix estimation

    eC.0 <- cluster.0[ind.0]
    eC.1 <- cluster.1[ind.1]
    cluster.df <- !(joint && clustered) && !(covariate.adjusted && clustered)
    k.side <- p + 1 + n.cov.i
    sigma.0.p <- rd2d_distance_vce(
      sigma0.half.p, sd.0, eC.0, dn.0, k.side, cluster.df, clusters
    )
    sigma.1.p <- rd2d_distance_vce(
      sigma1.half.p, sd.1, eC.1, dn.1, k.side, cluster.df, clusters
    )

    # sigma.0.p <-  t(sd.0 * sigma0.half.p) %*% (sd.0 * sigma0.half.p) * eN.0 / (eN.0 - (p+1)) * dn^2 # Meat matrices
    # sigma.1.p <-  t(sd.1 * sigma1.half.p) %*% (sd.1 * sigma1.half.p) * eN.1 / (eN.1 - (p+1)) * dn^2 # Meat matrices

    # Compute the coefficients for linear combination of (p+1)-th derivatives

    pmatrix.0 <- matrix(NA,nrow = eN.0, ncol = 1)
    pmatrix.1 <- matrix(NA, nrow = eN.1, ncol = 1)
    pmatrix.0[,1] <- (eX.0/dn.0)^(p+1) * eweight.0
    pmatrix.1[,1] <- (eX.1/dn.1)^(p+1) * eweight.1
    pmatrix.0 <- t(R.0.p) %*% pmatrix.0
    pmatrix.1 <- t(R.1.p) %*% pmatrix.1
    dmm <- p+2
    coeff.0 <- rep(0, dmm)
    coeff.1 <- rep(0, dmm)
    coeff.0[dmm] <- as.vector( t(as.matrix(vec)) %*% inv.gamma0.p %*% pmatrix.0 )
    coeff.1[dmm] <- as.vector( t(as.matrix(vec)) %*% inv.gamma1.p %*% pmatrix.1 )

    # Compute bias constants for lp using p-th order model

    B.p.0 <- coeff.0 %*% model.deriv.0$coef
    B.p.1 <- coeff.1 %*% model.deriv.1$coef

    # Compute regularization terms for lp using (p+1)-th order model

    R.0 <- coeff.0[dmm]^2 * model.deriv.0$se[dmm]^2
    R.1 <- coeff.1[dmm]^2 * model.deriv.1$se[dmm]^2

    if (verbose){
      print("model.deriv.0$coef = ")
      print(model.deriv.0$coef)
      print("model.deriv.1$coef = ")
      print(model.deriv.1$coef)
    }

    # Compute variance constants for lp using p-th order model

    # -- inv.gamma0 and sigma0.half using dn.0

    V.half.0 <- V.half.1 <- NULL
    if (clustered && joint) {
      scores.0 <- rd2d_cluster_sums(sigma0.half.p * sd.0, eC.0, clusters)
      scores.1 <- rd2d_cluster_sums(sigma1.half.p * sd.1, eC.1, clusters)
      V.half.0 <- as.vector(dn.0 * scores.0 %*% inv.gamma0.p %*% as.matrix(vec))
      V.half.1 <- as.vector(dn.1 * scores.1 %*% inv.gamma1.p %*% as.matrix(vec))
    }

    V.p.0 <- t(as.matrix(vec)) %*% inv.gamma0.p %*% sigma.0.p %*% inv.gamma0.p %*% as.matrix(vec)
    V.p.1 <- t(as.matrix(vec)) %*% inv.gamma1.p %*% sigma.1.p %*% inv.gamma1.p %*% as.matrix(vec)
    V.p.01 <- 0
    if ((joint || covariate.adjusted) && clustered) {
      if (joint) {
        scale <- rd2d_joint_scale(
          vce, eN.0 + eN.1, 2 * (p + 1) + n.cov.i,
          c(eC.0, eC.1), clustered
        )
        V.p.0 <- V.p.0 * scale^2
        V.p.1 <- V.p.1 * scale^2
        V.half.0 <- V.half.0 * scale
        V.half.1 <- V.half.1 * scale
      } else {
        scale.0 <- rd2d_joint_scale(vce, eN.0, p + 1 + n.cov.i, eC.0, clustered)
        scale.1 <- rd2d_joint_scale(vce, eN.1, p + 1 + n.cov.i, eC.1, clustered)
        V.p.0 <- V.p.0 * scale.0^2
        V.p.1 <- V.p.1 * scale.1^2
      }
    }
    if (clustered && joint) {
      V.p.01 <- as.numeric(crossprod(V.half.0, V.half.1))
    }

    # Optimal bandwidth for treatment effect using p-th order model

    V.diff <- V.p.0 + V.p.1 - 2 * V.p.01
    hn <- (2 * V.diff / ( (2 * p + 2) * ((B.p.0 -  B.p.1)^2 + scaleregul * R.0 + scaleregul * R.1)) )^(1/(2 * p + 4))

    if (verbose){
      print(paste("V.p.0 = ", V.p.0, ", V.p.1 = ", V.p.1, sep = ""))
      print("coeff.0 = "); print(coeff.0)
      print("coeff.1 = "); print(coeff.1)
    }
    results[i,] <- c(hn, hn, B.p.0, B.p.1, V.p.0, V.p.1, V.p.01,
                     R.0, R.1, eN.0, eN.1,
                     bw.min.0, bw.min.1, bw.max.0, bw.max.1)
  }

  return(results)
}

############################### rd2d_distance_fit ##################################

rd2d_distance_fit <- function(Y, distance, h, p, b, kernel, vce, bwcheck,
                              masspoints, cluster, cbands = TRUE,
                              masspoint.counts = NULL,
                              fitmethod = c("separate", "joint"),
                              covs.eff = NULL){

  fitmethod <- match.arg(fitmethod)

  distance_mat <- distance
  neval <- ncol(distance_mat)

  if (!is.null(b)){
    eval <- as.data.frame(b)
  } else {
    eval <- matrix(NA, nrow = neval, ncol = 2)
    eval <- as.data.frame(eval)
  }
  colnames(eval) <- c("x.1", "x.2")

  hgrid <- h[,1]
  hgrid.1 <- h[,2]
  neval <- ncol(distance_mat)
  d <- (distance_mat[,1] >= 0)
  N <- length(Y)
  N.1 <- sum(d)
  N.0 <- N - N.1
  Estimate=matrix(NA,neval,10)
  Estimate = as.data.frame(Estimate)
  colnames(Estimate)=c("b1", "b2", "h0", "h1", "N0","N1", "mu0","mu1","se0","se1")

  # to store matrices

  Indicators.0 <- list()
  Indicators.1 <- list()
  inv.designs.0 <- list()
  inv.designs.1 <- list()
  sig.halfs.0 <- list()
  sig.halfs.1 <- list()
  resd.0 <- list()
  resd.1 <- list()
  joint <- fitmethod == "joint"
  covariate.adjusted <- !is.null(covs.eff)
  covs.eff <- if (covariate.adjusted) as.matrix(covs.eff) else NULL
  n.cov <- if (covariate.adjusted) ncol(covs.eff) else 0
  clustered <- !is.null(cluster)
  clusters <- unique(cluster)

  # Check for mass points

  M.vec <- rep(N, neval)
  M.0.vec <- rep(N.0, neval)
  M.1.vec <- rep(N.1, neval)
  is_mass_point <- 0
  if (!is.null(masspoint.counts)) {
    M.vec <- masspoint.counts$M.vec
    M.0.vec <- masspoint.counts$M.0.vec
    M.1.vec <- masspoint.counts$M.1.vec
  } else if (masspoints == "check" | masspoints == "adjust"){
    for (j in 1:ncol(distance_mat)){
      distance <- distance_mat[,j]
      unique.counts <- rd2d_distance_unique_counts(distance)
      M <- unique.counts["M"]
      M.0 <- unique.counts["M.0"]
      M.1 <- unique.counts["M.1"]
      mass <- 1 - M / N
      M.vec[j] <- M
      M.0.vec[j] <- M.0
      M.1.vec[j] <- M.1
      if (mass >= 0.2){is_mass_point <- 1}
    }
    if (is_mass_point > 0){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & (masspoints == "check" | masspoints == "adjust")) bwcheck <- 50 + p + 1
    }
  }

  for (i in 1:neval) {

    b1 <- eval[i,1]
    b2 <- eval[i,2]
    distance <- distance_mat[,i]
    n.cov.i <- n.cov
    d <- (distance >= 0)
    distance <- abs(distance)
    side.0 <- !d
    side.1 <- d
    distance.0 <- distance[side.0]
    distance.1 <- distance[side.1]
    y.0 <- Y[side.0]
    y.1 <- Y[side.1]
    covs.0 <- if (covariate.adjusted) covs.eff[side.0, , drop = FALSE] else NULL
    covs.1 <- if (covariate.adjusted) covs.eff[side.1, , drop = FALSE] else NULL
    cluster.0 <- cluster[side.0]
    cluster.1 <- cluster[side.1]

    h.0 <- hgrid[i]
    h.1 <- hgrid.1[i]

    M.0 <- M.0.vec[i]
    M.1 <- M.1.vec[i]

    # bandwidth restrictions
    bw.min.0 <- NA
    bw.min.1 <- NA
    bw.max.0 <- NA
    bw.max.1 <- NA

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      bounds.0 <- rd2d_distance_bwcheck_limits(distance.0, bwcheck)
      bounds.1 <- rd2d_distance_bwcheck_limits(distance.1, bwcheck)
      bw.min.0   <- bounds.0["min"]
      bw.min.1   <- bounds.1["min"]
      bw.max.0   <- bounds.0["max"]
      bw.max.1   <- bounds.1["max"]
      h.0     <- max(h.0, bw.min.0)
      h.1     <- max(h.1, bw.min.1)
      h.0     <- min(h.0, bw.max.0)
      h.1     <- min(h.1, bw.max.1)
    }

    if (covariate.adjusted) {
      Y.gamma <- c(y.0, y.1)
      Z.gamma <- rbind(covs.0, covs.1)
      gamma <- rd2d_distance_covariate_gamma(
        Y.gamma, Z.gamma, c(-distance.0, distance.1), h.0, h.1, p, kernel
      )
      n.cov.i <- rd2d_covs_gamma_rank(gamma, colnames(covs.eff))
      Y.adj <- as.numeric(rd2d_distance_apply_covariate_adjustment(
        Y.gamma, Z.gamma, gamma
      ))
      n0 <- length(y.0)
      y.0 <- Y.adj[seq_len(n0)]
      y.1 <- Y.adj[n0 + seq_len(length(y.1))]
    }

    w.0   <- kernel_weight(distance.0/h.0, kernel = kernel)/c(h.0^2)
    w.1   <- kernel_weight(distance.1/h.1, kernel = kernel)/c(h.1^2)
    ind.0 <- (w.0 > 0)
    ind.1 <- (w.1 > 0)

    eN.0  <- sum(ind.0); eN.1 <- sum(ind.1)
    eY.0  <- y.0[ind.0]; eY.1  <- y.1[ind.1]
    u.0  <- distance.0[ind.0]; u.1  <- distance.1[ind.1]
    eweight.0 <- w.0[ind.0]; eweight.1 <- w.1[ind.1]
    eC.0 <- cluster.0[ind.0]; eC.1 <- cluster.1[ind.1]

    R.0 <- outer(u.0 / h.0, 0:p, `^`)
    R.1 <- outer(u.1 / h.1, 0:p, `^`)

    sqrtw_R.0 <- sqrt(eweight.0) * as.matrix(R.0)
    sqrtw_R.1 <- sqrt(eweight.1) * as.matrix(R.1)
    inv.gamma0 <- qrXXinv(sqrtw_R.0)
    inv.gamma1 <- qrXXinv(sqrtw_R.1)

    sigma0.half <- eweight.0 * as.matrix(R.0); sigma1.half <- eweight.1 * as.matrix(R.1)

    XtWY.0 <- t(matrix(eY.0 * unname(eweight.0), nrow = 1) %*% as.matrix(R.0)); XtWY.1 <-t(matrix(eY.1 * unname(eweight.1), nrow = 1) %*% as.matrix(R.1))
    hbeta.0 <- inv.gamma0 %*% XtWY.0; hbeta.1 <- inv.gamma1 %*% XtWY.1
    mu0 <- hbeta.0[1];  mu1 <- hbeta.1[1]
    res0 <- (eY.0 - as.matrix(R.0) %*% hbeta.0)[,1]
    res1 <- (eY.1 - as.matrix(R.1) %*% hbeta.1)[,1]

    # standard deviation
    if (vce=="hc0") {
      w.vce.0 = 1
      w.vce.1 = 1
    } else if (vce=="hc1") {
      if (clustered && (joint || covariate.adjusted)) {
        w.vce.0 = 1
        w.vce.1 = 1
      } else if (joint) {
        w.vce.0 = sqrt((eN.0 + eN.1) / (eN.0 + eN.1 - 2 * (p + 1) - n.cov.i))
        w.vce.1 = w.vce.0
      } else {
        w.vce.0 = sqrt(eN.0/(eN.0 - p - 1 - n.cov.i))
        w.vce.1 = sqrt(eN.1/(eN.1 - p - 1 - n.cov.i))
      }
    } else if (vce=="hc2") {
      hii.0 <- rowSums((sqrtw_R.0 %*% inv.gamma0) * sqrtw_R.0)
      w.vce.0 = sqrt(1/(1-hii.0))
      hii.1 <- rowSums((sqrtw_R.1 %*% inv.gamma1) * sqrtw_R.1)
      w.vce.1 = sqrt(1/(1-hii.1))
    } else if (vce == "hc3"){
      hii.0 <- rowSums((sqrtw_R.0 %*% inv.gamma0) * sqrtw_R.0)
      w.vce.0 = 1/(1-hii.0)
      hii.1 <- rowSums((sqrtw_R.1 %*% inv.gamma1) * sqrtw_R.1)
      w.vce.1 = 1/(1-hii.1)
    }

    res0 <- res0 * w.vce.0
    res1 <- res1 * w.vce.1

    ############################ start: get_cov_half ###########################
    if (is.null(cluster)){
      cov.half.const.0 <- sweep(sigma0.half, 1, res0, "*") %*% inv.gamma0
      cov.half.const.1 <- sweep(sigma1.half, 1, res1, "*") %*% inv.gamma1
    } else {
      k <- p + 1
      clusters <- unique(cluster)
      g <- length(clusters)
      n.0 <- length(eC.0); n.1 <- length(eC.1)
      g.0 <- length(unique(eC.0)); g.1 <- length(unique(eC.1))
      cluster.df <- !(joint && clustered) && !(covariate.adjusted && clustered)
      w.w.0 <- if (cluster.df) ((n.0 - 1) / (n.0 - k - n.cov.i)) * (g.0 / (g.0 - 1)) else 1
      w.w.1 <- if (cluster.df) ((n.1 - 1) / (n.1 - k - n.cov.i)) * (g.1 / (g.1 - 1)) else 1
      cov.half.const.0 <- matrix(0, nrow = g, ncol = k)
      cov.half.const.1 <- matrix(0, nrow = g, ncol = k)

      for (l in 1:g) {
        ind.vce.0 <- as.logical(eC.0 == clusters[l])
        ind.vce.1 <- as.logical(eC.1 == clusters[l])
        w_R_i.0 <- sigma0.half[ind.vce.0,,drop=FALSE]
        w_R_i.1 <- sigma1.half[ind.vce.1,,drop=FALSE]
        resd_i.0 <- matrix(res0[ind.vce.0], ncol = 1)
        resd_i.1 <- matrix(res1[ind.vce.1], ncol = 1)
        w_R_resd_i.0 <- t(crossprod(w_R_i.0, resd_i.0)) * sqrt(w.w.0)
        w_R_resd_i.1 <- t(crossprod(w_R_i.1, resd_i.1)) * sqrt(w.w.1)
        cov.half.const.0[l,] <- as.vector(w_R_resd_i.0)
        cov.half.const.1[l,] <- as.vector(w_R_resd_i.1)
      }
      cov.half.const.0 <- cov.half.const.0 %*% inv.gamma0
      cov.half.const.1 <- cov.half.const.1 %*% inv.gamma1
      if (joint || covariate.adjusted) {
        if (joint) {
          scale <- rd2d_joint_scale(
            vce, eN.0 + eN.1, 2 * k + n.cov.i, c(eC.0, eC.1), clustered
          )
          cov.half.const.0 <- cov.half.const.0 * scale
          cov.half.const.1 <- cov.half.const.1 * scale
        } else {
          scale.0 <- rd2d_joint_scale(vce, eN.0, k + n.cov.i, eC.0, clustered)
          scale.1 <- rd2d_joint_scale(vce, eN.1, k + n.cov.i, eC.1, clustered)
          cov.half.const.0 <- cov.half.const.0 * scale.0
          cov.half.const.1 <- cov.half.const.1 * scale.1
        }
      }
    }
    ############################ end: get_cov_half #############################

    cluster.df <- !(joint && clustered) && !(covariate.adjusted && clustered)
    sigma.0 <- rd2d_distance_vce(
      sigma0.half, res0, eC.0, h.0, p + 1 + n.cov.i, cluster.df, clusters
    )
    sigma.1 <- rd2d_distance_vce(
      sigma1.half, res1, eC.1, h.1, p + 1 + n.cov.i, cluster.df, clusters
    )

    cov.const.0 <- t(inv.gamma0) %*% sigma.0 %*% inv.gamma0
    cov.const.1 <- t(inv.gamma1) %*% sigma.1 %*% inv.gamma1
    if ((joint || covariate.adjusted) && clustered) {
      if (joint) {
        scale <- rd2d_joint_scale(
          vce, eN.0 + eN.1, 2 * (p + 1) + n.cov.i,
          c(eC.0, eC.1), clustered
        )
        cov.const.0 <- cov.const.0 * scale^2
        cov.const.1 <- cov.const.1 * scale^2
      } else {
        scale.0 <- rd2d_joint_scale(vce, eN.0, p + 1 + n.cov.i, eC.0, clustered)
        scale.1 <- rd2d_joint_scale(vce, eN.1, p + 1 + n.cov.i, eC.1, clustered)
        cov.const.0 <- cov.const.0 * scale.0^2
        cov.const.1 <- cov.const.1 * scale.1^2
      }
    }

    vec <- rep(0,p+1)
    vec[1] <- 1
    se0 <- matrix(vec, nrow = 1) %*% cov.const.0 %*% matrix(vec, ncol = 1) / (h.0^2)
    se0 <- sqrt(se0[1,1])
    se1 <- matrix(vec, nrow = 1) %*% cov.const.1 %*% matrix(vec, ncol = 1) / (h.1^2)
    se1 <- sqrt(se1[1,1])

    if (cbands){
      full.ind.0 <- rep(FALSE, N)
      full.ind.1 <- rep(FALSE, N)
      full.ind.0[which(side.0)] <- ind.0
      full.ind.1[which(side.1)] <- ind.1
      Indicators.0[[i]] <- full.ind.0
      Indicators.1[[i]] <- full.ind.1
      inv.designs.0[[i]] <- inv.gamma0
      inv.designs.1[[i]] <- inv.gamma1
      sig.halfs.0[[i]] <- cov.half.const.0
      sig.halfs.1[[i]] <- cov.half.const.1
      resd.0[[i]] <- res0
      resd.1[[i]] <- res1
    }

    Estimate[i,] <- c(b1, b2, h.0, h.1, eN.0, eN.1, mu0, mu1, se0, se1)
  }

  return(list("Estimate" = Estimate, "Indicators.0" = Indicators.0, "Indicators.1" = Indicators.1,
              "inv.designs.0" = inv.designs.0, "inv.designs.1" = inv.designs.1, "sig.halfs.0" = sig.halfs.0,
              "sig.halfs.1" = sig.halfs.1, "resd.0" = resd.0, "resd.1" = resd.1,
              "M.vec" = M.vec, "M.0.vec" = M.0.vec, "M.1.vec" = M.1.vec,
              "joint" = joint, "clustered" = clustered))
}

rd2d_distance_project_cov_halves <- function(halves, inds, hgrid, p,
                                             clustered = FALSE) {
  neval <- length(halves)
  if (neval == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  first <- which(vapply(halves, function(x) !is.null(x) && length(x) > 0, logical(1)))[1]
  if (is.na(first)) return(matrix(0, nrow = 0, ncol = neval))

  n.rows <- if (clustered) nrow(halves[[first]]) else length(inds[[first]])
  projected <- matrix(0, nrow = n.rows, ncol = neval)
  vec <- rep(0, p + 1)
  vec[1] <- 1

  for (i in seq_len(neval)) {
    half <- halves[[i]]
    if (is.null(half) || length(half) == 0) next

    h <- as.numeric(hgrid[i])
    score <- half %*% get_invH_dist(h, p) %*% matrix(vec, ncol = 1)
    score <- as.numeric(score)

    if (clustered) {
      projected[, i] <- score
    } else {
      projected[inds[[i]], i] <- score
    }
  }

  projected
}

rd2d_distance_project_fit_sides <- function(fit, p, clustered = FALSE) {
  list(
    `0` = rd2d_distance_project_cov_halves(
      fit$sig.halfs.0, fit$Indicators.0, fit$Estimate$h0, p, clustered
    ),
    `1` = rd2d_distance_project_cov_halves(
      fit$sig.halfs.1, fit$Indicators.1, fit$Estimate$h1, p, clustered
    )
  )
}

rd2d_distance_cov_from_projects <- function(project.a, project.b,
                                            side = c("both", "0", "1"),
                                            symmetric = FALSE,
                                            joint = FALSE,
                                            clustered = FALSE) {
  side <- match.arg(side)

  if (side == "0") {
    cov.out <- crossprod(project.a$`0`, project.b$`0`)
  } else if (side == "1") {
    cov.out <- crossprod(project.a$`1`, project.b$`1`)
  } else if (clustered && joint) {
    cov.out <- crossprod(
      project.a$`1` - project.a$`0`,
      project.b$`1` - project.b$`0`
    )
  } else {
    cov.out <- crossprod(project.a$`0`, project.b$`0`) +
      crossprod(project.a$`1`, project.b$`1`)
  }

  if (symmetric) cov.out <- (cov.out + t(cov.out)) / 2
  cov.out
}

rd2d_distance_cov_from_fits <- function(fit.a, fit.b, p, clustered = FALSE,
                                        side = c("both", "0", "1")) {
  side <- match.arg(side)
  same.fit <- identical(fit.a, fit.b)
  project.a <- rd2d_distance_project_fit_sides(fit.a, p, clustered)
  project.b <- if (same.fit) project.a else {
    rd2d_distance_project_fit_sides(fit.b, p, clustered)
  }
  cov.out <- rd2d_distance_cov_from_projects(
    project.a, project.b, side = side, symmetric = same.fit,
    joint = isTRUE(fit.a$joint), clustered = clustered
  )
  cov.out
}

rd2d_distance_fuzzy_cov_tables <- function(fit.itt, fit.fs, p, tau.itt,
                                           tau.fs, outputs, clustered = FALSE,
                                           joint = FALSE,
                                           denom.tol = sqrt(.Machine$double.eps)) {
  outputs <- unique(outputs)
  neval <- nrow(fit.itt$Estimate)
  covs <- lapply(outputs, function(x) matrix(NA_real_, nrow = neval, ncol = neval))
  names(covs) <- outputs

  valid <- is.finite(tau.itt) & is.finite(tau.fs) & abs(tau.fs) > denom.tol

  need.main <- "main" %in% outputs
  need.itt <- any(c("main", "itt") %in% outputs)
  need.fs <- any(c("main", "fs") %in% outputs)
  need.itt.side <- any(c("itt.0", "itt.1") %in% outputs)
  need.fs.side <- any(c("fs.0", "fs.1") %in% outputs)

  if (need.itt || need.itt.side) {
    project.itt <- rd2d_distance_project_fit_sides(fit.itt, p, clustered)
  }
  if (need.fs || need.fs.side) {
    project.fs <- rd2d_distance_project_fit_sides(fit.fs, p, clustered)
  }

  if (need.itt) {
    cov.itt <- rd2d_distance_cov_from_projects(
      project.itt, project.itt, side = "both", symmetric = TRUE,
      joint = joint, clustered = clustered
    )
  }
  if (need.fs) {
    cov.fs <- rd2d_distance_cov_from_projects(
      project.fs, project.fs, side = "both", symmetric = TRUE,
      joint = joint, clustered = clustered
    )
  }

  if ("itt" %in% outputs) covs$itt <- cov.itt
  if ("fs" %in% outputs) covs$fs <- cov.fs

  if (need.main) {
    cov.itt.fs <- rd2d_distance_cov_from_projects(
      project.itt, project.fs, side = "both", joint = joint,
      clustered = clustered
    )
    cov.fs.itt <- t(cov.itt.fs)

    grad.itt <- 1 / tau.fs
    grad.fs <- -tau.itt / tau.fs^2
    cov.main <- outer(grad.itt, grad.itt) * cov.itt +
      outer(grad.itt, grad.fs) * cov.itt.fs +
      outer(grad.fs, grad.itt) * cov.fs.itt +
      outer(grad.fs, grad.fs) * cov.fs
    cov.main[!outer(valid, valid)] <- NA_real_
    covs$main <- (cov.main + t(cov.main)) / 2
  }

  if ("itt.0" %in% outputs) {
    covs$itt.0 <- rd2d_distance_cov_from_projects(
      project.itt, project.itt, side = "0", symmetric = TRUE
    )
  }
  if ("itt.1" %in% outputs) {
    covs$itt.1 <- rd2d_distance_cov_from_projects(
      project.itt, project.itt, side = "1", symmetric = TRUE
    )
  }
  if ("fs.0" %in% outputs) {
    covs$fs.0 <- rd2d_distance_cov_from_projects(
      project.fs, project.fs, side = "0", symmetric = TRUE
    )
  }
  if ("fs.1" %in% outputs) {
    covs$fs.1 <- rd2d_distance_cov_from_projects(
      project.fs, project.fs, side = "1", symmetric = TRUE
    )
  }

  covs
}

########################### Standardization ####################################

# Standardize input data and point of evaluations.

lpgeo_bwselect_std <- function(x, eval){

  scale.1 <- sd(x$x.1)
  scale.2 <- sd(x$x.2)

  x.std <- x
  eval.std <- eval
  x.std$x.1 <- x.std$x.1 / scale.1
  x.std$x.2 <- x.std$x.2 / scale.2
  eval.std$x.1 <- eval.std$x.1 / scale.1
  eval.std$x.2 <- eval.std$x.2 / scale.2

  return(list(x.std = x.std, eval.std = eval.std, scale.1 = scale.1, scale.2 = scale.2))
}

###################################### Misc ####################################

pow <- function(vec,p){
  # if input is a vecotr
  x <- vec[1]
  y <- vec[2]
  # vec: vector
  # p: power to raise
  # d: dimension of vec
  res <- c(1)
  if (p>= 1){
    for (i in 1:p){
      for (j in 0:i){
        res <- c(res, x^(i-j) * y^(j))
      }
    }
  }
  return(res)
}

# pow.d
pow.d <- function(x,p){
  # if input is a numeric
  res <- rep(NA,p+1)
  res[1] <- 1
  if (p >= 1){
    for (i in 1:p){
      res[i+1] <- x^i
    }
  }
  return(res)
}

get_basis_v1 <- function(u,p) {
  u <- as.matrix(u)
  if (dim(u)[2] == 1){
    pow.temp <- function(v){ return(pow.d(v,p))}
  } else {
    pow.temp <- function(v){ return(pow(v,p)) }
  }
  tmp <- apply(u,c(1),pow.temp,simplify = FALSE)
  res <- data.frame(do.call(rbind,as.matrix(tmp)))
  return(res)
}

rd2d_distance_unique <- function(distance){

  ord <- order(distance)

  distance <- distance[ord]

  N <- length(distance)

  # if x has one or no element
  if (N == 0) return(list(unique = NULL, freq = c(), index = c()))
  if (N == 1) return(list(unique = distance, freq = 1, index = 1))

  # else
  uniqueIndex <- c(c(distance[2:N] != distance[1:(N-1)]), TRUE)
  unique <- distance[uniqueIndex]
  nUnique <- length(unique)

  # all are distinct
  if (nUnique == N) return(list(unique=unique, freq=rep(1,N), index=1:N))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=N, index=N))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, index=(1:N)[uniqueIndex]))
}

get_invH_dist <- function(h,p){
  result <- rep(NA, p+1)
  result[1] <- 1

  if (p >= 1){
    for (j in 1:p){
      result[j+1] <- 1/h^j
    }
  }
  result <- diag(result)
  return(result)
}
