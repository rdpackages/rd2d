################################################################################
#' @title Bandwidth Selection for Distance-Based Methods for Boundary Discontinuity Design
#' @description
#' \code{rdbw2d.distance} implements bandwidth selectors for distance-based local
#' polynomial boundary discontinuity (BD) estimation and inference.
#' The command targets level effects at two-dimensional boundary points.
#' See \href{https://arxiv.org/abs/2510.26051}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026)}
#' for methodological background.
#'
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param distance Signed distance scores; a numeric matrix or data frame of
#' dimension \eqn{N \times J}, where \eqn{N} is the sample size and \eqn{J}
#' is the number of evaluation points. Non-negative values identify
#' observations on the treated side and negative values identify observations
#' on the control side.
#' @param b Optional evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, dimension \eqn{J \times 2}.
#' @param p Polynomial order for point estimation. Default is \code{p = 1}.
#' @param kink.unknown Logical value or vector of length 2 controlling
#' unknown-kink bandwidth adjustments. A scalar \code{TRUE} is expanded to
#' \code{c(TRUE, TRUE)}, and a scalar \code{FALSE} is expanded to
#' \code{c(FALSE, FALSE)}. With the vector form, the first element controls
#' whether the point-estimation bandwidth is shrunk to the unknown-kink rate;
#' the second controls whether the inference bandwidth is further shrunk in
#' \code{\link{rd2d.distance}}. For bandwidth selection, only the first element
#' is used. Default is \code{c(FALSE, FALSE)}.
#' @param kink.position Optional boundary positions of known kink points.
#' Either a logical vector with one entry per boundary point, where \code{TRUE}
#' identifies a kink point, or an integer vector with indices between 1 and
#' the number of boundary points. Requires \code{b} so distances between
#' boundary points can be computed.
#' @param kernel Kernel function to use. Options are \code{"tri"} or
#' \code{"triangular"} (triangular, default), \code{"epa"} or
#' \code{"epanechnikov"} (Epanechnikov), \code{"uni"} or \code{"uniform"}
#' (uniform), and \code{"gau"} or \code{"gaussian"} (Gaussian).
#' @param bwselect Bandwidth selection strategy. Options:
#' \itemize{
#' \item \code{"mserd"}. One common MSE-optimal bandwidth selector for the boundary RD treatment effect estimator for each evaluation point (default).
#' \item \code{"cerrd"}. CER-optimal counterpart of \code{"mserd"}.
#' \item \code{"imserd"}. IMSE-optimal bandwidth selector for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"icerrd"}. CER-optimal counterpart of \code{"imserd"}.
#' \item \code{"msetwo"}. Two different MSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator for each evaluation point.
#' \item \code{"certwo"}. CER-optimal counterpart of \code{"msetwo"}.
#' \item \code{"imsetwo"}. Two IMSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"icertwo"}. CER-optimal counterpart of \code{"imsetwo"}.
#' }
#' @param bwparam Target parameter used for fuzzy automatic bandwidth
#' selection. Options are \code{"main"}, which selects bandwidths for the
#' linearized fuzzy Wald ratio, and \code{"itt"}, which selects bandwidths for
#' the reduced-form outcome. Ignored in sharp designs.
#' @param vce Variance-covariance estimator for standard errors.
#' Options:
#' \describe{
#'   \item{\code{"hc0"}}{Heteroskedasticity-robust variance estimator without small sample adjustment (White robust).}
#'   \item{\code{"hc1"}}{Heteroskedasticity-robust variance estimator with degrees-of-freedom correction. (default)}
#'   \item{\code{"hc2"}}{Heteroskedasticity-robust variance estimator using leverage adjustments.}
#'   \item{\code{"hc3"}}{More conservative heteroskedasticity-robust variance estimator (similar to jackknife correction).}
#' }
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used. Default is \code{20 + p + 1}.
#' @param masspoints Strategy for handling mass points in the running variable.
#' Options:
#' \describe{
#'   \item{\code{"check"}}{(default) Check for repeated values and adjust inference if needed.}
#'   \item{\code{"adjust"}}{Adjust bandwidths to guarantee a sufficient number of unique support points.}
#'   \item{\code{"off"}}{Ignore mass points completely.}
#' }
#' @param cluster Cluster ID variable used for cluster-robust variance estimation with degrees-of-freedom weights. Default is \code{cluster = NULL}.
#' @param scaleregul Scaling factor for the regularization term in bandwidth selection. Default is \code{1}.
#' @param cqt Constant controlling subsample fraction for initial bias estimation. Default is \code{0.5}.
#' @param fuzzy Optional treatment receipt/status variable used for fuzzy RD
#' bandwidth selection. If supplied, \code{bwparam} controls whether the
#' selector targets the fuzzy Wald ratio or the reduced-form outcome.
#'
#' @return An object of class \code{"rdbw2d.distance"}, containing:
#' \describe{
#'   \item{\code{bws}}{Data frame of optimal bandwidths for each evaluation point:
#'     \describe{
#'       \item{\code{b1}}{First coordinate of the evaluation point \eqn{b1}.}
#'       \item{\code{b2}}{Second coordinate of the evaluation point \eqn{b2}.}
#'       \item{\code{h0}}{Bandwidth for observations with negative signed distance.}
#'       \item{\code{h1}}{Bandwidth for observations with non-negative signed distance.}
#'     }
#'   }
#'   \item{\code{mseconsts}}{Data frame of intermediate bias and variance constants used for MSE/IMSE calculations, including \code{N.Co} and \code{N.Tr}, the effective sample sizes for observations with negative and non-negative signed distances, respectively.}
#'   \item{\code{opt}}{A list of options and settings used in estimation, including \code{p}, \code{kernel}, sample size \eqn{N}, and user-specified choices.}
#'   \item{\code{call}}{Matched function call.}
#' }
#'
#' @seealso \code{\link{rd2d.distance}}, \code{\link{rd2d}}, \code{\link{summary.rdbw2d.distance}}, \code{\link{print.rdbw2d.distance}}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @references
#' \itemize{
#' \item{Cattaneo, M. D., and Titiunik, R. (2022).
#' Regression Discontinuity Designs. \doi{10.1146/annurev-economics-051520-021409}.}
#' \item{\href{https://arxiv.org/abs/2510.26051}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods.}
#' \item{\href{https://arxiv.org/abs/2511.06474}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Boundary Discontinuity Designs: Theory and Practice.}
#' \item{\href{https://arxiv.org/abs/2505.07989}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2025).}
#' rd2d: Causal Inference in Boundary Discontinuity Designs.}
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 800
#'
#' # Generate running variables x1 and x2
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#'
#' # Define treatment assignment: treated if x1 >= 0
#' assignment <- as.numeric(x1 >= 0)
#'
#' # Generate outcome variable Y with some treatment effect
#' Y <- 3 + 2 * x1 + 1.5 * x2 + 1.5 * assignment + rnorm(n, sd = 0.5)
#'
#' # Define evaluation points (e.g., at the origin and another point)
#' eval <- data.frame(x.1 = c(0, 0), x.2 = c(0, 1))
#'
#' # Compute Euclidean distances to evaluation points
#' distance.a <- sqrt((x1 - eval$x.1[1])^2 + (x2 - eval$x.2[1])^2)
#' distance.b <- sqrt((x1 - eval$x.1[2])^2 + (x2 - eval$x.2[2])^2)
#'
#' # Combine distances into a matrix
#' distance <- as.data.frame(cbind(distance.a, distance.b))
#'
#' # Assign positive distances for treatment group, negative for control
#' assignment_expanded <- matrix(rep(2 * assignment - 1, times = ncol(distance)),
#'                               nrow = nrow(distance), ncol = ncol(distance))
#' distance <- distance * assignment_expanded
#'
#' # Run the rdbw2d.distance function
#' bws <- rdbw2d.distance(Y, distance = distance, b = eval, masspoints = "off", bwcheck = 10)
#'
#' # View the estimation results
#' print(bws)
#' summary(bws)
#'
#' # Fuzzy distance-based bandwidths
#' fuzzy <- as.numeric(runif(n) < ifelse(assignment == 1, 0.8, 0.2))
#' Y.fuzzy <- 3 + 2 * x1 + 1.5 * x2 + 1.5 * fuzzy + rnorm(n, sd = 0.5)
#' bws.fuzzy <- rdbw2d.distance(Y.fuzzy, distance = distance, b = eval, fuzzy = fuzzy,
#'                          bwparam = "main", masspoints = "off",
#'                          bwcheck = 10)
#' print(bws.fuzzy)
#' @export



rdbw2d.distance <- function(Y, distance, b = NULL, p = 1,
                   kink.unknown = c(FALSE, FALSE), kink.position = NULL,
                   kernel = c("tri","triangular", "epa","epanechnikov","uni","uniform","gau","gaussian"),
                   bwselect = c("mserd", "cerrd", "imserd", "icerrd",
                                "msetwo", "certwo", "imsetwo", "icertwo"),
                   bwparam = c("main", "itt"),
                   vce = c("hc1","hc0","hc2","hc3"),
                   bwcheck = 20 + p + 1, masspoints = c("check","adjust","off"),
                   cluster = NULL, scaleregul = 1, cqt = 0.5, fuzzy = NULL){

  # Input error handling

  bwselect <- match.arg(bwselect)
  bwselect.base <- rd2d_bwselect_base(bwselect)
  bwselect.cer <- rd2d_bwselect_is_cer(bwselect)
  bwselect.common <- rd2d_bwselect_is_common(bwselect)
  kernel <- match.arg(kernel)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  kink.unknown <- rd2d_validate_kink_unknown(kink.unknown)
  bwparam <- match.arg(bwparam)
  rot <- NULL
  is.fuzzy <- !is.null(fuzzy)
  if (!is.fuzzy) bwparam <- "main"
  distance_mat <- distance

  # Check Errors

  exit=0

  if (length(Y) != nrow(distance_mat)) {
    print("Y and rows of distance must have the same length")
    exit <- 1
  }

  if (is.fuzzy && length(fuzzy) != length(Y)) {
    print("fuzzy must have the same length as Y")
    exit <- 1
  }

  if (is.fuzzy && !(is.logical(fuzzy) || is.numeric(fuzzy))) {
    print("fuzzy must be a logical or numeric vector")
    exit <- 1
  }

  if (kernel!="gau" & kernel!="gaussian" & kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit <- 1
  }

  if (!is.null(b)){
    if (nrow(b) != ncol(distance_mat) || ncol(b) != 2){
      print("b must have 2 columns and the same number of rows as distance's number of columns")
      exit <- 1
    }
  }

  if (is.null(p)) p <- 1
  valid.p <- is.numeric(p) && length(p) == 1 &&
    is.finite(p) && p >= 0 &&
    abs(p - round(p)) < sqrt(.Machine$double.eps)
  if (!valid.p) {
    print("p must be a nonnegative integer")
    exit <- 1
  } else {
    p <- as.integer(round(p))
  }

  if (!is.null(cluster) && !(vce %in% c("hc0", "hc1"))) {
    warning("When cluster is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  # Data Cleaning

  neval <- ncol(distance_mat)
  kink.position <- rd2d_validate_kink_position(kink.position, neval, b)
  if (any(kink.position) && kink.unknown[1]) {
    stop("Use either kink.position or kink.unknown, not both.", call. = FALSE)
  }

  if (!is.null(b)){
    eval <- as.data.frame(b)
  } else {
    eval <- matrix(NA, nrow = neval, ncol = 2)
    eval <- as.data.frame(eval)
  }

  colnames(eval) <- c("x.1", "x.2")

  na.ok <- complete.cases(distance_mat) & complete.cases(Y)
  if (is.fuzzy) na.ok <- na.ok & complete.cases(fuzzy)
  distance_mat <- distance_mat[na.ok,,drop = FALSE]
  Y <- Y[na.ok]
  if (is.fuzzy) fuzzy <- as.numeric(fuzzy[na.ok])
  cluster <- cluster[na.ok]
  d <- (distance_mat[,1] >= 0)

  N <- length(Y)
  N.1 <- sum(d)
  N.0 <- N - N.1

  kernel   <- tolower(kernel)

  e_level <- matrix(0, nrow = neval, ncol = p+1)
  e_level[,1] <- 1

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Store variance and bias constants for IMSE

  bconst <- rep(NA, neval)
  vconst <- rep(NA, neval)

  # Check for mass points

  M.vec <- rep(N, neval)
  M.0.vec <- rep(N.0, neval)
  M.1.vec <- rep(N.1, neval)
  is_mass_point <- 0
  if (masspoints == "check" | masspoints == "adjust"){
    for (j in 1:ncol(distance_mat)){
      distance.score <- distance_mat[,j]
      unique.const <- rd2d_distance_unique(distance.score)
      unique <- unique.const$unique
      M.0 <- sum(unique < 0)
      M.1 <- sum(unique >= 0)
      M <- M.0 + M.1
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

  results <- rdbw2d_distance_bw(
    Y = Y, distance = distance_mat, p = p, kernel = kernel, target.vec = e_level, rot = rot,
    vce = vce, cluster = cluster, bwcheck = bwcheck, scaleregul = scaleregul,
    cqt = cqt, fuzzy = if (is.fuzzy) fuzzy else NULL, bwparam = bwparam
  )

  if (bwselect.base == "mserd"){
    hn.grid <- ( 2 * (results$v.0   + results$v.1) /
              ( (2 * p + 2) * ( (results$b.0 - results$b.1)^2 +
              scaleregul * results$r.0 + scaleregul * results$r.1) ) )^(1/(2 * p + 4))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid <- pmax(hn.grid, results$bw.min.0, results$bw.min.1)
      hn.grid <- pmin(hn.grid, results$bw.max.0, results$bw.max.1)
    }
    results$h.0 <- hn.grid
    results$h.1 <- hn.grid
  }

  if (bwselect.base == "msetwo"){
    hn.grid.0 <- ( 2 * results$v.0 /
                   ( (2 * p + 2) * ( results$b.0^2 + scaleregul * results$r.0) ) )^(1/(2 * p + 4))
    hn.grid.1 <- ( 2 * results$v.1 /
                   ( (2 * p + 2) * ( results$b.1^2 + scaleregul * results$r.1) ) )^(1/(2 * p + 4))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid.0 <- pmax(hn.grid.0, results$bw.min.0)
      hn.grid.0 <- pmin(hn.grid.0, results$bw.max.0)
      hn.grid.1 <- pmax(hn.grid.1, results$bw.min.1)
      hn.grid.1 <- pmin(hn.grid.1, results$bw.max.1)
    }
    results$h.0 <- hn.grid.0
    results$h.1 <- hn.grid.1
  }

  if (bwselect.base == "imserd"){
    V.V <- mean(results$v.0) + mean(results$v.1)
    B.B <- mean( (results$b.0 - results$b.1)^2 + scaleregul * results$r.0 + scaleregul * results$r.1)
    hIMSE <- (2 * V.V / ( (2 * p + 2)* B.B ) )^(1/(2 * p + 4))
    hn.grid <- rep(hIMSE, nrow(results))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid <- pmax(hn.grid, results$bw.min.0, results$bw.min.1)
      hn.grid <- pmin(hn.grid, results$bw.max.0, results$bw.max.1)
    }
    results$h.0 <- hn.grid
    results$h.1 <- hn.grid
  }

  if (bwselect.base == "imsetwo"){
    V.V.0 <- mean(results$v.0)
    V.V.1 <- mean(results$v.1)
    B.B.0 <- mean( results$b.0^2 + scaleregul * results$r.0)
    B.B.1 <- mean( results$b.1^2 + scaleregul * results$r.1)
    hIMSE.0 <- (2 * V.V.0 / ( (2 * p + 2)* B.B.0 ) )^(1/(2 * p + 4))
    hIMSE.1 <- (2 * V.V.1 / ( (2 * p + 2)* B.B.1 ) )^(1/(2 * p + 4))
    hn.grid.0 <- rep(hIMSE.0, nrow(results))
    hn.grid.1 <- rep(hIMSE.1, nrow(results))

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid.0 <- pmax(hn.grid.0, results$bw.min.0)
      hn.grid.0 <- pmin(hn.grid.0, results$bw.max.0)
      hn.grid.1 <- pmax(hn.grid.1, results$bw.min.1)
      hn.grid.1 <- pmin(hn.grid.1, results$bw.max.1)
    }
    results$h.0 <- hn.grid.0
    results$h.1 <- hn.grid.1
  }

  smooth.exp <- if (bwselect.cer) 1 / (p + 4) else 1 / (2 * p + 4)

  if (bwselect.cer) {
    if (bwselect.common) {
      results$h.0 <- results$h.0 * rd2d_cer_factor(M.vec, p)
      results$h.1 <- results$h.1 * rd2d_cer_factor(M.vec, p)
    } else {
      results$h.0 <- results$h.0 * rd2d_cer_factor(M.0.vec, p)
      results$h.1 <- results$h.1 * rd2d_cer_factor(M.1.vec, p)
    }
  }

  # Kink adjustments. Known-kink bandwidths use the adaptive selector based on
  # distance to the nearest kink point; unknown-kink bandwidths use the global
  # non-smooth rate.

  if (any(kink.position)) {
    distance.to.kink <- rd2d_distance_to_known_kink(eval, kink.position)
    if (bwselect.common) {
      h.rate <- results$h.0 * rd2d_bw_rate_factor(M.vec, smooth.exp, 1 / 4)
      h.bound <- pmax(h.rate, distance.to.kink)
      results$h.0 <- pmin(results$h.0, h.bound)
      results$h.1 <- pmin(results$h.1, h.bound)
    } else {
      h.rate.0 <- results$h.0 * rd2d_bw_rate_factor(M.0.vec, smooth.exp, 1 / 4)
      h.rate.1 <- results$h.1 * rd2d_bw_rate_factor(M.1.vec, smooth.exp, 1 / 4)
      results$h.0 <- pmin(results$h.0, pmax(h.rate.0, distance.to.kink))
      results$h.1 <- pmin(results$h.1, pmax(h.rate.1, distance.to.kink))
    }
  }

  if (kink.unknown[1]){
    if (bwselect.common){
      results$h.0 <- results$h.0 * rd2d_bw_rate_factor(M.vec, smooth.exp, 1 / 4)
      results$h.1 <- results$h.1 * rd2d_bw_rate_factor(M.vec, smooth.exp, 1 / 4)
    } else {
      results$h.0 <- results$h.0 * rd2d_bw_rate_factor(M.0.vec, smooth.exp, 1 / 4)
      results$h.1 <- results$h.1 * rd2d_bw_rate_factor(M.1.vec, smooth.exp, 1 / 4)
    }
  }

  # Outputs

  bws <- results[,c("h.0", "h.1")]
  bws <- cbind(eval, bws)
  colnames(bws) <- c("b1","b2","h0", "h1")
  rownames(bws) <- NULL
  rownames(results) <- NULL

  clustered <- !is.null(cluster)

  out        <- list(bws = bws, mseconsts = results,
                     opt = list(N=N, N.0 = N.0, N.1 = N.1, M = M.vec, M.0 = M.0.vec,
                                M.1 = M.1.vec, neval=neval, p=p, b = eval,
                                kernel=kernel.type,
                                kink.unknown = kink.unknown,
                                kink.position = kink.position,
                                bwselect=bwselect, bwparam = bwparam,
                                bwcheck = bwcheck,
                                cluster= cluster, clustered = clustered,
                                fuzzy = is.fuzzy,
                                vce = vce, masspoints = masspoints,
                                scaleregul = scaleregul, cqt = cqt))
  out$call   <- match.call()
  class(out) <- "rdbw2d.distance"
  return(out)
}

# Print method

################################################################################
#' Print Method for Bandwidth Selection (Distance-Based) in 2D Local Polynomial RD Design
#'
#' @description
#' Print method for displaying summary information from distance-based bandwidth selection in 2D local polynomial regression discontinuity (RD) designs, as produced by \code{\link{rdbw2d.distance}}.
#'
#' @param x An object of class \code{rdbw2d.distance}, returned by \code{\link{rdbw2d.distance}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return No return value. This function is called for its side effects: it prints summary information of \code{\link{rdbw2d.distance}}.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso \code{\link{rdbw2d.distance}} for distance-based bandwidth selection in 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rdbw2d.distance}}, \code{\link{summary.rdbw2d.distance}}.
#' @export

print.rdbw2d.distance <- function(x,...){
  cat("Call: rdbw2d.distance\n\n")

  cat(sprintf("Number of Obs.         %-10s\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("Unknown kink           %s\n", paste(tolower(x$opt$kink.unknown), collapse = ", ")))
  if (any(x$opt$kink.position)) {
    cat(sprintf("Known kink positions   %s\n", paste(which(x$opt$kink.position), collapse = ", ")))
  }
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  if (isTRUE(x$opt$fuzzy)) {
    cat(sprintf("Fuzzy                  on\n"))
    cat(sprintf("BW target              %s\n", x$opt$bwparam))
  }
  cat("\n")

  cat(sprintf("Number of Obs.         %-10s%-10s\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Order est. (p)         %-10s%-10s\n", x$opt$p, x$opt$p))
  cat("\n")
  #cat("Use summary(...) to show bandwidths.\n")
}

################################################################################
#' Summary Method for Bandwidth Selection in 2D Local Polynomial RD Design (Distance-Based)
#'
#' @description
#' Summarizes bandwidth selection results from a 2D local polynomial regression discontinuity (RD) design using distance-based methods, as returned by \code{\link{rdbw2d.distance}}.
#'
#' @param object An object of class \code{rdbw2d.distance}, returned by \code{\link{rdbw2d.distance}}.
#' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{subset}: Integer vector of indices of evaluation points to display. Defaults to all evaluation points.
#'     \item \code{sep}: Integer vector of length two. Controls spacing in the output.
#'       \code{sep[1]} controls spacing for the columns of evaluation points in the table.
#'       \code{sep[2]} controls spacing for the columns of bandwidths in the table.
#'       Default is \code{c(8, 14)}.
#'   }
#'
#' @return No return value. This function is called for its side effects: it prints a formatted summary of \code{\link{rdbw2d.distance}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso \code{\link{rdbw2d.distance}} for bandwidth selection using 2D local polynomial RD design with distance-based methods.
#'
#' Supported methods: \code{\link{print.rdbw2d.distance}}, \code{\link{summary.rdbw2d.distance}}.
#'
#' @export

summary.rdbw2d.distance <- function(object,...) {

  x <- object

  args <- list(...)

  if (is.null(args[['subset']])) {
    subset <- NULL
  } else {
    subset <- args[['subset']]
  }

  if (is.null(args[['sep']])) {
    sep <- c(8, 14)
  } else {
    sep <- args[['sep']]
  }

  cat("Call: rdbw2d.distance\n\n")

  cat(sprintf("Number of Obs.         %-10s\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("Unknown kink           %s\n", paste(tolower(x$opt$kink.unknown), collapse = ", ")))
  if (any(x$opt$kink.position)) {
    cat(sprintf("Known kink positions   %s\n", paste(which(x$opt$kink.position), collapse = ", ")))
  }
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  if (isTRUE(x$opt$fuzzy)) {
    cat(sprintf("Fuzzy                  on\n"))
    cat(sprintf("BW target              %s\n", x$opt$bwparam))
  }
  cat("\n")

  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat("\n")

  cat("Bandwidth Selection","\n")

  # Define column headers and their widths

  eval.specified <- !all(is.na(x$opt$b)) # TRUE is any entry is specified

  if (eval.specified){
    headers <- c("ID", "b1", "b2", "h0", "h1")
    col_widths <- c(4, sep[1], sep[1], sep[2], sep[2])
  } else {
    headers <- c("ID", "h0", "h1")
    col_widths <- c(4, sep[2], sep[2])
  }

  col_widths_default <- c(4, sep[1], sep[1], sep[2], sep[2])

  cat(strrep("=", sum(col_widths)), "\n")

  # Format headers using formatC for right alignment
  if (eval.specified){
    group_headers <- c(
      formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
      formatC("BW Control", width = col_widths_default[4], format = "s"),
      formatC("BW Treatment", width = col_widths_default[5], format = "s")
    )
  } else{
    group_headers <- c(
      formatC("", width = col_widths[1], format = "s"),
      formatC("BW Control", width = col_widths_default[4], format = "s"),
      formatC("BW Treatment", width = col_widths_default[5], format = "s")
    )
  }

  cat(paste(group_headers, collapse = ""), "\n")

  formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
  cat(paste(formatted_headers, collapse = ""), "\n")

  # Print separator line
  cat(strrep("=", sum(col_widths)), "\n")

  neval <- nrow(x$bws)
  if (is.null(subset)){
    subset <- seq_len(neval)
  } else{
    # input error handling
    if (!all(subset %in% seq_len(neval))) {
      warning("Invalid subset provided. Resetting to default: 1:neval")
      subset <- seq_len(neval)
    }
  }

  # Print each row of bandwidth estimates
  for (j in 1:nrow(x$bws)) {
    index <- formatC(j, width = col_widths_default[1], format = "d")
    bdy1 <- ifelse(is.na(x$bws[j, "b1"]),
                   formatC("NA", width = col_widths_default[2], format = "s"),
                   formatC(x$bws[j, "b1"], format = "f", digits = 3, width = col_widths_default[2]))
    bdy2 <- ifelse(is.na(x$bws[j, "b2"]),
                   formatC("NA", width = col_widths_default[3], format = "s"),
                   formatC(x$bws[j, "b2"], format = "f", digits = 3, width = col_widths_default[3]))
    control <- ifelse(is.na(x$bws[j, "h0"]),
                      formatC("NA", width = col_widths_default[4], format = "s"),
                      formatC(x$bws[j, "h0"], format = "f", digits = 3, width = col_widths_default[4]))
    treatment <- ifelse(is.na(x$bws[j, "h1"]),
                        formatC("NA", width = col_widths_default[5], format = "s"),
                        formatC(x$bws[j, "h1"], format = "f", digits = 3, width = col_widths_default[5]))

    # Combine formatted values and print the row
    if (eval.specified){
      row_vals <- c(index, bdy1, bdy2, control, treatment)
    } else {
      row_vals <- c(index, control, treatment)
    }

    if (j %in% subset) cat(paste(row_vals, collapse = ""), "\n")
  }

  # Print closing separator line
  cat(strrep("=", sum(col_widths)), "\n")
}
