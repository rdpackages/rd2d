################################################################################
#' @title Distance-Based Methods for Boundary Discontinuity Design
#'
#' @description
#' \code{rd2d.distance} implements distance-based local polynomial boundary
#' discontinuity (BD) point estimators with robust bias-corrected pointwise
#' confidence intervals and uniform confidence bands.
#' The command targets level effects at two-dimensional boundary points.
#' See \href{https://arxiv.org/abs/2510.26051}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026)}
#' for methodological background.
#'
#' Companion commands are: \code{rdbw2d.distance} for data-driven bandwidth selection.
#'
#' For other packages of RD designs, visit
#' <https://rdpackages.github.io/>
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param distance Signed distance scores; a numeric matrix or data frame of
#' dimension \eqn{N \times J}, where \eqn{N} is the sample size and \eqn{J}
#' is the number of evaluation points. Non-negative values identify
#' observations on the treated side and negative values identify observations
#' on the control side.
#' @param h Bandwidths. A positive scalar uses the same bandwidth for both
#' groups and all evaluation points. A matrix/data frame of size
#' \eqn{J \times 2} uses row-specific control and treated bandwidths. If not
#' specified, bandwidths are selected by \code{\link{rdbw2d.distance}}.
#' @param b Optional evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, dimension \eqn{J \times 2}.
#' @param p Polynomial order for point estimation. Default is \code{p = 1}.
#' @param q Polynomial order for bias-corrected estimation. Must satisfy
#' \eqn{q \geq p}. Default is \code{q = p + 1}, except when
#' \code{kink.unknown[1] = TRUE}, where the default is \code{q = p}.
#' @param kink.unknown Logical value or vector of length 2 controlling
#' unknown-kink bandwidth adjustments. A scalar \code{TRUE} is expanded to
#' \code{c(TRUE, TRUE)}, and a scalar \code{FALSE} is expanded to
#' \code{c(FALSE, FALSE)}. With the vector form, the first element controls
#' whether the point-estimation bandwidth is shrunk to the unknown-kink rate;
#' the second controls whether the inference bandwidth is further shrunk.
#' Default is \code{c(FALSE, FALSE)}. The second element can be \code{TRUE}
#' only when the first element is \code{TRUE}.
#' @param kink.position Optional boundary positions of known kink points.
#' Either a logical vector with one entry per boundary point, where \code{TRUE}
#' identifies a kink point, or an integer vector with indices between 1 and
#' the number of boundary points. Requires \code{b} so distances between
#' boundary points can be computed. This option applies only to automatic
#' bandwidth selection.
#' @param kernel Kernel function to use. Options are \code{"tri"} or
#' \code{"triangular"} (triangular, default), \code{"epa"} or
#' \code{"epanechnikov"} (Epanechnikov), \code{"uni"} or \code{"uniform"}
#' (uniform), and \code{"gau"} or \code{"gaussian"} (Gaussian).
#' @param level Nominal confidence level for intervals/bands, between 0 and 100 (default is 95).
#' @param cbands Logical. If \code{TRUE}, stores the covariance matrix needed
#' for uniform confidence bands and aggregate inference computed by
#' \code{summary(..., cbands = "main")}, \code{summary(..., WBATE = weights)},
#' or \code{summary(..., LBATE = TRUE)}. Default is \code{TRUE}.
#' @param side Type of confidence interval. Options: \code{"two"} (two-sided, default), \code{"left"} (left tail), or \code{"right"} (right tail).
#' @param repp Number of bootstrap repetitions used for critical value simulation. Default is \code{1000}.
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
#' \item \code{"user provided"}. User-provided bandwidths. If \code{h} is not \code{NULL}, then \code{bwselect} is overwritten to \code{"user provided"}.
#' }
#' @param bwparam Target parameter used for fuzzy automatic bandwidth
#' selection. Options are \code{"main"} (default), which selects bandwidths
#' for the linearized fuzzy Wald ratio, and \code{"itt"}, which selects
#' bandwidths for the reduced-form outcome. Ignored in sharp designs.
#' @param params.other Optional character vector requesting companion output
#' tables. In sharp designs, available values are \code{"main.0"} and
#' \code{"main.1"}. In fuzzy designs, available values are \code{"itt.0"},
#' \code{"itt.1"}, \code{"fs.0"}, and \code{"fs.1"}.
#' @param params.cov Optional character vector requesting covariance matrices
#' for aggregate inference or confidence bands. Available values are
#' \code{"main"} for sharp and fuzzy main effects, plus \code{"main.0"} and
#' \code{"main.1"} for sharp side-specific effects, and \code{"itt"},
#' \code{"itt.0"}, \code{"itt.1"}, \code{"fs"}, \code{"fs.0"}, and
#' \code{"fs.1"} for fuzzy reduced-form and first-stage effects. If
#' \code{cbands = TRUE}, \code{"main"} is added automatically.
#' @param vce Variance-covariance estimator for standard errors.
#' Options:
#' \describe{
#'   \item{\code{"hc0"}}{Heteroskedasticity-robust variance estimator without small sample adjustment (White robust).}
#'   \item{\code{"hc1"}}{Heteroskedasticity-robust variance estimator with degrees-of-freedom correction. (default)}
#'   \item{\code{"hc2"}}{Heteroskedasticity-robust variance estimator using leverage adjustments.}
#'   \item{\code{"hc3"}}{More conservative heteroskedasticity-robust variance estimator (similar to jackknife correction).}
#' }
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used. Default is \code{50 + p + 1}.
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
#' @param fuzzy Optional treatment receipt/status variable for fuzzy RD designs. If
#' supplied, \code{main} reports the fuzzy Wald treatment effect, \code{itt}
#' reports the reduced-form outcome discontinuity, and \code{fs} reports the
#' first-stage treatment receipt/status discontinuity.
#'
#' @return An object of class \code{"rd2d.distance"}, a list containing:
#' \describe{
#'   \item{\code{main}}{Data frame of point estimates, standard errors, confidence intervals, and bandwidths:
#'     \describe{
#'       \item{\code{b1}}{First coordinate of the evaluation point.}
#'       \item{\code{b2}}{Second coordinate of the evaluation point.}
#'       \item{\code{estimate.p}}{Point estimate with polynomial order \eqn{p}.}
#'       \item{\code{std.err.p}}{Standard error for \code{estimate.p}.}
#'       \item{\code{estimate.q}}{Bias-corrected estimate with polynomial order \eqn{q}.}
#'       \item{\code{std.err.q}}{Standard error for \code{estimate.q}.}
#'       \item{\code{t.value}}{t-statistic based on \eqn{\widehat{\tau}_{\text{distance},q}(\mathbf{b})}.}
#'       \item{\code{p.value}}{Two-sided p-value based on \eqn{T_{\text{distance},q}(\mathbf{b})}.}
#'       \item{\code{ci.lower}}{Lower bound of confidence interval.}
#'       \item{\code{ci.upper}}{Upper bound of confidence interval.}
#'       \item{\code{h0}}{Bandwidth used for the control group (negative signed distance).}
#'       \item{\code{h1}}{Bandwidth used for the treatment group (non-negative signed distance).}
#'       \item{\code{h0.rbc}}{Bandwidth used for control-side inference.}
#'       \item{\code{h1.rbc}}{Bandwidth used for treated-side inference.}
#'       \item{\code{N.Co}}{Effective sample size for the control side.}
#'       \item{\code{N.Tr}}{Effective sample size for the treatment side.}
#'     }
#'   }
#'   \item{\code{main.0}}{Sharp-design summary table for the control side only.}
#'   \item{\code{main.1}}{Sharp-design summary table for the treated side only.}
#'   \item{\code{itt}}{Fuzzy-design reduced-form outcome summary table.}
#'   \item{\code{itt.0}, \code{itt.1}}{Optional fuzzy-design side-specific
#'   reduced-form outcome summary tables requested through \code{params.other}.}
#'   \item{\code{fs}}{Fuzzy-design first-stage treatment receipt/status summary table.}
#'   \item{\code{fs.0}, \code{fs.1}}{Optional fuzzy-design side-specific
#'   first-stage treatment receipt/status summary tables requested through \code{params.other}.}
#'   \item{\code{bw}}{Bandwidth and effective-sample-size table.}
#'   \item{\code{tau.hat}}{Vector of point estimates \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'   \item{\code{tau.hat.q}}{Vector of bias-corrected estimates
#'   \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'   \item{\code{se.hat}}{Standard errors corresponding to \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'   \item{\code{se.hat.q}}{Standard errors corresponding to \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'   \item{\code{params.cov}}{List containing covariance matrices requested
#'   for aggregate or uniform inference.}
#'   \item{\code{cb}}{Pointwise confidence interval endpoints.}
#'   \item{\code{pvalues}}{Two-sided p-values based on bias-corrected
#'   estimates.}
#'   \item{\code{tvalues}}{t-statistics based on bias-corrected estimates.}
#'   \item{\code{tau.itt}, \code{tau.itt.q}, \code{tau.fs}, \code{tau.fs.q}}{Fuzzy
#'   reduced-form and first-stage point estimates with polynomial orders
#'   \eqn{p} and \eqn{q}; these are \code{NULL} in sharp designs.}
#'   \item{\code{rdmodel}}{Character label describing the fitted RD model.}
#'   \item{\code{call}}{Matched function call.}
#'   \item{\code{opt}}{A list of estimation options (e.g., \code{p}, \code{q}, \code{kernel}, \code{level}, etc.) and internal variables such as sample size \eqn{N}.}
#' }
#'
#' @seealso \code{\link{rdbw2d.distance}}, \code{\link{rd2d}}, \code{\link{print.rd2d.distance}}, \code{\link{summary.rd2d.distance}}
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
#' # Run the rd2d.distance function
#' result <- rd2d.distance(Y, distance = distance, b = eval, cbands = FALSE,
#'                     masspoints = "off", bwcheck = 10)
#'
#' # View the estimation results
#' print(result)
#' summary(result)
#'
#' # Fuzzy distance-based fit with a user-supplied bandwidth
#' fuzzy <- as.numeric(runif(n) < ifelse(assignment == 1, 0.8, 0.2))
#' Y.fuzzy <- 3 + 2 * x1 + 1.5 * x2 + 1.5 * fuzzy + rnorm(n, sd = 0.5)
#' fuzzy.result <- rd2d.distance(Y.fuzzy, distance = distance, h = 0.8, b = eval, fuzzy = fuzzy,
#'                           cbands = FALSE, masspoints = "off")
#' summary(fuzzy.result)
#' @export

rd2d.distance <- function(Y, distance, h = NULL, b = NULL, p = 1, q = NULL,
                      kink.unknown = c(FALSE, FALSE), kink.position = NULL,
                      kernel = c("tri","triangular", "epa","epanechnikov","uni","uniform","gau","gaussian"),
                      level = 95, cbands = TRUE, side = c("two", "left", "right"), repp = 1000,
                      bwselect = c("mserd", "cerrd", "imserd", "icerrd",
                                   "msetwo", "certwo", "imsetwo", "icertwo",
                                   "user provided"),
                      bwparam = c("main", "itt"),
                      params.other = NULL, params.cov = NULL,
                      vce = c("hc1","hc0","hc2","hc3"),
                      bwcheck = 50 + p + 1, masspoints = c("check","adjust","off"),
                      cluster = NULL, scaleregul = 1, cqt = 0.5, fuzzy = NULL){

  # Input error handling

  bwselect <- match.arg(bwselect)
  kernel <- match.arg(kernel)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  side <- match.arg(side)
  bwselect.base <- rd2d_bwselect_base(bwselect)
  bwselect.cer <- rd2d_bwselect_is_cer(bwselect)
  bwselect.common <- rd2d_bwselect_is_common(bwselect)
  kink.unknown <- rd2d_validate_kink_unknown(kink.unknown)
  bwparam <- match.arg(bwparam)
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
  q.user <- !is.null(q)

  valid.p <- is.numeric(p) && length(p) == 1 &&
    is.finite(p) && p >= 0 &&
    abs(p - round(p)) < sqrt(.Machine$double.eps)
  if (!valid.p) {
    print("p must be a nonnegative integer")
    exit <- 1
  } else {
    p <- as.integer(round(p))
  }

  if (!q.user) {
    q <- if (valid.p) {
      if (kink.unknown[1]) p else p + 1
    } else {
      NA_real_
    }
  }

  valid.q <- is.numeric(q) && length(q) == 1 &&
    is.finite(q) && q >= 0 &&
    abs(q - round(q)) < sqrt(.Machine$double.eps)
  if (!valid.q) {
    print("q must be a nonnegative integer")
    exit <- 1
  } else {
    q <- as.integer(round(q))
  }

  if (valid.p && valid.q && q < p) {
    print("q must be greater than or equal to p")
    exit <- 1
  }

  # level must be numeric in (0, 100)
  if (!is.numeric(level) || level <= 0 || level >= 100) {
    print("level must be a numeric value between 0 and 100")
    exit <- 1
  }

  # repp must be a positive integer
  if (!is.numeric(repp) || repp < 1 || repp != as.integer(repp)) {
    print("repp must be a positive integer")
    exit <- 1
  }

  # h must be either a positive scalar or a matrix/data.frame with same rows as b and 4 columns
  if (!is.null(h)) {
    if (length(h) == 1) {
      if (!is.numeric(h) || h <= 0) {
        print("If h is a scalar, it must be a positive numeric value")
        exit <- 1
      }
    } else if (!(is.matrix(h) || is.data.frame(h)) ||
               nrow(h) != ncol(distance_mat) || ncol(h) != 2) {
      print("If h is not a scalar, it must be a matrix or data frame with the same number of rows as b and 2 columns")
      exit <- 1
    }
  }

  if (is.null(h) & bwselect == "user provided"){
    exit <- 1
    print("Please provide bandwidths.")
  }

  allowed.params.other <- c(
    "itt.0", "itt.1", "fs.0", "fs.1", "main.0", "main.1"
  )
  allowed.params.cov <- c(
    "main", "main.0", "main.1", "itt", "itt.0", "itt.1",
    "fs", "fs.0", "fs.1"
  )
  if (is.null(params.other)) {
    params.other <- character(0)
  } else if (!is.character(params.other) || any(is.na(params.other))) {
    print("params.other must be NULL or a character vector")
    exit <- 1
  } else if (!all(params.other %in% allowed.params.other)) {
    print(paste("params.other must contain only:", paste(allowed.params.other, collapse = ", ")))
    exit <- 1
  } else {
    params.other <- unique(params.other)
  }

  if (is.null(params.cov)) {
    params.cov <- character(0)
  } else if (!is.character(params.cov) || any(is.na(params.cov))) {
    print("params.cov must be NULL or a character vector")
    exit <- 1
  } else if (!all(params.cov %in% allowed.params.cov)) {
    print(paste("params.cov must contain only:", paste(allowed.params.cov, collapse = ", ")))
    exit <- 1
  } else {
    params.cov <- unique(params.cov)
  }
  if (isTRUE(cbands)) params.cov <- unique(c(params.cov, "main"))

  if (exit == 0) {
    if (is.fuzzy) {
      invalid.design.params <- intersect(params.other, c("main.0", "main.1"))
      if (length(invalid.design.params) > 0) {
        print("main.0 and main.1 are available only for sharp rd2d.distance fits")
        exit <- 1
      }
      invalid.cov <- intersect(params.cov, c("main.0", "main.1"))
      if (length(invalid.cov) > 0) {
        print("main.0 and main.1 covariance matrices are available only for sharp rd2d.distance fits")
        exit <- 1
      }
      missing.side.cov <- setdiff(
        intersect(params.cov, c("itt.0", "itt.1", "fs.0", "fs.1")),
        params.other
      )
      if (length(missing.side.cov) > 0) {
        print(paste(
          "params.cov side outputs must also be requested in params.other:",
          paste(missing.side.cov, collapse = ", ")
        ))
        exit <- 1
      }
    } else {
      invalid.design.params <- setdiff(params.other, c("main.0", "main.1"))
      if (length(invalid.design.params) > 0) {
        print("itt.0, itt.1, fs.0, and fs.1 are available only for fuzzy rd2d.distance fits")
        exit <- 1
      }
      invalid.cov <- setdiff(params.cov, c("main", "main.0", "main.1"))
      if (length(invalid.cov) > 0) {
        print("itt, itt.0, itt.1, fs, fs.0, and fs.1 covariance matrices are available only for fuzzy rd2d.distance fits")
        exit <- 1
      }
    }
  }

  if (!is.null(cluster) && !(vce %in% c("hc0", "hc1"))) {
    warning("When cluster is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  ############################ Data Preparation ################################

  neval <- ncol(distance_mat)
  kink.position <- rd2d_validate_kink_position(kink.position, neval, b)
  if (any(kink.position) && kink.unknown[1]) {
    stop("Use either kink.position or kink.unknown, not both.", call. = FALSE)
  }
  if (!is.null(h) && any(kink.position)) {
    stop("kink.position applies only to automatic bandwidth selection; omit h to use it.", call. = FALSE)
  }
  if (!is.null(h) && any(kink.unknown)) {
    stop("kink.unknown applies only to automatic bandwidth selection; omit h to use it.", call. = FALSE)
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

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  ################################ Bandwidth ###################################

  if (is.null(h)){
    bws <- rdbw2d.distance(
      Y = Y, distance = distance_mat, b = b, p = p,
      kink.unknown = kink.unknown, kink.position = kink.position, kernel = kernel,
      bwselect = bwselect, bwparam = bwparam, vce = vce,
      bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
      scaleregul = scaleregul, cqt = cqt,
      fuzzy = if (is.fuzzy) fuzzy else NULL
    )
    bws <- bws$bws
    hgrid <- bws[,3]
    hgrid.1 <- bws[,4]
  } else {
    bwselect <- "user provided"
    bwselect.base <- "user provided"
    bwselect.cer <- FALSE
    bwselect.common <- FALSE
    # standardize bandwidth
    if (length(h) == 1){
      hgrid <- rep(h, neval)
      hgrid.1 <- rep(h, neval)
    } else {
      hgrid <- h[,1]
      hgrid.1 <- h[,2]
    }
  }

  hfull <- cbind(hgrid, hgrid.1)

  ###################### Point estimation and inference ########################

  need.p.cov.cache <- length(params.cov) > 0 && (q == p || (kink.unknown[2] && q > p))

  if (!is.fuzzy) {
    distancefit.p <- rd2d_distance_fit(
      Y = Y, distance = distance_mat, h = hfull, p = p, b = b, kernel = kernel,
      vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
      cbands = need.p.cov.cache
    )
    estimate.p <- distancefit.p$Estimate
    tau.hat.p <- estimate.p$mu1 - estimate.p$mu0
    se.hat.p <- sqrt(estimate.p$se0^2 + estimate.p$se1^2)
    eN0.p <- estimate.p$N0
    eN1.p <- estimate.p$N1
  } else {
    distancefit.Y.p <- rd2d_distance_fit(
      Y = Y, distance = distance_mat, h = hfull, p = p, b = b, kernel = kernel,
      vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
      cbands = TRUE
    )
    distancefit.fs.p <- rd2d_distance_fit(
      Y = fuzzy, distance = distance_mat, h = hfull, p = p, b = b, kernel = kernel,
      vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
      cbands = TRUE
    )
    estimate.itt.p <- distancefit.Y.p$Estimate
    estimate.fs.p <- distancefit.fs.p$Estimate
    tau.itt.p <- estimate.itt.p$mu1 - estimate.itt.p$mu0
    tau.fs.p <- estimate.fs.p$mu1 - estimate.fs.p$mu0
    se.itt.p <- sqrt(estimate.itt.p$se0^2 + estimate.itt.p$se1^2)
    se.fs.p <- sqrt(estimate.fs.p$se0^2 + estimate.fs.p$se1^2)
    denom.tol <- sqrt(.Machine$double.eps)
    valid.p <- is.finite(tau.fs.p) & abs(tau.fs.p) > denom.tol
    tau.hat.p <- rep(NA_real_, neval)
    tau.hat.p[valid.p] <- tau.itt.p[valid.p] / tau.fs.p[valid.p]
    cov.p <- rd2d_distance_fuzzy_cov_tables(
      distancefit.Y.p, distancefit.fs.p, p, tau.itt.p, tau.fs.p,
      outputs = "main", clustered = !is.null(cluster), denom.tol = denom.tol
    )
    se.hat.p <- sqrt(diag(cov.p$main))
    eN0.p <- estimate.itt.p$N0
    eN1.p <- estimate.itt.p$N1
  }

  fit.p.ref <- if (is.fuzzy) distancefit.Y.p else distancefit.p
  M.vec <- fit.p.ref$M.vec
  M.0.vec <- fit.p.ref$M.0.vec
  M.1.vec <- fit.p.ref$M.1.vec

  # inference bandwidth

  hfull.rbc <- hfull
  if (kink.unknown[2]){
    if (bwselect.common){
      hfull.rbc <- hfull * rd2d_bw_rate_factor(M.vec, 1 / 4, 1 / 3)
    } else {
      hfull.rbc[,1] <- hfull[,1] * rd2d_bw_rate_factor(M.0.vec, 1 / 4, 1 / 3)
      hfull.rbc[,2] <- hfull[,2] * rd2d_bw_rate_factor(M.1.vec, 1 / 4, 1 / 3)
    }
  }

  cov.tables <- list()
  clustered <- !is.null(cluster)
  reuse.q.fit <- q == p && identical(hfull.rbc, hfull)

  if (!is.fuzzy) {
    if (reuse.q.fit) {
      distancefit.q <- distancefit.p
    } else {
      distancefit.q <- rd2d_distance_fit(
        Y = Y, distance = distance_mat, h = hfull.rbc, p = q, b = b, kernel = kernel,
        vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
        cbands = length(params.cov) > 0
      )
    }
    estimate.q <- distancefit.q$Estimate
    tau.hat.q <- estimate.q$mu1 - estimate.q$mu0
    se.hat.q <- sqrt(estimate.q$se0^2 + estimate.q$se1^2)
    eN0.q <- estimate.q$N0
    eN1.q <- estimate.q$N1

    if (length(params.cov) > 0) {
      if ("main" %in% params.cov) {
        cov.tables$main <- rd2d_distance_cov_from_fits(
          distancefit.q, distancefit.q, q, clustered, side = "both"
        )
      }
      if ("main.0" %in% params.cov) {
        cov.tables$main.0 <- rd2d_distance_cov_from_fits(
          distancefit.q, distancefit.q, q, clustered, side = "0"
        )
      }
      if ("main.1" %in% params.cov) {
        cov.tables$main.1 <- rd2d_distance_cov_from_fits(
          distancefit.q, distancefit.q, q, clustered, side = "1"
        )
      }
    }
  } else {
    if (reuse.q.fit) {
      distancefit.Y.q <- distancefit.Y.p
      distancefit.fs.q <- distancefit.fs.p
    } else {
      distancefit.Y.q <- rd2d_distance_fit(
        Y = Y, distance = distance_mat, h = hfull.rbc, p = q, b = b, kernel = kernel,
        vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
        cbands = TRUE
      )
      distancefit.fs.q <- rd2d_distance_fit(
        Y = fuzzy, distance = distance_mat, h = hfull.rbc, p = q, b = b, kernel = kernel,
        vce = vce, bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
        cbands = TRUE
      )
    }
    estimate.itt.q <- distancefit.Y.q$Estimate
    estimate.fs.q <- distancefit.fs.q$Estimate
    tau.itt.q <- estimate.itt.q$mu1 - estimate.itt.q$mu0
    tau.fs.q <- estimate.fs.q$mu1 - estimate.fs.q$mu0
    se.itt.q <- sqrt(estimate.itt.q$se0^2 + estimate.itt.q$se1^2)
    se.fs.q <- sqrt(estimate.fs.q$se0^2 + estimate.fs.q$se1^2)

    valid.q <- is.finite(tau.fs.q) & abs(tau.fs.q) > denom.tol
    if (any(!valid.p | !valid.q)) {
      warning(
        paste(
          "Weak or zero first-stage fuzzy RD estimates detected;",
          "returning NA for affected fuzzy estimates."
        )
      )
    }

    tau.hat.q <- rep(NA_real_, neval)
    tau.hat.q[valid.q] <- tau.itt.q[valid.q] / tau.fs.q[valid.q]

    cov.outputs <- unique(c("main", params.cov))
    cov.q <- rd2d_distance_fuzzy_cov_tables(
      distancefit.Y.q, distancefit.fs.q, q, tau.itt.q, tau.fs.q,
      outputs = cov.outputs, clustered = clustered, denom.tol = denom.tol
    )
    se.hat.q <- sqrt(diag(cov.q$main))
    cov.tables <- cov.q[intersect(names(cov.q), params.cov)]
    eN0.q <- estimate.itt.q$N0
    eN1.q <- estimate.itt.q$N1
  }

  tvalues <- tau.hat.q/se.hat.q
  pvalues <- 2 * pnorm(abs(tvalues), lower.tail = FALSE)

  if (side == "two"){
    zval <- qnorm((level + 100)/ 200)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- tau.hat.q + zval * se.hat.q
  }
  if (side == "left"){
    zval <- qnorm(level / 100)
    CI.upper <- tau.hat.q + zval * se.hat.q
    CI.lower <- rep(-Inf, length(CI.upper))
  }
  if (side == "right"){
    zval <- qnorm(level / 100)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- rep(Inf, length(CI.lower))
  }

  cb.hat.q <- list(CI.l = CI.lower, CI.r = CI.upper)

  ############################### outputs ######################################

  main <- cbind(eval[,1], eval[,2], tau.hat.p, se.hat.p, tau.hat.q, se.hat.q, tvalues, pvalues,
                CI.lower, CI.upper, hfull[,1], hfull[,2],
                hfull.rbc[,1], hfull.rbc[,2], eN0.p, eN1.p)
  main <- as.data.frame(main)
  colnames(main) <- c("b1","b2","estimate.p","std.err.p","estimate.q","std.err.q",
                      "t.value", "p.value", "ci.lower","ci.upper",
                      "h0", "h1",
                      "h0.rbc", "h1.rbc", "N.Co", "N.Tr")
  rownames(main) <- NULL
  main.names <- colnames(main)

  make_result_table <- function(est.p, se.p, est.q, se.q,
                                h0 = hfull[, 1], h1 = hfull[, 2],
                                h0.rbc = hfull.rbc[, 1],
                                h1.rbc = hfull.rbc[, 2],
                                NCo = eN0.p, NTr = eN1.p) {
    t.comp <- est.q / se.q
    p.comp <- 2 * pnorm(abs(t.comp), lower.tail = FALSE)

    if (side == "two") {
      zval.comp <- qnorm((level + 100) / 200)
      ci.lower.comp <- est.q - zval.comp * se.q
      ci.upper.comp <- est.q + zval.comp * se.q
    }
    if (side == "left") {
      zval.comp <- qnorm(level / 100)
      ci.upper.comp <- est.q + zval.comp * se.q
      ci.lower.comp <- rep(-Inf, length(ci.upper.comp))
    }
    if (side == "right") {
      zval.comp <- qnorm(level / 100)
      ci.lower.comp <- est.q - zval.comp * se.q
      ci.upper.comp <- rep(Inf, length(ci.lower.comp))
    }

    out.comp <- data.frame(
      b1 = eval[, 1],
      b2 = eval[, 2],
      estimate.p = est.p,
      std.err.p = se.p,
      estimate.q = est.q,
      std.err.q = se.q,
      t.value = t.comp,
      p.value = p.comp,
      ci.lower = ci.lower.comp,
      ci.upper = ci.upper.comp,
      h0 = h0,
      h1 = h1,
      h0.rbc = h0.rbc,
      h1.rbc = h1.rbc,
      N.Co = NCo,
      N.Tr = NTr,
      check.names = FALSE
    )
    rownames(out.comp) <- NULL
    out.comp[, main.names]
  }

  if (!is.fuzzy) {
    main.0 <- make_result_table(
      estimate.p$mu0, estimate.p$se0, estimate.q$mu0, estimate.q$se0,
      h1 = rep(NA_real_, neval), h1.rbc = rep(NA_real_, neval),
      NTr = rep(NA_real_, neval)
    )
    main.1 <- make_result_table(
      estimate.p$mu1, estimate.p$se1, estimate.q$mu1, estimate.q$se1,
      h0 = rep(NA_real_, neval), h0.rbc = rep(NA_real_, neval),
      NCo = rep(NA_real_, neval)
    )
    itt <- NULL
    itt.0 <- NULL
    itt.1 <- NULL
    fs <- NULL
    fs.0 <- NULL
    fs.1 <- NULL
    tau.itt.p <- tau.itt.q <- tau.fs.p <- tau.fs.q <- NULL
  } else {
    itt <- make_result_table(tau.itt.p, se.itt.p, tau.itt.q, se.itt.q)
    fs <- make_result_table(tau.fs.p, se.fs.p, tau.fs.q, se.fs.q)

    itt.0 <- NA
    itt.1 <- NA
    fs.0 <- NA
    fs.1 <- NA

    if ("itt.0" %in% params.other) {
      itt.0 <- make_result_table(
        estimate.itt.p$mu0, estimate.itt.p$se0,
        estimate.itt.q$mu0, estimate.itt.q$se0,
        h1 = rep(NA_real_, neval), h1.rbc = rep(NA_real_, neval),
        NTr = rep(NA_real_, neval)
      )
    }
    if ("itt.1" %in% params.other) {
      itt.1 <- make_result_table(
        estimate.itt.p$mu1, estimate.itt.p$se1,
        estimate.itt.q$mu1, estimate.itt.q$se1,
        h0 = rep(NA_real_, neval), h0.rbc = rep(NA_real_, neval),
        NCo = rep(NA_real_, neval)
      )
    }
    if ("fs.0" %in% params.other) {
      fs.0 <- make_result_table(
        estimate.fs.p$mu0, estimate.fs.p$se0,
        estimate.fs.q$mu0, estimate.fs.q$se0,
        h1 = rep(NA_real_, neval), h1.rbc = rep(NA_real_, neval),
        NTr = rep(NA_real_, neval)
      )
    }
    if ("fs.1" %in% params.other) {
      fs.1 <- make_result_table(
        estimate.fs.p$mu1, estimate.fs.p$se1,
        estimate.fs.q$mu1, estimate.fs.q$se1,
        h0 = rep(NA_real_, neval), h0.rbc = rep(NA_real_, neval),
        NCo = rep(NA_real_, neval)
      )
    }
  }

  bw <- data.frame(
    b1 = eval[, 1],
    b2 = eval[, 2],
    h0 = hfull[, 1],
    h1 = hfull[, 2],
    N.Co = eN0.p,
    N.Tr = eN1.p,
    check.names = FALSE
  )
  rownames(bw) <- NULL

  rdmodel <- ifelse(is.fuzzy, "fuzzy rd2d.distance", "rd2d.distance")
  result.tables <- list(main = main, bw = bw)
  if (!is.fuzzy) {
    result.tables <- c(result.tables, list(main.0 = main.0, main.1 = main.1))
  } else {
    result.tables <- c(
      result.tables,
      list(
        itt = itt, itt.0 = itt.0, itt.1 = itt.1,
        fs = fs, fs.0 = fs.0, fs.1 = fs.1
      )
    )
  }

  out <- c(
    result.tables,
    list(
      opt=list(
        b = eval, p = p, q = q, kernel=kernel.type,
        kink.unknown = kink.unknown, kink.position = kink.position,
        N=N, N.0 = N.0, N.1 = N.1,
        M = M.vec, M.0 = M.0.vec, M.1 = M.1.vec, neval=neval,
        bwselect = bwselect, bwparam = bwparam, vce = vce,
        bwcheck = bwcheck, masspoints = masspoints, cluster = cluster,
        clustered = clustered, scaleregul = scaleregul, cqt = cqt,
        level = level, repp = repp, side = side, cbands = cbands,
        params.other = params.other, params.cov = params.cov,
        fuzzy = is.fuzzy,
        h0 = hfull[,1], h1 = hfull[,2],
        h0.rbc = hfull.rbc[,1], h1.rbc = hfull.rbc[,2],
        N.Co = eN0.p, N.Tr = eN1.p
      ),
      tau.hat = tau.hat.p, tau.hat.q = tau.hat.q,
      se.hat = se.hat.p, se.hat.q = se.hat.q,
      params.cov = cov.tables,
      cb = cb.hat.q, pvalues = pvalues, tvalues = tvalues,
      tau.itt = tau.itt.p, tau.itt.q = tau.itt.q,
      tau.fs = tau.fs.p, tau.fs.q = tau.fs.q,
      rdmodel = rdmodel
    )
  )
  out$call   <- match.call()
  class(out) <- "rd2d.distance"

  return(out)
}

################################################################################
#' Print Method for 2D Local Polynomial RD Estimation (Distance-Based)
#'
#' @description
#' Prints the results of a 2D local polynomial regression discontinuity (RD) estimation using distance-based evaluation, as obtained from \code{\link{rd2d.distance}}.
#'
#' @param x An object of class \code{rd2d.distance}, returned by \code{\link{rd2d.distance}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return
#' No return value. This function is called for its side effects: it prints the \code{\link{rd2d.distance}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso
#' \code{\link{rd2d.distance}} for estimation using distance-based methods in 2D local polynomial RD designs.
#'
#' Supported methods: \code{\link{print.rd2d.distance}}, \code{\link{summary.rd2d.distance}}.
#'
#' @export
#'
print.rd2d.distance <- function(x,...) {

  cat(paste(x$rdmodel, "\n", sep = ""))
  cat(paste("\n", sep = ""))

  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
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
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  cat("\n")
}

################################################################################
#' Summary Method for 2D Local Polynomial RD Estimation (Distance-Based)
#'
#' @description
#' Summarizes estimation and bandwidth results from a 2D local polynomial regression discontinuity (RD) design using distance-based methods, as returned by \code{\link{rd2d.distance}}.
#'
#' @param object An object of class \code{rd2d.distance}, returned by \code{\link{rd2d.distance}}.
#' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{cbands}: Character vector. Use \code{cbands = "main"} to
#'       display uniform confidence bands for the main distance-based estimates.
#'       Other stored outputs can be requested when their covariance matrices
#'       were stored through \code{params.cov}. The default displays pointwise
#'       confidence intervals.
#'     \item \code{WBATE}: Optional numeric weights for a weighted boundary
#'       average treatment effect row. The weights must match the full set of
#'       evaluation points and are normalized internally. The fitted object must
#'       contain the covariance matrix stored by \code{rd2d.distance(cbands = TRUE)}
#'       or requested through \code{params.cov}.
#'     \item \code{LBATE}: Logical. If \code{TRUE}, prints a largest boundary
#'       average treatment effect row. The fitted object must contain the
#'       covariance matrix stored by \code{rd2d.distance(cbands = TRUE)} or
#'       requested through \code{params.cov}.
#'     \item \code{subset}: Integer vector of indices of evaluation points to display.
#'       Defaults to all evaluation points.
#'     \item \code{output}: Character vector. Use \code{"main"} to display
#'       treatment effect estimates or \code{"bw"} to display bandwidth
#'       information. In sharp designs, \code{"main.0"} and \code{"main.1"}
#'       display side-specific estimates. In fuzzy designs, \code{"itt"} and
#'       \code{"fs"} display reduced-form and first-stage estimates, and
#'       \code{"itt.0"}, \code{"itt.1"}, \code{"fs.0"}, and \code{"fs.1"}
#'       display requested side-specific companion estimates.
#'     \item \code{sep}: Integer vector of length three. Controls spacing in the output.
#'       \code{sep[1]} controls spacing for the columns of boundary points, estimation,
#'       t-value, and p-value in the \code{"main"} table.
#'       \code{sep[2]} controls spacing for confidence intervals (or bands) in the \code{"main"} table.
#'       \code{sep[3]} controls spacing for the columns in the \code{"bw"} table.
#'       Default is \code{c(7, 17, 8)}.
#'   }
#'
#' @return Invisibly returns an object of class \code{"summary.rd2d.distance"}, a
#'   list with elements:
#'   \itemize{
#'     \item \code{tables}: named list of returned summary tables.
#'     \item \code{cbands}: named list of confidence-band endpoints for outputs
#'       requested through \code{cbands}.
#'     \item \code{outputs}: character vector of summarized outputs.
#'     \item \code{call}: matched summary call.
#'   }
#'   Requested WBATE and LBATE rows are appended to the returned estimation
#'   table(s), except for \code{output = "bw"}.
#'   The function is also called for its side effect of printing a formatted
#'   summary.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso \code{\link{rd2d.distance}} for estimation using distance-based 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rd2d.distance}}, \code{\link{summary.rd2d.distance}}.
#'
#' @export

summary.rd2d.distance <- function(object, ...) {

  x <- object

  args <- list(...)
  valid.summary.args <- c("cbands", "WBATE", "LBATE", "subset", "output", "sep")
  arg.names <- names(args)
  if (length(args) > 0 && (is.null(arg.names) || any(arg.names == ""))) {
    stop("All summary.rd2d.distance options must be named.", call. = FALSE)
  }
  invalid.args <- setdiff(arg.names, valid.summary.args)
  if (length(invalid.args) > 0) {
    stop(
      sprintf(
        "Unsupported summary.rd2d.distance option(s): %s.",
        paste(invalid.args, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (is.null(args[["cbands"]])) {
    cbands <- character(0)
  } else {
    cbands <- args[["cbands"]]
    if (!is.character(cbands) || any(is.na(cbands))) {
      stop("cbands must be NULL or a character vector.", call. = FALSE)
    }
    cbands <- unique(cbands)
  }
  WBATE <- args[["WBATE"]]
  LBATE <- isTRUE(args[["LBATE"]])

  if (is.null(args[["subset"]])) {
    subset <- NULL
  } else {
    subset <- args[["subset"]]
  }

  output.requested <- !is.null(args[["output"]])
  if (!output.requested) {
    output <- "main"
  } else {
    output <- args[["output"]]
  }

  if (is.null(args[["sep"]])) {
    sep <- c(7,17,8)
  } else {
    sep <- args[["sep"]]
  }

  valid.outputs.all <- c("main", "bw")
  if (isTRUE(x$opt$fuzzy)) {
    valid.outputs.all <- c(
      valid.outputs.all, "itt", "itt.0", "itt.1", "fs", "fs.0", "fs.1"
    )
  } else {
    valid.outputs.all <- c(valid.outputs.all, "main.0", "main.1")
  }
  valid.cband.outputs <- setdiff(valid.outputs.all, "bw")
  invalid.cbands <- setdiff(cbands, valid.cband.outputs)
  if (length(invalid.cbands) > 0) {
    stop(
      sprintf(
        "cbands contains unavailable output(s): %s.",
        paste(invalid.cbands, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  cat(paste(x$rdmodel, "\n", sep = ""))
  cat(paste("\n", sep = ""))

  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("Unknown kink           %s\n", paste(tolower(x$opt$kink.unknown), collapse = ", ")))
  if (any(x$opt$kink.position)) {
    cat(sprintf("Known kink positions   %s\n", paste(which(x$opt$kink.position), collapse = ", ")))
  }
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""), sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  if (isTRUE(x$opt$fuzzy)) {
    cat(sprintf("Fuzzy                  on\n"))
    cat(sprintf("BW target              %s\n", x$opt$bwparam))
  }
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  cat("\n")

  eval.specified <- !all(is.na(x$opt$b))

  covariance_available <- function(cov.mat, n) {
    is.matrix(cov.mat) &&
      identical(dim(cov.mat), c(n, n)) &&
      all(is.finite(cov.mat)) &&
      all(is.finite(diag(cov.mat))) &&
      all(diag(cov.mat) > 0)
  }

  print_table <- function(output) {
    if (!(output %in% valid.outputs.all)) {
      warning(
        paste(
          sprintf("output='%s' is not available for this rd2d.distance object.", output),
          "Resetting to output='main'."
        )
      )
      output <- "main"
    }

    if (!is.data.frame(x[[output]])) {
      warning(
        paste(
          sprintf("output='%s' was not computed for this rd2d.distance object.", output),
          "Request companion outputs with params.other in rd2d.distance().",
          "Resetting to output='main'."
        )
      )
      output <- "main"
    }

    WBATE.now <- WBATE
    LBATE.now <- LBATE
    if (output == "bw" && (!is.null(WBATE.now) || LBATE.now)) {
      warning("WBATE and LBATE are not available for output='bw'.")
      WBATE.now <- NULL
      LBATE.now <- FALSE
    }

    results <- x[[output]]
    table.title <- switch(
      output,
      main = ifelse(isTRUE(x$opt$fuzzy), "Treatment effect estimates.",
                    "Sharp main estimates."),
      main.0 = "Sharp control-side estimates.",
      main.1 = "Sharp treatment-side estimates.",
      itt = "Intention-to-treat estimates.",
      itt.0 = "Intention-to-treat control-side estimates.",
      itt.1 = "Intention-to-treat treatment-side estimates.",
      fs = "First-stage estimates.",
      fs.0 = "First-stage control-side estimates.",
      fs.1 = "First-stage treatment-side estimates.",
      bw = "Bandwidth information."
    )
    cat(sprintf("%s\n\n", table.title))

    if (output != "bw") {
      bands.requested <- output %in% cbands
      neval <- nrow(results)
      if (is.null(subset)) {
        subset.now <- seq_len(neval)
      } else if (!all(subset %in% seq_len(neval))) {
        warning("Invalid subset provided. Resetting to default: 1:neval")
        subset.now <- seq_len(neval)
      } else {
        subset.now <- subset
      }

      cov.table <- NULL
      if (!is.null(x$params.cov) && is.matrix(x$params.cov[[output]])) {
        cov.table <- x$params.cov[[output]]
      }
      if (!is.null(cov.table) &&
          !identical(dim(cov.table), c(nrow(results), nrow(results)))) {
        cov.table <- NULL
      }
      summary.side <- x$opt$side

      require_covariance <- function(reason) {
        if (!covariance_available(cov.table, nrow(results))) {
          verb <- if (identical(reason, "Uniform confidence bands")) {
            "require"
          } else {
            "requires"
          }
          stop(
            paste(
              sprintf(
                "%s %s a stored covariance matrix for output='%s'.",
                reason, verb, output
              ),
              sprintf("Rerun rd2d.distance(..., params.cov = \"%s\").", output)
            ),
            call. = FALSE
          )
        }
      }

      if (bands.requested) require_covariance("Uniform confidence bands")
      if (!is.null(WBATE.now)) require_covariance("WBATE inference")
      if (LBATE.now) require_covariance("LBATE inference")

      result.subset <- results[subset.now,,drop = FALSE]
      interval.lower <- results$ci.lower
      interval.upper <- results$ci.upper
      cbands.return <- NULL
      if (bands.requested) {
        cb.out <- rd2d_cb(
          results$estimate.q, cov.table, x$opt$repp, x$opt$side, x$opt$level
        )
        interval.lower <- cb.out$CB.l
        interval.upper <- cb.out$CB.r
        result.subset$cb.lower <- cb.out$CB.l[subset.now]
        result.subset$cb.upper <- cb.out$CB.r[subset.now]
        cbands.return <- result.subset[, c("cb.lower", "cb.upper"), drop = FALSE]
      }

      if (eval.specified) {
        headers <- c("ID", "b1", "b2", "Est.", "t", "P > |t|", sprintf("%d%% CI", x$opt$level))
        col_widths <- c(4, sep[1], sep[1], sep[1], sep[1], sep[1], sep[2])
      } else {
        headers <- c("ID", "Est.", "t", "P > |t|", sprintf("%d%% CI", x$opt$level))
        col_widths <- c(4, sep[1], sep[1], sep[1], sep[2])
      }
      if (bands.requested) {
        headers[length(headers)] <- sprintf("%d%% Unif. CB", x$opt$level)
      }

      rule.width <- sum(col_widths) + 2 * (length(headers) - 1)
      rule <- paste(rep("=", rule.width), collapse = "")
      cat(rule, "\n")
      formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
      cat(paste(formatted_headers, collapse = "  "), "\n")
      cat(rule, "\n")

      format_blank <- function(width) {
        formatC("", width = width, format = "s")
      }

      format_number <- function(value, width, digits = 4) {
        if (length(value) != 1 || is.na(value)) {
          format_blank(width)
        } else {
          formatC(value, format = "f", digits = digits, width = width)
        }
      }

      format_interval <- function(lower, upper, width) {
        if (length(lower) != 1 || length(upper) != 1 ||
            is.na(lower) || is.na(upper)) {
          format_blank(width)
        } else {
          formatC(
            paste0(
              "[",
              formatC(lower, format = "f", digits = 4),
              ", ",
              formatC(upper, format = "f", digits = 4),
              "]"
            ),
            width = width,
            format = "s"
          )
        }
      }

      aggregate_interval <- function(center, se, cval) {
        if (summary.side == "two") {
          return(c(center - cval * se, center + cval * se))
        }
        if (summary.side == "left") {
          return(c(-Inf, center + cval * se))
        }
        c(center - cval * se, Inf)
      }

      print_aggregate_row <- function(label, est, t.value = NA_real_,
                                      p.value = NA_real_,
                                      lower = NA_real_, upper = NA_real_) {
        if (eval.specified) {
          row_vals <- c(
            formatC(label, width = col_widths[1], format = "s"),
            format_blank(col_widths[2]),
            format_blank(col_widths[3]),
            format_number(est, col_widths[4]),
            format_number(t.value, col_widths[5]),
            format_number(p.value, col_widths[6]),
            format_interval(lower, upper, col_widths[7])
          )
        } else {
          row_vals <- c(
            formatC(label, width = col_widths[1], format = "s"),
            format_number(est, col_widths[2]),
            format_number(t.value, col_widths[3]),
            format_number(p.value, col_widths[4]),
            format_interval(lower, upper, col_widths[5])
          )
        }
        cat(paste(row_vals, collapse = "  "), "\n")
      }

      make_wbate <- function(weights, result.all, cov.mat) {
        if (is.null(weights)) return(NULL)
        if (!is.numeric(weights) || any(!is.finite(weights))) {
          warning("WBATE must be a finite numeric vector. Ignoring WBATE.")
          return(NULL)
        }
        if (length(weights) != nrow(result.all)) {
          warning("WBATE length must match the number of evaluation points.")
          return(NULL)
        }
        if (sum(weights) == 0) {
          warning("WBATE weights must have a nonzero sum. Ignoring WBATE.")
          return(NULL)
        }

        weights <- weights / sum(weights)
        est <- sum(weights * result.all[["estimate.p"]])
        out <- list(label = "WBATE", est = est)

        if (covariance_available(cov.mat, nrow(result.all)) &&
            all(is.finite(result.all[["estimate.q"]]))) {
          se <- sqrt(as.numeric(t(weights) %*% cov.mat %*% weights))
          if (is.finite(se) && se > 0) {
            center <- sum(weights * result.all[["estimate.q"]])
            t.value <- center / se
            p.value <- 2 * pnorm(abs(t.value), lower.tail = FALSE)
            cval <- if (summary.side == "two") {
              qnorm((x$opt$level + 100) / 200)
            } else {
              qnorm(x$opt$level / 100)
            }
            ci <- aggregate_interval(center, se, cval)
            out$center <- center
            out$se <- se
            out$t.value <- t.value
            out$p.value <- p.value
            out$lower <- ci[1]
            out$upper <- ci[2]
          }
        }

        out
      }

      make_lbate <- function(result.all, cov.mat) {
        if (!LBATE.now) return(NULL)
        est <- max(result.all[["estimate.p"]], na.rm = TRUE)
        center <- max(result.all[["estimate.q"]], na.rm = TRUE)
        out <- list(label = "LBATE", est = est, center = center)

        if (covariance_available(cov.mat, nrow(result.all)) &&
            all(is.finite(result.all[["estimate.q"]]))) {
          se <- sqrt(diag(cov.mat))
          cval <- rd2d_cval(cov.mat, rep = x$opt$repp, side = summary.side,
                            alpha = x$opt$level, lp = Inf)
          if (summary.side == "two") {
            out$lower <- max(result.all[["estimate.q"]] - cval * se)
            out$upper <- max(result.all[["estimate.q"]] + cval * se)
          } else if (summary.side == "left") {
            out$lower <- -Inf
            out$upper <- max(result.all[["estimate.q"]] + cval * se)
          } else {
            out$lower <- max(result.all[["estimate.q"]] - cval * se)
            out$upper <- Inf
          }
        }

        out
      }

      for (i in seq_len(neval)) {
        if (i %in% subset.now) {
          if (eval.specified) {
            row_vals <- c(
              formatC(i, width = col_widths[1], format = "d"),
              formatC(results$b1[i], format = "f", digits = 3, width = col_widths[2]),
              formatC(results$b2[i], format = "f", digits = 3, width = col_widths[3]),
              formatC(results$estimate.p[i], format = "f", digits = 4, width = col_widths[4]),
              formatC(results$t.value[i], format = "f", digits = 4, width = col_widths[5]),
              formatC(results$p.value[i], format = "f", digits = 4, width = col_widths[6]),
              formatC(
                paste0("[", formatC(interval.lower[i], format = "f", digits = 4),
                       ", ", formatC(interval.upper[i], format = "f", digits = 4), "]"),
                width = col_widths[7], format = "s"
              )
            )
          } else {
            row_vals <- c(
              formatC(i, width = col_widths[1], format = "d"),
              formatC(results$estimate.p[i], format = "f", digits = 4, width = col_widths[2]),
              formatC(results$t.value[i], format = "f", digits = 4, width = col_widths[3]),
              formatC(results$p.value[i], format = "f", digits = 4, width = col_widths[4]),
              formatC(
                paste0("[", formatC(interval.lower[i], format = "f", digits = 4),
                       ", ", formatC(interval.upper[i], format = "f", digits = 4), "]"),
                width = col_widths[5], format = "s"
              )
            )
          }
          cat(paste(row_vals, collapse = "  "), "\n")
        }
      }

      wbate.row <- make_wbate(WBATE.now, results, cov.table)
      lbate.row <- make_lbate(results, cov.table)
      aggregate.rows <- Filter(Negate(is.null), list(wbate.row, lbate.row))
      summary.table <- result.subset
      rownames(summary.table) <- as.character(subset.now)

      make_aggregate_return_row <- function(aggregate.row) {
        row <- as.list(rep(NA_real_, ncol(summary.table)))
        names(row) <- names(summary.table)
        row[["estimate.p"]] <- aggregate.row$est
        row[["estimate.q"]] <- if (!is.null(aggregate.row$center)) {
          aggregate.row$center
        } else {
          NA_real_
        }
        row[["std.err.q"]] <- if (!is.null(aggregate.row$se)) {
          aggregate.row$se
        } else {
          NA_real_
        }
        row[["t.value"]] <- if (!is.null(aggregate.row$t.value)) {
          aggregate.row$t.value
        } else {
          NA_real_
        }
        row[["p.value"]] <- if (!is.null(aggregate.row$p.value)) {
          aggregate.row$p.value
        } else {
          NA_real_
        }
        row$ci.lower <- if (!is.null(aggregate.row$lower)) {
          aggregate.row$lower
        } else {
          NA_real_
        }
        row$ci.upper <- if (!is.null(aggregate.row$upper)) {
          aggregate.row$upper
        } else {
          NA_real_
        }
        as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
      }

      if (length(aggregate.rows) > 0) {
        aggregate.table <- do.call(
          rbind,
          lapply(aggregate.rows, make_aggregate_return_row)
        )
        rownames(aggregate.table) <- vapply(
          aggregate.rows,
          function(aggregate.row) aggregate.row$label,
          character(1)
        )
        summary.table <- rbind(summary.table, aggregate.table)

        for (aggregate.row in aggregate.rows) {
          cat(paste(rep("-", rule.width), collapse = ""), "\n")
          print_aggregate_row(
            aggregate.row$label, aggregate.row$est,
            aggregate.row$t.value, aggregate.row$p.value,
            aggregate.row$lower, aggregate.row$upper
          )
        }
      }

      cat(rule, "\n")
      return(list(output = output, table = summary.table, cbands = cbands.return))
    }

    results <- x$bw
    neval <- nrow(results)
    if (is.null(subset)) {
      subset.now <- seq_len(neval)
    } else if (!all(subset %in% seq_len(neval))) {
      warning("Invalid subset provided. Resetting to default: 1:neval")
      subset.now <- seq_len(neval)
    } else {
      subset.now <- subset
    }

    if (eval.specified) {
      headers <- c("ID", "b1", "b2", "h0", "h1", "N.Co", "N.Tr")
      col_widths <- c(4, sep[3], sep[3], sep[3], sep[3], sep[3], sep[3])
    } else {
      headers <- c("ID", "h0", "h1", "N.Co", "N.Tr")
      col_widths <- c(4, sep[3], sep[3], sep[3], sep[3])
    }

    cat(strrep("=", sum(col_widths)), "\n")
    if (eval.specified) {
      group_headers <- c(
        formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
        formatC("      Bandwidths", width = col_widths[4] + col_widths[5], format = "s", flag = "-"),
        formatC("      Eff. N", width = col_widths[6] + col_widths[7], format = "s", flag = "-")
      )
    } else {
      group_headers <- c(
        formatC("  ", width = col_widths[1], format = "s", flag = "-"),
        formatC("   Bandwidths", width = col_widths[2] + col_widths[3], format = "s", flag = "-"),
        formatC("      Eff. N", width = col_widths[4] + col_widths[5], format = "s", flag = "-")
      )
    }
    cat(paste(group_headers, collapse = ""), "\n")

    formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
    cat(paste(formatted_headers, collapse = ""), "\n")
    cat(strrep("=", sum(col_widths)), "\n")

    for (j in seq_len(neval)) {
      if (j %in% subset.now) {
        if (eval.specified) {
          row_vals <- c(
            formatC(j, width = col_widths[1], format = "d"),
            formatC(results$b1[j], format = "f", digits = 3, width = col_widths[2]),
            formatC(results$b2[j], format = "f", digits = 3, width = col_widths[3]),
            formatC(results$h0[j], format = "f", digits = 3, width = col_widths[4]),
            formatC(results$h1[j], format = "f", digits = 3, width = col_widths[5]),
            formatC(results$N.Co[j], format = "d", width = col_widths[6]),
            formatC(results$N.Tr[j], format = "d", width = col_widths[7])
          )
        } else {
          row_vals <- c(
            formatC(j, width = col_widths[1], format = "d"),
            formatC(results$h0[j], format = "f", digits = 3, width = col_widths[2]),
            formatC(results$h1[j], format = "f", digits = 3, width = col_widths[3]),
            formatC(results$N.Co[j], format = "d", width = col_widths[4]),
            formatC(results$N.Tr[j], format = "d", width = col_widths[5])
          )
        }
        cat(paste(row_vals, collapse = ""), "\n")
      }
    }
    cat(strrep("=", sum(col_widths)), "\n")
    list(output = output, table = results[subset.now,,drop = FALSE], cbands = NULL)
  }

  if (!output.requested) {
    if (isTRUE(x$opt$fuzzy)) {
      summary.outputs <- c(
        "main", "itt", "fs",
        intersect(x$opt$params.other, c("itt.0", "itt.1", "fs.0", "fs.1"))
      )
    } else {
      summary.outputs <- c(
        "main",
        intersect(x$opt$params.other, c("main.0", "main.1"))
      )
    }
    summary.outputs <- summary.outputs[vapply(
      summary.outputs,
      function(out) is.data.frame(x[[out]]),
      logical(1)
    )]
  } else {
    summary.outputs <- output
  }
  summary.result <- list(
    tables = list(),
    cbands = list(),
    outputs = summary.outputs,
    call = match.call()
  )
  for (out in summary.outputs) {
    printed <- print_table(out)
    summary.result$tables[[printed$output]] <- printed$table
    if (!is.null(printed$cbands)) {
      summary.result$cbands[[printed$output]] <- printed$cbands
    }
  }
  summary.result$outputs <- names(summary.result$tables)
  class(summary.result) <- "summary.rd2d.distance"
  invisible(summary.result)
}
