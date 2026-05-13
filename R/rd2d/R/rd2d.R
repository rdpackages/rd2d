################################################################################
#' @title Location-Based Methods for Boundary Discontinuity Design
#'
#' @description
#' \code{rd2d} implements location-based local polynomial boundary
#' discontinuity (BD) point estimators with robust bias-corrected pointwise
#' confidence intervals and uniform confidence bands.
#' See \href{https://arxiv.org/abs/2505.05670}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026)}
#' for methodological background.
#'
#' Companion commands are: \code{rdbw2d} for data-driven bandwidth selection.
#'
#' For other packages of RD designs, visit
#' <https://rdpackages.github.io/>
#'
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where
#' \eqn{N} is the sample size.
#' @param X Bivariate running variable (a.k.a score variable); a numeric
#' matrix or data frame of dimension \eqn{N \times 2}, with each row
#' \eqn{\mathbf{X}_i = (X_{1i}, X_{2i})}.
#' @param assignment Treatment assignment indicator; a logical or binary vector
#' indicating assignment to the treated side.
#' @param b Evaluation points; a matrix or data frame specifying boundary
#' points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, of dimension
#' \eqn{J \times 2}.
#' @param h Bandwidths. Either a positive scalar (same bandwidth for all
#' dimensions and groups), or a matrix/data frame of size \eqn{J \times 4},
#' corresponding to \eqn{h_{\text{control},1}},
#' \eqn{h_{\text{control},2}}, \eqn{h_{\text{treated},1}},
#' \eqn{h_{\text{treated},2}} at each evaluation point. If not specified,
#' bandwidth \code{h} is computed by the companion command
#' \code{rdbw2d}; for fuzzy designs, \code{bwparam} determines whether
#' automatic bandwidths target the fuzzy Wald ratio or the reduced-form
#' outcome discontinuity. Default is \code{h = NULL}.
#' @param fuzzy Optional treatment receipt/status variable used for fuzzy RD
#' estimation. The default is \code{NULL}, which estimates the sharp RD design.
#' @param deriv The order of the derivatives of the regression functions to be
#' estimated; a nonnegative integer vector of length 2 specifying the number
#' of derivatives in each coordinate (e.g., \eqn{c(1,2)} corresponds to
#' \eqn{\partial_1 \partial_2^2}).
#' @param tangvec Tangent vectors; a matrix or data frame of dimension
#' \eqn{J \times 2} specifying directional derivatives. Overrides
#' \code{deriv} if provided.
#' @param p Polynomial order for point estimation (\eqn{p = 1} by default).
#' @param q Polynomial order for robust confidence interval construction. Must
#' satisfy \eqn{q \geq p}; default is \eqn{q = p + 1}.
#' @param kernel Kernel function to use. Options are \code{"tri"} or
#' \code{"triangular"} (triangular, default), \code{"epa"} or
#' \code{"epanechnikov"} (Epanechnikov), \code{"uni"} or \code{"uniform"}
#' (uniform), and \code{"gau"} or \code{"gaussian"} (Gaussian).
#' @param kernel_type Kernel structure. Either \code{"prod"} for product
#' kernels (default) or \code{"rad"} for radial kernels.
#' @param vce Variance-covariance estimation method. Options are:
#' \itemize{
#' \item \code{"hc0"}: heteroskedasticity-robust plug-in residual variance
#' estimator without small-sample adjustment.
#' \item \code{"hc1"}: heteroskedasticity-robust plug-in residual variance
#' estimator with HC1 small-sample adjustment (default).
#' \item \code{"hc2"}: heteroskedasticity-robust plug-in residual variance
#' estimator with HC2 adjustment.
#' \item \code{"hc3"}: heteroskedasticity-robust plug-in residual variance
#' estimator with HC3 adjustment.
#' }
#' Default is \code{"hc1"}.
#' @param masspoints Handling of mass points in the running variable. Options
#' are:
#' \itemize{
#' \item \code{"check"}: detects presence of mass points and reports the
#' number of unique observations (default).
#' \item \code{"adjust"}: adjusts preliminary bandwidths to ensure a minimum
#' number of unique observations within each side of the cutoff.
#' \item \code{"off"}: ignores presence of mass points.
#' }
#' @param cluster Cluster ID variable used for cluster-robust variance estimation
#' with degrees-of-freedom weights. Default is \code{cluster = NULL}.
#' @param level Nominal confidence level for intervals/bands, between 0 and
#' 100 (default is 95).
#' @param params.other Optional character vector requesting companion-side
#' result tables. Options are \code{"main.0"}, \code{"main.1"} for sharp
#' fits, and \code{"itt.0"}, \code{"itt.1"}, \code{"fs.0"},
#' \code{"fs.1"} for fuzzy fits. Unrequested companion tables are returned
#' as \code{NA}. Default is \code{NULL}.
#' @param params.cov Optional character vector requesting covariance matrices
#' to store for later use by \code{\link{summary.rd2d}}. Options are
#' \code{"main"}, \code{"main.0"}, \code{"main.1"}, \code{"itt"},
#' \code{"itt.0"}, \code{"itt.1"}, \code{"fs"}, \code{"fs.0"}, and
#' \code{"fs.1"}, subject to the sharp/fuzzy design and
#' \code{params.other}. Default is \code{NULL}, which computes no covariance
#' matrices.
#' @param side Type of confidence interval. Options: \code{"two"} (two-sided,
#' default), \code{"left"} (left tail), or \code{"right"} (right tail).
#' @param repp Number of repetitions for critical value simulation (used in
#' uniform confidence bands). Default is 1000.
#' @param bwselect Bandwidth selection strategy. Options:
#' \itemize{
#' \item \code{"mserd"}. One common MSE-optimal bandwidth selector for the
#' boundary RD treatment effect estimator for each evaluation point (default).
#' \item \code{"cerrd"}. CER-optimal counterpart of \code{"mserd"}.
#' \item \code{"imserd"}. IMSE-optimal bandwidth selector for the boundary RD
#' treatment effect estimator based on all evaluation points.
#' \item \code{"icerrd"}. CER-optimal counterpart of \code{"imserd"}.
#' \item \code{"msetwo"}. Two different MSE-optimal bandwidth selectors
#' (control and treatment) for the boundary RD treatment effect estimator for
#' each evaluation point.
#' \item \code{"certwo"}. CER-optimal counterpart of \code{"msetwo"}.
#' \item \code{"imsetwo"}. Two IMSE-optimal bandwidth selectors (control and
#' treatment) for the boundary RD treatment effect estimator based on all
#' evaluation points.
#' \item \code{"icertwo"}. CER-optimal counterpart of \code{"imsetwo"}.
#' \item \code{"user provided"}. User-provided bandwidths. If \code{h} is not
#' \code{NULL}, then \code{bwselect} is overwritten to \code{"user provided"}.
#' }
#' @param bwparam Target parameter used for fuzzy automatic bandwidth
#' selection. Options are \code{"main"} (default), which selects bandwidths for
#' the linearized fuzzy Wald ratio, and \code{"itt"}, which selects bandwidths
#' for the reduced-form outcome discontinuity. This option is ignored in sharp
#' designs and when \code{h} is supplied.
#' @param method Bandwidth selection method for bias estimator based on local
#' polynomials. Either \code{"dpi"} (default) for data-driven plug-in MSE
#' optimal bandwidth selector or \code{"rot"} for rule-of-thumb bandwidth
#' selector.
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth
#' used in the calculations is enlarged so that at least \code{bwcheck} unique
#' observations are used. Default is \code{50 + p + 1}.
#' @param scaleregul Scaling factor for the regularization term in bandwidth
#' selection. Default is 3.
#' @param scalebiascrct Scaling factor used for bias correction based on higher
#' order expansions. Default is 1.
#' @param stdvars Logical. If \code{TRUE}, the running variables \eqn{X_{1i}}
#' and \eqn{X_{2i}} are standardized before computing automatic bandwidths.
#' Default is \code{TRUE}. Standardization only affects automatic bandwidth
#' selection when bandwidths are not manually provided through \code{h}.
#'
#' @return An object of class \code{"rd2d"}, a list with components:
#' \describe{
#'   \item{\code{main}}{A data frame with point estimates, standard errors,
#'   t-statistics, p-values, confidence intervals, and
#'   bandwidths at each evaluation point.
#'     \describe{
#'       \item{\code{b1}, \code{b2}}{First and second coordinate of
#'       evaluation points \eqn{\mathbf{b} = (b_1,b_2)}.}
#'       \item{\code{estimate.p}}{Point estimate \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'       \item{\code{std.err.p}}{Standard error of
#'       \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'       \item{\code{estimate.q}}{Bias-corrected point estimate
#'       \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'       \item{\code{std.err.q}}{Standard error of the bias-corrected estimate
#'       \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'       \item{\code{t.value}, \code{p.value}}{t-statistic and p-value based
#'       on the bias-corrected estimate.}
#'       \item{\code{ci.lower}, \code{ci.upper}}{Pointwise confidence
#'       intervals.}
#'       \item{\code{h01}, \code{h02}, \code{h11}, \code{h12}}{Bandwidths
#'       used in each coordinate and group.}
#'       \item{\code{N.Co}, \code{N.Tr}}{Effective sample size for the
#'       control and treatment sides, respectively.}
#'     }
#'   }
#'   \item{\code{main.0}, \code{main.1}}{For sharp RD only, companion-side
#'   outcome tables requested through \code{params.other}; otherwise
#'   \code{NA}.}
#'   \item{\code{bw}}{Bandwidth and effective sample size table.}
#'   \item{\code{itt}}{For fuzzy RD only, same structure as \code{main} but for
#'   the reduced-form or intention-to-treat outcome discontinuity.}
#'   \item{\code{itt.0}, \code{itt.1}}{For fuzzy RD only, outcome summaries
#'   for the control and treated sides requested through
#'   \code{params.other}; otherwise \code{NA}.}
#'   \item{\code{fs}}{For fuzzy RD only, same structure as \code{main} but for
#'   the first-stage treatment receipt/status discontinuity.}
#'   \item{\code{fs.0}, \code{fs.1}}{For fuzzy RD only, first-stage summaries
#'   for the control and treated sides requested through
#'   \code{params.other}; otherwise \code{NA}.}
#'   \item{\code{tau.hat}}{Point estimates at each evaluation point.}
#'   \item{\code{tau.hat.q}}{Bias-corrected estimates at each evaluation
#'   point.}
#'   \item{\code{se.hat}}{Standard errors corresponding to \code{tau.hat}.}
#'   \item{\code{se.hat.q}}{Standard errors corresponding to
#'   \code{tau.hat.q}.}
#'   \item{\code{params.cov}}{List of covariance matrices for aggregate
#'   inference in \code{\link{summary.rd2d}}. Contains only entries requested
#'   through \code{params.cov}.}
#'   \item{\code{cb}}{List with pointwise confidence intervals.}
#'   \item{\code{pvalues}}{Two-sided p-values based on bias-corrected
#'   estimates.}
#'   \item{\code{tvalues}}{t-statistics based on bias-corrected estimates.}
#'   \item{\code{tau.itt}, \code{tau.itt.q}}{For fuzzy RD, reduced-form outcome
#'   estimates using polynomial orders \eqn{p} and \eqn{q}; otherwise
#'   \code{NA}.}
#'   \item{\code{tau.fs}, \code{tau.fs.q}}{For fuzzy RD, first-stage treatment receipt/status
#'   estimates using polynomial orders \eqn{p} and \eqn{q}; otherwise
#'   \code{NA}.}
#'   \item{\code{rdmodel}}{Character label describing the fitted RD model.}
#'   \item{\code{call}}{Matched function call.}
#'   \item{\code{opt}}{List of options used in the function call.}
#' }
#'
#' @seealso \code{\link{rdbw2d}}, \code{\link{print.rd2d}},
#' \code{\link{summary.rd2d}}
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
#' \item{\href{https://arxiv.org/abs/2505.05670}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods.}
#' \item{\href{https://arxiv.org/abs/2511.06474}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Boundary Discontinuity Designs: Theory and Practice.}
#' \item{\href{https://arxiv.org/abs/2505.07989}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2025).}
#' rd2d: Causal Inference in Boundary Discontinuity Designs.}
#' }
#'
#' @examples
#' # Simulated example
#' set.seed(123)
#' n <- 800
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' assignment <- as.numeric(X1 > 0)
#' Y <- 3 + 2 * X1 + 1.5 * X2 + assignment + rnorm(n)
#' X <- cbind(X1, X2)
#' b <- matrix(c(0, 0, 0, 1), ncol = 2)
#'
#' # Estimate treatment effect using rd2d
#' result <- rd2d(Y, X, assignment, b, params.cov = "main",
#'                masspoints = "off", bwcheck = 10, repp = 49)
#' print(result)
#' summary(result, cbands = "main")
#'
#' # Fuzzy RD example
#' fuzzy <- as.numeric(runif(n) < ifelse(assignment == 1, 0.8, 0.2))
#' Y.fuzzy <- 3 + 2 * X1 + 1.5 * X2 + 1.5 * fuzzy + rnorm(n)
#' result.fuzzy <- rd2d(Y.fuzzy, X, assignment, b, fuzzy = fuzzy,
#'                      bwparam = "main", masspoints = "off",
#'                      bwcheck = 10, repp = 49)
#' print(result.fuzzy)
#' summary(result.fuzzy, output = "itt")
#' @export
rd2d <- function(Y, X, assignment, b, h = NULL, deriv = c(0,0), tangvec = NULL,
                 p = 1, q = NULL,
                 kernel = c("tri", "triangular", "epa", "epanechnikov",
                            "uni", "uniform", "gau", "gaussian"),
                 kernel_type = c("prod","rad"),
                 vce = c("hc1","hc0","hc2","hc3"),
                 masspoints = c("check", "adjust", "off"),cluster = NULL,
                 level = 95, params.other = NULL, params.cov = NULL,
                 side = c("two", "left", "right"), repp = 1000,
                 bwselect = c("mserd", "cerrd", "imserd", "icerrd",
                              "msetwo", "certwo", "imsetwo", "icertwo",
                              "user provided"),
                 bwparam = c("main", "itt"),
                 method = c("dpi", "rot"), bwcheck = 50 + p + 1,
                 scaleregul = 3, scalebiascrct = 1, stdvars = TRUE,
                 fuzzy = NULL){

  ######################## Input error handling ################################

  kernel <- match.arg(kernel)
  kernel_type <- match.arg(kernel_type)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  side <- match.arg(side)
  bwselect <- match.arg(bwselect)
  bwparam <- match.arg(bwparam)
  method <- match.arg(method)

  d <- assignment # renaming the variable
  is.fuzzy <- !is.null(fuzzy)
  if (!is.fuzzy) bwparam <- "main"

  exit <- 0

  # X must be matrix/data.frame with 2 columns before checking nrow/ncol.
  valid.X <- is.matrix(X) || is.data.frame(X)
  if (!valid.X || ncol(X) != 2) {
    print("X must be a matrix or data frame with exactly 2 columns")
    exit <- 1
  } else if (length(Y) != length(d) || length(Y) != nrow(X)) {
    print("Y, assignment, and rows of X must have the same length")
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

  if (!is.null(cluster) && length(cluster) != length(Y)) {
    print("cluster must have the same length as Y")
    exit <- 1
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
    print(
      paste(
        "params.other must contain only:",
        paste(allowed.params.other, collapse = ", ")
      )
    )
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
    print(
      paste(
        "params.cov must contain only:",
        paste(allowed.params.cov, collapse = ", ")
      )
    )
    exit <- 1
  } else {
    params.cov <- unique(params.cov)
  }

  if (exit == 0) {
    if (is.fuzzy) {
      invalid.design.params <- intersect(params.other, c("main.0", "main.1"))
      if (length(invalid.design.params) > 0) {
        print("main.0 and main.1 are available only for sharp rd2d fits")
        exit <- 1
      }
    } else {
      invalid.design.params <- setdiff(params.other, c("main.0", "main.1"))
      if (length(invalid.design.params) > 0) {
        print("itt.0, itt.1, fs.0, and fs.1 are available only for fuzzy rd2d fits")
        exit <- 1
      }
    }
  }

  if (exit == 0 && length(params.cov) > 0) {
    if (is.fuzzy) {
      invalid.cov <- intersect(params.cov, c("main.0", "main.1"))
      if (length(invalid.cov) > 0) {
        print("main.0 and main.1 covariance matrices are available only for sharp rd2d fits")
        exit <- 1
      }
      missing.side.cov <- setdiff(
        intersect(params.cov, c("itt.0", "itt.1", "fs.0", "fs.1")),
        params.other
      )
      if (length(missing.side.cov) > 0) {
        print(
          paste(
            "params.cov side outputs must also be requested in params.other:",
            paste(missing.side.cov, collapse = ", ")
          )
        )
        exit <- 1
      }
    } else {
      invalid.cov <- setdiff(params.cov, c("main", "main.0", "main.1"))
      if (length(invalid.cov) > 0) {
        print("itt, itt.0, itt.1, fs, fs.0, and fs.1 covariance matrices are available only for fuzzy rd2d fits")
        exit <- 1
      }
      missing.side.cov <- setdiff(
        intersect(params.cov, c("main.0", "main.1")),
        params.other
      )
      if (length(missing.side.cov) > 0) {
        print(
          paste(
            "params.cov side outputs must also be requested in params.other:",
            paste(missing.side.cov, collapse = ", ")
          )
        )
        exit <- 1
      }
    }
  }

  # d must be logical or contain only 0 and 1
  if (!(is.logical(d) || all(d %in% c(0, 1)))) {
    print(
      "d must be a logical vector or a numeric vector containing only 0 and 1"
    )
    exit <- 1
  }

  # b must be a matrix with 2 columns
  if (!(is.matrix(b) || is.data.frame(b)) || ncol(b) != 2) {
    print("b must be a matrix with 2 columns")
    exit <- 1
  }

  valid.p <- is.numeric(p) && length(p) == 1 &&
    is.finite(p) && p >= 0 &&
    abs(p - round(p)) < sqrt(.Machine$double.eps)
  if (!valid.p) {
    print("p must be a nonnegative integer")
    exit <- 1
  } else {
    p <- as.integer(round(p))
  }

  if (is.null(q)) {
    q <- if (valid.p) p + 1 else NA_real_
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

  # h must be scalar or matrix/data.frame with same rows as b and 4 columns
  if (!is.null(h)) {
    if (length(h) == 1) {
      if (!is.numeric(h) || !is.finite(h) || h <= 0) {
        print("If h is a scalar, it must be a positive numeric value")
        exit <- 1
      }
    } else if (!(is.matrix(h) || is.data.frame(h)) ||
               nrow(h) != nrow(b) || ncol(h) != 4) {
      print(
        paste(
          "If h is not a scalar, it must be a matrix or data frame with",
          "the same number of rows as b and 4 columns"
        )
      )
      exit <- 1
    } else {
      h.numeric <- suppressWarnings(as.matrix(h) + 0)
      if (!is.numeric(h.numeric) || any(!is.finite(h.numeric)) ||
          any(h.numeric <= 0)) {
        print(
          paste(
            "If h is a matrix or data frame, all entries must be finite",
            "positive numeric values"
          )
        )
        exit <- 1
      }
    }
  }

  # deriv must be a valid two-dimensional multi-index with sum <= p
  valid.deriv <- is.numeric(deriv) && length(deriv) == 2 &&
    all(is.finite(deriv)) && all(deriv >= 0) &&
    all(abs(deriv - round(deriv)) < sqrt(.Machine$double.eps))
  if (!valid.deriv) {
    print("deriv must be a nonnegative integer vector of length 2")
    exit <- 1
  } else {
    deriv <- as.integer(round(deriv))
  }
  if (valid.deriv && valid.p && sum(deriv) > p) {
    print(
      "Sum of deriv components must be less than or equal to polynomial order p"
    )
    exit <- 1
  }

  # tangvec, if provided, must be matrix/data.frame matching b rows and cols
  if (!is.null(tangvec)) {
    if (!(is.matrix(tangvec) || is.data.frame(tangvec)) ||
        nrow(tangvec) != nrow(b) || ncol(tangvec) != 2) {
      warning(
        paste(
          "tangvec must be a matrix or data frame with same number of rows",
          "as b and 2 columns"
        )
      )
      exit <- 1
    }
    if (valid.p && p < 1) {
      print("tangvec requires polynomial order p >= 1")
      exit <- 1
    }
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

  if (valid.p && valid.q && q < p){
    print("Parameter q must be no smaller than p. Please provide valid inputs.")
    exit <- 1
  }

  if (exit == 0 && !is.null(bwcheck)) {
    if (!is.numeric(bwcheck) || length(bwcheck) != 1 ||
        !is.finite(bwcheck) || bwcheck < 1 || bwcheck != as.integer(bwcheck)) {
      print("bwcheck must be NULL or a positive integer")
      exit <- 1
    } else {
      bwcheck <- as.integer(bwcheck)
    }
  }

  if (!is.null(cluster) && !(vce %in% c("hc0", "hc1"))) {
    warning(
      paste(
        "When cluster is specified, vce must be 'hc0' or 'hc1'.",
        "Resetting vce to 'hc1'."
      )
    )
    vce <- "hc1"
  }

  if (exit>0) stop()

  ############################ Data preparation ################################

  dat <- cbind(X[,1], X[,2], Y, d)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("x.1", "x.2", "y", "d")
  if (is.fuzzy) dat$fuzzy <- as.numeric(fuzzy)
  eval <- as.data.frame(b)
  colnames(eval) <- c("x.1", "x.2")
  neval <- dim(eval)[1]
  na.ok <- complete.cases(dat$x.1) &
    complete.cases(dat$x.2) &
    complete.cases(dat$y) &
    complete.cases(dat$d)
  if (is.fuzzy) na.ok <- na.ok & complete.cases(dat$fuzzy)
  if (!is.null(cluster)) na.ok <- na.ok & complete.cases(cluster)
  dat <- dat[na.ok,]
  if (!is.null(cluster)) cluster <- cluster[na.ok]
  N <- dim(dat)[1]
  N.0 <- dim(dat[dat$d == 0,])[1]
  N.1 <- dim(dat[dat$d == 1,])[1]

  if (is.null(p))         p <- 1
  kernel   <- tolower(kernel)

  e_deriv <- matrix(
    0,
    nrow = neval,
    ncol = factorial(p+2)/(factorial(p) * factorial(2))
  )
  deriv.sum <- deriv[1] + deriv[2]
  if (deriv.sum >= 1){
    deriv.index <- (factorial(deriv.sum+1) /
      (factorial(deriv.sum-1)* 2)) + deriv[2] + 1
    e_deriv[, deriv.index] <- prod(factorial(deriv))
  } else {
    e_deriv[,1] <- 1
  }

  if (!is.null(tangvec)){
    warning("Tangvec provided. Ignore option deriv.")
    e_deriv <- matrix(
      0,
      nrow = neval,
      ncol = factorial(p+2)/(factorial(p) * factorial(2))
    )
    e_deriv[,2] <- tangvec[,1]
    e_deriv[,3] <- tangvec[,2]
    deriv <- c(1,0) # standardization for latter codes
    deriv.sum <- 1
  }

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Check for mass points

  M <- NULL; M.0 <- NULL; M.1 <- NULL
  unique <- NULL


  # Only check for mass points if inputting bivariate coordinates.
  if (masspoints == "check" | masspoints == "adjust"){
    unique.const <- rd2d_unique(dat)
    unique <- unique.const$unique
    M.0 <- dim(unique[unique$d == 0,])[1]
    M.1 <- dim(unique[unique$d == 1,])[1]
    M <- M.0 + M.1
    mass <- 1 - M / N
    if (mass >= 0.2){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) &
          (masspoints == "check" | masspoints == "adjust")) {
        bwcheck <- 50 + p + 1
      }
    }
  }

  if (!is.null(bwcheck)) {
    if (masspoints == "adjust") {
      available.0 <- M.0
      available.1 <- M.1
      unit.label <- "unique mass points"
    } else {
      available.0 <- N.0
      available.1 <- N.1
      unit.label <- "observations"
    }

    if (available.0 < bwcheck || available.1 < bwcheck) {
      stop(
        sprintf(
          paste(
            "bwcheck=%d requires at least %d %s on each side of the cutoff;",
            "found %d on the control side and %d on the treatment side.",
            "Decrease bwcheck or set bwcheck=NULL."
          ),
          bwcheck, bwcheck, unit.label, available.0, available.1
        ),
        call. = FALSE
      )
    }
  }

  # if (!is.null(cluster)){
  #   unique.0 <- unique(cluster[d == FALSE])
  #   unique.1 <- unique(cluster[d == TRUE])
  #   M.0 <- length(unique.0)
  #   M.1 <- length(unique.1)
  #   M <- M.0 + M.1
  # }

  ################################ Bandwidth ###################################

  if (is.null(h)){
    # if (bwselect == "mserd" | bwselect == "imserd"){
    #   bws <- (rdbw2d(
    #     dat, eval, o, p, deriv, kernel, bwselect, method, vce, bwcheck,
    #     masspoints, cluster, scaleregul, bydist, inputdist
    #   ))$bws
    #   hgrid <- bws$hMSE.0
    #   hgrid.1 <- NULL
    # } else if (bwselect == "msetwo" | bwselect == "imsetwo"){
    #   bws <- (rdbw2d(
    #     dat, eval, o, p, deriv, kernel, bwselect, method, vce, bwcheck,
    #     masspoints, cluster, scaleregul, bydist, inputdist
    #   ))$bws
    #   hgrid <- bws$hMSE.0
    #   hgrid.1 <- bws$hMSE.1
    # }
    bws <- rdbw2d(
      Y = dat$y, X = as.matrix(dat[, c("x.1", "x.2")]), assignment = dat$d,
      b = b, p = p,
      deriv = deriv, tangvec = tangvec, kernel = kernel,
      kernel_type = kernel_type, bwselect = bwselect, bwparam = bwparam,
      method = method,
      vce = vce, bwcheck = bwcheck, masspoints = masspoints,
      cluster = cluster, scaleregul = scaleregul, scalebiascrct = scalebiascrct,
      stdvars = stdvars, fuzzy = if (is.fuzzy) dat$fuzzy else NULL
    )
    bws <- bws$bws
    hgrid <- cbind(bws[,3],bws[,4])
    hgrid.1 <- cbind(bws[,5],bws[,6])
  } else {
    bwselect <- "user provided"
    # standardize bandwidth
    if (length(h) == 1){
      hgrid <- matrix(h, nrow = neval, ncol = 2)
      hgrid.1 <- matrix(h, nrow = neval, ncol = 2)
    } else {
      hgrid <- cbind(h[,1],h[,2])
      hgrid.1 <- cbind(h[,3],h[,4])
    }
  }

  bwcheck.bounds <- rd2d_bwcheck_bounds(
    dat, eval, bwcheck, masspoints, unique = unique, metric = kernel_type,
    scale = c(sd(dat$x.1), sd(dat$x.2))
  )

  ###################### Point estimation and inference ########################

  count.q <- factorial(q + 2)/(factorial(q) * 2)
  count.p <- factorial(p + 2)/(factorial(p) * 2)

  e_deriv.q <- matrix(0, nrow = neval, ncol = count.q)
  e_deriv.q[,c(1:count.p)] <- e_deriv

  if (!is.fuzzy){
    rdfit.p <- rd2d_fit(
      dat, eval, e_deriv, deriv, p, hgrid, hgrid.1, kernel, kernel_type, vce,
      masspoints, cluster, bwcheck, unique, bwcheck.bounds
    )
    tau.hat.p <- rdfit.p$mu.1 - rdfit.p$mu.0
    se.hat.p <- sqrt(rdfit.p$se.0^2 + rdfit.p$se.1^2)
    h.01.p <- rdfit.p$h.0.x
    h.02.p <- rdfit.p$h.0.y
    h.11.p <- rdfit.p$h.1.x
    h.12.p <- rdfit.p$h.1.y
    eN.0.p <- rdfit.p$eN.0
    eN.1.p <- rdfit.p$eN.1

    if (q == p) {
      rdfit.q <- rdfit.p
    } else {
      rdfit.q <- rd2d_fit(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel, kernel_type,
        vce, masspoints, cluster, bwcheck, unique, bwcheck.bounds
      )
    }
    tau.hat.q <- rdfit.q$mu.1 - rdfit.q$mu.0
    se.hat.q <- sqrt(rdfit.q$se.0^2 + rdfit.q$se.1^2)
    h.01.q <- rdfit.q$h.0.x
    h.02.q <- rdfit.q$h.0.y
    h.11.q <- rdfit.q$h.1.x
    h.12.q <- rdfit.q$h.1.y
    eN.0.q <- rdfit.q$eN.0
    eN.1.q <- rdfit.q$eN.1
  } else {
    dat.fs <- dat
    dat.fs$y <- dat$fuzzy

    rdfit.p.multi <- rd2d_fit_multi(
      dat, eval, e_deriv, deriv, p, hgrid, hgrid.1, kernel, kernel_type, vce,
      masspoints, cluster, bwcheck, unique, bwcheck.bounds
    )
    rdfit.itt.p <- rdfit.p.multi$y
    rdfit.fs.p <- rdfit.p.multi$fuzzy
    tau.itt.p <- rdfit.itt.p$mu.1 - rdfit.itt.p$mu.0
    tau.fs.p <- rdfit.fs.p$mu.1 - rdfit.fs.p$mu.0
    se.itt.p <- sqrt(rdfit.itt.p$se.0^2 + rdfit.itt.p$se.1^2)
    se.fs.p <- sqrt(rdfit.fs.p$se.0^2 + rdfit.fs.p$se.1^2)

    if (q == p) {
      rdfit.q.multi <- rdfit.p.multi
    } else {
      rdfit.q.multi <- rd2d_fit_multi(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel, kernel_type,
        vce, masspoints, cluster, bwcheck, unique, bwcheck.bounds
      )
    }
    rdfit.itt.q <- rdfit.q.multi$y
    rdfit.fs.q <- rdfit.q.multi$fuzzy
    tau.itt.q <- rdfit.itt.q$mu.1 - rdfit.itt.q$mu.0
    tau.fs.q <- rdfit.fs.q$mu.1 - rdfit.fs.q$mu.0
    se.itt.q <- sqrt(rdfit.itt.q$se.0^2 + rdfit.itt.q$se.1^2)
    se.fs.q <- sqrt(rdfit.fs.q$se.0^2 + rdfit.fs.q$se.1^2)

    denom.tol <- sqrt(.Machine$double.eps)
    valid.p <- is.finite(tau.fs.p) & abs(tau.fs.p) > denom.tol
    valid.q <- is.finite(tau.fs.q) & abs(tau.fs.q) > denom.tol
    if (any(!valid.p | !valid.q)) {
      warning(
        paste(
          "Weak or zero first-stage fuzzy RD estimates detected;",
          "returning NA for affected fuzzy estimates."
        )
      )
    }

    tau.hat.p <- rep(NA_real_, neval)
    tau.hat.q <- rep(NA_real_, neval)
    tau.hat.p[valid.p] <- tau.itt.p[valid.p] / tau.fs.p[valid.p]
    tau.hat.q[valid.q] <- tau.itt.q[valid.q] / tau.fs.q[valid.q]

    h.01.p <- rdfit.itt.p$h.0.x
    h.02.p <- rdfit.itt.p$h.0.y
    h.11.p <- rdfit.itt.p$h.1.x
    h.12.p <- rdfit.itt.p$h.1.y
    eN.0.p <- rdfit.itt.p$eN.0
    eN.1.p <- rdfit.itt.p$eN.1

    h.01.q <- rdfit.itt.q$h.0.x
    h.02.q <- rdfit.itt.q$h.0.y
    h.11.q <- rdfit.itt.q$h.1.x
    h.12.q <- rdfit.itt.q$h.1.y
    eN.0.q <- rdfit.itt.q$eN.0
    eN.1.q <- rdfit.itt.q$eN.1

    se.hat.p <- rdbw2d_se_fuzzy(
      dat, eval, e_deriv, deriv, p, hgrid, hgrid.1,
      kernel, kernel_type, vce, cluster, tau.itt.p, tau.fs.p, se.itt.p, se.fs.p,
      denom.tol = denom.tol
    )
    if (q == p) {
      se.hat.q <- se.hat.p
    } else {
      se.hat.q <- rdbw2d_se_fuzzy(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1,
        kernel, kernel_type, vce, cluster, tau.itt.q, tau.fs.q, se.itt.q, se.fs.q,
        denom.tol = denom.tol
      )
    }
  }

  tvalues <- tau.hat.q/se.hat.q
  pvalues <- 2 * pnorm(abs(tvalues),lower.tail = FALSE)

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

  # Covariance matrices are optional and stored only for later summary() use.

  cov.tables <- list()
  if (!is.fuzzy) {
    if ("main" %in% params.cov) {
      cov.tables$main <- rdbw2d_cov(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster
      )
    }
    if ("main.0" %in% params.cov) {
      cov.tables$main.0 <- rdbw2d_cov_side(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 0
      )
    }
    if ("main.1" %in% params.cov) {
      cov.tables$main.1 <- rdbw2d_cov_side(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 1
      )
    }
  } else if (length(params.cov) > 0) {
    cov.tables <- rdbw2d_cov_fuzzy_tables(
      dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1,
      kernel, kernel_type, vce, cluster, tau.itt.q, tau.fs.q,
      outputs = params.cov, denom.tol = denom.tol
    )
  } else {
    if ("main" %in% params.cov) {
      cov.tables$main <- rdbw2d_cov_fuzzy(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1,
        kernel, kernel_type, vce, cluster, tau.itt.q, tau.fs.q,
        denom.tol = denom.tol
      )
    }
    if ("itt" %in% params.cov) {
      cov.tables$itt <- rdbw2d_cov(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster
      )
    }
    if ("fs" %in% params.cov) {
      cov.tables$fs <- rdbw2d_cov(
        dat.fs, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster
      )
    }
    if ("itt.0" %in% params.cov) {
      cov.tables$itt.0 <- rdbw2d_cov_side(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 0
      )
    }
    if ("itt.1" %in% params.cov) {
      cov.tables$itt.1 <- rdbw2d_cov_side(
        dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 1
      )
    }
    if ("fs.0" %in% params.cov) {
      cov.tables$fs.0 <- rdbw2d_cov_side(
        dat.fs, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 0
      )
    }
    if ("fs.1" %in% params.cov) {
      cov.tables$fs.1 <- rdbw2d_cov_side(
        dat.fs, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel,
        kernel_type, vce, cluster, side = 1
      )
    }
  }

  cb.hat.q <- list(CI.l = CI.lower, CI.r = CI.upper)

  clustered <- !is.null(cluster)

  ################################## Output ####################################

  main <- cbind(
    b[,1], b[,2], tau.hat.p, se.hat.p, tau.hat.q, se.hat.q, tvalues,
    pvalues, CI.lower, CI.upper, hgrid[,1], hgrid[,2],
    hgrid.1[,1], hgrid.1[,2], eN.0.p, eN.1.p
  )
  main <- as.data.frame(main)
  colnames(main) <- c(
    "b1","b2","estimate.p","std.err.p","estimate.q","std.err.q",
    "t.value", "p.value",
    "ci.lower","ci.upper", "h01", "h02",
    "h11", "h12", "N.Co", "N.Tr"
  )
  main.names <- colnames(main)

  make_result_table <- function(est.p, se.p, est.q, se.q,
                                h01 = hgrid[, 1], h02 = hgrid[, 2],
                                h11 = hgrid.1[, 1], h12 = hgrid.1[, 2],
                                NCo = eN.0.p, NTr = eN.1.p) {
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
      b1 = b[, 1],
      b2 = b[, 2],
      estimate.p = est.p,
      std.err.p = se.p,
      estimate.q = est.q,
      std.err.q = se.q,
      t.value = t.comp,
      p.value = p.comp,
      ci.lower = ci.lower.comp,
      ci.upper = ci.upper.comp,
      h01 = h01,
      h02 = h02,
      h11 = h11,
      h12 = h12,
      N.Co = NCo,
      N.Tr = NTr,
      check.names = FALSE
    )
    out.comp[, main.names]
  }

  make_bw_table <- function() {
    data.frame(
      b1 = b[, 1],
      b2 = b[, 2],
      h01 = hgrid[, 1],
      h02 = hgrid[, 2],
      h11 = hgrid.1[, 1],
      h12 = hgrid.1[, 2],
      N.Co = eN.0.p,
      N.Tr = eN.1.p,
      check.names = FALSE
    )
  }

  bw <- make_bw_table()

  if (!is.fuzzy) {
    main.0 <- NA
    main.1 <- NA

    if ("main.0" %in% params.other) {
      main.0 <- make_result_table(
        rdfit.p$mu.0, rdfit.p$se.0, rdfit.q$mu.0, rdfit.q$se.0,
        h11 = rep(NA_real_, neval), h12 = rep(NA_real_, neval),
        NTr = rep(NA_real_, neval)
      )
    }

    if ("main.1" %in% params.other) {
      main.1 <- make_result_table(
        rdfit.p$mu.1, rdfit.p$se.1, rdfit.q$mu.1, rdfit.q$se.1,
        h01 = rep(NA_real_, neval), h02 = rep(NA_real_, neval),
        NCo = rep(NA_real_, neval)
      )
    }

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
        rdfit.itt.p$mu.0, rdfit.itt.p$se.0, rdfit.itt.q$mu.0, rdfit.itt.q$se.0,
        h11 = rep(NA_real_, neval), h12 = rep(NA_real_, neval),
        NTr = rep(NA_real_, neval)
      )
    }

    if ("itt.1" %in% params.other) {
      itt.1 <- make_result_table(
        rdfit.itt.p$mu.1, rdfit.itt.p$se.1, rdfit.itt.q$mu.1, rdfit.itt.q$se.1,
        h01 = rep(NA_real_, neval), h02 = rep(NA_real_, neval),
        NCo = rep(NA_real_, neval)
      )
    }

    if ("fs.0" %in% params.other) {
      fs.0 <- make_result_table(
        rdfit.fs.p$mu.0, rdfit.fs.p$se.0, rdfit.fs.q$mu.0, rdfit.fs.q$se.0,
        h11 = rep(NA_real_, neval), h12 = rep(NA_real_, neval),
        NTr = rep(NA_real_, neval)
      )
    }

    if ("fs.1" %in% params.other) {
      fs.1 <- make_result_table(
        rdfit.fs.p$mu.1, rdfit.fs.p$se.1, rdfit.fs.q$mu.1, rdfit.fs.q$se.1,
        h01 = rep(NA_real_, neval), h02 = rep(NA_real_, neval),
        NCo = rep(NA_real_, neval)
      )
    }
  }

  rdmodel <- ifelse(is.fuzzy, "fuzzy rd2d", "rd2d")

  result.tables <- list(
    main = main,
    bw = bw
  )
  if (!is.fuzzy) {
    result.tables <- c(
      result.tables,
      list(
        main.0 = main.0,
        main.1 = main.1
      )
    )
  } else {
    result.tables <- c(
      result.tables,
      list(
        itt = itt,
        itt.0 = itt.0,
        itt.1 = itt.1,
        fs = fs,
        fs.0 = fs.0,
        fs.1 = fs.1
      )
    )
  }

  out <- c(
    result.tables,
    list(
      opt = list(
        b = b, deriv = deriv, tangvec = tangvec, p = p, q = q,
        kernel = kernel.type, kernel_type = kernel_type, N = N, N.0 = N.0,
        N.1 = N.1, M = M, M.0 = M.0, M.1 = M.1, neval = neval,
        bwselect = bwselect, bwparam = bwparam, method = method,
        vce = vce, bwcheck = bwcheck,
        masspoints = masspoints, cluster = cluster, clustered = clustered,
        scaleregul = scaleregul, scalebiascrct = scalebiascrct,
        stdvars = stdvars, fuzzy = is.fuzzy, level = level, repp = repp,
        side = side, params.other = params.other, params.cov = params.cov,
        h01 = hgrid[,1], h02 = hgrid[,2],
        h11 = hgrid.1[,1], h12 = hgrid.1[,2],
        N.Co = eN.0.p, N.Tr = eN.1.p
      ),
      tau.hat = tau.hat.p,
      tau.hat.q = tau.hat.q,
      se.hat = se.hat.p,
      se.hat.q = se.hat.q,
      params.cov = cov.tables,
      cb = cb.hat.q,
      pvalues = pvalues,
      tvalues = tvalues,
      tau.itt = tau.itt.p,
      tau.itt.q = tau.itt.q,
      tau.fs = tau.fs.p,
      tau.fs.q = tau.fs.q,
      rdmodel = rdmodel
    )
  )
  out$call   <- match.call()
  class(out) <- "rd2d"

  return(out)
}

################################################################################
#' Print Method for 2D Local Polynomial RD Estimation
#' @description
#' Prints the results of a 2D local polynomial regression discontinuity (RD)
#' estimation, as obtained from \code{\link{rd2d}}.
#'
#' @param x An object of class \code{rd2d}, returned by \code{\link{rd2d}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return
#' No return value. This function is called for its side effects, which are to
#' print the \code{\link{rd2d}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso
#' \code{\link{rd2d}} for conducting 2D local polynomial RD estimation.
#'
#' Supported methods: \code{\link{print.rd2d}}, \code{\link{summary.rd2d}}.
#'
#' @export

print.rd2d <- function(x,...) {

  cat(paste(x$rdmodel, "\n", sep = ""))
  cat(paste("\n", sep = ""))

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  bw.type <- paste(
    paste(x$opt$bwselect, x$opt$method, sep = "-"),
    ifelse(x$opt$stdvars, "-std", ""),
    sep = ""
  )
  kernel.name <- paste(
    tolower(x$opt$kernel),
    x$opt$kernel_type,
    sep = "-"
  )
  vce.method <- paste(
    x$opt$vce,
    ifelse(x$opt$clustered, "-clustered", ""),
    sep = ""
  )
  cat(sprintf("BW type.               %s\n", bw.type))
  cat(sprintf("Kernel                 %s\n", kernel.name))
  cat(sprintf("VCE method             %s\n", vce.method))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  if (isTRUE(x$opt$fuzzy)) {
    cat(sprintf("Fuzzy                  on\n"))
    cat(sprintf("BW target              %s\n", x$opt$bwparam))
  }
  # cat(sprintf(
  #   "Standardization        %s\n",
  #   ifelse(x$opt$stdvars, "on", "off")
  # ))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf(
    "Estimand (deriv)       %-10d   %-10d\n",
    x$opt$deriv[1], x$opt$deriv[2]
  ))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")

}


################################################################################
#' Summary Method for 2D Local Polynomial RD Estimation
#'
#' @description
#' Summarizes estimation and bandwidth results from a 2D local polynomial
#' regression discontinuity (RD) design, as produced by \code{\link{rd2d}}.
#'
#' @param object An object of class \code{rd2d}, typically returned by
#' \code{\link{rd2d}}.
#' @param ... Optional named arguments. Unsupported option names produce an
#'   error. Supported options include:
#'   \itemize{
#'     \item \code{cbands}: Optional character vector requesting uniform
#'       confidence bands for displayed outputs. Options are \code{"main"},
#'       \code{"main.0"}, \code{"main.1"}, \code{"itt"}, \code{"itt.0"},
#'       \code{"itt.1"}, \code{"fs"}, \code{"fs.0"}, and \code{"fs.1"},
#'       subject to the sharp/fuzzy design and covariance matrices stored by
#'       \code{\link{rd2d}}.
#'     \item \code{WBATE}: Optional numeric weights for a weighted boundary
#'       average treatment effect row. The weights must match the full set of
#'       evaluation points for the selected output. The selected output must
#'       have a matching covariance matrix stored by \code{\link{rd2d}}.
#'     \item \code{LBATE}: Logical. If \code{TRUE}, prints a largest boundary
#'       average treatment effect row. The selected output must have a
#'       matching covariance matrix stored by \code{\link{rd2d}}.
#'     \item \code{subset}: Integer vector of indices of evaluation points to
#'       display. Aggregates and uniform-band critical values are still
#'       computed using all evaluation points. Defaults to all evaluation
#'       points.
#'     \item \code{output}: Character. Use \code{"main"} to display sharp or
#'       fuzzy main results; \code{"main.0"} or \code{"main.1"} to display
#'       sharp control- or treatment-side main results; \code{"itt"},
#'       \code{"itt.0"}, or \code{"itt.1"} to display fuzzy
#'       reduced-form or intention-to-treat results; \code{"fs"},
#'       \code{"fs.0"}, or \code{"fs.1"} to display fuzzy first-stage treatment receipt/status
#'       results; and
#'       \code{"bw"} to display bandwidth information. If omitted,
#'       \code{summary()} prints all default and requested estimate tables
#'       and does not print \code{"bw"}. Use \code{output = "bw"} to print
#'       bandwidth information.
#'     \item \code{sep}: Integer vector of length three. Controls spacing in
#'       the output.
#'       \code{sep[1]} controls spacing for the columns of bandwidths,
#'       estimation,
#'       t-statistic, and p-value in the \code{"main"} table.
#'       \code{sep[2]} controls spacing for the confidence interval
#'       (confidence bands)
#'       in the \code{"main"} table.
#'       \code{sep[3]} controls spacing for the columns in the
#'       \code{"bw"} table.
#'       Default is \code{c(7, 17, 8)}.
#'   }
#'
#' @return Invisibly returns a list with displayed tables and uniform
#' confidence bands requested through \code{cbands}. Each returned estimate
#' table has the same columns as the corresponding \code{\link{rd2d}} output,
#' with \code{cb.lower} and \code{cb.upper} added only when uniform bands are
#' requested. Requested WBATE and LBATE rows are appended to the corresponding
#' returned table.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com} \cr
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}
#'
#' @seealso \code{\link{rd2d}} for estimation using 2D local polynomial RD
#' design.
#'
#' Supported methods: \code{\link{print.rd2d}}, \code{\link{summary.rd2d}}.
#'
#' @export

summary.rd2d <- function(object, ...) {

  x <- object

  args <- list(...)
  valid.summary.args <- c("cbands", "WBATE", "LBATE", "subset", "output", "sep")
  arg.names <- names(args)
  if (length(args) > 0 && (is.null(arg.names) || any(arg.names == ""))) {
    stop("All summary.rd2d options must be named.", call. = FALSE)
  }
  invalid.args <- setdiff(arg.names, valid.summary.args)
  if (length(invalid.args) > 0) {
    stop(
      sprintf(
        "Unsupported summary.rd2d option(s): %s.",
        paste(invalid.args, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (is.null(args[['cbands']])) {
    cbands <- character(0)
  } else {
    cbands <- args[['cbands']]
    if (!is.character(cbands) || any(is.na(cbands))) {
      stop("cbands must be NULL or a character vector.", call. = FALSE)
    }
    cbands <- unique(cbands)
  }
  WBATE <- args[['WBATE']]
  LBATE <- isTRUE(args[['LBATE']])

  if (is.null(args[['subset']])) {
    subset <- NULL
  } else {
    subset <- args[['subset']]
  }

  output.requested <- !is.null(args[['output']])
  if (!output.requested) {
    output <- "main"
  } else {
    output <- args[['output']]
  }

  if (is.null(args[['sep']])) {
    sep <- c(7,17,8)
  } else {
    sep <- args[['sep']]
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

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  bw.type <- paste(
    paste(x$opt$bwselect, x$opt$method, sep = "-"),
    ifelse(x$opt$stdvars, "-std", ""),
    sep = ""
  )
  kernel.name <- paste(
    tolower(x$opt$kernel),
    x$opt$kernel_type,
    sep = "-"
  )
  vce.method <- paste(
    x$opt$vce,
    ifelse(x$opt$clustered, "-clustered", ""),
    sep = ""
  )
  cat(sprintf("BW type.               %s\n", bw.type))
  cat(sprintf("Kernel                 %s\n", kernel.name))
  cat(sprintf("VCE method             %s\n", vce.method))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  if (isTRUE(x$opt$fuzzy)) {
    cat(sprintf("Fuzzy                  on\n"))
    cat(sprintf("BW target              %s\n", x$opt$bwparam))
  }
  # cat(sprintf(
  #   "Standardization        %s\n",
  #   ifelse(x$opt$stdvars, "on", "off")
  # ))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf(
    "Estimand (deriv)       %-10d   %-10d\n",
    x$opt$deriv[1], x$opt$deriv[2]
  ))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")

  print_table <- function(output) {

  valid.outputs <- valid.outputs.all
  if (!(output %in% valid.outputs)) {
    warning(
      paste(
        sprintf("output='%s' is not available for this rd2d object.", output),
        "Resetting to output='main'."
      )
    )
    output <- "main"
  }

  if (!is.data.frame(x[[output]])) {
    warning(
      paste(
        sprintf("output='%s' was not computed for this rd2d object.", output),
        "Request companion outputs with params.other in rd2d().",
        "Resetting to output='main'."
      )
    )
    output <- "main"
  }

  if (output == "bw" && (!is.null(WBATE) || LBATE)) {
    warning("WBATE and LBATE are not available for output='bw'.")
    WBATE <- NULL
    LBATE <- FALSE
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

  if (output != "bw"){
    bands.requested <- output %in% cbands
    headers <- c(
      "ID", "b1", "b2", "estimate", "t.value", "p.value",
      sprintf("%d%% CI", x$opt$level)
    )
    if (bands.requested){
      headers[length(headers)] <- sprintf("%d%% Unif. CB", x$opt$level)
    }

    col_widths <- pmax(c(4, sep[1], sep[1], sep[1], sep[1], sep[1], sep[2]),
                       nchar(headers))

    # Format and print header row
    rule.width <- sum(col_widths) + 2 * (length(headers) - 1)
    rule <- paste(rep("=", rule.width), collapse = "")
    cat(rule, "\n")
    formatted_headers <- mapply(
      function(h, w) formatC(h, width = w, format = "s"),
      headers,
      col_widths
    )
    cat(paste(formatted_headers, collapse = "  "), "\n")
    cat(rule, "\n")

    neval <- nrow(results)
    if (is.null(subset)){
      subset <- seq_len(neval)
    } else{
      # input error handling
      if (!all(subset %in% seq_len(neval))) {
        warning("Invalid subset provided. Resetting to default: 1:neval")
        subset <- seq_len(neval)
      }
    }
    result.subset <- results[subset,,drop = FALSE]

    cov.table <- NULL
    if (!is.null(x$params.cov) && is.matrix(x$params.cov[[output]])) {
      cov.table <- x$params.cov[[output]]
    }
    if (!is.null(cov.table) &&
        !identical(dim(cov.table), c(nrow(results), nrow(results)))) {
      cov.table <- NULL
    }
    summary.side <- x$opt$side

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

    covariance_available <- function(cov.mat, n) {
      !is.null(cov.mat) &&
        identical(dim(cov.mat), c(n, n)) &&
        all(is.finite(cov.mat)) &&
        all(is.finite(diag(cov.mat))) &&
        all(diag(cov.mat) > 0)
    }

    require_covariance <- function(reason) {
      if (!covariance_available(cov.table, nrow(results))) {
        stop(
          sprintf(
            paste(
              "%s requires a stored covariance matrix for output='%s'.",
              "Rerun rd2d(..., params.cov = \"%s\")."
            ),
            reason, output, output
          ),
          call. = FALSE
        )
      }
    }

    if (bands.requested) require_covariance("Uniform confidence bands")
    if (!is.null(WBATE)) require_covariance("WBATE inference")
    if (LBATE) require_covariance("LBATE inference")

    interval.lower <- results$ci.lower
    interval.upper <- results$ci.upper
    cbands.return <- NULL
    if (bands.requested) {
      cb.out <- rd2d_cb(
        results[["estimate.q"]], cov.table, x$opt$repp, summary.side,
        x$opt$level
      )
      interval.lower <- cb.out$CB.l
      interval.upper <- cb.out$CB.r
      result.subset$cb.lower <- cb.out$CB.l[subset]
      result.subset$cb.upper <- cb.out$CB.r[subset]
      cbands.return <- result.subset[, c("cb.lower", "cb.upper"), drop = FALSE]
    }

    print_aggregate_row <- function(label, est, t.value = NA_real_,
                                    p.value = NA_real_,
                                    lower = NA_real_, upper = NA_real_) {
      row_vals <- c(
        formatC(label, width = col_widths[1], format = "s"),
        format_blank(col_widths[2]),
        format_blank(col_widths[3]),
        format_number(est, col_widths[4]),
        format_number(t.value, col_widths[5]),
        format_number(p.value, col_widths[6]),
        format_interval(lower, upper, col_widths[7])
      )
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
      if (!LBATE) return(NULL)
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

    # Print each row of the results
    for (i in 1:nrow(results)) {

      index_formatted <- formatC(i, width = col_widths[1], format = "d")

      b1_formatted <- ifelse(
        is.na(results$b1[i]),
        formatC("NA", width = col_widths[2], format = "s"),
        formatC(results$b1[i], format = "f", digits = 3, width = col_widths[2])
      )

      b2_formatted <- ifelse(
        is.na(results$b2[i]),
        formatC("NA", width = col_widths[3], format = "s"),
        formatC(results$b2[i], format = "f", digits = 3, width = col_widths[3])
      )

      coef_formatted <- ifelse(
        is.na(results[["estimate.p"]][i]),
        formatC("NA", width = col_widths[4], format = "s"),
        formatC(
          results[["estimate.p"]][i],
          format = "f",
          digits = 4,
          width = col_widths[4]
        )
      )

      tvalues_formatted <- ifelse(
        is.na(results[["t.value"]][i]),
        formatC("NA", width = col_widths[5], format = "s"),
        formatC(
          results[["t.value"]][i],
          format = "f",
          digits = 4,
          width = col_widths[5]
        )
      )

      pvalues_formatted <- ifelse(
        is.na(results[i, "p.value"]),
        formatC("NA", width = col_widths[6], format = "s"),
        formatC(
          results[i, "p.value"],
          format = "f",
          digits = 4,
          width = col_widths[6]
        )
      )

      ci.lower <- interval.lower[i]
      ci.upper <- interval.upper[i]
      ci_formatted <- ifelse(
        is.na(ci.lower) | is.na(ci.upper),
        formatC("NA", width = col_widths[7], format = "s"),
        formatC(
          paste0(
            "[",
            formatC(ci.lower, format = "f", digits = 4),
            ", ",
            formatC(ci.upper, format = "f", digits = 4),
            "]"
          ),
          width = col_widths[7],
          format = "s"
        )
      )

      # Print
      row_vals <- c(
        index_formatted,
        b1_formatted,
        b2_formatted,
        coef_formatted,
        tvalues_formatted,
        pvalues_formatted,
        ci_formatted
      )
      if (i %in% subset) cat(paste(row_vals, collapse = "  "), "\n")
    }

    wbate.row <- make_wbate(WBATE, results, cov.table)
    lbate.row <- make_lbate(results, cov.table)
    aggregate.rows <- Filter(Negate(is.null), list(wbate.row, lbate.row))
    summary.table <- result.subset
    rownames(summary.table) <- as.character(subset)
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
      row$ci.lower <- aggregate.row$lower
      row$ci.upper <- aggregate.row$upper
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
    }
    if (length(aggregate.rows) > 0) {
      for (aggregate.row in aggregate.rows) {
        cat(paste(rep("-", rule.width), collapse = ""), "\n")
        print_aggregate_row(
          aggregate.row$label,
          aggregate.row$est,
          aggregate.row$t.value,
          aggregate.row$p.value,
          aggregate.row$lower,
          aggregate.row$upper
        )
      }
    }
    cat(rule, "\n")
    return(
      list(
        output = output,
        table = summary.table,
        cbands = cbands.return
      )
    )
  } else {
    headers <- c("ID", "b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr")
    col_widths <- c(
      4, sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3]
    )

    cat(strrep("=", sum(col_widths)), "\n")
    group_headers <- c(
      formatC(
        "        Bdy Points",
        width = col_widths[1] + col_widths[2] + col_widths[3],
        format = "s",
        flag = "-"
      ),
      formatC(
        "      BW Control",
        width = col_widths[4] + col_widths[5],
        format = "s",
        flag = "-"
      ),
      formatC(
        "    BW Treatment",
        width = col_widths[6] + col_widths[7],
        format = "s",
        flag = "-"
      ),
      formatC(
        "       Eff. N",
        width = col_widths[8] + col_widths[9],
        format = "s",
        flag = "-"
      )
    )
    cat(paste(group_headers, collapse = ""), "\n")

    # Format and print header row
    formatted_headers <- mapply(
      function(h, w) formatC(h, width = w, format = "s"),
      headers,
      col_widths
    )
    cat(paste(formatted_headers, collapse = ""), "\n")
    cat(strrep("=", sum(col_widths)), "\n")

    neval <- nrow(results)
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
    for (j in 1:nrow(results)) {
      index <- formatC(j, width = col_widths[1], format = "d")
      bdy1 <- ifelse(
        is.na(results[j, "b1"]),
        formatC("NA", width = col_widths[2], format = "s"),
        formatC(
          results[j, "b1"],
          format = "f",
          digits = 3,
          width = col_widths[2]
        )
      )
      bdy2 <- ifelse(
        is.na(results[j, "b2"]),
        formatC("NA", width = col_widths[3], format = "s"),
        formatC(
          results[j, "b2"],
          format = "f",
          digits = 3,
          width = col_widths[3]
        )
      )
      control1 <- ifelse(
        is.na(results[j, "h01"]),
        formatC("NA", width = col_widths[4], format = "s"),
        formatC(
          results[j, "h01"],
          format = "f",
          digits = 3,
          width = col_widths[4]
        )
      )
      control2 <- ifelse(
        is.na(results[j, "h02"]),
        formatC("NA", width = col_widths[5], format = "s"),
        formatC(
          results[j, "h02"],
          format = "f",
          digits = 3,
          width = col_widths[5]
        )
      )
      treatment1 <- ifelse(
        is.na(results[j, "h11"]),
        formatC("NA", width = col_widths[6], format = "s"),
        formatC(
          results[j, "h11"],
          format = "f",
          digits = 3,
          width = col_widths[6]
        )
      )
      treatment2 <- ifelse(
        is.na(results[j, "h12"]),
        formatC("NA", width = col_widths[7], format = "s"),
        formatC(
          results[j, "h12"],
          format = "f",
          digits = 3,
          width = col_widths[7]
        )
      )
      NCo <- ifelse(
        is.na(results[j, "N.Co"]),
        formatC("NA", width = col_widths[8], format = "s"),
        formatC(results[j, "N.Co"], format = "d", width = col_widths[8])
      )
      NTr <- ifelse(
        is.na(results[j, "N.Tr"]),
        formatC("NA", width = col_widths[9], format = "s"),
        formatC(results[j, "N.Tr"], format = "d", width = col_widths[9])
      )
      # Combine formatted values and print the row
      row_vals <- c(
        index, bdy1, bdy2, control1, control2, treatment1, treatment2, NCo,
        NTr
      )
      if (j %in% subset) cat(paste(row_vals, collapse = ""), "\n")
    }
    cat(strrep("=", sum(col_widths)), "\n")
    return(
      list(
        output = output,
        table = results[subset,,drop = FALSE],
        cbands = NULL
      )
    )
  }

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

  class(summary.result) <- "summary.rd2d"
  invisible(summary.result)

}
