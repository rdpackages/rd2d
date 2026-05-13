################################################################################
#' @title rd2d: Estimation and Inference for Boundary Discontinuity Designs
#'
#' @description
#' \pkg{rd2d} implements pointwise and uniform estimation and inference
#' procedures for boundary discontinuity (BD) designs using local polynomial
#' methods. The package includes location-based and distance-based methods,
#' sharp and fuzzy designs, automatic bandwidth selection, pointwise confidence
#' intervals, and uniform confidence bands.
#' Distance-based methods in this package target level effects at
#' two-dimensional boundary points.
#'
#' Included functions are: \link{rd2d} for location-based estimation and
#' inference, \link{rdbw2d} for location-based bandwidth selection,
#' \link{rd2d.distance} for distance-based estimation and inference, and
#' \link{rdbw2d.distance} for distance-based bandwidth selection.
#'
#' \code{print()} and \code{summary()} methods are available for all four
#' functions.
#'
#' Related Stata, R, and Python packages useful for inference in RD designs are
#' described at:
#'
#' \href{https://rdpackages.github.io/}{https://rdpackages.github.io/}
#'
#' For an introduction to regression discontinuity designs, see
#' Cattaneo and Titiunik (2022, \doi{10.1146/annurev-economics-051520-021409})
#' and references therein.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com}.
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}.
#' Ruiqi Rae Yu, Princeton University. \email{raeyuuuu@gmail.com}.
#'
#' @references
#' \itemize{
#' \item{Cattaneo, M. D., and Titiunik, R. (2022).
#' Regression Discontinuity Designs. \doi{10.1146/annurev-economics-051520-021409}.}
#' \item{\href{https://arxiv.org/abs/2505.05670}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods.}
#' \item{\href{https://arxiv.org/abs/2510.26051}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods.}
#' \item{\href{https://arxiv.org/abs/2511.06474}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026).}
#' Boundary Discontinuity Designs: Theory and Practice.}
#' \item{\href{https://arxiv.org/abs/2505.07989}{Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2025).}
#' rd2d: Causal Inference in Boundary Discontinuity Designs.}
#' }
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats integrate
#' @importFrom stats optimize
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @importFrom stats as.formula
#' @importFrom stats complete.cases
#' @importFrom stats cov
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom stats predict
#' @importFrom stats density
#' @importFrom MASS mvrnorm
#' @importFrom MASS ginv
#' @importFrom expm sqrtm

#' @import ggplot2
#'
#' @aliases rd2d-package
"_PACKAGE"
