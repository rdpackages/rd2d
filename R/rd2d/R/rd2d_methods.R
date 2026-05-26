if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "estimate", "lower", "upper", "bandwidth", "series"))
}

rd2d_table <- function(x, output = "main") {
  if (!is.character(output) || length(output) != 1 || is.na(output)) {
    stop("output must be a single character string.", call. = FALSE)
  }
  tab <- x[[output]]
  if (!is.data.frame(tab)) {
    available <- names(x)[vapply(x, is.data.frame, logical(1))]
    stop(
      paste0(
        "output must name an available result table. Available outputs: ",
        paste(available, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  tab
}

rd2d_vcov_table <- function(x, output = "main") {
  if (!is.character(output) || length(output) != 1 || is.na(output)) {
    stop("output must be a single character string.", call. = FALSE)
  }
  cov <- x$params.cov[[output]]
  if (is.null(cov)) {
    stop(
      paste0(
        "vcov requires the requested covariance matrix. Re-run the fit with ",
        "params.cov including \"", output, "\"."
      ),
      call. = FALSE
    )
  }
  cov
}

rd2d_order_suffix <- function(order) {
  order <- match.arg(order, c("q", "p"))
  paste0(".", order)
}

rd2d_eval_axis <- function(tab) {
  n <- nrow(tab)
  if (!all(c("b1", "b2") %in% names(tab))) {
    return(list(x = seq_len(n), label = "Evaluation point"))
  }
  b1 <- tab$b1
  b2 <- tab$b2
  unique_finite <- function(x) unique(x[is.finite(x)])
  if (length(unique_finite(b1)) > 1 && length(unique_finite(b2)) <= 1) {
    list(x = b1, label = "b1")
  } else if (length(unique_finite(b2)) > 1 && length(unique_finite(b1)) <= 1) {
    list(x = b2, label = "b2")
  } else {
    list(x = seq_len(n), label = "Evaluation point")
  }
}

rd2d_confint_matrix <- function(object, output = "main", order = c("q", "p"),
                                level = 0.95, side = NULL) {
  suffix <- rd2d_order_suffix(order)
  tab <- rd2d_table(object, output)
  est.col <- paste0("estimate", suffix)
  se.col <- paste0("std.err", suffix)
  if (!all(c(est.col, se.col) %in% names(tab))) {
    stop("Selected output does not contain estimates and standard errors.",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1 ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be a numeric value between 0 and 1.", call. = FALSE)
  }
  if (is.null(side)) side <- object$opt$side
  if (is.null(side)) side <- "two"
  side <- match.arg(side, c("two", "left", "right"))

  estimate <- tab[[est.col]]
  se <- tab[[se.col]]
  if (side == "two") {
    alpha <- 1 - level
    zval <- qnorm(1 - alpha / 2)
    lower <- estimate - zval * se
    upper <- estimate + zval * se
    colnames.out <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2)), " %")
  } else if (side == "left") {
    zval <- qnorm(level)
    lower <- rep(-Inf, length(estimate))
    upper <- estimate + zval * se
    colnames.out <- c("lower", paste0(format(100 * level), " %"))
  } else {
    zval <- qnorm(level)
    lower <- estimate - zval * se
    upper <- rep(Inf, length(estimate))
    colnames.out <- c(paste0(format(100 * (1 - level)), " %"), "upper")
  }

  out <- cbind(lower, upper)
  colnames(out) <- colnames.out
  rownames(out) <- if (!is.null(rownames(tab))) rownames(tab) else seq_along(estimate)
  out
}

rd2d_coef <- function(object, output = "main", order = c("q", "p")) {
  suffix <- rd2d_order_suffix(order)
  tab <- rd2d_table(object, output)
  est.col <- paste0("estimate", suffix)
  if (!est.col %in% names(tab)) {
    stop("Selected output does not contain estimates.", call. = FALSE)
  }
  out <- tab[[est.col]]
  names(out) <- if (!is.null(rownames(tab))) rownames(tab) else seq_along(out)
  out
}

rd2d_plot_estimates <- function(x, output = "main", order = c("q", "p"),
                                ci = TRUE, level = NULL, draw = TRUE) {
  suffix <- rd2d_order_suffix(order)
  tab <- rd2d_table(x, output)
  est.col <- paste0("estimate", suffix)
  if (!est.col %in% names(tab)) {
    stop("Selected output does not contain estimates.", call. = FALSE)
  }
  axis <- rd2d_eval_axis(tab)
  plot.dat <- data.frame(
    x = axis$x,
    estimate = tab[[est.col]]
  )
  if (isTRUE(ci)) {
    if (is.null(level)) {
      level <- if (!is.null(x$opt$level)) x$opt$level / 100 else 0.95
    }
    ci.mat <- rd2d_confint_matrix(x, output = output, order = sub("^\\.", "", suffix),
                                  level = level)
    plot.dat$lower <- ci.mat[, 1]
    plot.dat$upper <- ci.mat[, 2]
  }

  p <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = x, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
                        color = "grey55") +
    ggplot2::labs(x = axis$label, y = "Estimate")
  if (nrow(plot.dat) > 1) {
    p <- p + ggplot2::geom_line(linewidth = 0.4, color = "#1f77b4")
  }
  if (isTRUE(ci)) {
    p <- p + ggplot2::geom_pointrange(
      ggplot2::aes(ymin = lower, ymax = upper),
      linewidth = 0.35, fatten = 1.8, color = "#1f77b4"
    )
  } else {
    p <- p + ggplot2::geom_point(size = 1.8, color = "#1f77b4")
  }
  p <- p + ggplot2::theme_minimal()
  if (isTRUE(draw)) print(p)
  invisible(p)
}

rd2d_plot_bandwidths <- function(x, distance = FALSE, draw = TRUE) {
  tab <- x$bws
  h.cols <- if (isTRUE(distance)) c("h0", "h1") else c("h01", "h02", "h11", "h12")
  h.cols <- intersect(h.cols, names(tab))
  if (length(h.cols) == 0) {
    stop("No bandwidth columns found to plot.", call. = FALSE)
  }
  labels <- if (isTRUE(distance)) {
    c(h0 = "Control", h1 = "Treatment")[h.cols]
  } else {
    c(
      h01 = "Control x1", h02 = "Control x2",
      h11 = "Treatment x1", h12 = "Treatment x2"
    )[h.cols]
  }
  axis <- rd2d_eval_axis(tab)
  plot.dat <- data.frame(
    x = rep(axis$x, times = length(h.cols)),
    bandwidth = unlist(tab[h.cols], use.names = FALSE),
    series = rep(unname(labels), each = nrow(tab))
  )
  p <- ggplot2::ggplot(
    plot.dat,
    ggplot2::aes(x = x, y = bandwidth, color = series)
  ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::labs(x = axis$label, y = "Bandwidth", color = NULL) +
    ggplot2::theme_minimal()
  if (isTRUE(draw)) print(p)
  invisible(p)
}

#' Extract rd2d Point Estimates
#'
#' @param object An object returned by \code{\link{rd2d}} or
#'   \code{\link{rd2d.distance}}.
#' @param output Result table to extract. Default is \code{"main"}.
#' @param order Polynomial order estimate to extract: \code{"q"} for the
#'   inference order or \code{"p"} for the estimation order.
#' @param ... Additional arguments ignored.
#'
#' @return A named numeric vector of estimates.
#' @export
coef.rd2d <- function(object, output = "main", order = c("q", "p"), ...) {
  rd2d_coef(object, output = output, order = order)
}

#' @rdname coef.rd2d
#' @export
coef.rd2d.distance <- function(object, output = "main", order = c("q", "p"), ...) {
  rd2d_coef(object, output = output, order = order)
}

#' Extract rd2d Covariance Matrices
#'
#' @param object An object returned by \code{\link{rd2d}} or
#'   \code{\link{rd2d.distance}}.
#' @param output Covariance output to extract. Default is \code{"main"}.
#' @param ... Additional arguments ignored.
#'
#' @return A covariance matrix.
#' @export
vcov.rd2d <- function(object, output = "main", ...) {
  rd2d_vcov_table(object, output = output)
}

#' @rdname vcov.rd2d
#' @export
vcov.rd2d.distance <- function(object, output = "main", ...) {
  rd2d_vcov_table(object, output = output)
}

#' Confidence Intervals for rd2d Estimates
#'
#' @param object An object returned by \code{\link{rd2d}} or
#'   \code{\link{rd2d.distance}}.
#' @param parm Ignored; included for compatibility with \code{\link{confint}}.
#' @param level Confidence level in \eqn{(0, 1)}. Default is \code{0.95}.
#' @param output Result table to use. Default is \code{"main"}.
#' @param order Polynomial order estimate to use: \code{"q"} or \code{"p"}.
#' @param side Confidence interval side. Defaults to the side stored in the fit.
#' @param ... Additional arguments ignored.
#'
#' @return A two-column matrix with lower and upper confidence limits.
#' @export
confint.rd2d <- function(object, parm, level = 0.95, output = "main",
                         order = c("q", "p"), side = NULL, ...) {
  rd2d_confint_matrix(
    object, output = output, order = order, level = level, side = side
  )
}

#' @rdname confint.rd2d
#' @export
confint.rd2d.distance <- function(object, parm, level = 0.95, output = "main",
                                  order = c("q", "p"), side = NULL, ...) {
  rd2d_confint_matrix(
    object, output = output, order = order, level = level, side = side
  )
}

#' Plot rd2d Estimates and Bandwidths
#'
#' @param x An object returned by \code{\link{rd2d}},
#'   \code{\link{rd2d.distance}}, \code{\link{rdbw2d}}, or
#'   \code{\link{rdbw2d.distance}}.
#' @param output Result table to plot for estimation objects. Default is
#'   \code{"main"}.
#' @param order Polynomial order estimate to plot for estimation objects:
#'   \code{"q"} or \code{"p"}.
#' @param ci Logical. If \code{TRUE}, include pointwise confidence intervals for
#'   estimation objects.
#' @param level Confidence level for plotted intervals. Defaults to the level
#'   stored in the fit, or \code{0.95} if unavailable.
#' @param draw Logical. If \code{TRUE}, print the ggplot object.
#' @param ... Additional arguments ignored.
#'
#' @return Invisibly, a \code{ggplot} object.
#' @export
plot.rd2d <- function(x, output = "main", order = c("q", "p"), ci = TRUE,
                      level = NULL, draw = TRUE, ...) {
  rd2d_plot_estimates(
    x, output = output, order = order, ci = ci, level = level, draw = draw
  )
}

#' @rdname plot.rd2d
#' @export
plot.rd2d.distance <- function(x, output = "main", order = c("q", "p"),
                               ci = TRUE, level = NULL, draw = TRUE, ...) {
  rd2d_plot_estimates(
    x, output = output, order = order, ci = ci, level = level, draw = draw
  )
}

#' @rdname plot.rd2d
#' @export
plot.rdbw2d <- function(x, draw = TRUE, ...) {
  rd2d_plot_bandwidths(x, distance = FALSE, draw = draw)
}

#' @rdname plot.rd2d
#' @export
plot.rdbw2d.distance <- function(x, draw = TRUE, ...) {
  rd2d_plot_bandwidths(x, distance = TRUE, draw = draw)
}
