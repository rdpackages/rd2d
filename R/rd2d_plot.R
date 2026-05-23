################################################################################
# RD2D Package
# Plot Illustration
################################################################################

rm(list = ls(all = TRUE))

library(rd2d)
library(ggplot2)

# Generate boundary evaluation points.
make_eval_grid <- function(neval = 40) {
  half <- ceiling(neval / 2)
  rbind(
    data.frame(
      x.1 = rep(0, half),
      x.2 = 40 - (seq_len(half) - 1) * 40 / half
    ),
    data.frame(
      x.1 = (seq_len(neval - half) - 1) * 56 / half,
      x.2 = rep(0, neval - half)
    )
  )
}

# Generate signed distances to each boundary evaluation point.
make_signed_distances <- function(X, eval, assignment) {
  distance <- sapply(seq_len(nrow(eval)), function(j) {
    sqrt((X$x.1 - eval$x.1[j])^2 + (X$x.2 - eval$x.2[j])^2)
  })
  distance * matrix(
    2 * assignment - 1,
    nrow = nrow(distance),
    ncol = ncol(distance)
  )
}

# Label each estimand for plots.
label_for_output <- function(output) {
  switch(
    output,
    main = "Fuzzy",
    itt = "ITT",
    fs = "First stage",
    itt.0 = "ITT, control side",
    output
  )
}

# Keep boundary rows separate from aggregate rows.
point_rows <- function(table) {
  table[!(rownames(table) %in% c("WBATE", "LBATE")), , drop = FALSE]
}

# Keep WBATE and LBATE rows for reference lines.
aggregate_rows <- function(table) {
  table[rownames(table) %in% c("WBATE", "LBATE"), , drop = FALSE]
}

# Compute summary tables without printing them.
summary_for_plot <- function(...) {
  invisible(utils::capture.output(out <- summary(...)))
  out
}

# Plot point estimates, CIs, and CBs.
make_inference_plot <- function(summary_object, output) {
  table <- summary_object$tables[[output]]
  points <- point_rows(table)
  aggregates <- aggregate_rows(table)

  plot_dat <- data.frame(
    index = seq_len(nrow(points)),
    estimate = points$estimate.p,
    ci.lower = points$ci.lower,
    ci.upper = points$ci.upper,
    cb.lower = points$cb.lower,
    cb.upper = points$cb.upper
  )

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = index, y = estimate)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cb.lower, ymax = cb.upper, fill = "95% CB"),
      alpha = 0.16
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci.lower, ymax = ci.upper, color = "95% CI"),
      width = 0,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = label_for_output(output)),
      size = 1.9
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "grey45"
    ) +
    ggplot2::geom_vline(
      xintercept = 21,
      linewidth = 0.35,
      color = "grey80"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Fuzzy" = "#1b6ca8",
        "ITT" = "#1b6ca8",
        "First stage" = "#1b6ca8",
        "ITT, control side" = "#1b6ca8",
        "95% CI" = "#1b6ca8"
      ),
      breaks = c(label_for_output(output), "95% CI"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("95% CB" = "#1b6ca8"),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Boundary evaluation point",
      y = label_for_output(output)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 5, 10, 15, 21, 25, 30, 35, 40)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  if (nrow(aggregates) > 0) {
    aggregate_dat <- data.frame(
      aggregate = rownames(aggregates),
      estimate = aggregates$estimate.p
    )
    p <- p +
      ggplot2::geom_hline(
        data = aggregate_dat,
        ggplot2::aes(yintercept = estimate, linetype = aggregate),
        color = "grey30",
        linewidth = 0.35
      ) +
      ggplot2::scale_linetype_manual(
        values = c(WBATE = "dashed", LBATE = "dotdash"),
        name = NULL
      )
  }

  p
}

# Plot point estimates over the boundary.
make_effect_heatmap <- function(summary_object, output) {
  table <- point_rows(summary_object$tables[[output]])
  plot_dat <- data.frame(
    b1 = table$b1,
    b2 = table$b2,
    index = seq_len(nrow(table)),
    estimate = table$estimate.p
  )

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = b1, y = b2)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = estimate),
      width = 2.7,
      height = 2.7,
      color = "white",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(ggplot2::aes(label = index), size = 2.6) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      name = label_for_output(output)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Score 1", y = "Score 2") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.82, 0.82),
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

# Plot p-values over the boundary.
make_pvalue_heatmap <- function(summary_object, output) {
  table <- point_rows(summary_object$tables[[output]])
  plot_dat <- data.frame(
    b1 = table$b1,
    b2 = table$b2,
    index = seq_len(nrow(table)),
    p.value = table$p.value
  )
  plot_dat$p.group <- cut(
    plot_dat$p.value,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.10, Inf),
    labels = c("0.000", "< 0.010", "< 0.050", "< 0.100", ">= 0.100")
  )

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = b1, y = b2)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = p.group),
      width = 2.7,
      height = 2.7,
      color = "white",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(ggplot2::aes(label = index), size = 2.6) +
    ggplot2::scale_fill_manual(
      values = c(
        "0.000" = "#b2182b",
        "< 0.010" = "#ef8a62",
        "< 0.050" = "#fddbc7",
        "< 0.100" = "#d1e5f0",
        ">= 0.100" = "#67a9cf"
      ),
      drop = FALSE,
      name = "p-value"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Score 1", y = "Score 2") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.82, 0.82),
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

# Plot distance-based point estimates, CIs, and CBs.
make_distance_inference_plot <- function(summary_object,
                                         output = "main",
                                         y_label = "distance-based effect",
                                         color = "#7b3294") {
  table <- point_rows(summary_object$tables[[output]])
  plot_dat <- data.frame(
    index = seq_len(nrow(table)),
    estimate = table$estimate.p,
    ci.lower = table$ci.lower,
    ci.upper = table$ci.upper,
    cb.lower = table$cb.lower,
    cb.upper = table$cb.upper
  )

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = index, y = estimate)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cb.lower, ymax = cb.upper),
      fill = color,
      alpha = 0.14
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci.lower, ymax = ci.upper),
      color = color,
      width = 0,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(color = color, size = 1.9) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "grey45"
    ) +
    ggplot2::geom_vline(
      xintercept = 21,
      linewidth = 0.35,
      color = "grey80"
    ) +
    ggplot2::labs(x = "Boundary evaluation point", y = y_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# Load simulated data.
dat <- read.csv("rd2d_data.csv")
X <- dat[, c("x.1", "x.2")]
Y <- dat$Y
A <- dat$assignment
D <- dat$fuzzy

# Set plot inputs.
neval <- as.integer(Sys.getenv("RD2D_ILLUSTRATION_NEVAL", "40"))
repp <- as.integer(Sys.getenv("RD2D_ILLUSTRATION_REPP", "499"))
eval <- make_eval_grid(neval)
distance <- make_signed_distances(X, eval, A)
wbate_weights <- rep(1, neval)

# Location-based fuzzy estimation.
fit_location <- rd2d(
  Y, X, A, eval,
  fuzzy = D,
  params.other = "itt.0",
  params.cov = c("main", "itt", "fs", "itt.0"),
  masspoints = "off"
)

# Distance-based sharp estimation.
fit_distance <- rd2d.distance(
  Y,
  distance = distance,
  b = eval,
  masspoints = "off",
  cbands = TRUE
)

# Distance-based fuzzy estimation.
fit_distance_fuzzy <- rd2d.distance(
  Y, distance = distance, b = eval,
  fuzzy = D, bwparam = "itt",
  params.cov = c("main", "itt", "fs"),
  masspoints = "off"
)

# Compute summaries used by the plots.
summaries <- list(
  main = summary_for_plot(
    fit_location, output = "main", cbands = "main",
    WBATE = wbate_weights, LBATE = TRUE, repp = repp
  ),
  itt = summary_for_plot(
    fit_location, output = "itt", cbands = "itt",
    WBATE = wbate_weights, LBATE = TRUE, repp = repp
  ),
  fs = summary_for_plot(
    fit_location, output = "fs", cbands = "fs",
    WBATE = wbate_weights, LBATE = TRUE, repp = repp
  ),
  itt.0 = summary_for_plot(
    fit_location, output = "itt.0", cbands = "itt.0",
    WBATE = wbate_weights, LBATE = TRUE, repp = repp
  ),
  distance_sharp = summary_for_plot(
    fit_distance, output = "main", cbands = "main", repp = repp
  ),
  distance_fuzzy = summary_for_plot(
    fit_distance_fuzzy, output = "main", cbands = "main",
    WBATE = wbate_weights, LBATE = TRUE, repp = repp
  )
)

# Build plots.
plots <- list()

for (output in c("main", "itt", "fs", "itt.0")) {
  prefix <- if (output == "main") "fuzzy" else gsub("\\.", "", output)
  plots[[paste0(prefix, "_inference")]] <- make_inference_plot(
    summaries[[output]], output
  )
  plots[[paste0(prefix, "_heatmap")]] <- make_effect_heatmap(
    summaries[[output]], output
  )
  plots[[paste0(prefix, "_pvalue_heatmap")]] <- make_pvalue_heatmap(
    summaries[[output]], output
  )
}

plots$distance_sharp_inference <- make_distance_inference_plot(
  summaries$distance_sharp,
  output = "main",
  y_label = "distance-based sharp effect",
  color = "#7b3294"
)

plots$distance_fuzzy_inference <- make_distance_inference_plot(
  summaries$distance_fuzzy,
  output = "main",
  y_label = "distance-based fuzzy effect",
  color = "#a6611a"
)

if (interactive()) {
  invisible(lapply(plots, print))
}
