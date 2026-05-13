################################################################################
# rd2d R Package
# Illustration: plots
################################################################################

rm(list = ls(all = TRUE))

# This script reads rd2d_illustration_results.rds and demonstrates how to build
# inference plots and heatmaps from the summary.rd2d return tables.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/"))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(dirname(frame_files[length(frame_files)]), winslash = "/"))
  }

  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required for this plotting illustration.", call. = FALSE)
}

results_file <- file.path(output_dir, "rd2d_illustration_results.rds")
if (!file.exists(results_file)) {
  stop("Run R/rd2d_illustration.R before R/rd2d_plot.R.", call. = FALSE)
}

obj <- readRDS(results_file)

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

point_rows <- function(table) {
  table[!(rownames(table) %in% c("WBATE", "LBATE")), , drop = FALSE]
}

aggregate_rows <- function(table) {
  table[rownames(table) %in% c("WBATE", "LBATE"), , drop = FALSE]
}

make_inference_plot <- function(summary_object, output, file) {
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
    ggplot2::geom_point(ggplot2::aes(color = label_for_output(output)), size = 1.9) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey45") +
    ggplot2::geom_vline(xintercept = 21, linewidth = 0.35, color = "grey80") +
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
    ggplot2::scale_fill_manual(values = c("95% CB" = "#1b6ca8"), name = NULL) +
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
      ggplot2::scale_linetype_manual(values = c(WBATE = "dashed", LBATE = "dotdash"),
                                     name = NULL)
  }

  ggplot2::ggsave(file.path(plot_dir, file), p, width = 6.5, height = 4.2)
  p
}

make_effect_heatmap <- function(summary_object, output, file) {
  table <- point_rows(summary_object$tables[[output]])
  plot_dat <- data.frame(
    b1 = table$b1,
    b2 = table$b2,
    index = seq_len(nrow(table)),
    estimate = table$estimate.p
  )

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = b1, y = b2)) +
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

  ggplot2::ggsave(file.path(plot_dir, file), p, width = 5.6, height = 4.8)
  p
}

make_pvalue_heatmap <- function(summary_object, output, file) {
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

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = b1, y = b2)) +
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

  ggplot2::ggsave(file.path(plot_dir, file), p, width = 5.6, height = 4.8)
  p
}

make_distance_inference_plot <- function(summary_object, output = "main", file,
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

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = index, y = estimate)) +
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
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey45") +
    ggplot2::geom_vline(xintercept = 21, linewidth = 0.35, color = "grey80") +
    ggplot2::labs(x = "Boundary evaluation point", y = y_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(file.path(plot_dir, file), p, width = 6.5, height = 4.2)
  p
}

plots <- list()

for (output in c("main", "itt", "fs", "itt.0")) {
  prefix <- if (output == "main") "fuzzy" else gsub("\\.", "", output)
  plots[[paste0(prefix, "_inference")]] <- make_inference_plot(
    obj$summaries[[output]],
    output,
    sprintf("rd2d_illustration_%s_inference.png", prefix)
  )
  plots[[paste0(prefix, "_heatmap")]] <- make_effect_heatmap(
    obj$summaries[[output]],
    output,
    sprintf("rd2d_illustration_%s_heatmap.png", prefix)
  )
  plots[[paste0(prefix, "_pvalue_heatmap")]] <- make_pvalue_heatmap(
    obj$summaries[[output]],
    output,
    sprintf("rd2d_illustration_%s_heatmap_pvalue.png", prefix)
  )
}

plots$distance_sharp_inference <- make_distance_inference_plot(
  obj$summaries$distance_sharp,
  output = "main",
  file = "rd2d_illustration_distance_sharp_inference.png",
  y_label = "distance-based sharp effect",
  color = "#7b3294"
)

plots$distance_fuzzy_inference <- make_distance_inference_plot(
  obj$summaries$distance_fuzzy,
  output = "main",
  file = "rd2d_illustration_distance_fuzzy_inference.png",
  y_label = "distance-based fuzzy effect",
  color = "#a6611a"
)

cat(sprintf("Illustration plots saved to: %s\n", plot_dir))
