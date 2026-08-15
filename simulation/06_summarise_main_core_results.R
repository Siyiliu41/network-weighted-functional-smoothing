# Summarise the completed, frozen 16-cell core simulation.
#
# This script reads the saved main result, writes analysis tables and figures,
# and never refits a model or changes the random-number-generator state.

main_summary_script_version <- "2026-08-14-core-summary-r2"

find_project_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    is_project_root <-
      dir.exists(file.path(here, "simulation")) &&
      file.exists(file.path(here, "netfunsmooth", "DESCRIPTION"))

    if (is_project_root) {
      return(here)
    }

    parent <- dirname(here)

    if (identical(parent, here)) {
      stop(
        "Could not find the project root. Expected a `simulation` directory ",
        "and `netfunsmooth/DESCRIPTION`."
      )
    }

    here <- parent
  }
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package `ggplot2` is required to generate Figure 1. ",
    "Install it with install.packages(\"ggplot2\")."
  )
}

suppressPackageStartupMessages(library(ggplot2))

project_root <- find_project_root()

main_result_path <- file.path(
  project_root,
  "simulation", "results", "core", "main", "main_core_result.rds"
)

analysis_root <- file.path(
  project_root,
  "simulation", "results", "core", "main", "analysis"
)

table_directory <- file.path(analysis_root, "tables")
figure_directory <- file.path(analysis_root, "figures")

if (!file.exists(main_result_path)) {
  stop("Cannot find completed main result: ", main_result_path)
}

main_result <- readRDS(main_result_path)

required_entries <- c(
  "config",
  "quality_checks",
  "replication_summary",
  "completed_at",
  "main_version"
)

if (!is.list(main_result) || !all(required_entries %in% names(main_result))) {
  stop("The main result has an incompatible structure.")
}

if (!isTRUE(all(main_result$quality_checks))) {
  stop(
    "The saved main result did not pass all quality checks; ",
    "pause interpretation."
  )
}

required_config_entries <- c(
  "n_replications",
  "primary_methods",
  "maximum_failure_rate"
)

if (!is.list(main_result$config) ||
    !all(required_config_entries %in% names(main_result$config))) {
  stop("The saved main-result configuration is incomplete.")
}

replication_summary <- main_result$replication_summary

required_columns <- c(
  "cell_id",
  "truth_structure",
  "alpha",
  "neighbourhood",
  "sigma",
  "replication",
  "method",
  "success",
  "warning_count",
  "node_averaged_ise",
  "median_node_ise",
  "worst_node_ise",
  "elapsed_seconds"
)

if (!is.data.frame(replication_summary) ||
    !all(required_columns %in% names(replication_summary))) {
  stop("`replication_summary` is missing required columns.")
}

expected_attempted <- main_result$config$n_replications

if (!identical(as.integer(expected_attempted), 200L)) {
  stop("This summary script is reserved for the frozen B = 200 main result.")
}

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

finite_mean <- function(x) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  mean(x)
}

finite_median <- function(x) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  stats::median(x)
}

finite_mcse <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  stats::sd(x) / sqrt(length(x))
}

make_method_summary <- function(data) {
  groups <- split(
    data,
    interaction(data$cell_id, data$method, drop = TRUE, lex.order = TRUE)
  )

  out <- lapply(groups, function(x) {
    successful <- x$success %in% TRUE & is.finite(x$node_averaged_ise)
    elapsed <- x$elapsed_seconds[successful]

    data.frame(
      cell_id = x$cell_id[[1L]],
      truth_structure = x$truth_structure[[1L]],
      alpha = x$alpha[[1L]],
      neighbourhood = x$neighbourhood[[1L]],
      sigma = x$sigma[[1L]],
      method = x$method[[1L]],
      attempted = nrow(x),
      successful = sum(successful),
      failure_rate = 1 - mean(successful),
      total_warnings = sum(x$warning_count, na.rm = TRUE),
      mise = finite_mean(x$node_averaged_ise[successful]),
      mise_mcse = finite_mcse(x$node_averaged_ise[successful]),
      median_aise = finite_median(x$node_averaged_ise[successful]),
      mean_median_node_ise = finite_mean(x$median_node_ise[successful]),
      median_median_node_ise = finite_median(x$median_node_ise[successful]),
      mean_worst_node_ise = finite_mean(x$worst_node_ise[successful]),
      median_worst_node_ise = finite_median(x$worst_node_ise[successful]),
      mean_elapsed_seconds = finite_mean(elapsed),
      median_elapsed_seconds = finite_median(elapsed),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  result[order(
    result$truth_structure,
    result$alpha,
    result$neighbourhood,
    result$sigma,
    result$method
  ), ]
}

make_network_paired_summary <- function(data) {
  cells <- sort(unique(data$cell_id))

  out <- lapply(cells, function(cell) {
    x <- data[data$cell_id == cell, , drop = FALSE]
    nodewise <- x[x$method == "nodewise", , drop = FALSE]
    network <- x[x$method == "network", , drop = FALSE]

    paired <- merge(
      nodewise[, c("replication", "node_averaged_ise", "success")],
      network[, c("replication", "node_averaged_ise", "success")],
      by = "replication",
      all = TRUE,
      sort = TRUE,
      suffixes = c("_nodewise", "_network")
    )

    successful <- paired$success_nodewise %in% TRUE &
      paired$success_network %in% TRUE &
      is.finite(paired$node_averaged_ise_nodewise) &
      is.finite(paired$node_averaged_ise_network)

    paired <- paired[successful, , drop = FALSE]

    delta <- paired$node_averaged_ise_nodewise -
      paired$node_averaged_ise_network

    design <- x[1L, c(
      "cell_id",
      "truth_structure",
      "alpha",
      "neighbourhood",
      "sigma"
    )]

    data.frame(
      design,
      attempted = nrow(nodewise),
      jointly_successful = length(delta),
      mean_paired_improvement = finite_mean(delta),
      paired_mcse = finite_mcse(delta),
      lower_two_mcse = finite_mean(delta) - 2 * finite_mcse(delta),
      upper_two_mcse = finite_mean(delta) + 2 * finite_mcse(delta),
      win_rate = if (length(delta)) mean(delta > 0) else NA_real_,
      mise_ratio_network_vs_nodewise = if (length(delta)) {
        sum(paired$node_averaged_ise_network) /
          sum(paired$node_averaged_ise_nodewise)
      } else {
        NA_real_
      },
      relative_improvement_percent = if (length(delta)) {
        100 * (
          1 -
            sum(paired$node_averaged_ise_network) /
            sum(paired$node_averaged_ise_nodewise)
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  result[order(
    result$truth_structure,
    result$alpha,
    result$neighbourhood,
    result$sigma
  ), ]
}

make_rook_queen_summary <- function(data) {
  data <- data[data$method == "network", , drop = FALSE]

  keys <- c("truth_structure", "alpha", "sigma")

  groups <- split(
    data,
    interaction(data[keys], drop = TRUE, lex.order = TRUE)
  )

  out <- lapply(groups, function(x) {
    rook <- x[x$neighbourhood == "rook", , drop = FALSE]
    queen <- x[x$neighbourhood == "queen", , drop = FALSE]

    paired <- merge(
      rook[, c("replication", "node_averaged_ise", "success")],
      queen[, c("replication", "node_averaged_ise", "success")],
      by = "replication",
      all = TRUE,
      sort = TRUE,
      suffixes = c("_rook", "_queen")
    )

    successful <- paired$success_rook %in% TRUE &
      paired$success_queen %in% TRUE &
      is.finite(paired$node_averaged_ise_rook) &
      is.finite(paired$node_averaged_ise_queen)

    paired <- paired[successful, , drop = FALSE]

    difference <- paired$node_averaged_ise_rook -
      paired$node_averaged_ise_queen

    data.frame(
      truth_structure = x$truth_structure[[1L]],
      alpha = x$alpha[[1L]],
      sigma = x$sigma[[1L]],
      attempted = nrow(rook),
      jointly_successful = length(difference),
      mean_rook_minus_queen_aise = finite_mean(difference),
      paired_mcse = finite_mcse(difference),
      queen_win_rate = if (length(difference)) {
        mean(difference > 0)
      } else {
        NA_real_
      },
      queen_to_rook_mise_ratio = if (length(difference)) {
        sum(paired$node_averaged_ise_queen) /
          sum(paired$node_averaged_ise_rook)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  result[order(result$truth_structure, result$alpha, result$sigma), ]
}

method_summary <- make_method_summary(replication_summary)

network_paired_summary <- make_network_paired_summary(replication_summary)

rook_queen_summary <- make_rook_queen_summary(replication_summary)

primary_method_rows <- method_summary$method %in%
  main_result$config$primary_methods

analysis_checks <- c(
  saved_main_quality_checks_passed = isTRUE(all(main_result$quality_checks)),
  expected_16_core_cells = length(unique(replication_summary$cell_id)) == 16L,
  expected_12800_method_rows = nrow(replication_summary) == 16L * 200L * 4L,
  no_duplicate_cell_replication_method = !anyDuplicated(
    replication_summary[, c("cell_id", "replication", "method")]
  ),
  all_primary_methods_present = all(
    main_result$config$primary_methods %in% unique(replication_summary$method)
  ),
  all_method_cells_have_200_attempts = all(method_summary$attempted == 200L),
  primary_failure_rates_at_most_5_percent = isTRUE(all(
    method_summary$failure_rate[primary_method_rows] <=
      main_result$config$maximum_failure_rate
  )),
  all_network_pairs_have_200_successes = all(
    network_paired_summary$jointly_successful == 200L
  ),
  all_rook_queen_pairs_have_200_successes = all(
    rook_queen_summary$jointly_successful == 200L
  )
)

if (!isTRUE(all(analysis_checks))) {
  stop(
    "Post-main summary audit failed: ",
    paste(names(analysis_checks)[!analysis_checks], collapse = ", "),
    "."
  )
}

utils::write.csv(
  method_summary,
  file.path(table_directory, "table_1_method_performance.csv"),
  row.names = FALSE
)

utils::write.csv(
  network_paired_summary,
  file.path(table_directory, "table_2_network_vs_nodewise.csv"),
  row.names = FALSE
)

utils::write.csv(
  rook_queen_summary,
  file.path(table_directory, "table_3_rook_vs_queen_network.csv"),
  row.names = FALSE
)

# Figure 1: relative MISE reduction of network smoothing versus nodewise
# smoothing.

figure_1 <- file.path(
  figure_directory,
  "figure_1_network_relative_improvement.png"
)

plot_data <- network_paired_summary

figure_1_checks <- c(
  expected_16_rows = nrow(plot_data) == 16L,
  all_network_pairs_successful = all(plot_data$jointly_successful == 200L),
  expected_truth_structures = identical(
    sort(unique(plot_data$truth_structure)),
    c("cluster", "coordinate")
  ),
  expected_noise_levels = identical(
    sort(unique(plot_data$sigma)),
    c(0.2, 0.5)
  ),
  expected_neighbourhoods = identical(
    sort(unique(plot_data$neighbourhood)),
    c("queen", "rook")
  ),
  finite_relative_improvements = all(
    is.finite(plot_data$relative_improvement_percent)
  ),
  values_fit_fixed_y_axis = all(
    plot_data$relative_improvement_percent >= -0.5 &
      plot_data$relative_improvement_percent <= 21.5
  )
)

if (!isTRUE(all(figure_1_checks))) {
  stop(
    "Figure 1 input check failed: ",
    paste(names(figure_1_checks)[!figure_1_checks], collapse = ", "),
    "."
  )
}

plot_data$truth_label <- factor(
  plot_data$truth_structure,
  levels = c("cluster", "coordinate"),
  labels = c("Cluster-structured truth", "Coordinate-smooth truth")
)

plot_data$sigma_label <- factor(
  plot_data$sigma,
  levels = c(0.2, 0.5),
  labels = c("Low noise: sigma = 0.2", "High noise: sigma = 0.5")
)

plot_data$graph_label <- factor(
  plot_data$neighbourhood,
  levels = c("rook", "queen"),
  labels = c("Rook graph", "Queen graph")
)

figure_1_plot <- ggplot(
  plot_data,
  aes(
    x = alpha,
    y = relative_improvement_percent,
    colour = graph_label,
    linetype = graph_label,
    shape = graph_label,
    group = graph_label
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "grey55",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(
    rows = vars(truth_label),
    cols = vars(sigma_label)
  ) +
  scale_x_continuous(
    breaks = c(0.5, 0.9),
    limits = c(0.45, 0.95)
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, by = 5),
    limits = c(-0.5, 21.5),
    labels = function(x) paste0(x, "%")
  ) +
  scale_colour_manual(
    values = c(
      "Rook graph" = "#0072B2",
      "Queen graph" = "#D55E00"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Rook graph" = "solid",
      "Queen graph" = "dashed"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Rook graph" = 16,
      "Queen graph" = 17
    )
  ) +
  labs(
    title = "Relative MISE reduction of network-weighted smoothing",
    x = expression(alpha),
    y = "Relative MISE reduction (%)",
    colour = "Supplied graph",
    linetype = "Supplied graph",
    shape = "Supplied graph"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

ggplot2::ggsave(
  filename = figure_1,
  plot = figure_1_plot,
  width = 7.3,
  height = 5.5,
  units = "in",
  dpi = 300
)

# Figure 2: absolute MISE of all four estimators in every core cell.

figure_2 <- file.path(
  figure_directory,
  "figure_2_method_mise_by_cell.png"
)

grDevices::png(
  figure_2,
  width = 2000,
  height = 1050,
  res = 180
)

old_par <- graphics::par(no.readonly = TRUE)

graphics::par(
  mar = c(10, 4.5, 2.5, 1.4),
  las = 2
)

cell_order <- unique(network_paired_summary$cell_id)

methods <- c("raw", "pooled", "nodewise", "network")

colours <- c(
  raw = "grey45",
  pooled = "#CC79A7",
  nodewise = "#009E73",
  network = "#0072B2"
)

bar_matrix <- vapply(methods, function(method) {
  x <- method_summary[
    method_summary$method == method,
    c("cell_id", "mise")
  ]

  x$mise[match(cell_order, x$cell_id)]
}, numeric(length(cell_order)))

bar_matrix <- t(bar_matrix)

graphics::barplot(
  bar_matrix,
  beside = TRUE,
  col = colours[methods],
  border = NA,
  names.arg = cell_order,
  ylab = "MISE",
  main = "Absolute reconstruction error by core cell"
)

graphics::legend(
  "topright",
  legend = methods,
  fill = colours[methods],
  bty = "n",
  horiz = TRUE
)

graphics::par(old_par)
grDevices::dev.off()

analysis_result <- list(
  analysis_version = main_summary_script_version,
  r_version = R.version.string,
  ggplot2_version = as.character(utils::packageVersion("ggplot2")),
  session_info = utils::sessionInfo(),
  input_path = main_result_path,
  input_completed_at = main_result$completed_at,
  input_main_version = main_result$main_version,
  analysis_checks = analysis_checks,
  figure_1_checks = figure_1_checks,
  method_summary = method_summary,
  network_paired_summary = network_paired_summary,
  rook_queen_summary = rook_queen_summary,
  table_paths = list.files(table_directory, full.names = TRUE),
  figure_paths = c(figure_1, figure_2),
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

saveRDS(
  analysis_result,
  file.path(analysis_root, "main_core_analysis_result.rds"),
  version = 3
)

cat("\nMain core-result analysis checks:\n")
print(analysis_checks)

cat("\nFigure 1 input checks:\n")
print(figure_1_checks)

cat("\nNetwork versus nodewise, all 16 cells:\n")
print(network_paired_summary)

cat("\nOutputs written to:\n", analysis_root, "\n", sep = "")