# Prespecified cluster-boundary difference-in-differences analysis.
#
# Run from the repository root after the B = 200 main run has completed:
#   source("simulation/07_cluster_boundary_did.R")
#
# This script only reads completed main checkpoints. It never refits a model,
# changes an RNG stream, or uses pilot / known-answer replications.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package `ggplot2` is required to generate the boundary figure. ",
    "Install it with install.packages(\"ggplot2\")."
  )
}

suppressPackageStartupMessages(library(ggplot2))

boundary_did_script_version <- "2026-08-14-cluster-boundary-did-r2"

main_result_path <- "simulation/results/core/main/main_core_result.rds"
analysis_root <- "simulation/results/core/main/analysis"
table_directory <- file.path(analysis_root, "tables")
figure_directory <- file.path(analysis_root, "figures")

if (!file.exists(main_result_path)) {
  stop("Cannot find completed main result: ", main_result_path)
}
main_result <- readRDS(main_result_path)

required_entries <- c(
  "config", "quality_checks", "replication_summary", "checkpoint_paths"
)
if (!is.list(main_result) || !all(required_entries %in% names(main_result))) {
  stop("The main result has an incompatible structure.")
}
if (!isTRUE(all(main_result$quality_checks))) {
  stop("The saved main result did not pass all quality checks; pause interpretation.")
}
if (!identical(as.integer(main_result$config$n_replications), 200L)) {
  stop("This analysis is reserved for the frozen B = 200 main result.")
}

checkpoint_paths <- main_result$checkpoint_paths
if (!is.character(checkpoint_paths) || length(checkpoint_paths) != 200L ||
    !all(file.exists(checkpoint_paths))) {
  stop("All 200 completed main checkpoints are required for boundary DiD analysis.")
}

# On the column-wise 5 x 5 lattice, nodes 11:15 are column 3 (boundary)
# and nodes 16:20 are column 4 (matched interior).
boundary_nodes <- as.character(11:15)
interior_nodes <- as.character(16:20)
expected_nodes <- as.character(seq_len(25L))

finite_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

get_node_ise <- function(cell_result, method) {
  value <- cell_result$node_ise[[method]]
  if (!is.numeric(value) || length(value) != 25L || any(!is.finite(value))) {
    stop("Invalid node-level ISE for method ", method, ".")
  }
  if (is.null(names(value)) || !setequal(names(value), expected_nodes)) {
    stop("Node-level ISE has an incompatible node ordering or naming scheme.")
  }
  value[expected_nodes]
}

extract_boundary_did <- function(checkpoint) {
  result <- checkpoint$result
  if (!is.list(result) || !is.list(result$cell_results) ||
      !is.data.frame(result$summary)) {
    stop("A main checkpoint has an incompatible result structure.")
  }

  cluster_design <- result$summary[
    result$summary$truth_structure == "cluster" &
      result$summary$method == "network",
    c("cell_id", "truth_structure", "alpha", "neighbourhood", "sigma"),
    drop = FALSE
  ]
  if (nrow(cluster_design) != 8L || anyDuplicated(cluster_design$cell_id)) {
    stop("Each checkpoint must contain exactly eight cluster-network cells.")
  }

  lapply(seq_len(nrow(cluster_design)), function(i) {
    design <- cluster_design[i, , drop = FALSE]
    cell_result <- result$cell_results[[design$cell_id]]
    if (is.null(cell_result) || !is.data.frame(cell_result$summary)) {
      stop("Missing detailed result for cell ", design$cell_id, ".")
    }

    cell_summary <- cell_result$summary
    nodewise_ok <- cell_summary$success[cell_summary$method == "nodewise"] %in% TRUE
    network_ok <- cell_summary$success[cell_summary$method == "network"] %in% TRUE
    jointly_successful <- length(nodewise_ok) == 1L &&
      length(network_ok) == 1L && nodewise_ok && network_ok

    out <- data.frame(
      replication = as.integer(result$replication),
      cell_id = design$cell_id,
      truth_structure = design$truth_structure,
      alpha = design$alpha,
      neighbourhood = design$neighbourhood,
      sigma = design$sigma,
      jointly_successful = jointly_successful,
      boundary_relative_ise = NA_real_,
      interior_relative_ise = NA_real_,
      boundary_did = NA_real_,
      stringsAsFactors = FALSE
    )
    if (!jointly_successful) return(out)

    # d_ibc = ISE_network - ISE_nodewise. Positive = relative network loss.
    relative_ise <- get_node_ise(cell_result, "network") -
      get_node_ise(cell_result, "nodewise")
    out$boundary_relative_ise <- mean(relative_ise[boundary_nodes])
    out$interior_relative_ise <- mean(relative_ise[interior_nodes])
    out$boundary_did <- out$boundary_relative_ise -
      out$interior_relative_ise
    out
  })
}

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

checkpoints <- lapply(checkpoint_paths, readRDS)
did_rows <- do.call(
  rbind,
  unlist(lapply(checkpoints, extract_boundary_did), recursive = FALSE)
)
rownames(did_rows) <- NULL

groups <- split(
  did_rows,
  interaction(
    did_rows[, c("cell_id", "alpha", "neighbourhood", "sigma")],
    drop = TRUE, lex.order = TRUE
  )
)
boundary_did_summary <- do.call(rbind, lapply(groups, function(x) {
  keep <- x$jointly_successful %in% TRUE & is.finite(x$boundary_did)
  did <- x$boundary_did[keep]
  mcse <- finite_mcse(did)
  data.frame(
    cell_id = x$cell_id[[1L]],
    truth_structure = "cluster",
    alpha = x$alpha[[1L]],
    neighbourhood = x$neighbourhood[[1L]],
    sigma = x$sigma[[1L]],
    attempted = nrow(x),
    jointly_successful = sum(keep),
    mean_boundary_relative_ise = mean(x$boundary_relative_ise[keep]),
    mean_interior_relative_ise = mean(x$interior_relative_ise[keep]),
    mean_boundary_did = mean(did),
    paired_mcse = mcse,
    lower_two_mcse = mean(did) - 2 * mcse,
    upper_two_mcse = mean(did) + 2 * mcse,
    boundary_loss_rate = mean(did > 0),
    stringsAsFactors = FALSE
  )
}))
rownames(boundary_did_summary) <- NULL
boundary_did_summary <- boundary_did_summary[
  order(
    boundary_did_summary$alpha,
    boundary_did_summary$neighbourhood,
    boundary_did_summary$sigma
  ),
]

analysis_checks <- c(
  saved_main_quality_checks_passed = isTRUE(all(main_result$quality_checks)),
  all_200_main_checkpoints_present = length(checkpoints) == 200L,
  expected_1600_cluster_replication_rows = nrow(did_rows) == 8L * 200L,
  no_duplicate_cell_replication = !anyDuplicated(
    did_rows[, c("cell_id", "replication")]
  ),
  all_boundary_cells_have_200_attempts =
    all(boundary_did_summary$attempted == 200L),
  all_boundary_cells_have_two_or_more_joint_successes =
    all(boundary_did_summary$jointly_successful >= 2L)
)
if (!all(analysis_checks)) {
  stop(
    "Cluster-boundary DiD audit failed: ",
    paste(names(analysis_checks)[!analysis_checks], collapse = ", "), "."
  )
}

utils::write.csv(
  boundary_did_summary,
  file.path(table_directory, "table_4_cluster_boundary_did.csv"),
  row.names = FALSE
)

# Figure 3: boundary difference-in-differences for cluster-structured truth.
# Positive values indicate attenuation of the network estimator's advantage
# at boundary nodes. Error bars show plus or minus two paired MCSEs.

figure_path <- file.path(
  figure_directory,
  "figure_3_cluster_boundary_did.png"
)

plot_data <- boundary_did_summary

boundary_figure_checks <- c(
  expected_eight_rows = nrow(plot_data) == 8L,
  all_cells_have_200_attempts = all(plot_data$attempted == 200L),
  all_pairs_successful = all(plot_data$jointly_successful == 200L),
  expected_noise_levels = identical(
    sort(unique(plot_data$sigma)),
    c(0.2, 0.5)
  ),
  expected_graphs = identical(
    sort(unique(plot_data$neighbourhood)),
    c("queen", "rook")
  ),
  finite_estimates = all(is.finite(plot_data$mean_boundary_did)),
  finite_intervals = all(
    is.finite(plot_data$lower_two_mcse) &
      is.finite(plot_data$upper_two_mcse)
  )
)

if (!isTRUE(all(boundary_figure_checks))) {
  stop(
    "Boundary figure input check failed: ",
    paste(
      names(boundary_figure_checks)[!boundary_figure_checks],
      collapse = ", "
    ),
    "."
  )
}

# Rescale solely for display: one plotted unit equals 0.001 Boundary DiD.
plot_data$boundary_did_scaled <- 1000 * plot_data$mean_boundary_did
plot_data$lower_scaled <- 1000 * plot_data$lower_two_mcse
plot_data$upper_scaled <- 1000 * plot_data$upper_two_mcse

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

# Small horizontal offsets keep Rook and Queen estimates distinguishable.
plot_data$alpha_position <- plot_data$alpha + ifelse(
  plot_data$neighbourhood == "rook",
  -0.025,
  0.025
)

y_range <- range(
  c(plot_data$lower_scaled, plot_data$upper_scaled, 0),
  finite = TRUE
)

y_padding <- max(0.15, 0.08 * diff(y_range))
y_limits <- y_range + c(-y_padding, y_padding)

boundary_figure <- ggplot(
  plot_data,
  aes(
    x = alpha_position,
    y = boundary_did_scaled,
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
  geom_errorbar(
    aes(ymin = lower_scaled, ymax = upper_scaled),
    width = 0.012,
    linewidth = 0.8
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(cols = vars(sigma_label)) +
  scale_x_continuous(
    breaks = c(0.5, 0.9),
    limits = c(0.42, 0.98)
  ) +
  scale_y_continuous(
    limits = y_limits,
    breaks = pretty(y_limits, n = 5)
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
    title = "Boundary attenuation of network smoothing",
    x = expression(alpha),
    y = expression("Boundary DiD (" %*% 10^{-3} * ")"),
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
  filename = figure_path,
  plot = boundary_figure,
  width = 7.3,
  height = 3.7,
  units = "in",
  dpi = 300
)

analysis_result <- list(
  analysis_version = boundary_did_script_version,
  input_path = main_result_path,
  input_completed_at = main_result$completed_at,
  analysis_checks = analysis_checks,
  boundary_figure_checks = boundary_figure_checks,
  boundary_nodes = boundary_nodes,
  interior_nodes = interior_nodes,
  did_rows = did_rows,
  boundary_did_summary = boundary_did_summary,
  table_path = file.path(table_directory, "table_4_cluster_boundary_did.csv"),
  figure_path = figure_path,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(
  analysis_result,
  file.path(analysis_root, "cluster_boundary_did_analysis_result.rds"),
  version = 3
)

cat("\nCluster-boundary DiD analysis checks:\n")
print(analysis_checks)
cat("\nCluster-boundary DiD (network relative to nodewise):\n")
print(boundary_did_summary)
cat("\nOutputs written to:\n", analysis_root, "\n", sep = "")
cat("\nBoundary figure input checks:\n")
print(boundary_figure_checks)