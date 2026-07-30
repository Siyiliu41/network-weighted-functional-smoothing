# Known-answer simulation
#
# Purpose:
# Verify that network-weighted smoothing improves over nodewise smoothing
# in a setting with a deliberately strong and smooth graph signal.
#
# This script must run from the repository root.

devtools::load_all(
  "netfunsmooth",
  quiet = TRUE
)

source("simulation/R/make_graph.R")
source("simulation/R/generate_truth.R")
source("simulation/R/simulate_data.R")
source("simulation/R/fit_baselines.R")
source("simulation/R/metrics.R")


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

config <- list(
  n_side = 5L,
  neighbourhood = "rook",
  n_grid = 30L,
  signal_scale = 1,
  sigma = 0.4,
  seed = 20260730,
  k_intercept = 10L,
  k_yindex = 10L,
  k_nodewise = 10L,
  k_pooled = 10L
)

t_grid <- seq(
  0,
  1,
  length.out = config$n_grid
)


# -------------------------------------------------------------------------
# Generate graph and true curves
# -------------------------------------------------------------------------

graph <- make_lattice_graph(
  n_side = config$n_side,
  neighbourhood = config$neighbourhood
)

truth_object <- generate_graph_smooth_truth(
  graph = graph,
  t_grid = t_grid,
  signal_scale = config$signal_scale
)


# -------------------------------------------------------------------------
# Known-answer diagnostics for the truth
# -------------------------------------------------------------------------

maximum_centring_error <- max(
  abs(colMeans(truth_object$deviations))
)

curve_distances <- as.matrix(
  stats::dist(truth_object$truth)
)

adjacency_matrix <- as.matrix(
  igraph::as_adjacency_matrix(
    graph,
    sparse = FALSE
  )
)

upper_triangle <- upper.tri(
  adjacency_matrix
)

mean_adjacent_distance <- mean(
  curve_distances[
    adjacency_matrix == 1 & upper_triangle
  ]
)

mean_nonadjacent_distance <- mean(
  curve_distances[
    adjacency_matrix == 0 & upper_triangle
  ]
)

shape_only_curves <- sweep(
  truth_object$truth,
  MARGIN = 1,
  STATS = rowMeans(truth_object$truth),
  FUN = "-"
)

maximum_shape_distance <- max(
  as.matrix(
    stats::dist(shape_only_curves)
  )
)

stopifnot(
  maximum_centring_error < 1e-12,
  mean_adjacent_distance < mean_nonadjacent_distance,
  maximum_shape_distance > 0.1
)


# -------------------------------------------------------------------------
# Add observation noise
# -------------------------------------------------------------------------

simulation_data <- simulate_observed_curves(
  truth_object = truth_object,
  sigma = config$sigma,
  seed = config$seed
)


# -------------------------------------------------------------------------
# Fit network-weighted smoother
# -------------------------------------------------------------------------

network_runtime <- system.time({
  network_fit <- netf_smooth(
    curves = simulation_data$curves,
    graph = simulation_data$graph,
    sandwich = "none",
    bs.int = list(
      bs = "ps",
      k = config$k_intercept,
      m = c(2, 1)
    ),
    bs.yindex = list(
      bs = "ps",
      k = config$k_yindex,
      m = c(2, 1)
    )
  )
})[["elapsed"]]

network_predictions <- predict(
  network_fit
)

network_prediction_matrix <- as.matrix(
  network_predictions,
  arg = simulation_data$t_grid
)


# -------------------------------------------------------------------------
# Fit nodewise smoother
# -------------------------------------------------------------------------

nodewise_runtime <- system.time({
  nodewise_fit <- fit_nodewise_smoothing(
    observed = simulation_data$observed,
    t_grid = simulation_data$t_grid,
    k = config$k_nodewise,
    keep_models = FALSE
  )
})[["elapsed"]]


# -------------------------------------------------------------------------
# Fit pooled smoother
# -------------------------------------------------------------------------

pooled_runtime <- system.time({
  pooled_fit <- fit_pooled_smoothing(
    observed = simulation_data$observed,
    t_grid = simulation_data$t_grid,
    k = config$k_pooled,
    keep_model = FALSE
  )
})[["elapsed"]]


# -------------------------------------------------------------------------
# Calculate reconstruction errors
# -------------------------------------------------------------------------

raw_node_ise <- calculate_node_ise(
  estimate = simulation_data$observed,
  truth = simulation_data$truth,
  t_grid = simulation_data$t_grid
)

pooled_node_ise <- calculate_node_ise(
  estimate = pooled_fit$predictions,
  truth = simulation_data$truth,
  t_grid = simulation_data$t_grid
)

nodewise_node_ise <- calculate_node_ise(
  estimate = nodewise_fit$predictions,
  truth = simulation_data$truth,
  t_grid = simulation_data$t_grid
)

network_node_ise <- calculate_node_ise(
  estimate = network_prediction_matrix,
  truth = simulation_data$truth,
  t_grid = simulation_data$t_grid
)

raw_summary <- summarise_ise(
  raw_node_ise
)

pooled_summary <- summarise_ise(
  pooled_node_ise
)

nodewise_summary <- summarise_ise(
  nodewise_node_ise
)

network_summary <- summarise_ise(
  network_node_ise
)


# -------------------------------------------------------------------------
# Aggregate results
# -------------------------------------------------------------------------

make_result_row <- function(
    method,
    ise_summary,
    elapsed_seconds) {
  
  data.frame(
    method = method,
    node_averaged_ise =
      ise_summary$node_averaged_ise,
    median_node_ise =
      ise_summary$median_node_ise,
    worst_node_ise =
      ise_summary$worst_node_ise,
    elapsed_seconds = elapsed_seconds,
    row.names = NULL
  )
}

performance <- rbind(
  make_result_row(
    "raw",
    raw_summary,
    NA_real_
  ),
  make_result_row(
    "pooled",
    pooled_summary,
    pooled_runtime
  ),
  make_result_row(
    "nodewise",
    nodewise_summary,
    nodewise_runtime
  ),
  make_result_row(
    "network",
    network_summary,
    network_runtime
  )
)

relative_improvement_percent <-
  100 * (
    nodewise_summary$node_averaged_ise -
      network_summary$node_averaged_ise
  ) /
  nodewise_summary$node_averaged_ise


# -------------------------------------------------------------------------
# Known-answer gate
# -------------------------------------------------------------------------

stopifnot(
  network_summary$node_averaged_ise <
    raw_summary$node_averaged_ise,
  
  network_summary$node_averaged_ise <
    nodewise_summary$node_averaged_ise,
  
  pooled_summary$node_averaged_ise >
    nodewise_summary$node_averaged_ise,
  
  inherits(network_fit, "netf_fit"),
  
  isTRUE(network_fit$model$converged),
  
  all(is.finite(network_prediction_matrix))
)


# -------------------------------------------------------------------------
# Save lightweight results
# -------------------------------------------------------------------------

known_answer_result <- list(
  config = config,
  diagnostics = list(
    maximum_centring_error =
      maximum_centring_error,
    mean_adjacent_distance =
      mean_adjacent_distance,
    mean_nonadjacent_distance =
      mean_nonadjacent_distance,
    maximum_shape_distance =
      maximum_shape_distance
  ),
  performance = performance,
  relative_improvement_percent =
    relative_improvement_percent,
  node_ise = list(
    raw = raw_node_ise,
    pooled = pooled_node_ise,
    nodewise = nodewise_node_ise,
    network = network_node_ise
  ),
  smoothing_parameters =
    network_fit$model$sp,
  session_info =
    utils::sessionInfo()
)

saveRDS(
  known_answer_result,
  file = "simulation/results/known_answer_result.rds"
)


# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

print(
  performance,
  row.names = FALSE
)

message(
  sprintf(
    paste0(
      "Known-answer check passed: network smoothing reduced ",
      "node-averaged ISE by %.2f%% relative to nodewise smoothing."
    ),
    relative_improvement_percent
  )
)