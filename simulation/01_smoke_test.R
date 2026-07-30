# Sequential smoke test
#
# Purpose:
# Run a small number of complete replications before implementing the full
# scenario grid or parallel execution.

devtools::load_all(
  "netfunsmooth",
  quiet = TRUE
)

source("simulation/R/make_graph.R")
source("simulation/R/generate_truth.R")
source("simulation/R/simulate_data.R")
source("simulation/R/fit_baselines.R")
source("simulation/R/metrics.R")
source("simulation/R/run_replication.R")


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

smoke_config <- list(
  n_side = 5L,
  neighbourhood = "rook",
  n_grid = 30L,
  signal_scale = 1,
  sigma = 0.4,
  base_seed = 20260730,
  k_intercept = 10L,
  k_yindex = 10L,
  k_nodewise = 10L,
  k_pooled = 10L
)

n_replications <- 3L


# -------------------------------------------------------------------------
# Run replications sequentially
# -------------------------------------------------------------------------

replication_results <- lapply(
  seq_len(n_replications),
  function(replication_id) {
    
    message(
      sprintf(
        "Running smoke-test replication %d of %d.",
        replication_id,
        n_replications
      )
    )
    
    run_replication(
      config = smoke_config,
      replication_id = replication_id,
      keep_predictions = FALSE
    )
  }
)


# -------------------------------------------------------------------------
# Combine method-level summaries
# -------------------------------------------------------------------------

smoke_summary <- do.call(
  rbind,
  lapply(
    replication_results,
    function(result) {
      result$summary
    }
  )
)

rownames(smoke_summary) <- NULL


# -------------------------------------------------------------------------
# Paired network-versus-nodewise comparison
# -------------------------------------------------------------------------

nodewise_rows <- smoke_summary[
  smoke_summary$method == "nodewise",
  c(
    "replication",
    "seed",
    "node_averaged_ise"
  )
]

names(nodewise_rows)[3] <-
  "nodewise_ise"

network_rows <- smoke_summary[
  smoke_summary$method == "network",
  c(
    "replication",
    "seed",
    "node_averaged_ise"
  )
]

names(network_rows)[3] <-
  "network_ise"

paired_comparison <- merge(
  nodewise_rows,
  network_rows,
  by = c("replication", "seed"),
  all = TRUE,
  sort = TRUE
)

paired_comparison$relative_improvement_percent <-
  100 * (
    paired_comparison$nodewise_ise -
      paired_comparison$network_ise
  ) /
  paired_comparison$nodewise_ise


# -------------------------------------------------------------------------
# Aggregate performance
# -------------------------------------------------------------------------

mean_ise_by_method <- stats::aggregate(
  node_averaged_ise ~ method,
  data = smoke_summary,
  FUN = mean,
  na.rm = TRUE
)

mean_runtime_by_method <- stats::aggregate(
  elapsed_seconds ~ method,
  data = smoke_summary[
    !is.na(smoke_summary$elapsed_seconds),
    ,
    drop = FALSE
  ],
  FUN = mean,
  na.rm = TRUE
)

aggregate_performance <- merge(
  mean_ise_by_method,
  mean_runtime_by_method,
  by = "method",
  all.x = TRUE,
  sort = FALSE
)

method_order <- c(
  "raw",
  "pooled",
  "nodewise",
  "network"
)

aggregate_performance <- aggregate_performance[
  match(
    method_order,
    aggregate_performance$method
  ),
  ,
  drop = FALSE
]

rownames(aggregate_performance) <- NULL


# -------------------------------------------------------------------------
# Save results before applying the smoke-test gate
# -------------------------------------------------------------------------

smoke_test_result <- list(
  config = smoke_config,
  n_replications = n_replications,
  summary = smoke_summary,
  paired_comparison = paired_comparison,
  aggregate_performance = aggregate_performance,
  network_smoothing_parameters = lapply(
    replication_results,
    function(result) {
      result$network_smoothing_parameters
    }
  ),
  session_info = utils::sessionInfo()
)

saveRDS(
  smoke_test_result,
  file = "simulation/results/smoke_test_result.rds"
)


# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

print(
  smoke_summary,
  row.names = FALSE
)

cat("\nAverage performance by method:\n")

print(
  aggregate_performance,
  row.names = FALSE
)

cat("\nPaired network-versus-nodewise comparison:\n")

print(
  paired_comparison,
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Smoke-test gate
# -------------------------------------------------------------------------

network_mean_ise <- mean(
  network_rows$network_ise
)

nodewise_mean_ise <- mean(
  nodewise_rows$nodewise_ise
)

mean_relative_improvement <- mean(
  paired_comparison$relative_improvement_percent
)

stopifnot(
  nrow(smoke_summary) ==
    n_replications * 4L,
  
  all(smoke_summary$success),
  
  sum(smoke_summary$warning_count) == 0L,
  
  all(is.na(smoke_summary$error_message)),
  
  all(
    smoke_summary$converged[
      smoke_summary$method == "network"
    ]
  ),
  
  network_mean_ise <
    nodewise_mean_ise,
  
  mean_relative_improvement > 0
)

message(
  sprintf(
    paste0(
      "Smoke test passed: mean network improvement over ",
      "nodewise smoothing was %.2f%% across %d replications."
    ),
    mean_relative_improvement,
    n_replications
  )
)