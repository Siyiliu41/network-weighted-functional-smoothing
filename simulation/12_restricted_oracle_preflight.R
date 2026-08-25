# One-dataset implementation check for the restricted oracle.
#
# This preflight is not the 50-ID oracle pilot and does not save an RDS file.
# Run from the repository root with:
#   source("simulation/12_restricted_oracle_preflight.R")
#   run_restricted_oracle_preflight()

devtools::load_all("netfunsmooth", quiet = TRUE)

oracle_preflight_sources <- c(
  "simulation/R/make_graph.R",
  "simulation/R/rng.R",
  "simulation/R/generate_truth.R",
  "simulation/R/simulate_data.R",
  "simulation/R/model_diagnostics.R",
  "simulation/R/metrics.R",
  "simulation/R/run_replication.R",
  "simulation/R/core_design.R",
  "simulation/R/restricted_oracle.R"
)

invisible(lapply(oracle_preflight_sources, source))

run_restricted_oracle_preflight <- function() {
  # Independent implementation-check seed; not part of the pilot/main IDs.
  preflight_stream <- make_replication_streams(
    base_seed = 20260824L,
    n_replications = 1L
  )[[1L]]
  
  generated <- generate_core_scenarios(
    replication_stream = preflight_stream
  )
  
  scenario <- generated$scenarios[["coordinate_a90_rook_s20"]]
  
  reml_run <- capture_method_run(
    function() {
      fit_oracle_network(
        simulation_data = scenario$simulated_data,
        graph = scenario$graph
      )
    }
  )
  
  if (!isTRUE(reml_run$success)) {
    stop("Restricted-oracle preflight failed: initial REML fit was unsuccessful.")
  }
  
  oracle <- run_restricted_oracle(
    simulation_data = scenario$simulated_data,
    graph = scenario$graph,
    reml_run = reml_run
  )
  
  if (!identical(
    names(oracle$smoothing_parameter_layout),
    c("intercept", "graph", "time")
  )) {
    stop("Restricted-oracle preflight failed: parameter layout is not intercept/graph/time.")
  }
  
  if (!isTRUE(oracle$success)) {
    stop("Restricted-oracle preflight failed: at least one required candidate fit failed.")
  }
  
  if (!is.finite(oracle$zero_offset_max_abs_difference) ||
      oracle$zero_offset_max_abs_difference > 1e-8) {
    stop("Restricted-oracle preflight failed: the (0, 0) candidate does not reproduce REML within 1e-8.")
  }
  
  candidate_count <- nrow(oracle$candidate_table)
  
  if (candidate_count < 25L || candidate_count > 33L) {
    stop("Restricted-oracle preflight failed: unexpected coarse-to-fine candidate count.")
  }
  
  message("Restricted-oracle preflight passed.")
  
  print(list(
    scenario_id = scenario$cell$cell_id,
    smoothing_parameter_layout = oracle$smoothing_parameter_layout,
    coarse_candidates = sum(
      oracle$candidate_table$stage == "coarse"
    ),
    total_candidates = candidate_count,
    zero_offset_max_abs_difference =
      oracle$zero_offset_max_abs_difference
  ))
  
  invisible(oracle)
}