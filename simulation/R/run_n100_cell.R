# Run one of the two fixed N = 100 scalability cells.
#
# Unlike run_core_cell(), this runner deliberately omits raw and pooled
# methods. It retains the same method-level failure and diagnostic structure.

run_n100_cell <- function(
    scenario,
    replication_id,
    seed,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    keep_predictions = FALSE) {
  
  required_scenario_entries <- c(
    "cell",
    "graph",
    "truth_object",
    "simulated_data"
  )
  
  if (!is.list(scenario) ||
      !all(required_scenario_entries %in% names(scenario))) {
    stop("`scenario` must be a fixed N = 100 scenario.")
  }
  
  if (!is.numeric(replication_id) ||
      length(replication_id) != 1L ||
      is.na(replication_id) ||
      replication_id < 1L ||
      replication_id %% 1L != 0L) {
    stop("`replication_id` must be a positive integer.")
  }
  
  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed)) {
    stop("`seed` must be one finite number.")
  }
  
  cell <- scenario$cell
  simulation_data <- scenario$simulated_data
  
  required_cell_columns <- c(
    "cell_id",
    "truth_structure",
    "alpha",
    "neighbourhood",
    "sigma"
  )
  
  if (!is.data.frame(cell) ||
      nrow(cell) != 1L ||
      !all(required_cell_columns %in% names(cell))) {
    stop("`scenario$cell` must describe one N = 100 cell.")
  }
  
  expected_cell_ids <- c(
    "n100_coordinate_a90_rook_s50",
    "n100_coordinate_a00_rook_s50"
  )
  
  if (!cell$cell_id %in% expected_cell_ids ||
      !identical(cell$truth_structure, "coordinate") ||
      !identical(cell$neighbourhood, "rook") ||
      !isTRUE(all.equal(cell$sigma, 0.5)) ||
      !isTRUE(cell$alpha %in% c(0, 0.9)) ||
      igraph::vcount(scenario$graph) != 100L ||
      igraph::ecount(scenario$graph) != 180L ||
      !is.matrix(simulation_data$observed) ||
      !identical(dim(simulation_data$observed), c(100L, 50L))) {
    stop("`scenario` differs from the prespecified N = 100 design.")
  }
  
  basis_dimensions <- c(
    k_intercept = k_intercept,
    k_deviation = k_deviation,
    k_nodewise = k_nodewise
  )
  
  if (anyNA(basis_dimensions) ||
      any(basis_dimensions < 3L) ||
      any(basis_dimensions %% 1L != 0L)) {
    stop("All basis dimensions must be integers of at least 3.")
  }
  
  replication_id <- as.integer(replication_id)
  seed <- as.integer(seed)
  
  nodewise_run <- capture_method_run(
    function() {
      fit <- fit_nodewise_smoothing(
        observed = simulation_data$observed,
        t_grid = simulation_data$t_grid,
        k = k_nodewise,
        keep_models = FALSE
      )
      
      converged <- fit$diagnostics$fit$converged
      all_converged <-
        length(converged) == nrow(simulation_data$observed) &&
        all(converged %in% TRUE)
      
      if (!all_converged) {
        stop("At least one nodewise smoothing model did not converge.")
      }
      
      if (any(!is.finite(fit$predictions))) {
        stop("The nodewise smoothing models produced non-finite predictions.")
      }
      
      node_ise <- calculate_node_ise(
        estimate = fit$predictions,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )
      
      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = all_converged,
        diagnostics = fit$diagnostics,
        predictions = if (keep_predictions) fit$predictions else NULL
      )
    }
  )
  
  network_run <- capture_method_run(
    function() {
      fit <- netf_smooth(
        curves = simulation_data$curves,
        graph = scenario$graph,
        sandwich = "none",
        bs.int = list(
          bs = "ps",
          k = k_intercept,
          m = c(2, 1)
        ),
        bs.yindex = list(
          bs = "ps",
          k = k_deviation,
          m = c(2, 1)
        )
      )
      
      diagnostics <- extract_gam_diagnostics(
        model = fit$model,
        fit_id = "network"
      )
      
      converged <- isTRUE(diagnostics$fit$converged[1L])
      if (!converged) {
        stop("The network-weighted model did not converge.")
      }
      
      prediction_matrix <- as.matrix(
        predict(fit),
        arg = simulation_data$t_grid
      )
      
      if (any(!is.finite(prediction_matrix))) {
        stop("The network-weighted model produced non-finite predictions.")
      }
      
      node_ise <- calculate_node_ise(
        estimate = prediction_matrix,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )
      
      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = converged,
        diagnostics = diagnostics,
        smoothing_parameters = fit$model$sp,
        predictions = if (keep_predictions) prediction_matrix else NULL
      )
    }
  )
  
  method_runs <- list(
    network = network_run,
    nodewise = nodewise_run
  )
  
  summary_table <- do.call(
    rbind,
    lapply(names(method_runs), function(method_name) {
      make_replication_summary_row(
        replication_id = replication_id,
        seed = seed,
        method = method_name,
        run_result = method_runs[[method_name]]
      )
    })
  )
  
  summary_table$cell_id <- cell$cell_id
  summary_table$truth_structure <- cell$truth_structure
  summary_table$alpha <- cell$alpha
  summary_table$neighbourhood <- cell$neighbourhood
  summary_table$sigma <- cell$sigma
  
  summary_table <- summary_table[
    ,
    c(
      "cell_id",
      "truth_structure",
      "alpha",
      "neighbourhood",
      "sigma",
      "replication",
      "seed",
      "method",
      "success",
      "converged",
      "warning_count",
      "warning_messages",
      "error_message",
      "node_averaged_ise",
      "median_node_ise",
      "worst_node_ise",
      "elapsed_seconds"
    )
  ]
  
  list(
    cell = cell,
    replication = replication_id,
    seed = seed,
    basis_dimensions = as.list(basis_dimensions),
    summary = summary_table,
    network_diagnostics = if (network_run$success) {
      network_run$value$diagnostics
    } else {
      NULL
    },
    nodewise_diagnostics = if (nodewise_run$success) {
      nodewise_run$value$diagnostics
    } else {
      NULL
    },
    predictions = if (keep_predictions) {
      lapply(method_runs, function(x) {
        if (x$success) x$value$predictions else NULL
      })
    } else {
      NULL
    }
  )
}