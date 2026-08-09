# Fit all four methods for one core simulation cell.
#
# This file does not modify the original run_replication() workflow.

run_core_cell <- function(
    scenario,
    replication_id,
    seed,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    k_pooled = 10L,
    keep_predictions = FALSE) {

  required_scenario_entries <- c(
    "cell",
    "graph",
    "truth_object",
    "simulated_data"
  )

  if (
    !is.list(scenario) ||
    !all(required_scenario_entries %in% names(scenario))
  ) {
    stop(
      "`scenario` must be one element returned by ",
      "`generate_core_scenarios()`."
    )
  }

  if (
    !is.numeric(replication_id) ||
    length(replication_id) != 1L ||
    is.na(replication_id) ||
    replication_id < 1 ||
    replication_id %% 1 != 0
  ) {
    stop("`replication_id` must be a positive integer.")
  }

  if (
    !is.numeric(seed) ||
    length(seed) != 1L ||
    is.na(seed) ||
    !is.finite(seed)
  ) {
    stop("`seed` must be one finite number.")
  }

  basis_dimensions <- c(
    k_intercept = k_intercept,
    k_deviation = k_deviation,
    k_nodewise = k_nodewise,
    k_pooled = k_pooled
  )

  if (
    anyNA(basis_dimensions) ||
    any(basis_dimensions < 3) ||
    any(basis_dimensions %% 1 != 0)
  ) {
    stop("All basis dimensions must be integers of at least 3.")
  }

  replication_id <- as.integer(replication_id)
  seed <- as.integer(seed)

  simulation_data <- scenario$simulated_data
  cell <- scenario$cell

  if (!is.data.frame(cell) || nrow(cell) != 1L) {
    stop("`scenario$cell` must be a one-row data frame.")
  }


  # Raw observations -----------------------------------------------------

  raw_run <- capture_method_run(
    function() {
      node_ise <- calculate_node_ise(
        estimate = simulation_data$observed,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )

      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = NA,
        predictions = if (keep_predictions) {
          simulation_data$observed
        } else {
          NULL
        }
      )
    }
  )

  # Raw observations require no fitting.
  raw_run$elapsed_seconds <- NA_real_


  # Pooled smoothing -----------------------------------------------------

  pooled_run <- capture_method_run(
    function() {
      fit <- fit_pooled_smoothing(
        observed = simulation_data$observed,
        t_grid = simulation_data$t_grid,
        k = k_pooled,
        keep_model = FALSE
      )

      node_ise <- calculate_node_ise(
        estimate = fit$predictions,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )

      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = isTRUE(
          fit$diagnostics$fit$converged[1L]
        ),
        diagnostics = fit$diagnostics,
        predictions = if (keep_predictions) {
          fit$predictions
        } else {
          NULL
        }
      )
    }
  )


  # Nodewise smoothing ---------------------------------------------------

  nodewise_run <- capture_method_run(
    function() {
      fit <- fit_nodewise_smoothing(
        observed = simulation_data$observed,
        t_grid = simulation_data$t_grid,
        k = k_nodewise,
        keep_models = FALSE
      )

      node_ise <- calculate_node_ise(
        estimate = fit$predictions,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )

      nodewise_convergence <-
        fit$diagnostics$fit$converged

      all_nodewise_models_converged <-
        length(nodewise_convergence) > 0L &&
        all(nodewise_convergence %in% TRUE)

      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = all_nodewise_models_converged,
        diagnostics = fit$diagnostics,
        predictions = if (keep_predictions) {
          fit$predictions
        } else {
          NULL
        }
      )
    }
  )


  # Network-weighted smoothing ------------------------------------------

  network_run <- capture_method_run(
    function() {
      fit <- netf_smooth(
        curves = simulation_data$curves,

        # Use the graph assigned to this estimator cell,
        # not the graph stored in the truth object.
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

      converged <- isTRUE(
        diagnostics$fit$converged[1L]
      )

      if (!converged) {
        stop(
          "The network-weighted model did not converge."
        )
      }

      prediction_tfd <- predict(fit)

      prediction_matrix <- as.matrix(
        prediction_tfd,
        arg = simulation_data$t_grid
      )

      if (any(!is.finite(prediction_matrix))) {
        stop(
          "The network-weighted model produced ",
          "non-finite predictions."
        )
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

        # Retained for compatibility with existing pilot summaries.
        smoothing_parameters = fit$model$sp,

        predictions = if (keep_predictions) {
          prediction_matrix
        } else {
          NULL
        }
      )
    }
  )


  # Combine method-level results ----------------------------------------

  method_runs <- list(
    raw = raw_run,
    pooled = pooled_run,
    nodewise = nodewise_run,
    network = network_run
  )

  summary_rows <- lapply(
    names(method_runs),
    function(method_name) {
      make_replication_summary_row(
        replication_id = replication_id,
        seed = seed,
        method = method_name,
        run_result = method_runs[[method_name]]
      )
    }
  )

  summary_table <- do.call(
    rbind,
    summary_rows
  )

  rownames(summary_table) <- NULL

  summary_table$cell_id <- cell$cell_id
  summary_table$truth_structure <-
    cell$truth_structure
  summary_table$alpha <- cell$alpha
  summary_table$neighbourhood <-
    cell$neighbourhood
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


  # Collect node-level ISE results.

  node_ise <- lapply(
    method_runs,
    function(run_result) {
      if (run_result$success) {
        run_result$value$node_ise
      } else {
        NULL
      }
    }
  )


  # Retain predictions only when requested.

  predictions <- if (keep_predictions) {
    lapply(
      method_runs,
      function(run_result) {
        if (run_result$success) {
          run_result$value$predictions
        } else {
          NULL
        }
      }
    )
  } else {
    NULL
  }


  # Retain the old network smoothing-parameter object for
  # compatibility with the existing pilot-summary code.

  network_smoothing_parameters <- if (
    network_run$success
  ) {
    network_run$value$smoothing_parameters
  } else {
    NULL
  }


  # Collect diagnostics from all fitted smoothing methods.
  # Raw has no diagnostics because it does not fit a GAM.

  model_diagnostics <- lapply(
    c("pooled", "nodewise", "network"),
    function(method_name) {
      run_result <- method_runs[[method_name]]

      if (run_result$success) {
        run_result$value$diagnostics
      } else {
        NULL
      }
    }
  )

  names(model_diagnostics) <- c(
    "pooled",
    "nodewise",
    "network"
  )


  # Return all results from this cell and replication.

  list(
    cell = cell,
    replication = replication_id,
    seed = seed,
    basis_dimensions = as.list(
      basis_dimensions
    ),
    summary = summary_table,
    node_ise = node_ise,
    predictions = predictions,
    model_diagnostics = model_diagnostics,
    network_smoothing_parameters =
      network_smoothing_parameters
  )
}