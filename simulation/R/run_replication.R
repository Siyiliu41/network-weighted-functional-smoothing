# Run one complete simulation replication.
#
# A replication uses the same simulated dataset for all methods, allowing
# paired performance comparisons.

capture_method_run <- function(run) {
  
  if (!is.function(run)) {
    stop("`run` must be a function.")
  }
  
  warning_messages <- character()
  
  start_time <- proc.time()[["elapsed"]]
  
  value <- tryCatch(
    withCallingHandlers(
      run(),
      warning = function(warning_condition) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(warning_condition)
        )
        
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) {
      error_condition
    }
  )
  
  elapsed_seconds <-
    proc.time()[["elapsed"]] - start_time
  
  if (inherits(value, "error")) {
    return(
      list(
        success = FALSE,
        value = NULL,
        warning_messages = unique(warning_messages),
        error_message = conditionMessage(value),
        elapsed_seconds = elapsed_seconds
      )
    )
  }
  
  list(
    success = TRUE,
    value = value,
    warning_messages = unique(warning_messages),
    error_message = NA_character_,
    elapsed_seconds = elapsed_seconds
  )
}


collapse_messages <- function(messages) {
  
  if (length(messages) == 0L) {
    return(NA_character_)
  }
  
  paste(
    unique(messages),
    collapse = " | "
  )
}


make_replication_summary_row <- function(
    replication_id,
    seed,
    method,
    run_result) {
  
  if (run_result$success) {
    ise_summary <- run_result$value$ise_summary
    converged <- run_result$value$converged
  } else {
    ise_summary <- data.frame(
      node_averaged_ise = NA_real_,
      median_node_ise = NA_real_,
      worst_node_ise = NA_real_
    )
    
    converged <- NA
  }
  
  data.frame(
    replication = replication_id,
    seed = seed,
    method = method,
    success = run_result$success,
    converged = converged,
    warning_count =
      length(run_result$warning_messages),
    warning_messages =
      collapse_messages(run_result$warning_messages),
    error_message =
      run_result$error_message,
    node_averaged_ise =
      ise_summary$node_averaged_ise,
    median_node_ise =
      ise_summary$median_node_ise,
    worst_node_ise =
      ise_summary$worst_node_ise,
    elapsed_seconds =
      run_result$elapsed_seconds,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


run_replication <- function(
    config,
    replication_id,
    keep_predictions = FALSE) {
  
  required_config <- c(
    "n_side",
    "neighbourhood",
    "n_grid",
    "signal_scale",
    "sigma",
    "base_seed",
    "k_intercept",
    "k_yindex",
    "k_nodewise",
    "k_pooled"
  )
  
  if (
    !is.list(config) ||
    !all(required_config %in% names(config))
  ) {
    stop(
      "`config` is missing one or more required entries: ",
      paste(required_config, collapse = ", "),
      "."
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
  
  replication_id <- as.integer(replication_id)
  
  seed <-
    as.integer(config$base_seed) +
    replication_id -
    1L
  
  t_grid <- seq(
    0,
    1,
    length.out = config$n_grid
  )
  
  graph <- make_lattice_graph(
    n_side = config$n_side,
    neighbourhood = config$neighbourhood
  )
  
  truth_object <- generate_graph_smooth_truth(
    graph = graph,
    t_grid = t_grid,
    signal_scale = config$signal_scale
  )
  
  simulation_data <- simulate_observed_curves(
    truth_object = truth_object,
    sigma = config$sigma,
    seed = seed
  )
  
  
  # Raw observations -------------------------------------------------------
  
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
  
  # No fitting is performed for the raw observations.
  raw_run$elapsed_seconds <- NA_real_
  
  
  # Pooled smoothing -------------------------------------------------------
  
  pooled_run <- capture_method_run(
    function() {
      
      fit <- fit_pooled_smoothing(
        observed = simulation_data$observed,
        t_grid = simulation_data$t_grid,
        k = config$k_pooled,
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
        converged = NA,
        predictions = if (keep_predictions) {
          fit$predictions
        } else {
          NULL
        }
      )
    }
  )
  
  
  # Nodewise smoothing -----------------------------------------------------
  
  nodewise_run <- capture_method_run(
    function() {
      
      fit <- fit_nodewise_smoothing(
        observed = simulation_data$observed,
        t_grid = simulation_data$t_grid,
        k = config$k_nodewise,
        keep_models = FALSE
      )
      
      node_ise <- calculate_node_ise(
        estimate = fit$predictions,
        truth = simulation_data$truth,
        t_grid = simulation_data$t_grid
      )
      
      list(
        ise_summary = summarise_ise(node_ise),
        node_ise = node_ise,
        converged = NA,
        predictions = if (keep_predictions) {
          fit$predictions
        } else {
          NULL
        }
      )
    }
  )
  
  
  # Network-weighted smoothing --------------------------------------------
  
  network_run <- capture_method_run(
    function() {
      
      fit <- netf_smooth(
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
      
      converged <- isTRUE(
        fit$model$converged
      )
      
      if (!converged) {
        stop("The network-weighted model did not converge.")
      }
      
      prediction_tfd <- predict(fit)
      
      prediction_matrix <- as.matrix(
        prediction_tfd,
        arg = simulation_data$t_grid
      )
      
      if (any(!is.finite(prediction_matrix))) {
        stop(
          "The network-weighted model produced non-finite predictions."
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
        smoothing_parameters = fit$model$sp,
        predictions = if (keep_predictions) {
          prediction_matrix
        } else {
          NULL
        }
      )
    }
  )
  
  
  # Combine method-level results ------------------------------------------
  
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
  
  node_ise <- lapply(
    method_runs,
    function(result) {
      if (result$success) {
        result$value$node_ise
      } else {
        NULL
      }
    }
  )
  
  predictions <- if (keep_predictions) {
    lapply(
      method_runs,
      function(result) {
        if (result$success) {
          result$value$predictions
        } else {
          NULL
        }
      }
    )
  } else {
    NULL
  }
  
  network_smoothing_parameters <-
    if (network_run$success) {
      network_run$value$smoothing_parameters
    } else {
      NULL
    }
  
  list(
    replication = replication_id,
    seed = seed,
    config = config,
    summary = summary_table,
    node_ise = node_ise,
    predictions = predictions,
    network_smoothing_parameters =
      network_smoothing_parameters
  )
}