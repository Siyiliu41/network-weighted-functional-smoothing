# Restricted oracle helpers for the two interaction smoothing parameters.
#
# The functional-intercept parameter remains fixed at its REML value. The
# oracle is a diagnostic: candidates are selected with the known true curves.

oracle_candidate_grid <- function(graph_offsets, time_offsets) {
  grid <- expand.grid(
    graph_offset = graph_offsets,
    time_offset = time_offsets,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  grid[order(grid$graph_offset, grid$time_offset), , drop = FALSE]
}


network_smoothing_parameter_layout <- function(network_fit) {
  model <- if (inherits(network_fit, "netf_fit")) {
    network_fit$model
  } else {
    network_fit
  }
  
  if (!inherits(model, "gam") || length(model$sp) != 3L ||
      is.null(model$smooth)) {
    stop("Restricted oracle requires a fitted GAM with exactly three smoothing parameters.")
  }
  
  labels <- vapply(model$smooth, function(smooth) smooth$label, character(1L))
  intercept_term <- which(grepl("^s\\(yindex\\.vec", labels))
  interaction_term <- which(grepl("^ti\\(node, ?yindex\\.vec", labels))
  
  if (length(intercept_term) != 1L || length(interaction_term) != 1L) {
    stop("Could not identify the functional intercept and node-by-time interaction.")
  }
  
  intercept_smooth <- model$smooth[[intercept_term]]
  interaction_smooth <- model$smooth[[interaction_term]]
  
  intercept_sp <- seq.int(intercept_smooth$first.sp, intercept_smooth$last.sp)
  interaction_sp <- seq.int(
    interaction_smooth$first.sp,
    interaction_smooth$last.sp
  )
  
  if (length(intercept_sp) != 1L || length(interaction_sp) != 2L ||
      is.null(interaction_smooth$margin)) {
    stop("The fitted model does not have the required one-plus-two penalty layout.")
  }
  
  margin_classes <- vapply(
    interaction_smooth$margin,
    function(margin) paste(class(margin), collapse = " "),
    character(1L)
  )
  
  graph_margin <- which(grepl("mrf", margin_classes))
  time_margin <- which(grepl("ps", margin_classes))
  
  if (length(graph_margin) != 1L || length(time_margin) != 1L) {
    stop("Could not identify the graph and time interaction penalties.")
  }
  
  layout <- c(
    intercept = intercept_sp,
    graph = interaction_sp[[graph_margin]],
    time = interaction_sp[[time_margin]]
  )
  
  if (!identical(sort(unname(layout)), seq_along(model$sp))) {
    stop("The extracted smoothing-parameter layout is inconsistent with model$sp.")
  }
  
  layout
}


fit_oracle_network <- function(
    simulation_data,
    graph,
    k_intercept = 10L,
    k_deviation = 15L,
    sp = NULL) {
  
  fit_arguments <- list(
    curves = simulation_data$curves,
    graph = graph,
    sandwich = "none",
    bs.int = list(bs = "ps", k = k_intercept, m = c(2, 1)),
    bs.yindex = list(bs = "ps", k = k_deviation, m = c(2, 1))
  )
  
  if (!is.null(sp)) {
    if (!is.numeric(sp) || length(sp) != 3L ||
        any(!is.finite(sp)) || any(sp <= 0)) {
      stop("Fixed smoothing parameters must be three positive finite values.")
    }
    
    fit_arguments$sp <- sp
  }
  
  fit <- do.call(netf_smooth, fit_arguments)
  
  diagnostics <- extract_gam_diagnostics(
    model = fit$model,
    fit_id = "network"
  )
  
  if (!isTRUE(diagnostics$fit$converged[[1L]])) {
    stop("The network-weighted model did not converge.")
  }
  
  predictions <- as.matrix(predict(fit), arg = simulation_data$t_grid)
  
  if (any(!is.finite(predictions))) {
    stop("The network-weighted model produced non-finite predictions.")
  }
  
  node_ise <- calculate_node_ise(
    estimate = predictions,
    truth = simulation_data$truth,
    t_grid = simulation_data$t_grid
  )
  
  list(
    fit = fit,
    predictions = predictions,
    node_ise = node_ise,
    ise_summary = summarise_ise(node_ise),
    smoothing_parameters = fit$model$sp,
    diagnostics = diagnostics
  )
}


make_oracle_candidate_row <- function(
    candidate_id,
    stage,
    graph_offset,
    time_offset,
    fixed_sp,
    run_result) {
  
  success <- isTRUE(run_result$success)
  estimate <- if (success) run_result$value$ise_summary else NULL
  
  data.frame(
    candidate_id = as.integer(candidate_id),
    stage = stage,
    graph_offset = graph_offset,
    time_offset = time_offset,
    fixed_sp_intercept = fixed_sp[[1L]],
    fixed_sp_graph = fixed_sp[[2L]],
    fixed_sp_time = fixed_sp[[3L]],
    success = success,
    warning_count = length(run_result$warning_messages),
    error_message = run_result$error_message,
    node_averaged_ise = if (success) estimate$node_averaged_ise else NA_real_,
    median_node_ise = if (success) estimate$median_node_ise else NA_real_,
    worst_node_ise = if (success) estimate$worst_node_ise else NA_real_,
    elapsed_seconds = run_result$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}


evaluate_oracle_candidates <- function(
    simulation_data,
    graph,
    reml_sp,
    sp_layout,
    candidates,
    stage,
    candidate_id_start = 1L) {
  
  runs <- vector("list", nrow(candidates))
  rows <- vector("list", nrow(candidates))
  
  for (candidate_index in seq_len(nrow(candidates))) {
    fixed_sp <- reml_sp
    
    fixed_sp[[sp_layout[["graph"]]]] <-
      reml_sp[[sp_layout[["graph"]]]] *
      exp(candidates$graph_offset[[candidate_index]])
    
    fixed_sp[[sp_layout[["time"]]]] <-
      reml_sp[[sp_layout[["time"]]]] *
      exp(candidates$time_offset[[candidate_index]])
    
    run_result <- capture_method_run(
      function() {
        fit_oracle_network(
          simulation_data = simulation_data,
          graph = graph,
          sp = fixed_sp
        )
      }
    )
    
    candidate_id <- candidate_id_start + candidate_index - 1L
    
    runs[[candidate_index]] <- run_result
    
    rows[[candidate_index]] <- make_oracle_candidate_row(
      candidate_id = candidate_id,
      stage = stage,
      graph_offset = candidates$graph_offset[[candidate_index]],
      time_offset = candidates$time_offset[[candidate_index]],
      fixed_sp = fixed_sp[unname(sp_layout)],
      run_result = run_result
    )
  }
  
  list(
    runs = runs,
    table = do.call(rbind, rows)
  )
}


select_oracle_candidate <- function(candidate_table) {
  successful <- candidate_table[
    candidate_table$success &
      is.finite(candidate_table$node_averaged_ise),
    ,
    drop = FALSE
  ]
  
  if (nrow(successful) == 0L) {
    return(NA_integer_)
  }
  
  minimum_ise <- min(successful$node_averaged_ise)
  tolerance <- 1e-12 * max(1, minimum_ise)
  
  tied <- successful[
    abs(successful$node_averaged_ise - minimum_ise) <= tolerance,
    ,
    drop = FALSE
  ]
  
  tied$distance_from_reml <- sqrt(
    tied$graph_offset^2 + tied$time_offset^2
  )
  
  tied <- tied[
    order(
      tied$distance_from_reml,
      tied$graph_offset,
      tied$time_offset
    ),
    ,
    drop = FALSE
  ]
  
  tied$candidate_id[[1L]]
}


run_restricted_oracle <- function(
    simulation_data,
    graph,
    reml_run = NULL,
    coarse_graph_offsets = c(-4, -2, 0, 2, 4),
    coarse_time_offsets = c(-4, -2, 0, 2, 4),
    refinement_offsets = c(-1, 0, 1),
    permitted_offset_range = c(-8, 8)) {
  
  if (is.null(reml_run)) {
    reml_run <- capture_method_run(
      function() {
        fit_oracle_network(
          simulation_data = simulation_data,
          graph = graph
        )
      }
    )
  }
  
  if (!isTRUE(reml_run$success)) {
    return(list(
      success = FALSE,
      completion_status = "reml-failed",
      reml_run = reml_run,
      candidate_table = data.frame()
    ))
  }
  
  reml_sp <- reml_run$value$smoothing_parameters
  
  if (!is.numeric(reml_sp) || length(reml_sp) != 3L ||
      any(!is.finite(reml_sp)) || any(reml_sp <= 0)) {
    stop("REML did not return three positive finite smoothing parameters.")
  }
  
  sp_layout <- network_smoothing_parameter_layout(reml_run$value$fit)
  
  coarse <- evaluate_oracle_candidates(
    simulation_data = simulation_data,
    graph = graph,
    reml_sp = reml_sp,
    sp_layout = sp_layout,
    candidates = oracle_candidate_grid(
      coarse_graph_offsets,
      coarse_time_offsets
    ),
    stage = "coarse"
  )
  
  coarse_winner_id <- select_oracle_candidate(coarse$table)
  
  if (is.na(coarse_winner_id)) {
    return(list(
      success = FALSE,
      completion_status = "no-successful-coarse-candidate",
      reml_run = reml_run,
      candidate_table = coarse$table,
      smoothing_parameter_layout = sp_layout
    ))
  }
  
  coarse_winner <- coarse$table[
    coarse$table$candidate_id == coarse_winner_id,
    ,
    drop = FALSE
  ]
  
  refinement_grid <- oracle_candidate_grid(
    coarse_winner$graph_offset + refinement_offsets,
    coarse_winner$time_offset + refinement_offsets
  )
  
  refinement_grid <- refinement_grid[
    refinement_grid$graph_offset >= permitted_offset_range[[1L]] &
      refinement_grid$graph_offset <= permitted_offset_range[[2L]] &
      refinement_grid$time_offset >= permitted_offset_range[[1L]] &
      refinement_grid$time_offset <= permitted_offset_range[[2L]],
    ,
    drop = FALSE
  ]
  
  refinement_key <- paste(
    refinement_grid$graph_offset,
    refinement_grid$time_offset
  )
  
  coarse_key <- paste(
    coarse$table$graph_offset,
    coarse$table$time_offset
  )
  
  refinement_grid <- refinement_grid[
    !refinement_key %in% coarse_key,
    ,
    drop = FALSE
  ]
  
  refinement <- if (nrow(refinement_grid) == 0L) {
    list(
      runs = list(),
      table = coarse$table[FALSE, , drop = FALSE]
    )
  } else {
    evaluate_oracle_candidates(
      simulation_data = simulation_data,
      graph = graph,
      reml_sp = reml_sp,
      sp_layout = sp_layout,
      candidates = refinement_grid,
      stage = "refinement",
      candidate_id_start = nrow(coarse$table) + 1L
    )
  }
  
  candidate_table <- rbind(coarse$table, refinement$table)
  all_runs <- c(coarse$runs, refinement$runs)
  
  selected_candidate_id <- select_oracle_candidate(candidate_table)
  selected_run <- all_runs[[selected_candidate_id]]
  
  zero_row <- candidate_table[
    candidate_table$graph_offset == 0 &
      candidate_table$time_offset == 0,
    ,
    drop = FALSE
  ]
  
  zero_run <- all_runs[[zero_row$candidate_id[[1L]]]]
  
  zero_difference <- if (isTRUE(zero_run$success)) {
    max(abs(
      zero_run$value$predictions -
        reml_run$value$predictions
    ))
  } else {
    NA_real_
  }
  
  all_required_candidates_successful <- all(candidate_table$success)
  
  list(
    success = all_required_candidates_successful,
    completion_status = if (all_required_candidates_successful) {
      "complete"
    } else {
      "incomplete-candidate-search"
    },
    reml_run = reml_run,
    candidate_table = candidate_table,
    coarse_winner_id = coarse_winner_id,
    selected_candidate_id = selected_candidate_id,
    oracle_run = selected_run,
    reml_oracle_gap = if (all_required_candidates_successful) {
      reml_run$value$ise_summary$node_averaged_ise -
        selected_run$value$ise_summary$node_averaged_ise
    } else {
      NA_real_
    },
    zero_offset_max_abs_difference = zero_difference,
    smoothing_parameter_layout = sp_layout
  )
}