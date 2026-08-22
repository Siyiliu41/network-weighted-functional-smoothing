# Fixed rewired-graph control.
#
# Run from the repository root:
# source("simulation/10_run_rewired_graph_control.R")
# run_rewired_graph_control("dry_run")
# run_rewired_graph_control("main")

rewired_control_script_version <- "2026-08-19-rewired-control-r1"

source("simulation/05_main_core_simulation.R")

rewired_control_cells <- data.frame(
  truth_structure = c("coordinate", "cluster"),
  alpha = c(0.9, 0.9),
  sigma = c(0.5, 0.5),
  stringsAsFactors = FALSE
)
rewired_control_cells$core_rook_cell_id <- sprintf(
  "%s_a%02d_rook_s%02d",
  rewired_control_cells$truth_structure,
  round(100 * rewired_control_cells$alpha),
  round(100 * rewired_control_cells$sigma)
)
rewired_control_cells$cell_id <- sub(
  "_rook_", "_rewired_", rewired_control_cells$core_rook_cell_id, fixed = TRUE
)

make_rewired_control_graph <- function() {
  rook <- make_lattice_graph(n_side = 5L, neighbourhood = "rook")
  rewired <- make_degree_preserving_rewire(
    graph = rook,
    seed = 20260819L,
    n_accepted_swaps = 500L,
    min_replaced_edge_proportion = 0.75
  )
  audit <- attr(rewired, "rewire_audit")
  adjacency <- as.matrix(igraph::as_adjacency_matrix(rewired, sparse = FALSE))
  dimnames(adjacency) <- list(igraph::V(rewired)$name, igraph::V(rewired)$name)
  fingerprint_file <- tempfile(pattern = "rewired_adjacency_", fileext = ".csv")
  on.exit(unlink(fingerprint_file), add = TRUE)
  utils::write.csv(adjacency, fingerprint_file, row.names = TRUE)
  audit$adjacency_md5 <- unname(tools::md5sum(fingerprint_file))
  list(rook = rook, rewired = rewired, adjacency = adjacency, audit = audit)
}

make_rewired_control_scenarios <- function(replication_stream, graph_bundle) {
  validate_rng_stream(replication_stream)
  t_grid <- seq(0, 1, length.out = 50L)
  draws <- generate_core_random_draws(
    replication_stream = replication_stream,
    n_nodes = 25L,
    n_time = length(t_grid),
    n_components = 3L
  )
  scenarios <- lapply(seq_len(nrow(rewired_control_cells)), function(i) {
    cell <- rewired_control_cells[i, , drop = FALSE]
    error_name <- sprintf(
      "error_%s_a%02d", cell$truth_structure, round(10 * cell$alpha)
    )
    truth_object <- generate_core_truth(
      graph = graph_bundle$rook,
      t_grid = t_grid,
      truth_structure = cell$truth_structure,
      alpha = cell$alpha,
      standard_normal_draws = draws$coefficient_draws[[cell$truth_structure]]
    )
    simulated_data <- simulate_observed_curves(
      truth_object = truth_object,
      sigma = cell$sigma,
      standard_normal_errors = draws$error_draws[[error_name]]
    )
    list(
      cell = data.frame(
        cell_id = cell$cell_id,
        truth_structure = cell$truth_structure,
        alpha = cell$alpha,
        neighbourhood = "rewired",
        sigma = cell$sigma,
        stringsAsFactors = FALSE
      ),
      graph = graph_bundle$rewired,
      truth_object = truth_object,
      simulated_data = simulated_data
    )
  })
  names(scenarios) <- rewired_control_cells$cell_id
  scenarios
}

run_rewired_control_replication <- function(
    replication_id, replication_stream, base_seed, config, graph_bundle,
    show_progress = TRUE) {
  scenarios <- make_rewired_control_scenarios(replication_stream, graph_bundle)
  cell_results <- lapply(names(scenarios), function(cell_id) {
    if (show_progress) message(sprintf("[rewired %03d] %s", replication_id, cell_id))
    run_core_cell(
      scenario = scenarios[[cell_id]],
      replication_id = replication_id,
      seed = as.integer(base_seed + 100L + replication_id - 1L),
      k_intercept = config$k_intercept,
      k_deviation = config$k_deviation,
      k_nodewise = config$k_nodewise,
      k_pooled = config$k_pooled,
      keep_predictions = FALSE
    )
  })
  names(cell_results) <- names(scenarios)
  summary_table <- do.call(rbind, lapply(cell_results, function(x) x$summary))
  rownames(summary_table) <- NULL
  if (nrow(summary_table) != 8L ||
      anyDuplicated(summary_table[, c("cell_id", "method")])) {
    stop("A rewired replication must contain two complete four-method cells.")
  }
  list(
    replication = as.integer(replication_id),
    replication_stream = replication_stream,
    summary = summary_table,
    cell_results = cell_results
  )
}

make_rewired_pairing_checks <- function(rewired_summary, main_summary) {
  checks <- lapply(seq_len(nrow(rewired_control_cells)), function(i) {
    cell <- rewired_control_cells[i, , drop = FALSE]
    rewire_rows <- rewired_summary[
      rewired_summary$cell_id == cell$cell_id &
        rewired_summary$method != "network", , drop = FALSE
    ]
    rook_rows <- main_summary[
      main_summary$cell_id == cell$core_rook_cell_id &
        main_summary$replication == rewire_rows$replication[[1L]] &
        main_summary$method != "network", , drop = FALSE
    ]
    rewire_rows <- rewire_rows[order(rewire_rows$method), , drop = FALSE]
    rook_rows <- rook_rows[order(rook_rows$method), , drop = FALSE]
    nrow(rewire_rows) == 3L && nrow(rook_rows) == 3L &&
      identical(rewire_rows$method, rook_rows$method) &&
      identical(rewire_rows$success, rook_rows$success) &&
      identical(rewire_rows$node_averaged_ise, rook_rows$node_averaged_ise) &&
      identical(rewire_rows$warning_count, rook_rows$warning_count)
  })
  c(non_network_results_match_frozen_rook = all(unlist(checks)))
}

make_rewired_paired_summary <- function(rewired_summary, main_summary) {
  output <- lapply(seq_len(nrow(rewired_control_cells)), function(i) {
    cell <- rewired_control_cells[i, , drop = FALSE]
    rewire <- rewired_summary[
      rewired_summary$cell_id == cell$cell_id &
        rewired_summary$method == "network", , drop = FALSE
    ]
    rook <- main_summary[
      main_summary$cell_id == cell$core_rook_cell_id &
        main_summary$method == "network", , drop = FALSE
    ]
    rewire <- rewire[order(rewire$replication), , drop = FALSE]
    rook <- rook[order(rook$replication), , drop = FALSE]
    if (!identical(rewire$replication, rook$replication)) {
      stop("Rewired and correct-rook replication identifiers do not align.")
    }
    joint <- rewire$success %in% TRUE & rook$success %in% TRUE &
      is.finite(rewire$node_averaged_ise) &
      is.finite(rook$node_averaged_ise)
    difference <- rewire$node_averaged_ise[joint] -
      rook$node_averaged_ise[joint]
    data.frame(
      truth_structure = cell$truth_structure,
      alpha = cell$alpha,
      sigma = cell$sigma,
      attempted = nrow(rewire),
      jointly_successful = length(difference),
      mean_rewired_minus_rook_aise = if (length(difference)) mean(difference) else NA_real_,
      paired_mcse = if (length(difference) >= 2L) stats::sd(difference) / sqrt(length(difference)) else NA_real_,
      rewired_to_rook_rmise = if (length(difference)) {
        sum(rewire$node_averaged_ise[joint]) / sum(rook$node_averaged_ise[joint])
      } else NA_real_,
      rewired_worse_win_rate = if (length(difference)) mean(difference > 0) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

run_rewired_graph_control <- function(run_mode = c("dry_run", "main")) {
  run_mode <- match.arg(run_mode)
  main_path <- "simulation/results/core/main/main_core_result.rds"
  if (!file.exists(main_path)) stop("The completed main result is required: ", main_path)
  main_result <- readRDS(main_path)
  if (!isTRUE(all(main_result$quality_checks))) stop("Main quality checks did not pass.")
  
  main_config <- main_result$config
  if (!identical(main_config$n_replications, 200L) ||
      !identical(main_config$stream_positions, c(first = 101L, last = 300L))) {
    stop("The saved main configuration is not the frozen B = 200 design.")
  }
  main_streams <- make_post_pilot_streams(main_config$base_seed, 300L)[101:300]
  graph_bundle <- make_rewired_control_graph()
  config <- list(
    control_version = rewired_control_script_version,
    control_type = "supplementary_fixed_rewired_graph",
    base_seed = main_config$base_seed,
    stream_positions = c(first = 101L, last = 300L),
    n_replications = 200L,
    cells = rewired_control_cells,
    k_intercept = main_config$k_intercept,
    k_deviation = main_config$k_deviation,
    k_nodewise = main_config$k_nodewise,
    k_pooled = main_config$k_pooled,
    maximum_failure_rate = main_config$maximum_failure_rate
  )
  code_files <- c(
    "simulation/10_run_rewired_graph_control.R",
    "simulation/05_main_core_simulation.R", "simulation/R/make_graph.R",
    "simulation/R/rng.R", "simulation/R/generate_truth.R",
    "simulation/R/simulate_data.R", "simulation/R/model_diagnostics.R",
    "simulation/R/fit_baselines.R", "simulation/R/metrics.R",
    "simulation/R/run_replication.R", "simulation/R/run_core_cell.R",
    "simulation/R/run_core_replication.R",
    sort(list.files("netfunsmooth/R", pattern = "[.]R$", full.names = TRUE))
  )
  code_fingerprint <- make_code_fingerprint(code_files)
  root <- "simulation/results/core/negative_controls/rewired"
  
  run_one <- function(replication_id) {
    result <- run_rewired_control_replication(
      replication_id, main_streams[[replication_id]], config$base_seed,
      config, graph_bundle
    )
    list(
      control_version = config$control_version,
      config = config,
      code_fingerprint = code_fingerprint,
      graph_audit = graph_bundle$audit,
      replication_id = as.integer(replication_id),
      replication_stream = main_streams[[replication_id]],
      result = result,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  }
  
  if (run_mode == "dry_run") {
    path <- file.path(root, "preflight", "dry_run_result.rds")
    if (file.exists(path)) {
      old <- readRDS(path)
      if (!identical(old$config, config) ||
          !identical(old$code_fingerprint, code_fingerprint) ||
          !identical(old$graph_audit, graph_bundle$audit)) {
        stop("Existing rewired preflight is incompatible.")
      }
      message("Validated existing rewired preflight: ", path)
      return(invisible(old))
    }
    record <- run_one(1L)
    checks <- c(
      complete_rows = nrow(record$result$summary) == 8L,
      all_methods_successful = all(record$result$summary$success),
      no_warnings = sum(record$result$summary$warning_count) == 0L,
      network_converged = all(record$result$summary$converged[
        record$result$summary$method == "network"
      ]),
      make_rewired_pairing_checks(record$result$summary, main_result$replication_summary)
    )
    record$checks <- checks
    record$graph_adjacency <- graph_bundle$adjacency
    save_rds_atomically(record, path)
    print(graph_bundle$audit)
    print(checks)
    if (!all(checks)) stop("Rewired preflight failed.")
    message("Preflight passed. Now call run_rewired_graph_control(\"main\").")
    return(invisible(record))
  }
  
  preflight_path <- file.path(root, "preflight", "dry_run_result.rds")
  if (!file.exists(preflight_path)) stop("Run and pass the rewired dry run first.")
  preflight <- readRDS(preflight_path)
  if (!identical(preflight$config, config) ||
      !identical(preflight$code_fingerprint, code_fingerprint) ||
      !identical(preflight$graph_audit, graph_bundle$audit) ||
      !isTRUE(all(preflight$checks))) {
    stop("The rewired preflight is missing, incompatible, or failed.")
  }
  
  directory <- file.path(root, "checkpoints")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  checkpoints <- vector("list", config$n_replications)
  for (replication_id in seq_len(config$n_replications)) {
    path <- file.path(directory, sprintf("replication_%03d.rds", replication_id))
    if (file.exists(path)) {
      checkpoint <- readRDS(path)
      if (!identical(checkpoint$config, config) ||
          !identical(checkpoint$code_fingerprint, code_fingerprint) ||
          !identical(checkpoint$graph_audit, graph_bundle$audit) ||
          !identical(checkpoint$replication_stream, main_streams[[replication_id]]) ||
          nrow(checkpoint$result$summary) != 8L) {
        stop("Existing rewired checkpoint is incompatible: ", path)
      }
      checkpoints[[replication_id]] <- checkpoint
      next
    }
    message(sprintf("[rewired %03d/%03d] running.", replication_id, config$n_replications))
    checkpoint <- tryCatch(run_one(replication_id), error = function(e) e)
    if (inherits(checkpoint, "error")) {
      stop("Rewired replication ", replication_id, " failed atomically: ",
           conditionMessage(checkpoint))
    }
    save_rds_atomically(checkpoint, path)
    checkpoints[[replication_id]] <- checkpoint
  }
  
  summary_table <- do.call(rbind, lapply(checkpoints, function(x) x$result$summary))
  rownames(summary_table) <- NULL
  if (nrow(summary_table) != config$n_replications * 8L ||
      anyDuplicated(summary_table[, c("cell_id", "replication", "method")])) {
    stop("Rewired aggregation failed its completeness or duplicate-row check.")
  }
  paired_summary <- make_rewired_paired_summary(summary_table, main_result$replication_summary)
  network_rows <- summary_table[summary_table$method == "network", , drop = FALSE]
  quality_checks <- c(
    complete = nrow(summary_table) == config$n_replications * 8L,
    network_cells_have_200_attempts =
      all(table(network_rows$cell_id) == config$n_replications),
    network_failure_rates_at_most_5_percent =
      all(tapply(network_rows$success %in% TRUE, network_rows$cell_id,
                 function(x) 1 - mean(x)) <= config$maximum_failure_rate),
    paired_cells_have_two_or_more_successes =
      all(paired_summary$jointly_successful >= 2L)
  )
  final_result <- list(
    control_version = config$control_version,
    config = config,
    code_fingerprint = code_fingerprint,
    graph_audit = graph_bundle$audit,
    graph_adjacency = graph_bundle$adjacency,
    quality_checks = quality_checks,
    paired_summary = paired_summary,
    replication_summary = summary_table,
    session_info = utils::sessionInfo(),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  result_path <- file.path(root, "rewired_graph_control_result.rds")
  save_rds_atomically(final_result, result_path)
  print(graph_bundle$audit)
  print(quality_checks)
  message("Rewired control saved to: ", result_path)
  if (!all(quality_checks)) stop("Rewired control completed but failed a quality check.")
  invisible(final_result)
}

message("Rewired control loaded. First run: run_rewired_graph_control(\"dry_run\").")