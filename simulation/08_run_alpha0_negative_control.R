# Prespecified uninformative-graph negative control (alpha = 0)
#
# Save this file as simulation/08_run_alpha0_negative_control.R and run from
# the repository root:
#   source("simulation/08_run_alpha0_negative_control.R")
#   run_alpha0_negative_control("dry_run")
#   run_alpha0_negative_control("main")
#
# This driver is deliberately separate from the 16-cell core simulation.  It
# reuses the 200 main replication streams, but writes only below
# simulation/results/core/negative_controls/alpha0/.  The two estimator-graph
# cells use exactly the same truth and observation-error draw within each
# replication.  The nodewise/raw/pooled results are consequently also checked
# for exact pairing across the two graph labels.

alpha0_control_script_version <- "2026-08-11-alpha0-negative-control-r1"

source("simulation/05_main_core_simulation.R")


make_alpha0_control_scenarios <- function(replication_stream) {
  validate_rng_stream(replication_stream)

  graphs <- list(
    rook = make_lattice_graph(n_side = 5L, neighbourhood = "rook"),
    queen = make_lattice_graph(n_side = 5L, neighbourhood = "queen")
  )
  t_grid <- seq(0, 1, length.out = 50L)
  random_draws <- generate_core_random_draws(
    replication_stream = replication_stream,
    n_nodes = 25L,
    n_time = length(t_grid),
    n_components = 3L
  )

  # alpha = 0 removes the graph-structured coefficient component.  The
  # coordinate component stream and coordinate_a50 error stream are reused
  # deterministically for every main replication identifier.
  truth_object <- generate_core_truth(
    graph = graphs$rook,
    t_grid = t_grid,
    truth_structure = "coordinate",
    alpha = 0,
    standard_normal_draws = random_draws$coefficient_draws$coordinate
  )
  simulated_data <- simulate_observed_curves(
    truth_object = truth_object,
    sigma = 0.5,
    standard_normal_errors = random_draws$error_draws$coordinate_a50
  )

  scenarios <- lapply(names(graphs), function(neighbourhood) {
    list(
      cell = data.frame(
        cell_id = sprintf("coordinate_a00_%s_s50", neighbourhood),
        truth_structure = "coordinate",
        alpha = 0,
        neighbourhood = neighbourhood,
        sigma = 0.5,
        stringsAsFactors = FALSE
      ),
      graph = graphs[[neighbourhood]],
      truth_object = truth_object,
      simulated_data = simulated_data
    )
  })
  names(scenarios) <- names(graphs)

  list(
    graphs = graphs,
    random_draws = random_draws,
    truth_object = truth_object,
    simulated_data = simulated_data,
    scenarios = scenarios
  )
}


run_alpha0_control_replication <- function(
    replication_id, replication_stream, base_seed, config, show_progress = TRUE) {
  scenarios <- make_alpha0_control_scenarios(replication_stream)
  cell_results <- lapply(names(scenarios$scenarios), function(neighbourhood) {
    if (show_progress) {
      message(sprintf("[alpha0 %03d] %s estimator graph", replication_id, neighbourhood))
    }
    run_core_cell(
      scenario = scenarios$scenarios[[neighbourhood]],
      replication_id = replication_id,
      seed = base_seed,
      k_intercept = config$k_intercept,
      k_deviation = config$k_deviation,
      k_nodewise = config$k_nodewise,
      k_pooled = config$k_pooled,
      keep_predictions = FALSE
    )
  })
  names(cell_results) <- names(scenarios$scenarios)

  summary_table <- do.call(rbind, lapply(cell_results, `[[`, "summary"))
  rownames(summary_table) <- NULL
  if (nrow(summary_table) != 8L || anyDuplicated(summary_table[, c("cell_id", "method")])) {
    stop("The alpha = 0 replication did not contain exactly two complete cells.")
  }

  list(
    replication = as.integer(replication_id),
    replication_stream = replication_stream,
    summary = summary_table,
    cell_results = cell_results
  )
}


validate_alpha0_checkpoint <- function(
    checkpoint, replication_id, replication_stream, config, code_fingerprint) {
  needed <- c("control_version", "config", "code_fingerprint", "replication_id",
              "replication_stream", "result")
  if (!is.list(checkpoint) || !all(needed %in% names(checkpoint))) {
    stop("alpha = 0 checkpoint ", replication_id, " has an incompatible structure.")
  }
  checks <- c(
    version = identical(checkpoint$control_version, config$control_version),
    config = identical(checkpoint$config, config),
    fingerprint = identical(checkpoint$code_fingerprint, code_fingerprint),
    replication_id = identical(checkpoint$replication_id, as.integer(replication_id)),
    stream = identical(checkpoint$replication_stream, replication_stream),
    complete_result = is.data.frame(checkpoint$result$summary) && nrow(checkpoint$result$summary) == 8L
  )
  if (!all(checks)) {
    stop("alpha = 0 checkpoint ", replication_id, " is incompatible: ",
         paste(names(checks)[!checks], collapse = ", "), ".")
  }
  invisible(checkpoint)
}


make_alpha0_pairing_checks <- function(summary_table) {
  non_network <- summary_table[summary_table$method != "network", , drop = FALSE]
  groups <- split(non_network, list(non_network$replication, non_network$method), drop = TRUE)
  exact_matches <- vapply(groups, function(x) {
    x <- x[order(x$neighbourhood), , drop = FALSE]
    nrow(x) == 2L &&
      identical(x$success[[1L]], x$success[[2L]]) &&
      identical(x$node_averaged_ise[[1L]], x$node_averaged_ise[[2L]]) &&
      identical(x$warning_count[[1L]], x$warning_count[[2L]])
  }, logical(1L))
  c(non_network_results_are_exactly_paired = all(exact_matches))
}


run_alpha0_negative_control <- function(run_mode = c("dry_run", "main")) {
  run_mode <- match.arg(run_mode)
  main_path <- "simulation/results/core/main/main_core_result.rds"
  if (!file.exists(main_path)) stop("The completed main result is required: ", main_path)
  main_result <- readRDS(main_path)
  if (!isTRUE(all(main_result$quality_checks))) {
    stop("The main simulation has not passed its quality checks. Pause the control run.")
  }

  config_main <- main_result$config
  if (!identical(config_main$n_replications, 200L) ||
      !identical(config_main$stream_positions, c(first = 101L, last = 300L))) {
    stop("The saved main configuration is not the frozen B = 200 design.")
  }
  main_streams <- make_post_pilot_streams(config_main$base_seed, 300L)[101:300]

  config <- list(
    control_version = alpha0_control_script_version,
    control_type = "prespecified_uninformative_graph",
    base_seed = config_main$base_seed,
    stream_positions = c(first = 101L, last = 300L),
    n_replications = 200L,
    truth_structure = "coordinate",
    alpha = 0,
    sigma = 0.5,
    estimator_graphs = c("rook", "queen"),
    k_intercept = config_main$k_intercept,
    k_deviation = config_main$k_deviation,
    k_nodewise = config_main$k_nodewise,
    k_pooled = config_main$k_pooled,
    maximum_failure_rate = config_main$maximum_failure_rate
  )
  code_files <- c(
    "simulation/08_run_alpha0_negative_control.R",
    "simulation/05_main_core_simulation.R", "simulation/R/make_graph.R",
    "simulation/R/rng.R", "simulation/R/generate_truth.R",
    "simulation/R/simulate_data.R", "simulation/R/model_diagnostics.R",
    "simulation/R/fit_baselines.R", "simulation/R/metrics.R",
    "simulation/R/run_replication.R", "simulation/R/run_core_cell.R",
    "simulation/R/run_core_replication.R",
    sort(list.files("netfunsmooth/R", pattern = "[.]R$", full.names = TRUE))
  )
  code_fingerprint <- make_code_fingerprint(code_files)
  root <- "simulation/results/core/negative_controls/alpha0"

  if (run_mode == "dry_run") {
    path <- file.path(root, "preflight", "dry_run_result.rds")
    if (file.exists(path)) {
      old <- readRDS(path)
      checks <- c(config = identical(old$config, config),
                  fingerprint = identical(old$code_fingerprint, code_fingerprint),
                  stream = identical(old$replication_stream, main_streams[[1L]]))
      if (!all(checks)) stop("Existing alpha = 0 preflight is incompatible: ",
                             paste(names(checks)[!checks], collapse = ", "), ".")
      message("Validated existing alpha = 0 preflight: ", path)
      return(invisible(old))
    }
    result <- run_alpha0_control_replication(1L, main_streams[[1L]], config$base_seed, config)
    checks <- c(
      complete_rows = nrow(result$summary) == 8L,
      all_methods_successful = all(result$summary$success),
      no_warnings = sum(result$summary$warning_count) == 0L,
      network_converged = all(result$summary$converged[result$summary$method == "network"]),
      pairing = make_alpha0_pairing_checks(result$summary)
    )
    record <- list(config = config, code_fingerprint = code_fingerprint,
                   replication_stream = main_streams[[1L]], result = result, checks = checks,
                   completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    save_rds_atomically(record, path)
    print(checks)
    if (!all(checks)) stop("alpha = 0 preflight failed. Inspect the saved record before B = 200.")
    message("alpha = 0 preflight passed. Now call run_alpha0_negative_control(\"main\").")
    return(invisible(record))
  }

  preflight_path <- file.path(root, "preflight", "dry_run_result.rds")
  if (!file.exists(preflight_path)) stop("Run and pass the alpha = 0 dry run first.")
  preflight <- readRDS(preflight_path)
  if (!identical(preflight$config, config) ||
      !identical(preflight$code_fingerprint, code_fingerprint) ||
      !isTRUE(all(preflight$checks))) {
    stop("The alpha = 0 preflight is missing, incompatible, or failed.")
  }

  checkpoint_directory <- file.path(root, "checkpoints")
  dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
  checkpoints <- vector("list", config$n_replications)
  for (replication_id in seq_len(config$n_replications)) {
    path <- file.path(checkpoint_directory, sprintf("replication_%03d.rds", replication_id))
    if (file.exists(path)) {
      checkpoint <- readRDS(path)
      validate_alpha0_checkpoint(checkpoint, replication_id, main_streams[[replication_id]], config, code_fingerprint)
      checkpoints[[replication_id]] <- checkpoint
      message(sprintf("[alpha0 %03d/%03d] validated existing checkpoint.", replication_id, config$n_replications))
      next
    }
    message(sprintf("[alpha0 %03d/%03d] running.", replication_id, config$n_replications))
    start <- proc.time()[["elapsed"]]
    result <- tryCatch(
      run_alpha0_control_replication(replication_id, main_streams[[replication_id]], config$base_seed, config),
      error = function(e) e
    )
    if (inherits(result, "error")) stop("alpha = 0 replication ", replication_id,
                                        " failed atomically: ", conditionMessage(result))
    checkpoint <- list(
      control_version = config$control_version, config = config,
      code_fingerprint = code_fingerprint, replication_id = as.integer(replication_id),
      replication_stream = main_streams[[replication_id]], result = result,
      elapsed_seconds = proc.time()[["elapsed"]] - start,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    save_rds_atomically(checkpoint, path)
    checkpoints[[replication_id]] <- checkpoint
  }

  summary_table <- do.call(rbind, lapply(checkpoints, function(x) x$result$summary))
  rownames(summary_table) <- NULL
  expected_rows <- config$n_replications * 2L * 4L
  method_summary <- make_main_method_summary(summary_table)
  paired_summary <- make_paired_cell_summary(summary_table)
  primary <- method_summary[method_summary$method %in% c("network", "nodewise"), , drop = FALSE]
  quality_checks <- c(
    expected_rows = nrow(summary_table) == expected_rows,
    no_duplicate_cell_replication_method = !anyDuplicated(summary_table[, c("cell_id", "replication", "method")]),
    all_primary_cells_have_200_attempts = all(primary$attempted == config$n_replications),
    primary_failure_rates_at_most_5_percent = all(primary$failure_rate <= config$maximum_failure_rate),
    all_network_nodewise_pairs_have_two_or_more_successes = all(paired_summary$jointly_successful >= 2L),
    pairing = make_alpha0_pairing_checks(summary_table)
  )
  result <- list(
    control_version = config$control_version, config = config,
    protocol_status = "prespecified_cell_run_after_the_main_core_block",
    code_fingerprint = code_fingerprint, quality_checks = quality_checks,
    method_summary = method_summary, paired_summary = paired_summary,
    replication_summary = summary_table,
    checkpoint_paths = file.path(checkpoint_directory, sprintf("replication_%03d.rds", 1:config$n_replications)),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  path <- file.path(root, "alpha0_control_result.rds")
  save_rds_atomically(result, path)
  cat("\nalpha = 0 negative-control quality checks:\n")
  print(quality_checks)
  cat("\nResult saved to:\n", path, "\n", sep = "")
  if (!all(quality_checks)) stop("The alpha = 0 control completed but failed a quality check.")
  invisible(result)
}


message("alpha = 0 control driver loaded. First run: run_alpha0_negative_control(\"dry_run\").")