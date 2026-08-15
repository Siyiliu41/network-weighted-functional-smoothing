# Frozen main core simulation driver
#
# Usage from the repository root:
#   source("simulation/05_main_core_simulation.R")
#   run_main_core_simulation("dry_run")
#   run_main_core_simulation("main")
#
# The dry run fits one complete 16-cell replication using stream position 101.
# It is an implementation preflight only and is not included in the main
# summaries. The main run then starts again at replication 1 with the same
# stream, as a deterministic re-computation of the preflight dataset.

main_core_script_version <- "2026-08-10-main-driver-r2-interface-compatible"

devtools::load_all("netfunsmooth", quiet = TRUE)

source("simulation/R/make_graph.R")
source("simulation/R/rng.R")
source("simulation/R/generate_truth.R")
source("simulation/R/simulate_data.R")
source("simulation/R/model_diagnostics.R")
source("simulation/R/fit_baselines.R")
source("simulation/R/metrics.R")
source("simulation/R/run_replication.R")
source("simulation/R/run_core_cell.R")
source("simulation/R/run_core_replication.R")


# Compatibility wrapper: `replication_stream` is the source of randomness.
# Older and newer `run_core_replication()` implementations differ only in
# whether they also accept bookkeeping arguments `base_seed` and `phase`.
run_main_replication <- function(
    replication_id,
    replication_stream,
    base_seed,
    config,
    show_progress = TRUE) {

  arguments <- list(
    replication_id = replication_id,
    replication_stream = replication_stream,
    k_intercept = config$k_intercept,
    k_deviation = config$k_deviation,
    k_nodewise = config$k_nodewise,
    k_pooled = config$k_pooled,
    keep_predictions = FALSE,
    show_progress = show_progress
  )

  supported_arguments <- names(formals(run_core_replication))

  arguments$seed <- as.integer(base_seed + 100L + replication_id - 1L)

  if ("base_seed" %in% supported_arguments) {
    arguments$base_seed <- base_seed
  }

  if ("phase" %in% supported_arguments) {
    arguments$phase <- "main"
  }

  do.call(run_core_replication, arguments)
}


normalised_file_md5 <- function(path) {
  if (!file.exists(path)) stop("Cannot fingerprint missing file: ", path)

  temporary_file <- tempfile(pattern = "normalised_", fileext = ".txt")
  on.exit(unlink(temporary_file), add = TRUE)

  normalised_text <- paste0(
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
    "\n"
  )
  writeBin(charToRaw(enc2utf8(normalised_text)), temporary_file)
  unname(tools::md5sum(temporary_file))
}


make_code_fingerprint <- function(paths) {
  paths <- sort(unique(paths))
  hashes <- vapply(paths, normalised_file_md5, character(1L))
  names(hashes) <- paths
  hashes
}


# This is deliberately local to the post-pilot drivers. The original pilot
# consumed stream positions 1--50 directly, and the formal sanity check used
# 51--100. The generic pilot-era ledger must not be used to relabel that
# history.
make_post_pilot_streams <- function(base_seed, n_streams) {
  if (
    !is.numeric(base_seed) || length(base_seed) != 1L ||
    is.na(base_seed) || !is.finite(base_seed) ||
    base_seed != as.integer(base_seed) ||
    !is.numeric(n_streams) || length(n_streams) != 1L ||
    is.na(n_streams) || n_streams < 1L ||
    n_streams != as.integer(n_streams)
  ) {
    stop("`base_seed` and `n_streams` must be finite positive integers.")
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(base_seed))
  next_stream <- get(".Random.seed", envir = .GlobalEnv)
  streams <- vector("list", as.integer(n_streams))

  for (i in seq_along(streams)) {
    streams[[i]] <- next_stream
    next_stream <- parallel::nextRNGStream(next_stream)
  }
  streams
}


save_rds_atomically <- function(object, path, overwrite = FALSE) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing file: ", path)
  }

  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"), tmpdir = directory,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, version = 3)

  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish file atomically: ", path)
  }
  invisible(path)
}


validate_main_checkpoint <- function(
    checkpoint, replication_id, replication_stream, config, code_fingerprint) {
  required_entries <- c(
    "main_version", "config", "code_fingerprint", "replication_id",
    "replication_stream", "result"
  )
  if (!is.list(checkpoint) || !all(required_entries %in% names(checkpoint))) {
    stop("Main checkpoint ", replication_id, " has an incompatible structure.")
  }

  checks <- c(
    version = identical(checkpoint$main_version, config$main_version),
    config = identical(checkpoint$config, config),
    fingerprint = identical(checkpoint$code_fingerprint, code_fingerprint),
    replication_id = identical(checkpoint$replication_id, as.integer(replication_id)),
    replication_stream = identical(checkpoint$replication_stream, replication_stream)
  )
  if (!all(checks)) {
    stop(
      "Main checkpoint ", replication_id, " is incompatible: ",
      paste(names(checks)[!checks], collapse = ", "), "."
    )
  }

  result <- checkpoint$result
  if (!is.list(result) || !is.data.frame(result$summary) ||
      nrow(result$summary) != 16L * 4L) {
    stop("Main checkpoint ", replication_id, " has an incomplete result.")
  }
  invisible(checkpoint)
}


make_main_method_summary <- function(summary_table) {
  groups <- split(summary_table, list(summary_table$cell_id, summary_table$method),
                  drop = TRUE)
  result <- lapply(groups, function(x) {
    successful <- x$success %in% TRUE & is.finite(x$node_averaged_ise)
    data.frame(
      cell_id = x$cell_id[[1L]],
      truth_structure = x$truth_structure[[1L]],
      alpha = x$alpha[[1L]],
      neighbourhood = x$neighbourhood[[1L]],
      sigma = x$sigma[[1L]],
      method = x$method[[1L]],
      attempted = nrow(x),
      successful = sum(successful),
      failure_rate = 1 - mean(successful),
      total_warnings = sum(x$warning_count),
      mean_aise = if (any(successful)) mean(x$node_averaged_ise[successful]) else NA_real_,
      median_aise = if (any(successful)) stats::median(x$node_averaged_ise[successful]) else NA_real_,
      mean_elapsed_seconds = mean(x$elapsed_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}


make_paired_cell_summary <- function(summary_table) {
  cells <- sort(unique(summary_table$cell_id))
  result <- lapply(cells, function(cell_id) {
    cell_rows <- summary_table[summary_table$cell_id == cell_id, , drop = FALSE]
    nodewise <- cell_rows[cell_rows$method == "nodewise", , drop = FALSE]
    network <- cell_rows[cell_rows$method == "network", , drop = FALSE]
    nodewise <- nodewise[order(nodewise$replication), , drop = FALSE]
    network <- network[order(network$replication), , drop = FALSE]

    if (!identical(nodewise$replication, network$replication)) {
      stop("Primary-method replication IDs do not align in cell ", cell_id)
    }
    joint_success <-
      nodewise$success %in% TRUE & network$success %in% TRUE &
      is.finite(nodewise$node_averaged_ise) &
      is.finite(network$node_averaged_ise)
    delta <- nodewise$node_averaged_ise[joint_success] -
      network$node_averaged_ise[joint_success]
    n_joint <- length(delta)

    data.frame(
      cell_id = cell_id,
      truth_structure = nodewise$truth_structure[[1L]],
      alpha = nodewise$alpha[[1L]],
      neighbourhood = nodewise$neighbourhood[[1L]],
      sigma = nodewise$sigma[[1L]],
      attempted = nrow(nodewise),
      jointly_successful = n_joint,
      mean_paired_delta = if (n_joint) mean(delta) else NA_real_,
      paired_mcse = if (n_joint >= 2L) stats::sd(delta) / sqrt(n_joint) else NA_real_,
      win_rate = if (n_joint) mean(delta > 0) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}


run_main_core_simulation <- function(run_mode = c("dry_run", "main")) {
  run_mode <- match.arg(run_mode)

  pilot_directory <- "simulation/results/core/pilot/checkpoints"
  formal_result_path <-
    "simulation/results/core/formal_known_answer/formal_known_answer_result.rds"
  if (!file.exists(formal_result_path)) {
    stop("The passed formal known-answer result is required before main simulation.")
  }
  formal_result <- readRDS(formal_result_path)
  if (!isTRUE(formal_result$formal_passed)) {
    stop("The formal known-answer check did not pass. Do not start main simulation.")
  }

  pilot_paths <- file.path(pilot_directory, sprintf("replication_%03d.rds", 1:50))
  if (!all(file.exists(pilot_paths))) {
    stop("All 50 pilot checkpoints are required to validate the main seed allocation.")
  }
  pilot_checkpoints <- lapply(pilot_paths, readRDS)
  base_seed <- pilot_checkpoints[[1L]]$config$base_seed
  if (is.null(base_seed) || any(!vapply(pilot_checkpoints, function(x) {
    identical(x$config$base_seed, base_seed)
  }, logical(1L)))) {
    stop("Pilot checkpoints do not contain one consistent frozen base seed.")
  }

  n_main <- 200L
  all_streams <- make_post_pilot_streams(base_seed, 100L + n_main)
  pilot_matches <- vapply(seq_len(50L), function(i) {
    identical(pilot_checkpoints[[i]]$replication_stream, all_streams[[i]])
  }, logical(1L))
  if (!all(pilot_matches)) {
    stop("Pilot seed-stream verification failed at replication(s): ",
         paste(which(!pilot_matches), collapse = ", "), ".")
  }

  formal_streams <- formal_result$seed_ledger$streams$sanity
  if (!is.list(formal_streams) || length(formal_streams) != 50L ||
      !all(vapply(seq_len(50L), function(i) {
        identical(formal_streams[[i]], all_streams[[50L + i]])
      }, logical(1L)))) {
    stop("Formal sanity seed-stream verification failed. Do not start main simulation.")
  }
  main_streams <- all_streams[101:(100L + n_main)]

  config <- list(
    main_version = main_core_script_version,
    base_seed = as.integer(base_seed),
    phase = "main",
    stream_positions = c(first = 101L, last = 300L),
    n_replications = n_main,
    n_side = 5L,
    n_grid = 50L,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    k_pooled = 10L,
    primary_methods = c("network", "nodewise"),
    maximum_failure_rate = 0.05
  )

  code_files <- c(
    "simulation/05_main_core_simulation.R", "simulation/R/make_graph.R",
    "simulation/R/rng.R", "simulation/R/generate_truth.R",
    "simulation/R/simulate_data.R", "simulation/R/model_diagnostics.R",
    "simulation/R/fit_baselines.R", "simulation/R/metrics.R",
    "simulation/R/run_replication.R", "simulation/R/run_core_cell.R",
    "simulation/R/run_core_replication.R",
    sort(list.files("netfunsmooth/R", pattern = "[.]R$", full.names = TRUE))
  )
  code_fingerprint <- make_code_fingerprint(code_files)
  root <- "simulation/results/core/main"

  if (run_mode == "dry_run") {
    dry_path <- file.path(root, "preflight", "dry_run_result.rds")
    if (file.exists(dry_path)) {
      previous <- readRDS(dry_path)
      checks <- c(
        config = identical(previous$config, config),
        fingerprint = identical(previous$code_fingerprint, code_fingerprint),
        stream = identical(previous$replication_stream, main_streams[[1L]])
      )
      if (!all(checks)) stop("Existing dry-run result is incompatible: ",
                             paste(names(checks)[!checks], collapse = ", "), ".")
      message("Validated existing dry run: ", dry_path)
      return(invisible(previous))
    }

    message("Main preflight: running one complete 16-cell replication (stream 101).")
    dry_result <- run_main_replication(
      replication_id = 1L,
      replication_stream = main_streams[[1L]],
      base_seed = base_seed,
      config = config,
      show_progress = TRUE
    )
    dry_summary <- dry_result$summary
    dry_checks <- c(
      complete_rows = nrow(dry_summary) == 16L * 4L,
      all_methods_successful = all(dry_summary$success),
      no_warnings = sum(dry_summary$warning_count) == 0L,
      network_converged = all(dry_summary$converged[dry_summary$method == "network"])
    )
    dry_record <- list(config = config, code_fingerprint = code_fingerprint,
                       replication_stream = main_streams[[1L]], result = dry_result,
                       checks = dry_checks, completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    save_rds_atomically(dry_record, dry_path)
    print(dry_checks)
    if (!all(dry_checks)) {
      stop("Main preflight failed. Do not start B = 200; inspect the saved record.")
    }
    message("Main preflight passed. Now call run_main_core_simulation(\"main\").")
    return(invisible(dry_record))
  }

  dry_path <- file.path(root, "preflight", "dry_run_result.rds")
  if (!file.exists(dry_path)) {
    stop("Run and pass run_main_core_simulation(\"dry_run\") before the main run.")
  }
  dry_record <- readRDS(dry_path)
  if (!identical(dry_record$config, config) ||
      !identical(dry_record$code_fingerprint, code_fingerprint) ||
      !isTRUE(all(dry_record$checks))) {
    stop("The dry-run record is missing, incompatible, or failed. Do not start main.")
  }

  checkpoint_directory <- file.path(root, "checkpoints")
  dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
  checkpoints <- vector("list", n_main)
  overall_start <- proc.time()[["elapsed"]]

  for (replication_id in seq_len(n_main)) {
    checkpoint_path <- file.path(checkpoint_directory,
                                 sprintf("replication_%03d.rds", replication_id))
    replication_stream <- main_streams[[replication_id]]
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      validate_main_checkpoint(checkpoint, replication_id, replication_stream,
                               config, code_fingerprint)
      checkpoints[[replication_id]] <- checkpoint
      message(sprintf("[Main %03d/%03d] validated existing checkpoint.",
                      replication_id, n_main))
      next
    }

    message(sprintf("[Main %03d/%03d] running.", replication_id, n_main))
    replication_start <- proc.time()[["elapsed"]]
    result <- tryCatch(
      run_main_replication(
        replication_id = replication_id,
        replication_stream = replication_stream,
        base_seed = base_seed,
        config = config,
        show_progress = TRUE
      ),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      stop("Replication ", replication_id, " did not finish atomically: ",
           conditionMessage(result), ". No checkpoint was written; resume will recompute it.")
    }

    checkpoint <- list(
      main_version = config$main_version, config = config,
      code_fingerprint = code_fingerprint, replication_id = as.integer(replication_id),
      replication_stream = replication_stream, result = result,
      elapsed_seconds = proc.time()[["elapsed"]] - replication_start,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    save_rds_atomically(checkpoint, checkpoint_path)
    checkpoints[[replication_id]] <- checkpoint
  }

  summary_table <- do.call(rbind, lapply(checkpoints, function(x) x$result$summary))
  rownames(summary_table) <- NULL
  expected_rows <- n_main * 16L * 4L
  if (nrow(summary_table) != expected_rows ||
      anyDuplicated(summary_table[, c("cell_id", "replication", "method")])) {
    stop("Main aggregation failed its completeness or duplicate-row check.")
  }

  method_summary <- make_main_method_summary(summary_table)
  paired_summary <- make_paired_cell_summary(summary_table)
  primary_summary <- method_summary[method_summary$method %in% config$primary_methods, , drop = FALSE]
  quality_checks <- c(
    complete = nrow(summary_table) == expected_rows,
    all_primary_cells_have_200_attempts = all(primary_summary$attempted == n_main),
    primary_failure_rates_at_most_5_percent = all(primary_summary$failure_rate <= config$maximum_failure_rate),
    all_paired_cells_have_two_or_more_successes = all(paired_summary$jointly_successful >= 2L)
  )

  main_result <- list(
    main_version = config$main_version, config = config,
    seed_allocation = list(pilot = 1:50, sanity = 51:100, main = 101:300),
    code_fingerprint = code_fingerprint, quality_checks = quality_checks,
    method_summary = method_summary, paired_summary = paired_summary,
    replication_summary = summary_table,
    checkpoint_paths = file.path(checkpoint_directory, sprintf("replication_%03d.rds", 1:n_main)),
    total_elapsed_seconds = proc.time()[["elapsed"]] - overall_start,
    session_info = utils::sessionInfo(),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  result_path <- file.path(root, "main_core_result.rds")
  save_rds_atomically(main_result, result_path)

  cat("\nMain core-simulation quality checks:\n")
  print(quality_checks)
  cat("\nResult saved to:\n", result_path, "\n", sep = "")
  if (!all(quality_checks)) {
    stop("Main simulation completed but failed a prespecified quality check. Results were saved; pause interpretation.")
  }
  message("Main core simulation completed successfully.")
  invisible(main_result)
}


message(
  "Main core driver loaded. First run: run_main_core_simulation(\"dry_run\")."
)