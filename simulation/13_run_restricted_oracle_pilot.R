# Restricted-oracle pilot: four fixed N = 25 cells, B = 50.
#
# Run from the repository root:
#   source("simulation/13_run_restricted_oracle_pilot.R")
#   run_restricted_oracle_pilot()
#
# The run is rep-major. Each checkpoint contains all four cells for one
# completed replication. Re-running resumes from the next replication.

restricted_oracle_pilot_version <- "2026-08-24-oracle-pilot-b50-r1"

devtools::load_all("netfunsmooth", quiet = TRUE)

oracle_source_files <- c(
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

invisible(lapply(oracle_source_files, source))


oracle_save_rds_atomically <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  
  saveRDS(object, temporary_path, version = 3)
  
  if (file.exists(path)) {
    file.remove(path)
  }
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not save checkpoint: ", path)
  }
  
  invisible(path)
}


oracle_code_fingerprint <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  
  if (length(missing_paths) > 0L) {
    stop(
      "Cannot fingerprint missing file(s): ",
      paste(missing_paths, collapse = ", ")
    )
  }
  
  hashes <- unname(tools::md5sum(paths))
  names(hashes) <- paths
  hashes
}


make_oracle_pilot_design <- function() {
  data.frame(
    cell_id = c(
      "coordinate_a90_rook_s20",
      "coordinate_a90_rook_s50",
      "cluster_a90_rook_s50",
      "coordinate_a00_rook_s50"
    ),
    truth_structure = c(
      "coordinate",
      "coordinate",
      "cluster",
      "coordinate"
    ),
    alpha = c(0.9, 0.9, 0.9, 0),
    neighbourhood = c("rook", "rook", "rook", "rook"),
    sigma = c(0.2, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
}


make_oracle_pilot_scenarios <- function(replication_stream) {
  generated <- generate_core_scenarios(
    replication_stream = replication_stream
  )
  
  required_core_cells <- c(
    "coordinate_a90_rook_s20",
    "coordinate_a90_rook_s50",
    "cluster_a90_rook_s50"
  )
  
  if (!all(required_core_cells %in% names(generated$scenarios))) {
    stop("The required oracle cells are missing from generate_core_scenarios().")
  }
  
  t_grid <- seq(0, 1, length.out = 50L)
  
  # alpha = 0 is outside the 16-cell core design. It uses the same coordinate
  # random-draw streams as the other coordinate settings, as in the existing
  # alpha = 0 negative-control script.
  alpha0_truth <- generate_core_truth(
    graph = generated$graphs$rook,
    t_grid = t_grid,
    truth_structure = "coordinate",
    alpha = 0,
    standard_normal_draws =
      generated$random_draws$coefficient_draws$coordinate
  )
  
  alpha0_data <- simulate_observed_curves(
    truth_object = alpha0_truth,
    sigma = 0.5,
    standard_normal_errors =
      generated$random_draws$error_draws$coordinate_a50
  )
  
  alpha0_scenario <- list(
    cell = data.frame(
      cell_id = "coordinate_a00_rook_s50",
      truth_structure = "coordinate",
      alpha = 0,
      neighbourhood = "rook",
      sigma = 0.5,
      stringsAsFactors = FALSE
    ),
    graph = generated$graphs$rook,
    truth_object = alpha0_truth,
    simulated_data = alpha0_data
  )
  
  scenarios <- c(
    generated$scenarios[required_core_cells],
    list(coordinate_a00_rook_s50 = alpha0_scenario)
  )
  
  expected_ids <- make_oracle_pilot_design()$cell_id
  
  if (!identical(names(scenarios), expected_ids)) {
    stop("The oracle-pilot scenario order differs from the frozen design.")
  }
  
  scenarios
}


run_oracle_pilot_cell <- function(scenario, replication_id) {
  cell_id <- scenario$cell$cell_id[[1L]]
  
  reml_run <- capture_method_run(
    function() {
      fit_oracle_network(
        simulation_data = scenario$simulated_data,
        graph = scenario$graph
      )
    }
  )
  
  oracle <- run_restricted_oracle(
    simulation_data = scenario$simulated_data,
    graph = scenario$graph,
    reml_run = reml_run
  )
  
  reml_ise <- if (isTRUE(reml_run$success)) {
    reml_run$value$ise_summary$node_averaged_ise[[1L]]
  } else {
    NA_real_
  }
  
  selected_candidate_id <- if (!is.null(oracle$selected_candidate_id)) {
    oracle$selected_candidate_id
  } else {
    NA_integer_
  }
  
  candidate_table <- oracle$candidate_table
  
  selected_candidate <- candidate_table[
    candidate_table$candidate_id == selected_candidate_id,
    ,
    drop = FALSE
  ]
  
  oracle_ise <- if (
    isTRUE(oracle$success) &&
    nrow(selected_candidate) == 1L
  ) {
    selected_candidate$node_averaged_ise[[1L]]
  } else {
    NA_real_
  }
  
  reml_summary <- data.frame(
    cell_id = cell_id,
    replication = as.integer(replication_id),
    success = isTRUE(reml_run$success),
    node_averaged_ise = reml_ise,
    elapsed_seconds = reml_run$elapsed_seconds,
    error_message = reml_run$error_message,
    stringsAsFactors = FALSE
  )
  
  oracle_summary <- data.frame(
    cell_id = cell_id,
    replication = as.integer(replication_id),
    reml_success = isTRUE(reml_run$success),
    oracle_success = isTRUE(oracle$success),
    completion_status = if (!is.null(oracle$completion_status)) {
      oracle$completion_status
    } else {
      NA_character_
    },
    selected_candidate_id = selected_candidate_id,
    coarse_winner_id = if (!is.null(oracle$coarse_winner_id)) {
      oracle$coarse_winner_id
    } else {
      NA_integer_
    },
    selected_graph_offset = if (nrow(selected_candidate) == 1L) {
      selected_candidate$graph_offset[[1L]]
    } else {
      NA_real_
    },
    selected_time_offset = if (nrow(selected_candidate) == 1L) {
      selected_candidate$time_offset[[1L]]
    } else {
      NA_real_
    },
    reml_node_averaged_ise = reml_ise,
    oracle_node_averaged_ise = oracle_ise,
    reml_oracle_ise_difference = reml_ise - oracle_ise,
    relative_reml_oracle_improvement = if (
      is.finite(reml_ise) &&
      reml_ise > 0 &&
      is.finite(oracle_ise)
    ) {
      (reml_ise - oracle_ise) / reml_ise
    } else {
      NA_real_
    },
    zero_offset_max_abs_difference =
      if (!is.null(oracle$zero_offset_max_abs_difference)) {
        oracle$zero_offset_max_abs_difference
      } else {
        NA_real_
      },
    stringsAsFactors = FALSE
  )
  
  if (nrow(candidate_table) > 0L) {
    candidate_table$cell_id <- cell_id
    candidate_table$replication <- as.integer(replication_id)
  }
  
  list(
    reml_summary = reml_summary,
    oracle_summary = oracle_summary,
    candidate_summary = candidate_table
  )
}


run_oracle_pilot_replication <- function(
    replication_id,
    replication_stream,
    show_progress = TRUE) {
  scenarios <- make_oracle_pilot_scenarios(replication_stream)
  
  cell_results <- lapply(names(scenarios), function(cell_id) {
    if (show_progress) {
      message(
        sprintf(
          "[Oracle %03d] %s",
          replication_id,
          cell_id
        )
      )
    }
    
    run_oracle_pilot_cell(
      scenario = scenarios[[cell_id]],
      replication_id = replication_id
    )
  })
  names(cell_results) <- names(scenarios)
  
  list(
    replication = as.integer(replication_id),
    replication_stream = replication_stream,
    reml_summary = do.call(
      rbind,
      lapply(cell_results, `[[`, "reml_summary")
    ),
    oracle_summary = do.call(
      rbind,
      lapply(cell_results, `[[`, "oracle_summary")
    ),
    candidate_summary = {
      tables <- Filter(
        function(x) nrow(x) > 0L,
        lapply(cell_results, `[[`, "candidate_summary")
      )
      
      if (length(tables) > 0L) do.call(rbind, tables) else data.frame()
    }
  )
}


validate_oracle_checkpoint <- function(
    checkpoint,
    replication_id,
    replication_stream,
    config,
    code_fingerprint) {
  required_entries <- c(
    "version",
    "config",
    "code_fingerprint",
    "replication_id",
    "replication_stream",
    "result"
  )
  
  if (!is.list(checkpoint) ||
      !all(required_entries %in% names(checkpoint))) {
    stop("Oracle checkpoint ", replication_id, " has an incompatible structure.")
  }
  
  checks <- c(
    version = identical(checkpoint$version, restricted_oracle_pilot_version),
    config = identical(checkpoint$config, config),
    fingerprint = identical(checkpoint$code_fingerprint, code_fingerprint),
    replication_id = identical(
      checkpoint$replication_id,
      as.integer(replication_id)
    ),
    stream = identical(checkpoint$replication_stream, replication_stream),
    complete_result =
      is.data.frame(checkpoint$result$reml_summary) &&
      nrow(checkpoint$result$reml_summary) == 4L &&
      is.data.frame(checkpoint$result$oracle_summary) &&
      nrow(checkpoint$result$oracle_summary) == 4L
  )
  
  if (!all(checks)) {
    stop(
      "Oracle checkpoint ", replication_id, " is incompatible: ",
      paste(names(checks)[!checks], collapse = ", "),
      "."
    )
  }
  
  invisible(checkpoint)
}


make_oracle_paired_summary <- function(oracle_summary) {
  cells <- split(oracle_summary, oracle_summary$cell_id)
  
  do.call(
    rbind,
    lapply(cells, function(cell_rows) {
      complete <- cell_rows$oracle_success &
        is.finite(cell_rows$reml_oracle_ise_difference)
      
      delta <- cell_rows$reml_oracle_ise_difference[complete]
      
      data.frame(
        cell_id = cell_rows$cell_id[[1L]],
        attempted = nrow(cell_rows),
        jointly_successful = length(delta),
        mean_paired_aisе_difference = if (length(delta)) mean(delta) else NA_real_,
        paired_mcse = if (length(delta) >= 2L) {
          stats::sd(delta) / sqrt(length(delta))
        } else {
          NA_real_
        },
        oracle_win_rate = if (length(delta)) mean(delta > 0) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )
}


make_oracle_boundary_audit <- function(
    oracle_summary,
    candidate_summary,
    cell_ids) {
  output <- lapply(cell_ids, function(cell_id) {
    cell_oracle <- oracle_summary[
      oracle_summary$cell_id == cell_id &
        !is.na(oracle_summary$coarse_winner_id),
      ,
      drop = FALSE
    ]
    
    winners <- lapply(seq_len(nrow(cell_oracle)), function(i) {
      candidate_summary[
        candidate_summary$cell_id == cell_id &
          candidate_summary$replication == cell_oracle$replication[[i]] &
          candidate_summary$candidate_id ==
          cell_oracle$coarse_winner_id[[i]],
        c("graph_offset", "time_offset"),
        drop = FALSE
      ]
    })
    winners <- Filter(function(x) nrow(x) == 1L, winners)
    
    if (length(winners) == 0L) {
      return(data.frame(
        cell_id = cell_id,
        successful_coarse_searches = 0L,
        graph_lower_boundary_frequency = NA_real_,
        graph_upper_boundary_frequency = NA_real_,
        time_lower_boundary_frequency = NA_real_,
        time_upper_boundary_frequency = NA_real_,
        range_inadequate = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    winners <- do.call(rbind, winners)
    
    frequencies <- c(
      mean(winners$graph_offset == -4),
      mean(winners$graph_offset == 4),
      mean(winners$time_offset == -4),
      mean(winners$time_offset == 4)
    )
    
    data.frame(
      cell_id = cell_id,
      successful_coarse_searches = nrow(winners),
      graph_lower_boundary_frequency = frequencies[[1L]],
      graph_upper_boundary_frequency = frequencies[[2L]],
      time_lower_boundary_frequency = frequencies[[3L]],
      time_upper_boundary_frequency = frequencies[[4L]],
      range_inadequate = any(frequencies > 0.05),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, output)
}


run_restricted_oracle_pilot <- function() {
  design <- make_oracle_pilot_design()
  
  config <- list(
    version = restricted_oracle_pilot_version,
    n_replications = 50L,
    base_seed = 20260915L,
    n_nodes = 25L,
    coarse_graph_offsets = c(-4, -2, 0, 2, 4),
    coarse_time_offsets = c(-4, -2, 0, 2, 4),
    cells = design
  )
  
  code_files <- c(
    "simulation/13_run_restricted_oracle_pilot.R",
    oracle_source_files,
    sort(list.files(
      "netfunsmooth/R",
      pattern = "[.]R$",
      full.names = TRUE
    ))
  )
  code_fingerprint <- oracle_code_fingerprint(code_files)
  
  root <- "simulation/results/restricted_oracle/pilot_50"
  checkpoint_directory <- file.path(root, "checkpoints")
  result_path <- file.path(root, "restricted_oracle_50id_result.rds")
  
  dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
  
  streams <- make_replication_streams(
    base_seed = config$base_seed,
    n_replications = config$n_replications
  )
  
  checkpoints <- vector("list", config$n_replications)
  
  for (replication_id in seq_len(config$n_replications)) {
    checkpoint_path <- file.path(
      checkpoint_directory,
      sprintf("replication_%03d.rds", replication_id)
    )
    
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      
      validate_oracle_checkpoint(
        checkpoint = checkpoint,
        replication_id = replication_id,
        replication_stream = streams[[replication_id]],
        config = config,
        code_fingerprint = code_fingerprint
      )
      
      checkpoints[[replication_id]] <- checkpoint
      
      message(
        sprintf(
          "[Oracle %03d/%03d] validated existing checkpoint.",
          replication_id,
          config$n_replications
        )
      )
      next
    }
    
    message(
      sprintf(
        "[Oracle %03d/%03d] running all four cells.",
        replication_id,
        config$n_replications
      )
    )
    
    start_time <- proc.time()[["elapsed"]]
    
    result <- tryCatch(
      run_oracle_pilot_replication(
        replication_id = replication_id,
        replication_stream = streams[[replication_id]]
      ),
      error = function(error_condition) error_condition
    )
    
    if (inherits(result, "error")) {
      stop(
        "Oracle replication ", replication_id,
        " did not finish atomically: ",
        conditionMessage(result),
        ". No checkpoint was written."
      )
    }
    
    checkpoint <- list(
      version = restricted_oracle_pilot_version,
      config = config,
      code_fingerprint = code_fingerprint,
      replication_id = as.integer(replication_id),
      replication_stream = streams[[replication_id]],
      result = result,
      elapsed_seconds = proc.time()[["elapsed"]] - start_time,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    
    oracle_save_rds_atomically(checkpoint, checkpoint_path)
    checkpoints[[replication_id]] <- checkpoint
    
    first_round_checks <- c(
      all_reml_fits_successful = all(result$reml_summary$success),
      all_oracle_searches_successful = all(result$oracle_summary$oracle_success),
      zero_offset_identity = all(
        is.finite(result$oracle_summary$zero_offset_max_abs_difference) &
          result$oracle_summary$zero_offset_max_abs_difference <= 1e-8
      )
    )
    
    if (!all(first_round_checks)) {
      print(first_round_checks)
      
      stop(
        "The completed replication was saved, but a first-round oracle check ",
        "failed. Pause and inspect this checkpoint before continuing."
      )
    }
  }
  
  reml_summary <- do.call(
    rbind,
    lapply(checkpoints, function(x) x$result$reml_summary)
  )
  oracle_summary <- do.call(
    rbind,
    lapply(checkpoints, function(x) x$result$oracle_summary)
  )
  candidate_summary <- do.call(
    rbind,
    lapply(checkpoints, function(x) x$result$candidate_summary)
  )
  
  paired_summary <- make_oracle_paired_summary(oracle_summary)
  
  boundary_audit <- make_oracle_boundary_audit(
    oracle_summary = oracle_summary,
    candidate_summary = candidate_summary,
    cell_ids = design$cell_id
  )
  
  quality_checks <- c(
    expected_reml_rows = nrow(reml_summary) == 4L * config$n_replications,
    expected_oracle_rows = nrow(oracle_summary) == 4L * config$n_replications,
    no_duplicate_oracle_rows = !anyDuplicated(
      oracle_summary[, c("cell_id", "replication")]
    ),
    all_reml_fits_successful = all(reml_summary$success),
    all_oracle_searches_successful = all(oracle_summary$oracle_success),
    zero_offset_identity = all(
      is.finite(oracle_summary$zero_offset_max_abs_difference) &
        oracle_summary$zero_offset_max_abs_difference <= 1e-8
    ),
    coarse_grid_not_inadequate = !any(boundary_audit$range_inadequate)
  )
  
  final_result <- list(
    version = restricted_oracle_pilot_version,
    config = config,
    code_fingerprint = code_fingerprint,
    quality_checks = quality_checks,
    paired_summary = paired_summary,
    reml_summary = reml_summary,
    oracle_summary = oracle_summary,
    candidate_summary = candidate_summary,
    coarse_boundary_audit = boundary_audit,
    checkpoint_paths = file.path(
      checkpoint_directory,
      sprintf("replication_%03d.rds", seq_len(config$n_replications))
    ),
    session_info = utils::sessionInfo(),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  
  oracle_save_rds_atomically(final_result, result_path)
  
  cat("\nRestricted-oracle B = 50 quality checks:\n")
  print(quality_checks)
  
  cat("\nPaired REML-minus-oracle AISE summary:\n")
  print(paired_summary)
  
  cat("\nCoarse-grid boundary audit:\n")
  print(boundary_audit)
  
  cat("\nResult saved to:\n", result_path, "\n", sep = "")
  
  if (!all(quality_checks)) {
    stop(
      "The B = 50 oracle pilot completed, but at least one quality check failed. ",
      "Results were saved; do not interpret them before inspection."
    )
  }
  
  invisible(final_result)
}


message(
  "Restricted-oracle B = 50 driver loaded. ",
  "Run: run_restricted_oracle_pilot()"
)