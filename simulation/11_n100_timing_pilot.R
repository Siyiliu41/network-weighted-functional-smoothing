# Five-replication timing pilot for the two fixed N = 100 cells.
#
# Usage:
#   source("simulation/11_n100_timing_pilot.R")
#   run_n100_timing_pilot()

n100_timing_script_version <- "2026-08-22-n100-timing-v1"

devtools::load_all("netfunsmooth", quiet = TRUE)

n100_sources <- c(
  "simulation/R/make_graph.R",
  "simulation/R/rng.R",
  "simulation/R/generate_truth.R",
  "simulation/R/simulate_data.R",
  "simulation/R/model_diagnostics.R",
  "simulation/R/fit_baselines.R",
  "simulation/R/metrics.R",
  "simulation/R/run_replication.R",
  "simulation/R/run_n100_cell.R",
  "simulation/R/run_n100_replication.R",
  "simulation/R/n100_utils.R"
)

invisible(lapply(n100_sources, source))


run_n100_timing_pilot <- function() {
  config <- list(
    version = n100_timing_script_version,
    phase = "timing_pilot",
    base_seed = 20260822L,
    n_timing_replications = 5L,
    n_side = 10L,
    n_grid = 50L,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    candidate_replication_counts = c(100L, 50L, 25L),
    serial_runtime_budget_hours = 8
  )
  
  code_files <- c(
    "simulation/11_n100_timing_pilot.R",
    n100_sources,
    sort(list.files(
      "netfunsmooth/R",
      pattern = "\\.R$",
      full.names = TRUE
    ))
  )
  
  code_fingerprint <- make_n100_code_fingerprint(code_files)
  
  result_root <- "simulation/results/n100/timing_pilot"
  checkpoint_directory <- file.path(result_root, "checkpoints")
  result_path <- file.path(
    result_root,
    "n100_timing_pilot_result.rds"
  )
  
  if (file.exists(result_path)) {
    previous <- readRDS(result_path)
    
    checks <- c(
      config = identical(previous$config, config),
      fingerprint = identical(
        previous$code_fingerprint,
        code_fingerprint
      )
    )
    
    if (!all(checks)) {
      stop(
        "Existing N = 100 timing result is incompatible: ",
        paste(names(checks)[!checks], collapse = ", "),
        "."
      )
    }
    
    message("Validated existing N = 100 timing-pilot result: ", result_path)
    return(invisible(previous))
  }
  
  replication_streams <- make_replication_streams(
    base_seed = config$base_seed,
    n_replications = config$n_timing_replications
  )
  
  checkpoints <- vector(
    "list",
    config$n_timing_replications
  )
  
  for (replication_id in seq_len(config$n_timing_replications)) {
    checkpoint_path <- file.path(
      checkpoint_directory,
      sprintf("replication_%03d.rds", replication_id)
    )
    
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      
      validate_n100_checkpoint(
        checkpoint = checkpoint,
        replication_id = replication_id,
        replication_stream = replication_streams[[replication_id]],
        config = config,
        code_fingerprint = code_fingerprint
      )
      
      checkpoints[[replication_id]] <- checkpoint
      
      message(sprintf(
        "[N100 timing %d/%d] validated existing checkpoint.",
        replication_id,
        config$n_timing_replications
      ))
      
      next
    }
    
    message(sprintf(
      "[N100 timing %d/%d] running.",
      replication_id,
      config$n_timing_replications
    ))
    
    started_at <- proc.time()[["elapsed"]]
    
    result <- tryCatch(
      run_n100_replication(
        replication_id = replication_id,
        replication_stream = replication_streams[[replication_id]],
        seed = config$base_seed + replication_id - 1L,
        k_intercept = config$k_intercept,
        k_deviation = config$k_deviation,
        k_nodewise = config$k_nodewise,
        keep_predictions = FALSE,
        show_progress = TRUE
      ),
      error = function(condition) condition
    )
    
    if (inherits(result, "error")) {
      stop(
        "Timing replication ",
        replication_id,
        " failed before an atomic checkpoint could be written: ",
        conditionMessage(result)
      )
    }
    
    checkpoint <- list(
      version = config$version,
      config = config,
      code_fingerprint = code_fingerprint,
      replication_id = as.integer(replication_id),
      replication_stream = replication_streams[[replication_id]],
      result = result,
      elapsed_seconds = proc.time()[["elapsed"]] - started_at,
      completed_at = format(
        Sys.time(),
        tz = "UTC",
        usetz = TRUE
      )
    )
    
    save_n100_rds_atomically(
      checkpoint,
      checkpoint_path
    )
    
    checkpoints[[replication_id]] <- checkpoint
  }
  
  timing_table <- make_n100_timing_table(checkpoints)
  
  replication_choice <- choose_n100_replication_count(
    timing_table = timing_table,
    budget_hours = config$serial_runtime_budget_hours,
    n_timing_replications = config$n_timing_replications
  )
  
  timing_result <- list(
    version = config$version,
    config = config,
    code_fingerprint = code_fingerprint,
    timing_table = timing_table,
    replication_choice = replication_choice,
    checkpoint_paths = file.path(
      checkpoint_directory,
      sprintf(
        "replication_%03d.rds",
        seq_len(config$n_timing_replications)
      )
    ),
    session_info = utils::sessionInfo(),
    completed_at = format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    )
  )
  
  save_n100_rds_atomically(
    timing_result,
    result_path
  )
  
  print(replication_choice)
  
  if (!replication_choice$computationally_feasible) {
    stop(
      "Even B = 25 is infeasible under the eight-hour serial budget. ",
      "Pause the N = 100 main block."
    )
  }
  
  message("N = 100 timing pilot completed: ", result_path)
  
  invisible(timing_result)
}