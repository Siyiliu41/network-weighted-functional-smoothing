# Sequential pilot simulation for the 16-cell core design.
#
# The pilot:
# - runs 50 Monte Carlo replications;
# - uses independent L'Ecuyer-CMRG replication streams;
# - saves one checkpoint per replication;
# - resumes from compatible checkpoints;
# - does not contribute to the final simulation estimates.


# -------------------------------------------------------------------------
# Load package and simulation functions
# -------------------------------------------------------------------------

devtools::load_all(
  "netfunsmooth",
  quiet = TRUE
)

simulation_sources <- c(
  "simulation/R/make_graph.R",
  "simulation/R/generate_truth.R",
  "simulation/R/simulate_data.R",
  "simulation/R/fit_baselines.R",
  "simulation/R/metrics.R",
  "simulation/R/rng.R",
  "simulation/R/core_design.R",
  "simulation/R/run_replication.R",
  "simulation/R/run_core_cell.R",
  "simulation/R/run_core_replication.R"
)

invisible(
  lapply(
    simulation_sources,
    source
  )
)


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

base_seed <- 20260805L
n_replications <- 50L

replication_seed_labels <- as.integer(
  base_seed + seq_len(n_replications) - 1L
)

code_files <- c(
  "simulation/02_pilot_simulation.R",
  simulation_sources,
  sort(
    list.files(
      "netfunsmooth/R",
      pattern = "\\.R$",
      full.names = TRUE
    )
  )
)

pilot_config <- list(
  pilot_version = "2026-08-06-v1",
  base_seed = base_seed,
  n_replications = n_replications,
  replication_seed_labels = replication_seed_labels,
  k_intercept = 10L,
  k_deviation = 15L,
  k_nodewise = 15L,
  k_pooled = 10L,
  code_md5 = tools::md5sum(code_files)
)

checkpoint_directory <-
  "simulation/results/core/pilot/checkpoints"

final_result_file <-
  "simulation/results/core/pilot/pilot_result.rds"

dir.create(
  checkpoint_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# Reproducible replication streams
# -------------------------------------------------------------------------

replication_streams <- make_replication_streams(
  base_seed = pilot_config$base_seed,
  n_replications = pilot_config$n_replications
)


# -------------------------------------------------------------------------
# Helper for safely writing checkpoints
# -------------------------------------------------------------------------

save_checkpoint <- function(object, path) {
  temporary_path <- tempfile(
    pattern = "pilot-checkpoint-",
    tmpdir = dirname(path),
    fileext = ".tmp"
  )

  on.exit(
    {
      if (file.exists(temporary_path)) {
        unlink(temporary_path)
      }
    },
    add = TRUE
  )

  saveRDS(
    object,
    file = temporary_path
  )

  if (!file.rename(temporary_path, path)) {
    stop(
      "Could not move the temporary checkpoint to: ",
      path
    )
  }

  invisible(path)
}


# -------------------------------------------------------------------------
# Run or restore replications
# -------------------------------------------------------------------------

replication_results <- vector(
  mode = "list",
  length = pilot_config$n_replications
)

pilot_start <- proc.time()[["elapsed"]]

for (
  replication_id in
  seq_len(pilot_config$n_replications)
) {
  checkpoint_file <- file.path(
    checkpoint_directory,
    sprintf(
      "replication_%03d.rds",
      replication_id
    )
  )

  if (file.exists(checkpoint_file)) {
    checkpoint <- tryCatch(
      readRDS(checkpoint_file),
      error = function(error) {
        stop(
          "Checkpoint could not be read: ",
          checkpoint_file,
          "\nDelete only this checkpoint and rerun the pilot.\n",
          "Original error: ",
          conditionMessage(error)
        )
      }
    )

    if (!identical(
      checkpoint$config,
      pilot_config
    )) {
      stop(
        "Checkpoint configuration or code fingerprint does not ",
        "match the current pilot:\n",
        checkpoint_file,
        "\nDo not combine results generated from different code versions."
      )
    }

    if (!identical(
      checkpoint$replication_stream,
      replication_streams[[replication_id]]
    )) {
      stop(
        "RNG stream mismatch in checkpoint:\n",
        checkpoint_file
      )
    }

    replication_results[[replication_id]] <-
      checkpoint$result

    message(
      sprintf(
        "[Pilot] Loaded replication %d/%d from checkpoint.",
        replication_id,
        pilot_config$n_replications
      )
    )

    next
  }

  message(
    sprintf(
      "[Pilot] Running replication %d/%d.",
      replication_id,
      pilot_config$n_replications
    )
  )

  replication_result <- run_core_replication(
    replication_id = replication_id,

    # This scalar is retained as a human-readable result label.
    # The actual random numbers come from replication_stream.
    seed = pilot_config$replication_seed_labels[
      replication_id
    ],

    replication_stream =
      replication_streams[[replication_id]],

    k_intercept = pilot_config$k_intercept,
    k_deviation = pilot_config$k_deviation,
    k_nodewise = pilot_config$k_nodewise,
    k_pooled = pilot_config$k_pooled,
    keep_predictions = FALSE,
    show_progress = TRUE
  )

  if (!identical(
    replication_result$replication_stream,
    replication_streams[[replication_id]]
  )) {
    stop(
      "The saved RNG stream does not match replication ",
      replication_id,
      "."
    )
  }

  checkpoint <- list(
    config = pilot_config,
    replication_id = replication_id,
    replication_stream =
      replication_streams[[replication_id]],
    result = replication_result
  )

  save_checkpoint(
    checkpoint,
    checkpoint_file
  )

  replication_results[[replication_id]] <-
    replication_result

  message(
    sprintf(
      "[Pilot] Saved checkpoint %d/%d.",
      replication_id,
      pilot_config$n_replications
    )
  )
}

pilot_elapsed_seconds <-
  proc.time()[["elapsed"]] - pilot_start


# -------------------------------------------------------------------------
# Combine method-level summaries
# -------------------------------------------------------------------------

pilot_summary <- do.call(
  rbind,
  lapply(
    replication_results,
    function(result) {
      result$summary
    }
  )
)

rownames(pilot_summary) <- NULL

expected_cells <- 16L
expected_methods <- 4L
expected_rows_per_replication <-
  expected_cells * expected_methods

network_rows <- pilot_summary[
  pilot_summary$method == "network",
  ,
  drop = FALSE
]


# -------------------------------------------------------------------------
# Pilot integrity checks
# -------------------------------------------------------------------------

rows_per_replication <- table(
  pilot_summary$replication
)

method_counts <- table(
  pilot_summary$method
)

streams_saved_correctly <- all(
  vapply(
    seq_along(replication_results),
    function(i) {
      identical(
        replication_results[[i]]$replication_stream,
        replication_streams[[i]]
      )
    },
    logical(1)
  )
)

pilot_checks <- list(
  expected_total_rows =
    nrow(pilot_summary) ==
    pilot_config$n_replications *
    expected_rows_per_replication,

  expected_number_of_replications =
    length(unique(pilot_summary$replication)) ==
    pilot_config$n_replications,

  expected_rows_per_replication =
    length(rows_per_replication) ==
    pilot_config$n_replications &&
    all(
      rows_per_replication ==
        expected_rows_per_replication
    ),

  expected_number_of_cells =
    length(unique(pilot_summary$cell_id)) ==
    expected_cells,

  expected_method_counts =
    all(
      c(
        "raw",
        "pooled",
        "nodewise",
        "network"
      ) %in% names(method_counts)
    ) &&
    all(
      method_counts[
        c(
          "raw",
          "pooled",
          "nodewise",
          "network"
        )
      ] ==
        pilot_config$n_replications *
        expected_cells
    ),

  all_fits_successful =
    isTRUE(all(pilot_summary$success)),

  no_warnings =
    sum(pilot_summary$warning_count) == 0L,

  no_error_messages =
    all(is.na(pilot_summary$error_message)),

  all_network_fits_converged =
    all(network_rows$converged %in% TRUE),

  all_ise_values_finite =
    all(
      is.finite(
        pilot_summary$node_averaged_ise
      )
    ),

  streams_saved_correctly =
    streams_saved_correctly
)

pilot_passed <- all(
  unlist(
    pilot_checks,
    use.names = FALSE
  )
)


# -------------------------------------------------------------------------
# Descriptive pilot summaries
# -------------------------------------------------------------------------

mean_ise_by_method <- stats::aggregate(
  node_averaged_ise ~ method,
  data = pilot_summary,
  FUN = mean,
  na.rm = TRUE
)

mean_runtime_by_method <- stats::aggregate(
  elapsed_seconds ~ method,
  data = pilot_summary[
    !is.na(pilot_summary$elapsed_seconds),
    ,
    drop = FALSE
  ],
  FUN = mean,
  na.rm = TRUE
)


# -------------------------------------------------------------------------
# Save combined pilot result
# -------------------------------------------------------------------------

pilot_result <- list(
  config = pilot_config,
  pilot_passed = pilot_passed,
  checks = pilot_checks,
  summary = pilot_summary,
  mean_ise_by_method = mean_ise_by_method,
  mean_runtime_by_method = mean_runtime_by_method,
  replication_elapsed_seconds = vapply(
    replication_results,
    function(result) {
      result$elapsed_seconds
    },
    numeric(1)
  ),
  total_elapsed_seconds = pilot_elapsed_seconds,
  checkpoint_directory = checkpoint_directory,
  session_info = utils::sessionInfo()
)

saveRDS(
  pilot_result,
  file = final_result_file
)


# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

cat("\nPilot integrity checks:\n")

print(
  unlist(pilot_checks)
)

cat("\nMean ISE by method:\n")

print(
  mean_ise_by_method,
  row.names = FALSE
)

cat("\nMean runtime by method:\n")

print(
  mean_runtime_by_method,
  row.names = FALSE
)

cat(
  sprintf(
    "\nTotal elapsed time: %.2f minutes\n",
    pilot_elapsed_seconds / 60
  )
)

cat(
  "\nCombined pilot result saved to:\n",
  final_result_file,
  "\n"
)

if (!pilot_passed) {
  stop(
    "The pilot completed and its results were saved, ",
    "but at least one integrity check failed."
  )
}

message("Pilot simulation completed successfully.")