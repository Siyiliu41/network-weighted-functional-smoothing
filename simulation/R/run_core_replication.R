# Run all 16 core-design cells for one Monte Carlo replication.

run_core_replication <- function(
    replication_id,
    seed,
    replication_stream,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    k_pooled = 10L,
    keep_predictions = FALSE,
    show_progress = TRUE) {

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

  validate_rng_stream(replication_stream)

  replication_id <- as.integer(replication_id)
  seed <- as.integer(seed)

  core_scenarios <- generate_core_scenarios(
    replication_stream = replication_stream
  )

  scenario_names <- names(core_scenarios$scenarios)

  cell_results <- vector(
    mode = "list",
    length = length(scenario_names)
  )

  names(cell_results) <- scenario_names

  replication_start <- proc.time()[["elapsed"]]

  for (i in seq_along(scenario_names)) {
    cell_id <- scenario_names[i]

    if (show_progress) {
      message(
        sprintf(
          "[Replication %d] Cell %d/%d: %s",
          replication_id,
          i,
          length(scenario_names),
          cell_id
        )
      )
    }

    cell_results[[cell_id]] <- run_core_cell(
      scenario = core_scenarios$scenarios[[cell_id]],
      replication_id = replication_id,
      seed = seed,
      k_intercept = k_intercept,
      k_deviation = k_deviation,
      k_nodewise = k_nodewise,
      k_pooled = k_pooled,
      keep_predictions = keep_predictions
    )
  }

  replication_elapsed <- (
    proc.time()[["elapsed"]] - replication_start
  )

  summary_table <- do.call(
    rbind,
    lapply(
      cell_results,
      function(cell_result) {
        cell_result$summary
      }
    )
  )

  rownames(summary_table) <- NULL

  expected_rows <- nrow(core_scenarios$design) * 4L

  if (nrow(summary_table) != expected_rows) {
    stop(
      "Unexpected number of summary rows: expected ",
      expected_rows,
      ", obtained ",
      nrow(summary_table),
      "."
    )
  }

  expected_methods <- c(
    "raw",
    "pooled",
    "nodewise",
    "network"
  )

  method_counts <- table(summary_table$method)

  if (
    !all(expected_methods %in% names(method_counts)) ||
    any(method_counts[expected_methods] != nrow(core_scenarios$design))
  ) {
    stop("Each method must occur once in every design cell.")
  }

  list(
    replication = replication_id,
    seed = seed,
    replication_stream = replication_stream,
    design = core_scenarios$design,
    basis_dimensions = list(
      k_intercept = as.integer(k_intercept),
      k_deviation = as.integer(k_deviation),
      k_nodewise = as.integer(k_nodewise),
      k_pooled = as.integer(k_pooled)
    ),
    summary = summary_table,
    cell_results = cell_results,
    elapsed_seconds = replication_elapsed
  )
}