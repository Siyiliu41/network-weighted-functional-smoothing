# Construct and run the two fixed N = 100 scalability cells.
#
# Both cells use the same L'Ecuyer replication stream, coefficient draws, and
# standard-normal observation-error matrix. Only alpha differs.

make_n100_design <- function() {
  data.frame(
    cell_id = c(
      "n100_coordinate_a90_rook_s50",
      "n100_coordinate_a00_rook_s50"
    ),
    truth_structure = c("coordinate", "coordinate"),
    alpha = c(0.9, 0),
    neighbourhood = c("rook", "rook"),
    sigma = c(0.5, 0.5),
    stringsAsFactors = FALSE
  )
}


validate_n100_design <- function(design) {
  expected <- make_n100_design()
  
  if (!is.data.frame(design) ||
      !identical(design, expected)) {
    stop("The N = 100 design differs from the two prespecified fixed cells.")
  }
  
  invisible(design)
}


make_n100_random_draws <- function(replication_stream) {
  validate_rng_stream(replication_stream)
  
  component_streams <- make_component_substreams(replication_stream)
  
  list(
    coefficient_draws = with_rng_stream(
      component_streams$truth_draws,
      matrix(
        stats::rnorm(100L * 3L),
        nrow = 100L,
        ncol = 3L
      )
    ),
    standard_normal_errors = with_rng_stream(
      component_streams$error_coordinate_a09,
      matrix(
        stats::rnorm(100L * 50L),
        nrow = 100L,
        ncol = 50L
      )
    ),
    component_streams = component_streams
  )
}


make_n100_scenario <- function(cell, graph, random_draws) {
  if (!is.data.frame(cell) || nrow(cell) != 1L) {
    stop("`cell` must be one row from `make_n100_design()`.")
  }
  
  required_draws <- c(
    "coefficient_draws",
    "standard_normal_errors"
  )
  
  if (!is.list(random_draws) ||
      !all(required_draws %in% names(random_draws))) {
    stop("`random_draws` is incomplete.")
  }
  
  truth_object <- generate_core_truth(
    graph = graph,
    t_grid = seq(0, 1, length.out = 50L),
    truth_structure = cell$truth_structure,
    alpha = cell$alpha,
    standard_normal_draws = random_draws$coefficient_draws,
    coefficient_scales = c(
      0.8 / sqrt(2),
      0.8 / sqrt(2),
      0.3
    )
  )
  
  simulated_data <- simulate_observed_curves(
    truth_object = truth_object,
    sigma = cell$sigma,
    standard_normal_errors = random_draws$standard_normal_errors
  )
  
  list(
    cell = cell,
    graph = graph,
    truth_object = truth_object,
    simulated_data = simulated_data
  )
}


run_n100_replication <- function(
    replication_id,
    replication_stream,
    seed,
    k_intercept = 10L,
    k_deviation = 15L,
    k_nodewise = 15L,
    keep_predictions = FALSE,
    show_progress = TRUE) {
  
  if (!is.numeric(replication_id) ||
      length(replication_id) != 1L ||
      is.na(replication_id) ||
      replication_id < 1L ||
      replication_id %% 1L != 0L) {
    stop("`replication_id` must be a positive integer.")
  }
  
  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed)) {
    stop("`seed` must be one finite number.")
  }
  
  validate_rng_stream(replication_stream)
  
  design <- make_n100_design()
  validate_n100_design(design)
  
  graph <- make_lattice_graph(
    n_side = 10L,
    neighbourhood = "rook"
  )
  
  if (igraph::vcount(graph) != 100L ||
      igraph::ecount(graph) != 180L) {
    stop("The N = 100 rook lattice must have 100 nodes and 180 edges.")
  }
  
  random_draws <- make_n100_random_draws(replication_stream)
  
  cell_results <- vector("list", nrow(design))
  names(cell_results) <- design$cell_id
  cell_elapsed_seconds <- numeric(nrow(design))
  names(cell_elapsed_seconds) <- design$cell_id
  
  for (cell_index in seq_len(nrow(design))) {
    cell <- design[cell_index, , drop = FALSE]
    
    if (show_progress) {
      message(sprintf(
        "[N100 replication %03d] %s",
        as.integer(replication_id),
        cell$cell_id
      ))
    }
    
    cell_start <- proc.time()[["elapsed"]]
    
    scenario <- make_n100_scenario(
      cell = cell,
      graph = graph,
      random_draws = random_draws
    )
    
    cell_results[[cell$cell_id]] <- run_n100_cell(
      scenario = scenario,
      replication_id = replication_id,
      seed = seed,
      k_intercept = k_intercept,
      k_deviation = k_deviation,
      k_nodewise = k_nodewise,
      keep_predictions = keep_predictions
    )
    
    cell_elapsed_seconds[[cell$cell_id]] <-
      proc.time()[["elapsed"]] - cell_start
  }
  
  summary_table <- do.call(
    rbind,
    lapply(cell_results, `[[`, "summary")
  )
  rownames(summary_table) <- NULL
  
  if (nrow(summary_table) != 4L ||
      anyDuplicated(
        summary_table[, c("cell_id", "replication", "method")]
      )) {
    stop("An N = 100 replication must contain two cells and two methods per cell.")
  }
  
  list(
    replication = as.integer(replication_id),
    seed = as.integer(seed),
    replication_stream = replication_stream,
    design = design,
    basis_dimensions = list(
      k_intercept = as.integer(k_intercept),
      k_deviation = as.integer(k_deviation),
      k_nodewise = as.integer(k_nodewise)
    ),
    summary = summary_table,
    cell_elapsed_seconds = cell_elapsed_seconds,
    cell_results = cell_results
  )
}