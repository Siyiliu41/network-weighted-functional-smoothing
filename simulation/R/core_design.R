# Core factorial design and common random numbers for the pilot simulation.

make_core_design <- function() {
  design <- expand.grid(
    truth_structure = c("coordinate", "cluster"),
    alpha = c(0.5, 0.9),
    neighbourhood = c("rook", "queen"),
    sigma = c(0.2, 0.5),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design$cell_id <- sprintf(
    "%s_a%02d_%s_s%02d",
    design$truth_structure,
    round(100 * design$alpha),
    design$neighbourhood,
    round(100 * design$sigma)
  )
  
  design <- design[
    order(
      design$truth_structure,
      design$alpha,
      design$neighbourhood,
      design$sigma
    ),
  ]
  
  rownames(design) <- NULL
  
  design[
    ,
    c(
      "cell_id",
      "truth_structure",
      "alpha",
      "neighbourhood",
      "sigma"
    )
  ]
}


generate_core_random_draws <- function(
    seed,
    n_nodes = 25L,
    n_time = 50L,
    n_components = 3L) {
  
  if (
    !is.numeric(seed) ||
    length(seed) != 1L ||
    is.na(seed) ||
    !is.finite(seed)
  ) {
    stop("`seed` must be one finite number.")
  }
  
  dimensions <- c(
    n_nodes = n_nodes,
    n_time = n_time,
    n_components = n_components
  )
  
  if (
    any(!is.numeric(dimensions)) ||
    anyNA(dimensions) ||
    any(dimensions < 1) ||
    any(dimensions != as.integer(dimensions))
  ) {
    stop(
      "`n_nodes`, `n_time`, and `n_components` ",
      "must be positive integers."
    )
  }
  
  set.seed(as.integer(seed))
  
  list(
    coefficient_draws = matrix(
      stats::rnorm(n_nodes * n_components),
      nrow = n_nodes,
      ncol = n_components
    ),
    error_draws = matrix(
      stats::rnorm(n_nodes * n_time),
      nrow = n_nodes,
      ncol = n_time
    )
  )
}


generate_core_scenarios <- function(
    seed,
    n_side = 5L,
    t_grid = seq(0, 1, length.out = 50L)) {
  
  design <- make_core_design()
  
  graphs <- list(
    rook = make_lattice_graph(
      n_side = n_side,
      neighbourhood = "rook"
    ),
    queen = make_lattice_graph(
      n_side = n_side,
      neighbourhood = "queen"
    )
  )
  
  node_names <- igraph::V(graphs$rook)$name
  
  if (!identical(node_names, igraph::V(graphs$queen)$name)) {
    stop("Rook and queen graphs must have identical node ordering.")
  }
  
  random_draws <- generate_core_random_draws(
    seed = seed,
    n_nodes = length(node_names),
    n_time = length(t_grid),
    n_components = 3L
  )
  
  truth_settings <- unique(
    design[, c("truth_structure", "alpha")]
  )
  
  truth_objects <- vector(
    mode = "list",
    length = nrow(truth_settings)
  )
  
  truth_names <- sprintf(
    "%s_a%02d",
    truth_settings$truth_structure,
    round(100 * truth_settings$alpha)
  )
  
  names(truth_objects) <- truth_names
  
  for (i in seq_len(nrow(truth_settings))) {
    truth_objects[[i]] <- generate_core_truth(
      # The truth depends on coordinates, not on graph adjacency.
      graph = graphs$rook,
      t_grid = t_grid,
      truth_structure = truth_settings$truth_structure[i],
      alpha = truth_settings$alpha[i],
      standard_normal_draws = random_draws$coefficient_draws
    )
  }
  
  scenarios <- vector(
    mode = "list",
    length = nrow(design)
  )
  
  names(scenarios) <- design$cell_id
  
  for (i in seq_len(nrow(design))) {
    truth_name <- sprintf(
      "%s_a%02d",
      design$truth_structure[i],
      round(100 * design$alpha[i])
    )
    
    truth_object <- truth_objects[[truth_name]]
    
    simulated_data <- simulate_observed_curves(
      truth_object = truth_object,
      sigma = design$sigma[i],
      standard_normal_errors = random_draws$error_draws
    )
    
    scenarios[[i]] <- list(
      cell = design[i, , drop = FALSE],
      graph = graphs[[design$neighbourhood[i]]],
      truth_object = truth_object,
      simulated_data = simulated_data
    )
  }
  
  list(
    seed = seed,
    design = design,
    graphs = graphs,
    random_draws = random_draws,
    truth_objects = truth_objects,
    scenarios = scenarios
  )
}