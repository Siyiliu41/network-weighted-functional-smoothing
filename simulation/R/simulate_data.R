# Add Gaussian observation noise and construct a tfd vector.

simulate_observed_curves <- function(
    truth_object,
    sigma = 0.4,
    seed = NULL) {
  
  required_components <- c(
    "truth",
    "t_grid",
    "graph"
  )
  
  if (
    !is.list(truth_object) ||
    !all(required_components %in% names(truth_object))
  ) {
    stop(
      "`truth_object` must contain `truth`, `t_grid`, and `graph`."
    )
  }
  
  true_curves <- truth_object$truth
  t_grid <- truth_object$t_grid
  
  if (
    !is.matrix(true_curves) ||
    ncol(true_curves) != length(t_grid)
  ) {
    stop(
      "`truth_object$truth` must be a matrix with one column ",
      "for every value in `t_grid`."
    )
  }
  
  if (
    !is.numeric(sigma) ||
    length(sigma) != 1L ||
    is.na(sigma) ||
    sigma < 0
  ) {
    stop("`sigma` must be a single non-negative number.")
  }
  
  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed)
    ) {
      stop("`seed` must be NULL or a single number.")
    }
    
    had_rng_state <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    
    if (had_rng_state) {
      previous_rng_state <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }
    
    on.exit(
      {
        if (had_rng_state) {
          assign(
            ".Random.seed",
            previous_rng_state,
            envir = .GlobalEnv
          )
        } else if (
          exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )
        ) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )
    
    set.seed(seed)
  }
  
  noise <- matrix(
    stats::rnorm(
      n = length(true_curves),
      mean = 0,
      sd = sigma
    ),
    nrow = nrow(true_curves),
    ncol = ncol(true_curves)
  )
  
  observed_curves <- true_curves + noise
  
  node_names <- rownames(true_curves)
  
  rownames(noise) <- node_names
  rownames(observed_curves) <- node_names
  
  curves_tfd <- tf::tfd(
    observed_curves,
    arg = t_grid
  )
  
  names(curves_tfd) <- node_names
  
  list(
    curves = curves_tfd,
    observed = observed_curves,
    truth = true_curves,
    noise = noise,
    t_grid = t_grid,
    graph = truth_object$graph,
    sigma = sigma,
    seed = seed
  )
}