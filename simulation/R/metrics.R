# Performance measures for reconstruction of node-specific functions.

trapezoid_weights <- function(t_grid) {
  
  if (
    !is.numeric(t_grid) ||
    length(t_grid) < 2L ||
    anyNA(t_grid) ||
    any(diff(t_grid) <= 0)
  ) {
    stop("`t_grid` must be a strictly increasing numeric vector.")
  }
  
  n_grid <- length(t_grid)
  
  weights <- numeric(n_grid)
  
  weights[1] <-
    (t_grid[2] - t_grid[1]) / 2
  
  weights[n_grid] <-
    (t_grid[n_grid] - t_grid[n_grid - 1L]) / 2
  
  if (n_grid > 2L) {
    weights[2:(n_grid - 1L)] <-
      (
        t_grid[3:n_grid] -
          t_grid[1:(n_grid - 2L)]
      ) / 2
  }
  
  weights
}


calculate_node_ise <- function(
    estimate,
    truth,
    t_grid) {
  
  if (!is.matrix(estimate) || !is.matrix(truth)) {
    stop("`estimate` and `truth` must both be matrices.")
  }
  
  if (!identical(dim(estimate), dim(truth))) {
    stop("`estimate` and `truth` must have identical dimensions.")
  }
  
  if (ncol(truth) != length(t_grid)) {
    stop(
      "The number of matrix columns must equal the length of `t_grid`."
    )
  }
  
  if (anyNA(estimate) || anyNA(truth)) {
    stop("`estimate` and `truth` must not contain missing values.")
  }
  
  truth_names <- rownames(truth)
  estimate_names <- rownames(estimate)
  
  if (
    !is.null(truth_names) &&
    !is.null(estimate_names)
  ) {
    if (!setequal(truth_names, estimate_names)) {
      stop(
        "The node names of `estimate` and `truth` do not agree."
      )
    }
    
    estimate <- estimate[
      truth_names,
      ,
      drop = FALSE
    ]
  }
  
  weights <- trapezoid_weights(t_grid)
  
  squared_error <- (estimate - truth)^2
  
  node_ise <- as.numeric(
    squared_error %*% weights
  )
  
  names(node_ise) <- truth_names
  
  node_ise
}


summarise_ise <- function(node_ise) {
  
  if (
    !is.numeric(node_ise) ||
    length(node_ise) == 0L ||
    anyNA(node_ise)
  ) {
    stop("`node_ise` must be a non-empty numeric vector.")
  }
  
  data.frame(
    node_averaged_ise = mean(node_ise),
    median_node_ise = stats::median(node_ise),
    worst_node_ise = max(node_ise)
  )
}