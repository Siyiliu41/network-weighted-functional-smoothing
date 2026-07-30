# Generate graph-smooth node-specific true functions.
#
# The truth uses three temporal basis functions. Their coefficients vary
# smoothly over the lattice coordinates, producing differences in phase,
# curvature, and local peak height rather than only vertical shifts.

generate_graph_smooth_truth <- function(
    graph,
    t_grid = seq(0, 1, length.out = 30L),
    signal_scale = 1) {
  
  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.")
  }
  
  if (
    !is.numeric(t_grid) ||
    length(t_grid) < 5L ||
    anyNA(t_grid) ||
    is.unsorted(t_grid, strictly = TRUE)
  ) {
    stop("`t_grid` must be a strictly increasing numeric vector.")
  }
  
  if (
    !is.numeric(signal_scale) ||
    length(signal_scale) != 1L ||
    is.na(signal_scale) ||
    signal_scale <= 0
  ) {
    stop("`signal_scale` must be a single positive number.")
  }
  
  node_names <- igraph::V(graph)$name
  row_coordinate <- igraph::vertex_attr(graph, "row")
  column_coordinate <- igraph::vertex_attr(graph, "column")
  
  if (
    is.null(node_names) ||
    is.null(row_coordinate) ||
    is.null(column_coordinate)
  ) {
    stop(
      "`graph` must contain node names and the vertex attributes ",
      "`row` and `column`."
    )
  }
  
  normalise_coordinate <- function(x) {
    x <- x - mean(x)
    
    maximum <- max(abs(x))
    
    if (maximum == 0) {
      return(x)
    }
    
    x / maximum
  }
  
  x_coordinate <- normalise_coordinate(column_coordinate)
  y_coordinate <- normalise_coordinate(row_coordinate)
  
  # Smooth coefficient fields over the graph.
  # The first two coefficients induce amplitude and phase differences.
  # The third coefficient controls a local peak.
  coefficients <- cbind(
    sine = 0.8 * x_coordinate,
    cosine = 0.8 * y_coordinate,
    peak = 0.6 * x_coordinate * y_coordinate
  )
  
  # Centre every coefficient field across nodes. This ensures that the
  # node-specific deviations sum to zero at every time point.
  coefficients <- sweep(
    coefficients,
    MARGIN = 2,
    STATS = colMeans(coefficients),
    FUN = "-"
  )
  
  temporal_basis <- cbind(
    sine = sin(2 * pi * t_grid),
    cosine = cos(2 * pi * t_grid),
    peak = exp(
      -0.5 * ((t_grid - 0.65) / 0.10)^2
    )
  )
  
  mean_curve <-
    0.4 +
    0.6 * sin(pi * t_grid) +
    0.2 * cos(2 * pi * t_grid)
  
  node_deviations <-
    signal_scale *
    coefficients %*% t(temporal_basis)
  
  true_curves <- sweep(
    node_deviations,
    MARGIN = 2,
    STATS = mean_curve,
    FUN = "+"
  )
  
  rownames(coefficients) <- node_names
  rownames(node_deviations) <- node_names
  rownames(true_curves) <- node_names
  
  colnames(temporal_basis) <- colnames(coefficients)
  
  list(
    truth = true_curves,
    mean_curve = mean_curve,
    deviations = node_deviations,
    coefficients = coefficients,
    temporal_basis = temporal_basis,
    t_grid = t_grid,
    graph = graph
  )
}