# Generate graph-smooth node-specific true functions.
#
# The truth uses three temporal basis functions. Their coefficients vary
# smoothly over the lattice coordinates, producing differences in phase,
# curvature, and local peak height rather than only vertical shifts.

# Functions for the core pilot data-generating mechanisms.

standardise_population <- function(x, label = "component") {
  x <- as.numeric(x)
  x <- x - mean(x)
  
  component_scale <- sqrt(mean(x^2))
  
  if (
    !is.finite(component_scale) ||
    component_scale <= sqrt(.Machine$double.eps)
  ) {
    stop(
      "Cannot standardise `",
      label,
      "`: its variance is zero."
    )
  }
  
  x / component_scale
}


make_temporal_basis <- function(t_grid) {
  if (
    !is.numeric(t_grid) ||
    length(t_grid) < 5L ||
    anyNA(t_grid) ||
    any(diff(t_grid) <= 0)
  ) {
    stop(
      "`t_grid` must be a strictly increasing numeric vector."
    )
  }
  
  cbind(
    sine = sin(2 * pi * t_grid),
    cosine = cos(2 * pi * t_grid),
    peak = exp(-100 * (t_grid - 0.65)^2)
  )
}


make_structured_fields <- function(
    graph,
    truth_structure = c("coordinate", "cluster")) {
  
  truth_structure <- match.arg(truth_structure)
  
  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.")
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
    x / max(abs(x))
  }
  
  x_coordinate <- normalise_coordinate(column_coordinate)
  y_coordinate <- normalise_coordinate(row_coordinate)
  
  coordinate_raw <- cbind(
    sine = 0.8 * x_coordinate,
    cosine = 0.8 * y_coordinate,
    peak = 0.6 * x_coordinate * y_coordinate
  )
  
  # These scales are also used for the cluster truth so that both truth
  # structures have comparable marginal coefficient magnitudes.
  centred_coordinate_raw <- sweep(
    coordinate_raw,
    MARGIN = 2,
    STATS = colMeans(coordinate_raw),
    FUN = "-"
  )
  
  coefficient_scales <- sqrt(
    colMeans(centred_coordinate_raw^2)
  )
  
  if (truth_structure == "coordinate") {
    raw_fields <- coordinate_raw
    cluster <- rep(NA_integer_, length(node_names))
  } else {
    # Columns 1--2 form cluster 1; columns 3--5 form cluster 2.
    cluster <- ifelse(column_coordinate <= 2L, 1L, 2L)
    
    contrasts <- rbind(
      c(1, 1, -1),
      c(-1, -1, 1)
    )
    
    colnames(contrasts) <- colnames(coordinate_raw)
    
    raw_fields <- contrasts[cluster, , drop = FALSE]
  }
  
  structured <- vapply(
    seq_len(ncol(raw_fields)),
    function(k) {
      standardise_population(
        raw_fields[, k],
        label = paste0("structured field ", k)
      )
    },
    numeric(nrow(raw_fields))
  )
  
  dimnames(structured) <- list(
    node_names,
    colnames(raw_fields)
  )
  
  list(
    structured = structured,
    raw_fields = raw_fields,
    coefficient_scales = coefficient_scales,
    cluster = cluster,
    coordinates = data.frame(
      node = node_names,
      row = row_coordinate,
      column = column_coordinate,
      stringsAsFactors = FALSE
    )
  )
}


make_unstructured_fields <- function(
    structured,
    standard_normal_draws) {
  
  if (!is.matrix(structured)) {
    stop("`structured` must be a matrix.")
  }
  
  if (
    !is.matrix(standard_normal_draws) ||
    !identical(dim(standard_normal_draws), dim(structured))
  ) {
    stop(
      "`standard_normal_draws` must be a matrix with the same ",
      "dimensions as `structured`."
    )
  }
  
  unstructured <- vapply(
    seq_len(ncol(structured)),
    function(k) {
      z <- standard_normal_draws[, k]
      structured_component <- structured[, k]
      
      # Remove the intercept.
      z <- z - mean(z)
      
      # Make the unstructured component orthogonal to the structured
      # component. This gives alpha an exact variance interpretation.
      z <- z -
        sum(z * structured_component) /
        sum(structured_component^2) *
        structured_component
      
      standardise_population(
        z,
        label = paste0("unstructured field ", k)
      )
    },
    numeric(nrow(structured))
  )
  
  dimnames(unstructured) <- dimnames(structured)
  
  unstructured
}


generate_core_truth <- function(
    graph,
    t_grid = seq(0, 1, length.out = 50L),
    truth_structure = c("coordinate", "cluster"),
    alpha = 0.9,
    standard_normal_draws) {
  
  truth_structure <- match.arg(truth_structure)
  
  if (
    !is.numeric(alpha) ||
    length(alpha) != 1L ||
    is.na(alpha) ||
    alpha < 0 ||
    alpha > 1
  ) {
    stop("`alpha` must be a single number in [0, 1].")
  }
  
  fields <- make_structured_fields(
    graph = graph,
    truth_structure = truth_structure
  )
  
  unstructured <- make_unstructured_fields(
    structured = fields$structured,
    standard_normal_draws = standard_normal_draws
  )
  
  coefficient_mixture <-
    sqrt(alpha) * fields$structured +
    sqrt(1 - alpha) * unstructured
  
  coefficients <- sweep(
    coefficient_mixture,
    MARGIN = 2,
    STATS = fields$coefficient_scales,
    FUN = "*"
  )
  
  temporal_basis <- make_temporal_basis(t_grid)
  
  mean_curve <-
    1 +
    0.5 * sin(2 * pi * t_grid) +
    0.25 * cos(4 * pi * t_grid)
  
  node_deviations <- coefficients %*% t(temporal_basis)
  
  true_curves <- sweep(
    node_deviations,
    MARGIN = 2,
    STATS = mean_curve,
    FUN = "+"
  )
  
  node_names <- igraph::V(graph)$name
  
  rownames(coefficients) <- node_names
  rownames(node_deviations) <- node_names
  rownames(true_curves) <- node_names
  
  list(
    truth = true_curves,
    mean_curve = mean_curve,
    deviations = node_deviations,
    coefficients = coefficients,
    structured_coefficients = fields$structured,
    unstructured_coefficients = unstructured,
    coefficient_scales = fields$coefficient_scales,
    temporal_basis = temporal_basis,
    t_grid = t_grid,
    graph = graph,
    truth_structure = truth_structure,
    alpha = alpha,
    cluster = fields$cluster,
    coordinates = fields$coordinates
  )
}

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