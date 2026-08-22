# Construct lattice graphs used in the simulation study.
#
# The graph nodes are ordered column-wise according to expand.grid().
# Node names are always "1", ..., "N", matching the curve names expected
# by netf_smooth().

make_lattice_graph <- function(
    n_side = 5L,
    neighbourhood = c("rook", "queen")) {
  
  neighbourhood <- match.arg(neighbourhood)
  
  if (
    !is.numeric(n_side) ||
    length(n_side) != 1L ||
    is.na(n_side) ||
    n_side < 2 ||
    n_side %% 1 != 0
  ) {
    stop("`n_side` must be a single integer greater than or equal to 2.")
  }
  
  n_side <- as.integer(n_side)
  
  coordinates <- expand.grid(
    row = seq_len(n_side),
    column = seq_len(n_side)
  )
  
  n_nodes <- nrow(coordinates)
  
  row_distance <- abs(
    outer(coordinates$row, coordinates$row, "-")
  )
  
  column_distance <- abs(
    outer(coordinates$column, coordinates$column, "-")
  )
  
  adjacency <- switch(
    neighbourhood,
    
    rook =
      row_distance + column_distance == 1,
    
    queen =
      row_distance <= 1 &
      column_distance <= 1 &
      row_distance + column_distance > 0
  )
  
  diag(adjacency) <- FALSE
  storage.mode(adjacency) <- "numeric"
  
  node_names <- as.character(seq_len(n_nodes))
  
  dimnames(adjacency) <- list(
    node_names,
    node_names
  )
  
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    diag = FALSE
  )
  
  graph <- igraph::set_vertex_attr(
    graph,
    name = "row",
    value = coordinates$row
  )
  
  graph <- igraph::set_vertex_attr(
    graph,
    name = "column",
    value = coordinates$column
  )
  
  graph
}

# Construct one fixed degree-preserving false graph for the rewired-graph
# control. Double-edge swaps preserve node IDs, degree sequence, edge count,
# and lattice-coordinate vertex attributes.

make_degree_preserving_rewire <- function(
    graph,
    seed = 20260819L,
    n_accepted_swaps = 500L,
    min_replaced_edge_proportion = 0.75) {
  
  if (!inherits(graph, "igraph") || igraph::is_directed(graph)) {
    stop("`graph` must be an undirected igraph object.")
  }
  
  if (
    !is.numeric(seed) ||
    length(seed) != 1L ||
    is.na(seed) ||
    seed != as.integer(seed)
  ) {
    stop("`seed` must be one integer.")
  }
  
  if (
    !is.numeric(n_accepted_swaps) ||
    length(n_accepted_swaps) != 1L ||
    is.na(n_accepted_swaps) ||
    n_accepted_swaps < 1L ||
    n_accepted_swaps != as.integer(n_accepted_swaps)
  ) {
    stop("`n_accepted_swaps` must be one positive integer.")
  }
  
  if (
    !is.numeric(min_replaced_edge_proportion) ||
    length(min_replaced_edge_proportion) != 1L ||
    is.na(min_replaced_edge_proportion) ||
    min_replaced_edge_proportion <= 0 ||
    min_replaced_edge_proportion >= 1
  ) {
    stop(
      "`min_replaced_edge_proportion` must lie strictly between zero and one."
    )
  }
  
  node_names <- igraph::V(graph)$name
  
  coordinates <- data.frame(
    row = igraph::vertex_attr(graph, "row"),
    column = igraph::vertex_attr(graph, "column")
  )
  
  if (
    is.null(node_names) ||
    anyNA(coordinates$row) ||
    anyNA(coordinates$column)
  ) {
    stop(
      "`graph` must have node names and `row`/`column` vertex attributes."
    )
  }
  
  original_adjacency <- as.matrix(
    igraph::as_adjacency_matrix(graph, sparse = FALSE)
  )
  
  storage.mode(original_adjacency) <- "integer"
  
  dimnames(original_adjacency) <- list(
    node_names,
    node_names
  )
  
  adjacency <- original_adjacency
  
  original_degree <- rowSums(original_adjacency)
  
  original_edge_count <- sum(
    original_adjacency[upper.tri(original_adjacency)]
  )
  
  # Preserve the caller's RNG state.
  old_kind <- RNGkind()
  
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  
  on.exit(
    {
      do.call(RNGkind, as.list(old_kind))
      
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (
        exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      ) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  
  set.seed(as.integer(seed))
  
  accepted_swaps <- 0L
  attempts <- 0L
  
  maximum_attempts <- as.integer(
    100L * n_accepted_swaps
  )
  
  while (
    accepted_swaps < n_accepted_swaps &&
    attempts < maximum_attempts
  ) {
    attempts <- attempts + 1L
    
    upper_coordinates <- which(
      upper.tri(adjacency),
      arr.ind = TRUE
    )
    
    edge_positions <- which(
      adjacency[upper.tri(adjacency)] == 1L
    )
    
    selected <- upper_coordinates[
      sample(edge_positions, size = 2L, replace = FALSE),
      ,
      drop = FALSE
    ]
    
    a <- selected[1L, 1L]
    b <- selected[1L, 2L]
    
    c <- selected[2L, 1L]
    d <- selected[2L, 2L]
    
    # Two valid degree-preserving double-edge swaps are possible.
    proposal <- if (stats::runif(1L) < 0.5) {
      c(a, d, c, b)
    } else {
      c(a, c, b, d)
    }
    
    u1 <- proposal[1L]
    v1 <- proposal[2L]
    
    u2 <- proposal[3L]
    v2 <- proposal[4L]
    
    # Reject loops, repeated vertices, and already-existing edges.
    if (
      length(unique(c(a, b, c, d))) < 4L ||
      u1 == v1 ||
      u2 == v2 ||
      adjacency[u1, v1] == 1L ||
      adjacency[u2, v2] == 1L
    ) {
      next
    }
    
    candidate <- adjacency
    
    candidate[a, b] <- candidate[b, a] <- 0L
    candidate[c, d] <- candidate[d, c] <- 0L
    
    candidate[u1, v1] <- candidate[v1, u1] <- 1L
    candidate[u2, v2] <- candidate[v2, u2] <- 1L
    
    candidate_graph <- igraph::graph_from_adjacency_matrix(
      candidate,
      mode = "undirected",
      diag = FALSE
    )
    
    # Retain a connected false graph.
    if (!igraph::is_connected(candidate_graph)) {
      next
    }
    
    adjacency <- candidate
    accepted_swaps <- accepted_swaps + 1L
  }
  
  retained_edge_count <- sum(
    (adjacency == 1L) &
      (original_adjacency == 1L)
  ) / 2L
  
  replaced_edge_proportion <- 1 -
    retained_edge_count / original_edge_count
  
  if (
    accepted_swaps < n_accepted_swaps ||
    replaced_edge_proportion < min_replaced_edge_proportion
  ) {
    stop(
      "Could not construct a sufficiently rewired connected graph."
    )
  }
  
  rewired <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    diag = FALSE
  )
  
  rewired <- igraph::set_vertex_attr(
    rewired,
    "name",
    value = node_names
  )
  
  rewired <- igraph::set_vertex_attr(
    rewired,
    "row",
    value = coordinates$row
  )
  
  rewired <- igraph::set_vertex_attr(
    rewired,
    "column",
    value = coordinates$column
  )
  
  if (
    !identical(
      as.integer(igraph::degree(rewired)),
      as.integer(original_degree)
    ) ||
    !igraph::is_connected(rewired)
  ) {
    stop(
      "Internal error: the rewired graph did not preserve graph invariants."
    )
  }
  
  attr(rewired, "rewire_audit") <- list(
    seed = as.integer(seed),
    accepted_swaps = accepted_swaps,
    attempts = attempts,
    original_edge_count = original_edge_count,
    retained_edge_count = retained_edge_count,
    replaced_edge_proportion = replaced_edge_proportion
  )
  
  rewired
}