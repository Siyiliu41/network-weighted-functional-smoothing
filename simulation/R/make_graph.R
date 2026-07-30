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
