#' Convert graph objects to neighbour lists
#'
#' Convert supported graph objects to neighbour lists compatible with
#' `mgcv::s(..., bs = "mrf")`.
#'
#' @param graph A graph object.
#' @param ... Additional arguments passed to methods.
#'
#' @returns A named list. Each element contains the neighbours of one node.
#'
#' @export
graph_to_nb <- function(graph, ...) {
  UseMethod("graph_to_nb")
}

#' @export
graph_to_nb.igraph <- function(graph, ...) {
  node_names <- igraph::V(graph)$name

  if (is.null(node_names)) {
    node_names <- as.character(seq_len(igraph::vcount(graph)))
    igraph::V(graph)$name <- node_names
  }

  nb_list <- igraph::adjacent_vertices(graph, igraph::V(graph))
  nb_list <- lapply(nb_list, names)
  names(nb_list) <- node_names

  nb_list
}

#' @export
graph_to_nb.matrix <- function(graph, ...) {
  checkmate::assert_matrix(graph, mode = "numeric")

  if (nrow(graph) != ncol(graph)) {
    cli::cli_abort("{.arg graph} must be a square adjacency matrix.")
  }

  if (anyNA(graph)) {
    cli::cli_abort("{.arg graph} must not contain missing values.")
  }

  if (!isSymmetric(graph)) {
    cli::cli_abort("{.arg graph} must be symmetric.")
  }

  diag(graph) <- 0

  node_names <- rownames(graph)

  if (is.null(node_names)) {
    node_names <- as.character(seq_len(nrow(graph)))
  }

  nb_list <- lapply(seq_len(nrow(graph)), function(i) {
    node_names[graph[i, ] != 0]
  })

  names(nb_list) <- node_names

  nb_list
}

#' @export
graph_to_nb.default <- function(graph, ...) {
  cli::cli_abort(
    "No {.fn graph_to_nb} method for objects of class {.cls {class(graph)}}."
  )
}
