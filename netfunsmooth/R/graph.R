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

  row_names <- rownames(graph)
  col_names <- colnames(graph)

  if (!is.null(row_names) && !is.null(col_names)) {
    if (anyNA(row_names) ||
        anyNA(col_names) ||
        any(!nzchar(row_names)) ||
        any(!nzchar(col_names))) {
      cli::cli_abort(
        "Row and column names of {.arg graph} must be non-empty."
      )
    }

    if (anyDuplicated(row_names) || anyDuplicated(col_names)) {
      cli::cli_abort(
        "Row and column names of {.arg graph} must be unique."
      )
    }

    if (!setequal(row_names, col_names)) {
      cli::cli_abort(
        "Row and column names of {.arg graph} must contain the same nodes."
      )
    }

    graph <- graph[
      ,
      match(row_names, col_names),
      drop = FALSE
    ]

    colnames(graph) <- row_names
  }
  if (!isSymmetric(graph, check.attributes = FALSE)) {
    cli::cli_abort("{.arg graph} must be symmetric.")
  }

  # Self-neighbours are not included.
  diag(graph) <- 0

  if (!is.null(row_names)) {
    node_names <- row_names
  } else if (!is.null(col_names)) {
    node_names <- col_names
  } else {
    node_names <- as.character(seq_len(nrow(graph)))
  }

  nb_list <- lapply(seq_len(nrow(graph)), function(i) {
    node_names[graph[i, ] != 0]
  })

  names(nb_list) <- node_names

  nb_list
}

#' @export
graph_to_nb.sf <- function(graph, queen = TRUE, ...) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    cli::cli_abort("Install {.pkg spdep} to use {.fn graph_to_nb} with {.cls sf} objects.")
  }

  nb <- spdep::poly2nb(graph, queen = queen, ...)

  nb_list <- lapply(nb, function(x) {
    if (identical(x, 0L)) {
      character()
    } else {
      as.character(x)
    }
  })

  names(nb_list) <- as.character(seq_along(nb_list))

  nb_list
}

#' @export
graph_to_nb.default <- function(graph, ...) {
  cli::cli_abort(
    "No {.fn graph_to_nb} method for objects of class {.cls {class(graph)}}."
  )
}
