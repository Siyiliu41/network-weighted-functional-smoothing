#' Convert graph objects to neighbour lists
#'
#' Convert supported graph objects to neighbour lists compatible with
#' `mgcv::s(..., bs = "mrf")`.
#'
#' @param graph A graph object.
#' @param id For an `sf` object, the name of a polygon-identifier column or a
#'   vector containing one unique identifier per polygon. Ignored for other
#'   graph classes.
#' @param ... Additional arguments passed to methods.
#'
#' @returns A named list. Each element contains the neighbours of one node.
#'
#' @export
graph_to_nb <- function(graph, id = NULL, ...) {
  UseMethod("graph_to_nb")
}

#' @export
graph_to_nb.igraph <- function(graph, id = NULL, ...) {
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
graph_to_nb.matrix <- function(graph, id = NULL, ...) {
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
graph_to_nb.sf <- function(graph, id = NULL, queen = TRUE, ...) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    cli::cli_abort("Install {.pkg spdep} to use {.fn graph_to_nb} with {.cls sf} objects.")
  }
  
  if (is.null(id)) {
    cli::cli_abort(
      "For an {.cls sf} graph, supply {.arg id} as a polygon-identifier column name or a vector of unique polygon identifiers."
    )
  }
  
  node_names <- if (
    is.character(id) &&
    length(id) == 1L &&
    id %in% names(graph)
  ) {
    as.character(graph[[id]])
  } else {
    as.character(id)
  }
  
  if (
    length(node_names) != nrow(graph) ||
    anyNA(node_names) ||
    any(!nzchar(node_names)) ||
    anyDuplicated(node_names)
  ) {
    cli::cli_abort(
      "{.arg id} must provide one unique, non-missing identifier per polygon."
    )
  }
  
  nb <- spdep::poly2nb(graph, queen = queen, ...)
  
  nb_list <- lapply(nb, function(x) {
    if (length(x) == 0L || (length(x) == 1L && identical(x, 0L))) {
      character()
    } else {
      unname(node_names[as.integer(x)])
    }
  })
  
  names(nb_list) <- node_names
  
  nb_list
}

#' @export
graph_to_nb.default <- function(graph, id = NULL, ...) {
  cli::cli_abort(
    "No {.fn graph_to_nb} method for objects of class {.cls {class(graph)}}."
  )
}