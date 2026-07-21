# R/netf_smooth.R

#' Network-weighted smoothing of functional data
#'
#' Estimates one smooth function for each node of a graph. The model consists
#' of a functional intercept representing the overall mean curve and
#' node-specific deviations that may vary over the functional domain.
#' A Markov random field penalty encourages curves at neighbouring nodes
#' to have similar shapes.
#'
#' If both the curves and graph nodes are named, curves are matched to nodes
#' by name and reordered to follow the graph-node order. An error is raised
#' if the names do not agree. Unnamed curves are matched to graph nodes
#' positionally.
#'
#' @param curves A [tf::tfd()] or [tf::tfb()] vector containing one function
#'   per graph node. All functions must share a common argument grid.
#' @param graph A graph representation supported by [graph_to_nb()].
#' @param sandwich Covariance estimator passed to [refund::pffr()].
#'   The default is `"none"` to make the result independent of changes in
#'   the default used by different versions of \pkg{refund}.
#' @param bs.int A list specifying the basis for the functional intercept;
#'   see [refund::pffr()]. If `NULL`, a P-spline basis is constructed
#'   automatically.
#' @param bs.yindex A list specifying the basis over the functional domain
#'   for the node-specific deviations; see [refund::pffr()]. If `NULL`,
#'   a P-spline basis with a grid-dependent basis dimension is used.
#' @param ... Additional arguments passed to [refund::pffr()].
#'
#' @return An object of class `netf_fit`.
#'
#' @details
#' The fitted model can be written as
#' \deqn{Y_i(t) = \mu(t) + g_i(t) + \epsilon_i(t),}
#' where \eqn{\mu(t)} is the overall mean function and \eqn{g_i(t)} is the
#' node-specific deviation. The deviations are smoothed jointly over the
#' graph and over the functional domain.
#'
#' The basis dimensions specified through `bs.int` and `bs.yindex` determine
#' the maximum complexity available to the mean function and the
#' node-specific deviations, respectively. Smoothness is additionally
#' controlled by smoothing parameters estimated during model fitting.
#'
#' @export
netf_smooth <- function(curves, graph, ...) {
  UseMethod("netf_smooth")
}

#' @export
netf_smooth.tfd <- function(curves, graph, ...) {
  fit_netf_smooth(curves = curves, graph = graph, ...)
}

#' @export
netf_smooth.tfb <- function(curves, graph, ...) {
  curve_names <- names(curves)
  curves <- tf::tfd(curves)
  names(curves) <- curve_names
  
  fit_netf_smooth(curves = curves, graph = graph, ...)
}

#' @export
netf_smooth.default <- function(curves, graph, ...) {
  cli::cli_abort(
    "No {.fn netf_smooth} method for objects of class {.cls {class(curves)}}."
  )
}

# refund::pffr() rebuilds the model formula internally via deparsing and
# re-parsing, which discards whatever environment the caller attached to
# it, so the `nb` argument can't be resolved through ordinary lexical
# scoping (a formula built with `env = environment()` fails with
# "object 'nb_list' not found"). Anything passed through `xt =` must
# therefore be reachable through an absolute name on the search path.
# Each fit registers its neighbour list under a fresh, unique key in a
# package-level registry and the formula retrieves it via
# asNamespace("netfunsmooth"), which resolves independent of environment
# (and, unlike `:::`, does not trigger an R CMD check NOTE). The key is
# unique per call so concurrent or nested fits never share (and can't
# clobber) each other's neighbour list, and the entry is removed once
# the fit completes. 
.nb_registry <- new.env(parent = emptyenv())

register_nb_list <- function(nb_list) {
  key <- basename(tempfile("nb"))
  assign(key, nb_list, envir = .nb_registry)
  key
}


build_pffr_formula <- function(key, k_node) {
  stats::as.formula(sprintf(
    'Y ~ s(node, bs = "mrf",
           xt = list(nb = get("%s", envir = get(".nb_registry", envir = asNamespace("netfunsmooth")))),
           k = %d)',
    key, k_node
  ))
}

fit_netf_smooth <- function(curves, graph,
                            sandwich = "none",
                            bs.int = NULL,
                            bs.yindex = NULL,
                            ...) {
  nb_list <- graph_to_nb(graph)
  n_nodes <- length(curves)

  if (length(nb_list) != n_nodes) {
    cli::cli_abort(
      "{.arg curves} must contain one curve for each graph node."
    )
  }

  if (is.null(names(nb_list))) {
    names(nb_list) <- as.character(seq_len(n_nodes))
  }
  node_names <- names(nb_list)

  # When both the curves and the graph nodes are named, match by name
  # instead of relying on positional order (curve i ~ node i), which is
  # too brittle for real data. Reorder the curves to the node order and
  # stop if the two name sets disagree.
  curve_names <- names(curves)
  if (!is.null(curve_names)) {
    if (anyDuplicated(curve_names) > 0) {
      cli::cli_abort("{.arg curves} must not contain duplicated names.")
    }
    if (!setequal(curve_names, node_names)) {
      unmatched_nodes <- setdiff(node_names, curve_names)
      unmatched_curves <- setdiff(curve_names, node_names)
      cli::cli_abort(c(
        "Names of {.arg curves} must match the graph node names.",
        "x" = if (length(unmatched_nodes) > 0) {
          "Node{?s} without a matching curve: {.val {unmatched_nodes}}."
        },
        "x" = if (length(unmatched_curves) > 0) {
          "Curve{?s} without a matching node: {.val {unmatched_curves}}."
        }
      ))
    }
    curves <- curves[node_names]
  }

  t_grid <- tf::tf_arg(curves)
  if (is.list(t_grid)) {
    cli::cli_abort(
      "{.arg curves} must be observed on a common grid; consider
       {.fn tf::tf_interpolate} first."
    )
  }
  n_t <- length(unique(t_grid))
  y_mat <- do.call(rbind, tf::tf_evaluate(curves, arg = t_grid))

  df_wide <- data.frame(
    node = factor(node_names, levels = node_names)
  )
  df_wide$Y <- I(y_mat)

  # The two basis specifications control the representational complexity along
  # the functional domain. `bs.int` applies to the overall mean function,
  # whereas `bs.yindex` applies to the node-specific deviations. The basis
  # dimensions are capped by the number of distinct grid points.
  #
  # These dimensions set an upper bound on model complexity; the effective
  # smoothness is additionally determined by the estimated smoothing
  # parameters. A larger default than in pffr() is used for `bs.yindex` to
  # accommodate moderately complex node-specific shape differences.
  if (is.null(bs.int)) {
    bs.int <- list(bs = "ps", k = min(20L, n_t), m = c(2, 1))
  }
  if (is.null(bs.yindex)) {
    bs.yindex <- list(bs = "ps", k = min(15L, n_t), m = c(2, 1))
  }

  key <- register_nb_list(nb_list)
  on.exit(rm(list = key, envir = .nb_registry), add = TRUE)

  form <- build_pffr_formula(key, k_node = n_nodes)

  fit <- refund::pffr(
    formula = form,
    yind = t_grid,
    data = df_wide,
    sandwich = sandwich,
    bs.int = bs.int,
    bs.yindex = bs.yindex,
    ...
  )

  new_netf_fit(
    model = fit,
    curves = curves,
    graph = graph,
    nb_list = nb_list,
    t_grid = t_grid,
    call = match.call()
  )
}
