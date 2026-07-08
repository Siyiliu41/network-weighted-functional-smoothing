# R/netf_smooth.R

#' Network-weighted smoothing for functional data
#'
#' Estimates node-specific smooth curves whose shapes are coupled over a
#' graph, via a functional intercept plus a node-by-t Markov random field
#' interaction fitted with [refund::pffr()].
#'
#' If `curves` is named, curves are matched to graph nodes by name and an
#' error is thrown when the names do not agree with the node names.
#' Unnamed curves are matched positionally (curve `i` corresponds to
#' node `i`).
#'
#' @param curves A `tfd` or `tfb` vector of curves, one function per graph node.
#' @param graph A graph representation supported by [graph_to_nb()].
#' @param sandwich Covariance type passed to [refund::pffr()]. Defaults to
#'   `"none"` so results do not depend on the installed refund version.
#' @param bs.int Basis specification (a list, see [refund::pffr()]) for the
#'   functional intercept, i.e. the mean curve over t. Together with
#'   `bs.yindex` this controls smoothing along t, in addition to the
#'   smoothing parameters over the graph.
#' @param bs.yindex Basis specification for the node-specific deviations
#'   over t. Defaults to a larger basis dimension than [refund::pffr()]'s
#'   own default (`k = 5`), which is too low to capture node-specific
#'   shape differences.
#' @param ... Additional arguments passed to [refund::pffr()].
#'
#' @returns An object of class `netf_fit`.
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
  fit_netf_smooth(curves = tf::tfd(curves), graph = graph, ...)
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
# the fit completes. The registry entry must never be called `nb`,
# because `nb` is also a function in mgcv.
.nb_registry <- new.env(parent = emptyenv())

register_nb_list <- function(nb_list) {
  key <- basename(tempfile("nb"))
  assign(key, nb_list, envir = .nb_registry)
  key
}

# The model is Y ~ 1 + s(node, bs = "mrf"): pffr expands this into a
# functional intercept over t plus a node x t tensor interaction
# (s(yindex.vec) + ti(node, yindex.vec, ...)), i.e. node-specific curves
# whose *shapes* vary over the graph. Do not wrap the smooth in c()
# (that would make it constant over t, reducing the model to a common
# shape plus node-specific vertical shifts) and do not remove the
# intercept.
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

  # The t-basis dimensions of the functional intercept (bs.int) and the
  # node-specific deviations (bs.yindex) are a third smoothing lever
  # besides the smoothing parameters over the graph and over t. pffr's
  # own bs.yindex default (k = 5) is far too low to let curve shapes
  # vary across nodes, so we raise it here; explicit user-supplied
  # values take precedence, and defaults are capped by the grid length.
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
