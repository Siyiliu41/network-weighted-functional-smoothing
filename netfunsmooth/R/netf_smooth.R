# R/netf_smooth.R

#' Network-weighted smoothing for functional data
#'
#' @param curves A `tfd` or `tfb` vector of curves, one function per graph node.
#' @param graph A graph representation supported by [graph_to_nb()].
#' @param sandwich Covariance type passed to [refund::pffr()]. Defaults to
#'   `"none"` so results do not depend on the installed refund version.
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

# refund::pffr() rebuilds the model formula internally and discards
# whatever environment the caller attached to it, so the `nb` argument
# can't be resolved through ordinary lexical scoping (a formula built
# with `env = environment()` fails with "object 'nb_list' not found").
# Each fit therefore registers its neighbour list under a fresh, unique
# key in a package-level registry and references it via `:::`, which
# resolves independent of environment. The key is unique per call so
# concurrent or nested fits never share (and can't clobber) each
# other's neighbour list, and the entry is removed once the fit
# completes.
.nb_registry <- new.env(parent = emptyenv())

register_nb_list <- function(nb_list) {
  key <- basename(tempfile("nb"))
  assign(key, nb_list, envir = .nb_registry)
  key
}

build_pffr_formula <- function(key) {
  stats::as.formula(
    sprintf(
      "Y ~ c(s(node, bs = 'mrf', xt = list(nb = netfunsmooth:::.nb_registry$%s)))",
      key
    )
  )
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

  # Default basis dimensions are capped by the grid length; explicit
  # user-supplied values take precedence.
  if (is.null(bs.int)) {
    bs.int <- list(bs = "ps", k = min(20L, n_t), m = c(2, 1))
  }
  if (is.null(bs.yindex)) {
    bs.yindex <- list(bs = "ps", k = min(5L, n_t), m = c(2, 1))
  }

  key <- register_nb_list(nb_list)
  on.exit(rm(list = key, envir = .nb_registry), add = TRUE)

  form <- build_pffr_formula(key)

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
