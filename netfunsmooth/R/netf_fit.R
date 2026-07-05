###############################################################
# netf_fit.R
#
# PURPOSE
# -------------------------------------------------------------
# EN:
# Defines the fitted model object returned by netf_smooth().
#
###############################################################


new_netf_fit <- function(
    model,
    curves,
    graph,
    nb_list,
    t_grid = NULL,
    call = NULL
) {
  structure(
    list(
      model = model,
      curves = curves,
      graph = graph,
      nb_list = nb_list,
      t_grid = t_grid,
      call = call
    ),
    class = "netf_fit"
  )
}


#' @export
print.netf_fit <- function(x, ...) {
  cat("<netf_fit>\n")
  cat("Network-weighted functional smoothing fit\n")
  cat("\n")
  cat("Number of curves:", length(x$curves), "\n")
  cat("Number of graph nodes:", length(x$nb_list), "\n")

  invisible(x)
}
