###############################################################
# predict.R
#
# PURPOSE
# -------------------------------------------------------------
# Defines predict.netf_fit().
#
###############################################################


#' Predict fitted curves from a network smoother
#'
#' @param object A `netf_fit` object.
#' @param newdata Optional new data passed to the underlying model.
#' @param ... Additional arguments passed to prediction methods.
#'
#' @returns A `tfd` object containing the fitted or predicted curves.
#'
#' @export
predict.netf_fit <- function(object, newdata = NULL, ...) {
  if (is.null(object$model)) {
    cli::cli_abort("{.arg object$model} is missing. Cannot compute predictions.")
  }

  pred <- if (is.null(newdata)) {
    stats::fitted(object$model)
  } else {
    stats::predict(object$model, newdata = newdata, ...)
  }

  pred <- as.matrix(pred)

  t_grid <- object$t_grid

  if (is.null(t_grid)) {
    t_grid <- seq_len(ncol(pred))
  }

  tf::tfd(pred, arg = t_grid)
}
