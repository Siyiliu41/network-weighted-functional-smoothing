# Extract compact diagnostics from an mgcv::gam model before the
# fitted model object is discarded.

extract_gam_diagnostics <- function(
    model,
    fit_id = NA_character_) {

  if (!inherits(model, "gam")) {
    stop("`model` must inherit from class <gam>.")
  }

  if (
    !is.character(fit_id) ||
    length(fit_id) != 1L
  ) {
    stop("`fit_id` must be one character value.")
  }


  # Fit- and optimiser-level diagnostics -------------------------------

  outer_info <- model$outer.info

  optimizer_convergence <- NA_character_
  optimizer_iterations <- NA_integer_
  maximum_absolute_gradient <- NA_real_

  if (!is.null(outer_info)) {
    if (!is.null(outer_info$conv)) {
      optimizer_convergence <- paste(
        outer_info$conv,
        collapse = " | "
      )
    }

    if (!is.null(outer_info$iter)) {
      optimizer_iterations <- as.integer(
        outer_info$iter[1L]
      )
    }

    if (
      !is.null(outer_info$grad) &&
      length(outer_info$grad) > 0L
    ) {
      maximum_absolute_gradient <- max(
        abs(as.numeric(outer_info$grad))
      )
    }
  }

  convergence_available <- !is.null(model$converged)

  fit_diagnostics <- data.frame(
    fit_id = fit_id,
    convergence_available = convergence_available,
    converged = if (convergence_available) {
      isTRUE(model$converged)
    } else {
      NA
    },
    optimizer = if (!is.null(model$optimizer)) {
      paste(model$optimizer, collapse = " / ")
    } else {
      NA_character_
    },
    optimizer_convergence = optimizer_convergence,
    optimizer_iterations = optimizer_iterations,
    maximum_absolute_gradient =
      maximum_absolute_gradient,
    stringsAsFactors = FALSE
  )


  # Penalty-level smoothing parameters ---------------------------------

  smoothing_parameters <- model$sp

  if (
    is.null(smoothing_parameters) ||
    length(smoothing_parameters) == 0L
  ) {
    penalty_diagnostics <- data.frame(
      fit_id = character(),
      penalty_id = integer(),
      parameter = character(),
      smoothing_parameter = numeric(),
      log_smoothing_parameter = numeric(),
      extreme_smoothing_parameter = logical(),
      stringsAsFactors = FALSE
    )
  } else {
    parameter_names <- names(smoothing_parameters)

    if (
      is.null(parameter_names) ||
      any(!nzchar(parameter_names))
    ) {
      parameter_names <- paste0(
        "sp_",
        seq_along(smoothing_parameters)
      )
    }

    log_smoothing_parameters <- log(
      as.numeric(smoothing_parameters)
    )

    penalty_diagnostics <- data.frame(
      fit_id = fit_id,
      penalty_id = seq_along(smoothing_parameters),
      parameter = parameter_names,
      smoothing_parameter =
        as.numeric(smoothing_parameters),
      log_smoothing_parameter =
        log_smoothing_parameters,
      extreme_smoothing_parameter =
        !is.finite(log_smoothing_parameters) |
        log_smoothing_parameters < -15 |
        log_smoothing_parameters > 15,
      stringsAsFactors = FALSE
    )
  }


  # Term-level EDF diagnostics -----------------------------------------

  smooth_terms <- model$smooth

  if (
    is.null(smooth_terms) ||
    length(smooth_terms) == 0L
  ) {
    term_diagnostics <- data.frame(
      fit_id = character(),
      term_id = integer(),
      term = character(),
      edf = numeric(),
      attainable_rank = integer(),
      edf_rank_ratio = numeric(),
      edf_near_rank = logical(),
      stringsAsFactors = FALSE
    )
  } else {
    term_rows <- lapply(
      seq_along(smooth_terms),
      function(term_index) {
        smooth_term <- smooth_terms[[term_index]]

        coefficient_indices <- seq.int(
          from = smooth_term$first.para,
          to = smooth_term$last.para
        )

        term_edf <- sum(
          model$edf[coefficient_indices]
        )

        # Number of coefficients remaining for this smooth after
        # mgcv's identifiability constraints.
        attainable_rank <- length(
          coefficient_indices
        )

        edf_rank_ratio <- term_edf / attainable_rank

        data.frame(
          fit_id = fit_id,
          term_id = term_index,
          term = smooth_term$label,
          edf = term_edf,
          attainable_rank = attainable_rank,
          edf_rank_ratio = edf_rank_ratio,
          edf_near_rank =
            is.finite(edf_rank_ratio) &&
            edf_rank_ratio >= 0.95,
          stringsAsFactors = FALSE
        )
      }
    )

    term_diagnostics <- do.call(
      rbind,
      term_rows
    )

    rownames(term_diagnostics) <- NULL
  }


  list(
    fit = fit_diagnostics,
    penalties = penalty_diagnostics,
    terms = term_diagnostics
  )
}