# Baseline methods for the simulation study.

fit_nodewise_smoothing <- function(
    observed,
    t_grid,
    k = 10L,
    keep_models = FALSE) {

  if (!is.matrix(observed)) {
    stop("`observed` must be a numeric matrix.")
  }

  if (
    !is.numeric(observed) ||
    anyNA(observed) ||
    any(!is.finite(observed))
  ) {
    stop("`observed` must contain finite numeric values.")
  }

  if (
    !is.numeric(t_grid) ||
    length(t_grid) != ncol(observed) ||
    anyNA(t_grid) ||
    any(diff(t_grid) <= 0)
  ) {
    stop(
      "`t_grid` must be strictly increasing and have one value ",
      "for every column of `observed`."
    )
  }

  if (
    !is.numeric(k) ||
    length(k) != 1L ||
    is.na(k) ||
    k < 3 ||
    k >= length(t_grid) ||
    k %% 1 != 0
  ) {
    stop(
      "`k` must be an integer greater than or equal to 3 ",
      "and smaller than the number of grid points."
    )
  }

  k <- as.integer(k)

  n_nodes <- nrow(observed)

  predictions <- matrix(
    NA_real_,
    nrow = n_nodes,
    ncol = length(t_grid)
  )

  models <- if (keep_models) {
    vector("list", n_nodes)
  } else {
    NULL
  }

  diagnostics_by_node <- vector(
    "list",
    n_nodes
  )

  node_ids <- rownames(observed)

  if (
    is.null(node_ids) ||
    length(node_ids) != n_nodes ||
    any(!nzchar(node_ids))
  ) {
    node_ids <- as.character(
      seq_len(n_nodes)
    )
  }

  model_data <- data.frame(
    t = t_grid,
    response = numeric(length(t_grid))
  )

  for (node_index in seq_len(n_nodes)) {

    model_data$response <- observed[node_index, ]

    model <- mgcv::gam(
      response ~ s(
        t,
        bs = "ps",
        k = k,
        m = c(2, 1)
      ),
      data = model_data,
      method = "REML"
    )

    diagnostics_by_node[[node_index]] <-
      extract_gam_diagnostics(
        model = model,
        fit_id = paste0(
          "node_",
          node_ids[node_index]
        )
      )

    predictions[node_index, ] <- stats::predict(
      model,
      newdata = data.frame(t = t_grid),
      type = "response"
    )

    if (keep_models) {
      models[[node_index]] <- model
    }
  }

  node_names <- rownames(observed)

  rownames(predictions) <- node_names
  colnames(predictions) <- colnames(observed)

  if (keep_models) {
    names(models) <- node_names
  }

  names(diagnostics_by_node) <- node_ids

  diagnostics <- list(
    fit = do.call(
      rbind,
      lapply(
        diagnostics_by_node,
        `[[`,
        "fit"
      )
    ),
    penalties = do.call(
      rbind,
      lapply(
        diagnostics_by_node,
        `[[`,
        "penalties"
      )
    ),
    terms = do.call(
      rbind,
      lapply(
        diagnostics_by_node,
        `[[`,
        "terms"
      )
    )
  )

  diagnostics <- lapply(
    diagnostics,
    function(diagnostic_table) {
      rownames(diagnostic_table) <- NULL
      diagnostic_table
    }
  )

  list(
    predictions = predictions,
    models = models,
    diagnostics = diagnostics,
    k = k
  )
}

fit_pooled_smoothing <- function(
    observed,
    t_grid,
    k = 10L,
    keep_model = FALSE) {

  if (!is.matrix(observed)) {
    stop("`observed` must be a numeric matrix.")
  }

  if (
    !is.numeric(observed) ||
    anyNA(observed) ||
    any(!is.finite(observed))
  ) {
    stop("`observed` must contain finite numeric values.")
  }

  if (
    !is.numeric(t_grid) ||
    length(t_grid) != ncol(observed) ||
    anyNA(t_grid) ||
    any(diff(t_grid) <= 0)
  ) {
    stop(
      "`t_grid` must be strictly increasing and have one value ",
      "for every column of `observed`."
    )
  }

  if (
    !is.numeric(k) ||
    length(k) != 1L ||
    is.na(k) ||
    k < 3 ||
    k >= length(t_grid) ||
    k %% 1 != 0
  ) {
    stop(
      "`k` must be an integer greater than or equal to 3 ",
      "and smaller than the number of grid points."
    )
  }

  k <- as.integer(k)

  n_nodes <- nrow(observed)

  pooled_data <- data.frame(
    t = rep(t_grid, times = n_nodes),
    response = as.vector(t(observed))
  )

  model <- mgcv::gam(
    response ~ s(
      t,
      bs = "ps",
      k = k,
      m = c(2, 1)
    ),
    data = pooled_data,
    method = "REML"
  )

  diagnostics <- extract_gam_diagnostics(
    model = model,
    fit_id = "pooled"
  )

  pooled_curve <- stats::predict(
    model,
    newdata = data.frame(t = t_grid),
    type = "response"
  )

  predictions <- matrix(
    rep(pooled_curve, each = n_nodes),
    nrow = n_nodes,
    ncol = length(t_grid)
  )

  rownames(predictions) <- rownames(observed)
  colnames(predictions) <- colnames(observed)

  list(
    predictions = predictions,
    model = if (keep_model) model else NULL,
    diagnostics = diagnostics,
    pooled_curve = pooled_curve,
    k = k
  )
}