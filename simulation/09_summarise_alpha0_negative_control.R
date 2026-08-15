# Summarise the completed alpha = 0 negative-control diagnostic.
#
# Run from the repository root after 08_run_alpha0_negative_control.R:
#   source("simulation/09_summarise_alpha0_negative_control.R")
#
# This script only reads the saved control result. It never fits a model or
# changes the random-number-generator state.

alpha0_summary_script_version <- "2026-08-15-alpha0-summary-r1"

find_project_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    is_project_root <-
      dir.exists(file.path(here, "simulation")) &&
      file.exists(file.path(here, "netfunsmooth", "DESCRIPTION"))

    if (is_project_root) {
      return(here)
    }

    parent <- dirname(here)

    if (identical(parent, here)) {
      stop(
        "Could not find the project root. Expected a `simulation` directory ",
        "and `netfunsmooth/DESCRIPTION`."
      )
    }

    here <- parent
  }
}

finite_mean <- function(x) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  mean(x)
}

finite_mcse <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  stats::sd(x) / sqrt(length(x))
}

project_root <- find_project_root()

control_root <- file.path(
  project_root,
  "simulation", "results", "core", "negative_controls", "alpha0"
)

control_result_path <- file.path(
  control_root,
  "alpha0_control_result.rds"
)

analysis_root <- file.path(control_root, "analysis")

table_directory <- file.path(analysis_root, "tables")

if (!file.exists(control_result_path)) {
  stop("Cannot find completed alpha = 0 result: ", control_result_path)
}

control_result <- readRDS(control_result_path)

required_entries <- c(
  "config",
  "quality_checks",
  "replication_summary",
  "completed_at"
)

if (!is.list(control_result) ||
    !all(required_entries %in% names(control_result))) {
  stop("The alpha = 0 result has an incompatible structure.")
}

if (!isTRUE(all(control_result$quality_checks))) {
  stop(
    "The saved alpha = 0 result did not pass all quality checks; ",
    "pause interpretation."
  )
}

required_config_entries <- c(
  "n_replications",
  "truth_structure",
  "alpha",
  "sigma",
  "estimator_graphs"
)

if (!is.list(control_result$config) ||
    !all(required_config_entries %in% names(control_result$config))) {
  stop("The saved alpha = 0 configuration is incomplete.")
}

replication_summary <- control_result$replication_summary

required_columns <- c(
  "cell_id",
  "truth_structure",
  "alpha",
  "neighbourhood",
  "sigma",
  "replication",
  "method",
  "success",
  "warning_count",
  "node_averaged_ise",
  "elapsed_seconds"
)

if (!is.data.frame(replication_summary) ||
    !all(required_columns %in% names(replication_summary))) {
  stop("`replication_summary` is missing required columns.")
}

if (!identical(as.integer(control_result$config$n_replications), 200L)) {
  stop("This summary script is reserved for the frozen B = 200 alpha = 0 result.")
}

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

make_method_summary <- function(data) {
  groups <- split(
    data,
    interaction(data$cell_id, data$method, drop = TRUE, lex.order = TRUE)
  )

  out <- lapply(groups, function(x) {
    successful <- x$success %in% TRUE & is.finite(x$node_averaged_ise)

    data.frame(
      cell_id = x$cell_id[[1L]],
      truth_structure = x$truth_structure[[1L]],
      alpha = x$alpha[[1L]],
      neighbourhood = x$neighbourhood[[1L]],
      sigma = x$sigma[[1L]],
      method = x$method[[1L]],
      attempted = nrow(x),
      successful = sum(successful),
      failure_rate = 1 - mean(successful),
      total_warnings = sum(x$warning_count, na.rm = TRUE),
      mise = finite_mean(x$node_averaged_ise[successful]),
      mise_mcse = finite_mcse(x$node_averaged_ise[successful]),
      mean_elapsed_seconds = finite_mean(x$elapsed_seconds[successful]),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  result[order(result$neighbourhood, result$method), , drop = FALSE]
}

make_network_paired_summary <- function(data) {
  graphs <- c("rook", "queen")

  out <- lapply(graphs, function(graph_name) {
    x <- data[data$neighbourhood == graph_name, , drop = FALSE]
    nodewise <- x[x$method == "nodewise", , drop = FALSE]
    network <- x[x$method == "network", , drop = FALSE]

    paired <- merge(
      nodewise[, c("replication", "node_averaged_ise", "success")],
      network[, c("replication", "node_averaged_ise", "success")],
      by = "replication",
      all = TRUE,
      sort = TRUE,
      suffixes = c("_nodewise", "_network")
    )

    jointly_successful <-
      paired$success_nodewise %in% TRUE &
      paired$success_network %in% TRUE &
      is.finite(paired$node_averaged_ise_nodewise) &
      is.finite(paired$node_averaged_ise_network)

    paired <- paired[jointly_successful, , drop = FALSE]

    delta <- paired$node_averaged_ise_nodewise -
      paired$node_averaged_ise_network

    data.frame(
      graph = graph_name,
      attempted = nrow(nodewise),
      jointly_successful = nrow(paired),
      paired_mise_nodewise = finite_mean(
        paired$node_averaged_ise_nodewise
      ),
      paired_mise_network = finite_mean(
        paired$node_averaged_ise_network
      ),
      mise_ratio_network_vs_nodewise = if (nrow(paired)) {
        sum(paired$node_averaged_ise_network) /
          sum(paired$node_averaged_ise_nodewise)
      } else {
        NA_real_
      },
      relative_improvement_percent = if (nrow(paired)) {
        100 * (
          1 -
            sum(paired$node_averaged_ise_network) /
            sum(paired$node_averaged_ise_nodewise)
        )
      } else {
        NA_real_
      },
      mean_paired_improvement = finite_mean(delta),
      paired_mcse = finite_mcse(delta),
      lower_two_mcse = finite_mean(delta) - 2 * finite_mcse(delta),
      upper_two_mcse = finite_mean(delta) + 2 * finite_mcse(delta),
      win_rate = if (length(delta)) mean(delta > 0) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL

  result[match(c("queen", "rook"), result$graph), , drop = FALSE]
}

check_non_network_pairing <- function(data) {
  non_network <- data[
    data$method %in% c("raw", "pooled", "nodewise"),
    ,
    drop = FALSE
  ]

  groups <- split(
    non_network,
    interaction(non_network$replication, non_network$method, drop = TRUE)
  )

  all(vapply(groups, function(x) {
    x <- x[order(x$neighbourhood), , drop = FALSE]

    nrow(x) == 2L &&
      identical(x$success[[1L]], x$success[[2L]]) &&
      identical(x$node_averaged_ise[[1L]], x$node_averaged_ise[[2L]]) &&
      identical(x$warning_count[[1L]], x$warning_count[[2L]])
  }, logical(1L)))
}

method_summary <- make_method_summary(replication_summary)
paired_summary <- make_network_paired_summary(replication_summary)

analysis_checks <- c(
  saved_control_quality_checks_passed = isTRUE(all(control_result$quality_checks)),
  expected_1600_method_rows = nrow(replication_summary) == 200L * 2L * 4L,
  expected_two_cells = length(unique(replication_summary$cell_id)) == 2L,
  expected_alpha_zero = all(replication_summary$alpha == 0),
  expected_coordinate_label = all(replication_summary$truth_structure == "coordinate"),
  expected_high_noise = all(replication_summary$sigma == 0.5),
  expected_rook_and_queen = identical(
    sort(unique(replication_summary$neighbourhood)),
    c("queen", "rook")
  ),
  no_duplicate_cell_replication_method = !anyDuplicated(
    replication_summary[, c("cell_id", "replication", "method")]
  ),
  all_method_cells_have_200_attempts = all(method_summary$attempted == 200L),
  all_network_pairs_have_200_successes = all(
    paired_summary$jointly_successful == 200L
  ),
  non_network_results_exactly_paired_across_graphs =
    check_non_network_pairing(replication_summary),
  finite_paired_estimates = all(vapply(
    paired_summary[, c(
      "paired_mise_nodewise",
      "paired_mise_network",
      "mise_ratio_network_vs_nodewise",
      "relative_improvement_percent",
      "mean_paired_improvement",
      "paired_mcse",
      "win_rate"
    )],
    function(x) all(is.finite(x)),
    logical(1L)
  ))
)

if (!isTRUE(all(analysis_checks))) {
  stop(
    "Alpha = 0 summary audit failed: ",
    paste(names(analysis_checks)[!analysis_checks], collapse = ", "),
    "."
  )
}

compact_table <- paired_summary[, c(
  "graph",
  "mise_ratio_network_vs_nodewise",
  "relative_improvement_percent",
  "win_rate"
)]

names(compact_table) <- c(
  "graph",
  "mise_ratio",
  "relative_improvement_percent",
  "win_rate"
)

utils::write.csv(
  method_summary,
  file.path(table_directory, "alpha0_method_performance.csv"),
  row.names = FALSE
)

utils::write.csv(
  paired_summary,
  file.path(table_directory, "alpha0_network_vs_nodewise_detail.csv"),
  row.names = FALSE
)

utils::write.csv(
  compact_table,
  file.path(table_directory, "table_5_alpha0_negative_control.csv"),
  row.names = FALSE
)

format_fixed <- function(x, digits) {
  formatC(x, format = "f", digits = digits)
}

latex_rows <- vapply(seq_len(nrow(compact_table)), function(i) {
  graph_label <- if (compact_table$graph[[i]] == "queen") "Queen" else "Rook"

  paste0(
    graph_label,
    " & ", format_fixed(compact_table$mise_ratio[[i]], 3L),
    " & ", format_fixed(compact_table$relative_improvement_percent[[i]], 2L),
    " & ", format_fixed(100 * compact_table$win_rate[[i]], 1L),
    " \\\\"
  )
}, character(1L))

latex_table <- c(
  "% Generated by simulation/09_summarise_alpha0_negative_control.R",
  "\\begin{tabular}{lrrr}",
  "  \\toprule",
  "  Graph & MISE ratio & Reduction (\\%) & Win rate (\\%) \\\\ ",
  "  \\midrule",
  paste0("  ", latex_rows),
  "  \\bottomrule",
  "\\end{tabular}"
)

latex_table_path <- file.path(
  table_directory,
  "table_5_alpha0_negative_control.tex"
)

writeLines(latex_table, latex_table_path, useBytes = TRUE)

analysis_result <- list(
  analysis_version = alpha0_summary_script_version,
  input_path = control_result_path,
  input_completed_at = control_result$completed_at,
  r_version = R.version.string,
  session_info = utils::sessionInfo(),
  analysis_checks = analysis_checks,
  method_summary = method_summary,
  paired_summary = paired_summary,
  compact_table = compact_table,
  table_paths = list.files(table_directory, full.names = TRUE),
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

saveRDS(
  analysis_result,
  file.path(analysis_root, "alpha0_negative_control_analysis_result.rds"),
  version = 3
)

cat("\nAlpha = 0 negative-control analysis checks:\n")
print(analysis_checks)

cat("\nNetwork versus nodewise, alpha = 0:\n")
print(paired_summary)

cat("\nOutputs written to:\n", analysis_root, "\n", sep = "")