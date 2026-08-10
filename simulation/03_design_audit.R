# Deterministic design audit for the core simulation.
#
# This script does not fit any smoothing model. It reconstructs the 50 pilot
# truth draws from their saved L'Ecuyer-CMRG streams and checks:
#
# 1. the aggregate relative unpenalised L2 projection error of the true
#    node-specific deviations for the fitted deviation basis;
# 2. the overall and deviation signal-to-noise ratios on the frozen
#    50-point observation grid;
# 3. the expected 6.25-fold SNR separation between sigma = 0.2 and 0.5;
# 4. pilot checkpoint compatibility and code fingerprints.
#
# Run from the repository root with:
# source("simulation/03_design_audit.R")


# -------------------------------------------------------------------------
# Preconditions and source files
# -------------------------------------------------------------------------

required_packages <- c("igraph", "mgcv")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before running the audit: ",
    paste(missing_packages, collapse = ", ")
  )
}

simulation_sources <- c(
  "simulation/R/make_graph.R",
  "simulation/R/generate_truth.R",
  "simulation/R/rng.R",
  "simulation/R/core_design.R"
)

required_input_files <- c(
  "Simulation_study_specification.md",
  "simulation/02_pilot_simulation.R",
  simulation_sources
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Run this script from the repository root. Missing file(s): ",
    paste(missing_input_files, collapse = ", ")
  )
}

invisible(lapply(simulation_sources, source))


# -------------------------------------------------------------------------
# Frozen audit configuration
# -------------------------------------------------------------------------

audit_config <- list(
  audit_version = "2026-08-10-v2",
  pilot_result_file =
    "simulation/results/core/pilot/pilot_result.rds",
  checkpoint_directory =
    "simulation/results/core/pilot/checkpoints",
  output_directory =
    "simulation/results/core/design_validation",
  expected_pilot_version = "2026-08-09-v2",
  expected_base_seed = 20260805L,
  n_replications = 50L,
  n_side = 5L,
  n_time = 50L,
  n_components = 3L,
  truth_structures = c("coordinate", "cluster"),
  alpha_values = c(0.5, 0.9),
  sigma_values = c(0.2, 0.5),
  k_deviation = 15L,
  quadrature_points = 4001L,
  rpe_threshold = 0.02,
  rpe_preferred_threshold = 0.01,
  expected_snr_ratio = (0.5^2) / (0.2^2),
  numerical_tolerance = 1e-10
)


# -------------------------------------------------------------------------
# Numerical helpers
# -------------------------------------------------------------------------

trapezoidal_weights <- function(t_grid) {
  if (
    !is.numeric(t_grid) ||
    length(t_grid) < 2L ||
    anyNA(t_grid) ||
    any(!is.finite(t_grid)) ||
    any(diff(t_grid) <= 0)
  ) {
    stop("`t_grid` must be a finite, strictly increasing numeric vector.")
  }
  
  intervals <- diff(t_grid)
  
  c(
    intervals[1L] / 2,
    (intervals[-length(intervals)] + intervals[-1L]) / 2,
    intervals[length(intervals)] / 2
  )
}


md5sum_with_lf_line_endings <- function(paths) {
  hashes <- vapply(
    paths,
    function(path) {
      file_size <- file.info(path)$size
      bytes <- readBin(path, what = "raw", n = file_size)
      
      if (length(bytes) >= 2L) {
        cr_positions <- which(bytes == as.raw(13L))
        cr_positions <- cr_positions[
          cr_positions < length(bytes)
        ]
        crlf_positions <- cr_positions[
          bytes[cr_positions + 1L] == as.raw(10L)
        ]
        
        if (length(crlf_positions) > 0L) {
          bytes <- bytes[-crlf_positions]
        }
      }
      
      temporary_file <- tempfile("design-audit-md5-")
      on.exit(unlink(temporary_file), add = TRUE)
      
      connection <- file(temporary_file, open = "wb")
      writeBin(bytes, connection)
      close(connection)
      
      unname(as.character(tools::md5sum(temporary_file)))
    },
    character(1)
  )
  
  names(hashes) <- paths
  hashes
}


make_deviation_basis <- function(t_grid, k_deviation) {
  smooth_specification <- mgcv::s(
    t,
    bs = "ps",
    k = k_deviation,
    m = c(2, 1)
  )
  
  smooth_object <- mgcv::smoothCon(
    smooth_specification,
    data = data.frame(t = t_grid),
    knots = NULL,
    absorb.cons = FALSE
  )[[1L]]
  
  basis <- smooth_object$X
  
  if (
    !is.matrix(basis) ||
    nrow(basis) != length(t_grid) ||
    ncol(basis) != k_deviation ||
    any(!is.finite(basis))
  ) {
    stop("Unexpected P-spline basis returned by `mgcv::smoothCon()`.")
  }
  
  basis
}


calculate_deviation_rpe <- function(deviations, basis, weights) {
  if (
    !is.matrix(deviations) ||
    !is.matrix(basis) ||
    ncol(deviations) != nrow(basis) ||
    length(weights) != nrow(basis)
  ) {
    stop("Incompatible deviations, basis, or quadrature weights.")
  }
  
  if (
    any(!is.finite(deviations)) ||
    any(!is.finite(basis)) ||
    any(!is.finite(weights)) ||
    any(weights <= 0)
  ) {
    stop("Projection inputs must contain finite values and positive weights.")
  }
  
  square_root_weights <- sqrt(weights)
  
  weighted_basis <- basis * square_root_weights
  weighted_deviations <- t(deviations) * square_root_weights
  
  basis_qr <- qr(weighted_basis, LAPACK = FALSE)
  
  if (basis_qr$rank != ncol(basis)) {
    stop("The deviation basis is not full column rank.")
  }
  
  coefficients <- qr.coef(
    basis_qr,
    weighted_deviations
  )
  
  projected_deviations <- t(
    basis %*% coefficients
  )
  
  squared_error <- (deviations - projected_deviations)^2
  
  numerator <- sum(
    sweep(
      squared_error,
      MARGIN = 2,
      STATS = weights,
      FUN = "*"
    )
  )
  
  denominator <- sum(
    sweep(
      deviations^2,
      MARGIN = 2,
      STATS = weights,
      FUN = "*"
    )
  )
  
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The integrated squared deviation norm must be positive.")
  }
  
  list(
    rpe_dev = sqrt(numerator / denominator),
    integrated_squared_error = numerator,
    integrated_squared_deviation = denominator,
    basis_rank = basis_qr$rank
  )
}


calculate_signal_measures <- function(truth_object, sigma) {
  true_curves <- truth_object$truth
  deviations <- truth_object$deviations
  
  if (
    !is.matrix(true_curves) ||
    !is.matrix(deviations) ||
    !identical(dim(true_curves), dim(deviations)) ||
    any(!is.finite(true_curves)) ||
    any(!is.finite(deviations))
  ) {
    stop("The truth object contains invalid truth or deviation matrices.")
  }
  
  if (
    !is.numeric(sigma) ||
    length(sigma) != 1L ||
    is.na(sigma) ||
    !is.finite(sigma) ||
    sigma <= 0
  ) {
    stop("`sigma` must be one finite positive number.")
  }
  
  grand_mean <- mean(true_curves)
  overall_signal_ms <- mean((true_curves - grand_mean)^2)
  deviation_signal_ms <- mean(deviations^2)
  
  data.frame(
    overall_signal_ms = overall_signal_ms,
    deviation_signal_ms = deviation_signal_ms,
    overall_signal_rms = sqrt(overall_signal_ms),
    deviation_signal_rms = sqrt(deviation_signal_ms),
    snr_overall = overall_signal_ms / sigma^2,
    snr_dev = deviation_signal_ms / sigma^2,
    stringsAsFactors = FALSE
  )
}


summarise_values <- function(data, group_columns, value_column) {
  split_key <- interaction(
    data[group_columns],
    drop = TRUE,
    lex.order = TRUE
  )
  
  groups <- split(data, split_key)
  
  rows <- lapply(
    groups,
    function(group) {
      values <- group[[value_column]]
      
      data.frame(
        group[1L, group_columns, drop = FALSE],
        n = length(values),
        minimum = min(values),
        median = stats::median(values),
        mean = mean(values),
        maximum = max(values),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
  
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


get_git_commit <- function() {
  commit <- tryCatch(
    system2(
      "git",
      c("rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(error) character()
  )
  
  if (length(commit) == 1L && nzchar(commit)) {
    commit
  } else {
    NA_character_
  }
}


# -------------------------------------------------------------------------
# Validate the completed pilot and its checkpoints
# -------------------------------------------------------------------------

if (!file.exists(audit_config$pilot_result_file)) {
  stop(
    "Pilot result not found: ",
    audit_config$pilot_result_file
  )
}

pilot_result <- readRDS(audit_config$pilot_result_file)
pilot_config <- pilot_result$config

expected_checkpoint_files <- file.path(
  audit_config$checkpoint_directory,
  sprintf(
    "replication_%03d.rds",
    seq_len(audit_config$n_replications)
  )
)

missing_checkpoint_files <- expected_checkpoint_files[
  !file.exists(expected_checkpoint_files)
]

if (length(missing_checkpoint_files) > 0L) {
  stop(
    "Missing pilot checkpoint(s): ",
    paste(basename(missing_checkpoint_files), collapse = ", ")
  )
}

pilot_config_checks <- c(
  pilot_version = identical(
    pilot_config$pilot_version,
    audit_config$expected_pilot_version
  ),
  base_seed = identical(
    pilot_config$base_seed,
    audit_config$expected_base_seed
  ),
  replication_count = identical(
    pilot_config$n_replications,
    audit_config$n_replications
  ),
  deviation_basis_dimension = identical(
    pilot_config$k_deviation,
    audit_config$k_deviation
  )
)

if (!all(pilot_config_checks)) {
  stop(
    "The saved pilot configuration does not match the frozen audit ",
    "configuration. Failed check(s): ",
    paste(names(pilot_config_checks)[!pilot_config_checks], collapse = ", ")
  )
}

saved_pilot_code_md5 <- pilot_config$code_md5

if (
  is.null(saved_pilot_code_md5) ||
  is.null(names(saved_pilot_code_md5)) ||
  any(!file.exists(names(saved_pilot_code_md5)))
) {
  stop("The saved pilot code fingerprint cannot be verified.")
}

current_pilot_code_md5 <- tools::md5sum(
  names(saved_pilot_code_md5)
)

current_pilot_code_md5_lf <- md5sum_with_lf_line_endings(
  names(saved_pilot_code_md5)
)

saved_pilot_code_md5_character <-
  as.character(saved_pilot_code_md5)

raw_fingerprint_matches <-
  as.character(current_pilot_code_md5) ==
  saved_pilot_code_md5_character

lf_fingerprint_matches <-
  as.character(current_pilot_code_md5_lf) ==
  saved_pilot_code_md5_character

file_fingerprint_matches <-
  raw_fingerprint_matches | lf_fingerprint_matches

pilot_code_matches <- all(file_fingerprint_matches)

line_ending_normalized_files <- names(saved_pilot_code_md5)[
  !raw_fingerprint_matches & lf_fingerprint_matches
]

if (!pilot_code_matches) {
  changed_files <- names(saved_pilot_code_md5)[
    !file_fingerprint_matches
  ]
  
  stop(
    "Current code does not match the code used for the pilot. ",
    "Changed file(s): ",
    paste(changed_files, collapse = ", ")
  )
}


# -------------------------------------------------------------------------
# Construct the common projection basis and audit the 50 truth draws
# -------------------------------------------------------------------------

observation_grid <- seq(
  0,
  1,
  length.out = audit_config$n_time
)

quadrature_grid <- seq(
  0,
  1,
  length.out = audit_config$quadrature_points
)

quadrature_weights <- trapezoidal_weights(
  quadrature_grid
)

deviation_basis <- make_deviation_basis(
  t_grid = quadrature_grid,
  k_deviation = audit_config$k_deviation
)

truth_graph <- make_lattice_graph(
  n_side = audit_config$n_side,
  neighbourhood = "rook"
)

rpe_rows <- vector(
  "list",
  audit_config$n_replications *
    length(audit_config$truth_structures) *
    length(audit_config$alpha_values)
)

snr_rows <- vector(
  "list",
  length(rpe_rows) * length(audit_config$sigma_values)
)

rpe_index <- 0L
snr_index <- 0L
checkpoint_config_matches <- logical(
  audit_config$n_replications
)

for (replication_id in seq_len(audit_config$n_replications)) {
  checkpoint <- readRDS(
    expected_checkpoint_files[replication_id]
  )
  
  checkpoint_config_matches[replication_id] <-
    identical(checkpoint$config, pilot_config) &&
    identical(checkpoint$replication_id, replication_id)
  
  if (!checkpoint_config_matches[replication_id]) {
    stop(
      "Incompatible pilot checkpoint: ",
      expected_checkpoint_files[replication_id]
    )
  }
  
  validate_rng_stream(checkpoint$replication_stream)
  
  random_draws <- generate_core_random_draws(
    replication_stream = checkpoint$replication_stream,
    n_nodes = audit_config$n_side^2,
    n_time = audit_config$n_time,
    n_components = audit_config$n_components
  )
  
  for (truth_structure in audit_config$truth_structures) {
    standard_normal_draws <-
      random_draws$coefficient_draws[[truth_structure]]
    
    for (alpha in audit_config$alpha_values) {
      truth_observed_grid <- generate_core_truth(
        graph = truth_graph,
        t_grid = observation_grid,
        truth_structure = truth_structure,
        alpha = alpha,
        standard_normal_draws = standard_normal_draws
      )
      
      truth_quadrature_grid <- generate_core_truth(
        graph = truth_graph,
        t_grid = quadrature_grid,
        truth_structure = truth_structure,
        alpha = alpha,
        standard_normal_draws = standard_normal_draws
      )
      
      if (!isTRUE(all.equal(
        truth_observed_grid$coefficients,
        truth_quadrature_grid$coefficients,
        tolerance = audit_config$numerical_tolerance,
        check.attributes = TRUE
      ))) {
        stop("Truth coefficients changed when only the evaluation grid changed.")
      }
      
      rpe <- calculate_deviation_rpe(
        deviations = truth_quadrature_grid$deviations,
        basis = deviation_basis,
        weights = quadrature_weights
      )
      
      rpe_index <- rpe_index + 1L
      
      rpe_rows[[rpe_index]] <- data.frame(
        replication = replication_id,
        truth_id = sprintf(
          "%s_a%02d",
          truth_structure,
          round(100 * alpha)
        ),
        truth_structure = truth_structure,
        alpha = alpha,
        k_deviation = audit_config$k_deviation,
        quadrature_points = audit_config$quadrature_points,
        rpe_dev = rpe$rpe_dev,
        integrated_squared_error =
          rpe$integrated_squared_error,
        integrated_squared_deviation =
          rpe$integrated_squared_deviation,
        basis_rank = rpe$basis_rank,
        centred_deviation_max_abs = max(
          abs(colMeans(truth_quadrature_grid$deviations))
        ),
        stringsAsFactors = FALSE
      )
      
      for (sigma in audit_config$sigma_values) {
        snr_index <- snr_index + 1L
        
        signal_measures <- calculate_signal_measures(
          truth_object = truth_observed_grid,
          sigma = sigma
        )
        
        snr_rows[[snr_index]] <- data.frame(
          replication = replication_id,
          truth_id = sprintf(
            "%s_a%02d",
            truth_structure,
            round(100 * alpha)
          ),
          truth_structure = truth_structure,
          alpha = alpha,
          sigma = sigma,
          signal_measures,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

rpe_detail <- do.call(rbind, rpe_rows)
snr_detail <- do.call(rbind, snr_rows)

rownames(rpe_detail) <- NULL
rownames(snr_detail) <- NULL


# -------------------------------------------------------------------------
# Summaries and prespecified pass checks
# -------------------------------------------------------------------------

rpe_summary <- summarise_values(
  data = rpe_detail,
  group_columns = c(
    "truth_id",
    "truth_structure",
    "alpha",
    "k_deviation"
  ),
  value_column = "rpe_dev"
)

snr_overall_summary <- summarise_values(
  data = snr_detail,
  group_columns = c(
    "truth_id",
    "truth_structure",
    "alpha",
    "sigma"
  ),
  value_column = "snr_overall"
)

snr_dev_summary <- summarise_values(
  data = snr_detail,
  group_columns = c(
    "truth_id",
    "truth_structure",
    "alpha",
    "sigma"
  ),
  value_column = "snr_dev"
)

snr_low <- snr_detail[
  snr_detail$sigma == min(audit_config$sigma_values),
  c(
    "replication",
    "truth_id",
    "truth_structure",
    "alpha",
    "snr_overall",
    "snr_dev"
  )
]

snr_high <- snr_detail[
  snr_detail$sigma == max(audit_config$sigma_values),
  c(
    "replication",
    "truth_id",
    "truth_structure",
    "alpha",
    "snr_overall",
    "snr_dev"
  )
]

snr_separation <- merge(
  snr_low,
  snr_high,
  by = c(
    "replication",
    "truth_id",
    "truth_structure",
    "alpha"
  ),
  suffixes = c("_sigma_0_2", "_sigma_0_5"),
  sort = TRUE
)

snr_separation$overall_ratio <-
  snr_separation$snr_overall_sigma_0_2 /
  snr_separation$snr_overall_sigma_0_5

snr_separation$deviation_ratio <-
  snr_separation$snr_dev_sigma_0_2 /
  snr_separation$snr_dev_sigma_0_5

snr_separation$overall_ratio_pass <-
  abs(
    snr_separation$overall_ratio -
      audit_config$expected_snr_ratio
  ) <= audit_config$numerical_tolerance

snr_separation$deviation_ratio_pass <-
  abs(
    snr_separation$deviation_ratio -
      audit_config$expected_snr_ratio
  ) <= audit_config$numerical_tolerance

audit_checks <- c(
  pilot_result_passed = isTRUE(pilot_result$pilot_passed),
  pilot_configuration_matches = all(pilot_config_checks),
  pilot_code_fingerprint_matches = pilot_code_matches,
  all_checkpoint_configurations_match =
    all(checkpoint_config_matches),
  expected_rpe_row_count =
    nrow(rpe_detail) == length(rpe_rows),
  expected_snr_row_count =
    nrow(snr_detail) == length(snr_rows),
  deviations_are_centred =
    max(rpe_detail$centred_deviation_max_abs) <=
    audit_config$numerical_tolerance,
  all_rpe_values_finite =
    all(is.finite(rpe_detail$rpe_dev)),
  rpe_dev_below_2_percent =
    max(rpe_detail$rpe_dev) < audit_config$rpe_threshold,
  all_snr_values_finite_and_positive = all(
    is.finite(snr_detail$snr_overall) &
      snr_detail$snr_overall > 0 &
      is.finite(snr_detail$snr_dev) &
      snr_detail$snr_dev > 0
  ),
  overall_snr_ratio_is_6_25 =
    all(snr_separation$overall_ratio_pass),
  deviation_snr_ratio_is_6_25 =
    all(snr_separation$deviation_ratio_pass)
)

rpe_preferred_check <-
  max(rpe_detail$rpe_dev) <=
  audit_config$rpe_preferred_threshold

audit_passed <- all(audit_checks)


# -------------------------------------------------------------------------
# Save an auditable result bundle
# -------------------------------------------------------------------------

dir.create(
  audit_config$output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

audit_files <- unique(
  c(
    "simulation/03_design_audit.R",
    required_input_files,
    names(saved_pilot_code_md5)
  )
)

audit_files <- audit_files[file.exists(audit_files)]
audit_code_md5 <- tools::md5sum(audit_files)

code_fingerprint <- data.frame(
  file = names(audit_code_md5),
  md5 = as.character(audit_code_md5),
  stringsAsFactors = FALSE
)

audit_result <- list(
  config = audit_config,
  audit_passed = audit_passed,
  checks = audit_checks,
  preferred_rpe_below_or_equal_1_percent =
    rpe_preferred_check,
  maximum_rpe_dev = max(rpe_detail$rpe_dev),
  expected_snr_ratio = audit_config$expected_snr_ratio,
  pilot_config = pilot_config,
  pilot_fingerprint_verification = list(
    saved_md5 = saved_pilot_code_md5,
    current_raw_md5 = current_pilot_code_md5,
    current_lf_normalized_md5 = current_pilot_code_md5_lf,
    line_ending_normalized_files = line_ending_normalized_files
  ),
  rpe_detail = rpe_detail,
  rpe_summary = rpe_summary,
  snr_detail = snr_detail,
  snr_overall_summary = snr_overall_summary,
  snr_dev_summary = snr_dev_summary,
  snr_separation = snr_separation,
  git_commit = get_git_commit(),
  code_md5 = audit_code_md5,
  session_info = utils::sessionInfo()
)

output_files <- list(
  result = file.path(
    audit_config$output_directory,
    "design_audit_result.rds"
  ),
  rpe_detail = file.path(
    audit_config$output_directory,
    "basis_rpe_detail.csv"
  ),
  rpe_summary = file.path(
    audit_config$output_directory,
    "basis_rpe_summary.csv"
  ),
  snr_detail = file.path(
    audit_config$output_directory,
    "snr_detail.csv"
  ),
  snr_overall_summary = file.path(
    audit_config$output_directory,
    "snr_overall_summary.csv"
  ),
  snr_dev_summary = file.path(
    audit_config$output_directory,
    "snr_deviation_summary.csv"
  ),
  snr_separation = file.path(
    audit_config$output_directory,
    "snr_separation.csv"
  ),
  code_fingerprint = file.path(
    audit_config$output_directory,
    "code_fingerprint.csv"
  )
)

saveRDS(audit_result, output_files$result)

utils::write.csv(
  rpe_detail,
  output_files$rpe_detail,
  row.names = FALSE
)

utils::write.csv(
  rpe_summary,
  output_files$rpe_summary,
  row.names = FALSE
)

utils::write.csv(
  snr_detail,
  output_files$snr_detail,
  row.names = FALSE
)

utils::write.csv(
  snr_overall_summary,
  output_files$snr_overall_summary,
  row.names = FALSE
)

utils::write.csv(
  snr_dev_summary,
  output_files$snr_dev_summary,
  row.names = FALSE
)

utils::write.csv(
  snr_separation,
  output_files$snr_separation,
  row.names = FALSE
)

utils::write.csv(
  code_fingerprint,
  output_files$code_fingerprint,
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

cat("\nDesign-audit checks:\n")
print(audit_checks)

if (length(line_ending_normalized_files) > 0L) {
  cat(
    paste0(
      "\nPilot fingerprint verified after CRLF-to-LF normalization for:\n",
      paste0("- ", line_ending_normalized_files, collapse = "\n"),
      "\n"
    )
  )
}

cat("\nDeviation-basis projection error by truth configuration:\n")
print(rpe_summary, row.names = FALSE)

cat(
  sprintf(
    paste0(
      "\nMaximum RPE_dev: %.6f (%.3f%%)\n",
      "Frozen threshold: < %.2f (%.0f%%)\n",
      "Preferred <= 1%% target met: %s\n"
    ),
    max(rpe_detail$rpe_dev),
    100 * max(rpe_detail$rpe_dev),
    audit_config$rpe_threshold,
    100 * audit_config$rpe_threshold,
    if (rpe_preferred_check) "TRUE" else "FALSE"
  )
)

cat("\nOverall SNR summary:\n")
print(snr_overall_summary, row.names = FALSE)

cat("\nDeviation SNR summary:\n")
print(snr_dev_summary, row.names = FALSE)

cat(
  sprintf(
    paste0(
      "\nExpected paired SNR ratio (sigma 0.2 / sigma 0.5): %.2f\n",
      "Observed overall ratio range: [%.6f, %.6f]\n",
      "Observed deviation ratio range: [%.6f, %.6f]\n"
    ),
    audit_config$expected_snr_ratio,
    min(snr_separation$overall_ratio),
    max(snr_separation$overall_ratio),
    min(snr_separation$deviation_ratio),
    max(snr_separation$deviation_ratio)
  )
)

cat(
  "\nAudit result saved to:\n",
  output_files$result,
  "\n"
)

if (!audit_passed) {
  stop(
    "The design audit completed and saved its outputs, ",
    "but at least one prespecified check failed."
  )
}

message("Deterministic design audit completed successfully.")