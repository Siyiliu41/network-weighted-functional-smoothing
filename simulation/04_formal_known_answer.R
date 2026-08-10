# Formal prespecified known-answer sanity check (Section 3.8)
#
# Fixed design:
# - coordinate-smooth truth
# - alpha = 1
# - rook estimator graph
# - sigma = 0.5
# - 5 x 5 lattice, N = 25
# - T = 50
# - B_sanity = 50 attempted replications
#
# Run this script from the repository root. Failed fits are recorded and are
# never replaced. Interrupted runs resume from validated replication-level
# checkpoints.

formal_known_answer_script_version <-
  "2026-08-10-no-make-seed-ledger"

message(
  "Formal known-answer script version: ",
  formal_known_answer_script_version
)

devtools::load_all(
  "netfunsmooth",
  quiet = TRUE
)

source("simulation/R/make_graph.R")
source("simulation/R/rng.R")
source("simulation/R/generate_truth.R")
source("simulation/R/simulate_data.R")
source("simulation/R/model_diagnostics.R")
source("simulation/R/fit_baselines.R")
source("simulation/R/metrics.R")
source("simulation/R/run_replication.R")
source("simulation/R/run_core_cell.R")


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

normalised_file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot fingerprint missing file: ", path)
  }
  
  lines <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  temporary_file <- tempfile(
    pattern = "normalised_",
    fileext = ".txt"
  )
  
  on.exit(
    unlink(temporary_file),
    add = TRUE
  )
  
  normalised_text <- paste0(
    paste(lines, collapse = "\n"),
    "\n"
  )
  
  writeBin(
    charToRaw(enc2utf8(normalised_text)),
    temporary_file
  )
  
  unname(tools::md5sum(temporary_file))
}


make_code_fingerprint <- function(paths) {
  paths <- sort(unique(paths))
  hashes <- vapply(
    paths,
    normalised_file_md5,
    character(1L)
  )
  names(hashes) <- paths
  hashes
}


make_corrected_seed_ledger <- function(
    base_seed,
    n_pilot,
    n_sanity) {
  
  if (
    !is.numeric(n_pilot) ||
    length(n_pilot) != 1L ||
    is.na(n_pilot) ||
    n_pilot < 1L ||
    n_pilot != as.integer(n_pilot) ||
    !is.numeric(n_sanity) ||
    length(n_sanity) != 1L ||
    is.na(n_sanity) ||
    n_sanity < 1L ||
    n_sanity != as.integer(n_sanity)
  ) {
    stop("`n_pilot` and `n_sanity` must be positive integers.")
  }
  
  n_pilot <- as.integer(n_pilot)
  n_sanity <- as.integer(n_sanity)
  
  # The completed pilot used make_replication_streams() directly, so it
  # consumed positions 1,...,n_pilot from the base seed. The formal sanity
  # check must therefore use the non-overlapping continuation. This records
  # the correction instead of pretending that the originally prespecified
  # phase labels were used by the pilot.
  all_streams <- make_replication_streams(
    base_seed = base_seed,
    n_replications = n_pilot + n_sanity
  )
  
  pilot_indices <- seq_len(n_pilot)
  sanity_indices <- n_pilot + seq_len(n_sanity)
  
  structure(
    list(
      rng_kind = "L'Ecuyer-CMRG",
      base_seed = as.integer(base_seed),
      prespecified_phase_order = c(
        "sanity",
        "pilot",
        "oracle_pilot",
        "timing_n100",
        "main",
        "n100"
      ),
      actual_allocation_order = c(
        "pilot_existing",
        "sanity_formal"
      ),
      allocation_correction = paste(
        "The completed pilot used stream positions 1 through",
        n_pilot,
        "directly from make_replication_streams().",
        "The formal sanity check therefore uses the non-overlapping",
        "continuation."
      ),
      stream_positions = list(
        pilot_existing = pilot_indices,
        sanity_formal = sanity_indices
      ),
      streams = list(
        pilot = all_streams[pilot_indices],
        sanity = all_streams[sanity_indices]
      )
    ),
    class = "corrected_simulation_seed_ledger"
  )
}


save_rds_atomically <- function(object, path) {
  directory <- dirname(path)
  
  if (!dir.exists(directory)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = directory,
    fileext = ".tmp"
  )
  
  on.exit(
    unlink(temporary_path),
    add = TRUE
  )
  
  saveRDS(
    object,
    file = temporary_path,
    version = 3
  )
  
  if (file.exists(path)) {
    stop(
      "Refusing to overwrite an existing checkpoint: ",
      path
    )
  }
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish checkpoint atomically: ", path)
  }
  
  invisible(path)
}


make_sanity_component_streams <- function(replication_stream) {
  validate_rng_stream(replication_stream)
  
  coordinate_unstructured <- replication_stream
  cluster_unstructured <- parallel::nextRNGSubStream(
    coordinate_unstructured
  )
  observation_error <- parallel::nextRNGSubStream(
    cluster_unstructured
  )
  
  list(
    coordinate_unstructured = coordinate_unstructured,
    cluster_unstructured = cluster_unstructured,
    observation_error = observation_error
  )
}


make_sanity_scenario <- function(
    replication_stream,
    config) {
  
  component_streams <- make_sanity_component_streams(
    replication_stream
  )
  
  graph <- make_lattice_graph(
    n_side = config$n_side,
    neighbourhood = config$neighbourhood
  )
  
  t_grid <- seq(
    0,
    1,
    length.out = config$n_grid
  )
  
  n_nodes <- igraph::vcount(graph)
  
  unstructured_draws <- with_rng_stream(
    component_streams$coordinate_unstructured,
    matrix(
      stats::rnorm(n_nodes * 3L),
      nrow = n_nodes,
      ncol = 3L
    )
  )
  
  standard_normal_errors <- with_rng_stream(
    component_streams$observation_error,
    matrix(
      stats::rnorm(n_nodes * config$n_grid),
      nrow = n_nodes,
      ncol = config$n_grid
    )
  )
  
  truth_object <- generate_core_truth(
    graph = graph,
    t_grid = t_grid,
    truth_structure = config$truth_structure,
    alpha = config$alpha,
    standard_normal_draws = unstructured_draws
  )
  
  simulated_data <- simulate_observed_curves(
    truth_object = truth_object,
    sigma = config$sigma,
    standard_normal_errors = standard_normal_errors
  )
  
  cell <- data.frame(
    cell_id = "sanity_coordinate_a100_rook_s50",
    truth_structure = config$truth_structure,
    alpha = config$alpha,
    neighbourhood = config$neighbourhood,
    sigma = config$sigma,
    stringsAsFactors = FALSE
  )
  
  list(
    scenario = list(
      cell = cell,
      graph = graph,
      truth_object = truth_object,
      simulated_data = simulated_data
    ),
    component_streams = component_streams
  )
}


validate_checkpoint <- function(
    checkpoint,
    replication_id,
    replication_stream,
    config,
    code_fingerprint) {
  
  required_entries <- c(
    "formal_version",
    "config",
    "code_fingerprint",
    "replication_id",
    "replication_stream",
    "component_streams",
    "result"
  )
  
  if (
    !is.list(checkpoint) ||
    !all(required_entries %in% names(checkpoint))
  ) {
    stop(
      "Checkpoint ",
      replication_id,
      " has an incompatible structure."
    )
  }
  
  checks <- c(
    formal_version = identical(
      checkpoint$formal_version,
      config$formal_version
    ),
    config = identical(
      checkpoint$config,
      config
    ),
    code_fingerprint = identical(
      checkpoint$code_fingerprint,
      code_fingerprint
    ),
    replication_id = identical(
      checkpoint$replication_id,
      as.integer(replication_id)
    ),
    replication_stream = identical(
      checkpoint$replication_stream,
      replication_stream
    )
  )
  
  if (!all(checks)) {
    stop(
      "Checkpoint ",
      replication_id,
      " is incompatible: ",
      paste(names(checks)[!checks], collapse = ", "),
      "."
    )
  }
  
  invisible(checkpoint)
}


# -------------------------------------------------------------------------
# Frozen configuration and seed ledger
# -------------------------------------------------------------------------

pilot_checkpoint_directory <-
  "simulation/results/core/pilot/checkpoints"

pilot_checkpoint_path <- file.path(
  pilot_checkpoint_directory,
  "replication_001.rds"
)

if (!file.exists(pilot_checkpoint_path)) {
  stop(
    "The first pilot checkpoint is required to recover and verify ",
    "the frozen base seed: ",
    pilot_checkpoint_path
  )
}

pilot_checkpoint <- readRDS(pilot_checkpoint_path)

if (
  is.null(pilot_checkpoint$config$base_seed) ||
  is.null(pilot_checkpoint$replication_stream)
) {
  stop(
    "The pilot checkpoint does not contain the frozen base seed ",
    "and replication stream."
  )
}

config <- list(
  formal_version = "2026-08-10-v2",
  base_seed = as.integer(
    pilot_checkpoint$config$base_seed
  ),
  pilot_n_replications = as.integer(
    pilot_checkpoint$config$n_replications
  ),
  n_replications = 50L,
  truth_structure = "coordinate",
  alpha = 1,
  neighbourhood = "rook",
  sigma = 0.5,
  n_side = 5L,
  n_grid = 50L,
  k_intercept = 10L,
  k_deviation = 15L,
  k_nodewise = 15L,
  k_pooled = 10L,
  maximum_failure_rate = 0.05,
  minimum_win_rate = 0.80,
  seed_allocation =
    "pilot positions 1:50; sanity positions 51:100"
)

if (config$pilot_n_replications != 50L) {
  stop(
    "The completed pilot does not contain the expected 50 ",
    "replications. Do not run the formal test."
  )
}

seed_ledger <- make_corrected_seed_ledger(
  base_seed = config$base_seed,
  n_pilot = config$pilot_n_replications,
  n_sanity = config$n_replications
)

pilot_checkpoint_paths <- file.path(
  pilot_checkpoint_directory,
  sprintf(
    "replication_%03d.rds",
    seq_len(config$pilot_n_replications)
  )
)

if (!all(file.exists(pilot_checkpoint_paths))) {
  stop(
    "All 50 pilot checkpoints are required to verify the corrected ",
    "seed allocation. Do not run the formal test."
  )
}

pilot_stream_matches <- vapply(
  seq_len(config$pilot_n_replications),
  function(replication_id) {
    checkpoint <- readRDS(
      pilot_checkpoint_paths[[replication_id]]
    )
    
    identical(
      checkpoint$config$base_seed,
      pilot_checkpoint$config$base_seed
    ) &&
      identical(
        checkpoint$replication_stream,
        seed_ledger$streams$pilot[[replication_id]]
      )
  },
  logical(1L)
)

if (!all(pilot_stream_matches)) {
  stop(
    "The reconstructed direct pilot streams do not reproduce pilot ",
    "checkpoint(s): ",
    paste(which(!pilot_stream_matches), collapse = ", "),
    ". Do not run the formal test."
  )
}

if (any(vapply(
  seed_ledger$streams$sanity,
  function(stream) {
    any(vapply(
      seed_ledger$streams$pilot,
      identical,
      logical(1L),
      y = stream
    ))
  },
  logical(1L)
))) {
  stop(
    "Pilot and formal-sanity streams overlap. Do not run the formal test."
  )
}

replication_seed_labels <-
  config$base_seed +
  config$pilot_n_replications +
  seq_len(config$n_replications) - 1L

message(
  "Verified all 50 pilot streams. Formal sanity uses the ",
  "non-overlapping continuation (stream positions 51-100)."
)

code_files <- c(
  "simulation/04_formal_known_answer.R",
  "simulation/R/make_graph.R",
  "simulation/R/rng.R",
  "simulation/R/generate_truth.R",
  "simulation/R/simulate_data.R",
  "simulation/R/model_diagnostics.R",
  "simulation/R/fit_baselines.R",
  "simulation/R/metrics.R",
  "simulation/R/run_replication.R",
  "simulation/R/run_core_cell.R",
  sort(list.files(
    "netfunsmooth/R",
    pattern = "[.]R$",
    full.names = TRUE
  ))
)

code_fingerprint <- make_code_fingerprint(code_files)

result_directory <-
  "simulation/results/core/formal_known_answer"

checkpoint_directory <- file.path(
  result_directory,
  "checkpoints"
)

dir.create(
  checkpoint_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# Run or resume the 50 attempted replications
# -------------------------------------------------------------------------

checkpoints <- vector(
  "list",
  config$n_replications
)

overall_start <- proc.time()[["elapsed"]]

for (replication_id in seq_len(config$n_replications)) {
  checkpoint_path <- file.path(
    checkpoint_directory,
    sprintf(
      "replication_%03d.rds",
      replication_id
    )
  )
  
  replication_stream <-
    seed_ledger$streams$sanity[[replication_id]]
  
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    
    validate_checkpoint(
      checkpoint = checkpoint,
      replication_id = replication_id,
      replication_stream = replication_stream,
      config = config,
      code_fingerprint = code_fingerprint
    )
    
    checkpoints[[replication_id]] <- checkpoint
    
    message(
      sprintf(
        "[Sanity %02d/%02d] validated existing checkpoint.",
        replication_id,
        config$n_replications
      )
    )
    
    next
  }
  
  message(
    sprintf(
      "[Sanity %02d/%02d] running.",
      replication_id,
      config$n_replications
    )
  )
  
  scenario_object <- make_sanity_scenario(
    replication_stream = replication_stream,
    config = config
  )
  
  replication_start <- proc.time()[["elapsed"]]
  
  cell_result <- run_core_cell(
    scenario = scenario_object$scenario,
    replication_id = replication_id,
    seed = replication_seed_labels[replication_id],
    k_intercept = config$k_intercept,
    k_deviation = config$k_deviation,
    k_nodewise = config$k_nodewise,
    k_pooled = config$k_pooled,
    keep_predictions = FALSE
  )
  
  replication_elapsed <-
    proc.time()[["elapsed"]] - replication_start
  
  checkpoint <- list(
    formal_version = config$formal_version,
    config = config,
    code_fingerprint = code_fingerprint,
    replication_id = as.integer(replication_id),
    replication_seed_label =
      replication_seed_labels[replication_id],
    replication_stream = replication_stream,
    component_streams =
      scenario_object$component_streams,
    result = cell_result,
    elapsed_seconds = replication_elapsed,
    completed_at = format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    )
  )
  
  save_rds_atomically(
    checkpoint,
    checkpoint_path
  )
  
  checkpoints[[replication_id]] <- checkpoint
}

total_elapsed_seconds <-
  proc.time()[["elapsed"]] - overall_start


# -------------------------------------------------------------------------
# Aggregate the prespecified paired comparison
# -------------------------------------------------------------------------

summary_table <- do.call(
  rbind,
  lapply(
    checkpoints,
    function(checkpoint) {
      checkpoint$result$summary
    }
  )
)

rownames(summary_table) <- NULL

expected_methods <- c(
  "raw",
  "pooled",
  "nodewise",
  "network"
)

method_counts <- table(summary_table$method)

if (
  nrow(summary_table) !=
  config$n_replications * length(expected_methods) ||
  !all(expected_methods %in% names(method_counts)) ||
  any(method_counts[expected_methods] != config$n_replications)
) {
  stop("The formal-test summary has unexpected method counts.")
}

if (anyDuplicated(
  summary_table[, c("replication", "method")]
)) {
  stop("Duplicate replication-method rows were found.")
}

nodewise_rows <- summary_table[
  summary_table$method == "nodewise",
  ,
  drop = FALSE
]

network_rows <- summary_table[
  summary_table$method == "network",
  ,
  drop = FALSE
]

nodewise_rows <- nodewise_rows[
  order(nodewise_rows$replication),
  ,
  drop = FALSE
]

network_rows <- network_rows[
  order(network_rows$replication),
  ,
  drop = FALSE
]

if (!identical(
  nodewise_rows$replication,
  network_rows$replication
)) {
  stop("Nodewise and network replication IDs are not aligned.")
}

nodewise_success <-
  nodewise_rows$success %in% TRUE &
  is.finite(nodewise_rows$node_averaged_ise)

network_success <-
  network_rows$success %in% TRUE &
  is.finite(network_rows$node_averaged_ise)

joint_success <-
  nodewise_success & network_success

paired_delta <-
  nodewise_rows$node_averaged_ise[joint_success] -
  network_rows$node_averaged_ise[joint_success]

n_joint <- length(paired_delta)

if (n_joint >= 2L) {
  mean_delta <- mean(paired_delta)
  paired_mcse <-
    stats::sd(paired_delta) / sqrt(n_joint)
  win_rate <- mean(paired_delta > 0)
} else {
  mean_delta <- NA_real_
  paired_mcse <- NA_real_
  win_rate <- NA_real_
}

nodewise_failure_rate <-
  1 - mean(nodewise_success)

network_failure_rate <-
  1 - mean(network_success)

criteria <- c(
  positive_mean_improvement =
    is.finite(mean_delta) && mean_delta > 0,
  improvement_exceeds_two_mcse =
    is.finite(mean_delta) &&
    is.finite(paired_mcse) &&
    mean_delta > 2 * paired_mcse,
  win_rate_at_least_80_percent =
    is.finite(win_rate) &&
    win_rate >= config$minimum_win_rate,
  nodewise_failure_rate_at_most_5_percent =
    nodewise_failure_rate <=
    config$maximum_failure_rate,
  network_failure_rate_at_most_5_percent =
    network_failure_rate <=
    config$maximum_failure_rate
)

formal_passed <- all(criteria)

paired_results <- data.frame(
  replication = nodewise_rows$replication,
  nodewise_success = nodewise_success,
  network_success = network_success,
  joint_success = joint_success,
  nodewise_aise = nodewise_rows$node_averaged_ise,
  network_aise = network_rows$node_averaged_ise,
  paired_delta = ifelse(
    joint_success,
    nodewise_rows$node_averaged_ise -
      network_rows$node_averaged_ise,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

formal_summary <- data.frame(
  attempted = config$n_replications,
  nodewise_successes = sum(nodewise_success),
  network_successes = sum(network_success),
  jointly_successful = n_joint,
  nodewise_failure_rate = nodewise_failure_rate,
  network_failure_rate = network_failure_rate,
  mean_paired_delta = mean_delta,
  paired_mcse = paired_mcse,
  two_paired_mcse = 2 * paired_mcse,
  win_rate = win_rate,
  formal_passed = formal_passed,
  row.names = NULL
)

formal_result <- list(
  formal_version = config$formal_version,
  config = config,
  seed_ledger = seed_ledger,
  code_fingerprint = code_fingerprint,
  criteria = criteria,
  formal_passed = formal_passed,
  formal_summary = formal_summary,
  paired_results = paired_results,
  method_summary = summary_table,
  checkpoint_paths = file.path(
    checkpoint_directory,
    sprintf(
      "replication_%03d.rds",
      seq_len(config$n_replications)
    )
  ),
  total_elapsed_seconds = total_elapsed_seconds,
  session_info = utils::sessionInfo(),
  completed_at = format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  )
)

result_path <- file.path(
  result_directory,
  "formal_known_answer_result.rds"
)

saveRDS(
  formal_result,
  file = result_path,
  version = 3
)


# -------------------------------------------------------------------------
# Report and enforce the gate
# -------------------------------------------------------------------------

cat("\nFormal known-answer summary:\n")
print(
  formal_summary,
  row.names = FALSE,
  digits = 6
)

cat("\nPrespecified criteria:\n")
print(criteria)

cat(
  "\nResult saved to:\n",
  result_path,
  "\n",
  sep = ""
)

if (!formal_passed) {
  stop(
    "The formal known-answer sanity check did not pass. ",
    "Pause the main simulation and inspect the implementation; ",
    "do not tune the DGP or pass criteria."
  )
}

message(
  "Formal known-answer sanity check passed."
)