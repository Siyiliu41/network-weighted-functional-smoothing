# Shared utilities for the fixed N = 100 scalability block.

normalised_file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot fingerprint missing file: ", path)
  }
  
  temporary_file <- tempfile(
    pattern = "normalised_",
    fileext = ".txt"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  
  normalised_text <- paste0(
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
    "\n"
  )
  
  writeBin(charToRaw(enc2utf8(normalised_text)), temporary_file)
  unname(tools::md5sum(temporary_file))
}


make_n100_code_fingerprint <- function(paths) {
  paths <- sort(unique(paths))
  
  hashes <- vapply(
    paths,
    normalised_file_md5,
    character(1L)
  )
  
  names(hashes) <- paths
  hashes
}


save_n100_rds_atomically <- function(object, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(path)) {
    stop("Refusing to overwrite existing file: ", path)
  }
  
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = directory,
    fileext = ".tmp"
  )
  
  on.exit(unlink(temporary_path), add = TRUE)
  
  saveRDS(object, temporary_path, version = 3)
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish file atomically: ", path)
  }
  
  invisible(path)
}


validate_n100_checkpoint <- function(
    checkpoint,
    replication_id,
    replication_stream,
    config,
    code_fingerprint) {
  
  required_entries <- c(
    "version",
    "config",
    "code_fingerprint",
    "replication_id",
    "replication_stream",
    "result"
  )
  
  if (!is.list(checkpoint) ||
      !all(required_entries %in% names(checkpoint))) {
    stop("N = 100 checkpoint has an incompatible structure.")
  }
  
  checks <- c(
    version = identical(checkpoint$version, config$version),
    config = identical(checkpoint$config, config),
    fingerprint = identical(checkpoint$code_fingerprint, code_fingerprint),
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
      "N = 100 checkpoint ",
      replication_id,
      " is incompatible: ",
      paste(names(checks)[!checks], collapse = ", "),
      "."
    )
  }
  
  result <- checkpoint$result
  
  if (!is.list(result) ||
      !is.data.frame(result$summary) ||
      nrow(result$summary) != 4L ||
      !is.numeric(result$cell_elapsed_seconds) ||
      length(result$cell_elapsed_seconds) != 2L) {
    stop("N = 100 checkpoint ", replication_id, " has an incomplete result.")
  }
  
  invisible(checkpoint)
}


make_n100_timing_table <- function(checkpoints) {
  expected_cells <- make_n100_design()$cell_id
  
  rows <- lapply(checkpoints, function(checkpoint) {
    result <- checkpoint$result
    summary <- result$summary
    
    lapply(expected_cells, function(cell_id) {
      cell_summary <- summary[
        summary$cell_id == cell_id,
        ,
        drop = FALSE
      ]
      
      network <- cell_summary[cell_summary$method == "network", , drop = FALSE]
      nodewise <- cell_summary[cell_summary$method == "nodewise", , drop = FALSE]
      
      if (nrow(network) != 1L || nrow(nodewise) != 1L) {
        stop("A timing replication is missing a primary-method record.")
      }
      
      data.frame(
        cell_id = cell_id,
        replication = result$replication,
        complete_elapsed_seconds =
          unname(result$cell_elapsed_seconds[[cell_id]]),
        network_success = network$success,
        nodewise_success = nodewise$success,
        stringsAsFactors = FALSE
      )
    })
  })
  
  timing_table <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(timing_table) <- NULL
  timing_table
}


choose_n100_replication_count <- function(
    timing_table,
    budget_hours = 8,
    n_timing_replications = 5L) {
  
  expected_cells <- make_n100_design()$cell_id
  
  if (!is.data.frame(timing_table) ||
      !all(c(
        "cell_id",
        "replication",
        "complete_elapsed_seconds"
      ) %in% names(timing_table))) {
    stop("The timing table is incomplete.")
  }
  
  if (!identical(as.integer(n_timing_replications), 5L) ||
      !identical(sort(unique(timing_table$cell_id)), sort(expected_cells)) ||
      anyDuplicated(timing_table[, c("cell_id", "replication")]) ||
      any(!is.finite(timing_table$complete_elapsed_seconds)) ||
      any(table(timing_table$cell_id) != n_timing_replications)) {
    stop(
      "The N = 100 timing pilot must contain exactly five finite timings ",
      "for each fixed cell."
    )
  }
  
  median_times <- stats::aggregate(
    complete_elapsed_seconds ~ cell_id,
    data = timing_table,
    FUN = stats::median
  )
  
  median_times <- median_times[
    match(expected_cells, median_times$cell_id),
    ,
    drop = FALSE
  ]
  
  candidates <- c(100L, 50L, 25L)
  projected_seconds <-
    candidates * sum(median_times$complete_elapsed_seconds)
  
  feasible <- candidates[
    projected_seconds <= budget_hours * 3600
  ]
  
  list(
    median_seconds_per_cell = median_times,
    projected_total_seconds = stats::setNames(
      projected_seconds,
      candidates
    ),
    attempted_replications_per_cell = if (length(feasible)) {
      max(feasible)
    } else {
      NA_integer_
    },
    computationally_feasible = length(feasible) > 0L
  )
}