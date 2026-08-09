# Reproducible random-number streams for the simulation study.
validate_rng_stream <- function(stream) {
  if (
    !is.integer(stream) ||
    length(stream) != 7L ||
    anyNA(stream)
  ) {
    stop(
      "`stream` must be a valid seven-integer ",
      "L'Ecuyer-CMRG RNG state."
    )
  }

  invisible(stream)
}

with_rng_stream <- function(stream, code) {
  validate_rng_stream(stream)

  old_kind <- RNGkind()
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit(
    {
      do.call(RNGkind, as.list(old_kind))

      if (had_seed) {
        assign(
          ".Random.seed",
          old_seed,
          envir = .GlobalEnv
        )
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        rm(
          ".Random.seed",
          envir = .GlobalEnv
        )
      }
    },
    add = TRUE
  )

  RNGkind("L'Ecuyer-CMRG")
  assign(
    ".Random.seed",
    stream,
    envir = .GlobalEnv
  )

  force(code)
}

make_initial_rng_stream <- function(base_seed) {
  if (
    !is.numeric(base_seed) ||
    length(base_seed) != 1L ||
    is.na(base_seed) ||
    !is.finite(base_seed)
  ) {
    stop("`base_seed` must be one finite number.")
  }

  old_kind <- RNGkind()
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit(
    {
      do.call(RNGkind, as.list(old_kind))

      if (had_seed) {
        assign(
          ".Random.seed",
          old_seed,
          envir = .GlobalEnv
        )
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        rm(
          ".Random.seed",
          envir = .GlobalEnv
        )
      }
    },
    add = TRUE
  )

  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(base_seed))

  .Random.seed
}

make_replication_streams <- function(
    base_seed,
    n_replications) {

  if (
    !is.numeric(n_replications) ||
    length(n_replications) != 1L ||
    is.na(n_replications) ||
    !is.finite(n_replications) ||
    n_replications < 1 ||
    n_replications %% 1 != 0
  ) {
    stop("`n_replications` must be a positive integer.")
  }

  n_replications <- as.integer(n_replications)

  streams <- vector(
    mode = "list",
    length = n_replications
  )

  streams[[1L]] <- make_initial_rng_stream(base_seed)

  if (n_replications > 1L) {
    for (i in 2:n_replications) {
      streams[[i]] <- parallel::nextRNGStream(
        streams[[i - 1L]]
      )
    }
  }

  names(streams) <- sprintf(
    "replication_%03d",
    seq_len(n_replications)
  )

  streams
}

make_component_substreams <- function(replication_stream) {
  validate_rng_stream(replication_stream)

  component_names <- c(
    "truth_draws",
    "error_coordinate_a05",
    "error_coordinate_a09",
    "error_cluster_a05",
    "error_cluster_a09"
  )

  streams <- vector(
    mode = "list",
    length = length(component_names)
  )

  streams[[1L]] <- replication_stream

  for (i in 2:length(streams)) {
    streams[[i]] <- parallel::nextRNGSubStream(
      streams[[i - 1L]]
    )
  }

  names(streams) <- component_names

  streams
}