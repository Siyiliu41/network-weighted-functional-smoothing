# Reproduce the Bavarian RKI COVID-19 application.
# See application/rki_covid/README.md for data provenance and outputs.

# Set this to the local clone of the package repository when `netfunsmooth`
# has not already been installed. Alternatively set the environment variable
# NETFUNSMOOTH_SOURCE_DIR before running the script.
netfunsmooth_source_dir <- Sys.getenv("NETFUNSMOOTH_SOURCE_DIR", unset = "")

analysis_start <- as.Date("2020-03-02")  # Monday
analysis_end <- as.Date("2023-03-05")    # Sunday
bkg_year <- "2020"
bkg_key_date <- "1231"

k_intercept <- 10L
k_deviation <- 15L
k_nodewise <- 15L

rki_url <- paste0(
  "https://media.githubusercontent.com/media/",
  "robert-koch-institut/",
  "SARS-CoV-2-Infektionen_in_Deutschland/main/",
  "Aktuell_Deutschland_SarsCov2_Infektionen.csv"
)

# ---- Paths and package checks -------------------------------------------

command_args <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_args, value = TRUE)

script_file <- if (length(file_argument) == 1L) {
  normalizePath(sub("^--file=", "", file_argument), winslash = "/")
} else {
  NA_character_
}

# When sourced interactively, the working directory must be the repository root.
repo_root <- if (!is.na(script_file)) {
  normalizePath(
    file.path(dirname(script_file), "..", ".."),
    winslash = "/"
  )
} else {
  normalizePath(getwd(), winslash = "/")
}

if (!file.exists(file.path(repo_root, "netfunsmooth", "DESCRIPTION"))) {
  stop(
    "`repo_root` must contain the `netfunsmooth` package directory: ",
    repo_root
  )
}

application_dir <- file.path(repo_root, "application", "rki_covid")
raw_dir <- file.path(application_dir, "data", "raw")
derived_dir <- file.path(application_dir, "data", "derived")
figure_dir <- file.path(application_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "data.table", "digest", "ffm", "ggplot2", "mgcv", "sf", "spdep", "tf"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (!nzchar(netfunsmooth_source_dir)) {
  netfunsmooth_source_dir <- file.path(repo_root, "netfunsmooth")
}

if (!dir.exists(netfunsmooth_source_dir)) {
  stop(
    "NETFUNSMOOTH_SOURCE_DIR does not exist: ",
    netfunsmooth_source_dir
  )
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Please install the `devtools` package first.")
}

devtools::load_all(netfunsmooth_source_dir, quiet = TRUE)

if (!requireNamespace("netfunsmooth", quietly = TRUE)) {
  stop("Could not load the local `netfunsmooth` package.")
}

library(data.table)

# ---- Helpers -------------------------------------------------------------

write_audit <- function(values, filename) {
  fwrite(
    data.table(metric = names(values), value = unname(values)),
    file.path(derived_dir, filename)
  )
}

fit_nodewise_smoothing <- function(observed, t_grid, k) {
  stopifnot(
    is.matrix(observed),
    is.numeric(observed),
    ncol(observed) == length(t_grid),
    all(is.finite(observed)),
    k >= 3L,
    k < length(t_grid)
  )
  
  fitted_values <- matrix(
    NA_real_, nrow = nrow(observed), ncol = ncol(observed),
    dimnames = dimnames(observed)
  )
  converged <- logical(nrow(observed))
  
  for (i in seq_len(nrow(observed))) {
    fit_i <- mgcv::gam(
      response ~ s(week, bs = "ps", k = k, m = c(2, 1)),
      data = data.frame(week = t_grid, response = observed[i, ]),
      method = "REML"
    )
    fitted_values[i, ] <- stats::predict(
      fit_i, newdata = data.frame(week = t_grid), type = "response"
    )
    converged[i] <- isTRUE(fit_i$converged)
  }
  
  list(fitted = fitted_values, converged = converged, k = k)
}

package_versions <- function(packages) {
  data.table(
    package = packages,
    version = vapply(packages, function(pkg) {
      as.character(utils::packageVersion(pkg))
    }, character(1))
  )
}

git_revision <- function(path) {
  if (!nzchar(path) || !dir.exists(path) ||
      !nzchar(Sys.which("git"))) {
    return(NA_character_)
  }
  result <- suppressWarnings(system2(
    "git", c("-C", path, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE
  ))
  if (length(result) == 1L && grepl("^[0-9a-f]{40}$", result)) result else NA_character_
}

# ---- RKI case data: frozen local snapshot and weekly aggregation --------

raw_file <- file.path(raw_dir, "rki_cases_download.csv")

if (!file.exists(raw_file)) {
  message("Downloading the current RKI CSV snapshot (approximately 400 MB).")
  temporary_file <- paste0(raw_file, ".part")
  on.exit(unlink(temporary_file), add = TRUE)
  utils::download.file(rki_url, destfile = temporary_file, mode = "wb", quiet = FALSE)
  if (file.info(temporary_file)$size < 100 * 1024^2) {
    stop("The RKI download is unexpectedly small and is not a valid CSV snapshot.")
  }
  if (!file.rename(temporary_file, raw_file)) {
    stop("Could not move the completed RKI download to: ", raw_file)
  }
}

raw_hash <- digest::digest(file = raw_file, algo = "sha256")
raw_info <- file.info(raw_file)
input_manifest <- data.table(
  item = "RKI case CSV",
  source_url = rki_url,
  local_file = normalizePath(raw_file, winslash = "/"),
  size_bytes = raw_info$size,
  modified_utc = format(raw_info$mtime, tz = "UTC", usetz = TRUE),
  sha256 = raw_hash,
  analysis_start = as.character(analysis_start),
  analysis_end = as.character(analysis_end),
  retrieved_or_checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
fwrite(input_manifest, file.path(derived_dir, "input_manifest.csv"))

rki_cases <- fread(
  raw_file,
  select = c("IdLandkreis", "Meldedatum", "AnzahlFall"),
  showProgress = TRUE
)
rki_cases[, Meldedatum := as.IDate(Meldedatum)]

case_audit <- c(
  n_rows = nrow(rki_cases),
  n_rki_units = uniqueN(rki_cases$IdLandkreis),
  missing_county_id = sum(is.na(rki_cases$IdLandkreis)),
  missing_report_date = sum(is.na(rki_cases$Meldedatum)),
  missing_case_count = sum(is.na(rki_cases$AnzahlFall)),
  first_reporting_date = as.character(min(rki_cases$Meldedatum)),
  last_reporting_date = as.character(max(rki_cases$Meldedatum))
)
write_audit(case_audit, "rki_case_audit.csv")

rki_window <- rki_cases[
  Meldedatum >= as.IDate(analysis_start) & Meldedatum <= as.IDate(analysis_end)
]
rki_window[, week_start := Meldedatum - ((as.integer(Meldedatum) + 3L) %% 7L)]

weekly_observed <- rki_window[
  , .(weekly_cases = sum(AnzahlFall)), by = .(IdLandkreis, week_start)
]

county_ids_411 <- sort(unique(rki_cases$IdLandkreis))
week_grid <- as.IDate(seq(analysis_start, analysis_end - 6L, by = "week"))
weekly_grid <- CJ(IdLandkreis = county_ids_411, week_start = week_grid, unique = TRUE)

weekly_cases <- merge(
  weekly_grid, weekly_observed,
  by = c("IdLandkreis", "week_start"), all.x = TRUE, sort = FALSE
)
weekly_cases[is.na(weekly_cases), weekly_cases := 0L]

# The RKI data report Berlin's twelve boroughs separately. Aggregate them to
# one Berlin unit (11000), matching the county-level BKG polygon layer.
weekly_cases[
  IdLandkreis >= 11001L & IdLandkreis <= 11012L,
  IdLandkreis := 11000L
]
weekly_cases <- weekly_cases[
  , .(weekly_cases = sum(weekly_cases)), by = .(IdLandkreis, week_start)
]
setorder(weekly_cases, IdLandkreis, week_start)
weekly_cases[, node_id := sprintf("%05d", IdLandkreis)]

weekly_audit <- c(
  n_nodes = uniqueN(weekly_cases$node_id),
  n_weeks = uniqueN(weekly_cases$week_start),
  expected_county_weeks = 400L * length(week_grid),
  observed_county_weeks = nrow(weekly_cases),
  duplicate_county_weeks = weekly_cases[, .N, by = .(node_id, week_start)][N > 1L, .N],
  missing_weekly_counts = sum(is.na(weekly_cases$weekly_cases)),
  negative_weekly_counts = sum(weekly_cases$weekly_cases < 0L),
  zero_county_weeks = sum(weekly_cases$weekly_cases == 0L)
)
write_audit(weekly_audit, "weekly_case_audit.csv")
stopifnot(
  weekly_audit["n_nodes"] == 400L,
  weekly_audit["expected_county_weeks"] == weekly_audit["observed_county_weeks"],
  weekly_audit["duplicate_county_weeks"] == 0L,
  weekly_audit["missing_weekly_counts"] == 0L,
  weekly_audit["negative_weekly_counts"] == 0L
)

# ---- County polygons, population denominators, and Queen graph ----------

kreise_sf <- ffm::bkg_admin_archive(
  level = "krs", scale = "250", key_date = bkg_key_date,
  year = bkg_year, timeout = 300
)
stopifnot(inherits(kreise_sf, "sf"))
names(kreise_sf) <- tolower(names(kreise_sf))

id_candidates <- c("ags", "ags_0", "ars", "ars_0")
id_field <- id_candidates[id_candidates %in% names(kreise_sf)][1L]
if (is.na(id_field)) {
  stop("No BKG administrative-code field found: ", paste(names(kreise_sf), collapse = ", "))
}

id_raw <- trimws(as.character(kreise_sf[[id_field]]))
if (anyNA(id_raw) || !all(grepl("^[0-9]+$", id_raw))) {
  stop("The BKG administrative-code field contains missing or non-numeric values.")
}
kreise_sf$node_id <- sprintf("%05d", as.integer(id_raw))

case_ids <- sort(unique(weekly_cases$node_id))

# The BKG layer contains 30 unpopulated areas sharing a county-code prefix
# with a regular Landkreis. They are not separate RKI observation units.
kreise_sf <- kreise_sf[
  kreise_sf$node_id %in% case_ids & !is.na(kreise_sf$ewz) & kreise_sf$ewz > 0,
]
kreise_sf <- kreise_sf[match(case_ids, kreise_sf$node_id), ]

polygon_audit <- c(
  n_matched_polygons = nrow(kreise_sf),
  n_unique_polygon_ids = uniqueN(kreise_sf$node_id),
  n_case_ids = length(case_ids),
  case_ids_without_polygons = length(setdiff(case_ids, kreise_sf$node_id)),
  duplicate_polygon_ids = sum(duplicated(kreise_sf$node_id)),
  missing_population = sum(is.na(kreise_sf$ewz)),
  nonpositive_population = sum(kreise_sf$ewz <= 0)
)
write_audit(polygon_audit, "polygon_audit.csv")
stopifnot(
  all(sf::st_is_valid(kreise_sf)),
  identical(kreise_sf$node_id, case_ids),
  polygon_audit["n_matched_polygons"] == 400L,
  polygon_audit["n_unique_polygon_ids"] == 400L,
  polygon_audit["case_ids_without_polygons"] == 0L,
  polygon_audit["duplicate_polygon_ids"] == 0L,
  polygon_audit["missing_population"] == 0L,
  polygon_audit["nonpositive_population"] == 0L
)

nb_queen <- spdep::poly2nb(kreise_sf, queen = TRUE, row.names = case_ids)
A_queen <- spdep::nb2mat(nb_queen, style = "B", zero.policy = TRUE)
rownames(A_queen) <- case_ids
colnames(A_queen) <- case_ids
degree <- spdep::card(nb_queen)
components <- spdep::n.comp.nb(nb_queen)

graph_audit <- c(
  n_nodes = length(case_ids),
  n_edges = sum(A_queen) / 2,
  min_degree = min(degree),
  median_degree = median(degree),
  mean_degree = mean(degree),
  max_degree = max(degree),
  n_isolated_nodes = sum(degree == 0L),
  n_connected_components = components$nc,
  adjacency_symmetric = isTRUE(all.equal(A_queen, t(A_queen)))
)
write_audit(graph_audit, "germany_queen_graph_audit.csv")
stopifnot(
  all(diag(A_queen) == 0),
  graph_audit["n_nodes"] == 400L,
  graph_audit["n_isolated_nodes"] == 0L,
  graph_audit["n_connected_components"] == 1L,
  graph_audit["adjacency_symmetric"] == 1
)

# ---- Incidence transformation and Bavaria subset ------------------------

population_lookup <- data.table(
  node_id = kreise_sf$node_id,
  population_2020 = as.numeric(kreise_sf$ewz)
)
weekly_analysis <- copy(weekly_cases)
weekly_analysis[
  population_lookup, on = .(node_id), population_2020 := i.population_2020
]
weekly_analysis[, weekly_incidence_100k := 100000 * weekly_cases / population_2020]
weekly_analysis[, log_weekly_incidence_100k := log1p(weekly_incidence_100k)]
setorder(weekly_analysis, node_id, week_start)

incidence_audit <- c(
  n_nodes = uniqueN(weekly_analysis$node_id),
  n_weeks = uniqueN(weekly_analysis$week_start),
  n_county_weeks = nrow(weekly_analysis),
  missing_incidence = sum(is.na(weekly_analysis$weekly_incidence_100k)),
  negative_incidence = sum(weekly_analysis$weekly_incidence_100k < 0),
  zero_incidence_weeks = sum(weekly_analysis$weekly_incidence_100k == 0),
  maximum_weekly_incidence_100k = max(weekly_analysis$weekly_incidence_100k),
  maximum_log_incidence = max(weekly_analysis$log_weekly_incidence_100k)
)
write_audit(incidence_audit, "incidence_audit.csv")
stopifnot(
  incidence_audit["n_nodes"] == 400L,
  incidence_audit["n_weeks"] == 157L,
  incidence_audit["missing_incidence"] == 0L,
  incidence_audit["negative_incidence"] == 0L
)

weekly_analysis[, week_id := as.character(week_start)]
incidence_wide <- dcast(
  weekly_analysis, node_id ~ week_id, value.var = "log_weekly_incidence_100k"
)
week_ids <- names(incidence_wide)[-1L]
log_incidence_matrix <- as.matrix(incidence_wide[, -1L])
storage.mode(log_incidence_matrix) <- "double"
rownames(log_incidence_matrix) <- incidence_wide$node_id
colnames(log_incidence_matrix) <- week_ids
stopifnot(identical(rownames(log_incidence_matrix), case_ids), all(is.finite(log_incidence_matrix)))

# Full Germany requires a multi-gigabyte allocation in the present pffr/mgcv
# implementation. The fitted application is therefore the induced, connected
# Bavaria subgraph (96 counties), not a silently truncated Germany-wide model.
bavaria_ids <- case_ids[substr(case_ids, 1L, 2L) == "09"]
bavaria_index <- match(bavaria_ids, case_ids)
rki_bavaria_matrix <- log_incidence_matrix[bavaria_index, , drop = FALSE]
A_bavaria <- A_queen[bavaria_index, bavaria_index, drop = FALSE]
kreise_bavaria_sf <- kreise_sf[bavaria_index, ]

nb_bavaria <- spdep::mat2listw(A_bavaria, style = "B", zero.policy = TRUE)$neighbours
bavaria_degree <- spdep::card(nb_bavaria)
bavaria_components <- spdep::n.comp.nb(nb_bavaria)
bavaria_graph_audit <- c(
  n_nodes = length(bavaria_ids),
  n_edges = sum(A_bavaria) / 2,
  min_degree = min(bavaria_degree),
  median_degree = median(bavaria_degree),
  mean_degree = mean(bavaria_degree),
  max_degree = max(bavaria_degree),
  n_isolated_nodes = sum(bavaria_degree == 0L),
  n_connected_components = bavaria_components$nc
)
write_audit(bavaria_graph_audit, "bavaria_queen_graph_audit.csv")
stopifnot(
  length(bavaria_ids) == 96L,
  identical(rownames(rki_bavaria_matrix), bavaria_ids),
  identical(rownames(A_bavaria), bavaria_ids),
  identical(colnames(A_bavaria), bavaria_ids),
  identical(kreise_bavaria_sf$node_id, bavaria_ids),
  bavaria_graph_audit["n_isolated_nodes"] == 0L,
  bavaria_graph_audit["n_connected_components"] == 1L
)

saveRDS(
  list(
    node_ids = bavaria_ids,
    week_ids = week_ids,
    observed_log_incidence = rki_bavaria_matrix,
    graph = A_bavaria,
    polygons = kreise_bavaria_sf,
    population_2020 = population_lookup[match(bavaria_ids, node_id)],
    configuration = list(
      analysis_start = analysis_start, analysis_end = analysis_end,
      bkg_year = bkg_year, bkg_key_date = bkg_key_date,
      outcome = "log1p(weekly cases per 100,000 residents)"
    )
  ),
  file.path(derived_dir, "bavaria_application_inputs.rds")
)

# ---- Nodewise and network fits ------------------------------------------

t_index <- seq_len(ncol(rki_bavaria_matrix))
rki_bavaria_curves <- tf::tfd(rki_bavaria_matrix, arg = t_index)
names(rki_bavaria_curves) <- bavaria_ids

nodewise_fit <- fit_nodewise_smoothing(
  observed = rki_bavaria_matrix, t_grid = t_index, k = k_nodewise
)
stopifnot(all(nodewise_fit$converged), all(is.finite(nodewise_fit$fitted)))

network_fit <- netfunsmooth::netf_smooth(
  curves = rki_bavaria_curves,
  graph = A_bavaria,
  sandwich = "none",
  bs.int = list(bs = "ps", k = k_intercept, m = c(2, 1)),
  bs.yindex = list(bs = "ps", k = k_deviation, m = c(2, 1))
)
stopifnot(isTRUE(network_fit$model$converged))

network_curves <- stats::predict(network_fit)
network_fitted <- as.matrix(network_curves, arg = t_index)
rownames(network_fitted) <- names(network_curves)
colnames(network_fitted) <- week_ids
stopifnot(
  identical(rownames(network_fitted), bavaria_ids),
  all(is.finite(network_fitted))
)

fit_audit <- c(
  n_nodes = nrow(rki_bavaria_matrix),
  n_weeks = ncol(rki_bavaria_matrix),
  nodewise_converged = sum(nodewise_fit$converged),
  network_converged = isTRUE(network_fit$model$converged),
  network_edf = sum(network_fit$model$edf),
  mean_abs_network_minus_nodewise = mean(abs(network_fitted - nodewise_fit$fitted)),
  max_abs_network_minus_nodewise = max(abs(network_fitted - nodewise_fit$fitted))
)
write_audit(fit_audit, "fit_audit_bavaria.csv")

# Keep the same compact result object used in the interactive analysis. This
# also means the plotting code below can be rerun without refitting, provided
# `rki_bavaria_results` has been retained in the R session.
rki_bavaria_results <- list(
  node_ids = bavaria_ids,
  week_ids = week_ids,
  observed_log_incidence = rki_bavaria_matrix,
  nodewise_fitted = nodewise_fit$fitted,
  network_fitted = network_fitted,
  graph = A_bavaria,
  polygons = kreise_bavaria_sf
)

# ---- Figures -------------------------------------------------------------

plot_family <- if (identical(.Platform$OS.type, "windows")) "Arial" else "sans"
publication_theme <- ggplot2::theme_minimal(base_size = 10, base_family = plot_family) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

map_difference <- data.frame(
  node_id = bavaria_ids,
  mean_abs_difference = rowMeans(abs(network_fitted - nodewise_fit$fitted))
)
kreise_bavaria_plot <- merge(
  kreise_bavaria_sf, map_difference, by = "node_id", all.x = TRUE, sort = FALSE
)

fig_map <- ggplot2::ggplot(kreise_bavaria_plot) +
  ggplot2::geom_sf(
    ggplot2::aes(fill = mean_abs_difference),
    colour = "white",
    linewidth = 0.08
  ) +
  ggplot2::scale_fill_viridis_c(
    option = "C",
    name = "Mean absolute difference\non log incidence scale"
  ) +
  ggplot2::labs(
    title = "Difference between network and nodewise fitted curves",
    subtitle = "Mean absolute difference across the 157 weekly observations"
  ) +
  ggplot2::theme_void(base_family = plot_family) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

selected_counties <- c("09162", "09564", "09172", "09663")
selected_names <- kreise_bavaria_sf$gen[match(selected_counties, bavaria_ids)]
if (anyNA(selected_names)) stop("A prespecified county is missing from the Bavaria subset.")
county_rows <- match(selected_counties, rki_bavaria_results$node_ids)

plot_grid <- expand.grid(
  node_id = selected_counties, week_index = seq_along(week_ids),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
plot_grid$county <- selected_names[match(plot_grid$node_id, selected_counties)]
plot_grid$week_start <- as.Date(week_ids[plot_grid$week_index])
plot_grid$observed <- as.vector(t(
  rki_bavaria_results$observed_log_incidence[county_rows, ]
))
plot_grid$nodewise <- as.vector(t(
  rki_bavaria_results$nodewise_fitted[county_rows, ]
))
plot_grid$network <- as.vector(t(
  rki_bavaria_results$network_fitted[county_rows, ]
))

plot_long <- rbind(
  data.frame(
    plot_grid[c("node_id", "county", "week_start")],
    series = "Observed weekly incidence",
    value = plot_grid$observed
  ),
  data.frame(
    plot_grid[c("node_id", "county", "week_start")],
    series = "Nodewise smoothing",
    value = plot_grid$nodewise
  ),
  data.frame(
    plot_grid[c("node_id", "county", "week_start")],
    series = "Network-weighted smoothing",
    value = plot_grid$network
  )
)

plot_long$series <- factor(
  plot_long$series,
  levels = c(
    "Observed weekly incidence",
    "Nodewise smoothing",
    "Network-weighted smoothing"
  )
)

fig_curves <- ggplot2::ggplot(
  plot_long,
  ggplot2::aes(
    x = week_start,
    y = value,
    colour = series,
    linetype = series,
    group = series
  )
) +
  ggplot2::geom_line(
    data = subset(plot_long, series == "Observed weekly incidence"),
    linewidth = 0.30,
    alpha = 0.60
  ) +
  ggplot2::geom_line(
    data = subset(plot_long, series == "Nodewise smoothing"),
    linewidth = 0.85
  ) +
  ggplot2::geom_line(
    data = subset(plot_long, series == "Network-weighted smoothing"),
    linewidth = 0.85
  ) +
  ggplot2::facet_wrap(~ county, ncol = 2) +
  ggplot2::scale_colour_manual(
    values = c(
      "Observed weekly incidence" = "grey45",
      "Nodewise smoothing" = "#D95F02",
      "Network-weighted smoothing" = "#1B9E77"
    ),
    breaks = c(
      "Observed weekly incidence",
      "Nodewise smoothing",
      "Network-weighted smoothing"
    ),
    name = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Observed weekly incidence" = "solid",
      "Nodewise smoothing" = "dashed",
      "Network-weighted smoothing" = "solid"
    ),
    breaks = c(
      "Observed weekly incidence",
      "Nodewise smoothing",
      "Network-weighted smoothing"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Weekly COVID-19 incidence curves in selected Bavarian counties",
    x = NULL,
    y = "log(1 + weekly reported cases per 100,000 residents)"
  ) +
  publication_theme +
  ggplot2::theme(
    legend.position = "bottom"
  )

if (!capabilities("cairo")) {
  warning("Cairo is unavailable; PDFs use the default graphics device and may have poorer font embedding.")
  pdf_device <- "pdf"
} else {
  pdf_device <- grDevices::cairo_pdf
}

ggplot2::ggsave(
  file.path(figure_dir, "bavaria_smoothing_difference_map.pdf"), fig_map,
  width = 7.2, height = 6.0, device = pdf_device
)
ggplot2::ggsave(
  file.path(figure_dir, "bavaria_curve_comparison.pdf"), fig_curves,
  width = 7.2, height = 5.6, device = pdf_device
)

# ---- Persist results and computational provenance ----------------------

results <- list(
  node_ids = bavaria_ids,
  week_ids = week_ids,
  observed_log_incidence = rki_bavaria_matrix,
  nodewise_fitted = nodewise_fit$fitted,
  network_fitted = network_fitted,
  nodewise_fit = nodewise_fit,
  network_fit = network_fit,
  graph = A_bavaria,
  polygons = kreise_bavaria_sf,
  audits = list(
    cases = case_audit, weekly = weekly_audit, polygons = polygon_audit,
    germany_graph = graph_audit, incidence = incidence_audit,
    bavaria_graph = bavaria_graph_audit, fit = fit_audit
  ),
  configuration = list(
    analysis_start = analysis_start, analysis_end = analysis_end,
    bkg_year = bkg_year, bkg_key_date = bkg_key_date,
    k_intercept = k_intercept, k_deviation = k_deviation,
    k_nodewise = k_nodewise, outcome = "log1p(weekly cases per 100,000 residents)"
  )
)
saveRDS(results, file.path(derived_dir, "bavaria_application_results.rds"))

all_packages <- c(required_packages, "netfunsmooth")
provenance <- list(
  run_finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  project_dir = repo_root,
  input_manifest = input_manifest,
  packages = package_versions(all_packages),
  r_version = R.version.string,
  platform = R.version$platform,
  netfunsmooth_source_dir = if (nzchar(netfunsmooth_source_dir)) normalizePath(netfunsmooth_source_dir, winslash = "/") else NA_character_,
  netfunsmooth_git_revision = git_revision(netfunsmooth_source_dir),
  session_info = utils::sessionInfo()
)
saveRDS(provenance, file.path(derived_dir, "computational_provenance.rds"))
fwrite(provenance$packages, file.path(derived_dir, "package_versions.csv"))

message("Reproduction complete. Key fit audit:")
print(fit_audit)
message("Outputs written below: ", application_dir)