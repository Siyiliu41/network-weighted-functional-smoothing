#!/usr/bin/env Rscript

# Reproducible DWD application for the Bachelor thesis
#
# Purpose
# -------
# Recreate the complete German DWD temperature application from the frozen
# station registry through download, preprocessing, graph construction,
# nodewise and network-weighted fitting, and all thesis outputs.
#
# Run from the repository root:
#   Rscript application/dwd_temperature/reproduce_dwd_temperature_application.R
#
# Required local project layout
# -----------------------------
#   netfunsmooth/                 package source directory
#   simulation/R/fit_baselines.R  nodewise comparator used in the simulation
#
# The script intentionally fixes the 25 stations selected for the thesis rather
# than repeating the exploratory station-selection step.  This makes the final
# application exactly reproducible even when the DWD station index changes.

options(stringsAsFactors = FALSE)

required_packages <- c(
  "devtools", "rdwd", "sf", "tf", "igraph", "ggplot2"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "), call. = FALSE
  )
}

find_project_root <- function(start = getwd()) {
  here <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "netfunsmooth")) &&
        file.exists(file.path(here, "simulation", "R", "fit_baselines.R"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop(
        "Could not find the project root. Run this script from the repository " ,
        "or keep both netfunsmooth/ and simulation/R/fit_baselines.R in it.",
        call. = FALSE
      )
    }
    here <- parent
  }
}

project_root <- find_project_root()
application_dir <- file.path(project_root, "application", "dwd_temperature")
data_dir <- file.path(application_dir, "data")
raw_dir <- file.path(data_dir, "raw")
derived_dir <- file.path(data_dir, "derived")
figure_dir <- file.path(application_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Frozen thesis registry. Station IDs are character strings to preserve leading
# zeroes. target_city controls the labels and the six curve panels.
station_registry <- data.frame(
  station_id = c(
    "04271", "04466", "01975", "00691", "01550", "00433", "02014",
    "03126", "02928", "01443", "06217", "03379", "01766", "15000",
    "01303", "15444", "04928", "05906", "05856", "01050", "02667",
    "01424", "04104", "03668", "15207"
  ),
  target_city = c(
    "Rostock", "Schleswig", "Hamburg", "Bremen", "Garmisch", "Berlin",
    "Hannover", "Magdeburg", "Leipzig", "Freiburg", "Saarbruecken",
    "Munich", "Muenster", "Aachen", "Essen", "Ulm", "Stuttgart",
    "Mannheim", "Passau", "Dresden", "Cologne", "Frankfurt", "Regensburg",
    "Nuremberg", "Kassel"
  )
)

analysis_start <- as.Date("2024-01-01")
analysis_end <- as.Date("2024-12-31")
knn_k <- 4L
selected_cities <- c("Garmisch", "Freiburg", "Frankfurt", "Berlin", "Hamburg", "Rostock")

write.csv(
  station_registry,
  file.path(derived_dir, "dwd_station_registry.csv"),
  row.names = FALSE
)

devtools::load_all(file.path(project_root, "netfunsmooth"), quiet = TRUE)
source(file.path(project_root, "simulation", "R", "fit_baselines.R"))

first_data_frame <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    idx <- which(vapply(x, is.data.frame, logical(1)))
    if (length(idx) > 0L) return(x[[idx[1L]]])
  }
  stop("rdwd::readDWD() did not return a data frame.", call. = FALSE)
}

download_and_read_station <- function(station_id) {
  # Historical files are required because the analysis year is 2024.
  link <- rdwd::selectDWD(
    id = as.integer(station_id), res = "daily", var = "kl", per = "historical",
    quiet = TRUE
  )
  if (length(link) != 1L) {
    stop("Expected one historical daily climate file for station ", station_id, ".", call. = FALSE)
  }
  first_data_frame(rdwd::dataDWD(
    link,
    dir = raw_dir,
    varnames = TRUE,
    quiet = TRUE
  ))
}

extract_temperature_series <- function(x, station_id) {
  date_column <- grep("^MESS_DATUM$", names(x), value = TRUE)
  temperature_column <- grep("^TMK$|TMK.Lufttemperatur", names(x), value = TRUE)

  if (length(date_column) != 1L || length(temperature_column) != 1L) {
    stop("Expected MESS_DATUM and TMK in station ", station_id, ".", call. = FALSE)
  }

  raw_date <- x[[date_column]]
  if (inherits(raw_date, "Date")) {
    date <- as.Date(raw_date)
  } else {
    date_text <- trimws(as.character(raw_date))
    date <- as.Date(date_text, format = "%Y%m%d")
    unresolved <- is.na(date)
    date[unresolved] <- as.Date(date_text[unresolved])
  }

  temperature <- suppressWarnings(
    as.numeric(as.character(x[[temperature_column]]))
  )
  temperature[temperature == -999] <- NA_real_

  out <- data.frame(date = date, temperature = temperature)
  out <- out[
    out$date >= analysis_start & out$date <= analysis_end,
    ,
    drop = FALSE
  ]

  expected_dates <- seq(analysis_start, analysis_end, by = "day")
  missing_days <- sum(!expected_dates %in% out$date)
  duplicate_days <- sum(duplicated(out$date))
  missing_temperatures <- sum(is.na(out$temperature))

  if (nrow(out) != length(expected_dates) ||
      duplicate_days > 0L ||
      missing_temperatures > 0L ||
      missing_days > 0L) {
    stop(
      "Station ", station_id, " failed the 2024 completeness check ",
      "(rows = ", nrow(out),
      ", missing days = ", missing_days,
      ", duplicate days = ", duplicate_days,
      ", missing TMK values = ", missing_temperatures, ").",
      call. = FALSE
    )
  }

  out[order(out$date), , drop = FALSE]
}

message("Downloading and reading the frozen DWD station registry ...")
station_series <- lapply(station_registry$station_id, function(id) {
  extract_temperature_series(download_and_read_station(id), id)
})
names(station_series) <- station_registry$station_id

dates <- station_series[[1L]]$date
stopifnot(all(vapply(station_series, function(x) identical(x$date, dates), logical(1))))
temperature_matrix <- do.call(rbind, lapply(station_series, `[[`, "temperature"))
rownames(temperature_matrix) <- station_registry$station_id
colnames(temperature_matrix) <- as.character(dates)

# Retrieve the station coordinates from the same DWD metadata service and match
# them by fixed station ID. This avoids hand-copied coordinates in the analysis.
data("metaIndex", package = "rdwd", envir = environment())
metadata_index <- get("metaIndex", envir = environment())
id_column <- grep("Stations_id|station_id|^id$", names(metadata_index), ignore.case = TRUE, value = TRUE)[1L]
lon_column <- grep("geoLaenge|longitude", names(metadata_index), ignore.case = TRUE, value = TRUE)[1L]
lat_column <- grep("geoBreite|latitude", names(metadata_index), ignore.case = TRUE, value = TRUE)[1L]
name_column <- grep("Stationsname|station_name", names(metadata_index), ignore.case = TRUE, value = TRUE)[1L]
if (anyNA(c(id_column, lon_column, lat_column, name_column))) {
  stop("The DWD metadata index has unexpected column names.", call. = FALSE)
}
metadata_index[[id_column]] <- sprintf("%05d", as.integer(metadata_index[[id_column]]))
station_metadata <- merge(
  station_registry, metadata_index,
  by.x = "station_id", by.y = id_column, all.x = TRUE, sort = FALSE
)
station_metadata <- station_metadata[match(station_registry$station_id, station_metadata$station_id), , drop = FALSE]
if (anyNA(station_metadata[[lon_column]]) || anyNA(station_metadata[[lat_column]])) {
  stop("At least one frozen station could not be matched in the DWD metadata index.", call. = FALSE)
}
station_metadata$Stationsname <- station_metadata[[name_column]]
station_metadata$geoLaenge <- as.numeric(station_metadata[[lon_column]])
station_metadata$geoBreite <- as.numeric(station_metadata[[lat_column]])
station_metadata$node_id <- station_metadata$station_id

write.csv(
  station_metadata,
  file.path(derived_dir, "dwd_station_metadata.csv"),
  row.names = FALSE
)
write.csv(
  temperature_matrix,
  file.path(derived_dir, "dwd_temperature_2024.csv")
)

# Build the symmetric k-nearest-neighbour graph in projected ETRS89 / LAEA Europe
# coordinates (EPSG:3035), not in longitude/latitude.
station_sf <- sf::st_as_sf(
  station_metadata, coords = c("geoLaenge", "geoBreite"), crs = 4326, remove = FALSE
)
station_xy_m <- sf::st_coordinates(sf::st_transform(station_sf, 3035))
distance_m <- as.matrix(stats::dist(station_xy_m))
diag(distance_m) <- Inf
directed_knn <- matrix(0L, nrow(distance_m), ncol(distance_m),
                       dimnames = list(station_registry$station_id, station_registry$station_id))
for (i in seq_len(nrow(directed_knn))) {
  directed_knn[i, order(distance_m[i, ])[seq_len(knn_k)]] <- 1L
}
adjacency_knn4 <- 1L * ((directed_knn + t(directed_knn)) > 0L)
diag(adjacency_knn4) <- 0L
stopifnot(identical(adjacency_knn4, t(adjacency_knn4)))

edge_pairs <- which(upper.tri(adjacency_knn4) & adjacency_knn4 == 1L, arr.ind = TRUE)
edge_lengths_km <- sqrt(rowSums((station_xy_m[edge_pairs[, 1L], , drop = FALSE] -
                                   station_xy_m[edge_pairs[, 2L], , drop = FALSE])^2)) / 1000
edge_summary <- data.frame(
  n_nodes = nrow(adjacency_knn4), n_edges = nrow(edge_pairs),
  n_components = length(igraph::components(igraph::graph_from_adjacency_matrix(adjacency_knn4))$csize),
  min_edge_km = min(edge_lengths_km), median_edge_km = stats::median(edge_lengths_km),
  max_edge_km = max(edge_lengths_km), min_degree = min(rowSums(adjacency_knn4)),
  max_degree = max(rowSums(adjacency_knn4))
)
print(round(edge_summary, 1))
write.csv(
  adjacency_knn4,
  file.path(derived_dir, "dwd_adjacency_knn4.csv")
)
write.csv(
  edge_summary,
  file.path(derived_dir, "dwd_knn4_graph_summary.csv"),
  row.names = FALSE
)

# Fit the prespecified primary analysis. These basis dimensions are identical to
# the final simulation settings; smoothing parameters are selected by REML.
t_grid <- seq(0, 1, length.out = ncol(temperature_matrix))
temperature_curves <- tf::tfd(temperature_matrix, arg = t_grid)
names(temperature_curves) <- station_registry$station_id

nodewise_fit <- fit_nodewise_smoothing(
  observed = temperature_matrix, t_grid = t_grid, k = 15L, keep_models = FALSE
)
network_fit <- netf_smooth(
  curves = temperature_curves, graph = adjacency_knn4, sandwich = "none",
  bs.int = list(bs = "ps", k = 10L, m = c(2, 1)),
  bs.yindex = list(bs = "ps", k = 15L, m = c(2, 1))
)
if (!isTRUE(network_fit$model$converged)) stop("The network model did not converge.", call. = FALSE)
network_predictions <- as.matrix(predict(network_fit))
rownames(network_predictions) <- station_registry$station_id
colnames(network_predictions) <- colnames(temperature_matrix)
stopifnot(identical(dim(network_predictions), dim(temperature_matrix)), all(is.finite(network_predictions)))

application_fits <- list(
  raw = temperature_matrix, nodewise = nodewise_fit$predictions, network = network_predictions,
  t_grid = t_grid, dates = dates, station_metadata = station_metadata,
  graph = list(adjacency = adjacency_knn4, k = knn_k, edge_summary = edge_summary),
  basis = list(intercept_k = 10L, deviation_k = 15L, nodewise_k = 15L),
  network_model = network_fit$model
)
saveRDS(
  application_fits,
  file.path(derived_dir, "dwd_temperature_2024_primary_fits.rds")
)
write.csv(
  as.data.frame(network_fit$model$sp),
  file.path(derived_dir, "dwd_network_smoothing_parameters.csv")
)

# Thesis Figure: graph. Coordinates are shown in longitude/latitude only for display.
selected_indices <- match(selected_cities, station_metadata$target_city)
png(file.path(figure_dir, "dwd_knn4_station_graph.png"), width = 2000, height = 1800, res = 220)
par(mar = c(4.5, 4.8, 2.8, 1.0))
plot(station_metadata$geoLaenge, station_metadata$geoBreite, type = "n", asp = 1.45,
     xlab = "Longitude (°E)", ylab = "Latitude (°N)",
     main = "DWD stations and symmetrised 4-nearest-neighbour graph")
for (e in seq_len(nrow(edge_pairs))) {
  i <- edge_pairs[e, 1L]; j <- edge_pairs[e, 2L]
  segments(station_metadata$geoLaenge[i], station_metadata$geoBreite[i],
           station_metadata$geoLaenge[j], station_metadata$geoBreite[j], col = "grey75", lwd = 1.2)
}
points(station_metadata$geoLaenge, station_metadata$geoBreite, pch = 21, bg = "#1F78B4", col = "white", cex = 1.45)
points(station_metadata$geoLaenge[selected_indices], station_metadata$geoBreite[selected_indices],
       pch = 21, bg = "#E7298A", col = "white", cex = 1.65)
text(station_metadata$geoLaenge[selected_indices], station_metadata$geoBreite[selected_indices],
     labels = station_metadata$Stationsname[selected_indices], pos = 3, cex = 0.78)
legend("bottomleft", legend = c("Station", "Station shown in Figure 1", "k = 4 graph edge"),
       pch = c(21, 21, NA), pt.bg = c("#1F78B4", "#E7298A", NA), pt.cex = c(1.3, 1.5, NA),
       lty = c(NA, NA, 1), lwd = c(NA, NA, 1.2), col = c("black", "black", "grey60"), bty = "n", cex = 0.85)
dev.off()
# Thesis Figure 3: observed and fitted temperature curves.

global_ylim <- range(
  temperature_matrix,
  nodewise_fit$predictions,
  network_predictions,
  finite = TRUE
)

english_ticks <- as.Date(c(
  "2024-01-01",
  "2024-04-01",
  "2024-07-01",
  "2024-10-01",
  "2025-01-01"
))

selected_station_names <- station_metadata$Stationsname[selected_indices]

panel_labels <- paste0(
  "(", letters[seq_along(selected_station_names)], ") ",
  selected_station_names
)

curve_plot_data <- do.call(
  rbind,
  lapply(seq_along(selected_indices), function(panel_number) {
    i <- selected_indices[panel_number]

    data.frame(
      date = rep(dates, times = 3L),
      temperature = c(
        temperature_matrix[i, ],
        nodewise_fit$predictions[i, ],
        network_predictions[i, ]
      ),
      method = rep(
        c("Observed", "Nodewise", "Network"),
        each = length(dates)
      ),
      station = panel_labels[panel_number]
    )
  })
)

curve_plot_data$method <- factor(
  curve_plot_data$method,
  levels = c("Observed", "Nodewise", "Network")
)

curve_plot_data$station <- factor(
  curve_plot_data$station,
  levels = panel_labels
)

curve_plot <- ggplot2::ggplot(
  curve_plot_data,
  ggplot2::aes(
    x = date,
    y = temperature,
    colour = method,
    linetype = method,
    group = method
  )
) +
  ggplot2::geom_line(
    data = subset(curve_plot_data, method == "Observed"),
    linewidth = 0.30,
    alpha = 0.60
  ) +
  ggplot2::geom_line(
    data = subset(curve_plot_data, method == "Nodewise"),
    linewidth = 0.85
  ) +
  ggplot2::geom_line(
    data = subset(curve_plot_data, method == "Network"),
    linewidth = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ station,
    ncol = 2
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Observed" = "grey45",
      "Nodewise" = "#D95F02",
      "Network" = "#1B9E77"
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Observed" = "solid",
      "Nodewise" = "dashed",
      "Network" = "solid"
    )
  ) +
  ggplot2::scale_x_date(
    breaks = english_ticks,
    labels = c(
      "Jan\n2024",
      "Apr\n2024",
      "Jul\n2024",
      "Oct\n2024",
      "Jan\n2025"
    ),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = global_ylim,
    expand = ggplot2::expansion(mult = c(0.03, 0.05))
  ) +
  ggplot2::labs(
    title = "DWD daily mean temperatures in 2024:\nselected station fits",
    x = NULL,
    y = "Daily mean temperature (°C)",
    colour = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 15,
      hjust = 0.5,
      margin = ggplot2::margin(b = 12)
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11,
      margin = ggplot2::margin(b = 6)
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 8)
    ),
    axis.text = ggplot2::element_text(size = 9),
    axis.ticks.length = grid::unit(2, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = ggplot2::element_text(size = 10),
    legend.key.width = grid::unit(1.6, "cm"),
    panel.spacing = grid::unit(1.0, "lines"),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 8,
      l = 10
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "dwd_curve_comparison_selected_stations.pdf"
  ),
  plot = curve_plot,
  width = 7.5,
  height = 10,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE
)

# Descriptive diagnostic only: no true temperature functions are observed.
pairwise_curve_rms <- function(curves, pairs) sqrt(rowMeans((curves[pairs[, 1L], ] - curves[pairs[, 2L], ])^2))
all_pairs <- which(upper.tri(adjacency_knn4), arr.ind = TRUE)
is_edge <- adjacency_knn4[all_pairs] == 1L
summarise_pair_similarity <- function(curves, method) {
  edge_rms <- pairwise_curve_rms(curves, all_pairs[is_edge, , drop = FALSE])
  nonedge_rms <- pairwise_curve_rms(curves, all_pairs[!is_edge, , drop = FALSE])
  data.frame(method = method, n_edges = length(edge_rms), mean_edge_rms = mean(edge_rms),
             median_edge_rms = median(edge_rms), mean_nonedge_rms = mean(nonedge_rms),
             median_nonedge_rms = median(nonedge_rms), edge_to_nonedge_ratio = mean(edge_rms) / mean(nonedge_rms))
}
pair_similarity <- rbind(summarise_pair_similarity(temperature_matrix, "Observed"),
                         summarise_pair_similarity(nodewise_fit$predictions, "Nodewise"),
                         summarise_pair_similarity(network_predictions, "Network"))
pair_similarity$relative_edge_rms_reduction_vs_nodewise <- c(NA_real_, 0,
                                                             100 * (1 - pair_similarity$mean_edge_rms[3L] / pair_similarity$mean_edge_rms[2L]))
write.csv(
  pair_similarity,
  file.path(derived_dir, "dwd_pair_similarity_diagnostic.csv"),
  row.names = FALSE
)

fit_diagnostics <- data.frame(
  station_id = station_metadata$station_id, station = station_metadata$Stationsname,
  target_city = station_metadata$target_city,
  rms_network_vs_nodewise = sqrt(rowMeans((network_predictions - nodewise_fit$predictions)^2))
)
write.csv(
  fit_diagnostics,
  file.path(derived_dir, "dwd_fit_diagnostics.csv"),
  row.names = FALSE
)

message("DWD application completed successfully. Outputs written to: ", application_dir)
pair_similarity_print <- pair_similarity
numeric_columns <- vapply(pair_similarity_print, is.numeric, logical(1))

pair_similarity_print[numeric_columns] <- lapply(
  pair_similarity_print[numeric_columns],
  round,
  digits = 3
)

print(pair_similarity_print)