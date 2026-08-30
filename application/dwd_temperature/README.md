# DWD temperature application

This directory contains the reproducible real-data application using daily
mean air temperatures from the Climate Data Center of the German Weather
Service (DWD). It illustrates network-weighted smoothing for functional data
observed at geographical point locations. Because the latent temperature
curves are unknown, the application is descriptive and is not an external
accuracy evaluation.

## Analysis design

- **Outcome:** daily mean air temperature (`TMK.Lufttemperatur`).
- **Time period:** 1 January to 31 December 2024 (366 daily observations).
- **Spatial units:** 25 prespecified DWD stations with complete records on the
  common 2024 date grid.
- **Graph:** symmetrised four-nearest-neighbour graph constructed from station
  coordinates projected to ETRS89 / LAEA Europe (EPSG:3035).
- **Models:** separate nodewise P-splines (`k = 15`) and the proposed network
  smoother (temporal intercept `k = 10`, node-specific deviation `k = 15`).

The station registry is fixed in the reproduction script so that rerunning the
application does not repeat an exploratory station-selection step. Station
metadata are matched by station ID only after restricting the DWD metadata
index to the same product used in the analysis:

```text
res = daily
var = kl
per = historical
hasfile = TRUE
```

The script requires exactly one matching metadata row for every station and
stops if a station is missing, duplicated, or lacks coordinates.

## Reproduce

From the repository root, run either:

```r
source("application/dwd_temperature/reproduce_dwd_temperature_application.R")
```

or from a shell:

```sh
Rscript application/dwd_temperature/reproduce_dwd_temperature_application.R
```

Required CRAN packages are `bit64`, `devtools`, `ggplot2`, `igraph`, `rdwd`,
`sf`, and `tf`. The script loads the local `netfunsmooth/` package and the
nodewise comparator in `simulation/R/fit_baselines.R` automatically. On the
first run, `rdwd` downloads the historical daily climate files for the 25
stations.

## Outputs

Generated files are written below this directory and are excluded from Git.

### Data and audit files

- `data/derived/dwd_station_registry.csv`: fixed station IDs and display
  locations;
- `data/derived/dwd_station_metadata.csv`: product-filtered DWD station
  metadata and coordinates;
- `data/derived/dwd_temperature_2024.csv`: complete 25 by 366 temperature
  matrix;
- `data/derived/dwd_adjacency_knn4.csv`: symmetrised four-nearest-neighbour
  adjacency matrix;
- `data/derived/dwd_knn4_graph_summary.csv`: node, edge, component, degree, and
  edge-distance checks;
- `data/derived/dwd_temperature_2024_primary_fits.rds`: observed data, fitted
  curves, graph, basis settings, and fitted network model;
- `data/derived/dwd_network_smoothing_parameters.csv`: fitted network-model
  smoothing parameters;
- `data/derived/dwd_pair_similarity_diagnostic.csv`: descriptive pairwise
  curve-difference diagnostic for graph edges and non-edges;
- `data/derived/dwd_fit_diagnostics.csv`: station-level RMS differences between
  the network and nodewise fitted curves.

### Figures

- `figures/dwd_knn4_station_graph.png`: station locations and the final graph;
- `figures/dwd_curve_comparison_selected_stations.pdf`: observed, nodewise,
  and network-weighted curves for six selected stations.

The raw downloaded DWD files are retained under `data/raw/`. The curated
outputs used for the thesis are also available under
`simulation/results/thesis_release/applications/dwd/`.

## Interpretation

The pairwise curve-difference output is a regularisation diagnostic, not an
accuracy measure. Smaller differences between fitted curves are an expected
consequence of smoothing and do not by themselves establish that the network
fit is closer to an unobserved true temperature function or that any reduction
is specific to the supplied graph topology.
