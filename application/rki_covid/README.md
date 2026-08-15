# RKI COVID-19 application

This directory contains the reproducible real-data application using weekly
COVID-19 case incidence in Bavarian counties. It is an illustration of the
network-weighted smoother on polygonal areal data; because no latent truth is
available, it is not an accuracy evaluation.

## Analysis design

- **Outcome:** `log(1 + weekly reported COVID-19 cases per 100,000 residents)`.
- **Time period:** 2 March 2020 to 5 March 2023 (157 Monday-starting weeks).
- **Spatial units:** the 96 Bavarian Landkreise and kreisfreie Staedte.
- **Graph:** Queen contiguity of BKG VG250-EW county polygons, with population
  denominators and boundaries dated 31 December 2020.
- **Models:** separate nodewise P-splines (`k = 15`) and the proposed network
  smoother (temporal intercept `k = 10`, node-specific deviation `k = 15`).

The script first audits the complete 400-unit German panel and graph. It then
fits the connected Bavarian induced subgraph. The restriction is explicit:
the 400-node tensor-product fit exceeded the available memory under the
current `pffr`/`mgcv` implementation.

## Reproduce

From the repository root, run:

```r
source("application/rki_covid/reproduce_rki_bavaria_application.R")
```

Required CRAN packages are `data.table`, `digest`, `devtools`, `ffm`,
`ggplot2`, `mgcv`, `sf`, `spdep`, and `tf`. The script loads the local
`netfunsmooth/` package automatically. On the first run it downloads the
current RKI case-data snapshot (approximately 400 MB).

## Outputs and reproducibility record

The script writes generated outputs below this directory:

- `data/raw/rki_cases_download.csv`: frozen local RKI input snapshot;
- `data/derived/input_manifest.csv`: source URL, file size, timestamp, and
  SHA-256 checksum for that snapshot;
- `data/derived/*_audit.csv`: data, polygon, graph, incidence, and fit audits;
- `data/derived/bavaria_application_inputs.rds` and
  `bavaria_application_results.rds`: generated inputs and fitted results;
- `data/derived/computational_provenance.rds` and `package_versions.csv`:
  R session and package-version record;
- `figures/bavaria_smoothing_difference_map.pdf` and
  `figures/bavaria_curve_comparison.pdf`: the generated figures.

The raw RKI file is excluded from Git because of its size. Do not delete it if
you need to reproduce the exact analysis later: the RKI repository rebuilds
the current CSV over time. The manifest checksum identifies the precise local
snapshot used for a run.