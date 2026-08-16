# Network-Weighted Functional Smoothing

Bachelor thesis project in Statistics and Data Science at LMU Munich.

This repository contains an R implementation and evaluation workflow for
network-weighted smoothing of functional data. The objective is to estimate one
smooth function for each observed graph node while borrowing information across
neighbouring nodes.

## Overview

For functional observations (Y_i(t)) at nodes (i=1,\ldots,N) of a known
undirected graph, the method combines:

* smoothness over the functional domain;
* node-specific temporal deviations; and
* graph-based regularisation that encourages neighbouring nodes to have
  similar estimated curves.

The estimator targets reconstruction at observed graph nodes. It is not a
method for predicting functions at unobserved spatial locations.

## Method

The implementation is provided by the `netfunsmooth` package in this
repository and uses the `mgcv` / `refund::pffr` infrastructure.

The fitted model includes a shared temporal smooth and a node-by-time
interaction with a Markov random-field (MRF) graph penalty. In the underlying
`pffr` representation, the interaction has the form

```r
ti(
  node,
  yindex.vec,
  bs = c("mrf", "ps"),
  xt = list(nb = nb_list)
)
```

where `nb_list` is the supplied graph neighbourhood list. Smoothing parameters
are selected by REML.

Functional observations are handled using the `tf` / tidyfun ecosystem.
Graphs can be supplied as adjacency matrices, `igraph` objects, or supported
spatial polygon objects through `graph_to_nb()`.

## Simulation Study

The simulation study follows the ADEMP framework.

The completed core experiment uses a (5 \times 5) lattice with 25 nodes and
50 time points. It varies:

* truth structure: coordinate-smooth or cluster-structured;
* structured-signal proportion: (\alpha = 0.5) or (0.9);
* supplied graph: rook or queen neighbourhood; and
* noise standard deviation: (\sigma = 0.2) or (0.5).

This gives 16 core design cells, each with 200 attempted replications.
Coordinate-smooth and cluster-structured truths are generated independently of
the supplied estimation graph, so the rook-versus-queen comparison changes the
estimator rather than the target curves.

The primary comparison is network-weighted smoothing versus separate nodewise
smoothing. Raw observations and pooled smoothing are included as diagnostic
benchmarks. Performance is assessed using integrated squared error, node-
averaged ISE, aggregate relative MISE, paired improvements, paired Monte Carlo
standard errors, win rates, and fit diagnostics.

The completed supplementary analyses include:

* a known-answer sanity check;
* a cluster-boundary difference-in-differences diagnostic; and
* an (\alpha=0) uninformative-graph diagnostic.


For the full prespecified design and execution record, see
[`Simulation_study_specification.md`](Simulation_study_specification.md).

## Real-Data Applications

Two reproducible application workflows are provided:

* `application/dwd_temperature/reproduce_dwd_temperature_application.R`
  Daily mean DWD temperatures for 25 German stations in 2024.

* `application/rki_covid/reproduce_rki_bavaria_application.R`
  Weekly COVID-19 incidence curves for Bavarian districts.

These applications illustrate fitting and descriptive behaviour on observed
functional data. Because the latent true curves are unknown, they do not
provide external accuracy validation or establish a topology-specific effect.

## Repository Structure

```text
netfunsmooth/                         Installable R package
  R/                                  Estimator, graph conversion, and prediction methods
  tests/                              Automated package tests
  DESCRIPTION                          Package metadata and dependencies

simulation/
  R/                                  Graph construction, DGP, fitting, metrics, and RNG helpers
  02_pilot_simulation.R               Pilot simulation
  03_design_audit.R                   Basis-projection design audit
  04_formal_known_answer.R            Known-answer sanity check
  05_main_core_simulation.R           Final 16-cell core simulation
  06_summarise_main_core_results.R    Core-result summaries and Figure 1
  07_cluster_boundary_did.R           Cluster-boundary diagnostic
  08_run_alpha0_negative_control.R    Alpha = 0 diagnostic
  09_summarise_alpha0_negative_control.R
                                      Alpha = 0 summary
  results/                            Generated outputs (ignored by Git)

application/
  dwd_temperature/
    reproduce_dwd_temperature_application.R
                                      DWD temperature application
    README.md                         Application-specific instructions
  rki_covid/
    reproduce_rki_bavaria_application.R
                                      RKI Bavaria COVID-19 application
    README.md                         Application-specific instructions

Simulation_study_specification.md     ADEMP specification and execution record
README.md                             Project overview and reproduction guide
.gitignore                            Excludes generated results and downloaded data
```


## Reproducibility

# Install development tools if needed
```r
install.packages("devtools")
```

# Install every package declared in netfunsmooth/DESCRIPTION
```r
devtools::install_deps(
  "netfunsmooth",
  dependencies = TRUE,
  upgrade = "never"
)
```

# Load the local development version
```r
devtools::load_all("netfunsmooth")
```

The main simulation and its summaries can be run in sequence from the
repository root:

```r
source("simulation/05_main_core_simulation.R")
source("simulation/06_summarise_main_core_results.R")
source("simulation/07_cluster_boundary_did.R")
source("simulation/08_run_alpha0_negative_control.R")
source("simulation/09_summarise_alpha0_negative_control.R")
```

The DWD and RKI applications can be reproduced with:

```r
source("application/dwd_temperature/reproduce_dwd_temperature_application.R")
source("application/rki_covid/reproduce_rki_bavaria_application.R")
```

The simulation scripts use a prespecified `L'Ecuyer-CMRG` random-number
stream and substream ledger, explicit quality checks, and checkpoint/resume
handling. Generated simulation results, downloaded raw data, derived data
objects, and application figures are intentionally excluded from version
control; they can be recreated by the committed scripts.

## Supervisor

Prof. Dr. Fabian Scheipl
LMU Munich
