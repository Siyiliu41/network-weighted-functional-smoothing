# Network-Weighted Functional Smoothing

Bachelor thesis project in Statistics and Data Science at LMU Munich.

This repository contains an R implementation and evaluation workflow for
network-weighted smoothing of functional data. The objective is to estimate one
smooth function for each observed graph node while borrowing information across
neighbouring nodes.

## Overview

For functional observations $Y_i(t)$ at nodes $i=1,\ldots,N$ of a known
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

The completed core experiment uses a $5 \times 5$ lattice with 25 nodes and
50 time points. It varies:

* truth structure: coordinate-smooth or cluster-structured;
* structured-signal proportion: $\alpha = 0.5$ or $0.9$;
* supplied graph: rook or queen neighbourhood; and
* noise standard deviation: $\sigma = 0.2$ or $0.5$.

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
* a cluster-boundary difference-in-differences diagnostic;
* an $\alpha = 0$ zero-structured-signal diagnostic;
* a fixed degree-preserving rewired-graph control;
* an $N = 100$ computational-scaling timing gate; and
* a restricted-oracle tuning diagnostic.


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
  00_known_answer.R                    Initial known-answer check
  01_smoke_test.R                      Initial smoke test
  02_pilot_simulation.R                Pilot simulation
  03_design_audit.R                    Basis-projection design audit
  04_formal_known_answer.R             Formal known-answer sanity check
  05_main_core_simulation.R            Final 16-cell core simulation
  06_summarise_main_core_results.R     Core-result summaries and figures
  07_cluster_boundary_did.R            Cluster-boundary diagnostic
  08_run_alpha0_negative_control.R     Alpha = 0 zero-structured-signal diagnostic
  09_summarise_alpha0_negative_control.R
                                      Alpha = 0 diagnostic summaries
  10_run_rewired_graph_control.R       Fixed rewired-graph control
  11_n100_timing_pilot.R               N = 100 timing-feasibility assessment
  12_restricted_oracle_preflight.R     Restricted-oracle implementation checks
  13_run_restricted_oracle_pilot.R     Restricted-oracle diagnostic (50 replications per cell)
  results/                            Generated outputs (generally ignored by Git)
    thesis_release/                   Curated thesis outputs tracked in Git
      README.md                       Release inventory and reproduction references
      tables/                         Simulation summary tables
      figures/                        Simulation figures
      audits/                         Design-validation records
      metadata/                       Compact result objects and provenance
      applications/
        dwd/                          Retained DWD application outputs
        rki/                          Retained RKI application outputs

application/
  dwd_temperature/
    reproduce_dwd_temperature_application.R
                                      DWD temperature application
    README.md                         Application-specific instructions
  rki_covid/
    reproduce_rki_bavaria_application.R
                                      RKI Bavaria COVID-19 application
    README.md                         Application-specific instructions
thesis/
    Bachelorarbeit_Siyi_Liu.pdf
    
Simulation_study_specification.md     ADEMP specification and execution record
README.md                             Project overview and reproduction guide
.gitignore                            Excludes generated results and downloaded data
```


# Reproducibility

### Install development tools if needed
```r
install.packages("devtools")
```

### Install every package declared in netfunsmooth/DESCRIPTION
```r
devtools::install_deps(
  "netfunsmooth",
  dependencies = TRUE,
  upgrade = "never"
)
```

### Load the local development version
```r
devtools::load_all("netfunsmooth")
```

The main simulation and its summaries can be run in sequence from the
repository root:

```r
source("simulation/02_pilot_simulation.R") 
source("simulation/03_design_audit.R")
source("simulation/04_formal_known_answer.R")

source("simulation/05_main_core_simulation.R")
run_main_core_simulation("dry_run")
run_main_core_simulation("main")

source("simulation/06_summarise_main_core_results.R")
source("simulation/07_cluster_boundary_did.R")

source("simulation/08_run_alpha0_negative_control.R")
run_alpha0_negative_control("dry_run")
run_alpha0_negative_control("main")

source("simulation/09_summarise_alpha0_negative_control.R")

source("simulation/10_run_rewired_graph_control.R")
run_rewired_graph_control("dry_run")
run_rewired_graph_control("main")

source("simulation/11_n100_timing_pilot.R")
run_n100_timing_pilot()

source("simulation/12_restricted_oracle_preflight.R")
run_restricted_oracle_preflight()

source("simulation/13_run_restricted_oracle_pilot.R")
run_restricted_oracle_pilot()
```
The dry_run performs one complete 16-cell preflight replication and verifies
the required pilot and known-answer records. The full main run uses 200
attempted replications per cell and may take substantial time.

Except for the curated thesis release, generated checkpoints, fitted objects,
tables, and figures are excluded from version control. Released results can
be inspected directly. Rerunning the summary scripts requires the preceding
analysis outputs in their working directories, which are not supplied as a
complete checkpoint bundle.

The N = 100 timing pilot saves its feasibility assessment and stops if
even 25 replications per cell exceed the prespecified eight-hour budget.
This stop is an intended feasibility decision, not necessarily a fitting error.

The DWD and RKI applications can be reproduced with:

```r
source("application/dwd_temperature/reproduce_dwd_temperature_application.R")
source("application/rki_covid/reproduce_rki_bavaria_application.R")
```

The simulation scripts use a prespecified `L'Ecuyer-CMRG` random-number
stream and substream ledger, explicit quality checks, and checkpoint/resume
handling. Selected derived outputs are retained in the thesis release;
raw data and large intermediate objects are excluded. Exact reproduction
of the observed-data results requires the original input snapshots and
compatible software versions. A fresh download from the continuously
updated RKI source may differ from the snapshot used in the thesis.

The applications require additional packages listed in their respective
README files; installing the package dependencies alone may not install
all application dependencies.

## Thesis Release

The committed
[`simulation/results/thesis_release/`](simulation/results/thesis_release/)
contains a curated collection of numerical summaries, figures, audit files,
and computational metadata supporting the thesis. These files can be inspected
directly without rerunning the analyses.

The release contains:

- `tables/`: simulation summary tables;
- `figures/`: simulation figures;
- `audits/`: design-validation and diagnostic records;
- `metadata/`: compact result objects and computational provenance;
- `applications/`: selected DWD and RKI application outputs.

See the [release README](simulation/results/thesis_release/README.md)
for the file inventory and reproduction-script references. Output filenames
retain their script-level numbering and need not match the final thesis
table and figure numbers.

Raw source data, per-replication checkpoints, and large intermediate objects
are not included. Although generated result directories are generally
ignored by Git, the curated files under `thesis_release/` are explicitly
tracked. The release allows inspection of the retained results but is not
a complete checkpoint bundle for rerunning every analysis.

## Thesis

[Read the thesis (PDF)](thesis/Bachelorarbeit_Siyi_Liu.pdf)

Siyi Liu, *Network-Weighted Smoothing of Functional Data*.
Bachelor's thesis, Statistics and Data Science, LMU Munich, 2026.
Supervisor: Prof. Dr. Fabian Scheipl.
