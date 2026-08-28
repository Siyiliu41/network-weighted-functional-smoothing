# Thesis reproducibility release



This directory is a curated, lightweight release of the numerical outputs,

figures, audit files, and computational metadata supporting the empirical

claims in the thesis. It is intended for inspection without rerunning the full

simulation study.



Raw source data, per-replication checkpoint files, and large intermediate

objects are intentionally excluded.



## Contents



### `tables/`



Final simulation tables cited in the thesis:



- `table_1_method_performance.csv`: primary per-cell comparison of the

  network-weighted and nodewise estimators.

- `table_2_network_vs_nodewise.csv`: aggregated network-versus-nodewise

  comparison.

- `table_3_rook_vs_queen_network.csv`: comparison of results under rook and

  queen neighbourhood graphs.

- `table_4_cluster_boundary_did.csv`: cluster-boundary difference-in-differences

  diagnostic.

- `table_5_alpha0_negative_control.csv`: negative-control results for

  zero graph-signal strength.

- `alpha0_network_vs_nodewise_detail.csv`: detailed results underlying the

  alpha-zero negative-control summary.



### `audits/`



Design-validation and diagnostic outputs:



- `basis_rpe_summary.csv` and `basis_rpe_detail.csv`: unpenalised temporal

  basis approximation checks.

- `snr_overall_summary.csv`, `snr_deviation_summary.csv`, `snr_detail.csv`,

  and `snr_separation.csv`: signal-to-noise diagnostics for the simulation

  design.

- `code_fingerprint.csv`: fingerprints of the code used when the design audit

  was created.



### `figures/`



Final simulation figures:



- `figure_1_network_relative_improvement.png`

- `figure_2_method_mise_by_cell.png`

- `figure_3_cluster_boundary_did.png`



### `metadata/`



Compact R objects supporting the reported numerical summaries and execution

environment:



- `design_audit_result.rds`

- `formal_known_answer_result.rds`

- `main_core_analysis_result.rds`

- `alpha0_negative_control_analysis_result.rds`

- `rewired_graph_control_result.rds`

- `n100_timing_pilot_result.rds`

- `restricted_oracle_50id_result.rds`

- `session_info.txt`: saved R session information from the retained simulation

  result objects.



### `applications/`



Audits, diagnostics, fitted-parameter summaries, and final figures for the

two descriptive applications. Source data are not included.



- `applications/dwd/`: German Weather Service temperature application.

- `applications/rki/`: Bavarian COVID-19 incidence application.



## Reproduction scripts



The main simulation workflow is implemented in:



- `simulation/03_design_audit.R`

- `simulation/04_formal_known_answer.R`

- `simulation/05_main_core_simulation.R`

- `simulation/06_summarise_main_core_results.R`

- `simulation/07_cluster_boundary_did.R`

- `simulation/08_run_alpha0_negative_control.R`

- `simulation/09_summarise_alpha0_negative_control.R`

- `simulation/10_run_rewired_graph_control.R`

- `simulation/11_n100_timing_pilot.R`

- `simulation/12_restricted_oracle_preflight.R`

- `simulation/13_run_restricted_oracle_pilot.R`



The scripts `simulation/00_known_answer.R`, `simulation/01_smoke_test.R`, and

`simulation/02_pilot_simulation.R` provide earlier validation and pilot runs.

Shared simulation functions are in `simulation/R/`.



The release contains final derived outputs only. Reproducing all simulations

from scratch requires the package dependencies and can take substantial

computing time.

