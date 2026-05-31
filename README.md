# Network-weighted Functional Smoothing

Bachelor thesis project at LMU.

## Project Description

This project investigates network-weighted smoothing of functional data
using tensor-product smoothers in mgcv. The goal is to estimate
node-specific smooth functions while accounting for both temporal
smoothness and graph-based dependence.

## Methods

- Functional data representation via tidyfun
- Graph structures via igraph / spdep
- Tensor-product smoothing with
  te(node, t, bs = c("mrf", "ps"))
- REML-based smoothing parameter selection

## Repository Structure

├── R/             # helper functions
├── simulation/    # simulation studies
├── application/   # real data applications
├── results/       # generated results
└── report/        # thesis materials

## Supervisor

Prof. Dr. Fabian Scheipl
LMU Munich
