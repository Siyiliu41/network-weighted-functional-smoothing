# Simulation Study of Network-Weighted Smoothing for Functional Data

## Specification According to the ADEMP Framework

**Status:** Prespecified before conducting and inspecting the main simulation results  
**Framework:** ADEMP framework proposed by Morris, White, and Crowther (2019)

## 1. Overview

The simulation study evaluates the statistical performance of the network-weighted smoother implemented in the R package `netfunsmooth`.

For every node $i=1,\ldots,N$ in a known graph, a noisy functional observation

$$
Y_i(t_j)=f_i(t_j)+\varepsilon_{ij}
$$

is available. The objective is to reconstruct the unknown smooth node-specific curves $f_i(t)$. The study focuses on whether simultaneous smoothing over the functional domain and the graph improves reconstruction compared with methods that do not use the graph structure.

The study follows the ADEMP framework:

1. **Aims**
2. **Data-generating mechanisms**
3. **Estimands**
4. **Methods**
5. **Performance measures**

---

## 2. Aims

### 2.1 Primary aim

The primary aim is to determine whether incorporating a known graph structure improves the reconstruction of node-specific functions.

The central research question is:

> Under which combinations of graph-signal strength, curve structure, neighbourhood density, and observation noise does network-weighted smoothing achieve a lower reconstruction error than nodewise smoothing without graph information?

### 2.2 Secondary aims

The simulation study additionally addresses the following questions:

1. Does the advantage of network-weighted smoothing increase with graph-signal strength?
2. Is graph information particularly beneficial when the observations are strongly contaminated by noise?
3. How does neighbourhood density affect estimation accuracy and potential oversmoothing?
4. Can the method reconstruct both graph-smooth and block- or cluster-structured curves?
5. Does the reconstruction error increase near boundaries between clusters?
6. How well does REML select the temporal and graph smoothing parameters compared with an oracle choice?
7. Does more accurate curve reconstruction improve a downstream clustering task?
8. How frequently do numerical warnings, convergence problems, or failed fits occur?

### 2.3 Expected qualitative results

The following qualitative results are expected:

- Under strong graph structure and high observation noise, network-weighted smoothing should outperform nodewise smoothing.
- Under weak or absent graph structure, the advantage should become smaller or disappear.
- Pooled mean smoothing should only be competitive when the node-specific curves are nearly identical.
- In block-structured scenarios, graph smoothing may reduce variance within clusters but oversmooth across genuine cluster boundaries.

These expectations guide the interpretation of the results but will not be used to select or exclude simulation outcomes.

---

## 3. Data-Generating Mechanisms

### 3.1 Graph and observation grid

The core simulation uses a regular $5\times5$ lattice with

$$
N=25
$$

nodes.

Each function is observed on an equally spaced grid

$$
t_j\in[0,1],
\qquad
j=1,\ldots,T,
\qquad
T=50.
$$

Within each simulation replication, all methods are applied to exactly the same simulated dataset. This allows paired comparisons between methods.

### 3.2 Neighbourhood density

Two neighbourhood structures are considered.

#### Sparse graph

Rook adjacency is used: two nodes are considered neighbours if their grid cells share an edge.

This graph has 40 edges and an average node degree of 3.2.

#### Dense graph

Queen adjacency is used: in addition to shared edges, grid cells that share a corner are also considered neighbours.

This graph has 72 edges and an average node degree of 5.76.

Because the scale of the graph Laplacian changes with the number of edges and node degrees, a scaled or normalised graph Laplacian is used when constructing and comparing graph signals. This prevents apparent density effects from being caused only by a change in the overall penalty scale.

### 3.3 Observation model

The observed functions are generated as

$$
Y_i(t_j)
=
f_i(t_j)+\varepsilon_{ij},
$$

where, in the core simulation,

$$
\varepsilon_{ij}
\overset{\mathrm{iid}}{\sim}
N(0,\sigma^2).
$$

Gaussian and temporally independent errors provide a deliberately simple reference setting. Temporally correlated or heteroscedastic errors may subsequently be considered in robustness analyses.

### 3.4 True node-specific curves

The true function at node $i$ is decomposed as

$$
f_i(t)=\mu(t)+\delta_i(t),
$$

where $\mu(t)$ is a common smooth mean function and $\delta_i(t)$ is a node-specific deviation.

For identifiability, the deviations are centred across nodes for every $t$:

$$
\sum_{i=1}^{N}\delta_i(t)=0.
$$

The node-specific deviations are constructed using multiple temporal basis functions:

$$
\delta_i(t)
=
\sum_{k=1}^{K}\theta_{ik}\phi_k(t).
$$

The core simulation uses $K=3$, for example:

$$
\phi_1(t)=\sin(2\pi t),
$$

$$
\phi_2(t)=\cos(2\pi t),
$$

and

$$
\phi_3(t)
=
\exp\left\{-100(t-0.65)^2\right\}.
$$

Varying multiple coefficients $\theta_{ik}$ allows the node-specific functions to differ in amplitude, phase, curvature, and local peaks.

This is essential because the simulation should test whether the estimator can recover differences in curve shape across the graph, rather than only vertical shifts or amplitude differences.

### 3.5 Graph-signal strength

Graph-signal strength describes how strongly the similarity between the true node-specific curves agrees with the specified graph structure.

For each temporal basis coefficient, define

$$
\boldsymbol{\theta}_k
=
\sqrt{\alpha}\,
\boldsymbol{\theta}^{(G)}_k
+
\sqrt{1-\alpha}\,
\boldsymbol{\theta}^{(I)}_k,
$$

where:

- $\boldsymbol{\theta}^{(G)}_k$ is a graph-compatible coefficient vector;
- $\boldsymbol{\theta}^{(I)}_k$ is a coefficient vector independent of the graph;
- $\alpha$ controls the graph-signal strength.

After combining the two components, the coefficients are centred and standardised so that their marginal signal variance remains comparable across scenarios.

The core simulation considers

$$
\alpha\in\{0.5,0.9\},
$$

representing moderate and strong graph signals.

A separate negative-control scenario uses

$$
\alpha=0.
$$

In this setting, the graph contains no information about the similarity of the true curves. This scenario evaluates whether the graph penalty causes harmful oversmoothing when the graph structure is irrelevant.

Graph smoothness is additionally described using normalised Laplacian energy:

$$
R_k
=
\frac{
\boldsymbol{\theta}_k^\top L_{\mathrm{norm}}
\boldsymbol{\theta}_k
}{
\boldsymbol{\theta}_k^\top\boldsymbol{\theta}_k
}.
$$

Small values of $R_k$ indicate a strong and smooth graph signal.

### 3.6 Structure of the true curves

Two main data-generating mechanisms are considered.

#### DGP 1: Graph-smooth truth

The graph-compatible coefficient vectors $\boldsymbol{\theta}^{(G)}_k$ are constructed from low-frequency eigenvectors of the normalised graph Laplacian.

As a result, the shapes of the functions change gradually across the graph. Neighbouring nodes have similar but not identical curves.

This scenario corresponds closely to the main assumption of the network-weighted estimator. It also serves as a sanity check of whether the graph penalty has its intended effect.

#### DGP 2: Block- or cluster-structured truth

The graph nodes are divided into several connected clusters. Nodes within a cluster have similar coefficients, whereas coefficients differ substantially between clusters.

The coefficients can be represented as

$$
\theta_{ik}
=
\gamma_{c(i),k}+u_{ik},
$$

where:

- $c(i)$ denotes the true cluster assignment of node $i$;
- $\gamma_{c(i),k}$ is a cluster-specific coefficient;
- $u_{ik}$ is a smaller node-specific deviation.

This produces similar curves within clusters and discontinuous changes at cluster boundaries.

The scenario investigates the trade-off between:

- variance reduction through information sharing within clusters; and
- potential oversmoothing across genuine cluster boundaries.

A node is classified as a boundary node if at least one of its neighbours belongs to a different true cluster.

### 3.7 Noise level

The core simulation uses two noise levels:

$$
\sigma\in\{0.2,0.5\}.
$$

These represent low and high observation noise. A pilot simulation will verify that the two levels produce meaningfully different levels of difficulty, both visually and in terms of the signal-to-noise ratio.

The empirical signal-to-noise ratio is documented as

$$
\operatorname{SNR}
=
\frac{
\frac{1}{NT}
\sum_{i=1}^{N}
\sum_{j=1}^{T}
\left(f_i(t_j)-\bar f\right)^2
}{
\sigma^2
}.
$$

If the selected values do not produce a meaningful distinction between easy and difficult scenarios, they may be adjusted based on the pilot simulation. Any adjustment must be documented before the main simulation results are inspected.

### 3.8 Core simulation scenarios

The core factorial design combines the following factors:

| Factor | Levels |
|---|---|
| Structure of the truth | Graph-smooth, cluster-structured |
| Graph-signal strength | $\alpha=0.5,0.9$ |
| Neighbourhood density | Sparse, dense |
| Noise level | $\sigma=0.2,0.5$ |

This results in

$$
2\times2\times2\times2=16
$$

core scenarios.

The scenarios with $\alpha=0$ are evaluated as additional negative controls. They need not be combined with every other design factor if computational resources are limited.

### 3.9 Sensitivity analyses

The following sensitivity analyses will be conducted for selected representative scenarios:

1. larger graphs, such as $N=100$ on a $10\times10$ lattice;
2. different temporal basis dimensions for the mean function and node-specific deviations;
3. comparison of a deliberately low basis dimension with sufficiently flexible specifications;
4. temporally correlated observation errors;
5. mildly misspecified graphs created by adding or removing selected edges.

The temporal basis dimension for the node-specific deviations should be larger than the low default value $k=5$.

A possible core specification is

$$
K_{\mathrm{int}}=K_{\mathrm{dev}}=10,
$$

with smaller and larger basis dimensions examined in selected sensitivity scenarios.

---

## 4. Estimands

### 4.1 Primary estimand

The primary estimand is the complete collection of true node-specific functions:

$$
\mathcal{F}
=
\left\{
f_1(t),\ldots,f_N(t)
\right\}.
$$

Each method produces estimates

$$
\widehat{\mathcal{F}}
=
\left\{
\widehat f_1(t),\ldots,\widehat f_N(t)
\right\}.
$$

### 4.2 Secondary estimands

Secondary estimands include:

- reconstruction accuracy at individual nodes;
- average reconstruction accuracy across all nodes;
- reconstruction accuracy at cluster boundaries;
- reconstruction accuracy at interior cluster nodes;
- the true cluster structure in the block-structured DGP;
- the reconstruction error attainable under an optimal choice of smoothing parameters.

### 4.3 Distinction from spatial prediction

The primary simulation study concerns the reconstruction of noisy functions at observed graph nodes.

This differs from functional kriging, where the usual objective is to predict an entire curve at an unobserved spatial location. A comparison with functional kriging is therefore treated as a separate node-holdout analysis.

---

## 5. Methods

### M1: Network-weighted smoothing with REML

The primary method under investigation is `netf_smooth()` from the `netfunsmooth` package.

The estimator uses a tensor-product penalised representation with:

- a temporal smoothness penalty;
- an MRF penalty based on the graph Laplacian;
- REML-based selection of the smoothing parameters.

The underlying `pffr` specification contains a functional mean and a node-by-time interaction:

```r
Y ~ s(
  node,
  bs = "mrf",
  xt = list(nb = nb_list),
  k = N
)
```

Internally, this is expanded into a structure of the form

```r
s(yindex.vec) + ti(node, yindex.vec, ...)
```

The specification uses:

- the standard functional intercept;
- no `c()` wrapper around the MRF term;
- no `-1`.

This allows the node-specific functions to differ in shape rather than only in level.

### M2: Network-weighted oracle smoothing

As a diagnostic reference, the same network-weighted model is fitted with smoothing parameters selected over a prespecified grid to minimise the true reconstruction error.

The oracle method is not available in real applications because it uses the unknown true curves. It is included only to answer the following questions:

- How well does REML select the smoothing parameters?
- Is poor performance caused by the model itself or by smoothing-parameter selection?
- How large is the difference between practically attainable and theoretically optimal performance?

Because of its computational cost, the oracle method may be restricted to selected scenarios.

### M3: Nodewise univariate smoothing

Each node-specific function is smoothed separately over $t$:

$$
Y_i(t)=s_i(t)+\varepsilon_i(t).
$$

This method does not use graph information but preserves node-specific variation.

It is the most important comparator because it directly assesses whether the graph penalty provides an improvement beyond ordinary temporal smoothing.

Where possible, M1 and M3 use comparable temporal bases and the same type of smoothing-parameter selection.

### M4: Pooled mean smoothing

Node identity is ignored, and a single common mean function is estimated:

$$
Y_i(t)=\mu(t)+\varepsilon_i(t).
$$

The same estimated function is then assigned to every node:

$$
\widehat f_i(t)=\widehat\mu(t).
$$

This method represents complete pooling. It may have low variance but cannot represent systematic differences between node-specific curves.

### M5: Functional kriging with `SpatFD`

Functional kriging is not included as a direct core comparator for reconstruction at observed nodes because its primary purpose is prediction at unobserved spatial locations.

Instead, a separate leave-one-node-out analysis is planned:

1. One node is removed from the training data.
2. The method is fitted using the remaining nodes.
3. The complete function at the omitted node is predicted.
4. The procedure is repeated over several or all nodes.

This comparison requires spatial coordinates for the graph nodes. It assesses spatial prediction and will therefore be reported separately from the primary reconstruction benchmark.

---

## 6. Performance Measures

### 6.1 Node-specific integrated squared error

For node $i$, method $m$, and simulation replication $b$, the integrated squared error is

$$
\operatorname{ISE}_{ibm}
=
\int_0^1
\left[
\widehat f_{ibm}(t)-f_{ib}(t)
\right]^2dt.
$$

On the discrete observation grid, the integral is approximated numerically, preferably using the trapezoidal rule.

### 6.2 Node-averaged ISE

Within each simulation replication, the ISE is averaged across nodes:

$$
\operatorname{AISE}_{bm}
=
\frac{1}{N}
\sum_{i=1}^{N}
\operatorname{ISE}_{ibm}.
$$

Across $B$ simulation replications, the estimated mean integrated squared error is

$$
\widehat{\operatorname{MISE}}_m
=
\frac{1}{B}
\sum_{b=1}^{B}
\operatorname{AISE}_{bm}.
$$

This is the primary absolute performance measure.

### 6.3 Relative improvement over nodewise smoothing

The main comparative measure is defined relative to M3:

$$
\operatorname{RI}_{bm}
=
100
\frac{
\operatorname{AISE}_{b,\mathrm{nodewise}}
-
\operatorname{AISE}_{bm}
}{
\operatorname{AISE}_{b,\mathrm{nodewise}}
}.
$$

Its interpretation is:

- $\operatorname{RI}>0$: method $m$ performs better than nodewise smoothing;
- $\operatorname{RI}=0$: no performance difference;
- $\operatorname{RI}<0$: method $m$ performs worse.

The MISE ratio is additionally reported:

$$
\operatorname{rMISE}_m
=
\frac{
\widehat{\operatorname{MISE}}_m
}{
\widehat{\operatorname{MISE}}_{\mathrm{nodewise}}
}.
$$

Values below 1 favour the method under consideration.

### 6.4 Reconstruction error at cluster boundaries

For the block-structured DGP, performance is separately evaluated for boundary and interior nodes:

$$
\operatorname{MISE}_{\mathrm{boundary}}
$$

and

$$
\operatorname{MISE}_{\mathrm{interior}}.
$$

The difference

$$
\operatorname{MISE}_{\mathrm{boundary}}
-
\operatorname{MISE}_{\mathrm{interior}}
$$

is used as an indicator of potential oversmoothing across genuine cluster boundaries.

### 6.5 Oracle gap

For the REML-based network estimator, the difference from oracle performance is

$$
\operatorname{OracleGap}
=
\widehat{\operatorname{MISE}}_{\mathrm{network,REML}}
-
\widehat{\operatorname{MISE}}_{\mathrm{network,oracle}}.
$$

The ratio between estimated and oracle smoothing parameters may additionally be examined.

### 6.6 Monte Carlo uncertainty

A Monte Carlo standard error is reported for every estimated performance measure.

For the estimated MISE, for example,

$$
\operatorname{MCSE}
=
\frac{
\operatorname{SD}
\left(
\operatorname{AISE}_{1m},\ldots,
\operatorname{AISE}_{Bm}
\right)
}{
\sqrt{B}
}.
$$

The simulation results are therefore reported together with their Monte Carlo uncertainty rather than only as point estimates.

### 6.7 Numerical stability and computational cost

The following quantities are additionally recorded for each method:

- mean and median computation time;
- proportion of fits producing warnings;
- proportion of failed fits;
- proportion of non-converged fits;
- number of extreme or implausible estimates, where applicable.

Failed fits are not silently removed. Performance measures are reported for successful fits, while the numbers of successful and unsuccessful fits are presented separately.

Paired comparisons are based on replications in which both methods being compared fitted successfully. If failure rates are non-negligible, the potential effect of this restriction will be examined.

---

## 7. Downstream Task: Clustering

For the cluster-structured DGP, a secondary analysis examines whether the estimated functions preserve the true cluster structure.

The procedure is:

1. Represent all estimated functions using a common basis or common functional principal components.
2. Apply the same clustering algorithm and feature representation to the results of every smoothing method.
3. Compare the estimated clusters with the known true cluster labels.

Performance is assessed using the Adjusted Rand Index:

$$
\operatorname{ARI}\in[-1,1].
$$

A higher ARI indicates that the estimated functions preserve the true cluster structure more accurately.

The number of clusters is treated as known in the core analysis. This ensures that the evaluation primarily reflects the quality of curve reconstruction rather than the separate problem of selecting the number of clusters.

An outlier-detection analysis would require an additional DGP with explicitly generated anomalous nodes and is therefore considered an optional extension.

---

## 8. Monte Carlo Design

### 8.1 Number of replications

A pilot simulation with

$$
B_{\mathrm{pilot}}=20\text{--}50
$$

replications will first be conducted.

The pilot simulation is used only to verify:

- correct data and curve generation;
- the intended graph dependence;
- numerical stability;
- computation time;
- meaningful separation of the noise levels;
- correct implementation of the performance measures.

The initial target for the core simulation is

$$
B=500
$$

replications per scenario.

The final number of replications will be based on the resulting Monte Carlo standard errors. If the Monte Carlo uncertainty remains too large at $B=500$, the number of replications will be increased.

### 8.2 Random numbers and reproducibility

A unique and reproducible random seed is assigned to every scenario and replication.

Within each replication, all methods use:

- the same true functions;
- the same realisation of the observation errors;
- the same observed dataset.

Scenario definition, data generation, model fitting, performance evaluation, and result aggregation are implemented as separate functions.

Intermediate results are stored by scenario so that interrupted simulation runs can be resumed without repeating completed scenarios.

---

## 9. Analysis and Presentation

### 9.1 Main tables

For every combination of DGP, graph-signal strength, neighbourhood density, and noise level, the following are reported:

- estimated MISE;
- MCSE of the estimated MISE;
- relative improvement over nodewise smoothing;
- MISE ratio;
- computation time;
- number of successful fits;
- warning and failure rates.

### 9.2 Main figures

The planned figures include:

1. an interaction plot of relative improvement by graph-signal strength and noise level;
2. a heatmap of the MISE ratio between network-weighted and nodewise smoothing;
3. boxplots of replication-specific AISE values;
4. separate reconstruction errors for boundary and interior nodes;
5. selected example curves showing the truth, noisy observations, and estimates from the comparison methods;
6. clustering performance based on the Adjusted Rand Index.

The figures will show both average performance and variability across simulation replications.

---

## 10. Decisions to Be Finalised Before the Main Simulation

The following decisions will be finalised before the main simulation is conducted:

1. the exact values of the two noise levels;
2. the exact construction and standardisation of the coefficients under both DGPs;
3. the temporal basis dimensions for the functional mean and node-by-time deviations;
4. the smoothing-parameter grid for the oracle method;
5. whether the oracle method is evaluated in all or only selected scenarios;
6. the final number of Monte Carlo replications;
7. the clustering algorithm and functional representation;
8. the exact scope of the separate `SpatFD` experiment;
9. the sensitivity scenarios selected for $N=100$.

These decisions will be documented before the main simulation results are inspected.

---

## 11. Prioritisation Under Computational Constraints

If computational time is limited, the analyses will be prioritised as follows:

1. **Core benchmark:** network-REML versus nodewise smoothing;
2. **Sanity check:** a strongly graph-smooth setting in which graph smoothing should clearly help;
3. **Negative control:** $\alpha=0$;
4. **Pooled baseline;**
5. **Cluster-boundary and clustering analysis;**
6. **Oracle diagnostics;**
7. **Larger graphs and basis-dimension sensitivity;**
8. **Leave-one-node-out comparison with `SpatFD`;**
9. **Further robustness analyses.**

This prioritisation ensures that the main research questions can still be answered if not all optional analyses can be completed.

---

## 12. Summary of the Core Simulation

| ADEMP component | Specification |
|---|---|
| Aim | Assess the benefit of graph information for reconstructing node-specific curves |
| Data | $N=25$, $5\times5$ lattice, $T=50$, Gaussian noise |
| Truth | Graph-smooth or cluster-structured, with genuine shape differences across nodes |
| Graph signal | $\alpha=0.5,0.9$, with $\alpha=0$ as a negative control |
| Graph density | Rook and queen adjacency |
| Noise level | $\sigma=0.2,0.5$, to be validated in the pilot simulation |
| Methods | Network-REML, network-oracle, nodewise smoothing, pooled smoothing |
| Primary estimand | True node-specific functions $f_i(t)$ |
| Primary measure | Node-averaged ISE and relative improvement over nodewise smoothing |
| Secondary measures | Boundary MISE, oracle gap, ARI, computation time, warning and failure rates |
| Core design | 16 scenarios, initially $B=500$ replications |
| `SpatFD` | Separate leave-one-node-out spatial prediction experiment |

## Reference

Morris, T. P., White, I. R., and Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38, 2074–2102. <https://doi.org/10.1002/sim.8086>
