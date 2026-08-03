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

> Under which combinations of structured-signal proportion, curve structure, neighbourhood density, and observation noise does network-weighted smoothing achieve a lower reconstruction error than nodewise smoothing without graph information?

### 2.2 Secondary aims

The simulation study additionally addresses the following questions:

1. Does the advantage of network-weighted smoothing increase with structured-signal proportion?
2. Is graph information particularly beneficial when the observations are strongly contaminated by noise?
3. How does neighbourhood density affect estimation accuracy and potential oversmoothing?
4. Can the method reconstruct both coordinate-smooth and block- or cluster-structured curves?
5. Does the reconstruction error increase near boundaries between clusters?
6. In selected representative cells, how closely does REML approach restricted oracle smoothing-parameter selection?
7. Does more accurate curve reconstruction improve a downstream clustering task?
8. How frequently do numerical warnings, convergence problems, or failed fits occur?
9. What is the reconstruction cost of supplying an actively misleading graph rather than the correct graph?

The effect of neighbourhood density is evaluated primarily through the direct paired rook-versus-queen difference in AISE within otherwise identical design cells.

The cost of graph misspecification is evaluated through the direct paired rewired-versus-correct-rook difference in AISE within otherwise identical false-graph control cells.

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
t_j\in[0,1], \qquad j=1,\ldots,T, \qquad T=50.
$$

The formal core simulation and the prespecified known-answer sanity check in Section 3.8 use $T=50$. Existing legacy known-answer and smoke-test scripts use $T=30$ only as software-validation checks; the formal simulation drivers explicitly use the prespecified 50-point grid.

Within each simulation replication, all methods are applied to exactly the same simulated dataset. This allows paired comparisons between methods.

### 3.2 Neighbourhood density

Two neighbourhood structures are considered.

#### Sparse graph

Rook adjacency is used: two nodes are considered neighbours if their grid cells share an edge.

This graph has 40 edges and an average node degree of 3.2.

#### Dense graph

Queen adjacency is used: in addition to shared edges, grid cells that share a corner are also considered neighbours.

This graph has 72 edges and an average node degree of 5.76.

Neighbourhood density is varied only through the graph supplied to the estimator. It is not used to generate the true coefficient fields.

For every combination of truth structure, signal strength, noise level, and replication, the true curves and observation errors are generated once and reused for both the rook and queen fits. Thus,

$$
f_i^{\mathrm{rook}}(t)=f_i^{\mathrm{queen}}(t)
$$

and both fits receive the same noisy observations. Only the neighbourhood graph passed to the network-weighted estimator changes.

Although the truth-generation function may use an `igraph` object as a container for node names and grid coordinates, its adjacency structure is not used to generate the truth. Consequently, the paired rook-versus-queen comparison isolates the effect of the neighbourhood structure used by the estimator.

### 3.3 Observation model

The observed functions are generated as

$$
Y_i(t_j) = f_i(t_j)+\varepsilon_{ij},
$$

where, in the core simulation,

$$
\varepsilon_{ij} \overset{\mathrm{iid}}{\sim} N(0,\sigma^2).
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
\delta_i(t) = \sum_{k=1}^{K}\theta_{ik}\phi_k(t).
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
\phi_3(t) = \exp\{-100(t-0.65)^2\}.
$$

Varying multiple coefficients $\theta_{ik}$ allows the node-specific functions to differ in amplitude, phase, curvature, and local peaks.

This is essential because the simulation should test whether the estimator can recover differences in curve shape across the graph, rather than only vertical shifts or amplitude differences.

### 3.5 Structured-signal strength and graph compatibility

The true coefficient fields are generated independently of the rook or queen graph supplied to the estimator. For both core DGPs, the coefficient vector of temporal basis component $k$ is constructed as

$$
\boldsymbol{\theta}_k =
s_k\left[
\sqrt{\alpha}\,\boldsymbol{\theta}^{(S)}_k +
\sqrt{1-\alpha}\,\boldsymbol{\theta}^{(U)}_k
\right],
$$

where:

- $\boldsymbol{\theta}^{(S)}_k$ is the DGP-specific structured component;
- $\boldsymbol{\theta}^{(U)}_k$ is an unstructured component generated independently of the fitted graph;
- $s_k$ is a prespecified coefficient-scale factor;
- $\alpha$ is the proportion of coefficient variance assigned to the structured component.

Before mixing, both components are centred, scaled to unit empirical variance, and made orthogonal:

$$
\mathbf{1}^{\top}\boldsymbol{\theta}^{(S)}_k = 0,
\qquad
\mathbf{1}^{\top}\boldsymbol{\theta}^{(U)}_k = 0,
$$

$$
\frac{1}{N}\left\lVert\boldsymbol{\theta}^{(S)}_k\right\rVert_2^2 = 1,
\qquad
\frac{1}{N}\left\lVert\boldsymbol{\theta}^{(U)}_k\right\rVert_2^2 = 1,
\qquad
\left(\boldsymbol{\theta}^{(S)}_k\right)^{\top}
\boldsymbol{\theta}^{(U)}_k = 0.
$$

Under this construction, $\alpha$ has an exact variance-partition interpretation while the marginal coefficient variance remains fixed across signal-strength scenarios.

The core simulation uses

$$
\alpha\in\{0.5,0.9\},
$$

representing moderate and strong structured signals. An additional negative control uses $\alpha=0$, for which the coefficient variation is unrelated to the grid-based structured component.

Graph compatibility is measured after truth generation for each graph $G$, using standard graph-Laplacian energy definitions (Chung, 1997). Let $A_G$ denote its adjacency matrix, $D_G$ its degree matrix, and

$$
L_G=D_G-A_G
$$

the unnormalised Laplacian used by the MRF penalty. The corresponding energy is

$$
R^{\mathrm{MRF}}_{k,G} =
\frac{
\boldsymbol{\theta}_k^{\top}L_G\boldsymbol{\theta}_k
}{
\boldsymbol{\theta}_k^{\top}\boldsymbol{\theta}_k
}.
$$

Because this is the Laplacian associated with the neighbourhood list used by `bs = "mrf"`, $R^{\mathrm{MRF}}_{k,G}$ is the primary penalty-alignment diagnostic.

For a degree-adjusted descriptive comparison between graph densities, the normalised Laplacian

$$
L_{\mathrm{norm},G} =
I-D_G^{-1/2}A_GD_G^{-1/2}
$$

and its corresponding energy $R^{\mathrm{norm}}_{k,G}$ are additionally reported. Neither Laplacian is used to generate the true curves.

### 3.6 Structure of the true curves

Two main data-generating mechanisms are considered.

#### DGP 1: Coordinate-smooth truth

The structured coefficient fields are deterministic smooth functions of the lattice coordinates and do not depend on the adjacency graph.

Let $x_i$ and $y_i$ denote the centred and range-normalised column and row coordinates of node $i$, respectively, so that both coordinates lie in $[-1,1]$. The raw structured fields are

$$
h_{i1}=0.8x_i,
\qquad
h_{i2}=0.8y_i,
\qquad
h_{i3}=0.6x_iy_i.
$$

These fields correspond to the coordinate-based construction used in `generate_truth.R`. They generate gradual changes in amplitude, phase, and local peak height over the lattice.

For each component $k$, the centred and standardised version of $\boldsymbol{h}_k$ defines $\boldsymbol{\theta}^{(S)}_k$. The scale factor $s_k$ is set to the across-node standard deviation of the corresponding raw field, using divisor $N$, so that the case $\alpha=1$ reproduces the original coordinate-field scale.

The unstructured component is obtained by generating a standard-normal vector, removing its projection onto the intercept and $\boldsymbol{\theta}^{(S)}_k$, and rescaling the residual vector to unit across-node variance. This construction ensures the centring, scaling, and orthogonality conditions in Section 3.5.

The same structured fields and the same unstructured-component realisation are used for the rook and queen fits within a replication. Thus, neighbourhood density changes only the graph supplied to the estimator.

The setting $\alpha=1$ is reserved for the known-answer sanity check in which graph smoothing should clearly help. It is not part of the 16-cell core factorial design.

#### DGP 2: Block- or cluster-structured truth

The $5\times5$ lattice is partitioned into two connected rectangular clusters according to the fixed label matrix

$$
C=
\begin{pmatrix}
1&1&2&2&2\\
1&1&2&2&2\\
1&1&2&2&2\\
1&1&2&2&2\\
1&1&2&2&2
\end{pmatrix}.
$$

Rows and columns of $C$ correspond to the lattice row and column coordinates. Cluster 1 contains 10 nodes and cluster 2 contains 15 nodes. Both clusters are connected, and their common boundary lies inside the grid between columns 2 and 3.

For the three temporal coefficient fields, the raw cluster-level contrasts are

$$
\boldsymbol{\Gamma}=
\begin{pmatrix}
1&1&-1\\
-1&-1&1
\end{pmatrix},
$$

where row $c$ contains the three coefficient contrasts for cluster $c$. For node $i$ and temporal component $k$,

$$
h_{ik}=\Gamma_{c(i),k}.
$$

For each component, the node-level vector obtained by repeating these cluster contrasts is centred and scaled to unit across-node variance. This defines $\boldsymbol{\theta}^{(S)}_k$ in Section 3.5.

The unstructured component $\boldsymbol{\theta}^{(U)}_k$ is generated from a standard-normal vector, residualised with respect to the intercept and $\boldsymbol{\theta}^{(S)}_k$, and scaled to unit across-node variance. The final coefficient vector is

$$
\boldsymbol{\theta}_k =
s_k\left[
\sqrt{\alpha}\,\boldsymbol{\theta}^{(S)}_k +
\sqrt{1-\alpha}\,\boldsymbol{\theta}^{(U)}_k
\right].
$$

The same coefficient-scale factors $s_k$ are used for both core DGPs. Thus, $\alpha$ has the same variance-partition interpretation in the coordinate-smooth and cluster-structured scenarios.

The cluster labels, structured contrasts, and unstructured-component realisations are held fixed between the paired rook and queen fits.

Under edge-sharing grid adjacency, the complete set of nodes adjacent to the true cluster boundary consists of nodes in columns 2 and 3. For the primary boundary analysis, the evaluation sets are

$$
\mathcal{B}_{\mathrm{matched}} =
\lbrace i : \mathrm{col}(i)=3 \rbrace,
\qquad
\mathcal{I}_{\mathrm{matched}} =
\lbrace i : \mathrm{col}(i)=4 \rbrace.
$$

Both sets contain five nodes from cluster 2. Nodes are matched by row, so every boundary node and its interior comparator have the same outer-grid status and the same node degree under both rook and queen adjacency. This prevents the boundary comparison from being confounded by cluster membership or node degree.

The primary boundary effect is evaluated as a difference-in-differences relative to nodewise smoothing, as defined in Section 6.6. Results based on all nodes adjacent to the cluster boundary may additionally be reported as a secondary descriptive analysis.

### 3.7 Noise level and signal-to-noise ratios

The core simulation uses two noise levels:

$$
\sigma\in\{0.2,0.5\}.
$$

Within a replication, the same standard-normal error array is multiplied by the two values of $\sigma$. This provides a paired comparison of the noise-level scenarios.

The two values represent low and high observation noise. Since information under independent Gaussian noise is proportional to $1/\sigma^2$, these settings differ by a factor of

$$
\frac{0.5^2}{0.2^2}=6.25
$$

in information per node.

Two signal-to-noise ratios are documented. First, define the grand mean of the true evaluations as

$$
\overline{f} =
\frac{1}{NT}
\sum_{i=1}^{N}\sum_{j=1}^{T}f_i(t_j).
$$

The overall signal-to-noise ratio is

$$
\mathrm{SNR}_{\mathrm{overall}} =
\frac{
\frac{1}{NT}
\sum_{i=1}^{N}\sum_{j=1}^{T}
\left[f_i(t_j)-\overline{f}\right]^2
}{
\sigma^2
}.
$$

Because every method reconstructs the common mean function $\mu(t)$, the contrast between methods is driven primarily by the node-specific deviations. Therefore, the deviation signal-to-noise ratio is additionally reported:

$$
\mathrm{SNR}_{\mathrm{dev}} =
\frac{
\frac{1}{NT}
\sum_{i=1}^{N}\sum_{j=1}^{T}
\delta_i(t_j)^2
}{
\sigma^2
}.
$$

The pilot verifies that the two core noise levels produce meaningfully different values of both SNR measures and visibly different reconstruction difficulty. Any adjustment is made and documented before the main simulation results are inspected.

The core design contains one noisy curve per node. Under independent errors on a common observation grid, averaging $n_i$ replicate curves at node $i$ is equivalent to observing one curve with noise standard deviation $\sigma/\sqrt{n_i}$. Therefore, the number of replicate curves per node is not varied as a separate factor. This equivalence, and its dependence on the independent common-grid error assumption, is stated as a limitation.

An additional $\sigma=0.1$ high-information stage is optional and is not part of the 16-cell core design.

### 3.8 Core simulation scenarios

The core factorial design combines the following factors:

| Factor | Levels |
|---|---|
| Structure of the truth | Coordinate-smooth, cluster-structured |
| Structured-signal proportion | $\alpha=0.5,0.9$ |
| Neighbourhood density | Sparse, dense |
| Noise level | $\sigma=0.2,0.5$ |

This results in

$$
2\times2\times2\times2=16
$$

core scenarios.

#### Prespecified known-answer sanity check

Before the main simulation is launched, a separate known-answer sanity check is conducted in a setting where the graph penalty should clearly improve reconstruction. Its fixed design is:

| Component | Setting |
|---|---|
| Truth structure | Coordinate-smooth |
| Structured-signal proportion | $\alpha=1$ |
| Estimator graph | Rook |
| Noise level | $\sigma=0.5$ |
| Lattice | $5\times5$, $N=25$ |
| Observation grid | $T=50$ |
| Replications | $B_{\mathrm{sanity}}=50$ |
| Primary comparison | Network-REML versus nodewise smoothing |

Both methods are fitted to the same simulated dataset within each replication. Define the replication-specific sanity-check improvement as

$$
\Delta^{(\mathrm{sanity})}_b =
\mathrm{AISE}_{b,\mathrm{nodewise}} -
\mathrm{AISE}_{b,\mathrm{network}}.
$$

Its estimated mean and paired Monte Carlo standard error are

$$
\widehat{\Delta}^{(\mathrm{sanity})} =
\frac{1}{B_{\mathrm{sanity}}}
\sum_{b=1}^{B_{\mathrm{sanity}}}
\Delta^{(\mathrm{sanity})}_b,
$$

and

$$
\mathrm{MCSE}
\left(
\widehat{\Delta}^{(\mathrm{sanity})}
\right) =
\frac{
\mathrm{SD}
\left(
\Delta^{(\mathrm{sanity})}_1,\ldots,
\Delta^{(\mathrm{sanity})}_{B_{\mathrm{sanity}}}
\right)
}{
\sqrt{B_{\mathrm{sanity}}}
}.
$$

The sanity check is passed only if all three prespecified criteria hold:

$$
\widehat{\Delta}^{(\mathrm{sanity})}>0,
$$

$$
\widehat{\Delta}^{(\mathrm{sanity})}>
2\,
\mathrm{MCSE}
\left(
\widehat{\Delta}^{(\mathrm{sanity})}
\right),
$$

and

$$
\frac{1}{B_{\mathrm{sanity}}}
\sum_{b=1}^{B_{\mathrm{sanity}}}
\mathbf{1}
\left(
\Delta^{(\mathrm{sanity})}_b>0
\right)
\geq0.80.
$$

These criteria require a positive average improvement, an improvement that is clear relative to Monte Carlo uncertainty, and improvement in at least 80% of replications. The check is diagnostic and is not included in the 16-cell core design or in the final core performance estimates.

If the check fails, the main simulation is paused. Graph-to-node alignment, truth generation, basis representation, penalty construction, smoothing-parameter extraction, prediction, and performance-measure implementation are then inspected. The DGP and pass criteria are not tuned in response to the observed sanity-check performance. Any implementation correction is documented, and the sanity check is rerun with a new prespecified seed set before the main simulation starts.

This formal $T=50$ check is distinct from the existing $T=30$ software smoke test, whose purpose is only to verify that the workflow runs.

#### Additional negative controls

Two distinct negative controls are included because an uninformative graph and a false graph represent different failure modes.

**Uninformative graph control.** The setting $\alpha=0$ removes the structured coefficient component. It is evaluated for the coordinate-based DGP at $\sigma=0.5$ using both rook and queen estimator graphs. These two cells test whether graph regularisation is harmful when the supplied graph contains no systematic information about curve similarity.

**False graph control.** A fixed rewired graph is constructed from the rook graph using degree-preserving double-edge swaps. The rewired graph must:

- retain the original node labels;
- retain the rook degree sequence and number of edges;
- remain simple and connected;
- replace at least 75% of the original rook edges.

The rewiring seed and final adjacency list are fixed and saved before the main simulation. The same true curves, observation errors, and method settings are used for the correct-rook and rewired-graph fits.

The false-graph control is evaluated at $\alpha=0.9$ and $\sigma=0.5$ for both the coordinate-smooth and cluster-structured DGPs. The corresponding correct-rook cells are already part of the core factorial design, so only the additional rewired network fits are required.

To quantify whether the rewired graph is actively incompatible with the truth, define the purely structured deviation curve for DGP $d$ as

$$
\delta^{(S)}_{i,d}(t) =
\sum_{k=1}^{3}s_{k,d}\theta^{(S)}_{ik,d}\phi_k(t).
$$

For graph $G$, its integrated unnormalised Laplacian energy is

$$
E_{d,G} =
\frac{
\sum_{(i,j)\in\mathcal{E}_G}
\int_0^1
\left[
\delta^{(S)}_{i,d}(t)-\delta^{(S)}_{j,d}(t)
\right]^2dt
}{
\sum_{i=1}^{N}
\int_0^1
\left[\delta^{(S)}_{i,d}(t)\right]^2dt
}.
$$

The false-graph incompatibility ratio is

$$
Q_d =
\frac{E_{d,\mathrm{rewired}}}{E_{d,\mathrm{rook}}}.
$$

In addition to the structural requirements above, the retained rewired graph must satisfy

$$
Q_d\geq1.5
$$

separately for both the coordinate-smooth and cluster-structured DGPs. Rewiring candidates are generated using a prespecified seed sequence, and the first graph satisfying all structural and energy criteria is retained. The criterion is evaluated using only the fixed structured truth, before unstructured components and observation errors are generated and without inspecting any estimator-performance results. The final seed, adjacency list, edge-overlap proportion, $E_{d,G}$ values, and $Q_d$ values are saved before the main simulation.

The $\alpha=0$ control therefore represents useless graph information, whereas the rewired control represents actively misleading graph information.

### 3.9 Fixed larger-graph cells and optional robustness analyses

Node count is examined using a $10\times10$ lattice with

$$
N=100.
$$

Two larger-graph cells are prespecified:

| Purpose | Truth structure | $\alpha$ | Graph | $\sigma$ |
|---|---|---:|---|---:|
| Best case | Coordinate-smooth | 0.9 | Rook | 0.5 |
| Negative control | Coordinate-based unstructured truth | 0 | Rook | 0.5 |

The coordinate fields, temporal basis, coefficient scaling, and observation grid are constructed according to the same rules as for $N=25$. Only the lattice size and node count change.

Before the main larger-graph runs, five timing replications are performed for each of the two cells. Let $\widetilde{t}_{100}$ denote the median elapsed time of a network fit across these timing runs.

The planned replication count is

$$
B_{100}=
\begin{cases}
100, & 200\widetilde{t}_{100}\leq8\text{ hours},\\
50, & 200\widetilde{t}_{100}>8\text{ hours}.
\end{cases}
$$

The timing rule is evaluated before reconstruction-performance results from the $N=100$ cells are inspected. Convergence, log smoothing parameters, EDF, computation time, paired AISE improvement, and its MCSE are reported for both cells.

The two $N=100$ cells are prioritised ahead of the oracle diagnostics because node count is a distinct sample-size dimension.

Temporally correlated errors and additional graph perturbations are optional robustness analyses and are conducted only after completion of the core design, negative controls, and the two fixed $N=100$ cells.

### 3.10 Basis dimensions and projection check

The core basis dimensions are prespecified as

$$
K_{\mathrm{int}}=10, \qquad K_{\mathrm{dev}}=15.
$$

The larger deviation basis is used because the local peak in $\phi_3(t)$ cannot be represented sufficiently accurately by a P-spline basis with $K_{\mathrm{dev}}=10$.

During the pilot simulation, basis adequacy is assessed using the relative unpenalised $L^2$ projection error of the true node-specific deviations:

$$
\mathrm{RPE}_{\mathrm{dev}} =
\sqrt{
\frac{
\sum_{i=1}^{N}\int_0^1\left[\delta_i(t)-\Pi_{K_{\mathrm{dev}}}\delta_i(t)\right]^2dt
}{
\sum_{i=1}^{N}\int_0^1\delta_i(t)^2dt
}
}.
$$

Here, $\Pi_{K_{\mathrm{dev}}}\delta_i$ denotes the best unpenalised approximation of $\delta_i$ in the deviation basis used by the fitted model. The core simulation proceeds only if $\mathrm{RPE}_{\mathrm{dev}}\leq0.02$. Any necessary increase in the basis dimension will be made and documented before the main simulation results are inspected.

---

## 4. Estimands

### 4.1 Primary estimand

The primary estimand is the complete collection of true node-specific functions:

$$
\mathcal{F} = \{ f_1(t),\ldots,f_N(t) \}.
$$

Each method produces estimates

$$
\widehat{\mathcal{F}} = \{ \widehat f_1(t),\ldots,\widehat f_N(t) \}.
$$

### 4.2 Secondary estimands

Secondary estimands include:

- reconstruction accuracy at individual nodes;
- average reconstruction accuracy across all nodes;
- the effect of neighbourhood density, evaluated through the direct paired rook-versus-queen difference in AISE within otherwise identical design cells;
- the cost of an actively misleading graph, evaluated through the direct paired rewired-versus-correct-rook difference in AISE within otherwise identical false-graph control cells;
- reconstruction accuracy at cluster boundaries;
- reconstruction accuracy at interior cluster nodes;
- the true cluster structure in the block-structured DGP;
- the reconstruction error attainable under an optimal choice of smoothing parameters.

### 4.3 Distinction from spatial prediction

The primary simulation study concerns the reconstruction of noisy functions at observed graph nodes.

This differs from functional kriging, where the usual objective is to predict an entire curve at an unobserved spatial location. A comparison with functional kriging is therefore treated as a separate node-holdout analysis.

---

## 5. Methods

### Common fitting and evaluation rules

All reconstruction methods receive exactly the same observed $N\times T$ data matrix within a simulation cell and are evaluated against the same true curves on the prespecified 50-point grid $t_1,\ldots,t_T$. No method-specific preprocessing, observation weights, or outcome-dependent basis changes are used.

The temporal smooths use cubic P-spline bases with first-order difference penalties, specified by `bs = "ps"` and `m = c(2, 1)`. The basis dimensions are fixed before the main simulation at

$$
K_{\mathrm{int}}=10, \qquad K_{\mathrm{dev}}=15,
$$

subject only to the prespecified projection check in Section 3.10. Smoothing parameters are estimated by REML unless a method is explicitly labelled as an oracle. The basis dimension, penalty type, fitting method, and prediction grid are not adjusted after inspecting replication-specific reconstruction errors.

### M0: Raw observations

The no-smoothing benchmark is defined on the observation grid by

$$
\widehat f_i(t_j)=Y_i(t_j),
\qquad i=1,\ldots,N,\quad j=1,\ldots,T.
$$

It provides a descriptive reference for the reduction in observation noise achieved by smoothing. It is not a fitted model, does not use graph information, and is not treated as a primary competitor.

### M1: Network-weighted smoothing with REML

The primary method under investigation is `netf_smooth()` from the `netfunsmooth` package. The functional-regression representation and smoothing-parameter estimation follow Greven and Scheipl (2017) and Wood, Pya, and Säfken (2016).

The estimator uses a tensor-product penalised representation with:

- a temporal smoothness penalty;
- an MRF penalty based on the graph Laplacian;
- REML-based selection of the smoothing parameters.

The model is fitted to all $NT$ observations simultaneously with a Gaussian identity-link model. The functional intercept uses

```r
bs.int = list(bs = "ps", k = 10, m = c(2, 1))
```

and the temporal marginal basis of the node-by-time interaction uses

```r
bs.yindex = list(bs = "ps", k = 15, m = c(2, 1))
```

The graph marginal has `k = N`. The call uses `algorithm = "gam"`, `method = "REML"`, `tensortype = "ti"`, and `sandwich = "none"`; these settings are held fixed across all core cells. The sandwich setting concerns covariance estimation and does not change the reconstruction target, but fixing it avoids method-dependent computational overhead in the simulation.

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

The fitted values are returned for all observed nodes on the same 50-point grid used to generate and evaluate the data. Rook, queen, and rewired versions differ only in the neighbourhood list supplied to the MRF marginal.

### M2: Restricted network-weighted oracle smoothing

The oracle is a diagnostic reference and is not part of the full core benchmark. It uses the true curves and is therefore unavailable in real applications.

The fitted network model contains three smoothing parameters:

1. $\lambda_{\mathrm{int}}$ for the functional intercept;
2. $\lambda_{\mathrm{graph}}$ for the graph direction of the node-by-time interaction;
3. $\lambda_{\mathrm{time}}$ for the time direction of the node-by-time interaction.

For each oracle dataset, $\lambda_{\mathrm{int}}$ is fixed at its REML estimate. A two-dimensional grid search is then performed over the two interaction smoothing parameters.

The initial grid is defined relative to the corresponding REML estimates:

$$
\lambda_{d}^{(g)} =
\widehat{\lambda}_{d,\mathrm{REML}}\exp(g),
\qquad
g\in\{-4,-2,0,2,4\},
$$

for $d\in\{\mathrm{graph},\mathrm{time}\}$. This gives 25 candidate interaction-parameter pairs and includes the REML pair.

For each candidate pair, the model is fitted with all three smoothing parameters fixed, and the candidate with the smallest true node-averaged ISE is selected.

During the pilot, the frequency with which the selected candidate lies on an outer grid edge is recorded. If this occurs in more than 5% of oracle pilot fits, the grid is expanded before the main simulation. The final grid is then fixed and documented before oracle performance results are inspected.

The oracle is restricted to the following four cells:

| Truth structure | $\alpha$ | Graph | $\sigma$ |
|---|---:|---|---:|
| Coordinate-smooth | 0.9 | Rook | 0.2 |
| Coordinate-smooth | 0.9 | Rook | 0.5 |
| Cluster-structured | 0.9 | Rook | 0.5 |
| Coordinate-based negative control | 0 | Rook | 0.5 |

For each selected cell, the oracle is evaluated on the first

$$
B_{\mathrm{oracle}}=100
$$

completed main-simulation replications. The REML and oracle fits therefore use identical datasets.

The oracle analysis distinguishes limitations of the model representation from limitations of REML smoothing-parameter selection while keeping its computational cost separate from the core benchmark.

### M3: Nodewise univariate smoothing

Each node-specific function is smoothed separately over $t$:

$$
Y_i(t)=s_i(t)+\varepsilon_i(t).
$$

This method does not use graph information but preserves node-specific variation.

It is the most important comparator because it directly assesses whether the graph penalty provides an improvement beyond ordinary temporal smoothing.

For each node $i$, the $T=50$ observations are fitted separately using

```r
mgcv::gam(
  y ~ s(t, bs = "ps", k = 15, m = c(2, 1)),
  family = gaussian(),
  method = "REML"
)
```

Thus, M3 uses the same P-spline type, basis dimension, difference-penalty order, Gaussian likelihood, REML criterion, time range, and prediction grid as the temporal marginal of the network interaction. Each node receives its own independently estimated smoothing parameter. The basis dimension is fixed across nodes and replications and is not tuned using true-curve reconstruction error.

### M4: Pooled mean smoothing

Node identity is ignored, and a single common mean function is estimated:

$$
Y_i(t)=\mu(t)+\varepsilon_i(t).
$$

The same estimated function is then assigned to every node:

$$
\widehat f_i(t)=\widehat\mu(t).
$$

The model is fitted once to all $NT$ scalar observations using

```r
mgcv::gam(
  y ~ s(t, bs = "ps", k = 10, m = c(2, 1)),
  family = gaussian(),
  method = "REML"
)
```

It therefore uses the same temporal basis specification as the functional intercept of M1. The resulting mean curve is predicted on the prespecified 50-point grid and copied to all nodes.

This method represents complete pooling. It may have low variance but cannot represent systematic differences between node-specific curves. It is therefore a diagnostic benchmark for the cost of ignoring node identity rather than a primary competitor to M1.

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
\mathrm{ISE}_{ibm} = \int_0^1 \left[ \widehat f_{ibm}(t)-f_{ib}(t) \right]^2dt.
$$

On the discrete observation grid, the integral is approximated numerically, preferably using the trapezoidal rule.

### 6.2 Node-averaged ISE

Within each simulation replication, the ISE is averaged across nodes:

$$
\mathrm{AISE}_{bm} = \frac{1}{N} \sum_{i=1}^{N} \mathrm{ISE}_{ibm}.
$$

Across $B$ simulation replications, the estimated mean integrated squared error is

$$
\widehat{\mathrm{MISE}}_m = \frac{1}{B} \sum_{b=1}^{B} \mathrm{AISE}_{bm}.
$$

This is the primary absolute performance measure.

### 6.3 Primary paired comparison and relative measures

Because every method is fitted to the same simulated dataset within replication $b$, the primary comparison with nodewise smoothing is paired.

For method $m$, define the replication-specific AISE improvement

$$
\Delta_{bm} =
\mathrm{AISE}_{b,\mathrm{nodewise}} -
\mathrm{AISE}_{bm}.
$$

Thus, $\Delta_{bm}>0$ indicates that method $m$ has a lower node-averaged reconstruction error than nodewise smoothing.

The estimated mean paired improvement is

$$
\widehat{\Delta}_m =
\frac{1}{B}\sum_{b=1}^{B}\Delta_{bm}.
$$

Its Monte Carlo standard error is

$$
\mathrm{MCSE}\left(\widehat{\Delta}_m\right) =
\frac{
\mathrm{SD}\left(\Delta_{1m},\ldots,\Delta_{Bm}\right)
}{
\sqrt{B}
}.
$$

This paired MCSE is the primary Monte Carlo uncertainty measure for the network-versus-nodewise comparison.

At the aggregate level, the MISE ratio is

$$
\mathrm{rMISE}_m =
\frac{
\widehat{\mathrm{MISE}}_m
}{
\widehat{\mathrm{MISE}}_{\mathrm{nodewise}}
}.
$$

Values below 1 favour method $m$.

For descriptive replication-level summaries, relative improvement is defined as

$$
\mathrm{RI}_{bm} =
100\frac{\Delta_{bm}}{\mathrm{AISE}_{b,\mathrm{nodewise}}}.
$$

Because this measure contains a random denominator and may produce extreme ratios, it is summarised using the median, interquartile range, and boxplots. It is not used as the primary inferential comparison and is not averaged to determine the number of replications.

### 6.4 Direct paired contrast for neighbourhood density

To directly quantify the effect of neighbourhood density, rook and queen network fits are compared within the same replication. Because both fits use the same true curves, unstructured-component realisation, observation errors, and observed data, this is a paired comparison.

Let $c$ denote a design cell defined by truth structure, structured-signal proportion, noise level, and graph size, but excluding the graph-density factor. For replication $b$ in cell $c$, define

$$
\Delta^{(\mathrm{density})}_{bc} =
\mathrm{AISE}_{bc,\mathrm{rook}} -
\mathrm{AISE}_{bc,\mathrm{queen}}.
$$

The cell-specific mean density effect is estimated by

$$
\widehat{\Delta}^{(\mathrm{density})}_{c} =
\frac{1}{B_c}
\sum_{b=1}^{B_c}
\Delta^{(\mathrm{density})}_{bc}.
$$

Its paired Monte Carlo standard error is

$$
\mathrm{MCSE}\left(\widehat{\Delta}^{(\mathrm{density})}_{c}\right) =
\frac{
\mathrm{SD}\left(
\Delta^{(\mathrm{density})}_{bc}: b=1,\ldots,B_c
\right)
}{
\sqrt{B_c}
}.
$$

The contrast is interpreted as follows:

- $\widehat{\Delta}^{(\mathrm{density})}_{c}>0$ indicates lower node-averaged reconstruction error under the queen graph;
- $\widehat{\Delta}^{(\mathrm{density})}_{c}<0$ indicates lower node-averaged reconstruction error under the rook graph;
- values close to zero indicate little evidence of an accuracy difference attributable to neighbourhood density.

The density contrast is reported separately for each design cell and is not averaged across heterogeneous combinations of truth structure, structured-signal proportion, noise level, or graph size. Nodewise, pooled, and raw estimators do not depend on the neighbourhood graph. They are therefore fitted once per simulated dataset and their results are reused in both graph-density comparisons.

### 6.5 Direct paired contrast for graph misspecification

To quantify the cost of supplying an actively misleading graph, the rewired-graph and correct-rook network fits are compared within the same false-graph control replication. Both fits use the same true curves, unstructured-component realisation, observation errors, observed data, basis dimensions, and fitting settings.

Let $d\in\{\mathrm{coordinate},\mathrm{cluster}\}$ denote the DGP. For replication $b$, define

$$
\Delta^{(\mathrm{false})}_{bd} =
\mathrm{AISE}_{bd,\mathrm{rewired}} -
\mathrm{AISE}_{bd,\mathrm{rook}}.
$$

The DGP-specific mean cost of graph misspecification is estimated by

$$
\widehat{\Delta}^{(\mathrm{false})}_{d} =
\frac{1}{B_d}
\sum_{b=1}^{B_d}
\Delta^{(\mathrm{false})}_{bd}.
$$

Its paired Monte Carlo standard error is

$$
\mathrm{MCSE}\left(\widehat{\Delta}^{(\mathrm{false})}_{d}\right) =
\frac{
\mathrm{SD}\left(
\Delta^{(\mathrm{false})}_{bd}: b=1,\ldots,B_d
\right)
}{
\sqrt{B_d}
}.
$$

The contrast is interpreted as follows:

- $\widehat{\Delta}^{(\mathrm{false})}_{d}>0$ indicates that the rewired graph produces greater node-averaged reconstruction error than the correct rook graph;
- $\widehat{\Delta}^{(\mathrm{false})}_{d}<0$ indicates that the rewired graph produces lower error in that DGP;
- values close to zero indicate little evidence of an accuracy cost attributable to graph misspecification.

The contrast is reported separately for the coordinate-smooth and cluster-structured DGPs and is not averaged across them.

### 6.6 Matched boundary difference-in-differences

For the cluster-structured DGP, the primary boundary analysis compares the matched boundary and interior sets $\mathcal{B}_{\mathrm{matched}}$ and $\mathcal{I}_{\mathrm{matched}}$ defined in Section 3.6.

For node $i$, replication $b$, and method $m$, first define the ISE difference relative to nodewise smoothing:

$$
d_{ibm} =
\mathrm{ISE}_{ibm} -
\mathrm{ISE}_{ib,\mathrm{nodewise}}.
$$

Positive values of $d_{ibm}$ indicate that method $m$ has a larger reconstruction error than nodewise smoothing at node $i$.

The replication-specific matched boundary effect is

$$
\mathrm{DiD}_{bm} =
\frac{1}{|\mathcal{B}_{\mathrm{matched}}|}
\sum_{i\in\mathcal{B}_{\mathrm{matched}}}d_{ibm} -
\frac{1}{|\mathcal{I}_{\mathrm{matched}}|}
\sum_{i\in\mathcal{I}_{\mathrm{matched}}}d_{ibm}.
$$

In the core $5\times5$ cluster DGP, $|\mathcal{B}_{\mathrm{matched}}|=|\mathcal{I}_{\mathrm{matched}}|=5$, and boundary and interior nodes are matched by lattice row.

The interpretation is:

- $\mathrm{DiD}_{bm}>0$: method $m$ incurs additional reconstruction error near the true cluster boundary, relative to both matched interior nodes and nodewise smoothing;
- $\mathrm{DiD}_{bm}=0$: there is no additional boundary-specific loss relative to nodewise smoothing;
- $\mathrm{DiD}_{bm}<0$: method $m$ performs relatively better at the boundary than at the matched interior nodes.

Across replications, the estimated mean boundary effect is

$$
\widehat{\mathrm{DiD}}_m =
\frac{1}{B}\sum_{b=1}^{B}\mathrm{DiD}_{bm}.
$$

Its Monte Carlo standard error is calculated from the replication-specific paired effects:

$$
\mathrm{MCSE}\left(\widehat{\mathrm{DiD}}_m\right) =
\frac{
\mathrm{SD}\left(\mathrm{DiD}_{1m},\ldots,\mathrm{DiD}_{Bm}\right)
}{
\sqrt{B}
}.
$$

The primary boundary analysis concerns the network-weighted estimator. Boundary and interior ISE values may additionally be reported separately as descriptive summaries.

### 6.7 Oracle gap

For each oracle replication $b$, define the paired oracle gap

$$
G_b =
\mathrm{AISE}_{b,\mathrm{network,REML}} -
\mathrm{AISE}_{b,\mathrm{network,oracle}}.
$$

Because the oracle grid includes the REML interaction-parameter pair, $G_b$ should be non-negative apart from numerical fitting differences.

The estimated mean oracle gap is

$$
\widehat{G} =
\frac{1}{B_{\mathrm{oracle}}}
\sum_{b=1}^{B_{\mathrm{oracle}}}G_b,
$$

with paired Monte Carlo standard error

$$
\mathrm{MCSE}(\widehat{G}) =
\frac{
\mathrm{SD}(G_1,\ldots,G_{B_{\mathrm{oracle}}})
}{
\sqrt{B_{\mathrm{oracle}}}
}.
$$

For each interaction direction, the difference between oracle-selected and REML log smoothing parameters is additionally summarised:

$$
\log\widehat{\lambda}_{d,\mathrm{oracle}} -
\log\widehat{\lambda}_{d,\mathrm{REML}},
\qquad
d\in\{\mathrm{graph},\mathrm{time}\}.
$$

Oracle-grid boundary-selection frequencies are reported separately. Smoothing-parameter ratios are not used as a primary performance measure.

### 6.8 Monte Carlo uncertainty

Monte Carlo standard errors are reported together with all estimated performance summaries using a formula appropriate to the corresponding statistic.

For the estimated MISE of method $m$, the marginal MCSE is

$$
\mathrm{MCSE}\left(\widehat{\mathrm{MISE}}_m\right) =
\frac{
\mathrm{SD}\left(\mathrm{AISE}_{1m},\ldots,\mathrm{AISE}_{Bm}\right)
}{
\sqrt{B}
}.
$$

For the primary network-versus-nodewise comparison, uncertainty is instead calculated from the paired differences $\Delta_{bm}$, as defined in Section 6.3. The two marginal MCSEs are not combined because this would ignore the within-replication correlation between methods.

The final number of simulation replications is determined from a prespecified MCSE target using the pilot estimate of the standard deviation of the paired differences. The exact rule is defined in Section 8.1.

### 6.9 Smoothing diagnostics, numerical stability, and computational cost

For every successful smooth-model fit, term-level and penalty-level diagnostics are stored before the fitted model object is discarded.

For the network-weighted model, the following three log smoothing parameters are recorded:

1. the smoothing parameter of the functional intercept;
2. the graph-direction smoothing parameter of the node-by-time interaction;
3. the time-direction smoothing parameter of the node-by-time interaction.

If the smoothing parameters are denoted by $\lambda_{m\ell}$, the stored values are

$$
\log\lambda_{m\ell},
\qquad
\ell=1,2,3.
$$

Effective degrees of freedom are recorded for every smooth term. For term $q$, let $r_q$ denote its maximum attainable rank after identifiability constraints. Basis saturation is flagged when

$$
\frac{\mathrm{EDF}_q}{r_q}\geq0.95.
$$

This rank-based definition is used instead of comparing EDF directly with the nominal basis dimension because centring and identifiability constraints may reduce the maximum attainable EDF.

A REML smoothing parameter is flagged as numerically extreme if its logarithm is non-finite or lies outside the prespecified diagnostic interval

$$
-15\leq\log\lambda_{m\ell}\leq15.
$$

This flag is treated as a numerical diagnostic rather than as an automatic fit failure. For oracle fits, a separate boundary flag records whether the selected value of either interaction smoothing parameter is located on the minimum or maximum edge of the prespecified two-dimensional grid.

For nodewise and pooled smoothers, their available log smoothing parameters, EDF values, and EDF-to-rank ratios are recorded analogously.

The following replication-level outcomes are additionally stored for every method:

- elapsed computation time;
- warning messages and warning count;
- convergence status and available optimiser diagnostics;
- fit failure and error message;
- non-finite predictions;
- numerically extreme smoothing-parameter flags;
- EDF-near-rank flags;
- oracle-grid boundary flags, where applicable.

For replication $b$, design cell $c$, and method $m$, define the fit-success indicator

$$
S_{bcm}=
\begin{cases}
1, & \text{if the fit returns without an error, does not report non-convergence, and produces finite predictions at all nodes and time points},\\
0, & \text{otherwise}.
\end{cases}
$$

When a fitting routine exposes an optimiser convergence indicator, that indicator must confirm convergence. A warning alone, an extreme smoothing-parameter flag, EDF saturation, or an oracle-grid boundary selection is retained as a diagnostic and does not by itself set $S_{bcm}=0$.

The analysis reports mean, median, and selected quantiles of computation time, log smoothing parameters, and EDF values, together with the proportions of warnings, failures, non-convergence, extreme smoothing parameters, EDF saturation, and oracle boundary selections.

Failed fits are not silently removed. Performance measures are reported for successful fits, while successful and unsuccessful fit counts are presented separately. A failed fit counts as a completed attempted fit: its error and diagnostic record are saved, and its seed is neither replaced nor rerun merely to obtain a successful result.

For a paired comparison of methods $m_1$ and $m_2$, define the number of jointly successful replications as

$$
n^{\mathrm{pair}}_{c,m_1,m_2} =
\sum_{b=1}^{B}
S_{bc m_1}S_{bc m_2}.
$$

The paired mean and its MCSE are calculated only from these jointly successful replications, with $n^{\mathrm{pair}}_{c,m_1,m_2}$ replacing $B$ in both the mean and the MCSE denominator. The attempted replication count, method-specific failure rates, jointly successful count, and characteristics of excluded datasets are reported alongside every affected comparison.

Marginal method-specific performance estimates analogously use that method's number of successful fits, rather than the attempted count $B$, in their mean and MCSE formulas.

The prespecified fit-failure threshold is 5% per method and design cell. If any fitted core method exceeds this threshold in any core cell during the pilot, the main simulation is paused while the implementation and numerical diagnostics are investigated. In the main simulation, failed fits are never replaced. If a final failure rate exceeds 5%, the affected paired performance result is explicitly labelled as conditional on joint fitting success and interpreted together with the failure analysis.

---

## 7. Optional Downstream Task: Clustering

For the cluster-structured DGP, a secondary analysis examines whether the estimated functions preserve the two true clusters. This analysis is optional and is conducted only after completion of the primary reconstruction benchmark.

For method $m$ and replication $b$, the estimated common mean is

$$
\widehat{\mu}_{bm}(t_j) =
\frac{1}{N}\sum_{i=1}^{N}\widehat{f}_{ibm}(t_j),
$$

and the estimated node-specific deviations are

$$
\widehat{\delta}_{ibm}(t_j) =
\widehat{f}_{ibm}(t_j) -
\widehat{\mu}_{bm}(t_j).
$$

Each estimated deviation is represented on the common evaluation grid using trapezoidal integration weights $w_j$. The clustering feature vector is

$$
\boldsymbol{z}_{ibm} =
\left(
\sqrt{w_1}\widehat{\delta}_{ibm}(t_1),
\ldots,
\sqrt{w_T}\widehat{\delta}_{ibm}(t_T)
\right)^{\top}.
$$

Euclidean distances between these vectors approximate the $L^2$ distances between the estimated deviation curves.

The same clustering procedure is used for every smoothing method:

1. apply $k$-means clustering to $\boldsymbol{z}_{ibm}$;
2. fix the number of clusters at the true value of two;
3. use 50 random initialisations with a deterministic replication-specific seed;
4. retain the solution with the smallest within-cluster sum of squares;
5. compare the estimated labels with the fixed true cluster labels.

If a method produces fewer than two distinct feature vectors, as may occur under complete pooling, the result is classified as a degenerate clustering and assigned an ARI of 0.

Clustering agreement is assessed using the Adjusted Rand Index of Hubert and Arabie (1985). The ARI equals 1 for identical partitions, has expected value approximately 0 under random agreement under its reference model, and may be negative when agreement is worse than expected by chance. No fixed lower bound of $-1$ is asserted.

ARI distributions, median ARI, and the frequency of degenerate clustering results are reported descriptively. Clustering performance does not determine the main simulation replication count.


---

## 8. Monte Carlo Design

### 8.1 Number of replications

A pilot simulation with

$$
B_{\mathrm{pilot}}=50
$$

attempted replications per core cell is conducted before the main simulation. Here and below, a replication count denotes the number of prespecified datasets on which fitting is attempted, not the number of successful fits. Replications are not added or regenerated to compensate for fitting failures.

The pilot is used to verify:

- correct data and curve generation;
- successful completion of the prespecified known-answer sanity check in Section 3.8;
- the intended structured and unstructured signal variation;
- basis adequacy through the projection-error criterion in Section 3.10;
- numerical stability and computation time;
- meaningful separation of the noise levels;
- correct implementation of all performance measures;
- the variability of the paired network-versus-nodewise AISE differences.

Let $c$ index a core simulation cell. The replication-count rule is based only on the primary network-REML-versus-nodewise comparison in the 16 core cells; negative controls, oracle fits, larger-graph cells, clustering, and other optional analyses do not determine the core replication count. Let

$$
n^{\mathrm{pilot,pair}}_c =
\sum_{b=1}^{B_{\mathrm{pilot}}}
S_{bc,\mathrm{network}}
S_{bc,\mathrm{nodewise}}
$$

denote the number of jointly successful pilot pairs in core cell $c$. From those pairs, let

$$
\overline{\Delta}^{\mathrm{pilot}}_c =
\frac{1}{n^{\mathrm{pilot,pair}}_c}
\sum_{b:\,S_{bc,\mathrm{network}}S_{bc,\mathrm{nodewise}}=1}
\Delta_{bc}
$$

denote the pilot mean paired improvement, and let $s^{\mathrm{pilot}}_{\Delta,c}$ denote the corresponding sample standard deviation.

The MCSE tolerance for cell $c$ is

$$
\tau_c =
\max\left(
0.10\left|\overline{\Delta}^{\mathrm{pilot}}_c\right|,
0.01\widehat{\mathrm{MISE}}^{\mathrm{pilot}}_{\mathrm{nodewise},c}
\right).
$$

The first term targets a paired MCSE no larger than 10% of the estimated paired effect. The second term provides an absolute fallback when the paired effect is close to zero in a core cell.

The required number of main-simulation replications for cell $c$ is estimated as

$$
B_c =
\left\lceil
\left(
\frac{s^{\mathrm{pilot}}_{\Delta,c}}{\tau_c}
\right)^2
\right\rceil.
$$

A common replication count is used for all core cells:

$$
B =
50\left\lceil
\frac{
\max\left(200,\max_{c}B_c\right)
}{50}
\right\rceil.
$$

Thus, the main simulation attempts at least 200 replications per cell and rounds the required number upward to a multiple of 50. Using a common $B$ preserves a balanced set of attempted datasets across cells. The realised jointly successful counts used in paired analyses may be smaller and are reported separately.

If the calculated $B$ is not computationally feasible, it is not silently truncated. The unmet MCSE target is documented, and optional analyses are reduced according to Section 11 before any change to the core replication count is considered.

Pilot replications are not included in the final performance estimates. If the DGP or fitting procedure is changed after the pilot, the relevant pilot checks and the replication calculation are repeated before the main simulation is started.

### 8.2 Random numbers, processing order, and reproducibility

Each replication $b$ is assigned a unique master seed. Deterministic substreams derived from this master seed are used for truth generation, unstructured coefficient components, observation errors, and any method-specific random operations. Consequently, results do not depend on the order in which cells or methods are executed.

Common random numbers are used for paired comparisons. Within replication $b$:

- within each DGP and graph size, the same centred and standardised structured and unstructured coefficient components are reused across all values of $\alpha$, with only the prespecified weights $\sqrt{\alpha}$ and $\sqrt{1-\alpha}$ changed;
- the same resulting truth is reused across the paired rook and queen cells;
- the same standard-normal error realisation is reused across noise levels and multiplied by the corresponding value of $\sigma$;
- every method fitted within a cell receives exactly the same true functions and observed dataset;
- the network and nodewise methods are therefore compared on identical data;
- raw, nodewise, and pooled estimators, which do not depend on the estimator graph, are computed once per simulated dataset and their results are reused in the paired rook-versus-queen analysis.

The simulation is processed in replication-major order:

1. generate all random components required for replication $b$;
2. run every prespecified core cell for replication $b$;
3. fit all methods and calculate all performance measures for those cells;
4. save a replication-level result object containing either a success record or a failure record for every attempted fit, and mark replication $b$ as complete;
5. proceed to replication $b+1$.

A replication is included in the aggregated results only after all required core fits have been attempted and a success or failure record has been saved for each one. If execution stops before this complete result object is saved, the interrupted replication is recomputed from its original master seed when the run resumes. By contrast, a recorded model-fitting failure is a completed attempted fit and does not cause the dataset to be regenerated, the seed to be changed, or the fit to be repeatedly retried.

Results are checkpointed by replication rather than only by scenario. A completion manifest records the identifiers of replication objects that have been completely saved, irrespective of whether individual fits succeeded. This ensures that an interrupted run leaves the same number of attempted replications in every core cell and remains directly evaluable.

Pilot and main-simulation seeds are stored separately. Scenario definitions, random-number generation, data generation, model fitting, performance evaluation, and result aggregation are implemented as separate functions.

---

## 9. Analysis and Presentation

### 9.1 Main tables

For every core combination of truth structure, structured-signal proportion, neighbourhood density, and noise level, the following are reported:

- estimated MISE and its marginal MCSE for every available method, with raw and pooled results labelled as diagnostic benchmarks;
- mean paired AISE improvement $\widehat{\Delta}_m$ and its paired MCSE, with network-REML versus nodewise identified as the primary method comparison;
- direct paired rook-versus-queen AISE difference $\widehat{\Delta}^{(\mathrm{density})}_{c}$ and its paired MCSE for each otherwise identical design cell;
- direct paired rewired-versus-correct-rook AISE difference $\widehat{\Delta}^{(\mathrm{false})}_{d}$ and its paired MCSE for each false-graph control DGP;
- aggregate rMISE relative to nodewise smoothing;
- median and interquartile range of replication-specific RI;
- mean, median, and selected quantiles of computation time;
- numbers of successful, failed, warning-producing, and non-converged fits;
- numbers of attempted and jointly successful replications for every paired comparison, with conditional-on-joint-success labels where the 5% failure threshold is exceeded;
- log smoothing parameters, EDF, EDF-to-rank ratios, and diagnostic-flag rates.

For cluster-structured cells, matched boundary DiD and its paired MCSE are additionally reported. Restricted oracle results and the two $N=100$ cells are presented in separate tables.

The known-answer sanity check is reported separately with its mean paired improvement, paired MCSE, proportion of replications with positive improvement, and pass/fail status. It is not pooled with the core simulation results.

### 9.2 Main figures

The planned figures include:

1. an interaction plot of paired AISE improvement by structured-signal proportion and noise level;
2. a plot of the paired rook-versus-queen AISE difference and its MCSE by design cell;
3. a plot of the paired rewired-versus-correct-rook AISE difference and its MCSE for each false-graph control DGP;
4. a heatmap of aggregate rMISE for network-weighted versus nodewise smoothing;
5. boxplots of replication-specific paired AISE differences;
6. matched boundary difference-in-differences for the cluster DGP;
7. selected example curves showing truth, observations, and estimates;
8. term-level log smoothing parameters and EDF diagnostics;
9. optional clustering performance based on the Adjusted Rand Index.

Figures distinguish the 16 core cells from negative controls, larger-graph cells, oracle diagnostics, and other optional analyses.

---

## 10. Pilot-Dependent Checks Before the Main Simulation

The core design choices are prespecified above. The following checks are completed using pilot or timing runs before the main performance results are inspected:

1. run the prespecified known-answer sanity check in Section 3.8 and verify that all three pass criteria are satisfied;
2. verify that the relative unpenalised projection error is at most 2%;
3. verify that $\sigma=0.2$ and $\sigma=0.5$ produce meaningfully distinct overall and deviation SNR values;
4. verify that the fixed rewired graph satisfies all structural requirements and $Q_d\geq1.5$ separately for both DGPs;
5. determine the common core replication count using the paired-MCSE rule in Section 8.1;
6. verify that the fit-failure rate of every fitted core method is at most 5% in every core cell; if not, pause and investigate before launching the main simulation;
7. check the oracle-grid boundary-selection frequency and expand the grid if required;
8. determine whether $B_{100}=100$ or $B_{100}=50$ using the prespecified timing rule;
9. test with a toy example whether the `SpatFD` interface can predict a complete curve at an omitted node;
10. decide whether optional clustering, `SpatFD`, and robustness analyses remain feasible after completing the required analyses.

Every resulting choice, diagnostic value, seed, and final scenario table is saved before the main simulation is launched. The main simulation does not begin unless the known-answer check passes. If a pilot check changes the DGP or fitting procedure, the affected pilot checks, including the known-answer check when relevant, are repeated.

---

## 11. Prioritisation Under Computational Constraints

If computational time is limited, analyses are prioritised as follows:

1. **Known-answer sanity check:** coordinate-smooth truth with $\alpha=1$;
2. **Core benchmark:** network-REML versus nodewise smoothing;
3. **Negative controls:** $\alpha=0$ and the fixed rewired graph;
4. **Pooled baseline;**
5. **Matched cluster-boundary analysis and optional clustering task;**
6. **Two fixed $N=100$ cells;**
7. **Restricted oracle diagnostics;**
8. **Optional robustness analyses;**
9. **Leave-one-node-out comparison with `SpatFD`.**

If further reductions are necessary, `SpatFD` is omitted first, followed by optional robustness analyses and then a reduction of the oracle scope. The known-answer sanity check, core benchmark, negative controls, and fixed $N=100$ cells are retained.

---

## 12. Summary of the Core Simulation

| ADEMP component | Specification |
|---|---|
| Aim | Assess the benefit and potential harm of graph information for reconstructing node-specific curves |
| Data | $N=25$, $5\times5$ lattice, $T=50$, one noisy curve per node, Gaussian independent errors |
| Truth | Coordinate-smooth or cluster-structured curves generated independently of the estimator graph |
| Structured signal | $\alpha=0.5,0.9$ with fixed marginal coefficient variance |
| Graph density | Rook and queen adjacency applied to the same truth and observations |
| Noise level | $\sigma=0.2,0.5$, verified using overall and deviation SNR |
| Primary comparison methods | Network-REML and nodewise smoothing |
| Diagnostic benchmarks | Raw observations and pooled mean smoothing |
| Fixed fitting settings | Cubic P-splines with first-order difference penalties; $K_{\mathrm{int}}=10$, $K_{\mathrm{dev}}=15$; Gaussian identity-link models; REML; common 50-point prediction grid |
| Primary estimand | True node-specific functions $f_i(t)$ |
| Primary measures | MISE and paired AISE improvement over nodewise smoothing with paired MCSE |
| Secondary measures | Direct paired rook-versus-queen and rewired-versus-correct-rook AISE contrasts, rMISE, descriptive RI, matched boundary DiD, runtime, smoothing diagnostics, warnings and failures |
| Core design | 16 factorial cells with a pilot-determined common $B\geq200$ attempted replications; failed fits are recorded and never replaced |
| Failure policy | Success requires no error, no reported non-convergence, and finite predictions; pilot failure threshold 5% per fitted core method and cell; paired MCSE uses the jointly successful count |
| Known-answer check | Coordinate-smooth truth, $\alpha=1$, rook graph, $\sigma=0.5$, $N=25$, $T=50$, and $B_{\mathrm{sanity}}=50$, with three prespecified pass criteria |
| Negative controls | $\alpha=0$ uninformative graph and fixed degree-preserving rewired graph with $Q_d\geq1.5$ for both DGPs |
| Larger graph | Two fixed $N=100$ cells: best case and negative control |
| Oracle | Restricted 25-point grid in four cells with $B_{\mathrm{oracle}}=100$ |
| Clustering | Optional secondary analysis using two-cluster $k$-means and ARI |
| `SpatFD` | Optional separate leave-one-node-out spatial-prediction experiment |

## References

Chung, F. R. K. (1997). *Spectral Graph Theory*. CBMS Regional Conference Series in Mathematics, Vol. 92. American Mathematical Society. <https://bookstore.ams.org/cbms-92>

Greven, S., and Scheipl, F. (2017). A general framework for functional regression modelling. *Statistical Modelling*, 17(1–2), 1–35. <https://doi.org/10.1177/1471082X16681317>

Hubert, L., and Arabie, P. (1985). Comparing partitions. *Journal of Classification*, 2, 193–218. <https://doi.org/10.1007/BF01908075>

Morris, T. P., White, I. R., and Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38, 2074–2102. <https://doi.org/10.1002/sim.8086>

Wood, S. N., Pya, N., and Säfken, B. (2016). Smoothing parameter and model selection for general smooth models. *Journal of the American Statistical Association*, 111, 1548–1575. <https://doi.org/10.1080/01621459.2016.1180986>
