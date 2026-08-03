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

The formal core simulation uses $T=50$. The existing known-answer and smoke-test scripts use $T=30$ only as software-validation checks; the formal simulation drivers explicitly use the prespecified 50-point grid.

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

The primary boundary effect is evaluated as a difference-in-differences relative to nodewise smoothing, as defined in Section 6.4. Results based on all nodes adjacent to the cluster boundary may additionally be reported as a secondary descriptive analysis.

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

Before the main simulation, the pilot verifies that the rewired graph has substantially larger unnormalised Laplacian energy than the correct rook graph for the fixed structured truth. If this condition is not met, the rewiring is repeated using the prespecified seed sequence and the first graph satisfying the criteria above is retained. This selection is completed without inspecting method-performance results.

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

The primary method under investigation is `netf_smooth()` from the `netfunsmooth` package. The functional-regression representation and smoothing-parameter estimation follow Greven and Scheipl (2017) and Wood, Pya, and Säfken (2016).

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

### 6.4 Matched boundary difference-in-differences

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

### 6.5 Oracle gap

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

### 6.6 Monte Carlo uncertainty

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

### 6.7 Smoothing diagnostics, numerical stability, and computational cost

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

The analysis reports mean, median, and selected quantiles of computation time, log smoothing parameters, and EDF values, together with the proportions of warnings, failures, non-convergence, extreme smoothing parameters, EDF saturation, and oracle boundary selections.

Failed fits are not silently removed. Performance measures are reported for successful fits, while successful and unsuccessful fit counts are presented separately.

Paired comparisons use replications in which both methods being compared fitted successfully. If failure rates are non-negligible, the number and characteristics of excluded paired comparisons are reported.

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

replications per core cell is conducted before the main simulation.

The pilot is used to verify:

- correct data and curve generation;
- the intended structured and unstructured signal variation;
- basis adequacy through the projection-error criterion in Section 3.10;
- numerical stability and computation time;
- meaningful separation of the noise levels;
- correct implementation of all performance measures;
- the variability of the paired network-versus-nodewise AISE differences.

Let $c$ index a core simulation cell. From the pilot replications, let

$$
\overline{\Delta}^{\mathrm{pilot}}_c =
\frac{1}{B_{\mathrm{pilot}}}
\sum_{b=1}^{B_{\mathrm{pilot}}}
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

The first term targets a paired MCSE no larger than 10% of the estimated paired effect. The second term provides an absolute fallback when the paired effect is close to zero, as expected in some negative-control cells.

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

Thus, the main simulation uses at least 200 replications per cell and rounds the required number upward to a multiple of 50. Using a common $B$ preserves balanced paired comparisons across cells.

Pilot replications are not included in the final performance estimates. If the DGP or fitting procedure is changed after the pilot, the relevant pilot checks and the replication calculation are repeated before the main simulation is started.

### 8.2 Random numbers, processing order, and reproducibility

Each replication $b$ is assigned a unique master seed. Deterministic substreams derived from this master seed are used for truth generation, unstructured coefficient components, observation errors, and any method-specific random operations. Consequently, results do not depend on the order in which cells or methods are executed.

Common random numbers are used for paired comparisons. Within replication $b$:

- the same structured and unstructured coefficient components are reused across the paired rook and queen cells;
- the same standard-normal error realisation is reused across noise levels and multiplied by the corresponding value of $\sigma$;
- every method fitted within a cell receives exactly the same true functions and observed dataset;
- the network and nodewise methods are therefore compared on identical data.

The simulation is processed in replication-major order:

1. generate all random components required for replication $b$;
2. run every prespecified core cell for replication $b$;
3. fit all methods and calculate all performance measures for those cells;
4. save a replication-level result object and mark replication $b$ as complete;
5. proceed to replication $b+1$.

A replication is included in the aggregated results only after all required core cells have been completed and saved. If execution stops during replication $b$, that incomplete replication is recomputed when the run resumes.

Results are checkpointed by replication rather than only by scenario. A completion manifest records the successfully saved replication identifiers. This ensures that an interrupted run leaves the same number of completed replications in every core cell and remains directly evaluable.

Pilot and main-simulation seeds are stored separately. Scenario definitions, random-number generation, data generation, model fitting, performance evaluation, and result aggregation are implemented as separate functions.

---

## 9. Analysis and Presentation

### 9.1 Main tables

For every core combination of truth structure, structured-signal proportion, neighbourhood density, and noise level, the following are reported:

- estimated MISE and its marginal MCSE for every method;
- mean paired AISE improvement $\widehat{\Delta}_m$ and its paired MCSE;
- aggregate rMISE relative to nodewise smoothing;
- median and interquartile range of replication-specific RI;
- mean, median, and selected quantiles of computation time;
- numbers of successful, failed, warning-producing, and non-converged fits;
- log smoothing parameters, EDF, EDF-to-rank ratios, and diagnostic-flag rates.

For cluster-structured cells, matched boundary DiD and its paired MCSE are additionally reported. Restricted oracle results and the two $N=100$ cells are presented in separate tables.

### 9.2 Main figures

The planned figures include:

1. an interaction plot of paired AISE improvement by structured-signal proportion and noise level;
2. a heatmap of aggregate rMISE for network-weighted versus nodewise smoothing;
3. boxplots of replication-specific paired AISE differences;
4. matched boundary difference-in-differences for the cluster DGP;
5. selected example curves showing truth, observations, and estimates;
6. term-level log smoothing parameters and EDF diagnostics;
7. optional clustering performance based on the Adjusted Rand Index.

Figures distinguish the 16 core cells from negative controls, larger-graph cells, oracle diagnostics, and other optional analyses.

---

## 10. Pilot-Dependent Checks Before the Main Simulation

The core design choices are prespecified above. The following checks are completed using pilot or timing runs before the main performance results are inspected:

1. verify that the relative unpenalised projection error is at most 2%;
2. verify that $\sigma=0.2$ and $\sigma=0.5$ produce meaningfully distinct overall and deviation SNR values;
3. verify the structural properties and Laplacian-energy contrast of the fixed rewired graph;
4. determine the common core replication count using the paired-MCSE rule in Section 8.1;
5. check the oracle-grid boundary-selection frequency and expand the grid if required;
6. determine whether $B_{100}=100$ or $B_{100}=50$ using the prespecified timing rule;
7. test with a toy example whether the `SpatFD` interface can predict a complete curve at an omitted node;
8. decide whether optional clustering, `SpatFD`, and robustness analyses remain feasible after completing the required analyses.

Every resulting choice, diagnostic value, seed, and final scenario table is saved before the main simulation is launched. If a pilot check changes the DGP or fitting procedure, the affected pilot checks are repeated.

---

## 11. Prioritisation Under Computational Constraints

If computational time is limited, analyses are prioritised as follows:

1. **Core benchmark:** network-REML versus nodewise smoothing;
2. **Known-answer sanity check:** coordinate-smooth truth with $\alpha=1$;
3. **Negative controls:** $\alpha=0$ and the fixed rewired graph;
4. **Pooled baseline;**
5. **Matched cluster-boundary analysis and optional clustering task;**
6. **Two fixed $N=100$ cells;**
7. **Restricted oracle diagnostics;**
8. **Optional robustness analyses;**
9. **Leave-one-node-out comparison with `SpatFD`.**

If further reductions are necessary, `SpatFD` is omitted first, followed by optional robustness analyses and then a reduction of the oracle scope. The core benchmark, negative controls, and fixed $N=100$ cells are retained.

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
| Core methods | Network-REML, nodewise smoothing, and pooled smoothing |
| Primary estimand | True node-specific functions $f_i(t)$ |
| Primary measures | MISE and paired AISE improvement over nodewise smoothing with paired MCSE |
| Secondary measures | rMISE, descriptive RI, matched boundary DiD, runtime, smoothing diagnostics, warnings and failures |
| Core design | 16 factorial cells with a pilot-determined common $B\geq200$ |
| Negative controls | $\alpha=0$ uninformative graph and fixed degree-preserving rewired graph |
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
