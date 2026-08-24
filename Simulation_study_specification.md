# Simulation Study of Network-Weighted Smoothing for Functional Data

## Specification According to the ADEMP Framework

**Status:** Frozen prespecification before the corresponding simulation results
are generated or inspected. The document defines the core experiment, the two
fixed $N=100$ scalability cells, and the restricted-oracle diagnostic. It does
not report simulation results or thesis conclusions.

**Framework:** ADEMP framework proposed by Morris, White, and Crowther (2019)

### Planned scope

| Component | Planned design |
|---|---|
| Known-answer check | Coordinate-smooth sanity setting, $B_{\mathrm{sanity}}=50$, network-REML versus nodewise. |
| Core factorial experiment | 16 $N=25$ cells; common pilot-determined $B\geq200$. |
| Methods | Raw, pooled, nodewise, and network-REML in the core experiment. |
| Cluster-boundary difference-in-differences | Prespecified for the eight cluster-structured cells. |
| $\alpha=0$ diagnostic | Coordinate-smooth negative control at $\sigma=0.5$ for rook and queen graphs. |
| Rewired-graph control | Two selected false-graph ablation cells with a fixed degree-preserving rewired rook graph. |
| $N=100$ scalability cells | Two fixed $10\times10$ rook-lattice cells, following a timing-based replication rule. |
| Restricted-oracle comparison | Coarse-to-fine diagnostic in four selected $N=25$ cells. |
| Clustering task | Optional secondary analysis; not part of the current execution priority. |

All components below are design specifications. Results, diagnostics, and
substantive conclusions are reported only in the thesis.

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
The primary aim is to evaluate whether the complete joint network-weighted estimator improves reconstruction relative to separate nodewise smoothing
under the examined simulation settings. The comparison includes a common functional mean, jointly selected temporal smoothing parameters, and graph regularisation. It therefore does not isolate the contribution of the
particular supplied graph edges.

The central research question is:

> Under which combinations of signal structure, structured-signal proportion,
> supplied graph, and observation noise does the complete joint estimator
> achieve lower reconstruction error than separate nodewise smoothing?

### 2.2 Secondary aims

The simulation additionally addresses the following questions:

1. How does the relative performance of the complete joint estimator vary with
   structured-signal proportion?
2. How does relative performance vary with observation noise?
3. How do rook and queen neighbourhood graphs compare when the simulated truth
   and observations are held fixed?
4. Does the estimator improve reconstruction for both coordinate-smooth and
   cluster-structured signals?
5. Is the relative advantage of the estimator attenuated near cluster
   boundaries?
6. How frequently do numerical warnings, convergence problems, or failed fits
   occur?
7. In selected high-structure, high-noise settings, does replacing the
   scientifically aligned rook edges by a fixed degree-preserving rewired graph
   increase reconstruction error?
8. How do runtime, numerical stability, and descriptive reconstruction
   accuracy change in the two fixed $N=100$ cells?
9. In four representative $N=25$ cells, how much of the REML reconstruction
   error can be reduced by restricted reselection of the two interaction
   smoothing parameters?

The rook-versus-queen comparison changes the graph supplied to the estimator
while holding the true curves and observations fixed. It is descriptive and
does not by itself identify a graph-topology-specific causal effect.

The fixed rewired-graph control assesses the effect of one actively misleading
graph in two selected settings. It does not estimate a general effect over a
distribution of rewired graphs. The $N=100$ and restricted-oracle blocks are
reported separately from the primary core benchmark.

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

Gaussian and temporally independent errors provide a deliberately simple reference setting. Extensions to other error structures are outside the formal simulation and are listed in Section 13.

### 3.4 True node-specific curves

The true function at node $i$ is decomposed as

$$
f_i(t)=\mu(t)+\delta_i(t),
$$

where the common smooth mean function is fixed as

$$
\mu(t)=1+0.5\sin(2\pi t)+0.25\cos(4\pi t),
$$

and $\delta_i(t)$ is a node-specific deviation. This same $\mu(t)$ is used in every DGP, graph size, and replication.

For identifiability, the deviations are centred across nodes for every $t$:

$$
\sum_{i=1}^{N}\delta_i(t)=0.
$$

The node-specific deviations are constructed using multiple temporal basis functions:

$$
\delta_i(t) = \sum_{k=1}^{K}\theta_{ik}\phi_k(t).
$$

The core simulation fixes $K=3$ and uses exactly

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

Both methods are fitted to the same simulated dataset within each replication. Let $S_{b,\mathrm{network}}^{(\mathrm{sanity})}$ and $S_{b,\mathrm{nodewise}}^{(\mathrm{sanity})}$ be their success indicators, and define

$$
\mathcal{S}^{(\mathrm{sanity})}=
\left\lbrace b:S_{b,\mathrm{network}}^{(\mathrm{sanity})}
S_{b,\mathrm{nodewise}}^{(\mathrm{sanity})}=1\right\rbrace,
\qquad
n^{(\mathrm{sanity,pair})}=|\mathcal{S}^{(\mathrm{sanity})}|.
$$

For $b\in\mathcal{S}^{(\mathrm{sanity})}$, define the replication-specific sanity-check improvement as

$$
\Delta^{(\mathrm{sanity})}_b =
\mathrm{AISE}_{b,\mathrm{nodewise}} -
\mathrm{AISE}_{b,\mathrm{network}}.
$$

Its estimated mean and paired Monte Carlo standard error are

$$
\widehat{\Delta}^{(\mathrm{sanity})} =
\frac{1}{n^{(\mathrm{sanity,pair})}}
\sum_{b\in\mathcal{S}^{(\mathrm{sanity})}}
\Delta^{(\mathrm{sanity})}_b,
$$

and

$$
\mathrm{MCSE}
\left(
\widehat{\Delta}^{(\mathrm{sanity})}
\right) =
\frac{
\mathrm{SD}\left\lbrace
\Delta^{(\mathrm{sanity})}_b:b\in\mathcal{S}^{(\mathrm{sanity})}
\right\rbrace
}{
\sqrt{n^{(\mathrm{sanity,pair})}}
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
\frac{1}{n^{(\mathrm{sanity,pair})}}
\sum_{b\in\mathcal{S}^{(\mathrm{sanity})}}
\mathbf{1}
\left(
\Delta^{(\mathrm{sanity})}_b>0
\right)
\geq0.80.
$$

These criteria require a positive average improvement, an improvement that is clear relative to Monte Carlo uncertainty, and improvement in at least 80% of jointly successful replications. The attempted count, both method-specific success counts, and the jointly successful count are reported. If either method fails in more than 5% of the 50 attempted replications, the check is not passed even if the three performance criteria hold. The check is diagnostic and is not included in the 16-cell core design or in the main performance estimates.

If the check fails, the main simulation is paused. Graph-to-node alignment, truth generation, basis representation, penalty construction, smoothing-parameter extraction, prediction, and performance-measure implementation are then inspected. The DGP and pass criteria are not tuned in response to the observed sanity-check performance. Any implementation correction is documented, and the sanity check is rerun with a new prespecified seed set before the main simulation starts.

This formal $T=50$ check is distinct from the existing $T=30$ software smoke test, whose purpose is only to verify that the workflow runs.

The $\alpha=0$ diagnostic uses the coordinate-smooth DGP at $\sigma=0.5$ with
both rook and queen estimator graphs. The deliberately
structured coefficient component receives zero weight, so the node-specific
deviations are not deliberately aligned with the supplied graph.

This diagnostic does not isolate a graph-topology-specific effect. Instead, it
assesses whether the comparison with separate nodewise smoothing can reflect
generic joint shrinkage in addition to graph-aligned borrowing.

#### Rewired-graph control

To complement the $\alpha=0$ diagnostic, a fixed false-graph ablation was
prespecified in two selected settings: coordinate-smooth and cluster-structured
truth, both with $\alpha=0.9$ and $\sigma=0.5$. In each replication, the
correct rook-graph and rewired-graph network fits received exactly the same
true curves, observation errors, observed data, basis dimensions, and fitting
settings. Only the graph supplied to the additional network fit changes.

The fixed rewired graph is generated from the $5\times5$ rook graph by
degree-preserving double-edge swaps before estimator performance is inspected.
It must remain connected and retain the original degree sequence. Its seed,
adjacency-list hash, edge-overlap proportion, and swap audit are stored with
the simulation output.

The control uses the same attempted replication identifiers and RNG streams as
the corresponding correct-rook core cells. It is interpreted as a targeted
graph-ablation design, not as an average effect over all possible misspecified
graphs.

### 3.9 Fixed larger-graph scalability cells

Node count is a distinct design dimension from the amount of information in a
curve. The extension therefore uses a regular $10\times10$ lattice with

$$
N=100.
$$

It contains exactly the following two fixed rook-graph cells:

| Purpose | Truth structure | $\alpha$ | $\sigma$ |
|---|---|---:|---:|
| High-structure scalability cell | Coordinate-smooth | 0.9 | 0.5 |
| Negative control | Coordinate-based unstructured truth | 0 | 0.5 |

Both cells use one noisy curve per node and the same common $T=50$ grid,
temporal basis functions, mean curve, Gaussian iid observation model, and
basis dimensions as the core experiment. Under the iid common-grid
assumption, $n_i$ replicate curves at a node are equivalent to one curve with
noise standard deviation $\sigma/\sqrt{n_i}$. Consequently, replicate count is
not varied as an additional factor: the existing $\sigma\in\{0.2,0.5\}$ core
levels already span a factor of $6.25$ in information per node. This
equivalence does not generally hold for correlated errors or irregular grids.

The coordinate fields are evaluated on row and column coordinates
range-normalised to $[-1,1]$. To avoid confounding node count with coefficient
magnitude, the $N=100$ cells use the fixed $N=25$ reference scales

$$
\left(s_1^{(25)},s_2^{(25)},s_3^{(25)}\right)=
\left(\frac{0.8}{\sqrt{2}},\frac{0.8}{\sqrt{2}},0.3\right).
$$

The structured and unstructured components are nevertheless centred,
orthogonalised, and standardised on the 100-node lattice before mixing. The
network marginal uses $k=N=100$; all other fitting settings remain fixed. Only
network-REML and nodewise smoothing are fitted in this block; raw and pooled
benchmarks are not fitted or reported.

Before any $N=100$ reconstruction result is inspected, five independent timing
replications are run for each cell. Let $\widetilde t_{100,1}$ and
$\widetilde t_{100,2}$ be the resulting median complete-replication times. The
attempted count per cell is fixed by

$$
B_{100}=\max\left\{b\in\{100,50,25\}:
b\left(\widetilde t_{100,1}+\widetilde t_{100,2}\right)\leq8\text{ hours}\right\}.
$$

If even $25$ attempts per cell exceed this budget, the block is recorded as
computationally infeasible and is paused rather than reduced ad hoc. Failed
fits are retained and never replaced. Results will report method-specific and
joint success counts, MISE and MCSE, paired network-versus-nodewise AISE
improvement and MCSE, runtime, warnings, and failures. The $N=25$ values are
only descriptive references: no paired or causal-style effect of node count is
estimated.


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

Here, $\Pi_{K_{\mathrm{dev}}}\delta_i$ denotes the best unpenalised approximation of $\delta_i$ in the deviation basis used by the fitted model. The target is a relative projection error below 2%, preferably close to or below 1%. The core simulation proceeds only if $\mathrm{RPE}_{\mathrm{dev}}<0.02$. Any necessary increase in the basis dimension will be made and documented before the main simulation results are inspected.

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
- reconstruction accuracy at cluster boundaries;
- reconstruction accuracy at matched interior cluster nodes; and
- numerical stability and computation time for the $N=25$ methods;
- descriptive reconstruction accuracy and computation cost in the two fixed
  $N=100$ cells; and
- the restricted oracle gap between network-REML and the selected
  fixed-parameter oracle candidate in the four oracle cells.

For the two false-graph control cells, a supplementary estimand is
the paired difference in node-averaged integrated squared error (AISE) between
the network estimator supplied with the fixed rewired graph and the same
estimator supplied with the correct rook graph.

The $N=100$ and restricted-oracle estimands are reported separately from the
primary core estimands.

### 4.3 Distinction from spatial prediction

The primary simulation study concerns the reconstruction of noisy functions at observed graph nodes.

This differs from functional kriging, where the usual objective is to predict an entire curve at an unobserved spatial location. Functional kriging, including a possible comparison with `SpatFD`, is therefore outside the formal simulation design and is mentioned only as future work in Section 13. A fair comparison would require a separate node-holdout spatial-prediction ADEMP design rather than the observed-node reconstruction estimand used here.

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

The fitted values are returned for all observed nodes on the same 50-point grid used to generate and evaluate the data. The core experiment varies only the rook and queen neighbourhood lists supplied to the MRF marginal.

### M2: Restricted network-weighted oracle smoothing

The oracle is a diagnostic reference, not a real-data method. It uses the true
curves and is therefore unavailable in applications. Its results are reported
separately from the core benchmark.

The fitted network model contains three smoothing parameters:

1. $\lambda_{\mathrm{int}}$ for the functional intercept;
2. $\lambda_{\mathrm{graph}}$ for the graph direction of the node-by-time interaction;
3. $\lambda_{\mathrm{time}}$ for the time direction of the node-by-time interaction.

For each oracle dataset, $\lambda_{\mathrm{int}}$ is fixed at its REML estimate. A prespecified two-stage search is then performed over the two interaction smoothing parameters. This is a restricted coarse-to-fine oracle rather than a global oracle.

The initial grid is defined relative to the corresponding REML estimates:

$$
\lambda_{d}^{(g)} =
\widehat{\lambda}_{d,\mathrm{REML}}\exp(g),
\qquad
g\in\{-4,-2,0,2,4\},
$$

for $d\in\{\mathrm{graph},\mathrm{time}\}$. This gives 25 coarse candidate interaction-parameter pairs and includes the REML pair.

For each coarse candidate pair, the model is fitted with all three smoothing parameters fixed. After selecting the coarse candidate with the smallest true node-averaged ISE, a local $3\times3$ refinement is evaluated by adding

$$
(h_{\mathrm{graph}},h_{\mathrm{time}})
\in\{-1,0,1\}^2
$$

to the two log-scale offsets of the coarse winner. Duplicate candidates already evaluated in the coarse stage are not refitted. With the initial coarse grid, at most 33 distinct parameter pairs are therefore evaluated per dataset. The final oracle candidate is the candidate with the smallest true node-averaged ISE among all successfully evaluated coarse and refinement candidates.

During the pilot, the oracle search is attempted on the first 50 prespecified pilot replication identifiers in each of the four oracle cells. These 50 attempted searches are not replaced after failure and do not enter the final oracle summaries. Lower- and upper-edge selection frequencies of the **coarse-stage winner** are calculated among complete successful pilot searches and recorded separately for the graph and time directions, together with the attempted and successful counts. If a boundary is selected in more than 5% of successful oracle pilot searches, the coarse grid is expanded only in that direction and on that side by one additional offset spaced two log units from the current boundary. The same 50 pilot identifiers are rerun on the expanded grid and the check is repeated after each expansion. No coarse-grid boundary may be extended beyond the prespecified offset range $[-8,8]$. If the corresponding boundary-selection frequency remains above 5% at this maximum range, the search is labelled `range-inadequate`; the range is not enlarged further in response to oracle performance results.

The final coarse grid and all boundary-check results are fixed and documented before the main oracle analysis. Local-refinement offsets are clipped to the final $[-8,8]$ search range. A main-simulation selection on the minimum or maximum permitted offset is retained as a boundary diagnostic and does not trigger further expansion.

If candidate AISE values are equal within the numerical tolerance

$$
10^{-12}\max\{1,\min(\mathrm{AISE})\},
$$

ties are broken deterministically by first choosing the candidate with the smallest Euclidean distance from the REML log-parameter pair and then, if necessary, by lexicographic order of the graph and time offsets.

Candidate-level errors, non-convergence, and non-finite predictions are recorded. An oracle search is classified as successful only if the original REML fit and every candidate required by its prespecified coarse-to-fine path are successful. An incomplete search is retained as an attempted oracle fit but is excluded from the paired oracle-gap estimate; its replication identifier is not replaced.

The oracle is restricted to the following four cells:

| Truth structure | $\alpha$ | Graph | $\sigma$ |
|---|---:|---|---:|
| Coordinate-smooth | 0.9 | Rook | 0.2 |
| Coordinate-smooth | 0.9 | Rook | 0.5 |
| Cluster-structured | 0.9 | Rook | 0.5 |
| Coordinate-based negative control | 0 | Rook | 0.5 |

For each selected cell, the oracle is attempted for the first

$$
B_{\mathrm{oracle}}=100
$$

prespecified main-simulation replication identifiers, $b=1,\ldots,100$, irrespective of fit success. Thus $B_{\mathrm{oracle}}$ denotes attempted oracle replications. REML and oracle fits use identical datasets, failed searches are not replaced by later replication identifiers, and the selection of replication identifiers cannot depend on fitting outcomes.

Conditional on the fixed basis representation and the REML-selected functional-intercept smoothing parameter, the oracle analysis measures how much reconstruction error can be reduced by restricted reselection of the two interaction smoothing parameters. It does not separate all limitations of the model representation from smoothing-parameter selection, and its computational cost is kept separate from the core benchmark.

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

## 6. Performance Measures

### 6.1 Node-specific integrated squared error

For node $i$, method $m$, simulation replication $b$, and cell $c$, the integrated squared error is

$$
\mathrm{ISE}_{ibcm} = \int_0^1 \left[ \widehat f_{ibcm}(t)-f_{ibc}(t) \right]^2dt.
$$

On the equally spaced observation and evaluation grid, the integral is approximated using the fixed trapezoidal weights

$$
w_1=w_T=\frac{1}{2(T-1)},
\qquad
w_j=\frac{1}{T-1},\quad j=2,\ldots,T-1,
$$

so that

$$
\mathrm{ISE}_{ibcm}=
\sum_{j=1}^{T}w_j
\left[\widehat f_{ibcm}(t_j)-f_{ibc}(t_j)\right]^2.
$$

The same weights are used in every numerical integral and in the clustering features in Section 7.

### 6.2 Node-averaged ISE

Within each simulation replication, the ISE is averaged across nodes:

$$
\mathrm{AISE}_{bcm} = \frac{1}{N} \sum_{i=1}^{N} \mathrm{ISE}_{ibcm}.
$$

Let

$$
\mathcal{S}_{cm}=\{b:S_{bcm}=1\},
\qquad
n_{cm}=|\mathcal{S}_{cm}|
$$

denote the successful replication set and count for method $m$ in cell $c$. The estimated mean integrated squared error is

$$
\widehat{\mathrm{MISE}}_{cm} =
\frac{1}{n_{cm}}
\sum_{b\in\mathcal{S}_{cm}}\mathrm{AISE}_{bcm}.
$$

This is the primary absolute performance measure.

Individual-node reconstruction is summarised within each successful replication by

$$
\mathrm{MedISE}_{bcm}=\mathrm{median}_{i=1,\ldots,N}
\mathrm{ISE}_{ibcm},
\qquad
\mathrm{WorstISE}_{bcm}=\mathrm{max}_{i=1,\ldots,N}\mathrm{ISE}_{ibcm}.
$$

For each method and cell, the Monte Carlo mean and median of both quantities are reported over $\mathcal{S}_{cm}$. These are descriptive secondary measures and do not determine $B$.

### 6.3 Primary paired comparison and relative measures

Because every method is fitted to the same simulated dataset within replication $b$, the primary comparison with nodewise smoothing is paired.

For method $m$ in cell $c$, define the replication-specific AISE improvement

$$
\Delta_{bcm} =
\mathrm{AISE}_{bc,\mathrm{nodewise}} -
\mathrm{AISE}_{bcm}.
$$

Thus, $\Delta_{bcm}>0$ indicates that method $m$ has a lower node-averaged reconstruction error than nodewise smoothing.

For cell $c$, define the jointly successful set and count

$$
\mathcal{S}^{(\mathrm{pair})}_{c,m}=
\left\lbrace b:S_{bc,m}S_{bc,\mathrm{nodewise}}=1\right\rbrace,
\qquad
n^{(\mathrm{pair})}_{c,m}=|\mathcal{S}^{(\mathrm{pair})}_{c,m}|.
$$

The estimated mean paired improvement is

$$
\widehat{\Delta}_{c,m} =
\frac{1}{n^{(\mathrm{pair})}_{c,m}}
\sum_{b\in\mathcal{S}^{(\mathrm{pair})}_{c,m}}\Delta_{bcm}.
$$

Its Monte Carlo standard error is

$$
\mathrm{MCSE}\left(\widehat{\Delta}_{c,m}\right) =
\frac{
\mathrm{SD}\left\lbrace \Delta_{bcm}:b\in\mathcal{S}^{(\mathrm{pair})}_{c,m}\right\rbrace
}{
\sqrt{n^{(\mathrm{pair})}_{c,m}}
}.
$$

This paired MCSE is the primary Monte Carlo uncertainty measure for the network-versus-nodewise comparison.

For a failure-consistent aggregate relative measure, both numerator and denominator are calculated on the same jointly successful set:

$$
\mathrm{rMISE}_{c,m} =
\frac{
\displaystyle
\sum_{b\in\mathcal{S}^{(\mathrm{pair})}_{c,m}}
\mathrm{AISE}_{bc,m}
}{
\displaystyle
\sum_{b\in\mathcal{S}^{(\mathrm{pair})}_{c,m}}
\mathrm{AISE}_{bc,\mathrm{nodewise}}
}.
$$

Values below 1 favour method $m$.

For descriptive replication-level summaries, relative improvement is defined as

$$
\mathrm{RI}_{bcm} =
100\frac{\Delta_{bcm}}{\mathrm{AISE}_{bc,\mathrm{nodewise}}}.
$$

Because this measure contains a random denominator and may produce extreme ratios, it is calculated only for $b\in\mathcal{S}^{(\mathrm{pair})}_{c,m}$ and summarised using the median, interquartile range, and boxplots. It is not used as the primary inferential comparison and is not averaged to determine the number of replications.

### 6.4 Direct paired contrast for neighbourhood density

To directly quantify the effect of neighbourhood density, rook and queen network fits are compared within the same replication. Because both fits use the same true curves, unstructured-component realisation, observation errors, and observed data, this is a paired comparison.

Let $c$ denote a design cell defined by truth structure, structured-signal proportion, noise level, and graph size, but excluding the graph-density factor. For replication $b$ in cell $c$, define

$$
\Delta^{(\mathrm{density})}_{bc} =
\mathrm{AISE}_{bc,\mathrm{rook}} -
\mathrm{AISE}_{bc,\mathrm{queen}}.
$$

Let

$$
\mathcal{S}^{(\mathrm{density})}_c=
\left\lbrace b:S_{bc,\mathrm{rook}}S_{bc,\mathrm{queen}}=1\right\rbrace,
\qquad
n^{(\mathrm{density,pair})}_c=|\mathcal{S}^{(\mathrm{density})}_c|.
$$

The cell-specific mean density effect is estimated by

$$
\widehat{\Delta}^{(\mathrm{density})}_{c} =
\frac{1}{n^{(\mathrm{density,pair})}_c}
\sum_{b\in\mathcal{S}^{(\mathrm{density})}_c}
\Delta^{(\mathrm{density})}_{bc}.
$$

Its paired Monte Carlo standard error is

$$
\mathrm{MCSE}\left(\widehat{\Delta}^{(\mathrm{density})}_{c}\right) =
\frac{
\mathrm{SD}\left\lbrace
\Delta^{(\mathrm{density})}_{bc}:b\in\mathcal{S}^{(\mathrm{density})}_c
\right\rbrace
}{
\sqrt{n^{(\mathrm{density,pair})}_c}
}.
$$

The contrast is interpreted as follows:

- $\widehat{\Delta}^{(\mathrm{density})}_{c}>0$ indicates lower node-averaged reconstruction error under the queen graph;
- $\widehat{\Delta}^{(\mathrm{density})}_{c}<0$ indicates lower node-averaged reconstruction error under the rook graph;
- values close to zero indicate little evidence of an accuracy difference attributable to neighbourhood density.

The density contrast is reported separately for each design cell and is not averaged across heterogeneous combinations of truth structure, structured-signal proportion, noise level, or graph size. Nodewise, pooled, and raw estimators do not depend on the neighbourhood graph. They are therefore fitted once per simulated dataset and their results are reused in both graph-density comparisons.

### 6.5 Paired rewired-versus-rook contrast

For each of the two false-graph control cells, define

$$
\Delta^{(\mathrm{rewired})}_{bd} =
\mathrm{AISE}_{bd,\mathrm{rewired}} -
\mathrm{AISE}_{bd,\mathrm{rook}},
$$

where $d$ indexes the selected truth structure. Positive values indicate that
the fixed rewired graph yields higher reconstruction error than the correct
rook graph. The mean paired difference and its paired MCSE are computed over
the jointly successful pairs. The failure-consistent relative measure is

$$
\mathrm{rMISE}^{(\mathrm{rewired/rook})}_{d} =
\frac{\sum_b \mathrm{AISE}_{bd,\mathrm{rewired}}}
{\sum_b \mathrm{AISE}_{bd,\mathrm{rook}}},
$$

using the same paired-success set in both sums. Values above one favour the
correct rook graph. The replication-level win rate is the proportion of pairs
with $\Delta^{(\mathrm{rewired})}_{bd}>0$.

This targeted control assesses whether the supplied graph edges contribute
beyond generic joint shrinkage in the two selected settings. It does not
support a universal claim about all graph perturbations because only one fixed
rewired graph and two selected cells are evaluated.


### 6.6 Matched boundary difference-in-differences

For the cluster-structured DGP, the primary boundary analysis compares the matched boundary and interior sets $\mathcal{B}_{\mathrm{matched}}$ and $\mathcal{I}_{\mathrm{matched}}$ defined in Section 3.6.

For node $i$, replication $b$, cluster cell $c$, and method $m$, first define the ISE difference relative to nodewise smoothing:

$$
d_{ibcm} =
\mathrm{ISE}_{ibcm} -
\mathrm{ISE}_{ibc,\mathrm{nodewise}}.
$$

Positive values of $d_{ibcm}$ indicate that method $m$ has a larger reconstruction error than nodewise smoothing at node $i$.

The replication-specific matched boundary effect is

$$
\mathrm{DiD}_{bcm} =
\frac{1}{|\mathcal{B}_{\mathrm{matched}}|}
\sum_{i\in\mathcal{B}_{\mathrm{matched}}}d_{ibcm} -
\frac{1}{|\mathcal{I}_{\mathrm{matched}}|}
\sum_{i\in\mathcal{I}_{\mathrm{matched}}}d_{ibcm}.
$$

In the core $5\times5$ cluster DGP, $|\mathcal{B}_{\mathrm{matched}}|=|\mathcal{I}_{\mathrm{matched}}|=5$, and boundary and interior nodes are matched by lattice row.

The interpretation is:

- $\mathrm{DiD}_{bcm}>0$: method $m$ incurs additional reconstruction error near the true cluster boundary, relative to both matched interior nodes and nodewise smoothing;
- $\mathrm{DiD}_{bcm}=0$: there is no additional boundary-specific loss relative to nodewise smoothing;
- $\mathrm{DiD}_{bcm}<0$: method $m$ performs relatively better at the boundary than at the matched interior nodes.

For each cluster cell $c$, let

$$
\mathcal{S}^{(\mathrm{boundary})}_{c,m}=
\left\lbrace b:S_{bc,m}S_{bc,\mathrm{nodewise}}=1\right\rbrace,
\qquad
n^{(\mathrm{boundary,pair})}_{c,m}=
|\mathcal{S}^{(\mathrm{boundary})}_{c,m}|.
$$

Across jointly successful replications, the estimated mean boundary effect is

$$
\widehat{\mathrm{DiD}}_{c,m} =
\frac{1}{n^{(\mathrm{boundary,pair})}_{c,m}}
\sum_{b\in\mathcal{S}^{(\mathrm{boundary})}_{c,m}}
\mathrm{DiD}_{bcm}.
$$

Its Monte Carlo standard error is calculated from the replication-specific paired effects:

$$
\mathrm{MCSE}\left(\widehat{\mathrm{DiD}}_{c,m}\right) =
\frac{
\mathrm{SD}\left\lbrace \mathrm{DiD}_{bcm}:
b\in\mathcal{S}^{(\mathrm{boundary})}_{c,m}\right\rbrace
}{
\sqrt{n^{(\mathrm{boundary,pair})}_{c,m}}
}.
$$

The primary boundary analysis concerns the network-weighted estimator. Boundary and interior ISE values may additionally be reported separately as descriptive summaries.


### 6.8 Monte Carlo uncertainty

Monte Carlo standard errors are reported together with all estimated performance summaries using a formula appropriate to the corresponding statistic.

For the estimated MISE of method $m$ in cell $c$, the marginal MCSE is

$$
\mathrm{MCSE}\left(\widehat{\mathrm{MISE}}_{cm}\right) =
\frac{
\mathrm{SD}\left\lbrace \mathrm{AISE}_{bcm}:b\in\mathcal{S}_{cm}\right\rbrace
}{
\sqrt{n_{cm}}
}.
$$

For the primary network-versus-nodewise comparison, uncertainty is instead calculated from the paired differences $\Delta_{bcm}$, as defined in Section 6.3. The two marginal MCSEs are not combined because this would ignore the within-replication correlation between methods.

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

For a paired comparison of methods $m_1$ and $m_2$, let $\mathcal{A}_c$ be the prespecified attempted replication-ID set for block or cell $c$. Define the number of jointly successful replications as

$$
n^{\mathrm{pair}}_{c,m_1,m_2} =
\sum_{b\in\mathcal{A}_c}
S_{bc,m_1}S_{bc,m_2}.
$$

The paired mean and its MCSE are calculated only from these jointly successful replications, with $n^{\mathrm{pair}}_{c,m_1,m_2}$ in both the mean and the MCSE denominator. The attempted replication count, method-specific failure rates, jointly successful count, and characteristics of excluded datasets are reported alongside every affected comparison.

Marginal method-specific performance estimates analogously use that method's number of successful fits, rather than $|\mathcal{A}_c|$, in their mean and MCSE formulas. A mean is reported only when its relevant success count is at least one, and an MCSE only when that count is at least two; otherwise the statistic is recorded as undefined together with the count and failure records.

The prespecified fit-failure threshold is 5% per method and design cell. If any fitted core method exceeds this threshold in any core cell during the pilot, the main simulation is paused while the implementation and numerical diagnostics are investigated. In the main simulation, failed fits are never replaced. If a final failure rate exceeds 5%, the affected paired performance result is explicitly labelled as conditional on joint fitting success and interpreted together with the failure analysis.


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

Let $c$ index a core simulation cell. The replication-count rule is based only
on the primary network-REML-versus-nodewise comparison in the 16 core cells;
the separate $\alpha=0$ diagnostic does not determine the core replication
count. Let

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
\Delta_{bc,\mathrm{network}}
$$

denote the pilot mean paired improvement, and let $s^{\mathrm{pilot}}_{\Delta,c}$ denote the corresponding sample standard deviation. The pilot nodewise MISE used below is calculated only over successful nodewise fits in cell $c$. If fewer than two jointly successful pairs are available, or if no successful nodewise fit is available, the replication-count calculation is undefined and the pilot is paused for investigation.

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

Pilot replications are not included in the final performance estimates. If the DGP or fitting procedure is changed after the pilot, the relevant pilot checks and the replication calculation are repeated before the main simulation is started.

### 8.2 Runner implementation contract

Scenario identifiers are generated from a frozen scenario table and have the canonical form

```text
{block}__n{N}__{truth}__a{alpha}__{graph}__s{sigma}
```

with decimal points encoded as `p`, for example `core__n25__coordinate__a0p9__rook__s0p5`. Replication identifiers are positive integers within a simulation phase. The pair `(phase, replication_id)` identifies a random-number stream. Each attempted fit has a unique `fit_id`: it equals the method ID for M0, M1, M3, and M4, and includes the two candidate offsets for an M2 candidate. The tuple `(phase, scenario_id, replication_id, fit_id)` is the primary key of a fit-result record.

| Block | Scenario outputs | Attempted IDs | Required methods or operations | Prespecified reuse |
|---|---:|---:|---|---|
| Known-answer sanity | 1 | 50 | M1, M3 | Same dataset for both methods |
| Core factorial | 16 | Common $B$ | M0, M1, M3, M4 | Eight independent datasets per ID; rook/queen share truth and observations; M0/M3/M4 computed once per dataset |
| Uninformative-graph control | 2 | Same $B$ IDs | M0, M1, M3, M4 | Rook/queen share one dataset per ID; graph-independent methods computed once |
| Rewired-graph control | 2 | 200 per selected cell | Additional network fits with fixed rewired graph | Same data and stream positions as corresponding correct-rook core fits; paired against frozen rook results |
| Oracle pilot | 4 | First 50 pilot IDs | M1 plus complete M2 search | Same dataset and REML start within each cell; no final-performance contribution |
| Oracle main | 4 | First 100 main IDs | Reused M1 plus complete M2 search | REML and all oracle candidates use the identical dataset within each replication |
| Larger graph timing | 2 | 5 timing IDs per cell | M1, M3 | Independent timing data; timing results fix $B_{100}$ before performance results are inspected |
| Larger graph main | 2 | Fixed $B_{100}\in\{100,50,25\}$ per cell | M1, M3 | Same data for network and nodewise within each cell; no M0 or M4 |
| Clustering | Not executed | Not applicable | Not applicable | Archived plan only |

Every attempted fit record stores at least: `phase`, `scenario_id`, `replication_id`, `dataset_id`, `method_id`, `fit_id`, graph identifier and adjacency hash, master-stream identifier, component-substream identifiers, start and end timestamps, elapsed time, success indicator, convergence status, warning count and messages, error class and message, prediction-finiteness flag, AISE and node-level summaries when available, log smoothing parameters, EDF and attainable ranks, diagnostic flags, software versions, and a hash of the frozen specification and scenario table. Candidate-level M2 records additionally store offsets, fixed smoothing parameters, AISE, success diagnostics, and the deterministic selection result.

The checkpoint key is `(phase, replication_id)`. A checkpoint is published atomically only after every record required by the frozen scenario table for that phase and replication has been written, including explicit failure records. A manifest row stores the checkpoint key, expected and observed record counts, file hash, completion time, and status. Aggregation uses the record primary key `(phase, scenario_id, replication_id, fit_id)`; exact duplicate rows with identical content hashes are collapsed deterministically, whereas conflicting duplicates cause aggregation to stop with an error. No result is selected by file modification time.

### 8.3 Random numbers, processing order, and reproducibility

The simulation uses `RNGkind("L'Ecuyer-CMRG")`. A seed ledger is generated
serially from a fixed base seed. The pilot, known-answer check, and core
experiment use disjoint stream ranges: positions 1--50 for the pilot,
51--100 for the known-answer check, and 101--300 for the core run. The
separate $\alpha=0$ diagnostic deliberately reuses the core stream positions
while applying its own fixed data-generating mechanism and is analysed
separately.

The rewired-graph control also reuses the core stream positions
101--300 for its corresponding correct-rook cells. This is intentional: it
holds the generated truth and observations fixed while changing only the graph
supplied to the additional network fit.

Within each replication stream, `parallel::nextRNGSubStream()` assigns component substreams in the fixed order `coordinate_unstructured`, `cluster_unstructured`, `observation_error`, `clustering`, and `method_auxiliary`. Components that require several draws use a documented fixed internal order by graph size, DGP, temporal component, and noise-free dataset identifier. No seed is derived from task order, worker number, elapsed time, or character hashing. Consequently, results do not depend on the order in which cells or methods are executed.

Common random numbers are used for paired comparisons. Within replication $b$:

- within each DGP and graph size, the same centred and standardised structured and unstructured coefficient components are reused across all values of $\alpha$, with only the prespecified weights $\sqrt{\alpha}$ and $\sqrt{1-\alpha}$ changed;
- the same resulting truth is reused across the paired rook and queen cells;
- the same standard-normal error realisation is reused across noise levels and multiplied by the corresponding value of $\sigma$;
- every method fitted within a cell receives exactly the same true functions and observed dataset;
- the network and nodewise methods are therefore compared on identical data;
- raw, nodewise, and pooled estimators, which do not depend on the estimator graph, are computed once per simulated dataset and their results are reused in the paired rook-versus-queen analysis.

The simulation is processed in replication-major order within every
simulation phase. For each main-phase replication, all required core cells are
attempted before the replication-level checkpoint is published.

For each replication:

1. generate all random components required for replication $b$;
2. run every prespecified cell required in the current phase for replication $b$;
3. fit all required methods and calculate all performance measures for those cells;
4. save a replication-level result object containing either a success record or a failure record for every attempted fit, and mark replication $b$ as complete;
5. proceed to replication $b+1$.

A replication is included in the aggregated results only after all fits required by the frozen scenario table for the current phase have been attempted and a success or failure record has been saved for each one. If execution stops before this complete result object is saved, the interrupted replication is recomputed from its original master seed when the run resumes. By contrast, a recorded model-fitting failure is a completed attempted fit and does not cause the dataset to be regenerated, the seed to be changed, or the fit to be repeatedly retried.

Results are checkpointed by replication rather than only by scenario. A completion manifest records the identifiers of replication objects that have been completely saved, irrespective of whether individual fits succeeded. This ensures that an interrupted run leaves the same number of attempted replications in every core cell and remains directly evaluable.

Pilot and main-simulation seeds are stored separately. Scenario definitions, random-number generation, data generation, model fitting, performance evaluation, and result aggregation are implemented as separate functions.

---

## 9. Analysis and Presentation

### 9.1 Estimand-to-output and DGP-validation tables

The reporting contract is fixed as follows:

| Target | Measure | Output |
|---|---|---|
| Absolute average reconstruction | MISE and marginal MCSE | Core method-by-cell table |
| Individual-node reconstruction | Monte Carlo mean and median of $\mathrm{MedISE}$ and $\mathrm{WorstISE}$ | Secondary core table |
| Network benefit over nodewise | Paired mean AISE improvement and paired MCSE | Primary comparison table and figure |
| Neighbourhood-density effect | Paired rook-minus-queen AISE difference and MCSE | Density table and figure |
| False-graph cost | Paired rewired-minus-rook AISE difference, paired MCSE, rMISE, and win rate | Supplementary two-cell control table |
| Relative reconstruction | Joint-success rMISE; median and IQR of RI | Secondary core table and heatmap |
| Cluster-boundary loss | Matched boundary DiD and paired MCSE | Cluster table and figure |
| Restricted tuning potential | Restricted oracle gap and paired MCSE | Separate oracle table, after the prespecified run is complete |
| Scalability | Runtime, fit diagnostics, MISE, and paired improvement | Separate $N=100$ table, after the prespecified run is complete |
| Cluster preservation | Not executed | No reported output |

A separate deterministic design audit is conducted before estimator-performance
summaries are inspected. It records `SNR_overall`, `SNR_dev`, and `RPE_dev`
for each distinct generated truth configuration.

### 9.2 Main performance tables

For every core combination of truth structure, structured-signal proportion,
neighbourhood density, and noise level, the following will be reported:

- estimated MISE and its marginal MCSE for every available method, with raw and pooled results labelled as diagnostic benchmarks;
- mean paired AISE improvement $\widehat{\Delta}_{c,m}$ and its paired MCSE, with network-REML versus nodewise identified as the primary method comparison;
- direct paired rook-versus-queen AISE difference $\widehat{\Delta}^{(\mathrm{density})}_{c}$ and its paired MCSE for each otherwise identical design cell;
- aggregate rMISE relative to nodewise smoothing;
- Monte Carlo mean and median of the replication-specific median-node and worst-node ISE;
- median and interquartile range of replication-specific RI;
- mean, median, and selected quantiles of computation time;
- numbers of successful, failed, warning-producing, and non-converged fits;
- numbers of attempted and jointly successful replications for every paired comparison, with conditional-on-joint-success labels where the 5% failure threshold is exceeded;
- retained fit diagnostics, warnings, and failures.

For cluster-structured cells, matched boundary DiD and its paired MCSE are
additionally reported. The rewired-graph control is reported in a separate
two-row paired-comparison table. The oracle and $N=100$ results are each
reported in a separate table and are not pooled with the core summaries.

The known-answer sanity check is reported separately with its attempted, method-specific successful, and jointly successful counts; mean paired improvement; paired MCSE; proportion of jointly successful replications with positive improvement; and pass/fail status. It is not pooled with the core simulation results.

### 9.3 Planned figures

The primary presentation includes a faceted figure of the relative MISE reduction of
network-weighted smoothing relative to nodewise smoothing across the 16 core
cells, and a figure of the matched boundary difference-in-differences for the
cluster-structured cells. The known-answer check, the separate $\alpha=0$
diagnostic, and the rewired-graph control are reported in tables. Any $N=100$
or oracle figure is labelled separately; no clustering figure is planned.

---

## 10. Pre-Main Checks

Before any estimator-performance result is inspected, the following checks are
completed and documented:

1. the known-answer sanity check satisfies all prespecified criteria;
2. the relative unpenalised $L^2$ projection error of the deviations is below
   $2\%$;
3. the two core noise levels produce meaningfully different overall and
   deviation SNR values;
4. the fixed rewired graph satisfies its connectivity and degree-preservation
   requirements;
5. the core replication count is fixed by the paired-MCSE rule;
6. the fit-failure rate of every fitted core method is at most $5\%$ in every
   core cell;
7. the oracle candidate with offsets $(0,0)$ reproduces the REML predictions to
   within $10^{-8}$;
8. the 50-ID oracle pilot fixes the directional coarse grid before the oracle
   main run; and
9. five timing replications per $N=100$ cell fix $B_{100}$ under the eight-hour
   budget rule.

The frozen scenario table, seed ledger, expected record counts, atomic
checkpoint publication, resume behaviour, and deterministic duplicate
detection are validated before the corresponding main run.

---

## 11. Prioritisation Under Computational Constraints

Analyses are prioritised as follows:

1. known-answer sanity check and core benchmark;
2. $\alpha=0$ diagnostic, pooled benchmark, cluster-boundary diagnostic, and
   rewired-graph control;
3. two fixed $N=100$ cells; and
4. restricted oracle diagnostic.

The $N=100$ block precedes the oracle because node count is a distinct design
dimension. If further reduction is necessary, the restricted oracle is reduced
first; the $N=100$ block is not replaced by an ad hoc smaller replication count
when its timing rule declares it infeasible. Clustering remains optional future
work.

---

## 12. Summary of the Simulation Design

| ADEMP component | Specification |
|---|---|
| Aim | Assess the reconstruction benefit of graph-informed smoothing for node-specific curves at observed nodes |
| Data | $N=25$, $5\times5$ lattice, $T=50$, one noisy curve per node, Gaussian independent errors |
| Truth | Fixed $\mu(t)$ and three fixed temporal basis functions; coordinate-smooth or cluster-structured coefficient fields generated independently of the estimator graph |
| Structured signal | $\alpha=0.5,0.9$ with fixed marginal coefficient variance |
| Graph density | Rook and queen adjacency applied to the same truth and observations |
| Noise level | $\sigma=0.2,0.5$; overall and deviation SNR are checked in the pilot |
| Primary comparison methods | Network-REML and nodewise smoothing |
| Diagnostic benchmarks | Raw observations and pooled mean smoothing |
| Fixed fitting settings | Cubic P-splines with first-order difference penalties; $K_{\mathrm{int}}=10$, $K_{\mathrm{dev}}=15$; Gaussian identity-link models; REML; common 50-point prediction grid |
| Primary estimand | True node-specific functions $f_i(t)$ |
| Primary measures | MISE and paired AISE improvement over nodewise smoothing with paired MCSE |
| Secondary measures | Median-node and worst-node ISE, descriptive rook-versus-queen comparison, joint-success rMISE, replication-wise win rate, matched boundary DiD, runtime, warnings, and failures |
| Core design | 16 factorial cells with a common pilot-determined $B\geq200$; failed fits are recorded and never replaced |
| Failure policy | Success requires no error, no reported non-convergence, and finite predictions; pilot failure threshold 5% per fitted core method and cell; paired MCSE uses the jointly successful count |
| Known-answer check | Coordinate-smooth truth, $\alpha=1$, rook graph, $\sigma=0.5$, $N=25$, $T=50$, and $B_{\mathrm{sanity}}=50$, with three prespecified pass criteria |
| Negative control | $\alpha=0$ uninformative-graph diagnostic at $\sigma=0.5$ for rook and queen graphs |
| False-graph ablation | One fixed connected degree-preserving rewired rook graph; coordinate-smooth and cluster-structured truth at $\alpha=0.9$, $\sigma=0.5$; paired against the correct rook graph |
| Larger graph | Two fixed $N=100$ rook-lattice cells; timing-based $B_{100}\in\{100,50,25\}$; reported separately |
| Oracle | Restricted coarse-to-fine diagnostic in four selected $N=25$ cells; reported separately |
| Clustering | Optional secondary analysis |
| Reproducibility | Frozen scenario table; `L'Ecuyer-CMRG` stream/substream ledger; replication-major atomic checkpoints; explicit failure records; hash-based deterministic deduplication |

## 13. Future Work

A future extension may study prediction of complete functions at unobserved graph nodes and compare an adapted network method with functional kriging, for example through `SpatFD`. Because this is a spatial-prediction problem rather than reconstruction at observed nodes, it requires a separate node-holdout ADEMP specification with its own data-generating assumptions, comparators, tuning rules, performance measures, replication count, and failure policy. It is not part of the present simulation study.

Temporally correlated or heteroscedastic errors, a $\sigma=0.1$
high-information setting, additional graph perturbations beyond the fixed
rewiring, and the clustering task remain future work. The optional
$\sigma=0.1$ setting is not part of the current design because the existing
$\sigma=0.2$ and $\sigma=0.5$ levels already span a factor of $6.25$ in
information per node under the iid common-grid error model.

## References

Chung, F. R. K. (1997). *Spectral Graph Theory*. CBMS Regional Conference Series in Mathematics, Vol. 92. American Mathematical Society. <https://bookstore.ams.org/cbms-92>

Greven, S., and Scheipl, F. (2017). A general framework for functional regression modelling. *Statistical Modelling*, 17(1–2), 1–35. <https://doi.org/10.1177/1471082X16681317>

Hubert, L., and Arabie, P. (1985). Comparing partitions. *Journal of Classification*, 2, 193–218. <https://doi.org/10.1007/BF01908075>

Morris, T. P., White, I. R., and Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38, 2074–2102. <https://doi.org/10.1002/sim.8086>

Wood, S. N., Pya, N., and Säfken, B. (2016). Smoothing parameter and model selection for general smooth models. *Journal of the American Statistical Association*, 111, 1548–1575. <https://doi.org/10.1080/01621459.2016.1180986>
