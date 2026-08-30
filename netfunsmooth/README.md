# netfunsmooth

`netfunsmooth` is the R package developed for the bachelor's thesis
*Network-Weighted Smoothing of Functional Data*. It provides a
`tf`/`tidyfun`-compatible interface for graph-regularised smoothing of curves
observed on a common grid.

## Main interface

```r
fit <- netf_smooth(
  curves,
  graph,
  graph_id = NULL,
  sandwich = "none",
  bs.int = NULL,
  bs.yindex = NULL
)

fitted_curves <- predict(fit)
```

`curves` must be a `tf::tfd` or `tf::tfb` vector containing one curve per
graph node. `netf_smooth()` returns an object of class `netf_fit`; its
underlying `refund::pffr` model is available as `fit$model`. `predict(fit)`
returns a `tf::tfd` vector on the observed common grid. When input curves are
named, their names are retained in the predictions after graph alignment.

The fit uses a functional additive mixed-model representation with an MRF
smooth over the node factor. The temporal intercept basis is controlled by
`bs.int`; the temporal marginal basis for node-specific deviations is
controlled by `bs.yindex`. Additional arguments are passed to
`refund::pffr()`.

## Graph inputs and identifier handling

`graph_to_nb()` converts supported graph representations to the named
neighbour list required by the MRF smooth.

- **Named adjacency matrix:** must be numeric, square, symmetric, and have
  matching unique row and column names. Column order is aligned to row names;
  these names determine the graph-node identifiers.
- **`igraph` object:** vertex names are used as graph-node identifiers. If no
  vertex names are present, sequential character identifiers are assigned.
- **`sf` polygon layer:** Queen contiguity is computed by default with
  `spdep::poly2nb()`. Supply `graph_id` to `netf_smooth()` as either the name
  of a polygon-identifier column or a vector containing one unique,
  non-missing, non-empty identifier per polygon in polygon row order.
  These identifiers are preserved as graph-node names. When calling
  `graph_to_nb()` directly, the corresponding argument is named `id`.

For all graph types, `netf_smooth()` verifies that the number of curves agrees
with the number of graph nodes. If curves are named, their names must agree
exactly with the derived graph-node identifiers, and the curves are reordered
to match the graph-node order. Unnamed curves are matched by position and
assigned the corresponding graph-node names; their input order must therefore
already match the graph-node order. Keep identifiers with leading zeroes as
character strings.

For example, if `polygons` is an `sf` layer with a unique character column
named `node_id`, and `curves` is named with the same identifiers:

```r
fit <- netf_smooth(
  curves = curves,
  graph = polygons,
  graph_id = "node_id"
)

nb <- graph_to_nb(polygons, id = "node_id", queen = TRUE)
```

## Package contents

- `R/`: implementation of graph conversion, fitting, and prediction;
- `man/`: generated R help pages;
- `tests/testthat/`: unit and integration tests;
- `DESCRIPTION`: package metadata and dependencies.

## Development use

From the repository root, load the package source with:

```r
devtools::load_all("netfunsmooth")
```

Run the package tests with:

```r
devtools::test("netfunsmooth")
```