make_test_curves <- function(n_curves = 3L, n_grid = 8L) {
  t_grid <- seq(0, 1, length.out = n_grid)

  values <- vapply(
    seq_len(n_curves),
    function(i) i * sin(2 * pi * t_grid),
    numeric(n_grid)
  )

  values <- t(values)

  tf::tfd(values, arg = t_grid)
}


make_test_graph <- function(n_nodes = 3L) {
  graph <- matrix(
    0,
    nrow = n_nodes,
    ncol = n_nodes
  )

  if (n_nodes > 1L) {
    for (i in seq_len(n_nodes - 1L)) {
      graph[i, i + 1L] <- 1
      graph[i + 1L, i] <- 1
    }
  }

  node_names <- LETTERS[seq_len(n_nodes)]
  dimnames(graph) <- list(node_names, node_names)

  graph
}


test_that("unsupported curve classes produce an informative error", {
  graph <- make_test_graph(3)

  expect_error(
    netf_smooth(
      curves = list(1, 2, 3),
      graph = graph
    ),
    "No `netf_smooth\\(\\)` method"
  )
})


test_that("tfd inputs are dispatched to the common fitting function", {
  curves <- make_test_curves(3)
  graph <- make_test_graph(3)

  names(curves) <- c("A", "B", "C")

  local_mocked_bindings(
    fit_netf_smooth = function(curves, graph, ...) {
      list(
        curves = curves,
        graph = graph,
        extra = list(...)
      )
    },
    .package = "netfunsmooth"
  )

  result <- netf_smooth(
    curves = curves,
    graph = graph,
    test_argument = 42
  )

  expect_s3_class(result$curves, "tfd")
  expect_identical(names(result$curves), c("A", "B", "C"))
  expect_identical(result$graph, graph)
  expect_identical(result$extra$test_argument, 42)
})


test_that("curve and graph node counts must agree", {
  curves <- make_test_curves(2)
  graph <- make_test_graph(3)

  expect_error(
    netf_smooth(
      curves = curves,
      graph = graph
    ),
    "one curve for each graph node"
  )
})


test_that("curve names are aligned with graph node names", {
  curves <- make_test_curves(3)
  names(curves) <- c("C", "A", "B")

  graph <- make_test_graph(3)

  local_mocked_bindings(
    register_nb_list = function(nb_list) "test-key",
    build_pffr_formula = function(key, k_node) stats::as.formula("Y ~ 1"),
    .package = "netfunsmooth"
  )

  # This test checks alignment indirectly before pffr is reached.
  aligned <- netfunsmooth:::align_curves_to_nodes(
    curves,
    node_names = c("A", "B", "C")
  )

  expect_identical(names(aligned), c("A", "B", "C"))
})
