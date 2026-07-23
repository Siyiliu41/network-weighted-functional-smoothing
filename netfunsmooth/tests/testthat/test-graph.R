test_that("igraph objects are converted to named neighbour lists", {
  graph <- igraph::make_ring(3)
  igraph::V(graph)$name <- c("A", "B", "C")

  result <- graph_to_nb(graph)

  expect_named(result, c("A", "B", "C"))
  expect_equal(result$A, c("B", "C"))
  expect_equal(result$B, c("A", "C"))
  expect_equal(result$C, c("A", "B"))
})


test_that("unnamed igraph vertices receive numeric names", {
  graph <- igraph::make_ring(3)

  result <- graph_to_nb(graph)

  expect_named(result, c("1", "2", "3"))
  expect_equal(length(result), 3L)
})


test_that("adjacency matrices are converted correctly", {
  graph <- matrix(
    c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("A", "B", "C"),
      c("A", "B", "C")
    )
  )

  result <- graph_to_nb(graph)

  expected <- list(
    A = "B",
    B = c("A", "C"),
    C = "B"
  )

  expect_equal(result, expected)
})


test_that("matrix columns are reordered by node name", {
  graph <- matrix(
    c(
      0, 0, 1,
      1, 1, 0,
      0, 0, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("A", "B", "C"),
      c("C", "A", "B")
    )
  )

  result <- graph_to_nb(graph)

  expected <- list(
    A = "B",
    B = c("A", "C"),
    C = "B"
  )

  expect_equal(result, expected)
})


test_that("matrix diagonal entries are ignored", {
  graph <- matrix(
    c(
      1, 1,
      1, 1
    ),
    nrow = 2,
    dimnames = list(
      c("A", "B"),
      c("A", "B")
    )
  )

  result <- graph_to_nb(graph)

  expect_equal(result, list(A = "B", B = "A"))
})


test_that("invalid adjacency matrices are rejected", {
  expect_error(
    graph_to_nb(matrix(1:6, nrow = 2)),
    "square adjacency matrix"
  )

  graph_with_na <- matrix(
    c(0, NA, NA, 0),
    nrow = 2
  )

  expect_error(
    graph_to_nb(graph_with_na),
    "must not contain missing values"
  )

  asymmetric_graph <- matrix(
    c(
      0, 1,
      0, 0
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_error(
    graph_to_nb(asymmetric_graph),
    "must be symmetric"
  )
})


test_that("matrix row and column names must describe the same nodes", {
  graph <- diag(3)
  rownames(graph) <- c("A", "B", "C")
  colnames(graph) <- c("A", "B", "D")

  expect_error(
    graph_to_nb(graph),
    "must contain the same nodes"
  )
})


test_that("unsupported graph classes produce an informative error", {
  expect_error(
    graph_to_nb(data.frame(x = 1:3)),
    "No `graph_to_nb\\(\\)` method"
  )
})
