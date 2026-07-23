test_that("unnamed curves are matched to nodes by position", {
  curves <- list(10, 20, 30)
  node_names <- c("A", "B", "C")

  result <- align_curves_to_nodes(curves, node_names)

  expect_equal(result, list(A = 10, B = 20, C = 30))
  expect_named(result, node_names)
})


test_that("named curves are reordered to match graph nodes", {
  curves <- list(
    A = 10,
    B = 20,
    C = 30
  )

  result <- align_curves_to_nodes(
    curves,
    node_names = c("C", "A", "B")
  )

  expect_named(result, c("C", "A", "B"))
  expect_equal(unname(result), list(30, 10, 20))
})


test_that("curves already in node order remain unchanged", {
  curves <- list(
    A = 10,
    B = 20,
    C = 30
  )

  result <- align_curves_to_nodes(
    curves,
    node_names = c("A", "B", "C")
  )

  expect_identical(result, curves)
})


test_that("empty curve names are rejected", {
  curves <- list(10, 20, 30)
  names(curves) <- c("A", "", "C")

  expect_error(
    align_curves_to_nodes(curves, c("A", "B", "C")),
    "every curve must have a non-empty name"
  )
})


test_that("missing curve names are rejected", {
  curves <- list(
    A = 10,
    B = 20,
    C = 30
  )
  names(curves)[2] <- NA_character_

  expect_error(
    align_curves_to_nodes(curves, c("A", "B", "C")),
    "every curve must have a non-empty name"
  )
})


test_that("duplicate curve names are rejected", {
  curves <- list(
    A = 10,
    A = 20,
    C = 30
  )

  expect_error(
    align_curves_to_nodes(curves, c("A", "B", "C")),
    "must be unique"
  )
})


test_that("missing and unknown nodes are reported", {
  curves <- list(
    A = 10,
    B = 20,
    D = 30
  )

  expect_error(
    align_curves_to_nodes(curves, c("A", "B", "C")),
    regexp = "Missing curves for graph nodes: C"
  )

  expect_error(
    align_curves_to_nodes(curves, c("A", "B", "C")),
    regexp = "Curves without matching graph nodes: D"
  )
})
