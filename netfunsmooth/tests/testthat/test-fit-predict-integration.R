test_that("predict.netf_fit preserves curve names", {
  t_grid <- seq(0, 1, length.out = 20)

  values <- rbind(
    A = sin(2 * pi * t_grid),
    B = 0.8 * sin(2 * pi * t_grid + 0.2),
    C = 1.2 * sin(2 * pi * t_grid - 0.2)
  )

  curves <- tf::tfd(values, arg = t_grid)

  graph <- igraph::make_ring(3)
  igraph::V(graph)$name <- c("A", "B", "C")

  fit <- netf_smooth(
    curves = curves,
    graph = graph,
    bs.int = list(
      bs = "ps",
      k = 6,
      m = c(2, 1)
    ),
    bs.yindex = list(
      bs = "ps",
      k = 6,
      m = c(2, 1)
    )
  )

  pred <- predict(fit)

  expect_s3_class(fit, "netf_fit")
  expect_s3_class(pred, "tfd")
  expect_identical(names(pred), c("A", "B", "C"))
  expect_length(pred, 3L)
  expect_equal(tf::tf_arg(pred), t_grid)

  pred_matrix <- do.call(
    rbind,
    tf::tf_evaluate(pred, arg = t_grid)
  )

  expect_equal(dim(pred_matrix), c(3L, 20L))
  expect_equal(
    unname(pred_matrix),
    unname(stats::fitted(fit$model))
  )
})
