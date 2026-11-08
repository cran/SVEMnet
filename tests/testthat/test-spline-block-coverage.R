test_that("unmeasured spline blocks do not discard covered recipes in the scorer", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(u = seq(0, 1, length.out = 7), v = seq(0, 1, length.out = 7))
  measured <- transform(grid[seq(1, 49, by = 3), ], day = "A", y = u + v^2)
  spec <- spline_terms(y ~ u + v, grid, blocking = "day")
  set.seed(83)
  fit <- svem_spline(spec, measured, nBoot = 2)
  pool <- transform(grid[1:2, ], day = c("A", "B"))
  scored <- .svem_spline_scores(data.frame(recipe_id = c("a", "b")), pool,
    list(y = fit), list(y = NULL), function(pred, parameters) pred,
    c(y = 1), "arith", 1e-6)
  expect_equal(scored$table$y_spline_pred[1], predict(fit, pool[1, ]))
  expect_true(is.na(scored$table$y_spline_pred[2]))
  expect_identical(scored$status$status, "partial_coverage")
  expect_match(scored$status$message, "measured block levels")
})
