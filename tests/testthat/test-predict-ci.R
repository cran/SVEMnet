skip_on_cran()
skip_if_not_installed("SVEMnet")

test_that("predict returns percentile CI when requested", {
  set.seed(1)
  n  <- 60
  X1 <- rnorm(n); X2 <- rnorm(n)
  y  <- 1 + 0.8*X1 - 0.4*X2 + rnorm(n, 0, 0.3)
  d  <- data.frame(y, X1, X2)

  fit <- SVEMnet::SVEMnet(y ~ X1 + X2, d, nBoot = 40, glmnet_alpha = 1)
  pr  <- predict(fit, d[1:5, ], debias = FALSE, agg = "mean", interval = TRUE, level = 0.95)

  expect_true(is.list(pr))
  expect_true(all(c("fit","lwr","upr") %in% names(pr)))
  expect_equal(length(pr$fit), 5)
  expect_equal(length(pr$lwr), 5)
  expect_equal(length(pr$upr), 5)
})
