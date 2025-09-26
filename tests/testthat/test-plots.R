skip_on_cran()
skip_if_not_installed("SVEMnet")

test_that("basic plot methods work", {
  set.seed(3)
  n <- 50
  X1 <- rnorm(n); X2 <- rnorm(n)
  y  <- 1 + 0.7*X1 - 0.3*X2 + rnorm(n, 0, 0.2)
  d  <- data.frame(y, X1, X2)
  fit <- SVEMnet::SVEMnet(y ~ X1 + X2, d, nBoot = 20, glmnet_alpha = 1)
  expect_silent(plot(fit))
})
