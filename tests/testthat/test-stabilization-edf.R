test_that("conditional lasso EDF equals projection rank under collinearity", {
  x <- seq(-1, 1, length.out = 20)
  X <- cbind(x, x, x^2)
  y <- 2 + 3*x + sin(seq_along(x))
  w <- seq_along(x)
  context <- .svem_edf_context(X, y, w)
  coefficients <- matrix(c(2, 1, 1, .5), ncol = 1)
  # Weighted, centered projection: duplicate columns span only two slopes.
  Z <- cbind(1, X) * sqrt(w)
  expected <- qr(Z)$rank
  for (gamma in c(.2, .6, 1)) {
    expect_equal(.svem_conditional_edf(context, coefficients, .1, 1, gamma),
                 as.numeric(expected))
  }
  expect_equal(.svem_solver_support(coefficients[, 1]), 4L)
})

test_that("explicit private tolerance retains the research convention", {
  beta <- c(3, 1, 1e-9, 0)
  expect_identical(.svem_active(beta), c(1L, 2L))
  expect_identical(.svem_active(beta, 1e-7), 1L)
  expect_identical(.svem_active(c(3, 1e5, .001), 1e-7), 1L)
  for (bad in list(-1, NA_real_, Inf, c(0, 1), TRUE)) {
    expect_error(.svem_active(beta, bad), "tol")
  }
  X <- cbind(1:12, (1:12)^2, sin(1:12))
  context <- .svem_edf_context(X, 2 + 1:12, rep(1, 12))
  expect_equal(.svem_conditional_edf(context, matrix(beta), .1, 1), 3)
  expect_equal(.svem_conditional_edf(context, matrix(beta), .1, 1, 1, 1e-7), 2)
})
