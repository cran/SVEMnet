# Tests for the experimental complexity = "edf" option (gamma-blended
# effective degrees of freedom in the wAIC/wBIC penalty).
skip_on_cran()
skip_if_not_installed("SVEMnet")

.edf_toy_data <- function(n = 30, seed = 101) {
  set.seed(seed)
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  y  <- 1 + 2 * X1 - 1.5 * X2 + 0.5 * X3 + 1.2 * X1 * X2 + rnorm(n, sd = 0.5)
  data.frame(y, X1, X2, X3)
}

test_that("complexity = 'support' is identical to the default", {
  dat <- .edf_toy_data()
  set.seed(7); fit_def <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 25)
  set.seed(7); fit_sup <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 25,
                                  complexity = "support")
  expect_identical(fit_def$parms, fit_sup$parms)
  expect_identical(fit_def$best_lambdas, fit_sup$best_lambdas)
  expect_identical(fit_sup$complexity, "support")
  expect_true(all(is.na(fit_sup$best_edfs)))
})

test_that("edf reduces exactly to support for pure lasso (alpha = 1)", {
  dat <- .edf_toy_data(seed = 102)
  # With only alpha = 1 in the grid, the blended edf equals the support count
  # at every path point, so the entire selection must be identical.
  set.seed(8); fit_sup <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 25,
                                  glmnet_alpha = 1, complexity = "support")
  set.seed(8); fit_edf <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 25,
                                  glmnet_alpha = 1, complexity = "edf")
  expect_identical(fit_sup$parms, fit_edf$parms)
  expect_identical(fit_sup$best_lambdas, fit_edf$best_lambdas)
  expect_identical(fit_sup$best_relax_gammas, fit_edf$best_relax_gammas)
  expect_identical(fit_edf$complexity, "edf")
})

test_that("edf runs end-to-end with alpha < 1 and stays within bounds", {
  dat <- .edf_toy_data(seed = 103)
  set.seed(9)
  fit <- SVEMnet(y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2), dat, nBoot = 40,
                 glmnet_alpha = c(0.5), complexity = "edf")

  expect_true(all(is.finite(fit$parms)))
  edfs <- fit$best_edfs
  expect_true(all(is.finite(edfs)))
  # blended edf includes the intercept and can never exceed the full
  # parameter count; it must be at least 1 (intercept only)
  expect_true(all(edfs >= 1 - 1e-8))
  expect_true(all(edfs <= fit$nparm + 1e-8))
  # with alpha = 0.5 the ridge component shrinks, so at least some selected
  # path points should have non-integer effective df (evidence the blended
  # trace was actually used rather than the count)
  expect_true(any(abs(edfs - round(edfs)) > 1e-6))
  # predictions work as usual
  p <- predict(fit, dat)
  expect_true(all(is.finite(p)))
})

test_that("edf is reproducible under set.seed and consumes the same RNG stream", {
  dat <- .edf_toy_data(seed = 104)
  set.seed(11); f1 <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 20,
                              complexity = "edf")
  set.seed(11); f2 <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 20,
                              complexity = "edf")
  expect_identical(f1$parms, f2$parms)

  # the edf computation is deterministic, so a support fit with the same seed
  # must see the identical bootstrap weights: per-replicate n_eff diagnostics
  # (a pure function of the weights) must match between the two options
  set.seed(11); f3 <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 20,
                              complexity = "support")
  expect_identical(summary(f1$diagnostics$n_eff_summary),
                   summary(f3$diagnostics$n_eff_summary))
})

test_that("binomial falls back to support with a warning", {
  set.seed(12)
  n <- 90
  db <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  db$y <- rbinom(n, 1, plogis(0.9 * db$x1 - 0.6 * db$x2))
  expect_warning(
    fit_b <- SVEMnet(y ~ x1 + x2, db, nBoot = 10, family = "binomial",
                     complexity = "edf"),
    "gaussian"
  )
  expect_identical(fit_b$complexity, "support")
})

test_that("whole-model test refuses to propagate complexity", {
  dat <- .edf_toy_data(seed = 105)
  expect_warning(
    svem_significance_test_parallel(
      y ~ X1 + X2 + X3, dat,
      nPoint = 100, nSVEM = 1, nPerm = 12, nBoot = 10,
      nCore = 1, seed = 5, verbose = FALSE,
      complexity = "edf"
    ),
    "support-count"
  )
})
