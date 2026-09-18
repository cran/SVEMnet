# Conditional Gaussian EDF regression tests.
skip_on_cran()

.edf_data <- function(seed = 18, n = 36) {
  set.seed(seed)
  x <- matrix(rnorm(n * 4), n, 4)
  data.frame(y = 8 + x %*% c(3, -2, 1, .5) + rnorm(n), x)
}

test_that("corrected production EDF is invariant to response units", {
  dat <- .edf_data(); big <- dat; big$y <- 10 * dat$y
  for (relaxed in c(FALSE, TRUE)) {
    set.seed(2); pa <- SVEMnet(y ~ ., dat, nBoot = 6, glmnet_alpha = .5,
                              complexity = "edf", relaxed = relaxed)
    set.seed(2); pb <- SVEMnet(y ~ ., big, nBoot = 6, glmnet_alpha = .5,
                              complexity = "edf", relaxed = relaxed)
    expect_equal(pa$coef_matrix, pb$coef_matrix / 10, tolerance = 1e-8)
    expect_equal(pa$best_edfs, pb$best_edfs, tolerance = 1e-8)
  }
})

test_that("weighted nonuniform and zero-penalty EDF matches conditional trace", {
  d <- .edf_data(n = 60); X <- as.matrix(d[-1]); y <- d$y
  w <- seq_len(60); w <- w * 60 / sum(w)
  pf <- c(0, .5, 1, 2)
  ctx <- .svem_edf_context(X, y, w, list(penalty.factor = pf))
  fit <- do.call(glmnet::glmnet, .glmnet_prepare_call_args(
    list(x = X, y = y, weights = w, alpha = .5, penalty.factor = pf,
         nlambda = 30, thresh = 1e-13), glmnet::glmnet))
  co <- as.matrix(stats::coef(fit)); j <- ncol(co) - 2L
  act <- .svem_active(co[, j]); G <- crossprod(ctx$z[, act, drop = FALSE])
  ridge <- 60 * fit$lambda[j] * .5 / ctx$y_sd
  expected <- 1 + sum(diag(solve(G + diag(ridge * ctx$penalty[act], length(act)), G)))
  actual <- .svem_conditional_edf(ctx, co[, j, drop = FALSE], fit$lambda[j], .5)
  expect_equal(actual, expected, tolerance = 1e-9)
  # Verify the scale term against glmnet coefficient equations, not just itself.
  xm <- colSums(X * w) / 60; xc <- sweep(X, 2, xm, "-")
  sx <- sqrt(colSums(xc^2 * w) / 60); z <- sweep(xc, 2, sx, "/")
  theta <- co[1 + act, j] * sx[act]
  rhs <- as.vector(crossprod(z[, act, drop = FALSE], w * (y - sum(w * y) / 60))) -
    60 * fit$lambda[j] * .5 * ctx$penalty[act] * sign(theta)
  expect_equal(as.vector((G + diag(ridge * ctx$penalty[act], length(act))) %*% theta), unname(rhs), tolerance = 2e-6)
  # Singular zero-ridge limit must return a rank, never 0/0.
  Z <- cbind(X[, 1], X[, 1]); colnames(Z) <- c("a", "b")
  cc <- .svem_edf_context(Z, y, w)
  expect_equal(.svem_conditional_edf(cc, matrix(c(0, 1, 1), 3), 0, .5), 2)
})

test_that("bounded fits and invalid configurations do not invent EDF", {
  dat <- .edf_data()
  expect_error(SVEMnet(y ~ ., dat, nBoot = 2, complexity = "edf", lower.limits = 0), "bounds")
  fit <- SVEMnet(y ~ ., dat, nBoot = 2, complexity = "support", lower.limits = 0)
  expect_true(all(is.na(fit$best_edfs)))
})
