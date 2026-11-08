test_that("each training-fitted mean and variance is evaluated with validation weights", {
  d <- data.frame(y = c(1, 4, 2, 4, 4, 6, 8, 1), x = rep(0:1, each = 4))
  set.seed(3)
  u <- runif(nrow(d))
  train <- -log(u)
  valid <- -log1p(-u)
  group_means <- vapply(split(seq_len(nrow(d)), d$x), function(i) {
    weighted.mean(d$y[i], train[i])
  }, numeric(1))
  intercept <- weighted.mean(d$y, train)
  residuals <- cbind(d$y - intercept, d$y - group_means[d$x + 1L])
  scale <- diff(range(d$y)) / 2
  sq <- (residuals / scale)^2
  t <- colSums(train * sq) / sum(train)
  v <- colSums(valid * sq) / sum(valid)
  scores <- log(t) + v / t
  expect_equal(.svem_forward_gaussian_score(t, v), scores)

  # Compute the former explicit objectives from the same two candidate fits.
  n <- nrow(d)
  vsse <- n * v * scale^2
  neff <- sum(valid)^2 / sum(valid^2)
  losses <- list(gaussian_likelihood = scores, wAIC = n * log(vsse / n) + 2 * 1:2,
                 wBIC = n * log(vsse / n) + log(neff) * 1:2, wSSE = vsse)
  expect_identical(which.min(scores), 2L)
  expect_identical(which.min(losses$wAIC), 1L)
  betas <- rbind(c(intercept, 0), c(group_means[1L], diff(group_means)))
  for (objective in names(losses)) {
    set.seed(3)
    fit <- svem_forward(y ~ x, d, nBoot = 1, objective = objective)
    winner <- which.min(losses[[objective]])
    expect_equal(fit$best_ks, winner)
    expect_equal(unname(fit$parms), unname(betas[winner, ]))
    expect_equal(fit$member_diagnostics$objective_value, unname(losses[[objective]][winner]))
    expect_equal(fit$member_diagnostics$validation_sse, unname(vsse[winner]))
    expect_equal(as.numeric(predict(fit, d)), unname(d$y - residuals[, winner]))
  }
})

test_that("forward routing and fitted outputs are stable across response units", {
  d <- data.frame(y = c(1, 4, 2, 4, 4, 6, 8, 1), x = rep(0:1, each = 4))
  fit <- function(data, ...) {
    set.seed(3)
    svem_forward(y ~ x, data, nBoot = 6, ...)
  }
  a <- fit(d)
  expect_identical(a, fit(d, objective = "auto"))
  b <- fit(d, objective = "gaussian_likelihood")
  expect_identical(b$objective_input, "gaussian_likelihood")
  expect_identical(a$objective_input, "auto")
  expect_identical(a$auto_decision, "gaussian_likelihood")
  expect_identical(a$coef_matrix, b$coef_matrix)
  expect_identical(a$member_diagnostics, b$member_diagnostics)
  expect_identical(a$diagnostics, b$diagnostics)
  expect_identical(predict(a, d), predict(b, d))
  for (k in c(2^-400, -2^400)) {
    f <- fit(transform(d, y = k * y))
    expect_identical(f$best_ks, a$best_ks)
    expect_equal(f$coef_matrix / k, a$coef_matrix)
    expect_equal(f$member_diagnostics$validation_sse / k^2, a$member_diagnostics$validation_sse)
    expect_equal(f$member_diagnostics$objective_value, a$member_diagnostics$objective_value)
    expect_equal(as.numeric(predict(f, d)) / k, as.numeric(predict(a, d)))
  }
  shifted <- fit(transform(d, y = y + 2^45))
  expect_identical(shifted$best_ks, a$best_ks)
  expect_equal(shifted$member_diagnostics$objective_value, a$member_diagnostics$objective_value)
  expect_equal(shifted$coef_matrix[, 2], a$coef_matrix[, 2])
})

test_that("whole-term Kish caps, constant means and interpolation have defined scores", {
  d <- data.frame(f = factor(rep(letters[1:3], each = 2)), y = c(0, 1, 4, 5, 8, 9))
  set.seed(11)
  d$z <- c(-1, 1, 1, -1, -1, 1)
  a <- svem_forward(y ~ f + z, d, nBoot = 1)
  expect_identical(a$member_diagnostics$kish_validation_n_eff_adm, 2)
  expect_identical(a$diagnostics$admissibility_truncation_rate, 1)
  expect_identical(a$best_ks, 1L)
  set.seed(11)
  b <- svem_forward(y ~ 1, d, nBoot = 1)
  expect_identical(a$member_diagnostics$objective_value, b$member_diagnostics$objective_value)
  expect_equal(unname(a$parms[1]), unname(b$parms))

  d <- data.frame(x = -3:3, y = 3.5)
  f <- svem_forward(y ~ x, d, nBoot = 3)
  expect_equal(unname(f$parms), c(3.5, 0))
  expect_identical(f$member_diagnostics$objective_value, rep(log(.Machine$double.eps^2), 3))
  expect_equal(f$member_diagnostics$validation_sse, rep(0, 3))
  floor <- .Machine$double.eps^2
  expect_equal(.svem_forward_gaussian_score(c(0, floor / 4), c(floor / 2, floor)),
               log(floor) + c(0.5, 1))
  for (noise in c(0, 1e-14)) {
    d$y <- 2 * d$x + noise * (-1)^seq_len(nrow(d))
    exact <- svem_forward(y ~ x, d, nBoot = 1, weight_scheme = "Identity")
    small <- svem_forward(y ~ x, transform(d, y = y * 2^-400), nBoot = 1, weight_scheme = "Identity")
    expect_identical(exact$best_ks, 2L)
    expect_identical(exact$best_ks, small$best_ks)
    expect_equal(small$member_diagnostics$objective_value, exact$member_diagnostics$objective_value)
    expect_true(is.finite(exact$member_diagnostics$objective_value))
    expect_equal(as.numeric(predict(exact, d)), d$y, tolerance = 1e-12)
  }
  # Identical training means and variances at both admissible points.
  tied <- svem_forward(y ~ x, data.frame(x = c(1,-1,1,-1), y = c(1,1,-1,-1)),
                       nBoot = 1, weight_scheme = "Identity")
  expect_equal(tied$diagnostics$path_length_max, 2)
  expect_identical(tied$best_ks, 1L)
})


test_that("the best validation log score can occur after a worse intermediate point", {
  z <- rep(c(-1, 0, 1), 4)
  q <- seq_along(z) - 6.5
  d <- data.frame(x = q + z/3, w = q - z/3, y = z + rep(c(-.2, .1, .1, -.1), 3))
  set.seed(2)
  u <- runif(nrow(d))
  tr <- -log(u)
  va <- -log1p(-u)
  X <- model.matrix(y ~ x + w, d)
  candidates <- list(1L, 1:2, c(1L, 3L), 1:3)
  fits <- lapply(candidates, function(j) lm.wfit(X[, j, drop = FALSE], d$y, tr))
  t <- vapply(fits, function(f) weighted.mean(f$residuals^2, tr), numeric(1))
  v <- vapply(fits, function(f) weighted.mean(f$residuals^2, va), numeric(1))
  scores <- log(t) + v/t
  expect_lt(t[2], t[3])  # x enters first by training RSS.
  expect_gt(scores[2], scores[1])
  expect_lt(scores[4], scores[1])
  set.seed(2)
  fit <- svem_forward(y ~ x + w, d, nBoot = 1)
  expect_equal(fit$best_ks, 3L)
  expect_equal(fit$parms, fits[[4]]$coefficients)
  expect_equal(fit$member_diagnostics$objective_value + 2 * log(diff(range(d$y)) / 2),
               scores[4])
  expect_error(SVEMnet(y ~ x + w, d, objective = "gaussian_likelihood"), "arg")
})
