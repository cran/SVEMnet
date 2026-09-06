# Tests for forward_aicc() and svem_forward() (not run on CRAN)
skip_on_cran()
skip_if_not_installed("SVEMnet")

gen_forward_df <- function(n = 40L, sd = 0.25, with_factor = FALSE, seed = 7L) {
  set.seed(seed)
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  y  <- 1 + 2 * X1 - 1.5 * X2 + 1.2 * X1 * X2 + rnorm(n, 0, sd)
  if (with_factor) {
    F <- factor(sample(c("a", "b", "c"), n, TRUE))
    y <- y + 2 * (F == "b") - 1.5 * (F == "c")
    data.frame(y, X1, X2, X3, F)
  } else {
    data.frame(y, X1, X2, X3)
  }
}

## ---------------------------------------------------------------------------
## forward_aicc(): deterministic forward selection
## ---------------------------------------------------------------------------

test_that("forward_aicc recovers strong true terms and matches lm refit", {
  dat <- gen_forward_df()
  fit <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat)

  expect_s3_class(fit, "svem_forward_ic")
  expect_s3_class(fit, "svem_model")
  expect_true(all(c("X1", "X2", "X1:X2") %in% fit$selected_terms))

  # Coefficients equal an OLS refit on the selected design columns,
  # compared column-by-column (order-preserving, catches misassignment)
  sel_cols <- names(fit$parms)[fit$parms != 0]
  sel_cols <- setdiff(sel_cols, "(Intercept)")
  Xsel <- fit$training_X[, sel_cols, drop = FALSE]
  lm_fit <- lm(dat$y ~ Xsel)
  expect_equal(unname(fit$parms["(Intercept)"]),
               unname(coef(lm_fit)[1]), tolerance = 1e-8)
  expect_equal(unname(fit$parms[sel_cols]),
               unname(coef(lm_fit)[-1]), tolerance = 1e-8)

  # Predictions from predict() match the stored coefficients
  preds <- predict(fit, dat)
  expect_equal(preds, unname(fit$y_pred), tolerance = 1e-10)
})

test_that("forward_aicc selection path improves the criterion monotonically", {
  dat <- gen_forward_df()
  fit <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat)
  path <- fit$selection_path
  expect_true(nrow(path) >= 2L)
  expect_true(all(diff(path$criterion) < 0))
  expect_true(all(path$improvement > 0))
  expect_true(all(diff(path$k) > 0))
  # Final refit criterion agrees with the last accepted step
  expect_equal(fit$criterion_value, path$criterion[nrow(path)],
               tolerance = 1e-8)
})

test_that("forward_aicc criterion values match the closed-form formulas", {
  dat <- gen_forward_df()
  n <- nrow(dat)
  for (crit in c("AICc", "AIC", "BIC")) {
    fit <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat, criterion = crit)
    sel_cols <- setdiff(names(fit$parms)[fit$parms != 0], "(Intercept)")
    Xsel <- fit$training_X[, sel_cols, drop = FALSE]
    rss <- sum(residuals(lm(dat$y ~ Xsel))^2)
    k <- length(sel_cols) + 1L  # fitted coefficients including the intercept
    K <- k + 1L                 # + residual variance
    expected <- switch(crit,
      AICc = n * log(rss / n) + 2 * K + 2 * K * (K + 1) / (n - K - 1),
      AIC  = n * log(rss / n) + 2 * K,
      BIC  = n * log(rss / n) + K * log(n)
    )
    expect_equal(fit$criterion_value, expected, tolerance = 1e-8,
                 label = paste0(crit, " criterion_value"))
  }
})

test_that("forward_aicc supports AIC and BIC criteria", {
  dat <- gen_forward_df()
  fit_aic <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat, criterion = "AIC")
  fit_bic <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat, criterion = "BIC")
  expect_identical(fit_aic$criterion, "AIC")
  expect_identical(fit_bic$criterion, "BIC")
  expect_true(all(c("X1", "X2") %in% fit_aic$selected_terms))
  expect_true(all(c("X1", "X2") %in% fit_bic$selected_terms))
  # BIC penalizes size at least as hard as AIC for n >= 8
  expect_lte(length(fit_bic$selected_terms), length(fit_aic$selected_terms))
})

test_that("forward_aicc is deterministic", {
  dat <- gen_forward_df()
  f1 <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat)
  f2 <- forward_aicc(y ~ (X1 + X2 + X3)^2, dat)
  expect_identical(f1$parms, f2$parms)
  expect_identical(f1$selection_path, f2$selection_path)
})

test_that("forward_aicc works with a bigexp_spec and response swap", {
  dat <- gen_forward_df()
  dat$y2 <- 0.5 - 2 * dat$X2 + rnorm(nrow(dat), 0, 0.2)
  spec <- bigexp_terms(y ~ X1 + X2 + X3, dat,
                       factorial_order = 2, polynomial_order = 2)

  fit <- forward_aicc(spec, dat)
  expect_true(fit$used_bigexp_spec)
  expect_true(all(c("X1", "X2") %in% fit$selected_terms))
  expect_length(predict(fit, dat), nrow(dat))

  fit2 <- forward_aicc(spec, dat, response = "y2")
  expect_true("X2" %in% fit2$selected_terms)
  # X2 slope should be near -2 on the swapped response
  expect_equal(unname(fit2$parms["X2"]), -2, tolerance = 0.2)
})

test_that("forward_aicc keeps factor contrast columns together", {
  dat <- gen_forward_df(n = 60L, with_factor = TRUE)
  fit <- forward_aicc(y ~ X1 + X2 + X3 + F, dat)
  expect_true("F" %in% fit$selected_terms)
  f_cols <- grep("^F", names(fit$parms), value = TRUE)
  expect_length(f_cols, 2L)
  # Both contrast columns enter together and are estimated (nonzero)
  expect_true(all(fit$parms[f_cols] != 0))
  # The factor enters as one grouped step: model size jumps by 2 columns
  idx <- which(fit$selection_path$term == "F")
  expect_length(idx, 1L)
  k_prev <- if (idx == 1L) 1L else fit$selection_path$k[idx - 1L]
  expect_identical(fit$selection_path$k[idx] - k_prev, 2L)
})

test_that("forward_aicc rejects unsupported inputs", {
  dat <- gen_forward_df()
  expect_error(forward_aicc(y ~ 0 + X1, dat), "intercept convention")
  datf <- dat; datf$yf <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  expect_error(forward_aicc(yf ~ X1, datf), "numeric response")
})

test_that("forward_aicc returns intercept-only with a warning when n is too small for AICc", {
  dat <- data.frame(y = c(1, 2, 3), X1 = c(0.1, 0.5, 0.9))
  expect_warning(fit <- forward_aicc(y ~ X1, dat), "undefined")
  expect_length(fit$selected_terms, 0L)
  expect_identical(unname(fit$parms["X1"]), 0)
})

test_that("forward_aicc handles a constant response as intercept-only", {
  dat <- data.frame(y = rep(2, 20), X1 = rnorm(20), X2 = rnorm(20))
  expect_warning(fit <- forward_aicc(y ~ X1 + X2, dat), "noise floor")
  expect_length(fit$selected_terms, 0L)
  expect_equal(unname(fit$parms["(Intercept)"]), 2, tolerance = 1e-10)
})

test_that("forward selection is not truncated by large response means", {
  set.seed(12)
  n <- 40
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1e8 + 2 * x1 + rnorm(n, 0, 0.1), x1, x2)
  fit <- forward_aicc(y ~ x1 + x2, dat)
  expect_true("x1" %in% fit$selected_terms)
  expect_equal(unname(fit$parms[["x1"]]), 2, tolerance = 0.1)

  # Mid-path variant: a small later term must still enter
  dat2 <- data.frame(y = 1e5 + 2 * x1 + 0.03 * x2 + rnorm(n, 0, 1e-4), x1, x2)
  fit2 <- forward_aicc(y ~ x1 + x2, dat2)
  expect_setequal(fit2$selected_terms, c("x1", "x2"))
  expect_equal(unname(fit2$parms[["x2"]]), 0.03, tolerance = 1e-3)

  set.seed(13)
  sf <- svem_forward(y ~ x1 + x2, dat, nBoot = 5)
  expect_gt(sf$diagnostics$selection_frequencies[["x1"]], 0.9)
})

## ---------------------------------------------------------------------------
## svem_forward(): SVEM ensemble with forward-selection base learners
## ---------------------------------------------------------------------------

test_that("svem_forward fits, predicts, and interoperates with svem_model methods", {
  dat <- gen_forward_df()
  set.seed(11)
  fit <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 15)

  expect_s3_class(fit, "svem_forward")
  expect_s3_class(fit, "svem_model")
  expect_identical(fit$family, "gaussian")
  expect_identical(fit$objective_used, "wAIC")
  expect_equal(nrow(fit$coef_matrix), 15L)
  expect_identical(colnames(fit$coef_matrix)[1], "(Intercept)")

  preds <- predict(fit, dat)
  expect_length(preds, nrow(dat))
  expect_true(all(is.finite(preds)))

  out <- predict(fit, dat, se.fit = TRUE, interval = TRUE, level = 0.9)
  expect_named(out, c("fit", "se.fit", "lwr", "upr"))
  expect_true(all(out$lwr <= out$upr))

  # debias machinery mirrors SVEMnet (nBoot >= 10)
  expect_s3_class(fit$debias_fit, "lm")
  preds_db <- predict(fit, dat, debias = TRUE)
  expect_length(preds_db, nrow(dat))

  cc <- coef(fit)
  expect_identical(cc, fit$parms)

  nz <- svem_nonzero(fit, plot = FALSE, print_table = FALSE)
  expect_true(is.data.frame(nz))

  # Strong terms should be selected in most bootstraps
  freq <- fit$diagnostics$selection_frequencies
  expect_true(freq[["X1"]] > 0.8)
  expect_true(freq[["X2"]] > 0.8)
})

test_that("svem_forward is reproducible under set.seed", {
  dat <- gen_forward_df()
  set.seed(5)
  f1 <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 6)
  set.seed(5)
  f2 <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 6)
  expect_identical(f1$parms, f2$parms)
  expect_identical(f1$coef_matrix, f2$coef_matrix)
})

test_that("svem_forward supports wBIC and wSSE objectives", {
  dat <- gen_forward_df()
  set.seed(3)
  fit_bic <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 6,
                          objective = "wBIC")
  expect_identical(fit_bic$objective_used, "wBIC")
  expect_true(all(is.finite(fit_bic$parms)))

  set.seed(3)
  fit_sse <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 6,
                          objective = "wSSE")
  expect_identical(fit_sse$objective_used, "wSSE")
  # wSSE paths are never truncated by the admissibility ceiling
  expect_identical(fit_sse$diagnostics$admissibility_truncation_rate, 0)
})

test_that("svem_forward Identity weights with nBoot = 1 is deterministic", {
  dat <- gen_forward_df()
  f1 <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 1,
                     weight_scheme = "Identity")
  f2 <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 1,
                     weight_scheme = "Identity")
  expect_identical(f1$parms, f2$parms)
  expect_true(all(c("X1", "X2") %in%
                    names(which(f1$diagnostics$selection_frequencies == 1))))
})

test_that("svem_forward FRW_plain runs", {
  dat <- gen_forward_df()
  set.seed(4)
  fit <- svem_forward(y ~ X1 + X2 + X3, dat, nBoot = 5,
                      weight_scheme = "FRW_plain")
  expect_equal(nrow(fit$coef_matrix), 5L)
  expect_true(all(is.finite(fit$parms)))
})

test_that("svem_forward works with a bigexp_spec, factors, and response swap", {
  dat <- gen_forward_df(n = 60L, with_factor = TRUE)
  dat$y2 <- 1 + 3 * dat$X1 + rnorm(nrow(dat), 0, 0.2)
  spec <- bigexp_terms(y ~ X1 + X2 + X3 + F, dat, factorial_order = 2)

  set.seed(21)
  fit <- svem_forward(spec, dat, nBoot = 8)
  expect_true(fit$used_bigexp_spec)
  expect_length(predict(fit, dat), nrow(dat))

  # Factor contrast columns always enter together within each bootstrap:
  # each row has either all F main-effect columns zero or all nonzero
  f_cols <- grep("^F", colnames(fit$coef_matrix), value = TRUE)
  f_main <- f_cols[!grepl(":", f_cols)]
  expect_length(f_main, 2L)
  cm <- fit$coef_matrix[, f_main, drop = FALSE]
  n_nonzero <- rowSums(cm != 0)
  expect_true(all(n_nonzero %in% c(0L, length(f_main))))

  set.seed(22)
  fit2 <- svem_forward(spec, dat, response = "y2", nBoot = 8)
  expect_equal(unname(fit2$parms["X1"]), 3, tolerance = 0.5)
})

test_that("svem_forward skips exactly collinear candidates instead of failing", {
  set.seed(9)
  n <- 40
  A <- runif(n, 0.1, 0.6)
  B <- runif(n, 0.1, 0.35)
  C <- 1 - A - B  # exact mixture collinearity with the intercept
  y <- 2 + 3 * A - 2 * B + rnorm(n, 0, 0.2)
  dat <- data.frame(y, A, B, C)

  set.seed(10)
  fit <- svem_forward(y ~ A + B + C, dat, nBoot = 6)
  expect_true(all(is.finite(fit$parms)))
  # At most two of the three mixture components can be active per bootstrap
  active <- rowSums(fit$coef_matrix[, c("A", "B", "C"), drop = FALSE] != 0)
  expect_true(all(active <= 2))
})

test_that("svem_forward handles a constant response as intercept-only", {
  dat <- data.frame(y = rep(5, 25), X1 = rnorm(25), X2 = rnorm(25))
  set.seed(2)
  fit <- svem_forward(y ~ X1 + X2, dat, nBoot = 4)
  expect_true(all(fit$coef_matrix[, c("X1", "X2")] == 0))
  expect_true(all(abs(fit$coef_matrix[, "(Intercept)"] - 5) < 1e-8))
})

test_that("svem_forward predictions flag unseen factor levels as NA", {
  dat <- gen_forward_df(n = 60L, with_factor = TRUE)
  set.seed(13)
  fit <- svem_forward(y ~ X1 + F, dat, nBoot = 5)
  nd <- dat[1:3, c("X1", "F")]
  nd$F <- factor(c("a", "b", "zz"), levels = c("a", "b", "zz"))
  expect_warning(preds <- predict(fit, nd), "unseen")
  expect_true(is.na(preds[3]))
  expect_true(all(is.finite(preds[1:2])))
})

test_that("svem_forward validates inputs like SVEMnet", {
  dat <- gen_forward_df()
  expect_error(svem_forward(y ~ 0 + X1, dat), "intercept convention")
  expect_error(svem_forward(y ~ X1, dat, nBoot = 0), "nBoot")
  datf <- dat; datf$yf <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  expect_error(svem_forward(yf ~ X1, datf), "numeric response")
})

test_that("svem_forward matches an independent single-bootstrap oracle", {
  dat <- gen_forward_df()
  n <- nrow(dat)
  set.seed(123)
  fit <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 1)

  # Reproduce the replicate independently: same RNG stream, lm.wfit-based
  # greedy path on training weights, wAIC scoring with validation weights
  # applied to residuals on the ORIGINAL rows.
  eps <- .Machine$double.eps
  set.seed(123)
  U <- pmin(pmax(runif(n), eps), 1 - eps)
  w_train <- -log(U); w_valid <- -log1p(-U)
  w_train <- w_train * (n / sum(w_train))
  w_valid <- w_valid * (n / sum(w_valid))
  n_eff_raw <- sum(w_valid)^2 / (sum(w_valid^2) + eps)
  n_eff_adm <- max(2, min(n, n_eff_raw))

  Xf <- model.matrix(y ~ (X1 + X2 + X3)^2, dat)
  y  <- dat$y
  cols <- 1L
  remaining <- 2:ncol(Xf)
  points <- list(cols)
  while (length(remaining) && (length(cols) - 1L) < n_eff_adm) {
    wrss <- vapply(remaining, function(cc) {
      f <- lm.wfit(Xf[, c(cols, cc), drop = FALSE], y, w_train)
      sum(w_train * f$residuals^2)
    }, numeric(1))
    best <- remaining[which.min(wrss)]
    cols <- c(cols, best)
    remaining <- setdiff(remaining, best)
    points[[length(points) + 1L]] <- cols
  }
  scores <- vapply(points, function(idx) {
    f  <- lm.wfit(Xf[, idx, drop = FALSE], y, w_train)
    r  <- as.vector(Xf[, idx, drop = FALSE] %*% f$coefficients) - y
    ss <- max(sum(w_valid * r^2), eps)
    k  <- length(idx)
    if ((k - 1L) >= n_eff_adm) return(Inf)
    sum(w_valid) * log(ss / sum(w_valid)) + 2 * k
  }, numeric(1))
  winner <- which.min(scores)
  beta <- setNames(numeric(ncol(Xf)), colnames(Xf))
  f_win <- lm.wfit(Xf[, points[[winner]], drop = FALSE], y, w_train)
  beta[points[[winner]]] <- f_win$coefficients

  expect_equal(unname(fit$coef_matrix[1, ]), unname(beta), tolerance = 1e-7)
})

test_that("admissibility ceiling truncates wAIC paths but never wSSE paths", {
  set.seed(99)
  n <- 15
  Xs <- matrix(rnorm(n * 6), n, 6)
  colnames(Xs) <- paste0("X", 1:6)
  dat <- data.frame(y = 1 + 2 * Xs[, 1] - 1.5 * Xs[, 2] + rnorm(n, 0, 0.3), Xs)
  fml <- y ~ (X1 + X2 + X3 + X4 + X5 + X6)^2  # 21 slope terms >> n_eff_adm

  set.seed(3)
  f_aic <- svem_forward(fml, dat, nBoot = 6, objective = "wAIC")
  set.seed(3)
  f_sse <- svem_forward(fml, dat, nBoot = 6, objective = "wSSE")

  expect_gt(f_aic$diagnostics$admissibility_truncation_rate, 0)
  expect_identical(f_sse$diagnostics$admissibility_truncation_rate, 0)
  expect_gte(f_sse$diagnostics$path_length_max,
             f_aic$diagnostics$path_length_max)
})

test_that("badly scaled predictors survive acceptance and refit consistently", {
  # A timestamp-scale column has projected-norm ratio ~5e-8 against the
  # intercept: above the acceptance gate (1e-8) but below base qr()'s default
  # rank tolerance (1e-7). Regression test: the refits must honor the gate.
  set.seed(8)
  n <- 40
  ts <- 1.7e9 + seq(0, 300, length.out = n)
  x2 <- rnorm(n)
  y  <- 3 + 0.05 * (ts - mean(ts)) + 2 * x2 + rnorm(n, 0, 0.2)
  dat <- data.frame(y, ts, x2)

  fit <- forward_aicc(y ~ ts + x2, dat)
  expect_true("ts" %in% fit$selected_terms)
  expect_true(fit$parms[["ts"]] != 0)
  # The selected model must beat the intercept-only AICc it started from
  rss0  <- sum((y - mean(y))^2)
  aicc0 <- n * log(rss0 / n) + 2 * 2 + 2 * 2 * 3 / (n - 3)
  expect_lt(fit$criterion_value, aicc0)
  expect_gt(cor(fit$y_pred, y)^2, 0.9)

  set.seed(9)
  sfit <- svem_forward(y ~ ts + x2, dat, nBoot = 6)
  expect_gt(sfit$diagnostics$selection_frequencies[["ts"]], 0.5)
  expect_true(sfit$parms[["ts"]] != 0)
})

test_that("small independent directions are checked by exact refitting", {
  # Variation is meaningful but its norm relative to the timestamp origin
  # is below the projection shortcut's 1e-8 gate. The exact QR solver still
  # resolves it; silently skipping the term loses almost all of the signal.
  set.seed(81)
  n <- 40L
  ts <- 1.7e9 + seq(0, 30, length.out = n)
  x <- ts - mean(ts)
  y <- 3 + 0.5 * x + rnorm(n, 0, 0.1)
  dat <- data.frame(y, ts, x)

  fit <- forward_aicc(y ~ ts, dat)
  reference <- stats::lm(y ~ x, dat)
  expect_identical(fit$selected_terms, "ts")
  expect_equal(unname(fit$y_pred), unname(stats::fitted(reference)),
               tolerance = 1e-5)
  expect_lt(fit$selection_path$rss[1], 1)

  # The same fallback applies independently in each weighted SVEM path;
  # centering an ordinary linear predictor provides a stable oracle without
  # changing its model space, weights, objective, or selected model size.
  set.seed(82)
  sfit <- svem_forward(y ~ ts, dat, nBoot = 6)
  set.seed(82)
  sref <- svem_forward(y ~ x, dat, nBoot = 6)
  expect_equal(sfit$best_ks, sref$best_ks)
  expect_equal(sfit$diagnostics$selection_frequencies[["ts"]], 1)
  expect_equal(unname(sfit$y_pred), unname(sref$y_pred), tolerance = 1e-5)
})

test_that("ambiguous multi-column terms retain whole-group rank checks", {
  set.seed(83)
  n <- 50L
  ts <- 1.7e9 + seq(0, 30, length.out = n)
  x <- ts - mean(ts)
  z <- rnorm(n)
  y <- 3 + 0.5 * x + 2 * z + rnorm(n, 0, 0.1)
  dat <- data.frame(y, ts, x, z)

  fit <- forward_aicc(y ~ cbind(ts, z), dat)
  reference <- stats::lm(y ~ x + z, dat)
  expect_identical(fit$selected_terms, "cbind(ts, z)")
  expect_equal(unname(fit$y_pred), unname(stats::fitted(reference)),
               tolerance = 1e-5)
  expect_equal(fit$selection_path$k, 3L)

  # An actually rank-deficient block must still be skipped as a whole,
  # even when one of its directions alone would predict the response.
  bad <- forward_aicc(y ~ cbind(ts, ts), dat)
  expect_length(bad$selected_terms, 0L)
  expect_equal(unname(bad$y_pred), rep(mean(y), n), tolerance = 1e-10)
})

test_that("near-perfect candidate fits are compared by exact RSS", {
  # The projected RSS update is cancellation noise for near-perfect fits;
  # the exact-refit guard must pick the candidate with the smaller true RSS.
  set.seed(1)
  n <- 20
  cB <- rnorm(n)
  y  <- cB + 4e-9 * rnorm(n)
  c2 <- y + 1e-9 * rnorm(n)
  dat <- data.frame(y, cB, c2)

  # After c2 enters, the fit is at the noise floor, so selection stops with
  # a warning while cB remains unevaluated.
  expect_warning(fit <- forward_aicc(y ~ cB + c2, dat), "noise floor")
  expect_identical(fit$selected_terms[1], "c2")
  expect_true(all(fit$selection_path$rss > 0))
})

test_that("svem_forward interoperates with svem_random_table_multi", {
  dat <- gen_forward_df()
  set.seed(31)
  fit <- svem_forward(y ~ (X1 + X2 + X3)^2, dat, nBoot = 6)
  tab <- svem_random_table_multi(list(fit), n = 25)
  expect_true(is.list(tab))
  expect_true(is.data.frame(tab$all))
  expect_equal(nrow(tab$all), 25L)
  expect_true("y_pred" %in% colnames(tab$all))
})
