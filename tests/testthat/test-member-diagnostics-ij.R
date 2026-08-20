skip_on_cran()

.diagnostic_data <- function(n = 24L, seed = 441) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- 2 + 1.4 * x1 - 0.8 * x2 + 0.7 * x1 * x2 + rnorm(n, sd = 0.4)
  data.frame(y, x1, x2)
}

test_that("member diagnostics are additive and do not alter the fit", {
  dat <- .diagnostic_data()
  set.seed(91)
  base <- SVEMnet(y ~ x1 * x2, dat, nBoot = 12,
                  glmnet_alpha = c(0.5, 1), relaxed = TRUE,
                  objective = "wBIC")
  set.seed(91)
  kept <- SVEMnet(y ~ x1 * x2, dat, nBoot = 12,
                  glmnet_alpha = c(0.5, 1), relaxed = TRUE,
                  objective = "wBIC",
                  store_member_weights = TRUE)

  expect_identical(base$coef_matrix, kept$coef_matrix)
  expect_identical(base$parms, kept$parms)
  expect_null(base$member_weights)
  expect_true(is.list(kept$member_weights))
  expect_equal(dim(kept$member_weights$uniforms), c(12, nrow(dat)))
  expect_equal(rowMeans(kept$member_weights$training), rep(1, 12), tolerance = 1e-12)
  expect_equal(rowMeans(kept$member_weights$validation), rep(1, 12), tolerance = 1e-12)
  expect_equal(
    kept$member_weights$training,
    -log(kept$member_weights$uniforms) /
      rowMeans(-log(kept$member_weights$uniforms)),
    tolerance = 1e-12
  )
  expect_equal(
    kept$member_weights$validation,
    -log1p(-kept$member_weights$uniforms) /
      rowMeans(-log1p(-kept$member_weights$uniforms)),
    tolerance = 1e-12
  )

  md <- kept$member_diagnostics
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), kept$nBoot)
  expect_true(all(c(
    "objective", "validation_sse", "objective_value", "support_size", "effective_df",
    "kish_validation_n_eff", "residual_scale_support",
    "admissibility_warning"
  ) %in% names(md)))
  expect_true(all(is.finite(md$validation_sse)))
  expect_true(all(md$objective == "wBIC"))
  expect_true(all(is.finite(md$effective_df)))
  expect_equal(md$effective_df, kept$best_edfs)
  expect_equal(md$residual_scale_support^2,
               md$validation_sse / md$residual_df_support,
               tolerance = 1e-12)
})

test_that("experimental IJ variance matches the stored-weight formula", {
  dat <- .diagnostic_data(seed = 442)
  set.seed(92)
  fit <- SVEMnet(y ~ x1 * x2, dat, nBoot = 16,
                 glmnet_alpha = 1, relaxed = FALSE,
                 store_member_weights = TRUE)
  nd <- dat[1:3, c("x1", "x2")]
  got <- svem_ij_variance(fit, nd)

  mm <- model.matrix(~ x1 * x2, nd)
  Tmat <- mm %*% t(fit$coef_matrix)
  Tc <- Tmat - rowMeans(Tmat)
  W <- fit$member_weights$training
  Wc <- sweep(W, 2, colMeans(W), "-")
  raw <- rowSums((Tc %*% Wc / fit$nBoot)^2)
  correction <- nrow(dat) / fit$nBoot^2 * rowSums(Tc^2)

  expect_equal(got$fit, unname(rowMeans(Tmat)), tolerance = 1e-12)
  expect_equal(got$variance_raw, unname(raw), tolerance = 1e-12)
  expect_equal(got$finite_B_correction, unname(correction), tolerance = 1e-12)
  expect_equal(got$variance, unname(pmax(raw - correction, 0)), tolerance = 1e-12)
  expect_false(any(got$calibrated))
  expect_match(attr(got, "caution"), "Not calibrated")

  uncorrected <- svem_ij_variance(fit, nd, finite_B = FALSE)
  expect_equal(uncorrected$variance, unname(raw), tolerance = 1e-12)
  expect_equal(uncorrected$finite_B_correction, rep(0, 3))
})

test_that("IJ requires stored random weights and Gaussian family", {
  dat <- .diagnostic_data(seed = 443)
  set.seed(93)
  no_weights <- SVEMnet(y ~ x1 + x2, dat, nBoot = 4)
  expect_error(svem_ij_variance(no_weights, dat[1, ]), "Stored member weights")

  set.seed(94)
  ident <- SVEMnet(y ~ x1 + x2, dat, nBoot = 4,
                   weight_scheme = "Identity", store_member_weights = TRUE)
  expect_error(svem_ij_variance(ident, dat[1, ]), "Identity")

  dat$binary <- as.integer(dat$x1 + 0.25 * dat$x2 > 0)
  set.seed(96)
  binomial_fit <- suppressWarnings(SVEMnet(
    binary ~ x1 + x2, dat, nBoot = 4, family = "binomial",
    relaxed = FALSE, store_member_weights = TRUE
  ))
  expect_error(
    svem_ij_variance(binomial_fit, dat[1, c("x1", "x2")]),
    "Gaussian"
  )
})

test_that("IJ rebuilds Gaussian factor model matrices", {
  dat <- .diagnostic_data(seed = 446)
  dat$group <- factor(rep(c("low", "middle", "high"), length.out = nrow(dat)))
  dat$y <- dat$y + c(low = -0.2, middle = 0, high = 0.3)[dat$group]
  set.seed(97)
  fit <- SVEMnet(y ~ x1 + group, dat, nBoot = 8,
                 store_member_weights = TRUE)
  got <- svem_ij_variance(fit, dat[1:3, c("x1", "group")])
  expect_equal(nrow(got), 3L)
  expect_true(all(is.finite(got$fit)))
  expect_false(any(got$calibrated))
})

test_that("forward members expose the same diagnostic contract", {
  dat <- .diagnostic_data(seed = 444)
  set.seed(95)
  fit <- svem_forward(y ~ x1 * x2, dat, nBoot = 10,
                      store_member_weights = TRUE)
  expect_equal(nrow(fit$member_diagnostics), fit$nBoot)
  expect_equal(fit$best_edfs, as.numeric(fit$best_ks))
  expect_true(all(fit$member_diagnostics$objective == fit$objective_used))
  expect_equal(dim(fit$member_weights$training), c(fit$nBoot, nrow(dat)))
  expect_true(all(is.finite(fit$member_diagnostics$validation_sse)))
  ij <- svem_ij_variance(fit, dat[1:2, c("x1", "x2")])
  expect_equal(nrow(ij), 2)
  expect_false(any(ij$calibrated))
})

test_that("store_member_weights validates a logical scalar", {
  dat <- .diagnostic_data(seed = 445)
  expect_error(
    SVEMnet(y ~ x1 + x2, dat, nBoot = 2, store_member_weights = 1),
    "TRUE or FALSE"
  )
  expect_error(
    svem_forward(y ~ x1 + x2, dat, nBoot = 2, store_member_weights = NA),
    "TRUE or FALSE"
  )
})
