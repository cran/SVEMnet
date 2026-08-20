# Regression tests for the 3.2.3 maintenance fixes.
skip_on_cran()
skip_if_not_installed("SVEMnet")

test_that("glmnet_with_cv 1se rule selects a lambda at least as large as the min rule", {
  set.seed(123)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1.5, -1, 0.5, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n))
  df <- data.frame(y = y, X)
  colnames(df) <- c("y", paste0("x", seq_len(p)))

  fit_min <- glmnet_with_cv(y ~ ., df, glmnet_alpha = 1, nfolds = 5,
                            repeats = 2, choose_rule = "min", seed = 11)
  fit_1se <- glmnet_with_cv(y ~ ., df, glmnet_alpha = 1, nfolds = 5,
                            repeats = 2, choose_rule = "1se", seed = 11)

  # 1se must be at least as regularized as min (lambda_1se >= lambda_min)
  expect_gte(fit_1se$best_lambda, fit_min$best_lambda)

  # index-level property inside the aggregated CV curves: glmnet lambda paths
  # are decreasing, so the 1se index must not exceed the min index
  for (df_cv in fit_1se$cv_summary) {
    expect_lte(df_cv$idx_1se[1L], df_cv$idx_min[1L])
    expect_gte(df_cv$lambda[df_cv$idx_1se[1L]], df_cv$lambda[df_cv$idx_min[1L]])
  }

  # 1se should never select a denser model than min
  k_min <- fit_min$diagnostics$k_final_no_intercept
  k_1se <- fit_1se$diagnostics$k_final_no_intercept
  expect_lte(k_1se, k_min)
})

test_that("glmnet_with_cv relaxed fit jointly selects lambda and gamma", {
  set.seed(456)
  n <- 90; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, -1.5, rep(0, p - 2))
  y <- as.numeric(X %*% beta + rnorm(n))
  df <- data.frame(y = y, X)
  colnames(df) <- c("y", paste0("x", seq_len(p)))

  fit_rel <- glmnet_with_cv(y ~ ., df, glmnet_alpha = 1, nfolds = 5,
                            repeats = 1, choose_rule = "min",
                            relaxed = TRUE, seed = 21)
  expect_true(is.finite(fit_rel$best_gamma))
  expect_true(fit_rel$best_gamma %in% c(0, 0.25, 0.5, 0.75, 1))
  expect_true(any(grepl("jointly selected", fit_rel$note, fixed = TRUE)))
})

test_that("single-level discrete numeric supports round-trip unchanged", {
  set.seed(31)
  n <- 40
  dat <- data.frame(
    y  = rnorm(n),
    X1 = runif(n),
    D  = rep(5, n)
  )

  spec <- bigexp_terms(
    y ~ X1 + D,
    data             = dat,
    factorial_order  = 1,
    polynomial_order = 1,
    discrete_numeric = list(D = 5),
    report           = FALSE
  )
  fit <- SVEMnet(spec, dat, nBoot = 15)

  tab <- svem_random_table_multi(fit, n = 200)
  expect_true(all(tab$data$D == 5))
})

test_that("bigexp_terms honors dot exclusions like y ~ . - X2", {
  set.seed(32)
  dat <- data.frame(
    y  = rnorm(20),
    X1 = rnorm(20),
    X2 = rnorm(20),
    X3 = rnorm(20)
  )
  spec <- bigexp_terms(y ~ . - X2, data = dat,
                       factorial_order = 1, polynomial_order = 1,
                       report = FALSE)
  expect_false("X2" %in% spec$vars)
  expect_true(all(c("X1", "X3") %in% spec$vars))
})

test_that("bigexp_terms rejects transformed RHS terms with a clear error", {
  dat <- data.frame(y = rnorm(20), X1 = runif(20, 1, 2), X2 = rnorm(20))
  expect_error(
    bigexp_terms(y ~ log(X1) + X2, data = dat, report = FALSE),
    "not found as columns"
  )
})

test_that("SVEMnet routes direct maxit but warns on reserved arguments", {
  set.seed(33)
  n <- 30
  dat <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))

  # maxit is routed via the version-compat layer without a warning
  expect_warning(
    SVEMnet(y ~ x1 + x2, dat, nBoot = 5, maxit = 2e5),
    regexp = NA
  )

  # reserved arguments are dropped with a warning
  expect_warning(
    SVEMnet(y ~ x1 + x2, dat, nBoot = 5, nlambda = 100),
    "reserved"
  )
})

test_that("svem_random_table_multi errors instead of inventing factor levels", {
  set.seed(34)
  n <- 40
  dat <- data.frame(
    y  = rnorm(n),
    X1 = runif(n),
    G  = factor(sample(c("A", "B"), n, TRUE))
  )
  fit <- SVEMnet(y ~ X1 + G, dat, nBoot = 10)

  # corrupt the recorded levels to simulate a legacy/broken object
  fit$sampling_schema$factor_levels$G <- character(0)
  fit$xlevels$G <- NULL

  expect_error(svem_random_table_multi(fit, n = 50), "No recorded levels")
})

test_that("SVEMnet and glmnet_with_cv ignore user offset/weights with a warning", {
  set.seed(41)
  n <- 40
  dat <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  off <- runif(n, 0, 3)
  w   <- runif(n, 0.5, 2)

  # SVEMnet: offset ignored -> identical fit to no-offset call with same seed
  set.seed(42); fit_off <- suppressWarnings(SVEMnet(y ~ x1 + x2, dat, nBoot = 15, offset = off))
  set.seed(42); fit_ref <- SVEMnet(y ~ x1 + x2, dat, nBoot = 15)
  expect_warning(SVEMnet(y ~ x1 + x2, dat, nBoot = 3, offset = off), "offset")
  expect_equal(fit_off$parms, fit_ref$parms)

  # glmnet_with_cv: weights and offset ignored -> identical to plain call
  expect_warning(
    fw <- glmnet_with_cv(y ~ x1 + x2, dat, glmnet_alpha = 1, nfolds = 5,
                         repeats = 1, seed = 7, weights = w),
    "weights"
  )
  f0 <- glmnet_with_cv(y ~ x1 + x2, dat, glmnet_alpha = 1, nfolds = 5,
                       repeats = 1, seed = 7)
  expect_equal(fw$parms, f0$parms)
  expect_warning(
    glmnet_with_cv(y ~ x1 + x2, dat, glmnet_alpha = 1, nfolds = 5,
                   repeats = 1, seed = 7, offset = off),
    "offset"
  )
})

test_that("misspelled glmnet arguments are caught with a warning", {
  set.seed(43)
  dat <- data.frame(y = rnorm(30), x1 = rnorm(30), x2 = rnorm(30))
  expect_warning(
    SVEMnet(y ~ x1 + x2, dat, nBoot = 3, bogus_arg = 5),
    "unrecognized"
  )
  expect_warning(
    glmnet_with_cv(y ~ x1 + x2, dat, glmnet_alpha = 1, nfolds = 5,
                   repeats = 1, seed = 3, bogus_arg = 5),
    "unrecognized"
  )
})

test_that("significance test validates mixture variables", {
  set.seed(44)
  d <- data.frame(A = runif(30, 0.1, 0.6), B = runif(30, 0.1, 0.6))
  d$C <- pmax(0.05, 1 - d$A - d$B)
  d$y <- 2 + 3 * d$A + 1.5 * d$B + 1.2 * d$C + rnorm(30, 0, 0.3)

  expect_error(
    svem_significance_test_parallel(
      y ~ A + B + C, d,
      mixture_groups = list(list(vars = c("A", "B", "Czz"),
                                 lower = c(0.05, 0.05, 0.05),
                                 upper = c(0.9, 0.9, 0.9), total = 1)),
      nPoint = 100, nSVEM = 1, nPerm = 20, nBoot = 5,
      nCore = 1, verbose = FALSE
    ),
    "Mixture variables not in model predictors"
  )
})

test_that("bigexp_terms rejects single-level categorical predictors clearly", {
  dat <- data.frame(y = rnorm(20), x1 = rnorm(20), G = factor(rep("A", 20)))
  expect_error(
    bigexp_terms(y ~ x1 + G, data = dat, report = FALSE),
    "fewer than 2 levels"
  )
})

test_that("selection label round-trips into svem_export_candidates_csv", {
  set.seed(45)
  st <- data.frame(x1 = runif(50), x2 = runif(50), score = rnorm(50))
  my_label <- "round7_explore"
  sel <- svem_select_from_score_table(
    score_table = st, target = "score", direction = "max",
    k = 2, top_type = "n", top = 10, label = my_label
  )
  expect_identical(sel$label, my_label)
  out <- svem_export_candidates_csv(sel, write_file = FALSE)
  expect_true(all(out$selection_label == my_label))
})

test_that("uncertainty_measure anchors are shared between score_table and original data", {
  set.seed(35)
  n <- 40
  dat <- data.frame(y = rnorm(n), x1 = runif(n), x2 = runif(n))
  fit <- SVEMnet(y ~ x1 + x2, dat, nBoot = 20)

  scored <- svem_score_random(
    objects = list(y = fit),
    goals   = list(y = list(goal = "max", weight = 1)),
    data    = dat,
    n       = 500,
    verbose = FALSE
  )

  st <- scored$score_table
  od <- scored$original_data_scored
  expect_true(is.data.frame(od))

  # An original-data row and a sampled row with (nearly) the same raw CI width
  # must have (nearly) the same weighted normalized width, because the anchors
  # are shared.
  w_st <- st$y_upr - st$y_lwr
  w_od <- od$y_upr - od$y_lwr
  j <- which.min(abs(w_st - w_od[1L]))
  if (abs(w_st[j] - w_od[1L]) < 1e-3 * max(abs(w_od[1L]), 1e-8)) {
    expect_equal(od$y_ciw_w[1L], st$y_ciw_w[j], tolerance = 1e-2)
  }
  expect_true(all(od$uncertainty_measure >= 0 & od$uncertainty_measure <= 1))
})
