skip_on_cran()

# These tests protect the fitting and prediction conventions documented in
# Karl (2026), CILS 271, 105660, without changing the WMT-specific SVEM-Lasso
# conventions of Karl (2024), CILS 249, 105122.

test_that("one-predictor designs fit without exposing the glmnet sentinel", {
  set.seed(3401)
  n <- 42
  d <- data.frame(x = seq(-1.5, 1.5, length.out = n))
  d$y <- 1 + 2 * d$x + rnorm(n, sd = 0.15)
  d$yb <- as.integer(d$x + rnorm(n, sd = 0.4) > 0)
  d$f <- factor(ifelse(d$x > 0, "high", "low"))

  fits <- list(
    gaussian = SVEMnet(y ~ x, d, nBoot = 4, glmnet_alpha = 1,
                       relaxed = FALSE, penalty.factor = 2),
    binomial = SVEMnet(yb ~ x, d, nBoot = 4, glmnet_alpha = 1,
                       relaxed = FALSE, family = "binomial"),
    factor = SVEMnet(y ~ f, d, nBoot = 4, glmnet_alpha = 1,
                     relaxed = FALSE)
  )

  for (fit in fits) {
    expect_false(any(grepl("svem_internal_zero_sentinel", names(fit$parms))))
    expect_false(any(grepl("svem_internal_zero_sentinel", colnames(fit$training_X))))
    expect_false(any(grepl("svem_internal_zero_sentinel", colnames(fit$coef_matrix))))
    expect_length(predict(fit, d), n)
  }

  cv_g <- glmnet_with_cv(y ~ x, d, glmnet_alpha = 1, nfolds = 5,
                         repeats = 1, seed = 3401, penalty.factor = 2)
  cv_b <- glmnet_with_cv(yb ~ x, d, glmnet_alpha = 1, nfolds = 5,
                         repeats = 1, seed = 3401, family = "binomial")
  cv_f <- glmnet_with_cv(y ~ f, d, glmnet_alpha = 1, nfolds = 5,
                         repeats = 1, seed = 3401)

  for (fit in list(cv_g, cv_b, cv_f)) {
    expect_false(any(grepl("svem_internal_zero_sentinel", names(fit$parms))))
    expect_false(any(grepl("svem_internal_zero_sentinel", colnames(fit$training_X))))
    expect_length(predict(fit, d), n)
  }
})

test_that("fitters enforce response and intercept conventions", {
  set.seed(3402)
  d <- data.frame(y = rnorm(24), z = rnorm(24), x = rnorm(24))

  expect_error(SVEMnet(cbind(y, z) ~ x, d), "univariate response")
  expect_error(glmnet_with_cv(cbind(y, z) ~ x, d), "univariate response")
  expect_error(SVEMnet(y ~ 0 + x, d), "intercept convention")
  expect_error(glmnet_with_cv(y ~ x - 1, d), "intercept convention")

  spec_no_intercept <- bigexp_terms(y ~ x, d, intercept = FALSE)
  expect_error(SVEMnet(spec_no_intercept, d), "intercept convention")

  intercept_fit <- SVEMnet(y ~ 1, d, nBoot = 3, relaxed = FALSE)
  intercept_cv <- glmnet_with_cv(y ~ 1, d, nfolds = 5, repeats = 1)
  expect_identical(names(intercept_fit$parms), "(Intercept)")
  expect_identical(names(intercept_cv$parms), "(Intercept)")
  expect_equal(intercept_fit$diagnostics$fallback_rate, 0)

  dc <- transform(d, y = 4)
  constant_fit <- SVEMnet(y ~ x, dc, nBoot = 3, relaxed = FALSE)
  constant_cv <- glmnet_with_cv(y ~ x, dc, nfolds = 5, repeats = 1)
  expect_equal(predict(constant_fit, dc), rep(4, nrow(dc)))
  expect_equal(predict(constant_cv, dc), rep(4, nrow(dc)))
  expect_equal(constant_fit$diagnostics$fallback_rate, 0)
  expect_identical(constant_cv$schema$feature_names,
                   colnames(constant_cv$training_X))

  db <- transform(d, y = 1L)
  one_class_fit <- SVEMnet(y ~ x, db, nBoot = 3, relaxed = FALSE,
                           family = "binomial")
  one_class_cv <- glmnet_with_cv(y ~ x, db, nfolds = 5, repeats = 1,
                                 family = "binomial")
  expect_true(all(predict(one_class_fit, db) > 0.99))
  expect_true(all(predict(one_class_cv, db) > 0.99))
})

test_that("stored formula environments support transforms without capturing data", {
  set.seed(3403)
  d <- data.frame(x = seq(0.2, 2, length.out = 32))
  d$y <- log(d$x) + rnorm(nrow(d), sd = 0.02)

  poly_fit <- SVEMnet(y ~ poly(x, 2) + scale(x), d, nBoot = 3,
                      glmnet_alpha = 1, relaxed = FALSE)
  before <- predict(poly_fit, d)
  path <- tempfile(fileext = ".rds")
  saveRDS(poly_fit, path)
  poly_fit <- readRDS(path)
  expect_equal(predict(poly_fit, d), before)

  namespaced_fit <- SVEMnet(y ~ stats::poly(x, 2), d, nBoot = 3,
                            glmnet_alpha = 1, relaxed = FALSE)
  namespaced_before <- predict(namespaced_fit, transform(d, x = x + 0.1))
  namespaced_path <- tempfile(fileext = ".rds")
  saveRDS(namespaced_fit, namespaced_path)
  namespaced_fit <- readRDS(namespaced_path)
  expect_equal(predict(namespaced_fit, transform(d, x = x + 0.1)),
               namespaced_before)

  local_fit <- local({
    shift <- 0.35
    unrelated_large_object <- matrix(0, 1000, 50)
    bend <- function(z) (z - shift)^2
    SVEMnet(y ~ bend(x), d, nBoot = 3, glmnet_alpha = 1, relaxed = FALSE)
  })
  term_env <- environment(local_fit$terms)
  expect_true(exists("bend", term_env, inherits = FALSE))
  expect_false(exists("unrelated_large_object", term_env, inherits = FALSE))
  expect_false(exists("unrelated_large_object",
                      environment(get("bend", term_env, inherits = FALSE)),
                      inherits = FALSE))

  local_before <- predict(local_fit, d)
  path_local <- tempfile(fileext = ".rds")
  saveRDS(local_fit, path_local)
  rm(local_fit)
  gc()
  local_fit <- readRDS(path_local)
  expect_equal(predict(local_fit, d), local_before)

  default_fit <- local({
    default_shift <- 0.4
    bend_default <- function(z, shift = default_shift) (z - shift)^2
    SVEMnet(y ~ bend_default(x), d, nBoot = 3,
            glmnet_alpha = 1, relaxed = FALSE)
  })
  default_before <- predict(default_fit, d)
  default_path <- tempfile(fileext = ".rds")
  saveRDS(default_fit, default_path)
  default_fit <- readRDS(default_path)
  expect_equal(predict(default_fit, d), default_before)

  cv_local <- local({
    offset <- 0.2
    bend_cv <- function(z) (z - offset)^2
    glmnet_with_cv(y ~ bend_cv(x), d, glmnet_alpha = 1,
                   nfolds = 5, repeats = 1, seed = 3403)
  })
  expect_length(predict(cv_local, d), nrow(d))
})

test_that("prediction enforces raw predictor classes consistently", {
  set.seed(3404)
  d <- data.frame(x = rnorm(30), f = factor(rep(c("a", "b"), 15)))
  d$y <- d$x + as.numeric(d$f) + rnorm(30, sd = 0.1)
  fit <- SVEMnet(y ~ x + f, d, nBoot = 3, glmnet_alpha = 1,
                 relaxed = FALSE)
  cv <- glmnet_with_cv(y ~ x + f, d, glmnet_alpha = 1,
                       nfolds = 5, repeats = 1, seed = 3404)

  numeric_integer <- transform(d, x = as.integer(round(10 * x)))
  expect_length(predict(fit, numeric_integer), nrow(d))
  expect_length(predict(cv, numeric_integer), nrow(d))

  factor_character <- transform(d, f = as.character(f))
  expect_length(predict(fit, factor_character), nrow(d))
  expect_length(predict(cv, factor_character), nrow(d))

  bad <- transform(d, x = as.character(x))
  expect_error(predict(fit, bad), "fitted as numeric")
  expect_error(predict(cv, bad), "fitted as numeric")
  expect_error(SVEMnet:::.svem_member_predictions(fit, bad), "fitted as numeric")

  ordered_data <- transform(d, f = ordered(f))
  ordered_fit <- SVEMnet(y ~ x + f, ordered_data, nBoot = 3,
                         glmnet_alpha = 1, relaxed = FALSE)
  expect_error(predict(ordered_fit, d), "fitted as ordered")

  empty_missing <- d[0, setdiff(names(d), "x"), drop = FALSE]
  expect_error(predict(fit, empty_missing), "missing required predictor.*x")
  expect_error(predict(cv, empty_missing), "missing required predictor.*x")
})

test_that("CV and random-table controls reject silent coercion", {
  set.seed(3408)
  d <- data.frame(x = rnorm(24), y = rnorm(24))

  expect_error(glmnet_with_cv(y ~ x, d, nfolds = 4.5), "nfolds")
  expect_error(glmnet_with_cv(y ~ x, d, nfolds = 4, seed = 1.5), "seed")
  expect_error(glmnet_with_cv(y ~ x, d, nfolds = 4,
                              glmnet_alpha = c(1, Inf)), "glmnet_alpha")
  expect_error(glmnet_with_cv(y ~ x, transform(d, x = Inf), nfolds = 4),
               "Non-finite")
  expect_error(SVEMnet(y ~ x, d, nBoot = 2, exclude = 1.5),
               "integer index 1")
  expect_error(glmnet_with_cv(y ~ x, d, nfolds = 4, exclude = 1.5),
               "integer index 1")
  expect_warning(
    alpha_fallback <- glmnet_with_cv(y ~ x, d, nfolds = 4, repeats = 1,
                                     relaxed = TRUE, glmnet_alpha = 0,
                                     seed = 3408),
    "Dropping alpha = 0"
  )
  expect_identical(alpha_fallback$glmnet_alpha, 1)

  fit <- SVEMnet(y ~ x, d, nBoot = 2, glmnet_alpha = 1, relaxed = FALSE)
  expect_error(svem_random_table_multi(fit, n = 2.5), "single integer")
  expect_error(svem_random_table_multi(fit, n = 2, range_tol = Inf),
               "single finite number")
  expect_error(svem_random_table_multi(fit, n = 2, debias = 1),
               "TRUE or FALSE")
})

test_that("empty predictions are shaped and uncertainty needs two members", {
  set.seed(3405)
  d <- data.frame(x = rnorm(25), y = rnorm(25))
  fit <- SVEMnet(y ~ x, d, nBoot = 2, glmnet_alpha = 1, relaxed = FALSE)
  empty <- d[0, , drop = FALSE]

  expect_identical(predict(fit, empty), numeric(0L))
  out <- predict(fit, empty, se.fit = TRUE, interval = TRUE)
  expect_named(out, c("fit", "se.fit", "lwr", "upr"))
  expect_true(all(lengths(out) == 0L))
  members_empty <- SVEMnet:::.svem_member_predictions(fit, empty)
  expect_equal(dim(members_empty$pred_mat), c(0L, 2L))

  cv <- glmnet_with_cv(y ~ x, d, glmnet_alpha = 1,
                       nfolds = 5, repeats = 1, seed = 3405)
  expect_identical(predict(cv, empty), numeric(0L))

  one_member <- fit
  one_member$coef_matrix <- fit$coef_matrix[1, , drop = FALSE]
  expect_error(predict(one_member, d, se.fit = TRUE), "at least two")
  expect_error(predict(one_member, d, interval = TRUE), "at least two")
  expect_length(predict(one_member, d), nrow(d))

  db <- transform(d, y = as.integer(x > 0))
  fit_b <- SVEMnet(y ~ x, db, nBoot = 2, glmnet_alpha = 1,
                   relaxed = FALSE, family = "binomial")
  expect_identical(predict(fit_b, db[0, , drop = FALSE], type = "class"),
                   integer(0L))
})

test_that("all non-structural bootstrap fallbacks are an error", {
  set.seed(3406)
  d <- data.frame(x = rnorm(20), y = rnorm(20))
  expect_error(
    SVEMnet(y ~ x, d, nBoot = 2, relaxed = FALSE, exclude = 1),
    "All bootstrap fits fell back.*Reasons"
  )
})

test_that("relaxed CV jointly selects lambda and gamma like cv.glmnet", {
  set.seed(1)
  n <- 60
  p <- 8
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  beta <- c(3, -2, 1.5, rep(0, p - 3L))
  y <- drop(x %*% beta + rnorm(n, sd = 2))
  d <- data.frame(y, x)
  gamma_grid <- c(0, 0.5, 1)

  set.seed(901)
  foldid <- sample(rep(seq_len(5), length.out = n))
  direct <- glmnet::cv.glmnet(
    x, y, alpha = 1, foldid = foldid, relax = TRUE,
    gamma = gamma_grid
  )

  fit_min <- glmnet_with_cv(
    y ~ ., d, glmnet_alpha = 1, nfolds = 5, repeats = 1,
    choose_rule = "min", seed = 901, relaxed = TRUE,
    relax_gamma = gamma_grid
  )
  fit_1se <- glmnet_with_cv(
    y ~ ., d, glmnet_alpha = 1, nfolds = 5, repeats = 1,
    choose_rule = "1se", seed = 901, relaxed = TRUE,
    relax_gamma = gamma_grid
  )

  expect_equal(fit_min$best_lambda, direct$relaxed$lambda.min)
  expect_equal(fit_min$best_gamma, direct$relaxed$gamma.min)
  expect_equal(fit_1se$best_lambda, direct$relaxed$lambda.1se)
  expect_equal(fit_1se$best_gamma, direct$relaxed$gamma.1se)
  direct_coef <- drop(as.matrix(stats::coef(
    direct$glmnet.fit,
    s = direct$relaxed$lambda.min,
    gamma = direct$relaxed$gamma.min
  )))
  expect_equal(unname(fit_min$parms), unname(direct_coef))
  expect_lt(fit_min$best_gamma, 1)
  expect_false(isTRUE(all.equal(fit_min$best_lambda, direct$lambda.min)))

  summary_min <- fit_min$cv_summary[["1"]]
  expect_setequal(unique(summary_min$gamma), gamma_grid)
  expect_true(all(summary_min$n_repeats == 1L))
  expect_equal(fit_min$diagnostics$selected_gamma, fit_min$best_gamma)
  expect_match(fit_min$note, "jointly selected")
})

test_that("repeated relaxed CV aggregates every gamma surface", {
  set.seed(3410)
  n <- 54
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  d$y <- 2 * d$x1 - d$x2 + rnorm(n, sd = 1.2)

  fit <- glmnet_with_cv(
    y ~ ., d, glmnet_alpha = 1, nfolds = 5, repeats = 2,
    seed = 3410, relaxed = TRUE, relax_gamma = c(0, 0.5, 1)
  )
  sm <- fit$cv_summary[["1"]]

  expect_setequal(unique(sm$gamma), c(0, 0.5, 1))
  expect_true(all(sm$n_repeats == 2L))
  expect_true(fit$best_gamma %in% c(0, 0.5, 1))
})
