# Regression tests for the 2026-07-31 bug-fix release (3.3.1).
# Each block corresponds to a finding in the review; see NEWS [3.3.1].

skip_on_cran()

# ---------------------------------------------------------------------------
# H1: SVEMnet(spec$formula, data) must carry the spec (blocking, discrete)
# ---------------------------------------------------------------------------
test_that("spec$formula carries the bigexp_spec attribute into SVEMnet()", {
  set.seed(11)
  n <- 40L
  d <- data.frame(
    y  = rnorm(n),
    X1 = runif(n),
    D  = sample(c(0.5, 1, 2, 4), n, TRUE),
    Op = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  spec <- bigexp_terms(
    y ~ X1 + D, data = d,
    factorial_order = 2, polynomial_order = 2,
    blocking = "Op",
    discrete_numeric = "D"
  )

  expect_s3_class(attr(spec$formula, "bigexp_spec"), "bigexp_spec")
  # no-response bigexp_formula() returns the attributed formula too
  expect_s3_class(attr(bigexp_formula(spec), "bigexp_spec"), "bigexp_spec")

  fit_attr <- SVEMnet(spec$formula, d, nBoot = 8)
  fit_resp <- SVEMnet(bigexp_formula(spec, "y"), d, nBoot = 8)

  expect_identical(fit_attr$sampling_schema$blocking,
                   fit_resp$sampling_schema$blocking)
  expect_true("Op" %in% fit_attr$sampling_schema$blocking)
  expect_true("D" %in% fit_attr$sampling_schema$discrete_numeric)

  # blocking is pinned and discrete support respected in the sampled table
  set.seed(5)
  tab <- svem_random_table_multi(list(y = fit_attr), n = 100)$data
  expect_identical(length(unique(as.character(tab$Op))), 1L)
  expect_true(all(tab$D %in% c(0.5, 1, 2, 4)))
})

# ---------------------------------------------------------------------------
# M5: expansion formulas must not capture the training data environment
# ---------------------------------------------------------------------------
test_that("spec$formula environment does not embed the training data", {
  set.seed(12)
  d <- data.frame(y = rnorm(500), X1 = runif(500), X2 = runif(500))
  spec <- bigexp_terms(y ~ X1 + X2, data = d, 2, 2)
  expect_identical(environment(spec$formula), baseenv())
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(spec, f)
  # a spec with an embedded 500-row frame serialized to ~100 KB+; without it
  # the whole spec is tiny
  expect_lt(file.size(f), 50000)
})

# ---------------------------------------------------------------------------
# H2: glmnet_with_cv() stores contrasts captured before intercept drop
# ---------------------------------------------------------------------------
test_that("glmnet_with_cv() records non-NULL contrasts for factor predictors", {
  set.seed(13)
  n <- 50L
  d <- data.frame(
    y = rnorm(n),
    f = factor(sample(c("a", "b", "c"), n, TRUE)),
    x = runif(n)
  )
  fit <- glmnet_with_cv(y ~ f + x, d, nfolds = 5)
  expect_false(is.null(fit$schema$contrasts))
  expect_true("f" %in% names(fit$schema$contrasts))
})

# ---------------------------------------------------------------------------
# H3: calibration plot draws a non-empty binned reliability curve
# ---------------------------------------------------------------------------
test_that("plot(type='calibration') has finite calibration-line data", {
  set.seed(14)
  n <- 120L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  eta <- 0.8 * d$x1 - 0.6 * d$x2
  d$y <- rbinom(n, 1, plogis(eta))
  fit <- SVEMnet(y ~ x1 + x2, d, nBoot = 10, family = "binomial")
  plt <- plot(fit, type = "calibration", bins = 5)
  # find the calibration-line layer (uses p_mean/y_mean aesthetics)
  layer_datas <- lapply(plt$layers, function(l) l$data)
  has_calib <- vapply(layer_datas, function(dd) {
    is.data.frame(dd) && all(c("p_mean", "y_mean", "n") %in% names(dd)) &&
      nrow(dd) > 0 && all(is.finite(dd$p_mean)) && all(dd$n > 0)
  }, logical(1L))
  expect_true(any(has_calib))
  # '...' overrides no longer collide with the raw-point layer defaults
  expect_s3_class(plot(fit, type = "calibration", alpha = 0.5), "ggplot")
})

# ---------------------------------------------------------------------------
# H4: transformed responses (log(y) ~ .) across the four entry points
# ---------------------------------------------------------------------------
test_that("SVEMnet() fits a transformed response", {
  d <- gen_toy_df(50L)
  d$y <- abs(d$y) + 0.5
  fit <- SVEMnet(log(y) ~ X1 + X2, d, nBoot = 8)
  expect_s3_class(fit, "svem_model")
  pr <- predict(fit, d)
  expect_true(all(is.finite(pr)))
  # response-as-predictor guard still fires (including via the transform);
  # suppressWarnings: base R warns about the dropped RHS response first
  expect_error(suppressWarnings(SVEMnet(y ~ y + X1, d)), "same name as the response")
  expect_error(suppressWarnings(SVEMnet(log(y) ~ y + X1, d)), "same name as the response")
})

test_that("svem_random_table_multi() handles transformed-response models", {
  d <- gen_toy_df(50L)
  d$y <- abs(d$y) + 0.5
  fit <- SVEMnet(log(y) ~ X1 + X2, d, nBoot = 8)
  set.seed(3)
  out <- svem_random_table_multi(list(fit), n = 25)
  expect_true("log(y)_pred" %in% colnames(out$pred))
  expect_true(all(is.finite(out$pred[["log(y)_pred"]])))
})

test_that("svem_wmt_multi() infers names for transformed-response formulas", {
  # name inference must not crash before any WMT runs; use a failing-fast
  # control so the test stays cheap (the WMT itself errors -> NULL entry)
  d <- gen_toy_df(30L)
  d$y <- abs(d$y) + 0.5
  expect_error(
    out <- svem_wmt_multi(
      formulas = list(log(y) ~ X1 + X2),
      data = d,
      plot = FALSE, verbose = FALSE,
      wmt_control = list(nPoint = 10, nSVEM = 2, nPerm = 2, nBoot = 6,
                         nCore = safe_ncores(), seed = 1, verbose = FALSE)
    ),
    NA
  )
  expect_identical(names(out$p_values), "log(y)")
})

test_that("significance test runs with a transformed response", {
  skip_on_ci()
  set.seed(15)
  n <- 40L
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- exp(1 + 2 * d$x1 + rnorm(n, 0, 0.2))
  res <- svem_significance_test_parallel(
    log(y) ~ x1 + x2, d,
    nPoint = 100, nSVEM = 2, nPerm = 25, nBoot = 10,
    nCore = safe_ncores(), seed = 7, verbose = FALSE
  )
  expect_true(is.finite(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------------
# M4: blocking variables pinned in the significance-test evaluation grid
# ---------------------------------------------------------------------------
test_that("significance test rejects blocking mixture vars and runs with blocking", {
  skip_on_ci()
  set.seed(16)
  n <- 40L
  d <- data.frame(
    y  = rnorm(n),
    X1 = runif(n),
    Op = factor(sample(c("A", "B"), n, TRUE))
  )
  spec <- bigexp_terms(y ~ X1, data = d, 2, 2, blocking = "Op")
  expect_error(
    svem_significance_test_parallel(
      y ~ ., d, spec = spec,
      mixture_groups = list(list(vars = c("Op", "X1"))),
      nPoint = 10, nSVEM = 2, nPerm = 20, nBoot = 6,
      nCore = safe_ncores(), seed = 1, verbose = FALSE
    ),
    "blocking"
  )
  res <- svem_significance_test_parallel(
    y ~ ., d, spec = spec,
    nPoint = 100, nSVEM = 2, nPerm = 30, nBoot = 10,
    nCore = safe_ncores(), seed = 2, verbose = FALSE
  )
  expect_true(is.finite(res$p_value))
})

# ---------------------------------------------------------------------------
# M8/M9: degenerate arguments and ordered factors in the significance test
# ---------------------------------------------------------------------------
test_that("nSVEM/nPerm < 1 error up front; ordered factors are categorical", {
  d <- gen_toy_df(30L)
  expect_error(
    svem_significance_test_parallel(y ~ X1, d, nPerm = 0, nCore = 1L),
    "nPerm"
  )
  expect_error(
    svem_significance_test_parallel(y ~ X1, d, nSVEM = 0, nCore = 1L),
    "nSVEM"
  )

  skip_on_ci()
  set.seed(17)
  n <- 40L
  d2 <- data.frame(
    x1   = runif(n),
    dose = factor(sample(c("low", "mid", "high"), n, TRUE),
                  levels = c("low", "mid", "high"), ordered = TRUE)
  )
  d2$y <- 2 * d2$x1 + 0.5 * as.integer(d2$dose) + rnorm(n, 0, 0.3)
  res <- svem_significance_test_parallel(
    y ~ x1 + dose, d2,
    nPoint = 100, nSVEM = 2, nPerm = 25, nBoot = 10,
    nCore = safe_ncores(), seed = 4, verbose = FALSE
  )
  expect_true(is.finite(res$p_value))
})

# ---------------------------------------------------------------------------
# M1: candidate selection when k >= top-set size (pam k < n requirement)
# ---------------------------------------------------------------------------
test_that("svem_select_from_score_table() returns all top rows when k >= top size", {
  set.seed(18)
  n <- 100L
  tab <- data.frame(
    X1 = runif(n), X2 = runif(n),
    score = runif(n)
  )
  attr(tab, "svem_predictor_cols") <- c("X1", "X2")

  # k equal to top-set size: previously crashed in cluster::pam()
  out <- svem_select_from_score_table(tab, target = "score",
                                      k = 5, top_type = "n", top = 5)
  expect_identical(nrow(out$candidates), 5L)

  # 1-row top set: previously crashed in daisy()
  expect_warning(
    out1 <- svem_select_from_score_table(tab, target = "score",
                                         k = 3, top_type = "n", top = 1),
    "top set size"
  )
  expect_identical(nrow(out1$candidates), 1L)

  # single-row table end to end
  tab1 <- tab[1L, , drop = FALSE]
  attr(tab1, "svem_predictor_cols") <- c("X1", "X2")
  expect_warning(
    out2 <- svem_select_from_score_table(tab1, target = "score", k = 2),
    "top set size"
  )
  expect_identical(nrow(out2$best), 1L)
})

# ---------------------------------------------------------------------------
# L6: explicitly supplied predictor_cols must exist
# ---------------------------------------------------------------------------
test_that("svem_select_from_score_table() stops on misspelled predictor_cols", {
  tab <- data.frame(X1 = runif(20), score = runif(20))
  expect_error(
    svem_select_from_score_table(tab, target = "score", k = 2,
                                 predictor_cols = c("Tempp", "Presure")),
    "predictor_cols not found"
  )
})

# ---------------------------------------------------------------------------
# M2: specs with no active bounds; M3: predictor-cols attribute content
# ---------------------------------------------------------------------------
test_that("scoring with inactive specs appends NA joint columns; attrs are clean", {
  set.seed(19)
  n <- 40L
  d <- data.frame(
    X1 = runif(n), X2 = runif(n)
  )
  d$Potency <- 2 * d$X1 + rnorm(n, 0, 0.2)
  d$Size    <- 1 - d$X2 + rnorm(n, 0, 0.2)
  d$RunID   <- sprintf("run_%02d", seq_len(n))

  fits <- list(
    Potency = SVEMnet(Potency ~ X1 + X2, d, nBoot = 8),
    Size    = SVEMnet(Size ~ X1 + X2, d, nBoot = 8)
  )
  set.seed(21)
  scored <- svem_score_random(
    fits,
    goals = list(Potency = list(goal = "max", weight = 0.5),
                 Size    = list(goal = "min", weight = 0.5)),
    n = 50, data = d, verbose = FALSE,
    specs = list(Potency = list(lower = NA))   # no active bound anywhere
  )

  expect_true(all(c("p_joint_mean", "joint_in_spec_point")
                  %in% names(scored$score_table)))
  expect_true(all(is.na(scored$score_table$p_joint_mean)))

  # M3: attribute restricted to design factors (no responses / IDs)
  pc <- attr(scored$original_data_scored, "svem_predictor_cols")
  expect_setequal(pc, c("X1", "X2"))
  expect_false(any(c("Potency", "Size", "RunID") %in% pc))

  # non-numeric bound gives the intended message (spec failures now warn);
  # capture_warnings: as.numeric("abc") also raises a coercion warning
  set.seed(22)
  w <- testthat::capture_warnings(
    svem_score_random(
      fits,
      goals = list(Potency = list(goal = "max", weight = 0.5),
                   Size    = list(goal = "min", weight = 0.5)),
      n = 20, verbose = FALSE,
      specs = list(Potency = list(lower = "abc"))
    )
  )
  expect_true(any(grepl("numeric scalar bounds", w)))
})

# ---------------------------------------------------------------------------
# L7/L8/L9: scoring warnings and anchor validation
# ---------------------------------------------------------------------------
test_that("scoring warns on stale scored columns and bad wmt names; anchors validated", {
  set.seed(20)
  n <- 40L
  d <- data.frame(X1 = runif(n), X2 = runif(n))
  d$Potency <- 2 * d$X1 + rnorm(n, 0, 0.2)
  fits <- list(Potency = SVEMnet(Potency ~ X1 + X2, d, nBoot = 8))
  goals <- list(Potency = list(goal = "max", weight = 1))

  # L7: stale scored column in 'data'
  d_stale <- d
  d_stale$Potency_pred <- 0
  set.seed(23)
  expect_warning(
    svem_score_random(fits, goals = goals, n = 20, data = d_stale,
                      verbose = FALSE),
    "previous scoring round"
  )

  # L8: misaligned wmt multiplier names
  set.seed(24)
  expect_warning(
    svem_score_random(
      fits, goals = goals, n = 20, verbose = FALSE,
      wmt = list(multipliers = c(Potency = 0.5, PDl = 0.9))
    ),
    "match no response"
  )

  # L9: degenerate explicit anchors error out
  set.seed(25)
  expect_error(
    svem_score_random(
      fits,
      goals = list(Potency = list(goal = "max", weight = 1,
                                  lower_acceptable = 5, upper_acceptable = 1)),
      n = 20, verbose = FALSE
    ),
    "degenerate"
  )
})

# ---------------------------------------------------------------------------
# M7: svem_wmt_multi bookkeeping (NULL retention + unique names)
# ---------------------------------------------------------------------------
test_that("svem_wmt_multi keeps NULL entries and uniquifies names safely", {
  d <- gen_toy_df(30L)
  # second formula references a nonexistent column -> that WMT fails -> NULL
  out <- suppressWarnings(svem_wmt_multi(
    formulas = list(y = y ~ X1 + X2, bad = nope ~ X1),
    data = d,
    plot = FALSE, verbose = FALSE,
    wmt_control = list(nPoint = 10, nSVEM = 2, nPerm = 2, nBoot = 6,
                       nCore = safe_ncores(), seed = 1, verbose = FALSE)
  ))
  expect_identical(length(out$wmt_objects), 2L)
  expect_null(out$wmt_objects[["bad"]])
  expect_identical(length(out$p_values), 2L)

  # duplicate names cannot collide with existing user names
  nms <- make.unique(c("y", "y_1", "y"), sep = "_")
  expect_identical(anyDuplicated(nms), 0L)
})

# ---------------------------------------------------------------------------
# M6/L10/L11/L12: bigexp_terms gating and validation
# ---------------------------------------------------------------------------
test_that("factorial_order = 1 yields no interaction terms even with polynomials", {
  d <- data.frame(y = rnorm(30), X1 = runif(30), X2 = runif(30))
  spec <- bigexp_terms(y ~ X1 + X2, data = d,
                       factorial_order = 1, polynomial_order = 2)
  expect_false(grepl(":", spec$rhs, fixed = TRUE))
  expect_true(grepl("I(X1^2)", spec$rhs, fixed = TRUE))

  # 3-way partial cubics need factorial_order >= 3
  d3 <- data.frame(y = rnorm(30), X1 = runif(30), X2 = runif(30), X3 = runif(30))
  spec2 <- bigexp_terms(y ~ X1 + X2 + X3, data = d3,
                        factorial_order = 2, polynomial_order = 2,
                        include_pc_3way = TRUE)
  expect_false(grepl("I(X1^2):X2:X3", spec2$rhs, fixed = TRUE))
  spec3 <- bigexp_terms(y ~ X1 + X2 + X3, data = d3,
                        factorial_order = 3, polynomial_order = 2,
                        include_pc_3way = TRUE)
  expect_true(grepl("I(X1^2):X2:X3", spec3$rhs, fixed = TRUE))
})

test_that("bigexp_terms validates flags and diagnoses non-syntactic names", {
  d <- data.frame(y = rnorm(30), X1 = runif(30), B = factor(sample(c("a","b"), 30, TRUE)))
  expect_error(bigexp_terms(y ~ X1, data = d, intercept = 1), "TRUE or FALSE")
  expect_error(bigexp_terms(y ~ X1, data = d, include_pc_2way = 1), "TRUE or FALSE")

  # dot-prefixed predictor name must not trigger the y ~ . branch
  d2 <- data.frame(y = rnorm(30), .a = runif(30), B = factor(sample(c("a","b"), 30, TRUE)),
                   check.names = FALSE)
  expect_error(
    bigexp_terms(y ~ .a + B, data = d2, blocking = "B"),
    "both in the model RHS and in 'blocking'"
  )

  # non-syntactic column names get the real diagnosis
  d3 <- data.frame(y = rnorm(30), `my var` = runif(30), check.names = FALSE)
  expect_error(bigexp_terms(y ~ ., data = d3), "not syntactically valid")
})

# ---------------------------------------------------------------------------
# L3: a predictor literally named "Intercept" is kept
# ---------------------------------------------------------------------------
test_that("a predictor named Intercept survives fitting", {
  set.seed(26)
  n <- 40L
  d <- data.frame(y = rnorm(n), Intercept = runif(n), X1 = runif(n))
  d$y <- d$y + 3 * d$Intercept
  fit <- SVEMnet(y ~ Intercept + X1, d, nBoot = 8)
  expect_true("Intercept" %in% colnames(fit$training_X))
})

# ---------------------------------------------------------------------------
# L5: cv_object retained only on explicit keep = TRUE
# ---------------------------------------------------------------------------
test_that("glmnet_with_cv keeps cv_object only when keep = TRUE", {
  set.seed(27)
  d <- data.frame(y = rnorm(60), x1 = runif(60), x2 = runif(60))
  fit_default <- glmnet_with_cv(y ~ x1 + x2, d, relaxed = TRUE, nfolds = 5)
  expect_null(fit_default$meta$cv_object)
  fit_keep <- glmnet_with_cv(y ~ x1 + x2, d, relaxed = TRUE, nfolds = 5, keep = TRUE)
  expect_false(is.null(fit_keep$meta$cv_object))
})
