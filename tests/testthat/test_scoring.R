skip_on_cran()
skip_if_not_installed("SVEMnet")

test_that("svem_score_random basic scoring without WMT works", {
  set.seed(1)
  n <- 80
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 20)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 20)

  objs <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  res <- svem_score_random(
    objects         = objs,
    goals           = goals,
    n               = 1000,
    level           = 0.95,
    combine         = "geom",
    numeric_sampler = "random",
    verbose         = FALSE
  )

  expect_type(res, "list")
  expect_true(all(c("score_table",
                    "original_data_scored",
                    "weights_original",
                    "weights_final",
                    "wmt_p_values",
                    "wmt_multipliers") %in% names(res)))

  st <- res$score_table
  expect_s3_class(st, "data.frame")
  expect_equal(nrow(st), 1000L)
  expect_true(all(c("y1_pred", "y2_pred", "score", "uncertainty_measure") %in% colnames(st)))

  # CI & desirability columns present
  expect_true(all(c("y1_des", "y2_des",
                    "y1_lwr", "y1_upr",
                    "y2_lwr", "y2_upr",
                    "y1_ciw_w", "y2_ciw_w") %in% colnames(st)))

  # Weights sane and no WMT reweighting when wmt = NULL
  expect_equal(sum(res$weights_original), 1)
  expect_equal(res$weights_original, res$weights_final)

  # WMT bits should be "inactive"
  if (!is.null(res$wmt_p_values)) {
    expect_true(all(is.na(res$wmt_p_values)))
  }
  if (!is.null(res$wmt_multipliers)) {
    expect_true(all(res$wmt_multipliers == 1))
  }
})


test_that("svem_score_random uses external WMT multipliers when provided", {
  set.seed(2)
  n <- 80
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 20)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 20)

  objs <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  # Fake WMT object, shaped like svem_wmt_multi() output
  wmt_fake <- list(
    p_values    = c(y1 = 1e-4, y2 = 0.2),
    multipliers = c(y1 = 5,    y2 = 1)
  )

  res <- svem_score_random(
    objects         = objs,
    goals           = goals,
    n               = 500,
    numeric_sampler = "random",
    wmt             = wmt_fake,
    verbose         = FALSE
  )

  # weights_original normalized
  expect_equal(sum(res$weights_original), 1)

  # weights_final should be proportional to weights_original * multipliers
  w_user <- res$weights_original[c("y1", "y2")]
  mult   <- wmt_fake$multipliers[c("y1", "y2")]
  exp_final_raw <- w_user * mult
  exp_final <- exp_final_raw / sum(exp_final_raw)

  expect_equal(res$weights_final[c("y1", "y2")],
               exp_final,
               tolerance = 1e-8)

  # wmt_score present and different from plain score (since multipliers differ)
  st <- res$score_table
  expect_true("wmt_score" %in% colnames(st))
  expect_false(all(st$score == st$wmt_score))

  # p-values and multipliers returned with compatible names
  expect_true(all(names(res$wmt_p_values) %in% names(wmt_fake$p_values)))
  expect_true(all(names(res$wmt_multipliers) %in% names(wmt_fake$multipliers)))
})

test_that("svem_select_from_score_table picks best and k candidates", {
  set.seed(3)
  n <- 80
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 15)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 15)

  objs <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  scored <- svem_score_random(
    objects         = objs,
    goals           = goals,
    n               = 2000,
    numeric_sampler = "random",
    verbose         = FALSE
  )

  st <- scored$score_table

  sel <- svem_select_from_score_table(
    score_table = st,
    target      = "score",
    direction   = "max",
    k           = 5,
    top_type    = "frac",
    top         = 0.10
  )

  expect_true(is.list(sel))
  expect_true(all(c("best", "candidates") %in% names(sel)))

  best <- sel$best
  cand <- sel$candidates

  expect_s3_class(best, "data.frame")
  expect_equal(nrow(best), 1L)
  idx_max <- which.max(st$score)
  expect_equal(best$score, st$score[idx_max])

  expect_s3_class(cand, "data.frame")
  expect_equal(nrow(cand), 5L)

  # Candidates should come from the input table (at least by score)
  expect_true(all(cand$score %in% st$score))

  # No obvious duplicates in candidates
  expect_equal(nrow(unique(cand)), nrow(cand))
})

test_that("svem_select_from_score_table can target uncertainty_measure", {
  set.seed(4)
  n <- 60
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 10)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 10)

  objs <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  scored <- svem_score_random(
    objects         = objs,
    goals           = goals,
    n               = 1000,
    numeric_sampler = "random",
    verbose         = FALSE
  )

  st <- scored$score_table

  sel_unc <- svem_select_from_score_table(
    score_table = st,
    target      = "uncertainty_measure",
    direction   = "max",
    k           = 3,
    top_type    = "frac",
    top         = 0.10
  )

  expect_s3_class(sel_unc$best, "data.frame")
  expect_equal(nrow(sel_unc$best), 1L)

  max_u <- max(st$uncertainty_measure, na.rm = TRUE)
  expect_equal(sel_unc$best$uncertainty_measure, max_u)
})

test_that("svem_score_random adds spec-based design-space columns", {
  skip_on_cran()
  skip_if_not_installed("SVEMnet")

  set.seed(5)
  n  <- 80
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 20)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 20)

  objs  <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  specs <- list(
    y1 = list(lower = 75),   # one-sided lower
    y2 = list(upper = 70)    # one-sided upper
  )

  res <- svem_score_random(
    objects         = objs,
    goals           = goals,
    data            = dat,
    n               = 500,
    numeric_sampler = "random",
    specs           = specs,
    verbose         = FALSE
  )

  st <- res$score_table

  # Per-response columns present
  expect_true(all(c("y1_p_in_spec_mean",
                    "y2_p_in_spec_mean",
                    "y1_in_spec_point",
                    "y2_in_spec_point") %in% names(st)))

  # Joint columns present and sensible
  expect_true(all(c("p_joint_mean", "joint_in_spec_point") %in% names(st)))
  expect_true(all(st$y1_in_spec_point %in% c(0L, 1L, NA_integer_)))
  expect_true(all(st$y2_in_spec_point %in% c(0L, 1L, NA_integer_)))

  # When both point indicators are 1, joint_in_spec_point must be 1
  idx_both_in <- which(st$y1_in_spec_point == 1L & st$y2_in_spec_point == 1L)
  if (length(idx_both_in)) {
    expect_true(all(st$joint_in_spec_point[idx_both_in] == 1L))
  }
})


test_that("svem_score_random returns scored original data when data is supplied", {
  skip_on_cran()
  skip_if_not_installed("SVEMnet")

  set.seed(6)
  n  <- 40
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  y1 <- 80 + 10 * X1 - 5 * X2 + rnorm(n, 0, 3)
  y2 <- 60 +  5 * X1 + 2 * X2 + rnorm(n, 0, 4)
  dat <- data.frame(y1 = y1, y2 = y2, X1 = X1, X2 = X2)

  fit_y1 <- SVEMnet(y1 ~ (X1 + X2)^2, dat, nBoot = 15)
  fit_y2 <- SVEMnet(y2 ~ (X1 + X2)^2, dat, nBoot = 15)

  objs  <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  res <- svem_score_random(
    objects = objs,
    goals   = goals,
    data    = dat,
    n       = 500,
    verbose = FALSE
  )

  od <- res$original_data_scored
  expect_s3_class(od, "data.frame")
  expect_equal(nrow(od), n)

  expect_true(all(c("y1_pred", "y2_pred",
                    "y1_des", "y2_des",
                    "score", "uncertainty_measure") %in% names(od)))
})


test_that("svem_score_random errors when WMT is requested for binomial responses", {
  skip_on_cran()
  skip_if_not_installed("SVEMnet")

  set.seed(7)
  n  <- 60
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  p  <- plogis(0.5 + X1 - 0.5 * X2)
  yb <- rbinom(n, size = 1, prob = p)

  dat <- data.frame(yb = yb, X1 = X1, X2 = X2)

  fit_bin <- SVEMnet(yb ~ (X1 + X2)^2, dat, family = "binomial", nBoot = 10)

  objs  <- list(yb = fit_bin)
  goals <- list(
    yb = list(goal = "max", weight = 1)
  )

  wmt_fake <- list(
    p_values    = c(yb = 0.01),
    multipliers = c(yb = 2)
  )

  expect_error(
    svem_score_random(
      objects = objs,
      goals   = goals,
      n       = 200,
      wmt     = wmt_fake,
      verbose = FALSE
    ),
    "binomial",
    ignore.case = TRUE
  )
})

test_that("svem_score_random works with all-categorical predictors", {
  skip_on_cran()
  skip_if_not_installed("SVEMnet")

  set.seed(8)
  n <- 80

  F1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  F2 <- factor(sample(c("low", "high"), n, replace = TRUE))
  F3 <- factor(sample(c("X", "Y"), n, replace = TRUE))

  # Simple signal that depends only on factors
  mu1 <- 80 +
    ifelse(F1 == "B", 5, 0) +
    ifelse(F2 == "high", 3, 0) +
    ifelse(F3 == "Y", 2, 0)

  mu2 <- 60 +
    ifelse(F1 == "C", -4, 0) +
    ifelse(F2 == "high", 2, 0)

  y1 <- mu1 + rnorm(n, 0, 3)
  y2 <- mu2 + rnorm(n, 0, 4)

  dat <- data.frame(y1 = y1, y2 = y2, F1 = F1, F2 = F2, F3 = F3)

  fit_y1 <- SVEMnet(y1 ~ F1 * F2 * F3, dat, nBoot = 15)
  fit_y2 <- SVEMnet(y2 ~ F1 * F2 * F3, dat, nBoot = 15)

  objs <- list(y1 = fit_y1, y2 = fit_y2)
  goals <- list(
    y1 = list(goal = "max", weight = 0.5),
    y2 = list(goal = "max", weight = 0.5)
  )

  res <- svem_score_random(
    objects         = objs,
    goals           = goals,
    data            = dat,
    n               = 500,
    level           = 0.95,
    # numeric_sampler is irrelevant here, but should not break
    numeric_sampler = "random",
    verbose         = FALSE
  )

  expect_type(res, "list")
  expect_true(all(c("score_table",
                    "original_data_scored",
                    "weights_original",
                    "weights_final",
                    "wmt_p_values",
                    "wmt_multipliers") %in% names(res)))

  st <- res$score_table
  expect_s3_class(st, "data.frame")
  expect_equal(nrow(st), 500L)

  # Core columns present (predictions are *_pred in the score table)
  expect_true(all(c("y1_pred", "y2_pred", "score", "uncertainty_measure") %in% colnames(st)))

  # CIs & desirabilities present
  expect_true(all(c("y1_des", "y2_des",
                    "y1_lwr", "y1_upr",
                    "y2_lwr", "y2_upr") %in% colnames(st)))

  # All predictors in score_table should still be factors
  for (f in c("F1", "F2", "F3")) {
    expect_true(f %in% names(st))
    expect_true(is.factor(st[[f]]))
  }

  # Original data scored and factor structure preserved
  od <- res$original_data_scored
  expect_s3_class(od, "data.frame")
  expect_equal(nrow(od), n)
  for (f in c("F1", "F2", "F3")) {
    expect_true(f %in% names(od))
    expect_true(is.factor(od[[f]]))
  }

  expect_true(all(c("y1_pred", "y2_pred",
                    "y1_des",  "y2_des",
                    "score",   "uncertainty_measure",
                    # new CI-related columns on original_data_scored:
                    "y1_lwr",  "y1_upr",  "y1_ciw_w",
                    "y2_lwr",  "y2_upr",  "y2_ciw_w") %in% names(od)))

  # CI columns are numeric and aligned with rows
  for (nm in c("y1_lwr", "y1_upr", "y1_ciw_w",
               "y2_lwr", "y2_upr", "y2_ciw_w")) {
    expect_type(od[[nm]], "double")
    expect_equal(length(od[[nm]]), n)
  }

})
