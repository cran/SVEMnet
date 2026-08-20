skip_on_cran()

test_that("binomial ROC AUC and AP are invariant to ordering within ties", {
  make_binomial_plot_object <- function(y, p) {
    structure(
      list(actual_y = y, y_pred = p),
      class = c("svem_binomial", "svem_model")
    )
  }

  y <- c(1, 0, 0, 1, 1, 0, 0, 1)
  p <- c(0.9, 0.9, 0.6, 0.6, 0.3, 0.3, 0.1, 0.1)
  within_tie_permutation <- c(2, 1, 4, 3, 6, 5, 8, 7)

  roc_1 <- plot(make_binomial_plot_object(y, p), type = "roc")
  roc_2 <- plot(make_binomial_plot_object(
    y[within_tie_permutation], p[within_tie_permutation]
  ), type = "roc")
  pr_1 <- plot(make_binomial_plot_object(y, p), type = "pr")
  pr_2 <- plot(make_binomial_plot_object(
    y[within_tie_permutation], p[within_tie_permutation]
  ), type = "pr")

  expect_equal(roc_1$data, roc_2$data)
  expect_identical(roc_1$labels$title, roc_2$labels$title)
  expect_equal(pr_1$data, pr_2$data)
  expect_identical(pr_1$labels$title, pr_2$labels$title)

  constant <- make_binomial_plot_object(c(1, 0, 1, 0, 0), rep(0.4, 5))
  expect_identical(plot(constant, type = "roc")$labels$title,
                   "ROC (AUC = 0.500)")
  expect_identical(plot(constant, type = "pr")$labels$title,
                   "Precision-Recall (AP = 0.400)")
})

test_that("plot controls reject fractional and non-scalar values", {
  binomial_object <- structure(
    list(actual_y = c(0, 1), y_pred = c(0.2, 0.8)),
    class = c("svem_binomial", "svem_model")
  )
  gaussian_object <- structure(
    list(actual_y = c(0, 1), y_pred = c(0.1, 0.9)),
    class = "svem_model"
  )

  expect_error(plot(binomial_object, bins = 2.5), "single integer")
  expect_error(plot(binomial_object, jitter_width = c(0.1, 0.2)),
               "single finite number")
  expect_error(plot(gaussian_object, plot_debiased = NA), "TRUE or FALSE")
})

test_that("scoring and Thompson count controls require exact integers", {
  fake <- structure(
    list(
      formula = y ~ x,
      family = "gaussian",
      coef_matrix = matrix(0, nrow = 2, ncol = 2)
    ),
    class = "svem_model"
  )
  goal <- list(y = list(goal = "max", weight = 1))

  expect_error(svem_score_random(fake, goal, n = 1.5, verbose = FALSE),
               "single integer")
  expect_error(svem_score_random(fake, goal, n = 2, verbose = 1),
               "TRUE or FALSE")
  expect_error(svem_thompson_batch(fake, goal, batch_size = 1.5,
                                   candidates = data.frame(x = 0),
                                   verbose = FALSE),
               "single integer")
  expect_error(svem_thompson_batch(fake, goal, batch_size = 1, n = 2.5,
                                   verbose = FALSE),
               "single integer")
})

test_that("goal aliases and desirability controls fail clearly when ambiguous", {
  fake <- structure(
    list(
      formula = y ~ x,
      family = "gaussian",
      coef_matrix = matrix(0, nrow = 2, ncol = 2)
    ),
    class = "svem_model"
  )
  objects <- list(alias = fake)

  fake_z <- fake
  fake_z$formula <- z ~ x
  expect_error(
    svem_score_random(
      list(z = fake, other = fake_z),
      goals = list(z = list(goal = "max", weight = 1),
                   other = list(goal = "max", weight = 1)),
      n = 2, verbose = FALSE
    ),
    "Ambiguous response alias"
  )
  expect_error(
    svem_thompson_batch(
      list(z = fake, other = fake_z),
      goals = list(z = list(goal = "max", weight = 1),
                   other = list(goal = "max", weight = 1)),
      batch_size = 1, candidates = data.frame(x = 0), verbose = FALSE
    ),
    "Ambiguous response alias"
  )

  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(alias = list(goal = "max", weight = 1),
                   y = list(goal = "min", weight = 1)),
      n = 2, verbose = FALSE
    )),
    "duplicate responses"
  )
  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(typo = list(goal = "max", weight = 1)),
      n = 2, verbose = FALSE
    )),
    "Unknown goal response"
  )
  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(alias = list(goal = "max", weight = 0)),
      n = 2, verbose = FALSE
    )),
    "strictly positive"
  )
  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(alias = list(goal = "target", weight = 1,
                                target = c(1, 2))),
      n = 2, verbose = FALSE
    )),
    "target.*single finite number"
  )
  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(alias = list(goal = "max", weight = 1,
                                target = c(1, 2))),
      n = 2, verbose = FALSE
    )),
    "target.*single finite number"
  )
  expect_error(
    suppressWarnings(svem_score_random(
      objects,
      goals = list(alias = list(goal = "target", weight = 1,
                                target = 1, tol_left = 0)),
      n = 2, verbose = FALSE
    )),
    "tol_left.*single finite number"
  )
  expect_error(
    suppressWarnings(svem_thompson_batch(
      objects,
      goals = list(alias = list(goal = "target", weight = 1,
                                target = 1, lower_acceptable = 1)),
      batch_size = 1, candidates = data.frame(x = 0), verbose = FALSE
    )),
    "lower_acceptable < target"
  )
})

test_that("transformed response labels work in scoring and Thompson proposals", {
  set.seed(3401)
  d <- data.frame(x1 = runif(36), x2 = runif(36))
  d$y <- exp(0.5 + 1.2 * d$x1 - 0.4 * d$x2 + rnorm(36, sd = 0.08))
  fit <- SVEMnet(log(y) ~ x1 + x2, d, nBoot = 6,
                 glmnet_alpha = 1, relaxed = FALSE)
  goals <- list(list(goal = "max", weight = 1))

  set.seed(3402)
  scored <- svem_score_random(list(fit), goals, n = 12,
                              numeric_sampler = "uniform", verbose = FALSE)
  expect_true(all(c("log(y)_pred", "log(y)_des") %in%
                    names(scored$score_table)))

  set.seed(3403)
  proposed <- svem_thompson_batch(
    list(fit), goals, batch_size = 1,
    candidates = d[1:8, c("x1", "x2")], verbose = FALSE
  )
  expect_true(all(c("log(y)_draw", "log(y)_des", "log(y)_pred") %in%
                    names(proposed$proposals)))
})

test_that("candidate selection ranks only finite target values", {
  tab <- data.frame(
    x = 1:5,
    score = c(Inf, 0.8, NA_real_, 0.5, -Inf)
  )
  attr(tab, "svem_predictor_cols") <- "x"

  selected <- svem_select_from_score_table(
    tab, target = "score", direction = "max",
    k = 2, top_type = "n", top = 2
  )
  expect_identical(selected$best$x, 2L)
  expect_setequal(selected$candidates$x, c(2L, 4L))

  direct <- SVEMnet:::svem_select_candidates(
    tab, by = "score", top_frac = 1, k = 2,
    predictor_cols = "x", direction = "max"
  )
  expect_setequal(direct, c(2L, 4L))

  all_bad <- transform(tab, score = c(NA_real_, Inf, -Inf, NaN, NA_real_))
  expect_error(
    svem_select_from_score_table(all_bad, target = "score", k = 0),
    "non-finite"
  )
  expect_error(
    SVEMnet:::svem_select_candidates(
      all_bad, by = "score", top_frac = 1, k = 1,
      predictor_cols = "x"
    ),
    "no finite values"
  )
  expect_error(
    svem_select_from_score_table(tab, target = "score", k = 1.5),
    "single integer"
  )
  expect_error(
    svem_select_from_score_table(tab, target = "score", k = 1,
                                 top_type = "n", top = 1.5),
    "single integer"
  )
})
