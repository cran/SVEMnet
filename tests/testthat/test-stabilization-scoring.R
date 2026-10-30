.stabilization_scoring_fit <- function(nBoot = 3L) {
  d <- expand.grid(x = seq(-1, 1, length.out = 8), z = c(-1, 1))
  d$y <- 2 + 1.5 * d$x - .4 * d$z + .03 * sin(seq_len(nrow(d)))
  set.seed(481)
  list(data = d, fit = SVEMnet(y ~ x + z, d, nBoot = nBoot,
                              glmnet_alpha = 1, relaxed = FALSE))
}

test_that("scoring refuses generated-name collisions before consuming RNG", {
  a <- .stabilization_scoring_fit()
  goals <- list(y = list(goal = "max", weight = 1))
  for (name in c("score", "y_des", "y_lwr", "p_joint_mean")) {
    d <- a$data
    names(d)[names(d) == "x"] <- name
    fit <- SVEMnet(reformulate(c(name, "z"), "y"), d, nBoot = 2,
                   glmnet_alpha = 1, relaxed = FALSE)
    before <- .Random.seed
    expect_error(svem_score_random(list(y = fit), goals, n = 4, verbose = FALSE),
                 "would overwrite input variable")
    expect_identical(.Random.seed, before)
  }
  d <- a$data
  names(d)[names(d) == "x"] <- "alias_pred"
  fit <- SVEMnet(y ~ alias_pred + z, d, nBoot = 2, relaxed = FALSE)
  expect_error(svem_score_random(list(alias = fit),
                                 list(alias = list(goal = "max", weight = 1)),
                                 n = 4, verbose = FALSE), "alias_pred")

  # A second model's observed response must not be mistaken for stale output.
  d <- a$data
  d$y_pred <- d$y + 1
  fit2 <- SVEMnet(y_pred ~ x + z, d, nBoot = 2, relaxed = FALSE)
  expect_error(svem_score_random(list(y = a$fit, y_pred = fit2),
                                 list(y = goals$y, y_pred = goals$y),
                                 data = d, n = 4, verbose = FALSE),
               "supplied data.*y_pred")
})

test_that("rescoring replaces old outputs and preserves measurements", {
  a <- .stabilization_scoring_fit()
  d <- a$data
  d$y_pred <- -999
  d$y_des <- -888
  d$y_in_spec_point <- 0L
  d$p_joint_mean <- -1
  d$wmt_score <- -2
  d$Notes <- "measured"
  expect_warning(out <- svem_score_random(
    list(y = a$fit), list(y = list(goal = "max", weight = 1)),
    data = d, n = 12, numeric_sampler = "uniform", verbose = FALSE,
    specs = list(y = list(lower = -10, upper = 10))
  ), "replacing generated column.*previous scoring round")
  tab <- out$original_data_scored
  expect_identical(tab$y, d$y)
  expect_identical(tab$Notes, d$Notes)
  expect_equal(tab$y_pred, as.numeric(predict(a$fit, d)))
  expect_true(all(tab$y_in_spec_point == 1L))
  expect_true(all(tab$p_joint_mean == 1))
  expect_false("wmt_score" %in% names(tab))
  expect_false(any(grepl("[.]1$", names(tab))))
  expect_identical(attr(tab, "svem_predictor_cols"), c("x", "z"))
})

test_that("duplicate specification aliases fail before sampling", {
  a <- .stabilization_scoring_fit()
  before <- .Random.seed
  expect_error(svem_score_random(
    list(alias = a$fit), list(alias = list(goal = "max", weight = 1)),
    n = 4, verbose = FALSE,
    specs = list(alias = list(lower = 0), y = list(upper = 3))
  ), "Spec names/aliases map to duplicate responses")
  expect_identical(.Random.seed, before)
})

test_that("logical and directly factored numerics retain their sampling support", {
  d <- expand.grid(x = seq(-1, 1, length.out = 6), flag = c(FALSE, TRUE),
                   level = c(1, 2, 4))
  d$y <- 2 + d$x + d$flag + d$level + .01 * sin(seq_len(nrow(d)))
  set.seed(483)
  fit <- SVEMnet(y ~ x + flag + factor(level), d, nBoot = 3,
                 glmnet_alpha = 1, relaxed = FALSE)
  tab <- svem_random_table_multi(fit, n = 20, numeric_sampler = "uniform")
  expect_type(tab$data$flag, "logical")
  expect_true(all(tab$data$level %in% c(1, 2, 4)))
  expect_true(all(is.finite(tab$pred$y_pred)))
  expect_equal(tab$pred$y_pred, as.numeric(predict(fit, tab$data)))

  spec <- bigexp_terms(y ~ x, d, blocking = "flag", report = FALSE)
  fit_spec <- SVEMnet(spec, d, nBoot = 2, relaxed = FALSE)
  sampled <- svem_random_table_multi(fit_spec, n = 6)
  expect_type(sampled$data$flag, "logical")
  expect_length(unique(sampled$data$flag), 1L)
  expect_true(all(is.finite(sampled$pred$y_pred)))
})

test_that("Thompson and selection metadata cannot overwrite design columns", {
  a <- .stabilization_scoring_fit()
  goals <- list(y = list(goal = "max", weight = 1))
  for (name in c("slot", "ts_score", "y_pred", "y_draw", "y_des")) {
    cand <- a$data[c("x", "z")]
    cand[[name]] <- seq_len(nrow(cand))
    before <- .Random.seed
    expect_error(svem_thompson_batch(list(y = a$fit), goals, batch_size = 2,
                                     candidates = cand, verbose = FALSE),
                 "would overwrite input variable")
    expect_identical(.Random.seed, before)
  }
  tab <- data.frame(Notes_from_SVEMnet = 1:3, score = c(.1, .8, .4))
  expect_error(svem_select_from_score_table(tab, k = 0,
               predictor_cols = "Notes_from_SVEMnet", label = "new"),
               "would overwrite input variable")
  for (name in c("candidate_type", "selection_label")) {
    tab <- data.frame(x = 1:3, score = c(.1, .8, .4))
    names(tab)[1] <- name
    sel <- svem_select_from_score_table(tab, k = 0, predictor_cols = name)
    expect_error(svem_export_candidates_csv(sel, write_file = FALSE),
                 "would overwrite input variable")
  }
})

test_that("unavailable intervals are identified and keep the documented fallback", {
  a <- .stabilization_scoring_fit(nBoot = 1)
  expect_warning(out <- svem_score_random(
    list(y = a$fit), list(y = list(goal = "max", weight = 1)),
    n = 8, data = a$data, verbose = FALSE
  ), "constant 0.5 times its user weight")
  expect_true(all(is.na(out$score_table$y_lwr)))
  expect_equal(out$score_table$uncertainty_measure, rep(.5, 8))
  expect_true(all(is.finite(out$original_data_scored$y_pred)))
})

test_that("unexpected interval failures no longer pass silently", {
  a <- .stabilization_scoring_fit()
  # Point prediction uses the intact aggregated coefficients; interval
  # prediction must report the inconsistent member-coefficient schema.
  colnames(a$fit$coef_matrix)[2L] <- "missing_member_column"
  expect_warning(out <- svem_score_random(
    list(y = a$fit), list(y = list(goal = "max", weight = 1)),
    n = 6, verbose = FALSE
  ), "interval prediction failed for response 'y'.*score_table")
  expect_true(all(is.finite(out$score_table$y_pred)))
  expect_true(all(is.na(out$score_table$y_lwr)))
  expect_equal(out$score_table$uncertainty_measure, rep(.5, 6))
})

test_that("intercept-only ensembles support random tables and specification scoring", {
  set.seed(492)
  fit <- SVEMnet(y ~ 1, data.frame(y = seq_len(12)), nBoot = 3)
  tab <- svem_random_table_multi(fit, n = 5)
  expect_equal(dim(tab$data), c(5L, 0L))
  expect_length(tab$pred$y_pred, 5L)
  expect_equal(tab$pred$y_pred, rep(mean(fit$coef_matrix[, 1]), 5))
  out <- svem_score_random(
    list(y = fit), list(y = list(goal = "max", weight = 1)), n = 5,
    specs = list(y = list(lower = 0, upper = 20)), verbose = FALSE
  )
  expect_equal(out$score_table$y_p_in_spec_mean, rep(1, 5))
  expect_true(all(out$score_table$y_in_spec_point))
})
