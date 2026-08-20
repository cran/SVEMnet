# Tests for svem_thompson_batch()
skip_on_cran()

fit_two_models <- function(with_factor = FALSE, nBoot = 25L) {
  df <- gen_toy_df(n = 50L, with_factor = with_factor)
  set.seed(11)
  df$y2 <- 3 - 1.5 * df$X1 + 0.8 * df$X2 + rnorm(nrow(df), 0, 0.2)
  rhs <- if (with_factor) "X1 + X2 + X3 + F" else "X1 + X2 + X3"
  set.seed(42)
  f1 <- SVEMnet(stats::as.formula(paste("y ~", rhs)), df, nBoot = nBoot)
  f2 <- SVEMnet(stats::as.formula(paste("y2 ~", rhs)), df, nBoot = nBoot)
  list(objects = list(y = f1, y2 = f2), data = df)
}

basic_goals <- list(
  y  = list(goal = "max", weight = 0.7),
  y2 = list(goal = "min", weight = 0.3)
)

test_that("svem_thompson_batch returns a well-formed batch", {
  fits <- fit_two_models()
  set.seed(101)
  prop <- svem_thompson_batch(fits$objects, basic_goals,
                              batch_size = 4, n = 400, verbose = FALSE)

  expect_s3_class(prop, "svem_thompson_batch")
  expect_equal(nrow(prop$proposals), 4L)
  expect_true(all(c("X1", "X2", "X3", "slot", "ts_score",
                    "y_draw", "y_des", "y_pred",
                    "y2_draw", "y2_des", "y2_pred") %in% names(prop$proposals)))
  expect_equal(prop$proposals$slot, 1:4)
  # no candidate proposed twice
  expect_equal(anyDuplicated(prop$candidate_index), 0L)
  # scores/desirabilities in [0, 1]
  expect_true(all(prop$proposals$ts_score >= 0 & prop$proposals$ts_score <= 1))
  expect_true(all(prop$proposals$y_des  >= 0 & prop$proposals$y_des  <= 1))
  expect_true(all(prop$proposals$y2_des >= 0 & prop$proposals$y2_des <= 1))
  # members are valid bootstrap indices
  expect_true(all(prop$members >= 1))
  expect_true(all(prop$members[, "y"] <= nrow(fits$objects$y$coef_matrix)))
  expect_identical(colnames(prop$members), c("y", "y2"))
  # weights normalized
  expect_equal(sum(prop$weights), 1)
  expect_equal(unname(prop$weights["y"]), 0.7)
  # ds_params recorded for anchor reuse
  expect_identical(prop$ds_params$y$type, "max")
  expect_identical(prop$ds_params$y2$type, "min")
  # proposals within sampled numeric ranges of the training data schema
  expect_true(all(prop$proposals$X1 >= min(fits$data$X1) - 1e-8))
  expect_true(all(prop$proposals$X1 <= max(fits$data$X1) + 1e-8))
})

test_that("draws are reproducible under set.seed", {
  fits <- fit_two_models()
  set.seed(7)
  p1 <- svem_thompson_batch(fits$objects, basic_goals,
                            batch_size = 3, n = 300, verbose = FALSE)
  set.seed(7)
  p2 <- svem_thompson_batch(fits$objects, basic_goals,
                            batch_size = 3, n = 300, verbose = FALSE)
  expect_identical(p1$members, p2$members)
  expect_equal(p1$proposals, p2$proposals)
})

test_that("user-supplied candidates are used verbatim", {
  fits <- fit_two_models()
  cand <- fits$data[, c("X1", "X2", "X3")]
  set.seed(5)
  prop <- svem_thompson_batch(fits$objects, basic_goals,
                              batch_size = 3, candidates = cand,
                              verbose = FALSE)
  expect_equal(nrow(prop$candidates), nrow(cand))
  # each proposal is an actual candidate row
  for (b in seq_len(3)) {
    i <- prop$candidate_index[b]
    expect_equal(unlist(prop$proposals[b, c("X1", "X2", "X3")]),
                 unlist(cand[i, ]), tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
})

test_that("mixture constraints are honored in generated candidates", {
  df <- gen_toy_df(n = 50L)
  set.seed(42)
  f1 <- SVEMnet(y ~ A + B + C + X1, df, nBoot = 20)
  mix <- list(list(vars = c("A", "B", "C"),
                   lower = c(0.1, 0.1, 0.1),
                   upper = c(0.8, 0.8, 0.8),
                   total = 1.0))
  set.seed(6)
  prop <- svem_thompson_batch(list(y = f1),
                              list(y = list(goal = "max", weight = 1)),
                              batch_size = 4, n = 300,
                              mixture_groups = mix, verbose = FALSE)
  sums <- prop$proposals$A + prop$proposals$B + prop$proposals$C
  expect_true(all(abs(sums - 1) < 1e-6))
  expect_true(all(prop$proposals$A >= 0.1 - 1e-8 & prop$proposals$A <= 0.8 + 1e-8))
})

test_that("factor predictors and target goals are supported", {
  fits <- fit_two_models(with_factor = TRUE)
  goals <- list(
    y  = list(goal = "target", weight = 0.5, target = 4, tol = 1),
    y2 = list(goal = "max", weight = 0.5,
              lower_acceptable = 1, upper_acceptable = 4)
  )
  set.seed(9)
  prop <- svem_thompson_batch(fits$objects, goals,
                              batch_size = 3, n = 300, verbose = FALSE)
  expect_true(all(as.character(prop$proposals$F) %in% c("lo", "hi")))
  expect_identical(prop$ds_params$y$type, "target")
  expect_equal(prop$ds_params$y$T0, 4)
  # explicit anchors are used verbatim
  expect_equal(prop$ds_params$y2$L, 1)
  expect_equal(prop$ds_params$y2$U, 4)
})

test_that("binomial responses are handled on the probability scale", {
  df <- gen_toy_df(n = 60L)
  set.seed(12)
  df$yb <- rbinom(nrow(df), 1, plogis(scale(df$y)))
  set.seed(42)
  fb <- SVEMnet(yb ~ X1 + X2 + X3, df, nBoot = 20, family = "binomial")
  set.seed(13)
  prop <- svem_thompson_batch(list(yb = fb),
                              list(yb = list(goal = "max", weight = 1,
                                             lower_acceptable = 0,
                                             upper_acceptable = 1)),
                              batch_size = 2, n = 200, verbose = FALSE)
  expect_true(all(prop$proposals$yb_draw >= 0 & prop$proposals$yb_draw <= 1))
  expect_true(all(prop$proposals$yb_pred >= 0 & prop$proposals$yb_pred <= 1))
})

test_that("input validation errors are clear", {
  fits <- fit_two_models()

  # missing weight
  expect_error(
    svem_thompson_batch(fits$objects,
                        list(y = list(goal = "max"),
                             y2 = list(goal = "min", weight = 1)),
                        batch_size = 2, n = 100, verbose = FALSE),
    "goal and weight"
  )
  # bad goal keyword
  expect_error(
    svem_thompson_batch(fits$objects,
                        list(y = list(goal = "biggest", weight = 1),
                             y2 = list(goal = "min", weight = 1)),
                        batch_size = 2, n = 100, verbose = FALSE),
    "'max', 'min', or 'target'"
  )
  # batch larger than candidate table
  cand <- fits$data[1:3, c("X1", "X2", "X3")]
  expect_error(
    svem_thompson_batch(fits$objects, basic_goals,
                        batch_size = 5, candidates = cand, verbose = FALSE),
    "exceeds the number of candidates"
  )
  # model without bootstrap members
  no_boot <- fits$objects
  no_boot$y$coef_matrix <- NULL
  expect_error(
    svem_thompson_batch(no_boot, basic_goals,
                        batch_size = 2, n = 100, verbose = FALSE),
    "coef_matrix"
  )
})

test_that("print method runs quietly and returns invisibly", {
  fits <- fit_two_models()
  set.seed(3)
  prop <- svem_thompson_batch(fits$objects, basic_goals,
                              batch_size = 2, n = 200, verbose = FALSE)
  expect_output(out <- withVisible(print(prop)), "Thompson-sampling batch")
  expect_false(out$visible)
})
