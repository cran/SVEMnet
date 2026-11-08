plan_fixture <- function(blocks = FALSE) {
  d <- data.frame(x = seq(-1, 1, length.out = 16))
  d$y <- 2 + d$x + .1 * sin(seq_len(nrow(d)))
  d$z <- 3 - d$x^2 + .1 * cos(seq_len(nrow(d)))
  if (blocks) {
    d <- rbind(transform(d, day = "A"), transform(d, day = "B", y = y + 100))
    d$day <- factor(d$day, levels = c("A", "B"))
    spec <- bigexp_terms(y ~ x, d, factorial_order = 1, polynomial_order = 1,
                         blocking = "day", report = FALSE)
    pool <- data.frame(x = seq(-1, 1, length.out = 41),
                       day = factor("A", levels = levels(d$day)))
  } else {
    spec <- bigexp_terms(y ~ x, d, factorial_order = 1, polynomial_order = 1,
                         report = FALSE)
    pool <- data.frame(x = seq(-1, 1, length.out = 41))
  }
  list(data = d, spec = spec, pool = pool,
       goals = list(y = list(goal = "max", weight = .6,
                             lower_acceptable = 0, upper_acceptable = 110),
                    z = list(goal = "min", weight = .4,
                             lower_acceptable = 1, upper_acceptable = 4)),
       primary = list(nBoot = 3, glmnet_alpha = 1, relaxed = FALSE, companions = FALSE),
       forward = list(nBoot = 3), score = list(verbose = FALSE, combine = "mean"))
}

test_that("the planner composes the ordinary fits and scorer in the same RNG order", {
  f <- plan_fixture()
  set.seed(21)
  plan <- svem_plan_runs(f$spec, f$data, f$goals, f$pool, budget = 5,
                         fit_control = f$primary, forward_control = f$forward,
                         score_control = f$score)
  seed <- .Random.seed
  set.seed(21)
  fits <- do.call(svem_fit_responses, c(list(spec = f$spec, data = f$data,
                                            responses = names(f$goals)), f$primary))
  forward <- setNames(lapply(names(f$goals), function(r) {
    do.call(svem_forward, c(list(formula = f$spec, data = f$data, response = r), f$forward))
  }), names(f$goals))
  scores <- do.call(svem_score_random, c(list(objects = fits, goals = f$goals,
    data = f$data, candidates = f$pool, budget = 5, forward = forward,
    allocation = "portfolio", incumbent = TRUE), f$score))
  expect_identical(.Random.seed, seed)
  expect_identical(names(plan$fits), names(f$goals))
  expect_identical(names(plan$forward), names(f$goals))
  for (r in names(f$goals)) {
    expect_identical(plan$fits[[r]]$coef_matrix, fits[[r]]$coef_matrix)
    expect_identical(plan$forward[[r]]$coef_matrix, forward[[r]]$coef_matrix)
  }
  expect_identical(plan$scores, scores)
  expect_identical(plan$worklist, scores$worklist)
  expect_equal(sum(!plan$worklist$new_recipe), 1)
  expect_equal(nrow(plan$worklist) + plan$shortfall, 5)

  before <- .Random.seed
  rescored <- svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = fits,
                             forward = forward, score_control = f$score)
  expect_identical(.Random.seed, before)
  expect_identical(rescored$worklist, plan$worklist)
  expect_output(print(plan), "5 of 5")
  expect_s3_class(plot(plan), "ggplot")
  expect_s3_class(plot(plan, type = "responses"), "ggplot")
  expect_equal(nrow(ggplot2::ggplot_build(plot(plan))$data[[2L]]), nrow(plan$worklist))
  expect_equal(nrow(ggplot2::ggplot_build(plot(plan, type = "responses"))$data[[1L]]),
               nrow(plan$worklist) * 2L * 2L)
})

test_that("the incumbent is selected by observed outcomes instead of fitted utility", {
  f <- plan_fixture()
  f$data$y[8] <- 8
  goals <- list(y = list(goal = "max", weight = 1,
                         lower_acceptable = 0, upper_acceptable = 10))
  set.seed(31)
  plan <- svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
                         fit_control = f$primary, forward_control = f$forward,
                         score_control = f$score)
  expect_identical(plan$incumbent$source_row, 8L)
  expect_equal(plan$incumbent$observed_responses$y, 8)
  expect_equal(plan$incumbent$observed_score, .8)
  expect_equal(plan$worklist$x, f$data$x[8])
  expect_identical(plan$worklist$role, "incumbent_repeat")
  expect_false(plan$worklist$new_recipe)
  expect_false(which.max(plan$scores$original_data_scored$y_pred) ==
                 plan$incumbent$source_row)
  expect_true(plan$worklist$recipe_id %in% plan$scores$score_table$recipe_id)

  bad <- f$data
  bad$y[8] <- NA_real_
  expect_error(svem_plan_runs(f$spec, bad, goals, f$pool, budget = 1,
    fits = plan$fits, forward = plan$forward, incumbent = 8L,
    score_control = f$score), "Cannot reuse supplied fits")
  expect_error(svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
    fits = plan$fits, forward = plan$forward, incumbent = 8L,
    score_control = c(f$score, list(increments = c(x = .05)))), "No eligible exact historical recipe")
  explicit <- svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
    fits = plan$fits, forward = plan$forward, incumbent = 1L, score_control = f$score)
  expect_identical(explicit$incumbent$source_row, 1L)
  expect_true(explicit$incumbent$explicit)

  empty_pool <- f$data["x"]
  expect_warning(empty <- svem_plan_runs(f$spec, f$data, goals, empty_pool,
    budget = 1, fits = plan$fits, forward = plan$forward, incumbent = FALSE,
    score_control = f$score), "short by 1")
  expect_equal(nrow(empty$worklist), 0L)
  expect_equal(empty$shortfall, 1)
  expect_s3_class(plot(empty, type = "responses"), "ggplot")
})

test_that("incumbent context is explicit across additive block shifts", {
  f <- plan_fixture(blocks = TRUE)
  goals <- list(y = list(goal = "max", weight = 1,
                         lower_acceptable = 0, upper_acceptable = 110))
  set.seed(32)
  base <- svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score)
  expected_a <- which(f$data$day == "A")[which.max(f$data$y[f$data$day == "A"])]
  expected_b <- which.max(f$data$y)
  expect_identical(base$incumbent$source_row, expected_a)
  chosen_b <- svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
    fits = base$fits, forward = base$forward, reference = list(day = "B"),
    score_control = f$score)
  expect_identical(chosen_b$incumbent$source_row, expected_b)
  expect_identical(as.character(chosen_b$incumbent$source_data$day), "B")
  expect_identical(as.character(chosen_b$worklist$day), "A")
  explicit <- svem_plan_runs(f$spec, f$data, goals, f$pool, budget = 1,
    fits = base$fits, forward = base$forward, incumbent = expected_b,
    score_control = f$score)
  expect_identical(explicit$incumbent$source_row, expected_b)
  expect_error(svem_plan_runs(f$spec, f$data, goals, f$pool, fits = base$fits,
    forward = base$forward, reference = list(day = "C"), score_control = f$score),
    "No historical preparations")
  expect_error(svem_plan_runs(f$spec, f$data, goals, f$pool, fits = base$fits,
    forward = base$forward, reference = list(day = "B"), incumbent = expected_b,
    score_control = f$score), "either a historical row")
})

test_that("forward desirability uses each scientific endpoint and shared anchors", {
  f <- plan_fixture()
  set.seed(33)
  plan <- svem_plan_runs(f$spec, f$data, f$goals, f$pool, budget = 5,
                         fit_control = f$primary, forward_control = f$forward,
                         score_control = f$score)
  tab <- plan$scores$score_table
  clip <- function(x) pmin(1, pmax(0, x))
  manual <- .6 * clip(tab$y_forward_pred / 110) +
    .4 * clip((4 - tab$z_forward_pred) / 3)
  expect_equal(tab$forward_score, manual, tolerance = 1e-14)
  expect_equal(plan$scores$weights_original, c(y = .6, z = .4))
  for (r in names(f$goals)) {
    expect_equal(tab[[paste0(r, "_forward_pred")]],
      as.numeric(predict(plan$forward[[r]], tab, debias = FALSE)), tolerance = 1e-14)
  }
  incomplete <- f$data
  best <- which.max(f$data$y)
  incomplete$z[best] <- Inf
  expect_error(svem_plan_runs(f$spec, incomplete, f$goals, f$pool,
    fits = plan$fits, forward = plan$forward, incumbent = best,
    score_control = f$score), "Cannot reuse supplied fits")
  no_weight <- f$goals
  no_weight$z$weight <- 0
  expect_error(svem_plan_runs(f$spec, incomplete, no_weight, f$pool,
    fits = plan$fits, forward = plan$forward, incumbent = best,
    score_control = f$score), "Cannot reuse supplied fits")
})

test_that("planner argument ownership and supplied response identities are validated", {
  f <- plan_fixture()
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    fit_control = list(data = f$data)), "cannot replace")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    forward_control = list(response = "y")), "cannot replace")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    score_control = list(budget = 1)), "cannot replace")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    fit_control = list(nBoot = 3, nBoot = 4)), "unique argument names")
  set.seed(34)
  base <- svem_plan_runs(f$spec, f$data, f$goals, f$pool, budget = 3,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score)
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = base$fits,
    forward = base$forward, fit_control = f$primary), "cannot modify supplied fits")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = base$fits,
    forward = base$forward, forward_control = f$forward), "cannot modify supplied fits")

  swapped <- base$forward[c("z", "y")]
  names(swapped) <- c("y", "z")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = base$fits,
    forward = swapped, score_control = f$score), "response|label")
  wrong_primary <- base$fits[c("z", "y")]
  names(wrong_primary) <- c("y", "z")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = wrong_primary,
    forward = base$forward, score_control = f$score), "response|label")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = base$fits,
    forward = base$fits, score_control = f$score), "forward|Gaussian")
})

test_that("unlabeled physical recipe identities survive candidate-pool reordering", {
  f <- plan_fixture()
  set.seed(37)
  base <- svem_plan_runs(f$spec, f$data, f$goals, f$pool, budget = 5,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score)
  reordered <- svem_plan_runs(f$spec, f$data, f$goals, f$pool[nrow(f$pool):1, , drop = FALSE],
    budget = 5, fits = base$fits, forward = base$forward, score_control = f$score)
  expect_identical(reordered$worklist$x, base$worklist$x)
  expect_identical(reordered$worklist$recipe_id, base$worklist$recipe_id)
  expect_identical(reordered$nominations$best_recipe_id, base$nominations$best_recipe_id)
  expect_identical(reordered$nominations$new_recipe_id, base$nominations$new_recipe_id)
})

.plan_augmented_day_fixture <- function() {
  f <- plan_fixture(blocks = TRUE)
  f$data$day <- as.character(f$data$day)
  added <- f$data[1:4, , drop = FALSE]
  added$day <- "C"
  added$y <- added$y + 50
  f$data <- rbind(f$data, added)
  f$data$day <- factor(f$data$day, levels = c("A", "B", "C"))
  f
}

test_that("a stale specification rejects new measured levels before fitting or RNG", {
  f <- .plan_augmented_day_fixture()
  fits_called <- 0L
  forbidden_fit <- function(...) {
    fits_called <<- fits_called + 1L
    stop("A fit was attempted before the measured-level check.")
  }
  local_mocked_bindings(svem_fit_responses = forbidden_fit,
                         svem_forward = forbidden_fit, .package = "SVEMnet")
  set.seed(831)
  before <- .Random.seed
  expect_true(all(is.finite(f$data$y[f$data$day == "C"])))
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score),
    "Unseen level.*day.*C.*rebuild.*augmented_data")
  expect_identical(fits_called, 0L)
  expect_identical(.Random.seed, before)

  # The same guard covers a newly measured ingredient, not just a block.
  old <- plan_fixture()
  old$data$ingredient <- rep(c("lipid-A", "lipid-B"), length.out = nrow(old$data))
  spec <- bigexp_terms(y ~ x + ingredient, old$data, report = FALSE)
  changed <- old$data
  changed$ingredient[1L] <- "lipid-C"
  expect_error(svem_plan_runs(spec, changed, old$goals, old$pool,
    fit_control = old$primary, forward_control = old$forward, score_control = old$score),
    "Unseen level.*ingredient.*lipid-C.*rebuild")
  expect_identical(fits_called, 0L)
  expect_identical(.Random.seed, before)
})

test_that("a refreshed specification retains every measured new-day preparation", {
  f <- .plan_augmented_day_fixture()
  refreshed <- bigexp_terms(y ~ x, f$data, factorial_order = 1,
                            polynomial_order = 1, blocking = "day", report = FALSE)
  set.seed(832)
  plan <- svem_plan_runs(refreshed, f$data, f$goals, f$pool, budget = 1L,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score)
  for (response in names(f$goals)) {
    expect_identical(plan$fits[[response]]$nobs, nrow(f$data))
    expect_identical(plan$forward[[response]]$nobs, nrow(f$data))
    expect_equal(as.numeric(plan$fits[[response]]$actual_y), f$data[[response]])
    expect_equal(as.numeric(plan$forward[[response]]$actual_y), f$data[[response]])
  }
  expect_equal(nrow(plan$worklist) + plan$shortfall, 1L)
})

test_that("the measured-level preflight permits endpoint-specific missing responses", {
  f <- .plan_augmented_day_fixture()
  f$data$y[2L] <- NA_real_
  f$data$z[3L] <- NA_real_
  refreshed <- bigexp_terms(y ~ x, f$data, factorial_order = 1,
                            polynomial_order = 1, blocking = "day", report = FALSE)
  set.seed(833)
  plan <- svem_plan_runs(refreshed, f$data, f$goals, f$pool, budget = 1L,
    fit_control = f$primary, forward_control = f$forward, score_control = f$score)
  for (response in names(f$goals)) {
    used <- f$data[[response]][!is.na(f$data[[response]])]
    expect_equal(plan$fits[[response]]$nobs, length(used))
    expect_equal(plan$forward[[response]]$nobs, length(used))
    expect_equal(as.numeric(plan$fits[[response]]$actual_y), used)
    expect_equal(as.numeric(plan$forward[[response]]$actual_y), used)
  }
  expect_true(all(is.finite(unlist(plan$incumbent$observed_responses))))
})
