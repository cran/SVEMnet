test_that("shared spline fits and planning recognize numeric nominal labels", {
  measured <- expand.grid(x = seq(0, 1, length.out = 7), grade = factor(c("1", "2")))
  measured$y <- measured$x + as.numeric(measured$grade) / 4
  spec <- suppressWarnings(bigexp_terms(y ~ x + grade, measured, factorial_order = 1,
    polynomial_order = 1, report = FALSE))
  grid <- expand.grid(x = seq(.01, .99, length.out = 7), grade = factor(c("1", "2")))
  numeric_labels <- function(data) transform(data, grade = as.numeric(as.character(grade)))
  fit <- function(data, candidates) {
    set.seed(17)
    svem_fit_responses(spec, data, "y", method = "spline", domain_grid = candidates,
      rank = 2, nBoot = 2)$y
  }
  reference <- fit(measured, grid)
  expect_identical(reference$spec$smoother, "linear limit")
  expect_identical(predict(fit(measured, numeric_labels(grid)), grid), predict(reference, grid))
  expect_identical(predict(fit(numeric_labels(measured), grid), grid), predict(reference, grid))
  unknown <- numeric_labels(measured)
  unknown$grade[1L] <- 3
  expect_error(fit(unknown, grid), "nominal levels")

  plan <- function(candidates) {
    set.seed(17)
    svem_plan_runs(spec, measured, list(y = list(goal = "max", weight = 1)), candidates,
      budget = 1L, incumbent = FALSE, priority = "spline", spline = TRUE,
      spline_control = list(rank = 2, nBoot = 2),
      fit_control = list(nBoot = 2, glmnet_alpha = 1, relaxed = FALSE, companions = FALSE),
      forward_control = list(nBoot = 2), score_control = list(verbose = FALSE))
  }
  baseline <- plan(grid)
  imported <- plan(numeric_labels(grid))
  expect_s3_class(imported$spline$y, "svem_spline")
  expect_identical(imported$worklist, baseline$worklist)
})

test_that("the distinct spline reference site cap covers both basis versions", {
  grid <- data.frame(x = seq(0, 1, length.out = 65537L))
  expect_error(spline_terms(y ~ x, grid, rank = 3), "65536 distinct coordinate sites")
  expect_error(spline_terms(y ~ x + grade, transform(grid, grade = "a"), rank = 3),
    "65536 distinct coordinate sites")
})

test_that("mixture spline charts can use extreme physical scales", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(a = seq(0, 1, length.out = 12), grade = c("a", "b"))
  grid$b <- 1 - grid$a
  grid$y <- sin(3 * grid$a) + .1 * (grid$grade == "b")
  mixture <- function(total) list(list(vars = c("a", "b"), upper = rep(total, 2), total = total))
  for (formula in list(y ~ a + b, y ~ a + b + grade)) {
    s <- spline_terms(formula, grid, mixture_groups = mixture(1))
    set.seed(11)
    fit <- svem_spline(s, grid, nBoot = 3)
    for (total in c(1e-200, 1e200)) {
      scaled <- transform(grid, a = a * total, b = b * total)
      scaled_spec <- spline_terms(formula, scaled, mixture_groups = mixture(total))
      expect_identical(scaled_spec$version, s$version)
      expect_identical(scaled_spec$dimension, s$dimension)
      set.seed(11)
      scaled_fit <- svem_spline(scaled_spec, scaled, nBoot = 3)
      expect_equal(predict(scaled_fit, scaled), predict(fit, grid), tolerance = 1e-7)
    }
  }
})

test_that("flat general splines retain incomplete coverage diagnostics", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(x = 0:5 / 5, grade = factor(c("a", "b")))
  measured <- transform(grid, y = ifelse(grade == "a", 3, NA_real_))
  spec <- bigexp_terms(y ~ x + grade, measured, report = FALSE)
  primary <- svem_fit_responses(spec, measured, "y", nBoot = 2, companions = FALSE)
  smooth <- svem_fit_responses(spec, measured, "y", method = "spline", domain_grid = grid, nBoot = 2)
  pool <- transform(grid[grid$x < 1, ], x = x + .01)
  score <- function(candidates) suppressWarnings(svem_score_random(primary,
    list(y = list(goal = "max", weight = 1)), data = measured, candidates = candidates,
    spline = smooth, budget = 1L, allocation = "portfolio", priority = "spline",
    incumbent = FALSE, n_optima = 1L, n_explore = 0L, verbose = FALSE))
  partial <- score(pool)
  missing <- score(pool[pool$grade == "b", ])
  expect_identical(subset(partial$engine_status, engine == "spline")$status, "partial_coverage")
  expect_identical(subset(missing$engine_status, engine == "spline")$status, "no_finite_predictions")
  expect_match(subset(partial$engine_status, engine == "spline")$message, "5 recipe\\(s\\) have no prediction")
  expect_match(subset(missing$engine_status, engine == "spline")$message, "5 recipe\\(s\\) have no prediction")
  expect_equal(nrow(partial$worklist), 0L)
  expect_equal(nrow(missing$worklist), 0L)
})
