# Exercise ordinary public entry points on the same measured responses and feasible recipes.
spline_plan_fixture <- function() {
  pool <- expand.grid(u = seq(0, 1, length.out = 7), v = seq(0, 1, length.out = 7))
  pool$recipe_id <- sprintf("P%02d", seq_len(nrow(pool)))
  index <- 1L + ((0:23 * 19L) %% nrow(pool))
  measured <- pool[index, ]
  measured$activity <- with(measured, 1 + u - .25 * v + .4 * u * v) + .03 * cos(index)
  measured$size <- with(measured, .4 + v + .2 * u^2) + .02 * sin(index)
  measured$size[c(5L, 9L)] <- NA_real_
  spec <- bigexp_terms(activity ~ u + v, measured, report = FALSE)
  goals <- list(activity = list(goal = "max", weight = 1), size = list(goal = "min", weight = 1))
  set.seed(7314)
  fits <- svem_fit_responses(spec, measured, names(goals), nBoot = 4, companions = FALSE)
  forward <- stats::setNames(lapply(names(goals), function(r)
    svem_forward(spec, measured, response = r, nBoot = 4)), names(goals))
  list(pool = pool, data = measured, spec = spec, goals = goals, fits = fits, forward = forward)
}

test_that("the optional spline and explicit exploration are distinct physical planner roles", {
  skip_if_not_installed("mgcv")
  f <- spline_plan_fixture()
  set.seed(7315)
  splines <- svem_fit_responses(f$spec, f$data, names(f$goals), method = "spline",
                               domain_grid = f$pool, nBoot = 4)
  expect_identical(names(splines), names(f$goals))
  expect_true(all(vapply(splines, inherits, logical(1), "svem_spline")))
  plan_call <- function(data, saved = splines)
    svem_plan_runs(f$spec, data, f$goals, f$pool, budget = 4L, fits = f$fits,
      forward = f$forward, spline = saved, n_optima = 1L, n_explore = 2L,
      priority = c("incumbent", "spline", "svemnet", "forward"),
      exploration_fraction = 1, score_control = list(verbose = FALSE))
  seed <- .Random.seed
  plan <- plan_call(f$data)
  expect_identical(.Random.seed, seed)
  work <- plan$worklist
  expect_equal(nrow(work), 4)
  expect_false(anyDuplicated(work$recipe_id) > 0)
  expect_equal(sum(work$role == "incumbent_repeat"), 1)
  expect_equal(sum(work$role == "optimum"), 1)
  expect_equal(sum(work$role == "high_spread_medoid"), 2)
  nomination <- subset(plan$nominations, engine == "spline")
  expect_equal(nrow(nomination), 1)
  expect_true(grepl("spline", work$nominated_by[nomination$slot], fixed = TRUE))
  for (r in names(f$goals)) {
    col <- paste0(r, "_spline_pred")
    expect_equal(work[[col]], as.numeric(predict(splines[[r]], newdata = work)))
  }
  export <- svem_export_candidates_csv(plan, write_file = FALSE)
  expect_identical(export$recipe_id, work$recipe_id)
  expect_true(all(c("activity_spline_pred", "size_spline_pred", "spline_score") %in% names(export)))

  # Existing package checks accept row permutations and per-endpoint complete cases.
  again <- plan_call(f$data[nrow(f$data):1, ])
  expect_identical(again$worklist$recipe_id, work$recipe_id)
  expect_identical(again$worklist$spline_score, work$spline_score)
  expect_identical(.Random.seed, seed)

  changed <- f$data
  changed$activity[1] <- changed$activity[1] + .25
  stale <- splines
  stale$activity <- svem_fit_responses(f$spec, changed, "activity", method = "spline",
                                      domain_grid = f$pool, nBoot = 3)[[1]]
  seed <- .Random.seed
  expect_error(plan_call(f$data, stale), "training|measured|refit|stale|supplied")
  expect_identical(.Random.seed, seed)
})

test_that("an omitted spline and omitted counts preserve the ordinary worklist", {
  f <- spline_plan_fixture()
  args <- list(spec = f$spec, data = f$data, goals = f$goals, candidates = f$pool,
               budget = 3L, fits = f$fits, forward = f$forward,
               score_control = list(verbose = FALSE))
  seed <- .Random.seed
  direct <- do.call(svem_plan_runs, args)
  explicit <- do.call(svem_plan_runs, c(args, list(spline = NULL, n_optima = NULL, n_explore = NULL)))
  expect_identical(explicit$worklist, direct$worklist)
  expect_identical(explicit$nominations, direct$nominations)
  expect_false(any(grepl("spline", names(direct$worklist), fixed = TRUE)))
  expect_identical(.Random.seed, seed)
})

test_that("the planner fits a nominal study factor through the optional general spline", {
  skip_if_not_installed("mgcv")
  pool <- expand.grid(rate = seq(0, 1, length.out = 13),
                      ingredient = factor(c("A", "B")))
  pool$recipe_id <- sprintf("C%02d", seq_len(nrow(pool)))
  d <- pool[c(1:7, 13:20, 26), ]
  d$response <- with(d, 2 * rate + .3 * (ingredient == "B")) +
    .03 * cos(seq_len(nrow(d)))
  spec <- bigexp_terms(response ~ rate + ingredient, d, report = FALSE)
  set.seed(7316)
  plan <- suppressWarnings(svem_plan_runs(spec, d, list(response = list(goal = "max", weight = 1)),
    pool, budget = 2, spline = TRUE, n_optima = 1, n_explore = 0,
    priority = c("incumbent", "spline", "svemnet"),
    fit_control = list(nBoot = 3, companions = FALSE), forward_control = list(nBoot = 3),
    score_control = list(verbose = FALSE)))
  row <- subset(plan$nominations, engine == "spline")
  expect_equal(nrow(row), 1)
  expect_equal(row$slot,2L)
  expect_identical(row$fit_status, "ok")
  expect_equal(plan$worklist$role, c("incumbent_repeat", "optimum"))
  expect_match(plan$worklist$nominated_by[2], "spline")
})

