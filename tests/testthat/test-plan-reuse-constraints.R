.plan_reuse_goal <- list(y = list(goal = "max", weight = 1))

.plan_reuse_call <- function(spec, data, pool, ...) {
  svem_plan_runs(spec, data, .plan_reuse_goal, pool, budget = 1L,
                 score_control = list(verbose = FALSE), ...)
}

test_that("reused primary and forward models respect the declared blocking roles", {
  d <- data.frame(x = rep(seq(0, 1, length.out = 12), 2),
                  day = rep(c("A", "B"), each = 12))
  d$y <- d$x + ifelse(d$day == "A", 10, 0)
  spec <- bigexp_terms(y ~ x, d, blocking = "day", polynomial_order = 1,
                       report = FALSE)
  pool <- data.frame(x = seq(0, 1, length.out = 17), day = "B")
  set.seed(71)
  without_day <- list(y = SVEMnet(y ~ x, d, nBoot = 3, relaxed = FALSE))
  unblocked_day <- list(y = svem_forward(y ~ x + day, d, nBoot = 3))
  before <- .Random.seed
  expect_error(.plan_reuse_call(spec, d, pool, fits = without_day),
               "Cannot reuse supplied fits.*blocking roles")
  expect_error(.plan_reuse_call(spec, d, pool, forward = unblocked_day),
               "Cannot reuse supplied forward.*blocking roles")
  expect_identical(.Random.seed, before)

  plan <- .plan_reuse_call(spec, d, pool,
    fit_control = list(nBoot = 3, relaxed = FALSE, companions = FALSE),
    forward_control = list(nBoot = 3))
  reused <- .plan_reuse_call(spec, d, pool, fits = plan$fits, forward = plan$forward)
  expect_identical(reused$incumbent$source_data$day, "B")
  expect_identical(as.character(reused$worklist$day), "B")
  expect_identical(reused$worklist, plan$worklist)
})

test_that("mixture constraints cannot be dropped when reusing fitted models", {
  d <- data.frame(a = seq(.1, .9, length.out = 12))
  d$b <- 1 - d$a
  d$y <- d$a
  wide <- list(list(vars = c("a", "b")))
  tight <- list(list(vars = c("a", "b"), lower = c(.4, .4), upper = c(.6, .6)))
  old <- bigexp_terms(y ~ a + b, d, mixture_groups = wide, report = FALSE)
  spec <- bigexp_terms(y ~ a + b, d, mixture_groups = tight, report = FALSE)
  pool <- data.frame(a = seq(.4, .6, length.out = 7))
  pool$b <- 1 - pool$a
  set.seed(72)
  old_primary <- list(y = SVEMnet(old, d, nBoot = 3, relaxed = FALSE))
  old_forward <- list(y = svem_forward(old, d, nBoot = 3))
  before <- .Random.seed
  expect_error(.plan_reuse_call(spec, d, pool, fits = old_primary),
               "Cannot reuse supplied fits.*mixture constraints")
  expect_error(.plan_reuse_call(spec, d, pool, forward = old_forward),
               "Cannot reuse supplied forward.*mixture constraints")
  expect_error(svem_plan_runs(spec, d, .plan_reuse_goal, pool, budget = 1L,
    score_control = list(mixture_groups = wide)),
    "score_control mixture_groups must agree")
  expect_identical(.Random.seed, before)

  # A pure external fit with no fitted geometry can use scoring-only constraints.
  pure <- list(y = SVEMnet(y ~ a + b, d, nBoot = 3, relaxed = FALSE))
  plain_forward <- list(y = svem_forward(y ~ a + b, d, nBoot = 3))
  plan <- .plan_reuse_call(spec, d, pool, fits = pure, forward = plain_forward)
  expect_true(all(unlist(plan$worklist[c("a", "b")]) >= .4))
  expect_true(all(unlist(plan$worklist[c("a", "b")]) <= .6))
  expect_equal(plan$worklist$a, max(d$a[d$a >= .4 & d$a <= .6]))
  expect_identical(plan$fits, pure)
  expect_identical(plan$forward, plain_forward)

  # Existing geometry supplied only to the primary fit still constrains reuse.
  plain <- bigexp_terms(y ~ a + b, d, report = FALSE)
  fitted <- .plan_reuse_call(plain, d, pool,
    fit_control = list(nBoot = 3, relaxed = FALSE, companions = FALSE,
                       mixture_groups = tight), forward_control = list(nBoot = 3))
  reused <- .plan_reuse_call(plain, d, pool, fits = fitted$fits, forward = fitted$forward)
  expect_identical(reused$worklist, fitted$worklist)
  expect_equal(reused$worklist$a, plan$worklist$a)
})

test_that("discrete deployment support is checked independently of training rows", {
  d <- data.frame(x = rep(1:3, each = 6),
                  y = rep(c(1, 10, 3), each = 6) + sin(1:18) / 10)
  spec <- bigexp_terms(y ~ x, d, discrete_numeric = list(x = c(1, 3)),
                       polynomial_order = 1, report = FALSE)
  pool <- data.frame(x = c(1, 3))
  set.seed(73)
  continuous <- list(y = SVEMnet(y ~ x, d, nBoot = 3, relaxed = FALSE))
  before <- .Random.seed
  expect_error(.plan_reuse_call(spec, d, pool, fits = continuous),
               "Cannot reuse supplied fits.*discrete support for 'x'")
  expect_identical(.Random.seed, before)

  plan <- .plan_reuse_call(spec, d, pool,
    fit_control = list(nBoot = 3, relaxed = FALSE, companions = FALSE),
    forward_control = list(nBoot = 3))
  before <- .Random.seed
  reused <- .plan_reuse_call(spec, d[nrow(d):1, ], pool, fits = plan$fits,
                             forward = plan$forward)
  expect_identical(reused$worklist$x, 3)
  expect_identical(reused$worklist$x, plan$worklist$x)
  expect_identical(.Random.seed, before)

  for (field in c("discrete_values", "discrete_numeric")) {
    as_legacy <- function(fit) {
      fit$sampling_schema[[field]] <- fit$sampling_schema$discrete_levels
      fit$sampling_schema$discrete_levels <- NULL
      fit
    }
    legacy <- .plan_reuse_call(spec, d, pool, fits = lapply(plan$fits, as_legacy),
                              forward = lapply(plan$forward, as_legacy))
    expect_identical(legacy$worklist, plan$worklist)
  }
  expect_identical(.Random.seed, before)

  changed <- d
  changed$x[1] <- 2
  expect_error(.plan_reuse_call(spec, changed, pool, fits = plan$fits),
               "Cannot reuse supplied fits.*outcomes or design settings differ")
  categorical <- SVEMnet(y ~ factor(x), d, nBoot = 3, relaxed = FALSE)
  expect_error(.svem_plan_training_rows(categorical, rbind(d, data.frame(x = 4, y = 3))),
               "New measured discrete predictor settings")
})

test_that("numeric labels for a declared nominal block stay nominal in the plan", {
  d <- expand.grid(x = seq(0, 1, length.out = 8), day = factor(c("1", "2")))
  d$y <- d$x + ifelse(d$day == "1", 10, 0)
  expect_warning(spec <- bigexp_terms(y ~ x, d, blocking = "day",
    polynomial_order = 1, report = FALSE), "treat it as .*categorical")
  d$day <- as.integer(as.character(d$day))
  pool <- data.frame(x = seq(0, 1, length.out = 17), day = 2)
  set.seed(74)
  plan <- .plan_reuse_call(spec, d, pool,
    fit_control = list(nBoot = 3, relaxed = FALSE, companions = FALSE),
    forward_control = list(nBoot = 3))
  expect_identical(as.character(plan$worklist$day), "2")
  expect_identical(as.character(plan$incumbent$source_data$day), "2")
  expect_identical(unname(plan$fits$y$sampling_schema$var_classes["day"]), "factor")
  reused <- .plan_reuse_call(spec, d, pool, fits = plan$fits, forward = plan$forward)
  expect_identical(reused$worklist, plan$worklist)
})
