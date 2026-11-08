test_that("primary fit controls cannot dispatch the spline engine", {
  dat <- data.frame(x = seq(-1, 1, length.out = 8))
  dat$y <- 1 + dat$x
  spec <- bigexp_terms(y ~ x, dat, report = FALSE)
  goals <- list(y = list(goal = "max", weight = 1))
  expect_error(svem_plan_runs(spec, dat, goals, dat["x"],
    budget = 1, incumbent = FALSE, fit_control = list(method = "spline")),
    "fit_control cannot replace: method")
  expect_error(svem_plan_runs(spec, dat, goals, dat["x"],
    budget = 1, incumbent = FALSE, fit_control = list(domain_grid = dat["x"])),
    "fit_control cannot replace: domain_grid")
})


test_that("a complete model reuse request explains the additional spline role", {
  dat <- data.frame(x = seq(-1, 1, length.out = 8))
  dat$y <- 1 + dat$x
  spec <- bigexp_terms(y ~ x, dat, report = FALSE)
  expect_error(svem_plan_runs(spec, dat, list(y = list(goal = "max", weight = 1)),
    dat["x"], budget = 6, fits = list(), forward = list(), n_optima = 5),
    "Pass spline = previous_plan")
})

test_that("direct spline scoring respects continuous companion block declarations", {
  spline <- list(y = structure(list(), class = "svem_spline"))
  primary <- list(y = list(companions = structure(list(),
    input_map = list(continuous_blocks = "day"))))
  expect_error(.svem_spline_check_schema(spline, primary), "continuous companion block")
})

test_that("automatic spline mixtures explain the shared specification before fitting", {
  dat <- data.frame(a = seq(0, 1, length.out = 11))
  dat$b <- 1 - dat$a
  dat$y <- dat$a
  spec <- bigexp_terms(y ~ a + b, dat, report = FALSE)
  events <- character()
  local_mocked_bindings(svem_fit_responses = function(...) {
    events <<- c(events, "primary")
    stop("primary boundary")
  }, .package = "SVEMnet")
  expect_error(withCallingHandlers(svem_plan_runs(spec, dat,
    list(y = list(goal = "max", weight = 1)), dat[c("a", "b")],
    budget = 3, spline = TRUE,
    fit_control = list(mixture_groups = list(list(vars = c("a", "b"), total = 1)))),
    warning = function(w) {
      if (grepl("Declare mixture_groups", conditionMessage(w), fixed = TRUE)) {
        events <<- c(events, "mixture guidance")
        invokeRestart("muffleWarning")
      }
    }), "primary boundary")
  expect_identical(events, c("mixture guidance", "primary"))
})
