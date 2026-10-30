.capture_stabilization_wmt_grid <- function(formula, data, ...) {
  captured <- NULL
  wmt <- svem_significance_test_parallel
  local_env <- new.env(parent = environment(wmt))
  local_env$.svem_wmt_lapply <- function(..., T_data, f_use, y_mean, data_fit) {
    captured <<- list(grid = T_data, formula = f_use, mean = y_mean,
                      data = data_fit)
    stop("grid captured")
  }
  environment(wmt) <- local_env
  expect_error(wmt(formula, data, nPoint = 24, nSVEM = 1, nPerm = 20,
                   nBoot = 2, seed = 484, verbose = FALSE, ...), "grid captured")
  captured
}

test_that("WMT honors attached specs and preserves transformed responses", {
  d <- expand.grid(x = seq(-1, 1, length.out = 6), dose = c(1, 2, 4),
                   block = factor(c("a", "b")))
  d$y <- exp(2 + d$x + .1 * d$dose + (d$block == "b"))
  spec <- bigexp_terms(y ~ x + dose, d, blocking = "block",
                       discrete_numeric = "dose", report = FALSE)
  attached <- .capture_stabilization_wmt_grid(spec$formula, d)
  explicit <- .capture_stabilization_wmt_grid(y ~ x, d, spec = spec)
  expect_identical(attached$grid, explicit$grid)
  expect_true(all(attached$grid$dose %in% c(1, 2, 4)))
  expect_length(unique(attached$grid$block), 1L)

  f <- spec$formula
  f[[2L]] <- quote(log(y))
  transformed <- .capture_stabilization_wmt_grid(f, d)
  explicit_transformed <- .capture_stabilization_wmt_grid(log(y) ~ x, d, spec = spec)
  expect_identical(transformed$formula[[2L]], quote(log(y)))
  expect_identical(explicit_transformed$formula[[2L]], quote(log(y)))
  expect_equal(transformed$mean, mean(log(d$y)))
  expect_equal(explicit_transformed$mean, mean(log(d$y)))
  expect_identical(transformed$grid, explicit_transformed$grid)
})

test_that("WMT grids exclude formula constants and respect direct factor supports", {
  d <- expand.grid(x = seq(-1, 1, length.out = 6), dose = c(1, 2, 4))
  d$y <- 1 + d$x + d$dose
  center <- .25
  captured <- .capture_stabilization_wmt_grid(y ~ I(x - center) + factor(dose), d)
  expect_identical(names(captured$grid), c("x", "dose"))
  expect_true(all(captured$grid$dose %in% c(1, 2, 4)))
  expect_error(svem_significance_test_parallel(y ~ x + offset(dose), d,
               nPerm = 20, nPoint = 4, nSVEM = 1, nBoot = 2, verbose = FALSE),
               "offset")
})

test_that("failed whole-model responses warn even when progress is quiet", {
  d <- data.frame(x = 1:6, y = 2:7)
  expect_warning(out <- svem_wmt_multi(
    list(missing = absent ~ x), d, plot = FALSE, verbose = FALSE,
    wmt_control = list(nPerm = 20, nPoint = 4, nSVEM = 1, nBoot = 2)
  ), "failed for response 'missing'.*multiplier 1")
  expect_named(out$wmt_objects, "missing")
  expect_null(out$wmt_objects$missing)
  expect_true(is.na(out$p_values[["missing"]]))
  expect_equal(out$multipliers[["missing"]], 1)
})

.capture_stabilization_wmt_multi <- function(formulas, data, control = list()) {
  captured <- list()
  notices <- character()
  wmt <- svem_significance_test_parallel
  local_env <- new.env(parent = environment(wmt))
  local_env$.svem_wmt_lapply <- function(..., T_data, f_use, resp_name, y_mean) {
    captured[[length(captured) + 1L]] <<- list(
      grid = T_data, formula = f_use, response = resp_name, mean = y_mean)
    stop("dispatch captured", call. = FALSE)
  }
  environment(wmt) <- local_env
  wrapper <- svem_wmt_multi
  wrapper_env <- new.env(parent = environment(wrapper))
  wrapper_env$svem_significance_test_parallel <- wmt
  environment(wrapper) <- wrapper_env
  out <- withCallingHandlers(wrapper(
    formulas, data, plot = FALSE, verbose = FALSE,
    wmt_control = c(control, list(nPerm = 20, nPoint = 8, nSVEM = 1,
                                  nBoot = 2, verbose = FALSE, seed = 484))
  ), warning = function(w) {
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(out = out, captured = captured, notices = notices)
}

test_that("WMT wrapper keeps formula responses with attached or shared specs", {
  d <- data.frame(x = seq(-1, 1, length.out = 24), z = rep(c(0, 1), 12))
  d$y1 <- exp(1 + d$x + .1 * d$z)
  d$y2 <- exp(2 - d$x - .2 * d$z)
  spec <- bigexp_terms(y1 ~ x + z, d, factorial_order = 1,
                       polynomial_order = 1, report = FALSE)
  forms <- list(y1 = bigexp_formula(spec, "y1"),
                y2 = bigexp_formula(spec, "y2"))
  cases <- list(
    .capture_stabilization_wmt_multi(forms, d, list(response = "y1")),
    .capture_stabilization_wmt_multi(list(y1 = y1 ~ x, y2 = y2 ~ x), d,
                                     list(spec = spec, response = "y1"))
  )
  for (got in cases) {
    expect_named(got$out$p_values, c("y1", "y2"))
    expect_identical(vapply(got$captured, `[[`, "", "response"), c("y1", "y2"))
    expect_identical(lapply(got$captured, function(x) x$formula[[2L]]),
                     list(quote(y1), quote(y2)))
    expect_equal(vapply(got$captured, `[[`, 0, "mean"), c(mean(d$y1), mean(d$y2)))
    expect_true(all(vapply(got$captured, function(x)
      identical(names(x$grid), c("x", "z")), logical(1))))
    expect_equal(sum(grepl("ignoring.*response", got$notices)), 1)
    expect_false("response" %in% names(got$out$wmt_control))
  }
  expect_identical(cases[[2]]$out$wmt_control$spec, spec)

  transformed <- forms[[1]]
  transformed[[2L]] <- quote(log(y1))
  single <- .capture_stabilization_wmt_multi(
    list(log_y1 = transformed), d, list(response = "y1"))
  expect_named(single$out$p_values, "log_y1")
  expect_identical(single$captured[[1]]$formula[[2L]], quote(log(y1)))
  expect_identical(single$captured[[1]]$response, "y1")
  expect_equal(single$captured[[1]]$mean, mean(log(d$y1)))
  expect_equal(sum(grepl("ignoring.*response", single$notices)), 1)

  duplicate <- .capture_stabilization_wmt_multi(
    forms, d, list(response = "y1", response = "y2",
                    formula = y2 ~ z, formula = y1 ~ z,
                    data = NULL, data = NULL,
                    mixture_groups = NULL, mixture_groups = NULL))
  expect_identical(vapply(duplicate$captured, `[[`, "", "response"), c("y1", "y2"))
  expect_false(any(c("response", "formula", "data", "mixture_groups") %in%
                     names(duplicate$out$wmt_control)))
  expect_equal(sum(grepl("ignoring the following entries", duplicate$notices)), 1)

  partial <- .capture_stabilization_wmt_multi(forms, d, list(resp = "y1"))
  expect_identical(vapply(partial$captured, `[[`, "", "response"), c("y1", "y2"))
  expect_true(any(grepl("Ignoring unrecognized argument.*resp", partial$notices)))

  ordinary <- .capture_stabilization_wmt_multi(list(y1 = y1 ~ x, y2 = y2 ~ x), d)
  expect_identical(vapply(ordinary$captured, `[[`, "", "response"), c("y1", "y2"))
  expect_equal(vapply(ordinary$captured, `[[`, 0, "mean"), c(mean(d$y1), mean(d$y2)))
  expect_true(all(vapply(ordinary$captured, function(x)
    identical(names(x$grid), "x"), logical(1))))
  expect_false(any(grepl("ignoring|unrecognized", ordinary$notices)))
})

test_that("WMT controls are validated before dispatch or RNG use", {
  withr::local_seed(919)
  d <- data.frame(x = 1:6, y = 2:7)
  calls <- 0L
  wrapper <- svem_wmt_multi
  local_env <- new.env(parent = environment(wrapper))
  local_env$svem_significance_test_parallel <- function(...) {
    calls <<- calls + 1L
    stats::runif(1)
    list(p_value = .5)
  }
  environment(wrapper) <- local_env
  malformed <- list(list(5), list(seed = 1, 5),
                    stats::setNames(list(1), NA_character_),
                    list(seed = 1, seed = 2))
  patterns <- c(rep("fully named list", 3), "duplicated entries: seed")
  for (i in seq_along(malformed)) {
    before <- .Random.seed
    expect_error(wrapper(y ~ x, d, wmt_control = malformed[[i]],
                         plot = FALSE, verbose = FALSE), patterns[[i]])
    expect_identical(.Random.seed, before)
    expect_identical(calls, 0L)
  }
  expect_no_warning(out <- wrapper(y ~ x, d, wmt_control = list(),
                                   plot = FALSE, verbose = FALSE))
  expect_identical(calls, 1L)
  expect_identical(out$wmt_control, list())
  expect_equal(unname(out$p_values), .5)
})
