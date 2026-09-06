testthat::skip_on_cran()

# Regression invariants follow Karl (2024), Chemometrics and Intelligent
# Laboratory Systems 249, 105122 (doi:10.1016/j.chemolab.2024.105122), while
# the package-facing defaults and workflows follow Karl (2026), CILS 271,
# 105660 (doi:10.1016/j.chemolab.2026.105660).

test_that("WMT validates scalar controls without truncation", {
  d <- data.frame(y = seq_len(12), x = seq_len(12))
  call_wmt <- function(...) {
    do.call(
      svem_significance_test_parallel,
      utils::modifyList(
        list(
          formula = y ~ x,
          data = d,
          nPoint = 20,
          nSVEM = 1,
          nPerm = 20,
          nBoot = 2,
          nCore = 1,
          verbose = FALSE
        ),
        list(...)
      )
    )
  }

  expect_error(call_wmt(nPoint = 2.5), "nPoint")
  expect_error(call_wmt(nSVEM = c(1, 2)), "nSVEM")
  expect_error(call_wmt(nPerm = 4), "nPerm")
  expect_error(call_wmt(nBoot = 1), "nBoot")
  expect_error(call_wmt(nCore = 0), "nCore")
  expect_error(call_wmt(seed = -1), "seed")
  expect_error(call_wmt(percent = 100), "percent")
  expect_error(call_wmt(formula = y ~ x - 1), "intercept convention")
  expect_warning(
    expect_error(call_wmt(nPerm = 5, family = "binomial"), "binomial"),
    "nPerm < 20"
  )
})

test_that("WMT retries deterministically while holding a permutation fixed", {
  attempts <- 0L
  seen_y <- list()
  fit_fun <- function(formula, data, ...) {
    attempts <<- attempts + 1L
    seen_y[[attempts]] <<- data$y
    if (attempts < 3L) stop("deliberate retry")
    structure(list(), class = "wmt_stub")
  }
  predict_fun <- function(object, newdata, ...) {
    list(fit = seq_len(nrow(newdata)), se.fit = rep(2, nrow(newdata)))
  }
  d <- data.frame(y = seq_len(10), x = seq_len(10))

  ans <- SVEMnet:::.svem_wmt_fit_task(
    task = list(type = "permutation", index = 1L, seed = 123L),
    f_use = y ~ x,
    data_fit = d,
    resp_name = "y",
    nBoot = 2L,
    glmnet_alpha = 1,
    weight_scheme = "SVEM",
    objective = "wAIC",
    relaxed = FALSE,
    dots = list(),
    T_data = data.frame(x = 1:4),
    y_mean = mean(d$y),
    rng_sample_kind = "Rejection",
    fit_fun = fit_fun,
    predict_fun = predict_fun
  )

  expect_true(ans$ok)
  expect_identical(ans$retries, 2L)
  expect_length(ans$retry_failures, 2L)
  expect_identical(seen_y[[1L]], seen_y[[2L]])
  expect_identical(seen_y[[2L]], seen_y[[3L]])
  expect_false(identical(seen_y[[1L]], d$y))
})

test_that("WMT rejects degenerate member uncertainty and never drops slots", {
  fit_fun <- function(...) structure(list(), class = "wmt_stub")
  predict_fun <- function(object, newdata, ...) {
    list(fit = rep(1, nrow(newdata)), se.fit = rep(0, nrow(newdata)))
  }
  d <- data.frame(y = seq_len(8), x = seq_len(8))
  ans <- SVEMnet:::.svem_wmt_fit_task(
    task = list(type = "original", index = 1L, seed = 19L),
    f_use = y ~ x,
    data_fit = d,
    resp_name = "y",
    nBoot = 2L,
    glmnet_alpha = 1,
    weight_scheme = "SVEM",
    objective = "wAIC",
    relaxed = FALSE,
    dots = list(),
    T_data = data.frame(x = 1:3),
    y_mean = mean(d$y),
    rng_sample_kind = "Rejection",
    fit_fun = fit_fun,
    predict_fun = predict_fun
  )

  expect_false(ans$ok)
  expect_length(ans$retry_failures, 3L)
  expect_true(all(grepl("member SD", ans$retry_failures, fixed = TRUE)))
  expect_error(
    SVEMnet:::.svem_wmt_abort_failed_slots(list(ans), "original-data"),
    "no fits were dropped"
  )
})

test_that("WMT post-processing removes degenerate columns and uses strict retention", {
  M_pi <- cbind(
    c(-1, 1, 0, 0),
    c(0, 0, -1, 1),
    rep(5, 4)
  )
  M_Y <- matrix(c(0.5, -0.5, 5), nrow = 1L)

  expect_warning(
    reduced <- SVEMnet:::.svem_wmt_reduce_surfaces(M_Y, M_pi, percent = 50),
    "zero-variance"
  )
  expect_identical(reduced$removed_grid_columns, 3L)
  expect_identical(reduced$retained_components, 2L)
  expect_length(reduced$d_Y, 1L)
  expect_length(reduced$d_pi_Y, 4L)
  expect_identical(
    SVEMnet:::.svem_wmt_component_count(c(1, 1), percent = 50),
    2L
  )

  expect_warning(
    expect_error(
      SVEMnet:::.svem_wmt_reduce_surfaces(
        matrix(c(1, 1), nrow = 1L),
        matrix(1, nrow = 5L, ncol = 2L),
        percent = 90
      ),
      "All permutation-grid columns"
    ),
    "Removed 2"
  )
})

test_that("SHASHo right tails remain positive beyond CDF subtraction precision", {
  p <- SVEMnet:::.svem_wmt_shasho_upper_tail(
    q = 2, mu = 0, sigma = 1, nu = 0, tau = 2
  )
  transformed <- sinh(2 * asinh(2))
  expect_gt(p, 0)
  expect_equal(p, stats::pnorm(transformed, lower.tail = FALSE))
  expect_identical(
    1 - stats::pnorm(transformed),
    0
  )
})

test_that("legacy WMT seed schedule is stable", {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) .Random.seed else NULL
  on.exit({
    suppressWarnings(RNGkind(
      kind = old_kind[1L], normal.kind = old_kind[2L],
      sample.kind = old_kind[3L]
    ))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  legacy <- SVEMnet:::.svem_wmt_make_iter_seeds(
    5L, base_seed = 123L, sample_kind = "Rounding"
  )
  modern <- SVEMnet:::.svem_wmt_make_iter_seeds(
    5L, base_seed = 123L, sample_kind = "Rejection"
  )
  expect_identical(
    legacy,
    c(357285908L, 837187176L, 1631958398L, 1339477030L, 2000937474L)
  )
  expect_false(identical(modern, legacy))
})

test_that("private parallel WMT work leaves foreach registration unchanged", {
  skip_if_not_installed("foreach")
  before <- c(
    name = foreach::getDoParName(),
    workers = as.character(foreach::getDoParWorkers())
  )
  ans <- SVEMnet:::.svem_wmt_lapply(
    tasks = as.list(1:2),
    fun = identity,
    nCore = 2L,
    contrasts_opts = getOption("contrasts"),
    sample_kind = "Rejection"
  )
  after <- c(
    name = foreach::getDoParName(),
    workers = as.character(foreach::getDoParWorkers())
  )

  expect_identical(unlist(ans, use.names = FALSE), 1:2)
  expect_identical(after, before)
})

test_that("bounded WMT integration restores RNG and reports diagnostics", {
  set.seed(3400)
  n <- 36L
  d <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    x3 = rnorm(n)
  )
  d$y <- 1 + 2 * d$x1 - d$x2 + 0.5 * d$x3 + rnorm(n, sd = 0.35)

  set.seed(91)
  before_kind <- RNGkind()
  before_seed <- .Random.seed
  res <- svem_significance_test_parallel(
    y ~ x1 + x2 + x3,
    d,
    nPoint = 40,
    nSVEM = 1,
    nPerm = 20,
    nBoot = 8,
    nCore = 1,
    seed = 404,
    rng_sample_kind = "Rejection",
    verbose = FALSE
  )
  res_parallel <- svem_significance_test_parallel(
    y ~ x1 + x2 + x3,
    d,
    nPoint = 40,
    nSVEM = 1,
    nPerm = 20,
    nBoot = 8,
    nCore = 2,
    seed = 404,
    rng_sample_kind = "Rejection",
    verbose = FALSE
  )

  expect_s3_class(res, "svem_significance_test")
  expect_length(res$d_Y, 1L)
  expect_length(res$d_pi_Y, 20L)
  expect_true(is.finite(res$p_value))
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
  expect_identical(res$diagnostics$rng_sample_kind, "Rejection")
  expect_identical(res$diagnostics$retry_counts$total, 0L)
  expect_true(res$diagnostics$shasho_converged)
  expect_length(res$diagnostics$shasho_warnings, 0L)
  expect_identical(res_parallel$d_Y, res$d_Y)
  expect_identical(res_parallel$d_pi_Y, res$d_pi_Y)
  expect_identical(res_parallel$p_values, res$p_values)
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)
})

test_that("WMT restores RNG on an error exit", {
  d <- data.frame(y = 1:12, x = 1:12)
  set.seed(777)
  before_kind <- RNGkind()
  before_seed <- .Random.seed

  expect_error(
    svem_significance_test_parallel(
      y ~ x, d,
      mixture_groups = list(list(vars = "not_in_model")),
      nPoint = 20, nSVEM = 1, nPerm = 20, nBoot = 2,
      nCore = 1, seed = 9, verbose = FALSE
    ),
    "not in model predictors"
  )
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)
})

test_that("WMT evaluation grid reproduces raw predictor classes (logical, ordered)", {
  skip_on_ci()
  set.seed(23)
  n <- 40L
  d <- data.frame(
    x1   = runif(n),
    on   = sample(c(TRUE, FALSE), n, TRUE),
    dose = factor(sample(c("low", "mid", "high"), n, TRUE),
                  levels = c("low", "mid", "high"), ordered = TRUE)
  )
  d$y <- 2 * d$x1 + 0.7 * d$on + 0.5 * as.integer(d$dose) + rnorm(n, 0, 0.3)

  # Both the logical and the ordered predictor previously failed the
  # predict-time newdata class validator: the grid sampler rebuilt them as
  # plain factors, so every original-data fit slot aborted.
  res <- svem_significance_test_parallel(
    y ~ x1 + on + dose, d,
    nPoint = 60, nSVEM = 2, nPerm = 20, nBoot = 8,
    nCore = safe_ncores(), seed = 6, verbose = FALSE
  )
  expect_s3_class(res, "svem_significance_test")
  expect_true(is.finite(res$p_value))
})
