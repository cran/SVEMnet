test_that("random prediction tables respect fixed and boundary mixture components", {
  d <- expand.grid(A = c(.1, .2, .4), B = c(.1, .3, .5), x = c(-1, 1))
  d$C <- 1 - d$A - d$B
  d$y <- 1 + 2 * d$A - d$B + d$x + sin(seq_len(nrow(d))) / 10
  set.seed(906)
  fit <- SVEMnet(y ~ A + B + C + x, d, nBoot = 4, glmnet_alpha = 1)
  sample_grid <- function(lower, upper, total = 1) {
    svem_random_table_multi(
      fit, n = 12, numeric_sampler = "uniform",
      mixture_groups = list(list(vars = c("A", "B", "C"),
                                 lower = lower, upper = upper, total = total))
    )
  }

  fixed <- sample_grid(c(.2, 0, 0), c(.2, 1, 1))
  expect_equal(fixed$data$A, rep(.2, 12))
  expect_equal(rowSums(fixed$data[c("A", "B", "C")]), rep(1, 12))
  expect_true(all(as.matrix(fixed$data[c("B", "C")]) >= 0))
  expect_true(all(is.finite(as.matrix(fixed$pred))))

  upper <- sample_grid(c(0, 0, 0), c(.2, .3, .5))
  lower <- sample_grid(c(.2, .3, .5), c(1, 1, 1))
  all_fixed <- sample_grid(c(.2, .3, .5), c(.2, .3, .5))
  one_free <- sample_grid(c(.2, .3, 0), c(.2, .3, 1))
  expected <- data.frame(A = rep(.2, 12), B = rep(.3, 12), C = rep(.5, 12))
  for (out in list(upper, lower, all_fixed, one_free)) {
    expect_equal(out$data[c("A", "B", "C")], expected)
    expect_true(all(is.finite(as.matrix(out$pred))))
  }

  expect_error(sample_grid(c(.8, .1, .1), c(.5, .9, .9)), "lower <= upper")
  expect_error(sample_grid(c(NA, 0, 0), c(1, 1, 1)), "finite numeric vectors")
  expect_error(sample_grid(c(0, 0, 0), c(1, Inf, 1)), "finite numeric vectors")
  expect_error(sample_grid(c(0, 0, 0), c(1, 1, 1), c(1, 2)), "single finite number")
  expect_error(sample_grid(c(.5, .5, .5), c(1, 1, 1)), "Infeasible mixture")
})

test_that("whole-model evaluation grids respect fixed and boundary mixture components", {
  d <- expand.grid(A = c(.1, .2, .4), B = c(.1, .3, .5), x = c(-1, 1))
  d$C <- 1 - d$A - d$B
  d$y <- 1 + 2 * d$A - d$B + d$x
  captured <- NULL
  # Exercise the public validation and grid construction, stopping before
  # costly SVEM refits; this regression concerns only the evaluation grid.
  wmt <- svem_significance_test_parallel
  wmt_env <- new.env(parent = environment(wmt))
  wmt_env$.svem_wmt_lapply <- function(..., T_data) {
    captured <<- T_data
    stop("evaluation grid captured")
  }
  environment(wmt) <- wmt_env
  sample_grid <- function(lower, upper, total = 1) {
    wmt(
      y ~ A + B + C + x, d,
      mixture_groups = list(list(vars = c("A", "B", "C"),
                                 lower = lower, upper = upper, total = total)),
      nPoint = 12, nSVEM = 1, nPerm = 20, nBoot = 2,
      nCore = 1, seed = 906, verbose = FALSE
    )
  }

  expect_error(sample_grid(c(.2, 0, 0), c(.2, 1, 1)), "evaluation grid captured")
  expect_equal(captured$A, rep(.2, 12))
  expect_equal(rowSums(captured[c("A", "B", "C")]), rep(1, 12))
  expect_true(all(as.matrix(captured[c("B", "C")]) >= 0))

  expected <- data.frame(A = rep(.2, 12), B = rep(.3, 12), C = rep(.5, 12))
  expect_error(sample_grid(c(0, 0, 0), c(.2, .3, .5)), "evaluation grid captured")
  expect_equal(captured[c("A", "B", "C")], expected)
  expect_error(sample_grid(c(.2, .3, .5), c(1, 1, 1)), "evaluation grid captured")
  expect_equal(captured[c("A", "B", "C")], expected)
  expect_error(sample_grid(c(.2, .3, .5), c(.2, .3, .5)), "evaluation grid captured")
  expect_equal(captured[c("A", "B", "C")], expected)
  expect_error(sample_grid(c(.2, .3, 0), c(.2, .3, 1)), "evaluation grid captured")
  expect_equal(captured[c("A", "B", "C")], expected)

  expect_error(sample_grid(c(.8, .1, .1), c(.5, .9, .9)), "lower <= upper")
  expect_error(sample_grid(c(NA, 0, 0), c(1, 1, 1)), "finite numeric vectors")
  expect_error(sample_grid(c(0, 0, 0), c(1, Inf, 1)), "finite numeric vectors")
  expect_error(sample_grid(c(0, 0, 0), c(1, 1, 1), NA_real_), "single finite number")
  expect_error(sample_grid(c(.5, .5, .5), c(1, 1, 1)), "Infeasible mixture")
})
