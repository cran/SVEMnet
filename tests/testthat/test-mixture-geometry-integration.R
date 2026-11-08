geometry_fixture <- function() {
  d <- expand.grid(A = c(.2, .4, .7), S = c(2, 4, 7), z = c(-1, 1))
  d$B <- 1 - d$A
  d$T <- 10 - d$S
  d$y <- 2 * d$z + d$S + .1 * sin(seq_len(nrow(d)))
  groups <- list(list(vars = c("A", "B")),
                 list(vars = c("S", "T"), lower = c(1, 1),
                      upper = c(9, 9), total = 10L))
  spec <- bigexp_terms(y ~ A + B + S + T + z, d, report = FALSE,
                       factorial_order = 1, polynomial_order = 1,
                       include_pc_2way = FALSE, mixture_groups = groups)
  list(data = d, groups = groups, spec = spec)
}

test_that("both forward fits inherit and check the declared recipe domain", {
  x <- geometry_fixture()
  set.seed(103)
  fits <- list(forward_aicc(x$spec, x$data),
               svem_forward(bigexp_formula(x$spec, "y"), x$data, nBoot = 4))
  for (fit in fits) {
    expect_identical(fit$sampling_schema$mixture_groups, x$spec$settings$mixture_groups)
    grid <- svem_random_table_multi(fit, n = 12, numeric_sampler = "uniform")$data
    expect_equal(grid$A + grid$B, rep(1, 12))
    expect_equal(grid$S + grid$T, rep(10, 12))
    expect_error(svem_random_table_multi(fit, n = 12, mixture_groups = x$groups[1]),
                 "must agree")
  }
  bad <- x$data
  bad$T[2] <- bad$T[2] - 1
  expect_error(forward_aicc(x$spec, bad), "sum to 10")
  expect_error(svem_forward(x$spec, bad, nBoot = 2), "sum to 10")
})

test_that("whole-model quadrature inherits both groups and fits only SVEM", {
  x <- geometry_fixture()
  captured <- NULL
  wmt <- svem_significance_test_parallel
  wmt_env <- new.env(parent = environment(wmt))
  wmt_env$.svem_wmt_lapply <- function(..., T_data, dots) {
    captured <<- list(grid = T_data, dots = dots)
    stop("captured grid")
  }
  environment(wmt) <- wmt_env
  run <- function(f, ...) {
    wmt(f, x$data, nPoint = 12, nSVEM = 1, nPerm = 20, nBoot = 2,
        nCore = 1, seed = 123, verbose = FALSE, ...)
  }
  f <- bigexp_formula(x$spec, "y")
  expect_error(run(f), "captured grid")
  inherited <- captured$grid
  expect_equal(inherited$A + inherited$B, rep(1, 12))
  expect_equal(inherited$S + inherited$T, rep(10, 12))
  expect_error(run(y ~ 1, spec = x$spec, mixture_groups = x$groups), "captured grid")
  expect_equal(captured$grid, inherited)
  expect_error(run(f, mixture_groups = x$groups[1]), "conflicts")
  expect_warning(expect_error(run(f, companions = TRUE, companion_blocks = "factor",
                                  companion_continuous_blocks = "Ambient"),
                              "captured grid"), "Ignoring companion options")
  expect_false(any(c("companions", "companion_blocks", "companion_continuous_blocks")
                   %in% names(captured$dots)))
})

test_that("supplied Thompson pools honor declared physical constraints", {
  x <- geometry_fixture()
  set.seed(103)
  fit <- svem_forward(x$spec, x$data, nBoot = 4)
  propose <- function(pool, ...) {
    svem_thompson_batch(fit, goals = list(y = list(goal = "max", weight = 1)),
                        batch_size = 1, candidates = pool, verbose = FALSE, ...)
  }
  pool <- x$data[setdiff(names(x$data), "y")]
  expect_s3_class(propose(pool), "svem_thompson_batch")
  expect_error(propose(pool, mixture_groups = x$groups[1]), "must agree")
  expect_error(propose(transform(pool, T = T + 1)), "bounds or total")
  expect_error(propose(transform(pool, S = 0, T = 10)), "bounds or total")
})

test_that("a CV fit can predict a prepared formula without impersonating a companion", {
  x <- geometry_fixture()
  fit <- glmnet_with_cv(bigexp_formula(x$spec, "y"), x$data, glmnet_alpha = 1,
                        nfolds = 3, repeats = 1, seed = 103)
  expect_equal(predict_cv(fit, x$data), predict(fit, x$data))
  expect_error(predict(fit, x$data, engine = "gaupro"), "without an engine")
})
