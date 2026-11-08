test_that("legal limits do not rescale the reference geometry", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(a = seq(.2, .8, length.out = 7), b = seq(.35, .65, length.out = 7))
  narrow <- spline_terms(y ~ a + b, grid)
  bounds <- rbind(c(a = 0, b = -2), c(a = 1, b = 2))
  wide <- spline_terms(y ~ a + b, grid, bounds = bounds)
  expect_identical(wide$lower, narrow$lower)
  expect_identical(wide$widths, narrow$widths)
  expect_identical(wide$reference_X, narrow$reference_X)
  expect_identical(wide$penalty, narrow$penalty)
  expect_equal(wide$legal_lower, c(a = 0, b = -2))
  d <- transform(grid[seq(1, 49, by = 2), ], y = a + sin(4*b))
  set.seed(32); f <- do.call(svem_spline, list(spec = wide, data = d, nBoot = 3))
  set.seed(32); g <- svem_spline(narrow, d, nBoot = 3)
  expect_identical(f$coef_matrix, g$coef_matrix)
  expect_null(f[["call", exact = TRUE]])
  expect_null(f$spec$smooth[["X", exact = TRUE]])
  extra <- data.frame(a = c(.1, .9), b = c(-.1, .85))
  expect_true(all(.svem_spline_legal(wide, extra)))
  expect_false(any(.svem_spline_legal(narrow, extra)))
  expect_true(all(is.finite(predict(f, extra))))
  expect_equal(predict(unserialize(serialize(f, NULL)), extra), predict(f, extra))
  expect_error(predict(g, extra), "bounds")
  old <- narrow
  old$legal_lower <- old$legal_upper <- NULL
  expect_identical(.svem_spline_legal(old, rbind(grid, extra)), .svem_spline_legal(narrow, rbind(grid, extra)))
  expect_equal(.svem_spline_surface(old, grid[3:8, ]), narrow$reference_X[3:8, ], tolerance = 1e-10)
})

test_that("physical bounds and mixture closure use the candidate tolerance", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(A = seq(.1,.3,length.out=6), B = seq(.2,.4,length.out=6))
  grid$C <- 1 - grid$A - grid$B
  groups <- list(list(vars = names(grid), lower = c(0,0,0), upper = c(.4,1,1), total = 1))
  limits <- rbind(c(A=0,B=0,C=0), c(A=1,B=1,C=1))
  s <- spline_terms(y ~ A+B+C, grid, bounds = limits, mixture_groups = groups)
  a <- c(.4+5e-9, .4+2e-8, .2, .2, NA_real_)
  probe <- data.frame(A=a,B=.3,C=1-a-.3+c(0,0,4e-9,2e-8,0))
  expect_identical(.svem_spline_legal(s, probe), c(TRUE,FALSE,TRUE,FALSE,FALSE))
  expect_silent(.svem_check_candidate_mixtures(probe[c(1,3),], .svem_mixture_groups(groups), 1e-8))
  expect_equal(ncol(.svem_spline_coordinates(s, probe[c(1,3),])),2)
  expect_error(.svem_spline_coordinates(s,probe[2,,drop=FALSE]),"mixture")
})

test_that("whitened modes preserve the original curvature form", {
  skip_if_not_installed("mgcv")
  grid <- expand.grid(u = 0:8/8, v = 0:6/6)
  s <- spline_terms(y ~ u + v, grid)
  tmat <- s$transform
  raw_penalty <- s$smooth$S[[1L]]
  transformed <- crossprod(tmat, raw_penalty %*% tmat)
  values <- eigen((transformed+t(transformed))/2, symmetric=TRUE, only.values=TRUE)$values
  positive <- values[values > max(values)*sqrt(.Machine$double.eps)]
  lowest <- min(positive)
  # For these d=2 unpenalized polynomial directions, the independent raw
  # curvature is zero; the fitted smoother explicitly adds null_weight.
  curvature <- s$penalty
  curvature[length(curvature)-0:1] <- 0
  expect_length(positive, length(curvature)-2L)
  expect_equal(transformed/lowest, diag(curvature), tolerance=1e-7)
  expect_equal(crossprod(s$reference_X)/nrow(grid), diag(ncol(s$reference_X)), tolerance=1e-7)
  b <- cos(seq_along(curvature))
  original_coef <- drop(tmat %*% b)
  expect_equal(drop(crossprod(original_coef,raw_penalty %*% original_coef))/lowest,
               sum(curvature*b^2), tolerance=1e-7)
})

test_that("numeric block IDs retain exact values in scoring and reuse", {
  skip_if_not_installed("mgcv")
  grid <- data.frame(x=0:7/7)
  s <- spline_terms(y ~ x, grid, blocking="day")
  k <- 1 + 2*.Machine$double.eps
  ids <- rep(c(1,k,1e5), each=6)
  d <- data.frame(x=rep(grid$x,length.out=length(ids)), day=ids, y=rep(c(-1,3,1),each=6))
  f <- svem_spline(s,d,nBoot=3)
  expect_identical(f$block_levels$day,c(1,k,1e5))
  expect_equal(predict(f,transform(grid, day=100000L)), predict(f,transform(grid,day=1e5)))
  expect_error(predict(f,transform(grid,day=1+.Machine$double.eps)),"Unknown")
  expect_error(predict(f,transform(grid,day="100000")),"numeric spline block")
  wrong <- d
  wrong$day[1] <- k
  expect_error(.svem_spline_check_training(f,wrong),"differ")
  expect_silent(.svem_spline_check_training(f,d[nrow(d):1,]))
  expect_warning(svem_spline(s,transform(d,day=seq_len(nrow(d))),nBoot=1),
                 "unreplicated numeric block")
})

test_that("response-free knot subsampling preserves the caller's RNG state", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  grid <- data.frame(x=seq(0,1,length.out=2001))
  set.seed(6283)
  before <- .Random.seed
  s <- spline_terms(y ~ x, grid, rank=5)
  expect_identical(.Random.seed,before)
  expect_equal(s$rank,5)
  expect_equal(nrow(s$reference_X),2001)
})
