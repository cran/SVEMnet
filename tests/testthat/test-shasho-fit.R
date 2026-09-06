test_that("SHASHo density and tails agree with frozen upstream reference values", {
  # Captured from gamlss.dist 6.1-1; no upstream package is needed at test time.
  ref <- read.csv(test_path("fixtures", "shasho-6.1-1.csv"))
  log_density <- with(ref, SVEMnet:::.svem_wmt_shasho_log_density(x, mu, sigma, nu, tau))
  tail <- with(ref, SVEMnet:::.svem_wmt_shasho_upper_tail(x, mu, sigma, nu, tau))
  expect_equal(log_density, ref$log_density, tolerance = 1e-12)
  expect_equal(tail, 1 - ref$cdf, tolerance = 1e-12)

  x <- seq(-5, 5, length.out = 101)
  expect_equal(SVEMnet:::.svem_wmt_shasho_log_density(x, 2, 3, 0, 1),
               dnorm(x, 2, 3, log = TRUE), tolerance = 1e-13)
  area <- integrate(function(x) exp(SVEMnet:::.svem_wmt_shasho_log_density(
    x, -.4, .7, .6, .8)), -Inf, Inf)$value
  expect_equal(area, 1, tolerance = 1e-7)
  extreme <- SVEMnet:::.svem_wmt_shasho_upper_tail(2, 0, 1, 0, 2)
  expect_gt(extreme, 0)
  expect_equal(extreme, pnorm(sinh(2 * asinh(2)), lower.tail = FALSE))
})

test_that("SHASHo analytic scores agree with independent finite differences", {
  x <- seq(-2, 3, length.out = 21)
  for (par in list(c(.3, -.4, .7, .2), c(-.4, .3, -.2, -.5))) {
    numeric_score <- vapply(seq_along(par), function(i) {
      up <- down <- par
      up[i] <- up[i] + 1e-6
      down[i] <- down[i] - 1e-6
      (SVEMnet:::.svem_wmt_shasho_objective(up, x) -
         SVEMnet:::.svem_wmt_shasho_objective(down, x)) / 2e-6
    }, numeric(1L))
    expect_equal(SVEMnet:::.svem_wmt_shasho_objective(par, x, TRUE),
                 numeric_score, tolerance = 1e-6)
  }
})

test_that("SHASHo fitting recovers a known distribution and preserves location-scale changes", {
  # Deterministic quantiles avoid random recovery failures in a short CRAN test.
  p <- (seq_len(400) - .5) / 400
  x <- 2 + 3 * sinh((asinh(qnorm(p)) + .6) / .8)
  set.seed(42)
  before <- .Random.seed
  fit <- SVEMnet:::.svem_wmt_fit_shasho(x)
  expect_identical(.Random.seed, before)
  expect_s3_class(fit, "svem_shasho_fit")
  expect_true(fit$converged)
  expect_lt(fit$gradient_max, 1e-4)
  fitted_cdf <- 1 - do.call(SVEMnet:::.svem_wmt_shasho_upper_tail,
                            c(list(q = x), as.list(fit$parameters)))
  expect_lt(max(abs(fitted_cdf - p)), .005)
  shifted <- SVEMnet:::.svem_wmt_fit_shasho(100 + 7 * x)
  expect_equal(unname(shifted$parameters["mu"]),
               100 + 7 * unname(fit$parameters["mu"]), tolerance = 1e-5)
  expect_equal(unname(shifted$parameters["sigma"]),
               7 * unname(fit$parameters["sigma"]), tolerance = 1e-5)
  expect_equal(shifted$parameters[c("nu", "tau")],
               fit$parameters[c("nu", "tau")], tolerance = 1e-5)
  expect_equal(shifted$loglik, fit$loglik - length(x) * log(7), tolerance = 1e-8)
  expect_equal(unname(coef(fit, what = "mu")), unname(fit$parameters["mu"]))
  expect_equal(unname(coef(fit, what = "sigma")), log(unname(fit$parameters["sigma"])))
  expect_equal(unname(coef(fit, what = "tau")), log(unname(fit$parameters["tau"])))
})

test_that("invalid SHASHo inputs fail explicitly", {
  for (x in list(1:4, rep(1, 20), c(1:6, NA), c(1:6, Inf))) {
    expect_error(SVEMnet:::.svem_wmt_fit_shasho(x), "five distinct finite")
  }
  expect_error(SVEMnet:::.svem_wmt_shasho_upper_tail(1, 0, 0, 0, 1),
               "positive sigma and tau")
})
