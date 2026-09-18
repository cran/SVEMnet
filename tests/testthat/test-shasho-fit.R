test_that("SHASHo density and tails follow the normal transformation", {
  # Jones and Pewsey (2009): inverse normal transformation and its Jacobian.
  # The accepted manuscript uses epsilon = -nu in the package convention.
  p <- c(.0001, .01, .1, .25, .5, .75, .9, .99, .9999)
  r <- qnorm(p)
  for (theta in list(c(-.4, .7, .6, .8), c(2, 3, -.8, 1.7), c(0, 1, 0, 1))) {
    mu <- theta[1]; sigma <- theta[2]; nu <- theta[3]; tau <- theta[4]
    w <- (asinh(r) + nu) / tau
    x <- mu + sigma * sinh(w)
    inverse_jacobian <- sigma * cosh(w) / (tau * sqrt(1 + r^2) * dnorm(r))
    expect_equal(SVEMnet:::.svem_wmt_shasho_log_density(x, mu, sigma, nu, tau),
                 -log(inverse_jacobian), tolerance = 1e-11)
    expect_equal(SVEMnet:::.svem_wmt_shasho_upper_tail(x, mu, sigma, nu, tau),
                 1 - p, tolerance = 1e-12)
  }
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
  fitted_density <- do.call(SVEMnet:::.svem_wmt_shasho_log_density,
                             c(list(x = x), as.list(fit$parameters)))
  expect_equal(fit$loglik, sum(fitted_density), tolerance = 1e-8)
  expect_gte(fit$loglik,
             sum(SVEMnet:::.svem_wmt_shasho_log_density(x, 2, 3, .6, .8)) - 1e-8)
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
  for (x in list(1:4, rep(1, 20), c(1:6, NA), c(1:6, Inf), letters[1:6])) {
    expect_error(SVEMnet:::.svem_wmt_fit_shasho(x), "finite|numeric|distinct")
  }
  expect_error(SVEMnet:::.svem_wmt_shasho_upper_tail(1, 0, 0, 0, 1),
               "positive|sigma")
})

test_that("clustered permutation distances admit a stationary numerical fit", {
  skip_on_cran()
  # Captured from SVEMnet's logical/ordered-predictor WMT regression case.
  # Closely spaced lower distances challenge finite-difference diagnostics.
  x <- c(0.70223517305708572, 0.7000286199695539, 0.70036136437865693,
         0.72515472220679389, 2.6068503979483215, 1.0776373783520117,
         1.6447773267594754, 4.0402039782416903, 3.7250244284547196,
         3.0485624581368729, 0.77016262842733751, 3.139462828740351,
         0.75396591159814719, 1.5944190609787379, 1.1044195576598059,
         0.72985735848223909, 1.5420688360859978, 1.1911278350418719,
         1.1328990718151832, 1.8581570871848749)
  fit <- SVEMnet:::.svem_wmt_fit_shasho(x)
  expect_true(fit$converged)
  expect_lte(fit$gradient_max, 1e-4)
  expect_true(all(is.finite(fit$parameters)))
  log_density <- do.call(SVEMnet:::.svem_wmt_shasho_log_density,
                         c(list(x = x), as.list(fit$parameters)))
  expect_equal(fit$loglik, sum(log_density), tolerance = 1e-8)
  tail <- do.call(SVEMnet:::.svem_wmt_shasho_upper_tail,
                  c(list(q = sort(x)), as.list(fit$parameters)))
  expect_true(all(is.finite(tail) & tail >= 0 & tail <= 1))
  expect_true(all(diff(tail) <= 0))
})
