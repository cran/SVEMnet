center_toy_data <- function() {
  set.seed(3600)
  n <- 48L
  d <- data.frame(x = runif(n, 10, 14), z = runif(n, -5, -1),
                  m = runif(n, .1, .9),
                  g = factor(rep(c("a", "b", "c"), length.out = n)))
  d$y <- 2 + .8 * d$x + 2 * (d$x - 12)^2 - .5 * d$z +
    .3 * d$m + (d$g == "b") + rnorm(n, sd = .15)
  d$y2 <- d$y + .5 * d$z
  d
}

test_that("polynomial centering is opt-in and preserves the default spec", {
  d <- center_toy_data()
  s <- bigexp_terms(y ~ x + z + m + g, d, report = FALSE)
  off <- bigexp_terms(y ~ x + z + m + g, d,
                      center_polynomials = FALSE, report = FALSE)
  expect_identical(off, s)
  expect_null(s$settings$polynomial_centers)
  expect_match(s$rhs, "I(x^2)", fixed = TRUE)
  expect_match(s$rhs, "I(x^3)", fixed = TRUE)
  expect_identical(names(formals(bigexp_terms))[1:14],
    c("formula", "data", "factorial_order", "polynomial_order",
      "include_pc_2way", "include_pc_3way", "intercept", "blocking",
      "discrete_numeric", "audit", "audit_numeric_rate", "audit_unique_ratio",
      "audit_min_n", "report"))
})

test_that("fixed polynomial means exclude mixtures and leave raw terms and groups intact", {
  d <- center_toy_data()
  before <- d
  old_contrasts <- getOption("contrasts")
  on.exit(options(contrasts = old_contrasts), add = TRUE)
  for (contrast in c("contr.treatment", "contr.sum")) {
    options(contrasts = c(contrast, "contr.poly"))
    raw <- bigexp_terms(y ~ x + z + m + g, d, 3, 3,
                        include_pc_3way = TRUE, report = FALSE)
    centered <- bigexp_terms(y ~ x + z + m + g, d, 3, 3,
      include_pc_3way = TRUE, center_polynomials = TRUE,
      mixture_vars = "m", report = FALSE)
    X0 <- model.matrix(raw$formula, d)
    X <- model.matrix(centered$formula, d)
    expect_equal(centered$settings$polynomial_centers, c(x = mean(d$x), z = mean(d$z)))
    expect_identical(centered$settings$mixture_vars, "m")
    expect_identical(dim(X), dim(X0))
    expect_identical(attr(X, "assign"), attr(X0, "assign"))
    expect_identical(attr(X, "contrasts"), attr(X0, "contrasts"))
    raw_cols <- !grepl("I(", colnames(X0), fixed = TRUE)
    expect_identical(X[, raw_cols], X0[, raw_cols])
    expect_identical(unname(X[, "I(m^2)"]), d$m^2)
    expect_identical(unname(X[, "I(m^3)"]), d$m^3)
    expect_identical(unname(X[, "x:I(m^2)"]), d$x * d$m^2)
    expect_identical(centered$num_range, raw$num_range)
    expect_identical(centered$vars, raw$vars)
    expect_identical(centered$is_cat, raw$is_cat)
  }
  expect_identical(d, before)
})

test_that("powered factors use full-precision frozen centers for new rows and responses", {
  d <- center_toy_data()
  d$x <- d$x + pi
  d$y[c(1, 4)] <- NA_real_
  d$z[3] <- NA_real_
  old_outdec <- getOption("OutDec")
  on.exit(options(OutDec = old_outdec), add = TRUE)
  options(OutDec = ",")
  s <- bigexp_terms(y ~ x + z + m, d, 2, 3,
                    center_polynomials = TRUE, mixture_vars = "m", report = FALSE)
  centers <- s$settings$polynomial_centers
  expect_identical(centers[["x"]], mean(d$x))
  expect_identical(centers[["z"]], mean(d$z[is.finite(d$z)]))
  expect_identical(environment(s$formula), baseenv())
  restored <- unserialize(serialize(s, NULL))
  nd <- d[c(5, 8, 11), ]
  nd$x <- c(1, 30, 45)
  prep <- bigexp_prepare(restored, nd)
  expect_identical(prep$data, nd)
  X <- model.matrix(bigexp_formula(restored, "y2"), prep$data)
  powered_x <- which(attr(X, "assign") == 4L)
  # Main effects x/z/m are followed by the first quadratic x term.
  expect_length(powered_x, 1L)
  expect_identical(unname(X[, powered_x]), (nd$x - centers[["x"]])^2)
  expected <- stats::as.formula(paste("y2 ~", s$rhs), env = baseenv())
  expect_identical(X, model.matrix(expected, nd))
  one_at_a_time <- do.call(rbind, lapply(seq_len(nrow(nd)), function(i) {
    model.matrix(restored$formula, nd[i, , drop = FALSE])
  }))
  expect_equal(unname(X), unname(one_at_a_time), ignore_attr = TRUE)
  s_other <- bigexp_terms(y2 ~ x + z + m, d, 2, 3,
                          center_polynomials = TRUE, mixture_vars = "m", report = FALSE)
  expect_identical(s_other$rhs, s$rhs)
  expect_identical(s_other$settings$polynomial_centers, centers)
})

test_that("centering respects blocking, discrete metadata, degree one, and exclusions", {
  d <- center_toy_data()
  d$block <- rep(1:3, length.out = nrow(d))
  s <- bigexp_terms(y ~ x + z + m, d, 2, 2, blocking = "block",
    discrete_numeric = "m", center_polynomials = TRUE, mixture_vars = "m", report = FALSE)
  expect_identical(names(s$settings$polynomial_centers), c("x", "z"))
  expect_identical(s$settings$discrete_levels$m, sort(unique(d$m)))
  only_mix <- bigexp_terms(y ~ m, d, 1, 2, center_polynomials = TRUE,
                           mixture_vars = "m", report = FALSE)
  raw_mix <- bigexp_terms(y ~ m, d, 1, 2, report = FALSE)
  expect_identical(only_mix$rhs, raw_mix$rhs)
  linear <- bigexp_terms(y ~ x + z, d, 2, 1,
                         center_polynomials = TRUE, report = FALSE)
  raw_linear <- bigexp_terms(y ~ x + z, d, 2, 1, report = FALSE)
  expect_identical(linear$rhs, raw_linear$rhs)
  expect_length(linear$settings$polynomial_centers, 0L)
  for (bad in list(NA, 1, c(TRUE, FALSE))) {
    expect_error(bigexp_terms(y ~ x, d, center_polynomials = bad, report = FALSE),
                 "center_polynomials")
  }
  for (bad in list(NA_character_, "", 1, "missing", "g", "block")) {
    expect_error(bigexp_terms(y ~ x + g, d, blocking = "block",
      mixture_vars = bad, report = FALSE), "mixture_vars")
  }
})

test_that("centered specs flow through fitters and prediction without fitter changes", {
  skip_on_cran()
  d <- center_toy_data()
  s <- bigexp_terms(y ~ x + z + g, d, 1, 2,
                    center_polynomials = TRUE, report = FALSE)
  plain <- as.formula(paste("y ~", s$rhs), env = baseenv())
  for (fit_fun in list(
    function(f) { set.seed(910); SVEMnet(f, d, nBoot = 6, glmnet_alpha = 1, relaxed = FALSE) },
    function(f) { set.seed(911); svem_forward(f, d, nBoot = 6) },
    function(f) forward_aicc(f, d),
    function(f) glmnet_with_cv(f, d, nfolds = 3, repeats = 1, glmnet_alpha = 1, seed = 912)
  )) {
    fitted <- fit_fun(s$formula)
    oracle <- fit_fun(plain)
    nd <- d[c(2, 5, 17), , drop = FALSE]
    expect_equal(predict(fitted, nd), predict(oracle, nd), tolerance = 1e-12)
    restored <- unserialize(serialize(fitted, NULL))
    expect_identical(predict(fitted, nd), predict(restored, nd))
  }
  set.seed(911)
  fit <- svem_forward(s, d, nBoot = 6)
  set.seed(914)
  table <- svem_random_table_multi(list(fit), n = 10)
  expect_true(all(table$data$x >= min(d$x) & table$data$x <= max(d$x)))
  expect_equal(table$pred$y_pred, as.numeric(predict(fit, table$data, debias = FALSE)))
})
