.callback_stabilization_data <- function(decimal = TRUE) {
  withr::local_seed(125)
  d <- data.frame(x = rnorm(40), dose = rep(if (decimal) c(.1, .3) else c(0, 2), 20),
                  z = rnorm(40))
  d$y <- 2 + 3 * d$x + rnorm(40, sd = .3)
  d$yb <- rbinom(40, 1, plogis(.3 + .7 * d$x))
  d
}

.expect_callback_member_equality <- function(a, b, data) {
  for (field in c("coef_matrix", "best_alphas", "best_lambdas",
                  "best_relax_gammas", "member_diagnostics")) {
    expect_identical(a[[field]], b[[field]], info = field)
  }
  expect_identical(predict(a, data), predict(b, data))
}

test_that("callback exclusions match numeric exclusions with constant masks", {
  withr::local_seed(125)
  cases <- data.frame(decimal = c(TRUE, FALSE, TRUE, FALSE),
                      relaxed = c(TRUE, FALSE, FALSE, TRUE),
                      complexity = c("edf", "edf", "support", "support"),
                      alpha = c(.5, 1, 1, .5))
  for (i in seq_len(nrow(cases))) {
    d <- .callback_stabilization_data(cases$decimal[i])
    f <- if (cases$decimal[i]) y ~ x + I((dose - .2)^2) + z else
      y ~ x + I((dose - 1)^2) + z
    args <- list(formula = f, data = d, nBoot = 2, glmnet_alpha = cases$alpha[i],
                  relaxed = cases$relaxed[i], complexity = cases$complexity[i])
    callback <- function(x, y, weights, ...) 3L
    set.seed(125); a <- do.call(SVEMnet, c(args, list(exclude = callback)))
    set.seed(125); b <- do.call(SVEMnet, c(args, list(exclude = 3L)))
    .expect_callback_member_equality(a, b, d)
    expect_true(all(a$coef_matrix[, 3L] == 0))
    expect_identical(a$diagnostics$numerically_constant_columns,
                     colnames(a$coef_matrix)[3L])
  }
  for (relaxed in c(FALSE, TRUE)) {
    d <- .callback_stabilization_data()
    args <- list(formula = yb ~ x + I((dose - .2)^2) + z, data = d,
                  family = "binomial", nBoot = 2, glmnet_alpha = c(.5, 1),
                  relaxed = relaxed)
    set.seed(125); a <- do.call(SVEMnet, c(args, list(exclude = function(...) 3L)))
    set.seed(125); b <- do.call(SVEMnet, c(args, list(exclude = 3L)))
    .expect_callback_member_equality(a, b, d)
  }
})

test_that("CV callbacks retain original fold inputs and numeric path selections", {
  d <- .callback_stabilization_data()
  for (family in c("gaussian", "binomial")) for (relaxed in c(FALSE, TRUE)) {
    f <- if (family == "gaussian") y ~ x + I((dose - .2)^2) + z else
      yb ~ x + I((dose - .2)^2) + z
    X <- model.matrix(f, d)[, -1L, drop = FALSE]
    captured <- list()
    callback <- function(x, y, weights, ...) {
      captured[[length(captured) + 1L]] <<- list(x = x, y = y, weights = weights)
      c(FALSE, FALSE, TRUE)
    }
    args <- list(formula = f, data = d, family = family, glmnet_alpha = .5,
                  nfolds = 5, repeats = 1, seed = 125, relaxed = relaxed)
    a <- do.call(glmnet_with_cv, c(args, list(exclude = callback)))
    b <- do.call(glmnet_with_cv, c(args, list(exclude = 3L)))
    expect_identical(coef(a), coef(b))
    expect_identical(predict(a, d), predict(b, d))
    for (field in c("best_alpha", "best_lambda", "best_gamma")) {
      expect_identical(a[[field]], b[[field]])
    }
    expect_true(all(vapply(captured, function(z) {
      rows <- match(z$x[, 1L], d$x)
      expected_y <- if (family == "gaussian") d$y[rows] else
        cbind(1 - d$yb[rows], d$yb[rows])
      identical(unname(z$x), unname(X[rows, , drop = FALSE])) &&
        identical(colnames(z$x), colnames(X)) &&
        isTRUE(all.equal(unname(z$y), unname(expected_y), check.attributes = FALSE)) &&
        length(z$weights) == length(rows) && all(z$weights == 1)
    }, logical(1))))
    nrows <- vapply(captured, function(z) nrow(z$x), integer(1))
    expect_true(nrow(d) %in% nrows && any(nrows < nrow(d)))
  }
})

test_that("one-column callbacks never see or expose the glmnet sentinel", {
  withr::local_seed(125)
  d <- .callback_stabilization_data()
  captured <- list()
  callback <- function(x, y, weights, ...) {
    captured[[length(captured) + 1L]] <<- list(x = x, y = y, weights = weights)
    integer()
  }
  set.seed(125); a <- SVEMnet(y ~ x, d, nBoot = 2, glmnet_alpha = .5,
                             complexity = "edf", exclude = callback, penalty.factor = .7)
  set.seed(125); b <- SVEMnet(y ~ x, d, nBoot = 2, glmnet_alpha = .5,
                             complexity = "edf", penalty.factor = .7)
  .expect_callback_member_equality(a, b, d)
  expect_identical(colnames(a$coef_matrix), c("(Intercept)", "x"))
  expect_true(all(vapply(captured, function(z)
    identical(colnames(z$x), "x") && identical(unname(z$x[, 1L]), d$x) &&
      identical(as.numeric(z$y), d$y) && all(z$weights > 0) &&
      abs(sum(z$weights) - nrow(d)) < 1e-10, logical(1))))

  captured <- list()
  a <- glmnet_with_cv(y ~ x, d, glmnet_alpha = .5, nfolds = 5, repeats = 1,
                      seed = 125, exclude = callback)
  b <- glmnet_with_cv(y ~ x, d, glmnet_alpha = .5, nfolds = 5, repeats = 1, seed = 125)
  expect_identical(coef(a), coef(b))
  expect_identical(predict(a, d), predict(b, d))
  expect_identical(a$best_lambda, b$best_lambda)
  expect_true(all(vapply(captured, function(z) identical(colnames(z$x), "x"), logical(1))))

  d$x <- 1
  never <- function(...) stop("callback must not run")
  a <- SVEMnet(y ~ x, d, nBoot = 2, exclude = never)
  b <- glmnet_with_cv(y ~ x, d, nfolds = 5, repeats = 1, exclude = never)
  expect_true(all(a$coef_matrix[, "x"] == 0))
  expect_equal(unname(coef(b)["x"]), 0)
})

test_that("callback compatibility validates original indices and penalty inputs", {
  X <- cbind(x = 1:8, constant = rep(1, 8), z = 8:1)
  invalid <- list(4, 0, -.5, 1.5, Inf, NA_real_, "3", list(3),
                   c(FALSE, NA, TRUE), TRUE)
  for (value in invalid) {
    callback <- local({ answer <- value; function(...) answer })
    design <- SVEMnet:::.svem_prepare_glmnet_design(X, exclude = callback)
    expect_error(design$exclude(design$x, 1:8, rep(1, 8)), "callback must return indices")
  }
  for (value in list(NULL, integer(), numeric(), c(FALSE, FALSE, TRUE), 3L)) {
    callback <- local({ answer <- value; function(...) answer })
    design <- SVEMnet:::.svem_prepare_glmnet_design(X, exclude = callback)
    expected <- if (length(value)) 3L else integer()
    expect_identical(design$exclude(design$x, 1:8, rep(1, 8)), expected)
  }
  one <- SVEMnet:::.svem_prepare_glmnet_design(X[, 1L, drop = FALSE],
                                              exclude = function(...) 2L)
  expect_error(one$exclude(one$x, 1:8, rep(1, 8)), "indices in 1..1", fixed = TRUE)
  expect_error(SVEMnet:::.svem_prepare_glmnet_design(
    X[, 1L, drop = FALSE], args = list(penalty.factor = 1 + 1i),
    exclude = function(...) NULL), "one positive finite number")
  for (pf in list(1, c(1, NA, 1), c(1, NaN, 1), c("1", "1", "1"))) {
    expect_error(SVEMnet:::.svem_prepare_glmnet_design(
      X, args = list(penalty.factor = pf), exclude = function(...) 3L),
      "penalty.factor.*must be numeric")
  }
  ordinary <- cbind(x = 1:8, z = 8:1)
  callback <- function(...) integer()
  design <- SVEMnet:::.svem_prepare_glmnet_design(ordinary, exclude = callback)
  expect_identical(design$exclude, callback)
  expect_null(design$args$penalty.factor)
})

test_that("callback fixed exclusions preserve EDF and penalty normalization", {
  withr::local_seed(125)
  d <- .callback_stabilization_data()
  f <- y ~ x + I((dose - .2)^2) + z
  X <- model.matrix(f, d)[, -1L, drop = FALSE]
  callback <- function(...) 3L
  for (pf in list(c(.5, 17, 4), c(1, -2, 1))) {
    set.seed(125); a <- SVEMnet(f, d, nBoot = 2, glmnet_alpha = .5,
                               penalty.factor = pf, exclude = callback)
    set.seed(125); b <- SVEMnet(f, d, nBoot = 2, glmnet_alpha = .5,
                               penalty.factor = pf, exclude = 3L)
    expect_identical(a$coef_matrix, b$coef_matrix)
    expect_identical(a$best_lambdas, b$best_lambdas)
    expect_identical(predict(a, d), predict(b, d))
  }
  context <- SVEMnet:::.svem_edf_context(X, d$y, rep(1, nrow(d)),
    list(penalty.factor = c(1, Inf, 1), exclude = callback))
  expect_null(context$reason)
  expect_identical(context$penalty, rep(1, 3))
  context <- SVEMnet:::.svem_edf_context(X, d$y, rep(1, nrow(d)),
    list(penalty.factor = c(2, Inf, 1), exclude = callback))
  expect_match(context$reason, "nonuniform penalties requires a fixed exclude vector")
})
