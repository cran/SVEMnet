test_that("nested direct factor conversions are frozen before evaluation", {
  set.seed(6601)
  d <- data.frame(x = rep(1:3, each = 12), y = rep(c(1, 3, 2), each = 12) + rnorm(36, sd = .1))
  forms <- list(y ~ relevel(factor(x), ref = "2"),
                y ~ relevel(as.factor(x), ref = "2"),
                y ~ relevel(factor(x, labels = c("a", "b", "c")), ref = "b"),
                y ~ factor(x, c(3, 2, 1), labels = c("a", "b", "c")))
  for (f in forms) {
    fit <- SVEMnet(f, d, nBoot = 2, relaxed = FALSE, store_member_weights = TRUE)
    one <- d[1, , drop = FALSE]
    expect_equal(predict(fit, one), predict(fit, d)[1])
    expect_equal(rowMeans(SVEMnet:::.svem_member_predictions(fit, one)$pred_mat), predict(fit, one))
    expect_equal(svem_ij_variance(fit, one)$fit, predict(fit, one))
    original <- model.matrix(f, d)[, -1L, drop = FALSE]
    expect_identical(as.numeric(fit$training_X), as.numeric(original))
    expect_warning(unseen <- predict(fit, data.frame(x = 4)), "unseen or missing")
    expect_true(is.na(unseen))
  }
  cv <- glmnet_with_cv(forms[[1]], d, nfolds = 3, repeats = 1)
  expect_equal(predict(cv, d[1, , drop = FALSE]), predict(cv, d)[1])
  fit <- SVEMnet(forms[[1]], d, nBoot = 2, relaxed = FALSE)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_equal(predict(readRDS(path), d[1, , drop = FALSE]), predict(fit, d)[1])

  omitted <- d; omitted$y[omitted$x == 2] <- NA_real_
  expect_warning(fit <- SVEMnet(forms[[1]], omitted, nBoot = 2, relaxed = FALSE), "Removed unobserved")
  expect_true(is.finite(predict(fit, d[1, , drop = FALSE])))
  expect_warning(unknown <- predict(fit, data.frame(x = 2)), "unseen or missing")
  expect_true(is.na(unknown))
  omitted$g <- factor(omitted$x)
  expect_warning(raw_relevel <- SVEMnet(y ~ relevel(g, ref = "2"), omitted,
                                       nBoot = 2, relaxed = FALSE), "Removed unobserved")
  expect_true(is.finite(predict(raw_relevel, data.frame(g = "1"))))
  expect_warning(unknown <- predict(raw_relevel, data.frame(g = "2")), "unseen or missing")
  expect_true(is.na(unknown))

  local_formula <- local({factor <- function(x) x^2; y ~ factor(x)})
  tt <- attr(model.frame(local_formula, d), "terms")
  expect_identical(SVEMnet:::.svem_freeze_factor_predvars(tt, d), tt)
  expect_identical(SVEMnet:::.svem_direct_factor_support(tt, d), list())
  tt_null <- terms(y ~ factor(x, levels = NULL), data = d)
  attr(tt_null, "predvars") <- attr(tt_null, "variables")
  expect_identical(SVEMnet:::.svem_freeze_factor_predvars(tt_null, d), tt_null)

  excluded <- SVEMnet(y ~ factor(x, exclude = 3), d, nBoot = 2, relaxed = FALSE)
  expect_identical(excluded$sampling_schema$discrete_levels$x, c(1L, 2L))
  sampled <- svem_random_table_multi(excluded, n = 8)
  expect_true(all(sampled$data$x %in% 1:2))
  expect_true(all(is.finite(sampled$pred$y_pred)))

  dna <- d; dna$x[1] <- NA_real_
  nafit <- SVEMnet(y ~ factor(x, exclude = NULL), dna, nBoot = 2, relaxed = FALSE)
  expect_true(is.finite(predict(nafit, dna[1, , drop = FALSE])))
})

test_that("plain formulas drop unused levels after missing-row removal", {
  d <- data.frame(F = factor(rep(c("a", "b", "c"), each = 10)), y = rep(c(0, 1, 5), each = 10))
  d$y[d$F == "a"] <- NA_real_
  for (fun in list(function(f, d) SVEMnet(f, d, nBoot = 2, relaxed = FALSE),
                   function(f, d) glmnet_with_cv(f, d, nfolds = 3, repeats = 1),
                   forward_aicc, function(f, d) svem_forward(f, d, nBoot = 2))) {
    expect_warning(fit <- fun(y ~ F, d), "Removed unobserved")
    expect_identical(fit$xlevels$F, c("b", "c"))
    expect_warning(value <- predict(fit, data.frame(F = "a")), "unseen or missing")
    expect_true(is.na(value))
  }
  expect_error(SVEMnet(y ~ C(F, contr.sum), d, nBoot = 2), "explicit contrasts")
})

test_that("calibration remains stable at large prediction origins", {
  set.seed(6602)
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = rnorm(30))
  d$y <- 1e10 + 3 * d$x + rnorm(30, sd = .2)
  tolerance <- 8 * .Machine$double.eps * 1e10
  for (fun in list(function() SVEMnet(y ~ x + z, d, nBoot = 12, relaxed = FALSE),
                   function() svem_forward(y ~ x + z, d, nBoot = 12),
                   function() glmnet_with_cv(y ~ x + z, d, nfolds = 3, repeats = 1))) {
    fit <- fun()
    expect_identical(fit$debias_fit$rank, 2L)
    expect_true(sd(fit$y_pred_debiased) > 1)
    expect_lte(max(abs(predict(fit, d, debias = TRUE) - as.numeric(fit$y_pred_debiased))),
               tolerance)
    expect_identical(environment(fit$debias_fit$terms), baseenv())
    path <- tempfile(fileext = ".rds")
    saveRDS(fit, path)
    expect_equal(predict(readRDS(path), d, debias = TRUE), predict(fit, d, debias = TRUE))
    unlink(path)
  }
  X <- cbind(x = seq(-1, 1, length.out = 30))
  pred <- 2 + 3 * X[, 1]; y <- pred + sin(seq_along(pred))
  ordinary <- lm(y ~ pred)
  fit <- SVEMnet:::.svem_fit_calibration(y, pred, X, c(2, 3))$fit
  expect_identical(unname(coef(fit)), unname(coef(ordinary)))
  near <- SVEMnet:::.svem_fit_calibration(y, 1e10 + X[, 1] * 1e-6, X, c(1e10, 1e-6))
  expect_null(near$fit)
  expect_match(near$reason, "rounding-error")
})

test_that("rounding-constant columns are excluded with names and schema retained", {
  set.seed(6603)
  d <- data.frame(A = rep(c(.1, .3), each = 9), B = rep(c(.5, 1, 1.5), 6))
  d$y <- 10 + 20 * (d$A - .2) + 3 * (d$B - 1) + rnorm(18, sd = .2)
  f <- y ~ A + B + I((A - .2)^2)
  for (fun in list(function() SVEMnet(f, d, nBoot = 3, glmnet_alpha = c(.5, 1), relaxed = FALSE),
                   function() glmnet_with_cv(f, d, glmnet_alpha = c(.5, 1), nfolds = 3, repeats = 1))) {
    fit <- fun()
    expect_identical(fit$diagnostics$numerically_constant_columns, "I((A - 0.2)^2)")
    expect_identical(unname(fit$parms[["I((A - 0.2)^2)"]]), 0)
    expect_true(max(abs(fit$parms)) < 1e3)
    pred <- predict(fit, data.frame(A = .2, B = 1))
    expect_true(pred > min(d$y) && pred < max(d$y))
  }
  X <- cbind(decimal = (d$A - .2)^2, ordinary = d$B)
  for (scale in c(1e-6, 1, 1e6)) {
    expect_identical(SVEMnet:::.svem_numerically_constant_columns(X * scale), 1L)
  }
  expect_length(SVEMnet:::.svem_numerically_constant_columns(cbind(timestamp = 1e12 + 0:19)), 0)
  expect_length(SVEMnet:::.svem_numerically_constant_columns(cbind(square = d$A^2)), 0)
  expect_identical(SVEMnet:::.svem_numerically_constant_columns(cbind(constant = rep(.3, 20))), 1L)

  d$C <- rnorm(nrow(d))
  for (exclude in list(2L, function(x, y, weights, ...) 2L,
                       function(x, y, weights, ...) c(FALSE, TRUE, FALSE, FALSE))) {
    fit <- SVEMnet(update(f, . ~ . + C), d, nBoot = 2, relaxed = FALSE, exclude = exclude)
    expect_identical(unname(fit$parms[c("B", "I((A - 0.2)^2)")]), c(0, 0))
  }
  dc <- data.frame(y = 1:20, x = rep(c(.3, .30000000000000004), 10))
  for (fit in list(SVEMnet(y ~ x, dc, nBoot = 2),
                   glmnet_with_cv(y ~ x, dc), svem_forward(y ~ x, dc, nBoot = 2))) {
    expect_identical(unname(fit$parms[["x"]]), 0)
    expect_identical(fit$diagnostics$numerically_constant_columns, "x")
  }
})
