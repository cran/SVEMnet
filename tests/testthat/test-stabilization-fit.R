test_that("formula offsets fail explicitly at every fitting front end", {
  d <- data.frame(y = 1:20, x = seq(-1, 1, length.out = 20), z = 1:20)
  for (fun in list(SVEMnet, glmnet_with_cv, forward_aicc, svem_forward,
                   bigexp_terms)) {
    expect_error(fun(y ~ x + offset(z), d), "offset\\(\\)")
  }
  fit <- local({
    offset <- 0.2
    bend <- function(x) (x - offset)^2
    SVEMnet(y ~ bend(x), d, nBoot = 2, relaxed = FALSE)
  })
  expect_true(all(is.finite(predict(fit, d))))
})

test_that("evaluated categorical frames are frozen for every prediction consumer", {
  set.seed(36031)
  d <- data.frame(x = rep(1:3, each = 12))
  d$y <- d$x + rnorm(nrow(d), sd = 0.1)
  fit <- SVEMnet(y ~ factor(x), d, nBoot = 3, glmnet_alpha = 1,
                 relaxed = FALSE, store_member_weights = TRUE)
  cv <- glmnet_with_cv(y ~ factor(x), d, nfolds = 4, repeats = 1, seed = 3)
  fw <- forward_aicc(y ~ factor(x), d)
  sf <- svem_forward(y ~ factor(x), d, nBoot = 3)
  one <- d[13, , drop = FALSE]
  for (object in list(fit, cv, fw, sf)) {
    expect_equal(predict(object, one), predict(object, d)[13])
    expect_warning(bad <- predict(object, data.frame(x = 4L)), "unseen or missing")
    expect_true(is.na(bad))
  }
  expect_equal(rowMeans(SVEMnet:::.svem_member_predictions(fit, one)$pred_mat),
               predict(fit, one))
  expect_equal(svem_ij_variance(fit, one)$fit, predict(fit, one))
  expect_equal(predict(fit, one, interval = TRUE)$fit, predict(fit, one))

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_equal(predict(readRDS(path), one), predict(fit, one))
  old <- fit
  old$schema$model_frame <- NULL
  saveRDS(old, path)
  expect_equal(predict(readRDS(path), d), predict(fit, d))

  dc <- transform(d, a = ifelse(x == 1, "a", "b"), b = as.character(x))
  fc <- SVEMnet(y ~ paste(a, b), dc, nBoot = 2, relaxed = FALSE)
  expect_equal(predict(fc, dc[13, , drop = FALSE]), predict(fc, dc)[13])
  fl <- SVEMnet(y ~ I(x > 1), d, nBoot = 2, relaxed = FALSE)
  expect_equal(predict(fl, one), predict(fl, d)[13])
  explicit <- SVEMnet(y ~ factor(x, levels = c(3, 2, 1)), d,
                      nBoot = 2, relaxed = FALSE)
  expect_equal(predict(explicit, one), predict(explicit, d)[13])
  df <- transform(d, x = factor(x))
  refit <- SVEMnet(y ~ relevel(x, ref = "2"), df, nBoot = 2,
                   relaxed = FALSE)
  expect_equal(predict(refit, df[1, , drop = FALSE]), predict(refit, df)[1])
  ordered <- SVEMnet(y ~ ordered(x), d, nBoot = 2, relaxed = FALSE)
  expect_equal(predict(ordered, one), predict(ordered, d)[13])
})

test_that("raw ranges and direct categorical supports survive transformations", {
  set.seed(36032)
  d <- data.frame(x = seq(0.2, 2, length.out = 30), y = rnorm(30))
  d$y[30] <- NA_real_
  fit <- local({
    center <- 0.3
    SVEMnet(y ~ I(x - center), d, nBoot = 2, relaxed = FALSE)
  })
  expect_identical(fit$sampling_schema$predictor_vars, "x")
  expect_equal(unname(fit$sampling_schema$num_ranges[, "x"]), range(d$x[-30]))
  tab <- svem_random_table_multi(list(fit), n = 8)$all
  expect_true(all(is.finite(tab$y_pred)))
  for (fun in list(forward_aicc, function(f, d) svem_forward(f, d, nBoot = 2))) {
    object <- fun(y ~ I(x^2), d)
    expect_equal(unname(object$sampling_schema$num_ranges[, "x"]), range(d$x[-30]))
  }

  dc <- data.frame(x = rep(1:3, each = 10), y = rnorm(30))
  fc <- SVEMnet(y ~ as.factor(x), dc, nBoot = 2, relaxed = FALSE)
  expect_equal(fc$sampling_schema$discrete_levels$x, 1:3)
  sampled <- svem_random_table_multi(list(fc), n = 15)$all
  expect_true(all(sampled$x %in% 1:3))
  expect_true(all(is.finite(sampled$y_pred)))
  tt <- terms(y ~ cut(x, breaks = c(-Inf, 1.5, Inf)), data = dc)
  expect_identical(SVEMnet:::.svem_direct_factor_support(tt, dc), list())
})

test_that("spec logical and ordered predictors accept raw and prepared data", {
  set.seed(36033)
  d <- data.frame(x = runif(36), flag = rep(c(FALSE, TRUE), 18),
                  grade = ordered(rep(c("low", "mid", "high"), 12),
                                  levels = c("low", "mid", "high")))
  d$y <- d$x + d$flag + as.integer(d$grade) + rnorm(36, sd = 0.1)
  spec <- bigexp_terms(y ~ x + flag + grade, d, factorial_order = 1,
                       polynomial_order = 1, report = FALSE)
  for (fun in list(function(s, d) SVEMnet(s, d, nBoot = 2, relaxed = FALSE),
                   function(s, d) svem_forward(s, d, nBoot = 2))) {
    fit <- fun(spec, d)
    expect_identical(unname(fit$sampling_schema$var_classes[c("flag", "grade")]),
                     c("logical", "ordered"))
    expect_equal(predict(fit, d), predict(fit, bigexp_prepare(spec, d)$data))
    tab <- svem_random_table_multi(list(fit), n = 8)$all
    expect_true(is.logical(tab$flag))
    expect_true(is.ordered(tab$grade))
    expect_true(all(is.finite(predict(fit, tab))))
  }
  plain <- SVEMnet(y ~ x + flag, d, nBoot = 2, relaxed = FALSE)
  expect_null(plain$xlevels$flag)
  expect_identical(plain$sampling_schema$factor_levels$flag, c("FALSE", "TRUE"))
  expect_true(is.logical(svem_random_table_multi(list(plain), n = 8)$data$flag))

  blocked_spec <- bigexp_terms(y ~ x, d, blocking = c("flag", "grade"),
                               report = FALSE)
  for (fun in list(function(s, d) SVEMnet(s, d, nBoot = 2, relaxed = FALSE),
                   function(s, d) svem_forward(s, d, nBoot = 2))) {
    fit <- fun(blocked_spec, d)
    expect_equal(predict(fit, d), predict(fit, bigexp_prepare(blocked_spec, d)$data))
    tab <- svem_random_table_multi(fit, n = 8)$data
    expect_true(is.logical(tab$flag))
    expect_true(is.ordered(tab$grade))
    expect_identical(unique(as.character(tab$grade)), "low")
    expect_true(all(is.finite(predict(fit, tab))))
  }
})

test_that("internal support and Gaussian criteria are invariant to units", {
  set.seed(5105)
  d <- as.data.frame(matrix(runif(24 * 5, -1, 1), 24, 5))
  names(d) <- paste0("x", 1:5)
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(24, sd = 0.5)
  f <- y ~ (x1 + x2 + x3 + x4 + x5)^2 + I(x1^2) + I(x2^2)
  dx <- d; dx$x2 <- dx$x2 * 1e8
  dy <- d; dy$y <- dy$y * 1e-12
  for (complexity in c("support", "edf")) {
    getfit <- function(data) {
      set.seed(10)
      SVEMnet(f, data, nBoot = 2, glmnet_alpha = 0.5,
               relaxed = FALSE, complexity = complexity)
    }
    a <- getfit(d); b <- getfit(dx); c <- getfit(dy)
    expect_identical(a$member_diagnostics$support_size, b$member_diagnostics$support_size)
    expect_identical(a$member_diagnostics$support_size, c$member_diagnostics$support_size)
    expect_equal(a$y_pred, b$y_pred, tolerance = 1e-8)
    expect_equal(a$y_pred, c$y_pred / 1e-12, tolerance = 1e-8)
    expect_equal(a$best_edfs, c$best_edfs, tolerance = 1e-8)
  }
  set.seed(10); a <- svem_forward(f, d, nBoot = 2)
  set.seed(10); b <- svem_forward(f, dy, nBoot = 2)
  expect_identical(a$member_diagnostics$support_size, b$member_diagnostics$support_size)
  expect_equal(a$y_pred, b$y_pred / 1e-12, tolerance = 1e-8)
  tiny <- SVEMnet(y ~ 1, data.frame(y = c(0, 0)), nBoot = 2,
                  objective = "wSSE")
  expect_identical(tiny$member_diagnostics$objective_value, c(0, 0))
})

test_that("intercept-only member predictions and alias diagnostics are usable", {
  d <- data.frame(y = seq_len(20), F = factor(rep(c("b", "c"), each = 10),
                                            levels = c("a", "b", "c")))
  fit <- SVEMnet(y ~ 1, d, nBoot = 3)
  expect_equal(rowMeans(SVEMnet:::.svem_member_predictions(fit, d)$pred_mat),
               predict(fit, d))
  expect_warning(fw <- forward_aicc(y ~ F, d), "Removed unobserved")
  expect_identical(fw$selected_terms, "F")
  expect_warning(sf <- svem_forward(y ~ F, d, nBoot = 2), "Removed unobserved")
  expect_true(all(sf$member_diagnostics$support_size == 2L))
  expect_identical(sf$sampling_schema$factor_levels$F, c("b", "c"))
  expect_identical(sf$xlevels$F, c("b", "c"))

  full <- data.frame(y = seq_len(30), F = factor(rep(c("a", "b", "c"), each = 10)))
  spec <- bigexp_terms(y ~ F, full, report = FALSE)
  subset <- full[full$F != "a", ]
  locked <- SVEMnet(spec, subset, nBoot = 2, relaxed = FALSE)
  expect_identical(locked$sampling_schema$factor_levels$F, c("a", "b", "c"))
  expect_warning(locked_forward <- svem_forward(spec, subset, nBoot = 2),
                 "partially aliased")
  expect_identical(locked_forward$sampling_schema$factor_levels$F, c("a", "b", "c"))
})
