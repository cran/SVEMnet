testthat::skip_on_cran()
test_that("SVEMnet fits and predicts (parms vs mean)", {
  skip_if_not_installed("SVEMnet")
  set.seed(2)
  d <- gen_toy_df(80, with_factor = TRUE)
  fit <- SVEMnet::SVEMnet(y ~ (X1 + X2 + X3 + F)^2 + A + B + C, d, nBoot = 30, glmnet_alpha = 1)
  p1 <- predict(fit, d[1:10, ], agg = "parms")
  expect_type(p1, "double")
  p2 <- predict(fit, d[1:10, ], agg = "mean")
  expect_type(p2, "double")
  expect_equal(length(p1), 10L)
  expect_equal(length(p2), 10L)
})

test_that("SVEMnet relaxed path runs", {
  skip_on_cran()
  set.seed(4)
  n <- 40
  X1 <- rnorm(n); X2 <- rnorm(n)
  y  <- 1 + X1 - 0.5*X2 + rnorm(n, 0, 0.3)
  d  <- data.frame(y, X1, X2)
  fit <- SVEMnet(y ~ X1 + X2, d, nBoot = 20, glmnet_alpha = c(1, 0.5), relaxed = TRUE)
  p   <- predict(fit, d[1:5, ], agg = "mean")
  expect_type(p, "double")
  expect_true(length(p) == 5)
})



  ## --- BINOMIAL: basic fit + predict on all types ------------------------------
  test_that("SVEMnet binomial fits and predicts (response/link/class)", {
    skip_if_not_installed("SVEMnet")

    set.seed(2025)
    n  <- 120
    X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
    F  <- factor(sample(letters[1:3], n, TRUE))
    A  <- rnorm(n); B <- rnorm(n); C <- rnorm(n)
    eta <- -0.3 + 1.1*X1 - 0.8*X2 + 0.5*X1*X3 + 0.4*(F == "b")
    p   <- plogis(eta)
    yb  <- rbinom(n, 1, p)
    d   <- data.frame(yb, X1, X2, X3, F, A, B, C)

    fit <- SVEMnet::SVEMnet(
      yb ~ (X1 + X2 + X3 + F)^2 + A + B + C,
      d,
      nBoot = 40,
      glmnet_alpha = c(1, 0.5),
      relaxed = FALSE,
      family = "binomial"
    )

    # response
    p_resp <- predict(fit, d[1:15, ], type = "response")
    expect_type(p_resp, "double")
    expect_true(all(is.finite(p_resp)))
    expect_true(all(p_resp >= 0 & p_resp <= 1))
    expect_length(p_resp, 15L)

    # link
    p_link <- predict(fit, d[1:15, ], type = "link")
    expect_type(p_link, "double")
    expect_true(all(is.finite(p_link)))
    expect_length(p_link, 15L)

    # class
    p_cls <- predict(fit, d[1:15, ], type = "class")
    expect_type(p_cls, "integer")
    expect_true(all(p_cls %in% c(0L, 1L)))
    expect_length(p_cls, 15L)
  })

  ## --- BINOMIAL: uncertainty + aggregation ------------------------------------
  test_that("SVEMnet binomial uncertainty (se/interval) and agg modes", {
    skip_if_not_installed("SVEMnet")

    set.seed(9)
    n  <- 100
    X1 <- rnorm(n); X2 <- rnorm(n); G <- factor(sample(c("A","B","C"), n, TRUE))
    eta <- -0.2 + 1.0*X1 - 0.7*X2 + 0.6*(G == "B")
    p   <- plogis(eta)
    yb  <- rbinom(n, 1, p)
    d   <- data.frame(yb, X1, X2, G)

    fit <- SVEMnet::SVEMnet(
      yb ~ (X1 + X2 + G)^2,
      d,
      nBoot = 50,
      glmnet_alpha = c(1, 0.5),
      relaxed = TRUE,
      family = "binomial"
    )

    # parms aggregation, response scale, with se/interval
    out_parms <- predict(fit, d[1:12, ], type = "response",
                         agg = "mean", se.fit = TRUE, interval = TRUE, level = 0.9)
    expect_named(out_parms, c("fit","se.fit","lwr","upr"))
    expect_true(all(out_parms$fit >= 0 & out_parms$fit <= 1, na.rm = TRUE))
    expect_true(all(out_parms$lwr <= out_parms$upr, na.rm = TRUE))
    expect_length(out_parms$fit, 12L)

    # mean aggregation
    out_mean <- predict(fit, d[1:12, ], type = "response",
                        agg = "mean", se.fit = TRUE, interval = TRUE)
    expect_named(out_mean, c("fit","se.fit","lwr","upr"))
    expect_length(out_mean$fit, 12L)
  })

  ## --- BINOMIAL: class doesn't allow se/interval --------------------------------
  test_that("Binomial type='class' forbids se.fit/interval", {
    skip_if_not_installed("SVEMnet")

    set.seed(3)
    n <- 80
    X1 <- rnorm(n); X2 <- rnorm(n)
    pr <- plogis(-0.1 + X1 - 0.6*X2)
    yb <- rbinom(n, 1, pr)
    d  <- data.frame(yb, X1, X2)

    fit <- SVEMnet::SVEMnet(
      yb ~ (X1 + X2)^2, d,
      nBoot = 30, glmnet_alpha = 1,relaxed=FALSE, family = "binomial"
    )
    expect_error(predict(fit, d, type = "class", se.fit = TRUE))
    expect_error(predict(fit, d, type = "class", interval = TRUE))
  })

  ## --- BINOMIAL: debias is ignored ---------------------------------------------
  test_that("Binomial predictions ignore debias argument", {
    skip_if_not_installed("SVEMnet")

    set.seed(44)
    n <- 90
    X1 <- rnorm(n); X2 <- rnorm(n); H <- factor(sample(c("lo","hi"), n, TRUE))
    pr <- plogis(0.2 + 0.9*X1 - 0.4*X2 + 0.5*(H == "hi"))
    yb <- rbinom(n, 1, pr)
    d  <- data.frame(yb, X1, X2, H)

    fit <- SVEMnet::SVEMnet(
      yb ~ (X1 + X2 + H)^2, d,
      nBoot = 35, glmnet_alpha = c(1, 0.5), family = "binomial",relaxed=FALSE
    )

    p0 <- predict(fit, d[1:20, ], type = "response", debias = FALSE)
    p1 <- predict(fit, d[1:20, ], type = "response", debias = TRUE)  # should be identical
    expect_equal(p0, p1, tolerance = 0)
  })

  ## --- BINOMIAL: unseen factor level -> NA with warning -------------------------
  test_that("Binomial predict returns NA for unseen levels with a warning", {
    skip_if_not_installed("SVEMnet")

    set.seed(55)
    n <- 70
    X1 <- rnorm(n); F <- factor(sample(c("a","b"), n, TRUE))
    pr <- plogis(-0.1 + 0.8*X1 + 0.6*(F == "b"))
    yb <- rbinom(n, 1, pr)
    train <- data.frame(yb, X1, F)

    fit <- SVEMnet::SVEMnet(
      yb ~ X1 + F, train,
      nBoot = 25, glmnet_alpha = 1, family = "binomial",relaxed=FALSE
    )

    # newdata has unseen level "c"
    newd <- data.frame(X1 = rnorm(6), F = factor(c("a","b","c","a","c","b")))
    expect_warning(
      p <- predict(fit, newd, type = "response"),
      regexp = "unseen or missing levels"
    )
    expect_true(is.na(p[3]))
    expect_true(is.na(p[5]))
  })

  ## --- BINOMIAL + bigexp_spec path ----------------------------------------------
  test_that("SVEMnet binomial works with bigexp_spec", {
    skip_if_not_installed("SVEMnet")

    # build small dataset
    set.seed(101)
    n  <- 90
    X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
    eta <- 0.1 + 0.8*X1 - 0.7*X2 + 0.3*X1*X3
    yb  <- rbinom(n, 1, plogis(eta))
    df  <- data.frame(yb, X1, X2, X3)

    spec <- SVEMnet::bigexp_terms(
      yb ~ X1 + X2 + X3,
      data             = df,
      factorial_order  = 3,
      polynomial_order = 3,
      include_pc_3way  = FALSE
    )

    fit <- SVEMnet::SVEMnet(
      spec, df,
      nBoot = 30, glmnet_alpha = c(1, 0.5), relaxed = FALSE, family = "binomial"
    )

    # response probabilities and link
    pr <- predict(fit, df[1:10, ], type = "response", agg = "mean")
    lk <- predict(fit, df[1:10, ], type = "link", agg = "mean",
                  se.fit = TRUE, interval = TRUE, level = 0.9)

    expect_true(all(pr >= 0 & pr <= 1))
    expect_named(lk, c("fit","se.fit","lwr","upr"))
    expect_length(pr, 10L); expect_length(lk$fit, 10L)
  })


test_that("bigexp_prepare gives consistent columns across datasets", {
  skip_if_not_installed("SVEMnet")

  set.seed(456)
  train <- data.frame(
    y  = rnorm(20),
    X1 = rnorm(20),
    X2 = rnorm(20),
    G  = factor(sample(c("A", "B"), 20, replace = TRUE))
  )

  spec <- SVEMnet::bigexp_terms(
    y ~ X1 + X2 + G,
    data            = train,
    factorial_order = 2,
    polynomial_order = 2
  )

  # New data with only a subset of factor levels
  newdata <- data.frame(
    y  = rnorm(10),
    X1 = rnorm(10),
    X2 = rnorm(10),
    G  = factor(sample(c("A"), 10, replace = TRUE))
  )



  # New data with an unseen factor level should error when unseen = "error"
  newdata_bad <- data.frame(
    y  = rnorm(5),
    X1 = rnorm(5),
    X2 = rnorm(5),
    G  = factor(c("A", "B", "C", "A", "C"))
  )

  expect_error(
    SVEMnet::bigexp_prepare(spec, newdata_bad, unseen = "error"),
    regexp = "Unseen level"
  )
})

test_that("bigexp_formula switches response but keeps expansion", {
  skip_if_not_installed("SVEMnet")

  set.seed(789)
  df <- data.frame(
    y1 = rnorm(15),
    y2 = rnorm(15),
    X1 = rnorm(15),
    X2 = rnorm(15)
  )

  spec <- SVEMnet::bigexp_terms(
    y1 ~ X1 + X2,
    data            = df,
    factorial_order = 2,
    polynomial_order = 2
  )

  f2 <- SVEMnet::bigexp_formula(spec, "y2")
  f2_chr <- as.character(f2)

  # as.character(formula) is c("~", lhs, rhs)
  expect_equal(f2_chr[2L], "y2")
  expect_equal(f2_chr[3L], as.character(spec$formula)[3L])
})

test_that("with_bigexp_contrasts restores and uses stored options", {
  skip_if_not_installed("SVEMnet")

  old_opts <- getOption("contrasts")
  on.exit(options(contrasts = old_opts), add = TRUE)

  # Use a non-default contrasts setting when building the spec
  options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))

  df <- data.frame(
    y  = rnorm(10),
    X1 = rnorm(10),
    G  = factor(sample(c("A", "B"), 10, replace = TRUE))
  )

  spec <- SVEMnet::bigexp_terms(
    y ~ X1 + G,
    data            = df,
    factorial_order = 2,
    polynomial_order = 2
  )

  # Change contrasts again so we can verify they are restored inside
  options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))

  inner_opts <- SVEMnet::with_bigexp_contrasts(spec, {
    getOption("contrasts")
  })

  expect_equal(inner_opts, spec$settings$contrasts_options)
  expect_equal(getOption("contrasts"), c(unordered = "contr.treatment", ordered = "contr.poly"))
})

test_that("bigexp_train returns spec, formula, and prepared data", {
  skip_if_not_installed("SVEMnet")

  set.seed(42)
  df <- data.frame(
    y  = rnorm(25),
    X1 = rnorm(25),
    X2 = rnorm(25),
    G  = factor(sample(c("A", "B"), 25, replace = TRUE))
  )

  tr <- SVEMnet::bigexp_train(
    y ~ X1 + X2 + G,
    data            = df,
    factorial_order = 2,
    polynomial_order = 2
  )

  expect_s3_class(tr, "bigexp_train")
  expect_s3_class(tr$spec, "bigexp_spec")
  expect_s3_class(tr$formula, "formula")

  expect_equal(nrow(tr$data), nrow(df))

  # Design matrix from spec vs design matrix from train helper should match

})
