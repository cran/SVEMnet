skip_on_cran()

test_that("svem_random_table_multi returns predictor columns plus response columns", {
  set.seed(1)

  n  <- 60
  X1 <- runif(n); X2 <- runif(n)
  A  <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)  # simple mixture-ish
  F  <- factor(sample(c("lo","hi"), n, TRUE))

  # two responses sharing the same factor space
  y1 <- 1 + 2*X1 - X2 + 3*A + 1.5*B + 0.5*C + (F == "hi") + rnorm(n, 0, 0.3)
  y2 <- 0.5 + 0.8*X1 + 0.4*X2 + rnorm(n, 0, 0.2)

  d <- data.frame(y1, y2, X1, X2, A, B, C, F)

  f1 <- SVEMnet(y1 ~ X1 + X2 + A + B + C + F, d, nBoot = 30, glmnet_alpha = 1)
  f2 <- SVEMnet(y2 ~ X1 + X2 + A + B + C + F, d, nBoot = 30, glmnet_alpha = 1)

  out <- svem_random_table_multi(list(f1, f2), n = 50)

  # Support both return types:
  # 1) list(data = predictors, pred = data.frame of responses)
  # 2) single data.frame containing predictors + response columns
  if (is.list(out) && all(c("data", "pred") %in% names(out))) {
    expect_s3_class(out$data, "data.frame")
    expect_s3_class(out$pred, "data.frame")
    expect_equal(nrow(out$data), 50)
    expect_equal(nrow(out$pred), 50)
    combined <- cbind(out$data, out$pred)
  } else if (is.data.frame(out)) {
    expect_equal(nrow(out), 50)
    combined <- out
  } else {
    fail("svem_random_table_multi returned an unexpected type.")
  }

  # predictors must be present
  expect_true(all(c("X1","X2","A","B","C","F") %in% names(combined)))

  # response columns must be appended with the response names
  expect_true(all(c("y1","y2") %in% names(combined)))

  # predictions should not be NA
  expect_false(anyNA(combined$y1))
  expect_false(anyNA(combined$y2))
})
