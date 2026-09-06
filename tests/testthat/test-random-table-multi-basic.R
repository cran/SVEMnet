skip_on_cran()
skip_if_not_installed("SVEMnet")

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

  expect_type(out, "list")
  expect_true(all(c("data", "pred", "all") %in% names(out)))
  expect_s3_class(out$data, "data.frame")
  expect_s3_class(out$pred, "data.frame")
  expect_equal(nrow(out$data), 50L)
  expect_named(out$pred, c("y1_pred", "y2_pred"))
  expect_equal(out$all, cbind(out$data, out$pred))
  expect_true(all(c("X1", "X2", "A", "B", "C", "F") %in% names(out$data)))
  expect_true(all(is.finite(as.matrix(out$pred))))
  expect_equal(out$pred$y1_pred, as.numeric(predict(f1, out$data, debias = FALSE)))
  expect_equal(out$pred$y2_pred, as.numeric(predict(f2, out$data, debias = FALSE)))
})
