skip_on_cran()
skip_if_not_installed("SVEMnet")

test_that("svem_random_table_multi validates schemas and predicts", {
  set.seed(7)
  n <- 40
  X1 <- runif(n); X2 <- runif(n)
  A <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)
  F <- factor(sample(c("lo","hi"), n, TRUE))
  y <- 1 + 2*X1 - X2 + 3*A + 1.5*B + 0.5*C + (F=="hi") + rnorm(n, 0, 0.3)
  z <- 0.5 + 0.8*X1 + 0.4*X2 + rnorm(n, 0, 0.2)
  d <- data.frame(y, z, X1, X2, A, B, C, F)

  f1 <- SVEMnet::SVEMnet(y ~ X1 + X2 + A + B + C + F, d, nBoot = 30, glmnet_alpha = 1)
  f2 <- SVEMnet::SVEMnet(z ~ X1 + X2 + A + B + C + F, d, nBoot = 30, glmnet_alpha = 1)

  res <- SVEMnet::svem_random_table_multi(list(f1, f2), n = 50)

  expect_true(is.list(res))
  expect_true(all(c("data","pred","all") %in% names(res)))
  expect_true(all(c("X1","X2","A","B","C","F") %in% names(res$data)))
  expect_true(all(c("y","z") %in% names(res$pred) | all(c("y.pred","z.pred") %in% names(res$pred))))
  expect_equal(nrow(res$data), 50)
  expect_equal(nrow(res$pred), 50)
  expect_equal(nrow(res$all), 50)
})
