# Helper utilities and toy data generators for SVEMnet tests
# (Kept lightweight; all test files call skip_on_cran() at the top.)
suppressPackageStartupMessages({
  library(testthat)
})

gen_toy_df <- function(n = 60L, with_factor = FALSE) {
  set.seed(1)
  X1 <- runif(n); X2 <- runif(n); X3 <- rnorm(n)
  A <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)
  if (with_factor) {
    F <- factor(sample(c("lo","hi"), n, TRUE))
    y <- 1 + 2*X1 - X2 + 0.5*X3 + 3*A + 1.5*B + 0.5*C + (F == "hi") + rnorm(n, 0, 0.3)
    data.frame(y, X1, X2, X3, A, B, C, F)
  } else {
    y <- 1 + 2*X1 - X2 + 0.5*X3 + 3*A + 1.5*B + 0.5*C + rnorm(n, 0, 0.3)
    data.frame(y, X1, X2, X3, A, B, C)
  }
}
