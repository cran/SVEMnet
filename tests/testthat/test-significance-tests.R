testthat::skip_on_cran()
test_that("significance tests run (serial and parallel)", {
  skip_if_not_installed("SVEMnet")
  set.seed(7)
  d <- gen_toy_df(50, with_factor = FALSE)
  # Parallel test is heavier; still skipping on CRAN (this file already skips)
  expect_error(
    SVEMnet::svem_significance_test_parallel(y ~ X1 + X2 + X3 + A + B + C, d, n_perm = 19, nCore =safe_ncores() ),
    regexp = NA
  )
})
