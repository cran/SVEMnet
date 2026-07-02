testthat::skip_on_cran()
test_that("significance tests run (serial and parallel)", {
  skip_if_not_installed("SVEMnet")
  set.seed(7)
  d <- gen_toy_df(50, with_factor = FALSE)
  # Parallel test is heavier; still skipping on CRAN (this file already skips)
  expect_error(
    SVEMnet::svem_significance_test_parallel(
      y ~ X1 + X2 + X3 + A + B + C, d,
      nPoint  = 400,
      nSVEM   = 3,
      nPerm   = 30,
      nBoot   = 30,
      nCore   = safe_ncores(),
      seed    = 99,
      verbose = FALSE
    ),
    regexp = NA
  )
})

test_that("significance test works with a single original fit (nSVEM = 1)", {
  skip_if_not_installed("SVEMnet")
  set.seed(8)
  d <- gen_toy_df(50, with_factor = FALSE)
  res <- svem_significance_test_parallel(
    y ~ X1 + X2 + X3 + A + B + C, d,
    nPoint  = 300,
    nSVEM   = 1,
    nPerm   = 30,
    nBoot   = 25,
    nCore   = safe_ncores(),
    seed    = 42,
    verbose = FALSE
  )
  expect_s3_class(res, "svem_significance_test")
  expect_length(res$p_values, 1L)
  expect_true(is.finite(res$p_value))
})
