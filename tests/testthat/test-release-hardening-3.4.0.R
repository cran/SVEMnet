skip_on_cran()

test_that("citation version cannot drift from package metadata", {
  ver <- as.character(utils::packageVersion("SVEMnet"))
  cit <- utils::citation("SVEMnet")
  expect_true(any(grepl(paste("R package version", ver), format(cit),
                        fixed = TRUE)))
})

test_that("bigexp count arguments reject fractional values", {
  dat <- data.frame(y = 1:6, x = seq_len(6))

  expect_error(
    bigexp_terms(y ~ x, dat, factorial_order = 1.5, report = FALSE),
    "integer"
  )
  expect_error(
    bigexp_terms(y ~ x, dat, polynomial_order = 2.5, report = FALSE),
    "integer"
  )
  expect_error(
    bigexp_terms(y ~ x, dat, audit_min_n = 2.5, report = FALSE),
    "integer"
  )
})
