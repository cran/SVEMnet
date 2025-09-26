# tests/testthat/helper-cores.R
safe_ncores <- function() {
  nc <- parallel::detectCores(logical = TRUE)
  nc <- if (is.finite(nc)) nc else 1L
  max(1L, min(2L, nc))  # CRAN-safe: <= 2
}
