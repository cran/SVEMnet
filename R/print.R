#' Print Method for SVEM Significance Test
#'
#' Prints the median p-value from an object of class \code{svem_significance_test}.
#'
#' @param x An object of class \code{svem_significance_test}.
#' @param ... Additional arguments (unused).
#' @export
#' @method print svem_significance_test
print.svem_significance_test <- function(x, ...) {
  cat("SVEM Significance Test p-value:\n")
  print(x$p_value)
  invisible(x)
}
