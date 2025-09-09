#' Coefficient Nonzero Percentages from an SVEM Model
#'
#' Calculates the percentage of bootstrap iterations in which each coefficient
#' (excluding the intercept) is nonzero, using a small tolerance.
#'
#' @param object An object of class \code{svem_model}.
#' @param tol Numeric tolerance for "nonzero" (default \code{1e-7}).
#' @param ... Unused.
#' @return Invisibly returns a data frame with variables and percentages.
#' @import ggplot2
#' @export
#' @method coef svem_model
coef.svem_model <- function(object, tol = 1e-7, ...) {
  if (!"coef_matrix" %in% names(object)) {
    stop("The provided model does not contain 'coef_matrix'. Is it a valid 'svem_model'?")
  }

  cm <- object$coef_matrix
  if ("(Intercept)" %in% colnames(cm)) {
    cm <- cm[, colnames(cm) != "(Intercept)", drop = FALSE]
  }

  if (ncol(cm) == 0) {
    message("No non-intercept coefficients found.")
    return(invisible(data.frame()))
  }

  pct <- as.data.frame(100 * colMeans(abs(cm) > tol))
  colnames(pct)[1] <- "Percent of Bootstraps Nonzero"
  pct$Variable <- rownames(pct)
  pct <- pct[order(-pct[[1]]), , drop = FALSE]
  rownames(pct) <- pct$Variable
  print(pct[, "Percent of Bootstraps Nonzero", drop = FALSE])

  pct$Order <- seq_len(nrow(pct))
  gp <- ggplot2::ggplot(pct, ggplot2::aes(x = Order, y = `Percent of Bootstraps Nonzero`)) +
    ggplot2::geom_line(group = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = Variable), vjust = -0.5, size = 3) +
    ggplot2::scale_x_continuous(breaks = pct$Order, labels = pct$Variable) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Coefficient Percent Nonzero (by bootstrap)",
      x = "Variables",
      y = "Percent of Bootstraps Nonzero"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title  = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
    )

  print(gp)
  invisible(pct)
}
