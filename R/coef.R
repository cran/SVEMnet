#' Plot Coefficient Nonzero Percentages from a SVEMnet Model
#'
#' This function calculates the percentage of bootstrap iterations in which each coefficient
#' is nonzero.
#'
#' @param object An object of class \code{svem_model} returned by the \code{SVEMnet} function.
#' @param ... other arguments to pass.
#'
#' @return Invisibly returns a data frame containing the percentage of bootstraps
#'         where each coefficient is nonzero.
#'
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_minimal labs theme element_text
#' @importFrom stats median
#' @export
coef.svem_model <- function(object, ...) {
  # Check for coef_matrix
  if (!"coef_matrix" %in% names(object)) {
    stop("The provided model does not contain a 'coef_matrix'. Ensure it is a valid 'svem_model' object.")
  }

  coef_matrix <- object$coef_matrix

  # Exclude the intercept if present
  if ("(Intercept)" %in% colnames(coef_matrix)) {
    coef_matrix <- coef_matrix[, colnames(coef_matrix) != "(Intercept)", drop = FALSE]
  }

  if (ncol(coef_matrix) == 0) {
    message("No non-intercept coefficients found.")
    return(invisible(data.frame()))
  }

  # Percent (0-100), not fraction
  coef_pct <- as.data.frame(100 * colMeans(coef_matrix != 0))
  colnames(coef_pct)[1] <- "Percent of Bootstraps Nonzero"

  # Order by descending percentage
  coef_pct <- coef_pct[order(-coef_pct[, 1]), , drop = FALSE]

  # Print the data frame
  print(coef_pct)

  # Prep for plotting
  coef_pct$Variable <- rownames(coef_pct)
  coef_pct$Order <- seq_len(nrow(coef_pct))

  # Plot (namespaced calls to avoid import issues)
  plot <- ggplot2::ggplot(coef_pct, ggplot2::aes(x = Order, y = `Percent of Bootstraps Nonzero`)) +
    ggplot2::geom_line(group = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = Variable), vjust = -0.5, size = 3) +
    ggplot2::scale_x_continuous(breaks = coef_pct$Order, labels = coef_pct$Variable) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Coefficient Percentages with Labels",
      x = "Variables",
      y = "Percent of Bootstraps Nonzero"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
    )

  print(plot)
  invisible(coef_pct)
}
