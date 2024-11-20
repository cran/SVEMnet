#' Plot Method for SVEM Models
#'
#' Plots actual versus predicted values for an \code{svem_model} using \code{ggplot2}.
#'
#' @param x An object of class \code{svem_model}.
#' @param plot_debiased Logical; if \code{TRUE}, includes debiased predictions if available (default is \code{TRUE}).
#' @param ... Additional arguments passed to \code{ggplot2} functions.
#' @return A \code{ggplot} object showing actual versus predicted values.
#' @details
#' This function creates an actual vs. predicted plot for the SVEM model. If \code{plot_debiased} is \code{TRUE} and debiased predictions are available, it includes them in the plot.
#'
#' **Plot Features:**
#' \itemize{
#'   \item **Actual vs. Predicted Points:** Plots the actual response values against the predicted values from the SVEM model.
#'   \item **Debiased Predictions:** If available and \code{plot_debiased} is \code{TRUE}, debiased predictions are included.
#'   \item **Ideal Fit Line:** A dashed line representing perfect prediction (slope = 1, intercept = 0) is included for reference.
#' }
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @import ggplot2
#' @importFrom stats model.response

#' @export
plot.svem_model <- function(x, plot_debiased = TRUE, ...) {
  actual_y <- x$actual_y
  y_pred <- x$y_pred

  # Create a data frame for plotting
  plot_data <- data.frame(
    Actual = actual_y,
    Predicted = y_pred,
    Type = "Predictions"
  )

  # Add debiased predictions if requested and available
  if (plot_debiased && !is.null(x$y_pred_debiased)) {
    plot_data_debiased <- data.frame(
      Actual = actual_y,
      Predicted = x$y_pred_debiased,
      Type = "Debiased Predictions"
    )
    plot_data <- rbind(plot_data, plot_data_debiased)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = .data$Actual, y = .data$Predicted, color = .data$Type, shape = .data$Type)) +
    geom_point(...) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 1) +
    labs(title = paste("Actual vs Predicted Values -", class(x)[2]),
         x = "Actual y",
         y = "Predicted y") +
    theme_minimal() +
    scale_color_manual(values = c("Predictions" = "blue", "Debiased Predictions" = "red")) +
    scale_shape_manual(values = c("Predictions" = 16, "Debiased Predictions" = 17)) +
    coord_fixed()

  print(p)
}

