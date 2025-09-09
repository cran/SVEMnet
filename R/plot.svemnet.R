#' Plot Method for SVEM Models
#'
#' Plots actual versus predicted values for an \code{svem_model}.
#'
#' @param x An object of class \code{svem_model}.
#' @param plot_debiased Logical; if \code{TRUE}, includes debiased predictions if available.
#' @param ... Additional aesthetics passed to \code{geom_point()}.
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom utils tail
#' @export
#' @method plot svem_model
plot.svem_model <- function(x, plot_debiased = FALSE, ...) {
  actual_y <- x$actual_y
  y_pred   <- x$y_pred

  plot_data <- data.frame(
    Actual    = actual_y,
    Predicted = y_pred,
    Type      = "Predictions"
  )

  if (plot_debiased && !is.null(x$y_pred_debiased)) {
    plot_data <- rbind(
      plot_data,
      data.frame(Actual = actual_y, Predicted = x$y_pred_debiased, Type = "Debiased Predictions")
    )
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Actual, y = Predicted, color = Type, shape = Type)) +
    ggplot2::geom_point(...) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = paste("Actual vs Predicted -", utils::tail(class(x), 1)),
      x = "Actual y",
      y = "Predicted y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("Predictions" = "blue", "Debiased Predictions" = "red")) +
    ggplot2::scale_shape_manual(values = c("Predictions" = 16, "Debiased Predictions" = 17)) +
    ggplot2::coord_fixed()

  p
}
