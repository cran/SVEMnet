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
coef.svem_model <- function(object,...) {
  # Check if the model has a coef_matrix
  if (!"coef_matrix" %in% names(object)) {
    stop("The provided model does not contain a 'coef_matrix'. Ensure it is a valid 'svem_model' object.")
  }

  # Extract the coefficient matrix
  coef_matrix <- object$coef_matrix

  #  exclude the intercept if present
  intercept_present <- "Intercept" %in% colnames(coef_matrix)
  if (intercept_present) {
    coef_matrix <- coef_matrix[, !colnames(coef_matrix) %in% "Intercept", drop = FALSE]
  }

  # Calculate the percentage of bootstraps where coefficients are nonzero
  coef_pct <- as.data.frame(colMeans(coef_matrix != 0))

  # Rename the column for clarity
  colnames(coef_pct)[1] <- "Percent of Bootstraps Nonzero"

  # Order the coefficients by descending percentage
  coef_pct <- coef_pct[order(-coef_pct[,1]), , drop = FALSE]

  # Print the resulting data frame
  print(coef_pct)

  # Prepare data for plotting
  coef_pct$Variable <- rownames(coef_pct)
  coef_pct$Order <- seq_len(nrow(coef_pct))  # Create an order variable for the x-axis

  # Create the ggplot2 line and points plot
  plot <- ggplot(coef_pct, aes(x = Order, y = `Percent of Bootstraps Nonzero`)) +
    geom_line(group = 1, color = "steelblue") +  # Connect points with lines
    geom_point(color = "steelblue", size = 3) +  # Add points
    geom_text(aes(label = Variable), vjust = -0.5, size = 3) +  # Add labels above points
    scale_x_continuous(breaks = coef_pct$Order, labels = coef_pct$Variable) +  # Custom x-axis labels
    theme_minimal() +
    labs(
      title = "Coefficient Percentages with Labels",
      x = "Variables",
      y = "Percent of Bootstraps Nonzero"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  # Display the plot
  print(plot)

  # Invisibly return the data frame
  invisible(coef_pct)
}
