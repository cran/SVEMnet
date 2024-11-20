#' Plot SVEM Significance Test Results for Multiple Responses
#'
#' Plots the Mahalanobis distances for the original and permuted data from multiple SVEM significance test results.
#'
#' @param ... One or more objects of class \code{svem_significance_test}, which are the outputs from \code{\link{svem_significance_test}}.
#' @param labels Optional character vector of labels for the responses. If not provided, the function uses the response variable names.
#' @return A \code{ggplot} object showing the distributions of Mahalanobis distances for all responses.
#' @details
#' This function creates a combined plot of the Mahalanobis distances (\code{d_Y} and \code{d_pi_Y}) for the original and permuted data from multiple SVEM significance test results. It groups the data by response and source type, displaying original and permutation distances side by side for each response.
#'
#' **Usage Notes:**
#' \itemize{
#'   \item Use this function to compare the significance test results across multiple responses.
#'   \item The plot shows original and permutation distances next to each other for each response.
#' }
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @import ggplot2
#' @export
plot.svem_significance_test <- function(..., labels = NULL) {
  # Collect the result objects into a list
  results_list <- list(...)

  # Check that all inputs are of class 'svem_significance_test'
  if (!all(sapply(results_list, function(x) inherits(x, "svem_significance_test")))) {
    stop("All inputs must be objects of class 'svem_significance_test'.")
  }

  # If labels are provided, ensure they are unique; otherwise, generate unique labels
  if (is.null(labels)) {
    labels <- sapply(seq_along(results_list), function(i) {
      res <- results_list[[i]]
      if (!is.null(res$data_d$Response[1])) {
        response_name <- as.character(res$data_d$Response[1])
        # Append index to response name to ensure uniqueness
        paste0(response_name, "_", i)
      } else {
        paste0("Response_", i)
      }
    })
  } else {
    # Ensure labels are unique
    labels <- make.unique(labels)
  }

  # Combine data from all results
  plot_data <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    res <- results_list[[i]]
    if (!is.null(res$data_d)) {
      data_d <- res$data_d
      data_d$Response <- labels[i]  # Update the Response column with unique labels
      data_d
    } else {
      NULL
    }
  }))

  # Create a combined grouping variable
  plot_data$Group <- paste(plot_data$Response, plot_data$Source_Type, sep = "_")

  # Adjust the levels to control the order on the x-axis
  unique_groups <- unique(plot_data$Group)
  plot_data$Group <- factor(plot_data$Group, levels = unique_groups)

  # Create the plot with coord_cartesian to include y = 0
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = D, fill = Source_Type)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.2), alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Mahalanobis Distances for Original and Permuted Data",
      x = "Response and Source Type",
      y = "Mahalanobis Distance (D)",
      fill = "Source Type"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +  # Ensures y-axis includes 0 without removing data
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # Print the plot
  print(p)
}
