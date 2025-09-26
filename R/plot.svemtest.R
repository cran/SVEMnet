#' Plot SVEM Significance Test Results for Multiple Responses
#'
#' Plots the Mahalanobis-like distances for original and permuted data from
#' one or more SVEM significance test results.
#'
#' @param x An object of class \code{svem_significance_test}.
#' @param ... Optional additional \code{svem_significance_test} objects to include.
#' @param labels Optional character vector of labels for the responses.
#'   If not provided, the function uses the response names and ensures uniqueness.
#' @return A \code{ggplot} object showing the distributions of distances.
#' @details
#' If additional \code{svem_significance_test} objects are provided via \code{...},
#' they will be combined into a single plot alongside \code{x}.
#'
#' @import ggplot2
#' @name plot.svem_significance_test
#' @aliases plot_svem_significance_tests
#' @export
#' @method plot svem_significance_test
plot.svem_significance_test <- function(x, ..., labels = NULL) {
  # collect all inputs (x plus ...)
  results_list <- c(list(x), list(...))
  if (!length(results_list)) stop("Provide at least one svem_significance_test object.")
  if (!all(vapply(results_list, function(z) inherits(z, "svem_significance_test"), logical(1)))) {
    stop("All inputs must be objects of class 'svem_significance_test'.")
  }

  # build labels
  if (is.null(labels)) {
    labels <- vapply(seq_along(results_list), function(i) {
      res <- results_list[[i]]
      nm  <- if (!is.null(res$data_d) && length(res$data_d$Response)) as.character(res$data_d$Response[1]) else "Response"
      paste0(nm, "_", i)
    }, character(1))
  } else {
    labels <- make.unique(labels)
    if (length(labels) != length(results_list)) {
      stop("Length of 'labels' must match the number of results.")
    }
  }

  # combine
  pieces <- lapply(seq_along(results_list), function(i) {
    res <- results_list[[i]]
    dd  <- res$data_d
    if (is.null(dd) || !all(c("D","Source_Type","Response") %in% names(dd))) return(NULL)
    dd$Response <- labels[i]
    dd
  })
  plot_data <- do.call(rbind, pieces)
  if (is.null(plot_data) || !nrow(plot_data)) stop("No plottable data found in the provided objects.")

  plot_data$Source_Type <- factor(plot_data$Source_Type, levels = c("Original","Permutation"))
  plot_data$Group <- interaction(plot_data$Response, plot_data$Source_Type, sep = " / ", drop = TRUE)
  plot_data$Group <- factor(plot_data$Group, levels = unique(plot_data$Group))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = D, fill = Source_Type)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.2), alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Distances for Original vs Permutation by Response",
      x = "Response / Source",
      y = "Distance (D)",
      fill = "Source"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  p
}
