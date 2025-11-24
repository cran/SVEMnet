#' Plot SVEM significance test results for one or more responses
#'
#' Plots the Mahalanobis-like distances for original and permuted data from
#' one or more SVEM significance test results returned by
#' \code{svem_significance_test_parallel()}.
#'
#' @param x An object of class \code{svem_significance_test}.
#' @param ... Optional additional \code{svem_significance_test} objects to include
#'   in the same plot.
#' @param labels Optional character vector of labels for the responses. If not
#'   provided, the function uses inferred response names (from \code{data_d$Response}
#'   or \code{x$response}) and ensures uniqueness.
#'
#' @return A \code{ggplot2} object showing the distributions of distances for
#'   original vs. permuted data, grouped by response.
#'
#' @details
#' If additional \code{svem_significance_test} objects are provided via \code{...},
#' their distance tables (\code{$data_d}) are stacked and plotted together
#' using a shared x-axis grouping of \code{"Response / Source"} and a fill
#' aesthetic indicating \code{"Original"} vs \code{"Permutation"}.
#'
#' @import ggplot2
#' @name plot.svem_significance_test
#' @export
#' @method plot svem_significance_test
plot.svem_significance_test <- function(x, ..., labels = NULL) {
  # Collect all inputs (x plus ...)
  results_list <- c(list(x), list(...))
  if (!length(results_list)) stop("Provide at least one svem_significance_test object.")
  if (!all(vapply(results_list, function(z) inherits(z, "svem_significance_test"), logical(1)))) {
    stop("All inputs must be objects of class 'svem_significance_test'.")
  }

  # Build labels
  if (is.null(labels)) {
    labels <- vapply(seq_along(results_list), function(i) {
      res <- results_list[[i]]
      nm <- tryCatch({
        if (!is.null(res$data_d) && "Response" %in% names(res$data_d) && length(res$data_d$Response)) {
          as.character(res$data_d$Response[1])
        } else if (!is.null(res$response)) {
          as.character(res$response)
        } else {
          "Response"
        }
      }, error = function(e) "Response")
      paste0(nm, "_", i)
    }, character(1))
  } else {
    labels <- make.unique(labels)
    if (length(labels) != length(results_list)) {
      stop("Length of 'labels' must match the number of results.")
    }
  }

  # Combine plottable pieces
  pieces <- lapply(seq_along(results_list), function(i) {
    res <- results_list[[i]]
    dd  <- res$data_d
    if (is.null(dd) || !all(c("D", "Source_Type", "Response") %in% names(dd))) return(NULL)

    out <- dd[, c("D", "Source_Type", "Response"), drop = FALSE]
    out$D <- as.numeric(out$D)
    out$Source_Type <- as.character(out$Source_Type)
    out$Response <- labels[i]
    out
  })

  plot_data <- do.call(rbind, pieces)
  if (is.null(plot_data) || !nrow(plot_data)) stop("No plottable data found in the provided objects.")

  # Clean types / factors and drop non-finite distances
  plot_data <- plot_data[is.finite(plot_data$D), , drop = FALSE]
  plot_data$Source_Type <- factor(plot_data$Source_Type, levels = c("Original", "Permutation"))
  plot_data$Response    <- factor(plot_data$Response, levels = unique(plot_data$Response))
  plot_data$Group <- interaction(plot_data$Response, plot_data$Source_Type, sep = " / ", drop = TRUE)
  plot_data$Group <- factor(plot_data$Group, levels = unique(plot_data$Group))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = D, fill = Source_Type)) +
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
}
