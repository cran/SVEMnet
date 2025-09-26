#' Coefficient Nonzero Percentages from an SVEM Model
#'
#' Calculates the percentage of bootstrap iterations in which each coefficient
#' (excluding the intercept) is nonzero, using a small tolerance.
#'
#' @param object An object of class \code{svem_model}.
#' @param tol Numeric tolerance for "nonzero" (default \code{1e-7}).
#' @param plot Logical; if \code{TRUE}, draws a quick ggplot summary (default \code{TRUE}).
#' @param print_table Logical; if \code{TRUE}, prints a compact table (default \code{TRUE}).
#' @param ... Unused.
#'
#' @return Invisibly returns a data frame with columns:
#'   \code{Variable}, \code{Percent of Bootstraps Nonzero}.
#' @import ggplot2
#' @export
#' @method coef svem_model
coef.svem_model <- function(object, tol = 1e-7, plot = TRUE, print_table = TRUE, ...) {
  if (!is.list(object) || is.null(object$coef_matrix)) {
    stop("The provided model does not contain 'coef_matrix'. Is it a valid 'svem_model'?")
  }
  if (inherits(object, "svem_cv")) {
    stop("coef() for svem_cv objects is not defined (no bootstrap matrix).")
  }


  cm <- object$coef_matrix
  if (!is.matrix(cm)) cm <- as.matrix(cm)

  # Drop rows with any non-finite values to keep colMeans stable
  if (nrow(cm) > 0L) {
    good_rows <- rowSums(is.finite(cm)) == ncol(cm)
    if (any(!good_rows)) cm <- cm[good_rows, , drop = FALSE]
  }

  # Remove intercept column if present
  if (!is.null(colnames(cm)) && "(Intercept)" %in% colnames(cm)) {
    cm <- cm[, setdiff(colnames(cm), "(Intercept)"), drop = FALSE]
  }

  # Nothing to summarize?
  if (ncol(cm) == 0L) {
    message("No non-intercept coefficients found.")
    out <- data.frame(
      Variable = character(0),
      `Percent of Bootstraps Nonzero` = numeric(0),
      check.names = FALSE
    )
    if (isTRUE(print_table)) print(out)
    return(invisible(out))
  }

  # Compute % nonzero per coefficient (tolerant to NAs)
  pct_vec <- 100 * colMeans(abs(cm) > tol, na.rm = TRUE)

  # Build tidy table with stable names (spaces kept; check.names=FALSE)
  out <- data.frame(
    Variable = names(pct_vec),
    `Percent of Bootstraps Nonzero` = as.numeric(pct_vec),
    check.names = FALSE
  )
  out <- out[order(-out$`Percent of Bootstraps Nonzero`), , drop = FALSE]
  rownames(out) <- out$Variable

  # Print a compact one-column view if requested
  if (isTRUE(print_table)) {
    print(out[, "Percent of Bootstraps Nonzero", drop = FALSE])
  }

  # Optional plot (skip silently if ggplot2 is unavailable)
  if (isTRUE(plot)) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      out$Order <- seq_len(nrow(out))
      gp <- ggplot2::ggplot(out, ggplot2::aes(x = Order, y = `Percent of Bootstraps Nonzero`)) +
        ggplot2::geom_line(group = 1) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::geom_text(ggplot2::aes(label = Variable), vjust = -0.6, size = 3) +
        ggplot2::scale_x_continuous(breaks = out$Order, labels = out$Variable) +
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
    } else {
      message("Package 'ggplot2' not available; skipping plot.")
    }
  }

  invisible(out)
}
