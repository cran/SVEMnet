  #' Coefficients for SVEM Models
  #'
  #' Returns averaged coefficients from an \code{svem_model}.
  #'
  #' For Gaussian models, you can optionally return the debiased coefficients
  #' (if available) via \code{debiased = TRUE}. For Binomial models, \code{debiased}
  #' is ignored and the averaged coefficients are returned.
  #'
  #' @param object An object of class \code{svem_model}.
  #' @param debiased Logical; if \code{TRUE} and available (Gaussian fits),
  #'   return \code{parms_debiased} instead of \code{parms}. Default \code{FALSE}.
  #' @param ... Unused (for S3 compatibility).
  #'
  #' @return A named numeric vector of coefficients (including intercept).
  #'
  #' @seealso \code{\link{svem_nonzero}} for bootstrap nonzero percentages and a quick plot.
  #'
  #' @examples
  #' \donttest{
  #'   set.seed(1)
  #'   n  <- 200
  #'   x1 <- rnorm(n)
  #'   x2 <- rnorm(n)
  #'   eps <- rnorm(n, sd = 0.3)
  #'   y_g <- 1 + 2*x1 - 0.5*x2 + eps
  #'   dat_g <- data.frame(y_g, x1, x2)
  #'
  #'   # Small nBoot to keep runtime light in examples
  #'   fit_g <- SVEMnet(y_g ~ x1 + x2, data = dat_g, nBoot = 30, relaxed = TRUE)
  #'
  #'   # Averaged coefficients
  #'   cc <- coef(fit_g)
  #'   head(cc)
  #'
  #'   # Debiased (only if available for Gaussian fits)
  #'   ccd <- coef(fit_g, debiased = TRUE)
  #'   head(ccd)
  #'
  #'   # Binomial example (0/1 outcome)
  #'   set.seed(2)
  #'   n  <- 250
  #'   x1 <- rnorm(n)
  #'   x2 <- rnorm(n)
  #'   eta <- -0.4 + 1.1*x1 - 0.7*x2
  #'   p   <- 1/(1+exp(-eta))
  #'   y_b <- rbinom(n, 1, p)
  #'   dat_b <- data.frame(y_b, x1, x2)
  #'
  #'   fit_b <- SVEMnet(y_b ~ x1 + x2, data = dat_b, family = "binomial",
  #'                    nBoot = 30, relaxed = TRUE)
  #'
  #'   # Averaged coefficients (binomial)
  #'   coef(fit_b)
  #' }
  #'
  #' @export
  #' @method coef svem_model
  coef.svem_model <- function(object, debiased = FALSE, ...) {
    stopifnot(is.list(object), !is.null(object$parms))
    if (isTRUE(debiased) && !is.null(object$parms_debiased)) {
      return(object$parms_debiased)
    }
    object$parms
  }

#' Coefficient Nonzero Percentages (SVEM)
#'
#' Calculates the percentage of bootstrap iterations in which each coefficient
#' (excluding the intercept) is nonzero, using a small tolerance. Optionally
#' draws a quick \pkg{ggplot2} summary and/or prints a compact table.
#'
#' This function summarizes variable selection stability across SVEM bootstrap
#' refits. It expects \code{object$coef_matrix} to contain the per-bootstrap
#' coefficients (including an intercept column).
#'
#' @param object An object of class \code{svem_model}.
#' @param tol Numeric tolerance for "nonzero" (default \code{1e-7}).
#' @param plot Logical; if \code{TRUE}, draws a quick ggplot summary (default \code{TRUE}).
#' @param print_table Logical; if \code{TRUE}, prints a compact table (default \code{TRUE}).
#' @param ... Unused.
#'
#' @return Invisibly returns a data frame with columns:
#'   \itemize{
#'     \item \code{Variable}
#'     \item \code{Percent of Bootstraps Nonzero}
#'   }
#'
#' @seealso \code{\link[=coef.svem_model]{coef.svem_model}} for averaged (optionally debiased) coefficients.
#'
#' @examples
#' \donttest{
#'   ## ---------- Gaussian demo ----------
#'   set.seed(10)
#'   n  <- 220
#'   x1 <- rnorm(n)
#'   x2 <- rnorm(n)
#'   x3 <- rnorm(n)
#'   eps <- rnorm(n, sd = 0.4)
#'   y   <- 0.7 + 1.5*x1 - 0.8*x2 + 0.05*x3 + eps
#'   dat <- data.frame(y, x1, x2, x3)
#'
#'   fit <- SVEMnet(y ~ (x1 + x2 + x3)^2, data = dat, nBoot = 40, relaxed = TRUE)
#'
#'   # Table + plot
#'   nz <- svem_nonzero(fit, tol = 1e-7, plot = TRUE, print_table = TRUE)
#'   head(nz)
#'
#'   ## ---------- Binomial demo ----------
#'   set.seed(11)
#'   n  <- 260
#'   x1 <- rnorm(n)
#'   x2 <- rnorm(n)
#'   x3 <- rnorm(n)
#'   lp <- -0.3 + 0.9*x1 - 0.6*x2 + 0.2*x3
#'   p  <- 1/(1+exp(-lp))
#'   y  <- rbinom(n, 1, p)
#'   dat_b <- data.frame(y, x1, x2, x3)
#'
#'   fit_b <- SVEMnet(y ~ x1 + x2 + x3, data = dat_b,
#'                    family = "binomial", nBoot = 40, relaxed = TRUE)
#'
#'   # Plot optional; still summarizes the bootstrap selection frequencies
#'   svem_nonzero(fit_b, plot = TRUE, print_table = TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
svem_nonzero <- function(object, tol = 1e-7, plot = TRUE, print_table = TRUE, ...) {
  if (!is.list(object) || is.null(object$coef_matrix)) {
    stop("The provided model does not contain 'coef_matrix'. Is it a valid 'svem_model'?")
  }
  if (inherits(object, "svem_cv")) {
    stop("svem_nonzero() is not defined for 'svem_cv' objects (no bootstrap matrix).")
  }

  cm <- object$coef_matrix
  if (!is.matrix(cm)) cm <- as.matrix(cm)

  # Drop rows with any non-finite values to keep summaries stable
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

  # Compute % nonzero per coefficient
  pct_vec <- 100 * colMeans(abs(cm) > tol, na.rm = TRUE)

  out <- data.frame(
    Variable = names(pct_vec),
    `Percent of Bootstraps Nonzero` = as.numeric(pct_vec),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$`Percent of Bootstraps Nonzero`), , drop = FALSE]
  rownames(out) <- out$Variable

  if (isTRUE(print_table)) {
    print(out[, "Percent of Bootstraps Nonzero", drop = FALSE])
  }

  if (isTRUE(plot)) {
    # Safe helper columns (avoid spaces in aes mappings)
    out$Order          <- seq_len(nrow(out))
    out$PercentNonzero <- out[["Percent of Bootstraps Nonzero"]]

    gp <- ggplot2::ggplot(out, ggplot2::aes(x = .data$Order, y = .data$PercentNonzero)) +
      ggplot2::geom_line(ggplot2::aes(group = 1)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Variable), vjust = -0.6, size = 3) +
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
  }

  invisible(out)
}
