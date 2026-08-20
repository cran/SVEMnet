#' Coefficients for SVEM Models
#'
#' Extracts averaged coefficients from an \code{svem_model} fitted by
#' \code{\link{SVEMnet}}.
#'
#' For Gaussian fits, you can optionally request debiased coefficients (if they
#' were computed and stored) via \code{debiased = TRUE}. In that case, the
#' function returns \code{object$parms_debiased}. If debiased coefficients are
#' not available, or if \code{debiased = FALSE}, the function returns
#' \code{object$parms}, which are the ensemble-averaged coefficients across
#' bootstrap members.
#'
#' For Binomial models, \code{debiased} is ignored and the averaged
#' coefficients in \code{object$parms} are returned.
#'
#' @param object An object of class \code{svem_model}, typically returned by
#'   \code{\link{SVEMnet}}.
#' @param debiased Logical; if \code{TRUE} and debiased coefficients are
#'   available (Gaussian fits with \code{parms_debiased}), return those instead
#'   of \code{parms}. Default is \code{FALSE}.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A named numeric vector of coefficients (including the intercept).
#'
#' @details
#' This is a lightweight accessor around the stored components of an
#' \code{svem_model}:
#' \itemize{
#'   \item \code{parms}: ensemble-averaged coefficients over bootstrap members,
#'         on the model's link scale;
#'   \item \code{parms_debiased}: optional debiased coefficients (Gaussian only),
#'         if requested at fit time.
#' }
#' Passing \code{debiased = TRUE} has no effect if \code{parms_debiased} is
#' \code{NULL}.
#'
#' @seealso \code{\link{svem_nonzero}} for bootstrap nonzero percentages and a
#'   quick stability plot.
#'
#' @examples
#'   set.seed(1)
#'   n  <- 50
#'   x1 <- rnorm(n)
#'   x2 <- rnorm(n)
#'   eps <- rnorm(n, sd = 0.3)
#'   y_g <- 1 + 2*x1 - 0.5*x2 + eps
#'   dat_g <- data.frame(y_g, x1, x2)
#'
#'   fit_g <- SVEMnet(y_g ~ x1 + x2, data = dat_g,
#'                    nBoot = 10, glmnet_alpha = 1, relaxed = FALSE)
#'
#'   # Ensemble-averaged coefficients
#'   cc <- coef(fit_g)
#'   head(cc)
#'
#'   # Debiased (only if available for Gaussian fits)
#'   ccd <- coef(fit_g, debiased = TRUE)
#'   head(ccd)
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
#' Summarizes variable-selection stability across SVEM bootstrap refits by
#' computing the percentage of bootstrap iterations in which each coefficient
#' (excluding the intercept) is nonzero, using a small tolerance. Optionally
#' produces a quick \pkg{ggplot2} summary and/or prints a compact table.
#'
#' This function expects \code{object$coef_matrix} to contain the per-bootstrap
#' coefficients (including an intercept column), which \code{\link{SVEMnet}}
#' always stores in its fitted objects.
#' Rows correspond to bootstrap fits; columns correspond to coefficients.
#'
#' @param object An object of class \code{svem_model} with a non-empty
#'   \code{$coef_matrix} component. \code{svem_nonzero()} is not defined for
#'   \code{svem_cv} objects.
#' @param tol Numeric tolerance for "nonzero". Coefficients with
#'   \code{|beta| > tol} are counted as nonzero. Default is \code{1e-7}.
#' @param plot Logical; if \code{TRUE}, draws a quick \pkg{ggplot2} summary plot
#'   of the nonzero percentages (default \code{TRUE}).
#' @param print_table Logical; if \code{TRUE}, prints a compact table of
#'   nonzero percentages to the console (default \code{TRUE}).
#' @param ... Unused; included for future extension.
#'
#' @return
#' Invisibly returns a data frame with columns:
#' \itemize{
#'   \item \code{Variable}: coefficient name (excluding the intercept).
#'   \item \code{Percent of Bootstraps Nonzero}: percentage (0–100) of bootstrap
#'         fits in which \code{|beta| > tol}.
#' }
#'
#' If no non-intercept coefficients are found (for example, if only the
#' intercept is present), an empty data frame is returned and a message is
#' issued.
#'
#' @details
#' Internally, \code{svem_nonzero()}:
#' \itemize{
#'   \item checks for and drops rows of \code{coef_matrix} that contain any
#'         non-finite values, to keep summaries stable;
#'   \item drops an \code{"(Intercept)"} column if present;
#'   \item computes \code{100 * mean(|beta_j| > tol)} across bootstrap rows
#'         for each remaining coefficient.
#' }
#'
#' The plot is a simple line + point chart with labels, ordered by decreasing
#' nonzero percentage. It is intended as a quick diagnostic; for publication
#' graphics, you may want to customize the output data frame with your own
#' plotting code.
#'
#' @seealso \code{\link[=coef.svem_model]{coef.svem_model}} for ensemble-averaged
#'   (optionally debiased) coefficients.
#'
#' @examples
#'   set.seed(10)
#'   n  <- 60
#'   x1 <- rnorm(n)
#'   x2 <- rnorm(n)
#'   x3 <- rnorm(n)
#'   eps <- rnorm(n, sd = 0.4)
#'   y   <- 0.7 + 1.5*x1 - 0.8*x2 + 0.05*x3 + eps
#'   dat <- data.frame(y, x1, x2, x3)
#'
#'   fit <- SVEMnet(y ~ x1 + x2 + x3, data = dat,
#'                  nBoot = 10, glmnet_alpha = 1, relaxed = FALSE)
#'
#'   nz <- svem_nonzero(fit, tol = 1e-7, plot = FALSE, print_table = TRUE)
#'   head(nz)
#'
#' @import ggplot2
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

    gp <- ggplot2::ggplot(out, ggplot2::aes(x = Order, y = PercentNonzero)) +
      ggplot2::geom_line(ggplot2::aes(group = 1)) +
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
  }

  invisible(out)
}
