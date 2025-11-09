#' Plot Method for SVEM Models (Gaussian / Generic)
#'
#' Plots actual versus predicted values for an \code{svem_model}. This is the
#' default plot for models fit with \code{SVEMnet(..., family = "gaussian")} and
#' any other non-binomial models that share the \code{svem_model} class.
#'
#' Points show fitted values against observed responses; the dashed line is the
#' 45-degree identity. If available and requested, debiased predictions are included
#' as a second series.
#'
#' @param x An object of class \code{svem_model}.
#' @param plot_debiased Logical; if \code{TRUE}, include debiased predictions
#'   (when available) as an additional series. Default \code{FALSE}.
#' @param ... Additional aesthetics passed to \code{ggplot2::geom_point()}.
#' @details
#' This method assumes the fitted object stores the training response in
#' \code{$actual_y} and in-sample predictions in \code{$y_pred}, as produced
#' by \code{SVEMnet()} and \code{glmnet_with_cv()}.

#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#'   ## --- Gaussian example: simulate, fit, and plot --------------------------
#'   set.seed(2026)
#'   n  <- 300
#'   X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#'   eps <- rnorm(n, sd = 0.4)
#'   y_g <- 1.2 + 2*X1 - 0.7*X2 + 0.3*X3 + 1.1*(X1*X2) + 0.8*(X1^2) + eps
#'   dat_g <- data.frame(y_g, X1, X2, X3)
#'
#'   fit_g <- SVEMnet(
#'     y_g ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2),
#'     data          = dat_g,
#'     family        = "gaussian",
#'     glmnet_alpha  = c(1, 0.5),
#'     nBoot         = 60,
#'     objective     = "auto",
#'     weight_scheme = "SVEM",
#'     relaxed       = TRUE
#'   )
#'
#'   # Actual vs predicted (with and without debias overlay)
#'   plot(fit_g, plot_debiased = FALSE)
#'   plot(fit_g, plot_debiased = TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom rlang .data
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

  if (isTRUE(plot_debiased) && !is.null(x$y_pred_debiased)) {
    plot_data <- rbind(
      plot_data,
      data.frame(
        Actual    = actual_y,
        Predicted = x$y_pred_debiased,
        Type      = "Debiased Predictions"
      )
    )
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Actual, y = .data$Predicted,
                 color = .data$Type, shape = .data$Type)
  ) +
    ggplot2::geom_point(...) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         color = "black", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = paste("Actual vs Predicted -", utils::tail(class(x), 1)),
      x = "Actual y",
      y = "Predicted y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c(
      "Predictions" = "blue",
      "Debiased Predictions" = "red"
    )) +
    ggplot2::scale_shape_manual(values = c(
      "Predictions" = 16,
      "Debiased Predictions" = 17
    )) +
    ggplot2::coord_fixed()
}

#' Plot Method for SVEM Binomial Models
#'
#' Diagnostics for \code{svem_binomial} fits from \code{SVEMnet(..., family = "binomial")}.
#' Produces one of:
#' \itemize{
#'   \item \code{type = "calibration"}: Reliability curve (binned average predicted probability vs observed rate),
#'         with jittered raw points for context.
#'   \item \code{type = "roc"}: ROC curve with trapezoidal AUC in the title.
#'   \item \code{type = "pr"}: Precision–Recall curve with step-wise Average Precision (AP).
#' }
#'
#' For ROC/PR, simple one-class guards are used (returns a diagonal ROC and trivial PR).
#' The function assumes binomial models store \code{x$y_pred} on the \emph{probability} scale.
#'
#' @param x An object of class \code{svem_binomial}.
#' @param type One of \code{"calibration"}, \code{"roc"}, or \code{"pr"} (default \code{"calibration"}).
#' @param bins Integer number of equal-frequency bins for calibration (default \code{10}).
#' @param jitter_width Vertical jitter amplitude for raw points in calibration (default \code{0.05}).
#' @param ... Additional aesthetics passed to \code{ggplot2::geom_line()} or \code{ggplot2::geom_point()}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#'   ## --- Binomial example: simulate, fit, and plot --------------------------
#'   set.seed(2025)
#'   n  <- 600
#'   x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
#'   eta    <- -0.4 + 1.1*x1 - 0.8*x2 + 0.5*x3
#'   p_true <- plogis(eta)
#'   y      <- rbinom(n, 1, p_true)
#'   dat_b  <- data.frame(y, x1, x2, x3)
#'
#'   fit_b <- SVEMnet(
#'     y ~ x1 + x2 + x3 + I(x1^2) + (x1 + x2 + x3)^2,
#'     data          = dat_b,
#'     family        = "binomial",
#'     glmnet_alpha  = c(1, 0.5),
#'     nBoot         = 60,
#'     objective     = "auto",
#'     weight_scheme = "SVEM",
#'     relaxed       = TRUE
#'   )
#'
#'   # Calibration / ROC / PR
#'   plot(fit_b, type = "calibration", bins = 12)
#'   plot(fit_b, type = "roc")
#'   plot(fit_b, type = "pr")
#' }
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stats quantile aggregate runif
#' @export
#' @method plot svem_binomial
plot.svem_binomial <- function(x,
                               type = c("calibration", "roc", "pr"),
                               bins = 10,
                               jitter_width = 0.05,
                               ...) {
  type <- match.arg(type)
  stopifnot(is.numeric(x$actual_y), all(x$actual_y %in% c(0, 1)))
  y <- as.integer(x$actual_y)
  p <- as.numeric(x$y_pred)   # stored on probability scale for binomial in SVEMnet()
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]

  # minimal helpers (no extra packages)
  .roc_curve <- function(prob, yy) {
    o <- order(prob, decreasing = TRUE)
    yy <- yy[o]
    tp <- cumsum(yy == 1)
    fp <- cumsum(yy == 0)
    P  <- sum(yy == 1); N <- sum(yy == 0)
    if (P == 0 || N == 0) return(data.frame(fpr = c(0, 1), tpr = c(0, 1)))
    data.frame(fpr = fp / N, tpr = tp / P)
  }
  .auc_trap <- function(curve) {
    x <- curve$fpr; y <- curve$tpr
    o <- order(x, y)
    x <- x[o]; y <- y[o]
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }
  .pr_curve <- function(prob, yy) {
    o <- order(prob, decreasing = TRUE)
    yy <- yy[o]
    tp <- cumsum(yy == 1)
    fp <- cumsum(yy == 0)
    P  <- sum(yy == 1)
    precision <- tp / pmax(1, tp + fp)
    recall    <- if (P == 0) rep(0, length(tp)) else tp / P
    data.frame(recall = c(0, recall), precision = c(1, precision))
  }
  .ap_step <- function(curve) {
    r <- curve$recall; p <- curve$precision
    sum(diff(r) * utils::tail(p, -1))
  }

  if (type == "calibration") {
    # Equal-frequency bin edges; make sure endpoints cover [min(p), max(p)]
    q <- stats::quantile(p, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE, type = 8)
    q[1] <- min(q[1], min(p))
    q[length(q)] <- max(q[length(q)], max(p))
    # Remove duplicate edges if p has many ties
    q <- unique(as.numeric(q))
    if (length(q) < 2L) {
      # fallback to pretty breakpoints on [0,1]
      q <- seq(0, 1, length.out = max(2, bins + 1))
    }

    b <- cut(p, breaks = q, include.lowest = TRUE, right = TRUE, dig.lab = 12)

    # Aggregate means & bin sizes; drop empty/NA bins before plotting
    agg <- stats::aggregate(
      x = data.frame(p = p, y = y),
      by = list(bin = b),
      FUN = function(z) c(mean = mean(z), n = sum(is.finite(z)))
    )
    p_mean <- sapply(agg$p, `[`, 1)
    y_mean <- sapply(agg$y, `[`, 1)
    n_bin  <- sapply(agg$p, `[`, 2)
    calib  <- data.frame(p_mean = p_mean, y_mean = y_mean, n = n_bin)
    calib  <- calib[is.finite(calib$p_mean) & is.finite(calib$y_mean) & calib$n > 0, , drop = FALSE]

    # Jitter raw points in y, then squish into [0,1] to avoid "Removed rows" warnings
    pts <- data.frame(p = p, y = y)
    if (jitter_width > 0) {
      jy <- y + stats::runif(length(y), min = -jitter_width, max = jitter_width)
      jy <- pmin(1, pmax(0, jy)) # squish to [0,1]
      pts$y <- jy
    }

    plt <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x = .data$p, y = .data$y),
        alpha = 0.25, size = 1.2, na.rm = TRUE, ...
      ) +
      ggplot2::geom_line(
        data = calib,
        ggplot2::aes(x = .data$p_mean, y = .data$y_mean),
        linewidth = 1, na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = calib,
        ggplot2::aes(x = .data$p_mean, y = .data$y_mean, size = .data$n),
        alpha = 0.85, na.rm = TRUE
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = "Calibration plot (SVEM binomial)",
        x = "Predicted probability",
        y = "Observed event rate",
        size = "Bin size"
      ) +
      ggplot2::theme_minimal()

    return(plt)
  }

  if (type == "roc") {
    roc <- .roc_curve(p, y)
    auc <- .auc_trap(roc)
    plt <- ggplot2::ggplot(roc, ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
      ggplot2::geom_line(...) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = sprintf("ROC (AUC = %.3f)", auc),
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal()
    return(plt)
  }

  if (type == "pr") {
    pr <- .pr_curve(p, y)
    ap <- .ap_step(pr)
    plt <- ggplot2::ggplot(pr, ggplot2::aes(x = .data$recall, y = .data$precision)) +
      ggplot2::geom_line(...) +
      ggplot2::labs(
        title = sprintf("Precision-Recall (AP = %.3f)", ap),
        x = "Recall",
        y = "Precision"
      ) +
      ggplot2::theme_minimal()
    return(plt)
  }
}
