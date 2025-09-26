
#' Percentile Confidence Intervals for SVEM Predictions
#'
#' Computes predictions and percentile confidence intervals from bootstrap
#' member predictions for a fitted svem_model.
#'
#' @param object A svem_model.
#' @param newdata Data frame of predictors.
#' @param level Confidence level (default 0.95).
#' @param debias Logical; if TRUE, apply calibration to each member prediction (default FALSE).
#' @return A data.frame with columns fit, lwr, upr.
#' @export
predict_with_ci <- function(object, newdata, level = 0.95, debias = FALSE) {
  if (!inherits(object, "svem_model")) stop("object must be a svem_model.")
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.")
  if (is.null(object$coef_matrix) || !is.matrix(object$coef_matrix)) {
    stop("Percentile CI requires stored bootstrap member coefficients (coef_matrix).")
  }
  # Build design like predict.svem_model does
  terms_obj <- stats::delete.response(object$terms)
  environment(terms_obj) <- baseenv()
  xlev <- if (!is.null(object$xlevels) && is.list(object$xlevels)) object$xlevels else list()
  nd <- newdata
  if (length(xlev)) {
    for (v in names(xlev)) {
      if (v %in% names(nd)) {
        if (is.factor(nd[[v]])) {
          nd[[v]] <- factor(as.character(nd[[v]]), levels = xlev[[v]])
        } else if (is.character(nd[[v]])) {
          nd[[v]] <- factor(nd[[v]], levels = xlev[[v]])
        } else {
          nd[[v]] <- factor(as.character(nd[[v]]), levels = xlev[[v]])
        }
      }
    }
  }
  ctr <- if (!is.null(object$contrasts)) object$contrasts else NULL
  mf  <- stats::model.frame(terms_obj, data = nd, na.action = stats::na.pass)
  mm  <- stats::model.matrix(terms_obj, data = mf, contrasts.arg = ctr)
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) mm <- mm[, -int_col, drop = FALSE]
  bad_rows <- rowSums(is.na(mm)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows), " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm[is.na(mm)] <- 0
  }
  train_cols <- if (!is.null(object$training_X)) {
    colnames(object$training_X)
  } else if (!is.null(object$schema$feature_names)) {
    object$schema$feature_names
  } else {
    colnames(mm)
  }
  mm_use <- matrix(0, nrow(mm), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  }

  coef_matrix <- object$coef_matrix
  m  <- nrow(mm_use)
  nb <- nrow(coef_matrix)
  intercepts <- coef_matrix[, 1]
  betas      <- coef_matrix[, -1, drop = FALSE]
  pred_mat   <- mm_use %*% t(betas) + matrix(intercepts, nrow = m, ncol = nb, byrow = TRUE)

  if (debias && !is.null(object$debias_fit)) {
    ab <- stats::coef(object$debias_fit)
    if (length(ab) >= 2 && is.finite(ab[2])) {
      pred_mat <- ab[1] + ab[2] * pred_mat
    }
  }
  fit <- rowMeans(pred_mat)
  alpha <- (1 - level) / 2
  lwr <- apply(pred_mat, 1, stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE, type = 7)
  upr <- apply(pred_mat, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE, type = 7)
  if (any(bad_rows)) {
    fit[bad_rows] <- NA_real_; lwr[bad_rows] <- NA_real_; upr[bad_rows] <- NA_real_
  }
  data.frame(fit = as.numeric(fit), lwr = as.numeric(lwr), upr = as.numeric(upr))
}
