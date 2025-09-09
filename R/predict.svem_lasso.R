#' Predict Method for SVEM Models
#'
#' Generates predictions from a fitted svem_model.
#'
#' @param object An object of class svem_model (created by SVEMnet()).
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; default is FALSE. If TRUE, apply the linear calibration
#'   fit (y ~ y_pred) learned during training when available.
#' @param se.fit Logical; if TRUE, returns standard errors (default FALSE).
#' @param agg Aggregation method for ensemble predictions. One of "parms" (default)
#'   or "mean". "parms" uses the aggregated coefficients stored in object$parms
#'   (or parms_debiased if debias=TRUE). "mean" averages predictions from individual
#'   bootstrap members equally and optionally applies the debias calibration.
#' @param ... Additional arguments (currently unused).
#'
#' @importFrom stats terms reformulate na.pass sd delete.response coef model.frame model.matrix
#' @return A numeric vector of predictions, or a list with components \code{fit} and
#'   \code{se.fit} when \code{se.fit = TRUE}.
#'
#' @details
#' The function uses the training terms, factor levels (xlevels), and contrasts
#' saved by SVEMnet(). The terms object environment is set to baseenv() to avoid
#' unexpected lookup of objects in the original environment.
#'
#' Column handling is strict but robust:
#' \itemize{
#'   \item We require that the set of columns produced by \code{model.matrix} on
#'   \code{newdata} matches the set used in training.
#'   \item If \code{model.matrix} drops some training columns (e.g., a factor
#'   collapses to a single level in \code{newdata}), we add those missing columns
#'   as zeros so coefficients align.
#'   \item We then reorder columns to the exact training order before multiplying.
#' }
#'
#' Rows in \code{newdata} that contain unseen factor levels will yield NA predictions
#' (and NA standard errors if \code{se.fit=TRUE}); a warning is issued indicating how
#' many rows were affected.
#'
#' When \code{agg="mean"}, the predictive SE is the bootstrap spread (row-wise sd of
#' member predictions). If \code{debias=TRUE} and the calibration slope is available and
#' finite, SEs are scaled by abs(slope).
#'
#' @examples
#' set.seed(1)
#' n  <- 40
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' y  <- 1 + 0.8*X1 - 0.5*X2 + 0.2*X3 + rnorm(n, 0, 0.3)
#' dat <- data.frame(y, X1, X2, X3)
#' fit <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 30, relaxed = TRUE)
#' pred <- predict(fit, dat)
#' head(pred)
#'
#' # With SEs and debias
#' out <- predict(fit, dat, debias = TRUE, se.fit = TRUE, agg = "mean")
#' str(out)
#'
#' @export
#' @method predict svem_model
predict.svem_model <- function(object, newdata, debias = FALSE, se.fit = FALSE,
                               agg = c("parms", "mean"), ...) {
  agg <- match.arg(agg)

  if (!is.data.frame(newdata)) stop("newdata must be a data frame.")
  if (is.null(object$terms) || is.null(object$parms)) {
    stop("The fitted object is missing required components (terms/parms).")
  }

  # Use training terms but nuke environment to be safe
  terms_obj <- stats::delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  # Harmonize factor columns to training xlevels BEFORE model.frame()
  xlev <- if (!is.null(object$xlevels) && is.list(object$xlevels)) object$xlevels else list()
  if (length(xlev)) {
    for (v in names(xlev)) {
      if (v %in% names(newdata)) {
        if (is.factor(newdata[[v]])) {
          newdata[[v]] <- factor(as.character(newdata[[v]]), levels = xlev[[v]])
        } else if (is.character(newdata[[v]])) {
          newdata[[v]] <- factor(newdata[[v]], levels = xlev[[v]])
        } else {
          newdata[[v]] <- factor(as.character(newdata[[v]]), levels = xlev[[v]])
        }
      }
    }
  }

  # Build model frame/matrix (no xlev arg now; we've already coerced)
  ctr <- if (!is.null(object$contrasts)) object$contrasts else NULL
  mf  <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm  <- stats::model.matrix(terms_obj, data = mf, contrasts.arg = ctr)

  # Remove intercept column
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) mm <- mm[, -int_col, drop = FALSE]

  # Rows with NA columns (typically unseen levels)
  bad_rows <- rowSums(is.na(mm)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows),
            " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm[is.na(mm)] <- 0
  }

  # Determine training column order to mirror
  train_cols <- if (!is.null(object$training_X)) {
    colnames(object$training_X)
  } else if (!is.null(object$schema$feature_names)) {
    object$schema$feature_names
  } else {
    colnames(mm)  # last resort: no reordering possible
  }

  # Build a design matrix in training order (zero-fill missing)
  mm_use <- matrix(0, nrow(mm), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  }

  # Coefficients (includes "(Intercept)")
  coefs <- object$parms
  beta  <- stats::setNames(numeric(length(train_cols) + 1L),
                           c("(Intercept)", train_cols))
  common_beta <- intersect(names(coefs), names(beta))
  beta[common_beta] <- coefs[common_beta]

  # If we don't have member coefs, clamp features
  have_members <- !is.null(object$coef_matrix) && is.matrix(object$coef_matrix)

  if (se.fit || agg == "mean") {
    if (!have_members) {
      stop("se.fit=TRUE or agg='mean' requires bootstrap member predictions (coef_matrix), ",
           "which glmnet_with_cv() objects do not have.")
    }
  }

  # Predict via aggregate coefficients ("parms") — works for both fit types
  preds_parms <- as.numeric(beta[1L] + mm_use %*% beta[-1L])

  # Optionally compute member-based predictions (SVEMnet fits only)
  if (se.fit || agg == "mean") {
    coef_matrix <- object$coef_matrix  # nBoot x (p+1)
    m  <- nrow(mm_use)
    nb <- nrow(coef_matrix)
    intercepts <- coef_matrix[, 1]
    betas      <- coef_matrix[, -1, drop = FALSE]
    pred_mat   <- mm_use %*% t(betas) + matrix(intercepts, nrow = m, ncol = nb, byrow = TRUE)

    # mean aggregation (optionally debias via calibration)
    preds_mean <- rowMeans(pred_mat)
    if (debias && !is.null(object$debias_fit)) {
      ab <- stats::coef(object$debias_fit)
      if (length(ab) >= 2 && is.finite(ab[2])) preds_mean <- ab[1] + ab[2] * preds_mean
    }

    # choose output path
    preds_out <- if (agg == "parms") preds_parms else preds_mean

    if (se.fit) {
      se <- apply(pred_mat, 1, stats::sd)
      if (debias && !is.null(object$debias_fit)) {
        ab <- stats::coef(object$debias_fit)
        if (length(ab) >= 2 && is.finite(ab[2])) se <- abs(ab[2]) * se
      }
      if (any(bad_rows)) { preds_out[bad_rows] <- NA_real_; se[bad_rows] <- NA_real_ }
      return(list(fit = preds_out, se.fit = se))
    } else {
      if (any(bad_rows)) preds_out[bad_rows] <- NA_real_
      return(preds_out)
    }
  } else {
    # parms-only path
    if (debias && !is.null(object$debias_fit)) {
      ok <- !is.na(preds_parms)
      if (any(ok)) {
        preds_parms[ok] <- as.numeric(stats::predict(object$debias_fit,
                                                     newdata = data.frame(y_pred = preds_parms[ok])))
      }
    }
    if (any(bad_rows)) preds_parms[bad_rows] <- NA_real_
    return(preds_parms)
  }
}
