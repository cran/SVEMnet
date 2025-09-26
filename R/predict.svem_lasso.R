#' Predict Method for SVEM Models
#'
#' Generates predictions from a fitted svem_model, with optional bootstrap
#' standard errors and percentile confidence intervals.
#'
#' @param object An object of class svem_model (created by SVEMnet()).
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; default is FALSE. If TRUE, apply the linear calibration
#'   fit (y ~ y_pred) learned during training when available.
#' @param se.fit Logical; if TRUE, returns standard errors based on the bootstrap
#'   spread when member predictions are available (default FALSE).
#' @param interval Logical; if TRUE, returns percentile confidence limits based
#'   on bootstrap predictions when available (default FALSE).
#' @param level Confidence level for the percentile interval. Default 0.95.
#' @param agg Aggregation method for ensemble predictions. One of "parms" (default)
#'   or "mean". "parms" uses the aggregated coefficients stored in object$parms
#'   (or parms_debiased if debias=TRUE). "mean" averages predictions from individual
#'   bootstrap members equally and optionally applies the debias calibration.
#' @param ... Additional arguments (currently unused).
#'
#' @importFrom stats terms reformulate na.pass sd delete.response coef model.frame model.matrix quantile
#' @return A numeric vector of predictions, or a list with components as follows:
#'   - fit: predictions
#'   - se.fit: bootstrap standard errors when se.fit = TRUE
#'   - lwr, upr: percentile confidence limits when interval = TRUE
#'
#' @details
#' The function uses the training terms, factor levels (xlevels), and contrasts
#' saved by SVEMnet(). The terms object environment is set to baseenv() to avoid
#' unexpected lookup of objects in the original environment.
#'
#' Column handling follows these rules:
#'   - The set of columns produced by model.matrix on newdata is aligned to those
#'     used in training.
#'   - Any training columns dropped by model.matrix are added back as zeros.
#'   - Columns are reordered to the exact training order before multiplication.
#'
#' Rows in newdata that contain unseen factor levels will yield NA predictions
#' (and NA standard errors and intervals if requested); a warning is issued indicating
#' how many rows were affected.
#'
#' When agg = "mean", se.fit is the row-wise sd of member predictions. If debias = TRUE
#' and a finite calibration slope is available, both member predictions and their sd
#' are transformed by the calibration before aggregation. Percentile intervals are
#' computed from the transformed member predictions.
#'
#' @examples
#' set.seed(1)
#' n  <- 40
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' y  <- 1 + 0.8*X1 - 0.5*X2 + 0.2*X3 + rnorm(n, 0, 0.3)
#' dat <- data.frame(y, X1, X2, X3)
#' fit <- SVEMnet(y ~ (X1 + X2 + X3)^2, dat, nBoot = 30, relaxed = TRUE)
#'
#' # Mean aggregation with SEs and 95 percent CIs
#' out <- predict(fit, dat, debias = TRUE, se.fit = TRUE, interval = TRUE, agg = "mean")
#' str(out)
#'
#' @export
#' @method predict svem_model
predict.svem_model <- function(object, newdata, debias = FALSE, se.fit = FALSE,
                               interval = FALSE, level = 0.95,
                               agg = c("parms", "mean"), ...) {
  agg <- match.arg(agg)

  if (!is.data.frame(newdata)) stop("newdata must be a data frame.")
  if (is.null(object$terms) || is.null(object$parms)) {
    stop("The fitted object is missing required components (terms/parms).")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0,1).")
  }

  # Use training terms but set a safe environment
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

  # Build model frame/matrix
  ctr <- if (!is.null(object$contrasts)) object$contrasts else NULL
  if (!is.null(ctr) && !is.list(ctr)) ctr <- as.list(ctr)
  mf  <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm  <- stats::model.matrix(terms_obj, data = mf, contrasts.arg = ctr)

  # Remove intercept column
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) mm <- mm[, -int_col, drop = FALSE]

  # Rows with non-finite entries (typically unseen levels producing NA)
  bad_rows <- rowSums(!is.finite(mm)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows),
            " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm[!is.finite(mm)] <- 0
  }

  # Determine training column order to mirror
  train_cols <- if (!is.null(object$training_X)) {
    colnames(object$training_X)
  } else if (!is.null(object$schema$feature_names)) {
    object$schema$feature_names
  } else {
    colnames(mm)  # last resort
  }

  # Informative message if newdata had extra columns not seen in training
  extra_cols <- setdiff(colnames(mm), train_cols)
  if (length(extra_cols)) {
    message("Ignoring ", length(extra_cols), " column(s) not present at training: ",
            paste(utils::head(extra_cols, 5L), collapse = ", "),
            if (length(extra_cols) > 5L) " ..." else "")
  }

  # Build a design matrix in training order (zero-fill missing)
  mm_use <- matrix(0, nrow(mm), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  }
  storage.mode(mm_use) <- "double"

  # Aggregate-coefficient prediction ("parms")
  # Use parms_debiased when debias=TRUE and available, else parms
  coefs <- if (debias && !is.null(object$parms_debiased)) object$parms_debiased else object$parms
  beta  <- stats::setNames(numeric(length(train_cols) + 1L),
                           c("(Intercept)", train_cols))
  common_beta <- intersect(names(coefs), names(beta))
  beta[common_beta] <- coefs[common_beta]
  preds_parms <- as.numeric(beta[1L] + mm_use %*% beta[-1L])

  # Member predictions available?
  have_members <- !is.null(object$coef_matrix) && is.matrix(object$coef_matrix)
  if ((se.fit || interval || agg == "mean") && !have_members) {
    stop("se.fit = TRUE, interval = TRUE, or agg = 'mean' requires bootstrap member predictions (coef_matrix). ",
         "Re-fit with SVEMnet(..., nBoot >= 2) to populate coef_matrix.")
  }

  # If we need member predictions, compute them
  if (se.fit || interval || agg == "mean") {
    coef_matrix <- object$coef_matrix  # nBoot x (p+1)
    m  <- nrow(mm_use)
    nb <- nrow(coef_matrix)
    intercepts <- coef_matrix[, 1]
    betas      <- coef_matrix[, -1, drop = FALSE]
    pred_mat   <- mm_use %*% t(betas) + matrix(intercepts, nrow = m, ncol = nb, byrow = TRUE)

    # If debiasing, apply to each member, not just the mean
    if (debias && !is.null(object$debias_fit)) {
      ab <- stats::coef(object$debias_fit)
      if (length(ab) >= 2 && all(is.finite(ab[1:2]))) {
        pred_mat <- ab[1] + ab[2] * pred_mat
      }
    }

    # Mean aggregation
    preds_mean <- rowMeans(pred_mat)

    # SE from member spread
    if (se.fit) {
      se <- apply(pred_mat, 1, stats::sd)
    }

    # Percentile intervals
    if (interval) {
      alpha <- 1 - level
      qs <- t(apply(pred_mat, 1, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE))
      lwr <- qs[, 1]
      upr <- qs[, 2]
    }

    # Choose base fit to return
    preds_out <- if (agg == "parms") preds_parms else preds_mean

    # Mask bad rows
    if (any(bad_rows)) {
      preds_out[bad_rows] <- NA_real_
      if (se.fit) se[bad_rows] <- NA_real_
      if (interval) { lwr[bad_rows] <- NA_real_; upr[bad_rows] <- NA_real_ }
    }

    # Assemble return
    if (se.fit && interval) {
      return(list(fit = preds_out, se.fit = se, lwr = lwr, upr = upr))
    } else if (se.fit) {
      return(list(fit = preds_out, se.fit = se))
    } else if (interval) {
      return(list(fit = preds_out, lwr = lwr, upr = upr))
    } else {
      return(preds_out)
    }
  }

  # parms-only path (no member use)
  if (interval) {
    stop("interval = TRUE requires bootstrap member predictions; agg = 'parms' without coef_matrix cannot produce intervals. ",
         "Re-fit with SVEMnet(..., nBoot >= 2) to populate coef_matrix.")
  }
  if (se.fit) {
    stop("se.fit = TRUE requires bootstrap member predictions; agg = 'parms' without coef_matrix cannot produce SEs. ",
         "Re-fit with SVEMnet(..., nBoot >= 2) to populate coef_matrix.")
  }
  if (any(bad_rows)) preds_parms[bad_rows] <- NA_real_
  preds_parms
}
