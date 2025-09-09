#' Predict for svem_cv objects (and convenience wrapper)
#'
#' Generates predictions from a fitted object returned by \code{glmnet_with_cv()}.
#'
#' The design matrix for \code{newdata} is rebuilt using the training \code{terms}
#' (with environment set to \code{baseenv()}), along with the saved factor
#' \code{xlevels} and \code{contrasts} (stored in \code{object$schema}). Columns are
#' aligned robustly to the training order:
#' \itemize{
#'   \item Any training columns that \code{model.matrix()} drops for \code{newdata}
#'         (e.g., a factor collapses to a single level) are added back as zero columns.
#'   \item Columns are reordered to exactly match the training order.
#'   \item Rows with unseen factor levels are warned about and return \code{NA}.
#' }
#'
#' If \code{debias = TRUE} and a calibration fit \code{lm(y ~ y_pred)} exists with a
#' finite slope, predictions are transformed by \code{a + b * pred}.
#'
#' @name predict_cv
#' @aliases predict_cv predict.svem_cv
#'
#' @param object A list returned by \code{glmnet_with_cv()} (class \code{svem_cv}).
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; if TRUE and a debiasing fit is available, apply it.
#' @param strict Logical; if TRUE, require exact column-name match with training design
#'   (including intercept position) after alignment. Default \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predictions.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[,1] - 0.5 * X[,2] + rnorm(n)
#' df_ex <- data.frame(y = as.numeric(y), X)
#' colnames(df_ex) <- c("y", paste0("x", 1:p))
#' fit <- glmnet_with_cv(y ~ ., df_ex, glmnet_alpha = 1, nfolds = 5, repeats = 2, seed = 9)
#' preds_raw <- predict_cv(fit, df_ex)
#' preds_db  <- predict_cv(fit, df_ex, debias = TRUE)
#' cor(preds_raw, df_ex$y)
#'
#' @importFrom stats delete.response model.frame model.matrix na.pass predict setNames coef
#' @export
predict_cv <- function(object, newdata, debias = FALSE, strict = FALSE, ...) {
  # Delegate to the S3 method to keep a single code path.
  predict(object, newdata = newdata, debias = debias, strict = strict, ...)
}

#' @rdname predict_cv
#' @export
#' @method predict svem_cv
predict.svem_cv <- function(object, newdata, debias = FALSE, strict = FALSE, ...) {
  stopifnot(is.list(object), is.data.frame(newdata))
  if (is.null(object$terms) || is.null(object$parms) || is.null(object$training_X)) {
    stop("The object is missing required components (terms/parms/training_X).")
  }

  # Training terms (environment nuked for safety)
  terms_obj <- stats::delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  # Harmonize factor/character levels to training xlevels
  xlev <- if (!is.null(object$xlevels) && is.list(object$xlevels)) object$xlevels else NULL
  if (!is.null(xlev) && length(xlev)) {
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

  # Build model frame/matrix with saved contrasts; allow NAs (unseen levels become NA)
  ctr <- if (!is.null(object$contrasts)) object$contrasts else NULL
  mf_new <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm_new <- stats::model.matrix(terms_obj, data = mf_new, contrasts.arg = ctr)

  # Drop intercept (training_X is stored without intercept)
  int_col <- which(colnames(mm_new) == "(Intercept)")
  if (length(int_col)) mm_new <- mm_new[, -int_col, drop = FALSE]

  # Identify rows with any NA columns (typically from unseen levels)
  bad_rows <- rowSums(is.na(mm_new)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows),
            " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm_new[is.na(mm_new)] <- 0
  }

  # Strict check if requested
  if (isTRUE(strict)) {
    tr <- c("(Intercept)", colnames(object$training_X))
    te <- c("(Intercept)", colnames(mm_new))
    if (!identical(tr, te)) {
      stop("Column names in newdata do not match training design (strict=TRUE).")
    }
  }

  # Align coefficients to test design by name
  coefs <- object$parms                       # includes "(Intercept)"
  beta  <- stats::setNames(numeric(ncol(mm_new) + 1L),
                           c("(Intercept)", colnames(mm_new)))
  common <- intersect(names(coefs), names(beta))
  beta[common] <- coefs[common]

  # Predict
  preds <- as.numeric(beta[1L] + mm_new %*% beta[-1L])
  if (any(bad_rows)) preds[bad_rows] <- NA_real_

  # Optional debias calibration
  if (debias && !is.null(object$debias_fit)) {
    ok <- !is.na(preds)
    if (any(ok)) {
      preds[ok] <- as.numeric(stats::predict(object$debias_fit,
                                             newdata = data.frame(y_pred = preds[ok])))
    }
  }

  preds
}
