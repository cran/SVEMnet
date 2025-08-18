#' Predict Method for glmnet_with_cv Objects
#'
#' Generates predictions from a fitted object returned by \code{glmnet_with_cv()}.
#'
#' @param object A list returned by \code{glmnet_with_cv()}.
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; if \code{TRUE} and a debiasing fit is available, apply it (default \code{FALSE}).
#' @param strict Logical; if \code{TRUE}, require exact column-name match with training design (default \code{FALSE}).
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predictions.
#'
#' @details Columns are aligned by name. With \code{strict=TRUE}, a mismatch errors.
#'
#' @importFrom stats delete.response model.frame model.matrix na.pass predict setNames
#' @export

# Robust prediction for glmnet_with_cv() objects
predict_cv <- function(object, newdata, debias = FALSE, strict = FALSE, ...) {
  stopifnot(is.list(object), is.data.frame(newdata), !is.null(object$terms), !is.null(object$parms))

  # Build test model matrix using the training terms (no response)
  terms_obj <- delete.response(object$terms)
  # don't nuke its environment; it contains contrasts info

  mf_new <- model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm_new <- model.matrix(terms_obj, data = mf_new)

  # Mark rows with unseen levels / NA columns produced by model.matrix
  bad_rows <- rowSums(is.na(mm_new)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows), " row(s) in newdata contain unseen levels; returning NA for those rows.")
    # replace NAs with 0 so matrix multiply works; we'll drop them after
    mm_new[is.na(mm_new)] <- 0
  }

  # Coefficients (should include "(Intercept)")
  coefs <- object$parms

  # Optional strict check (train/test columns must match)
  if (isTRUE(strict)) {
    tr_cols <- colnames(object$training_X)
    # training_X was built without the intercept; add it for a fair check
    tr_cols <- c("(Intercept)", tr_cols)
    if (!identical(tr_cols, colnames(mm_new))) {
      stop("Column names in newdata do not match those used in training (strict=TRUE).")
    }
  }

  # Align by names (robust): build a full beta vector over test columns
  beta <- setNames(numeric(ncol(mm_new)), colnames(mm_new))
  common <- intersect(names(coefs), names(beta))
  beta[common] <- coefs[common]

  preds <- as.numeric(mm_new %*% beta)
  if (any(bad_rows)) preds[bad_rows] <- NA_real_

  if (debias && !is.null(object$debias_fit)) {
    # Keep NA rows as NA through debiasing
    ok <- !is.na(preds)
    preds[ok] <- as.numeric(predict(object$debias_fit, newdata = data.frame(y_pred = preds[ok])))
  }

  preds
}
