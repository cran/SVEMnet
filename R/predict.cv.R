#' Predict from glmnet_with_cv Fits (svem_cv Objects)
#'
#' Generate predictions from a fitted object returned by
#' \code{glmnet_with_cv()} (class \code{"svem_cv"}).
#'
#' The design matrix for \code{newdata} is rebuilt using the stored training
#' \code{terms} (with environment set to \code{baseenv()}), together with the
#' saved factor \code{xlevels} and \code{contrasts} (cached in
#' \code{object$schema}). Columns are then aligned back to the training
#' design in a robust way:
#' \itemize{
#'   \item Any training columns that \code{model.matrix()} drops for
#'         \code{newdata} (for example, a factor collapsing to a single level)
#'         are added back as zero columns.
#'   \item Columns are reordered to exactly match the training order.
#'   \item Rows containing unseen factor/character levels are warned about and
#'         their predictions are set to \code{NA}.
#' }
#'
#' For Gaussian fits (\code{family = "gaussian"}), the returned values are on
#' the original response (identity-link) scale. For binomial fits
#' (\code{family = "binomial"}), the returned values are probabilities in
#' \code{[0,1]} (logit-link inverted via \code{plogis()}).
#'
#' If \code{debias = TRUE} and a calibration model \code{lm(y ~ y_pred)} is
#' present with a finite slope, predictions are adjusted via
#' \code{a + b * pred}. Debiasing is only fitted and used for Gaussian models;
#' for binomial models the \code{debias} argument has no effect.
#'
#' \code{predict_cv()} is a small convenience wrapper that simply calls the
#' underlying S3 method \code{predict.svem_cv()}, keeping a single code path
#' for prediction from \code{glmnet_with_cv()} objects.
#'
#' @name predict_cv
#' @aliases predict_cv predict.svem_cv
#'
#' @param object A fitted object returned by \code{glmnet_with_cv()}
#'   (class \code{"svem_cv"}).
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; if \code{TRUE} and a debiasing fit is available,
#'   apply it. Has an effect only for Gaussian models where
#'   \code{debias_fit} is present.
#' @param strict Logical; if \code{TRUE}, enforce a strict column-name match
#'   between the aligned design for \code{newdata} and the training design
#'   (including the intercept position). Default \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predictions on the response scale:
#'   numeric fitted values for Gaussian models; probabilities in \code{[0,1]}
#'   for binomial models. Rows with unseen factor/character levels return
#'   \code{NA}.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[, 1] - 0.5 * X[, 2] + rnorm(n)
#' df_ex <- data.frame(y = as.numeric(y), X)
#' colnames(df_ex) <- c("y", paste0("x", 1:p))
#'
#' fit <- glmnet_with_cv(
#'   y ~ ., df_ex,
#'   glmnet_alpha = 1,
#'   nfolds = 5,
#'   repeats = 2,
#'   seed = 9,
#'   family = "gaussian"
#' )
#'
#' preds_raw <- predict_cv(fit, df_ex)
#' preds_db  <- predict_cv(fit, df_ex, debias = TRUE)
#' cor(preds_raw, df_ex$y)
#'
#' # Binomial example (probability predictions on [0,1] scale)
#' set.seed(2)
#' n2 <- 80; p2 <- 4
#' X2 <- matrix(rnorm(n2 * p2), n2, p2)
#' eta2 <- X2[, 1] - 0.8 * X2[, 2]
#' pr2 <- plogis(eta2)
#' y2 <- rbinom(n2, size = 1, prob = pr2)
#' df_bin <- data.frame(y = y2, X2)
#' colnames(df_bin) <- c("y", paste0("x", 1:p2))
#'
#' fit_bin <- glmnet_with_cv(
#'   y ~ ., df_bin,
#'   glmnet_alpha = c(0.5, 1),
#'   nfolds = 5,
#'   repeats = 2,
#'   seed = 11,
#'   family = "binomial"
#' )
#'
#' prob_hat <- predict_cv(fit_bin, df_bin)
#' summary(prob_hat)
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

  # Determine family name (gaussian / binomial) if available
  fam <- NULL
  if (!is.null(object$family)) {
    fam <- object$family
  } else if (!is.null(object$meta) && !is.null(object$meta$family)) {
    fam <- object$meta$family
  }
  if (inherits(fam, "family")) {
    fam_name <- fam$family
  } else if (is.character(fam) && length(fam) > 0L) {
    fam_name <- fam[1L]
  } else {
    fam_name <- NA_character_
  }
  is_gaussian <- identical(fam_name, "gaussian")
  is_binomial <- identical(fam_name, "binomial")

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
  if (!is.null(ctr)) {
    if (!is.list(ctr)) ctr <- as.list(ctr)
    ctr <- lapply(ctr, function(ci) {
      if (is.character(ci)) get(ci, mode = "function") else ci
    })
  }

  mf_new <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm_new <- stats::model.matrix(terms_obj, data = mf_new, contrasts.arg = ctr)

  # remove intercept: training_X is stored without it
  int_col <- which(colnames(mm_new) == "(Intercept)")
  if (length(int_col)) {
    mm_new <- mm_new[, -int_col, drop = FALSE]
  }

  # unseen levels -> NA columns
  bad_rows <- rowSums(!is.finite(mm_new)) > 0L
  if (any(bad_rows)) {
    warning(
      sum(bad_rows),
      " row(s) in newdata contain unseen or missing levels; returning NA for those rows."
    )
    mm_new[!is.finite(mm_new)] <- 0
  }

  # Mirror training columns exactly; zero-fill missing, reorder
  train_cols <- colnames(object$training_X)
  mm_use <- matrix(0, nrow(mm_new), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm_new), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm_new[, common_cols, drop = FALSE]
  }
  storage.mode(mm_use) <- "double"

  # Strict check if requested
  if (isTRUE(strict)) {
    tr <- c("(Intercept)", colnames(object$training_X))
    te <- c("(Intercept)", colnames(mm_new))
    if (!identical(tr, te)) {
      missing <- setdiff(tr, te)
      extra   <- setdiff(te, tr)
      stop(sprintf(
        "newdata columns do not match training design.\nMissing: %s\nExtra: %s",
        paste(missing, collapse = ", "),
        paste(extra,   collapse = ", ")
      ))
    }
  }

  # Align coefficients to test design by name
  coefs <- object$parms
  beta  <- stats::setNames(
    numeric(ncol(mm_use) + 1L),
    c("(Intercept)", colnames(mm_use))
  )
  common <- intersect(names(coefs), names(beta))
  beta[common] <- coefs[common]

  # Linear predictor
  lp <- as.numeric(beta[1L] + mm_use %*% beta[-1L])
  if (any(bad_rows)) lp[bad_rows] <- NA_real_

  # Map to response scale according to family
  preds <- if (is_binomial) {
    # For binomial, coefficients are on the logit scale; return probabilities
    plogis(lp)
  } else {
    # For gaussian or unknown family, return linear predictor as-is
    lp
  }

  # Optional debias calibration (gaussian only in practice)
  if (debias && is_gaussian && !is.null(object$debias_fit)) {
    ok <- !is.na(preds)
    if (any(ok)) {
      preds[ok] <- as.numeric(
        stats::predict(object$debias_fit, newdata = data.frame(y_pred = preds[ok]))
      )
    }
  }

  preds
}
