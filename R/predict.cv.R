#' Predict for \code{svem_cv} objects (and convenience wrapper)
#'
#' Generate predictions from a fitted object returned by
#' \code{glmnet_with_cv()}.
#'
#' The design matrix for \code{newdata} is rebuilt using the training
#' \code{terms} (with environment set to \code{baseenv()}), along with the
#' saved factor \code{xlevels} and \code{contrasts} (stored on the object and
#' cached in \code{object$schema}). Columns are aligned robustly to the
#' training order:
#' \itemize{
#'   \item Any training columns that \code{model.matrix()} drops for
#'         \code{newdata} (e.g., a factor collapses to a single level) are
#'         added back as zero columns.
#'   \item Columns are reordered to exactly match the training order.
#'   \item Rows with unseen factor levels are warned about and return
#'         \code{NA} predictions.
#' }
#'
#' For Gaussian fits (\code{family = "gaussian"}), the returned predictions
#' are on the original response (identity-link) scale. For binomial fits
#' (\code{family = "binomial"}), the returned predictions are probabilities
#' in \code{[0,1]} (logit-link inverted via \code{plogis}).
#'
#' If \code{debias = TRUE} and a calibration fit \code{lm(y ~ y_pred)}
#' exists with a finite slope, predictions are linearly transformed as
#' \code{a + b * pred}. Debiasing is only fitted and used for Gaussian
#' models; for binomial models the \code{debias} argument has no effect.
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
#' @param strict Logical; if \code{TRUE}, require an exact column-name match
#'   with the training design (including intercept position) after alignment.
#'   Default \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predictions on the response scale:
#'   numeric fitted values for Gaussian models; probabilities in \code{[0,1]}
#'   for binomial models.
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
  mf_new <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  mm_new <- stats::model.matrix(terms_obj, data = mf_new, contrasts.arg = ctr)

  # Drop intercept (training_X is stored without intercept)
  int_col <- which(colnames(mm_new) == "(Intercept)")
  if (length(int_col)) {
    mm_new <- mm_new[, -int_col, drop = FALSE]
  }

  # Identify rows with any NA columns (typically from unseen levels)
  bad_rows <- rowSums(is.na(mm_new)) > 0L
  if (any(bad_rows)) {
    warning(
      sum(bad_rows),
      " row(s) in newdata contain unseen or missing levels; returning NA for those rows."
    )
    mm_new[is.na(mm_new)] <- 0
  }

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
    numeric(ncol(mm_new) + 1L),
    c("(Intercept)", colnames(mm_new))
  )
  common <- intersect(names(coefs), names(beta))
  beta[common] <- coefs[common]

  # Linear predictor
  lp <- as.numeric(beta[1L] + mm_new %*% beta[-1L])
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
