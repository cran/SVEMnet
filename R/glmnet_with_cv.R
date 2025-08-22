#' Fit a glmnet Model with Cross-Validation
#'
#' A wrapper function for \code{\link[glmnet]{cv.glmnet}} that takes input arguments in a manner similar to \code{\link{SVEMnet}}.
#' This function searches over multiple \code{alpha} values by running \code{cv.glmnet()} for each provided \code{alpha}, and then
#' selects the combination of \code{alpha} and \code{lambda} with the best cross-validation performance.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param glmnet_alpha Elastic Net mixing parameter(s) (default is \code{c(1)}).
#' If multiple values are provided, \code{cv.glmnet} is run for each \code{alpha}, and the model with the lowest cross-validation
#' error is selected.
#' @param standardize Logical flag passed to \code{\link[glmnet]{glmnet}}. If \code{TRUE} (default), each variable is standardized
#' before model fitting.
#' @param nfolds Number of cross-validation folds (default is \code{10}).
#' @param ... Additional arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#'
#' @details
#' This function uses \code{\link[glmnet]{cv.glmnet}} to fit a generalized linear model with elastic net regularization,
#' performing k-fold cross-validation to select the regularization parameter \code{lambda}. If multiple \code{alpha} values are
#' provided, it selects the best-performing \code{alpha}-\code{lambda} pair based on the minimal cross-validation error.
#'
#' After fitting, the function calculates a debiasing linear model (if possible). This is done by regressing the actual responses
#' on the fitted values obtained from the selected model. The resulting linear model is stored in \code{debias_fit}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{parms}: Coefficients from the selected \code{cv.glmnet} model at \code{lambda.min}.
#'   \item \code{debias_fit}: A linear model of the form \code{y ~ y_pred} used for debiasing (if applicable).
#'   \item \code{glmnet_alpha}: The vector of \code{alpha} values considered.
#'   \item \code{best_alpha}: The selected \code{alpha} value that gave the best cross-validation result.
#'   \item \code{best_lambda}: The \code{lambda} value chosen by cross-validation at the selected \code{alpha}.
#'   \item \code{actual_y}: The response vector used in the model.
#'   \item \code{training_X}: The predictor matrix used in the model.
#'   \item \code{y_pred}: The fitted values from the final model (no debiasing).
#'   \item \code{y_pred_debiased}: Debiased fitted values if \code{debias_fit} is available.
#'   \item \code{formula}: The formula used for model fitting.
#'   \item \code{terms}: The terms object extracted from the model frame.
#' }
#'
#' @seealso \code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}, \code{\link{SVEMnet}}
#'
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent.
#' \emph{Journal of Statistical Software}, 33(1), 1-22. \doi{10.18637/jss.v033.i01}
#'
#' @examples
#' set.seed(0)
#' n <- 50
#' X1 <- runif(n)
#' X2 <- runif(n)
#' y <- 1 + 2*X1 + 3*X2 + rnorm(n)
#' data <- data.frame(y, X1, X2)
#'
#' model_cv <- glmnet_with_cv(y ~ X1 + X2, data = data, glmnet_alpha = c(1))
#' predictions <- predict_cv(model_cv, data)
#'
#' @importFrom stats model.frame model.response model.matrix lm predict var coef
#' @importFrom glmnet cv.glmnet
#' @export
glmnet_with_cv <- function(formula, data,
                           glmnet_alpha = c(0, 0.5, 1),
                           standardize = TRUE,
                           nfolds = 10,
                           ...){

  # -----------------------------
  # Build model frame & matrix
  # -----------------------------
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  if (nrow(mf) < nrow(as.data.frame(data))) {
    warning("Dropped ", nrow(as.data.frame(data)) - nrow(mf), " row(s) due to missing values.")
  }
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)

  # Remove intercept column for glmnet (it will add its own)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col)) {
    X <- X[, -intercept_col, drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  # -----------------------------
  # Near-zero-variance filter
  # -----------------------------
  drop_msg <- NULL
  # if (p > 0) {
  #   sds <- apply(X, 2, stats::sd)
  #   keep <- is.finite(sds) & (sds > 1e-8)
  #   if (!all(keep)) {
  #     drop_msg <- sprintf("Dropping %d near-zero-variance column(s).", sum(!keep))
  #     X <- X[, keep, drop = FALSE]
  #   }
  # }

  # If everything died, fit intercept-only
  if (ncol(X) == 0) {
    intercept <- mean(y)
    parms <- c("(Intercept)" = intercept)
    y_pred <- rep(intercept, n)
    debias_fit <- NULL
    y_pred_debiased <- NULL
    if (length(y) >= 10 && stats::var(y_pred) > 0) {
      debias_fit <- stats::lm(y ~ y_pred)
      y_pred_debiased <- stats::predict(debias_fit)
    }
    return(list(
      parms = parms,
      debias_fit = debias_fit,
      glmnet_alpha = glmnet_alpha,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      actual_y = y,
      training_X = X,
      y_pred = y_pred,
      y_pred_debiased = y_pred_debiased,
      formula = formula,
      terms = attr(mf, "terms"),
      note = c("intercept_only", drop_msg)
    ))
  }

  # -----------------------------
  # Safe nfolds (ensure ≥3 obs/fold)
  # -----------------------------
  max_folds   <- max(3L, floor(n / 3))
  nfolds_eff  <- min(max(5L, min(nfolds, n)), max_folds)
  foldid      <- sample(rep(seq_len(nfolds_eff), length.out = n))

  # Wide lambda path is helpful in narrow regimes
  cv_args <- list(
    standardize      = standardize,
    foldid           = foldid,
    family           = "gaussian",
    type.measure     = "mse",
    lambda.min.ratio = 1e-5
  )
  dots <- list(...)
  if (length(dots)) cv_args[names(dots)] <- dots

  # -----------------------------
  # Sanitize alpha vector
  # -----------------------------
  glmnet_alpha <- unique(as.numeric(glmnet_alpha))
  glmnet_alpha <- glmnet_alpha[is.finite(glmnet_alpha)]
  glmnet_alpha <- glmnet_alpha[glmnet_alpha >= 0 & glmnet_alpha <= 1]
  if (length(glmnet_alpha) == 0L) glmnet_alpha <- 1

  # Helper to try multiple alphas and pick the best cv model
  .fit_over_alphas <- function(alpha_vec) {
    best_model <- NULL; best_alpha <- NA_real_; best_cvm <- Inf
    for (alpha_val in alpha_vec) {
      fit_cv <- tryCatch({
        do.call(glmnet::cv.glmnet,
                c(list(x = X, y = y, alpha = alpha_val), cv_args))
      }, error = function(e) NULL)

      if (!is.null(fit_cv)) {
        min_cvm <- suppressWarnings(min(fit_cv$cvm))
        if (is.finite(min_cvm) && min_cvm < best_cvm) {
          best_cvm   <- min_cvm
          best_alpha <- alpha_val
          best_model <- fit_cv
        }
      }
    }
    list(model = best_model, alpha = best_alpha)
  }

  # -----------------------------
  # Try cv.glmnet over alphas
  # -----------------------------
  best <- .fit_over_alphas(glmnet_alpha)

  # If all failed and ridge was requested, retry without alpha=0
  if (is.null(best$model) && any(glmnet_alpha == 0)) {
    warning("cv.glmnet failed for all alphas; retrying without alpha=0 (ridge).")
    best <- .fit_over_alphas(setdiff(glmnet_alpha, 0))
  }

  # -----------------------------
  # Ridge fallback (manual CV)
  # -----------------------------
  if (is.null(best$model)) {
    warning("All cv.glmnet attempts failed; switching to ridge fallback with manual CV.")
    lam_seq <- 10 ^ seq(3, -5, length.out = 120)

    fit_ridge <- tryCatch(
      glmnet::glmnet(X, y, alpha = 0, lambda = lam_seq,
                     standardize = standardize, family = "gaussian"),
      error = function(e) NULL
    )
    if (is.null(fit_ridge)) {
      # last-resort intercept-only
      intercept <- mean(y)
      parms <- c("(Intercept)" = intercept)
      return(list(
        parms = parms,
        debias_fit = NULL,
        glmnet_alpha = glmnet_alpha,
        best_alpha = 0,
        best_lambda = NA_real_,
        actual_y = y,
        training_X = X,
        y_pred = rep(intercept, n),
        y_pred_debiased = NULL,
        formula = formula,
        terms = attr(mf, "terms"),
        note = c("ridge_fallback_failed", drop_msg)
      ))
    }

    # manual CV over lam_seq
    mse_by_lambda <- rep(NA_real_, length(lam_seq))
    for (j in seq_along(lam_seq)) {
      lam <- lam_seq[j]
      preds <- rep(NA_real_, n)
      for (fold in seq_len(nfolds_eff)) {
        tr_idx <- which(foldid != fold)
        te_idx <- which(foldid == fold)
        fit_j <- tryCatch(
          glmnet::glmnet(X[tr_idx, , drop = FALSE], y[tr_idx],
                         alpha = 0, lambda = lam, standardize = standardize,
                         family = "gaussian"),
          error = function(e) NULL
        )
        if (is.null(fit_j)) next
        preds[te_idx] <- drop(predict(fit_j, newx = X[te_idx, , drop = FALSE], s = lam))
      }
      mse_by_lambda[j] <- mean((preds - y)^2, na.rm = TRUE)
    }
    j_best <- which.min(mse_by_lambda)
    best_lambda <- lam_seq[j_best]

    # final ridge fit on full data
    final_ridge <- glmnet::glmnet(X, y, alpha = 0, lambda = best_lambda,
                                  standardize = standardize, family = "gaussian")

    # Robust coef extraction (avoid exact=TRUE path)
    coef_mat <- as.matrix(stats::coef(final_ridge, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)

    y_pred <- drop(predict(final_ridge, newx = X, s = best_lambda))

    debias_fit <- NULL
    y_pred_debiased <- NULL
    if (length(y) >= 10 && stats::var(y_pred) > 0) {
      debias_fit <- stats::lm(y ~ y_pred)
      y_pred_debiased <- stats::predict(debias_fit)
    }

    return(list(
      parms = best_coefs,
      debias_fit = debias_fit,
      glmnet_alpha = glmnet_alpha,
      best_alpha = 0,
      best_lambda = best_lambda,
      actual_y = y,
      training_X = X,
      y_pred = y_pred,
      y_pred_debiased = y_pred_debiased,
      formula = formula,
      terms = attr(mf, "terms"),
      note = c("ridge_fallback", drop_msg)
    ))
  }

  # -----------------------------
  # Success path: best cv.glmnet
  # -----------------------------
  best_model  <- best$model
  best_alpha  <- best$alpha
  best_lambda <- best_model$lambda.min

  # Coefficients (robust)
  coef_mat <- tryCatch(
    as.matrix(stats::coef(best_model, s = "lambda.min")),
    error = function(e) {
      # fallback via underlying fit if needed
      as.matrix(stats::coef(best_model$glmnet.fit, s = best_lambda))
    }
  )
  best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)

  # Training predictions for debias check
  y_pred <- drop(predict(best_model, newx = X, s = "lambda.min"))

  debias_fit <- NULL
  y_pred_debiased <- NULL
  if (length(y) >= 10 && stats::var(y_pred) > 0) {
    debias_fit <- stats::lm(y ~ y_pred)
    y_pred_debiased <- stats::predict(debias_fit)
  }

  list(
    parms = best_coefs,
    debias_fit = debias_fit,
    glmnet_alpha = glmnet_alpha,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    actual_y = y,
    training_X = X,
    y_pred = y_pred,
    y_pred_debiased = y_pred_debiased,
    formula = formula,
    terms = attr(mf, "terms"),
    note = drop_msg
  )
}
