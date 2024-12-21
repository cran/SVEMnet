#' Fit a glmnet Model with Cross-Validation
#'
#' A wrapper function for \code{\link[glmnet]{cv.glmnet}} that takes input arguments in a manner similar to \code{\link{SVEMnet}}.
#' This function searches over multiple \code{alpha} values by running \code{cv.glmnet()} for each provided \code{alpha}, and then
#' selects the combination of \code{alpha} and \code{lambda} with the best cross-validation performance.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param glmnet_alpha Elastic Net mixing parameter(s) (default is \code{c(0, 0.5, 1)}).
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
#' model_cv <- glmnet_with_cv(y ~ X1 + X2, data = data, glmnet_alpha = c(0,0.5,1))
#' predictions <- predict_cv(model_cv, data)
#'
#' @importFrom stats model.frame model.response model.matrix lm predict var
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @export
glmnet_with_cv <- function(formula, data, glmnet_alpha = c(0, 0.5, 1),
                           standardize = TRUE, nfolds = 10, ...) {
  # Create model frame and design matrix
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)

  # Remove intercept column if present
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X <- X[, -intercept_col, drop = FALSE]
  }

  # Precompute a consistent foldid
  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X)))

  # Loop over alphas and perform cv.glmnet
  best_cvm <- Inf
  best_alpha <- NA
  best_model <- NULL

  for (alpha_val in glmnet_alpha) {
    fit_cv <- tryCatch({
      cv.glmnet(
        X, y,
        alpha = alpha_val,
        standardize = standardize,
        foldid = foldid,
        ...
      )
    }, error = function(e) {
      warning(paste("Error in cv.glmnet for alpha", alpha_val, ":", e$message))
      return(NULL)
    })

    if (!is.null(fit_cv)) {
      # Identify minimal cv error for this alpha
      min_cvm <- min(fit_cv$cvm)
      if (min_cvm < best_cvm) {
        best_cvm <- min_cvm
        best_alpha <- alpha_val
        best_model <- fit_cv
      }
    }
  }

  if (is.null(best_model)) {
    stop("All attempts to fit cv.glmnet failed.")
  }

  best_lambda <- best_model$lambda.min
  best_coefs <- as.vector(coef(best_model, s = best_lambda))
  names(best_coefs) <- rownames(coef(best_model, s = best_lambda))

  # Predictions on training data
  y_pred <- drop(predict(best_model, newx = X, s = best_lambda))

  # Debiasing step (if possible)
  debias_fit <- NULL
  y_pred_debiased <- NULL
  if (length(y) >= 10 && var(y_pred) > 0) {
    debias_fit <- lm(y ~ y_pred)
    y_pred_debiased <- predict(debias_fit)
  }

  # Return model object
  result <- list(
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
    terms = attr(mf, "terms")
  )

  return(result)
}




#' Predict Function for glmnet_with_cv Models
#'
#' Generate predictions from the model fitted by \code{\link{glmnet_with_cv}}. This function accepts
#' new data and returns predictions, optionally debiased if a debiasing linear model was fit.
#'
#' @param object An object returned by \code{\link{glmnet_with_cv}}.
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; if \code{TRUE}, applies the debiasing linear model stored in
#' \code{object$debias_fit} (if available). Default is \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @details
#' Predictions are computed by forming the model matrix from \code{newdata} using the stored \code{formula} and \code{terms}
#' in the fitted model object. The coefficients used are those stored in \code{parms}. If \code{debias=TRUE} and a
#' \code{debias_fit} linear model is available, predictions are adjusted by that model.
#'
#' @return A numeric vector of predictions.
#'
#' @seealso \code{\link{glmnet_with_cv}}, \code{\link{SVEMnet}}, \code{\link{predict.svem_model}}
#'
#' @examples
#' set.seed(0)
#' n <- 50
#' X1 <- runif(n)
#' X2 <- runif(n)
#' y <- 1 + 2*X1 + 3*X2 + rnorm(n)
#' data <- data.frame(y, X1, X2)
#' model_cv <- glmnet_with_cv(y ~ X1 + X2, data = data, glmnet_alpha = c(0,0.5,1))
#' predictions <- predict_cv(model_cv, data)
#' predictions_debiased <- predict_cv(model_cv, data, debias = TRUE)
#'
#' @importFrom stats model.frame model.matrix na.pass coef
#' @export
predict_cv <- function(object, newdata, debias = FALSE, ...) {
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.")
  }

  # Use the stored terms object but remove the response variable
  terms_obj <- delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  mf <- model.frame(terms_obj, data = newdata, na.action = na.pass)
  X_new <- model.matrix(terms_obj, data = mf)

  # Remove intercept column if present
  intercept_col <- which(colnames(X_new) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X_new <- X_new[, -intercept_col, drop = FALSE]
  }

  # Check that newdata has same predictors
  training_colnames <- colnames(object$training_X)
  if (!all(training_colnames == colnames(X_new))) {
    stop("Column names in newdata do not match those in the training data.")
  }

  # Extract coefficients
  coefs <- object$parms
  intercept <- coefs[1]
  betas <- coefs[-1]

  # Compute predictions
  predictions <- drop(X_new %*% betas + intercept)

  # Debias if requested and available
  if (debias && !is.null(object$debias_fit)) {
    debias_coef <- coef(object$debias_fit)
    if (length(debias_coef) == 2) {
      a <- debias_coef[1]
      b <- debias_coef[2]
      predictions <- a + b * predictions
    } else {
      predictions <- rep(debias_coef[1], length(predictions))
    }
  }

  return(predictions)
}
