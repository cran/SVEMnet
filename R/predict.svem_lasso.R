#' Predict Method for SVEM Models
#'
#' Generates predictions from a fitted \code{svem_model}.
#'
#' @param object An object of class \code{svem_model}.
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; default is \code{FALSE}.
#' @param se.fit Logical; if \code{TRUE}, returns standard errors (default is \code{FALSE}).
#' @param agg Aggregation method for ensemble predictions. One of \code{"parms"} (default)
#'   or \code{"mean"}. \code{"parms"} uses the aggregated coefficients stored in
#'   \code{object$parms} (or \code{parms_debiased} if \code{debias=TRUE}). \code{"mean"}
#'   averages predictions from individual bootstrap members equally.
#' @param ... Additional arguments (currently unused).
#' @importFrom stats terms reformulate na.pass sd delete.response coef
#' @return Predictions or a list containing predictions and standard errors.
#' @details
#' A debiased fit is available (along with the standard fit). This is provided to
#' allow the user to match the output of JMP. The debiasing coefficients are always
#' calculated by \code{SVEMnet()}, and the predict() function determines whether the
#' raw or debiased predictions are returned via the \code{debias} argument (default
#' \code{FALSE}). When \code{se.fit=TRUE} and \code{debias=TRUE}, the reported SE is
#' the bootstrap spread scaled by \eqn{|b|} from the calibration \eqn{y \sim y_{pred}}.
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing
#' the structure of some of the code and the roxygen documentation. The code for the
#' significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @export
predict.svem_model <- function(object, newdata, debias = FALSE, se.fit = FALSE,
                               agg = c("parms", "mean"), ...) {
  agg <- match.arg(agg)

  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.")
  }

  if (is.null(object$training_X) || is.null(object$coef_matrix)) {
    stop("The fitted object lacks required components (training_X/coef_matrix). Was it created by SVEMnet()?")
  }

  terms_obj <- stats::delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  mf <- model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  X_new <- model.matrix(terms_obj, data = mf)

  # Remove intercept column (training_X stored without intercept)
  intercept_col <- which(colnames(X_new) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X_new <- X_new[, -intercept_col, drop = FALSE]
  }

  # Warn about unseen levels (rows with NAs after model.matrix)
  bad_rows <- rowSums(is.na(X_new)) > 0
  if (any(bad_rows)) {
    warning(sum(bad_rows), " row(s) in newdata contain unseen levels; returning NA predictions for those rows.")
  }

  # Column set must match training
  training_colnames <- colnames(object$training_X)
  if (!all(training_colnames == colnames(X_new))) {
    stop("Column names in newdata do not match those in the training data.")
  }

  # Member predictions matrix (used for se.fit and for agg = "mean")
  coef_matrix <- object$coef_matrix  # nBoot x (p + 1)
  nBoot <- nrow(coef_matrix)
  m <- nrow(X_new)

  intercepts <- coef_matrix[, 1]
  betas <- coef_matrix[, -1, drop = FALSE]
  predictions_matrix <- X_new %*% t(betas) +
    matrix(intercepts, nrow = m, ncol = nBoot, byrow = TRUE)

  # Aggregated coefficients path
  coefs_agg <- if (debias && !is.null(object$parms_debiased)) {
    object$parms_debiased
  } else {
    object$parms
  }

  if (!is.numeric(coefs_agg) || length(coefs_agg) != (ncol(X_new) + 1L)) {
    stop("Aggregated coefficients in object are inconsistent with the design matrix.")
  }

  # Choose aggregation
  if (agg == "parms") {
    predictions_mean <- as.vector(X_new %*% coefs_agg[-1L] + coefs_agg[1L])
  } else { # agg == "mean"
    predictions_mean <- rowMeans(predictions_matrix, na.rm = FALSE)
    if (debias && !is.null(object$debias_fit)) {
      # debiasing applies on the aggregated prediction, not per-member
      debias_coef <- stats::coef(object$debias_fit)
      if (length(debias_coef) == 2) {
        a <- debias_coef[1]; b <- debias_coef[2]
        predictions_mean <- a + b * predictions_mean
      } else {
        predictions_mean <- rep(debias_coef[1], m)
      }
    }
  }

  # Impose NA for bad rows (unseen levels)
  if (any(bad_rows)) predictions_mean[bad_rows] <- NA_real_

  if (se.fit) {
    # Bootstrap-based predictive SE from member spread
    predictions_se <- apply(predictions_matrix, 1, stats::sd, na.rm = TRUE)

    # If debiasing is applied, scale SE by |b|
    if (debias && !is.null(object$debias_fit)) {
      debias_coef <- stats::coef(object$debias_fit)
      if (length(debias_coef) >= 2 && is.finite(debias_coef[2])) {
        predictions_se <- abs(debias_coef[2]) * predictions_se
      }
    }

    if (any(bad_rows)) predictions_se[bad_rows] <- NA_real_
    return(list(fit = predictions_mean, se.fit = predictions_se))
  } else {
    return(predictions_mean)
  }
}
