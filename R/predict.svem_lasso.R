#' Predict Method for SVEM Models
#'
#' Generates predictions from a fitted \code{svem_model}.
#'
#' @param object An object of class \code{svem_model}.
#' @param newdata A data frame of new predictor values.
#' @param debias Logical; default is \code{FALSE}.
#' @param se.fit Logical; if \code{TRUE}, returns standard errors (default is \code{FALSE}).
#' @param ... Additional arguments.
#' @importFrom stats terms reformulate na.pass sd
#' @return Predictions or a list containing predictions and standard errors.
#' @details
#'
#'A debiased fit is available (along with the standard fit). This is provided to allow the user to match the output of JMP.\\ https://www.jmp.com/support/help/en/18.1/?utm_source=help&utm_medium=redirect#page/jmp/overview-of-selfvalidated-ensemble-models.shtml. The debiasing coefficients are always calculated by SVEMnet(), and the predict() function determines whether the raw or debiased predictions are returned via the \code{debias} argument. Default is \code{FALSE} based on performance on unpublished simulation studies.
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @export
predict.svem_model <- function(object, newdata, debias = FALSE, se.fit = FALSE, ...) {
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame.")
  }

  terms_obj <- delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  mf <- model.frame(terms_obj, data = newdata, na.action = na.pass)
  X_new <- model.matrix(terms_obj, data = mf)

  # Remove intercept column
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

  coef_matrix <- object$coef_matrix  # nBoot x (p + 1)
  intercepts <- coef_matrix[, 1]
  betas <- coef_matrix[, -1, drop = FALSE]

  nBoot <- nrow(coef_matrix)
  m <- nrow(X_new)

  predictions_matrix <- X_new %*% t(betas) + matrix(intercepts, nrow = m, ncol = nBoot, byrow = TRUE)
  predictions_mean <- rowMeans(predictions_matrix, na.rm = FALSE)

  if (debias && !is.null(object$debias_fit)) {
    debias_coef <- coef(object$debias_fit)
    if (length(debias_coef) == 2) {
      a <- debias_coef[1]; b <- debias_coef[2]
      predictions_mean <- a + b * predictions_mean
    } else {
      predictions_mean <- rep(debias_coef[1], m)
    }
  }

  if (se.fit) {
    predictions_se <- apply(predictions_matrix, 1, sd, na.rm = TRUE)
    return(list(fit = predictions_mean, se.fit = predictions_se))
  } else {
    return(predictions_mean)
  }
}

