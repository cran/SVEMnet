#' Experimental FRW Infinitesimal-Jackknife Variance Diagnostic
#'
#' Estimate the sampling-variance component of the raw SVEM ensemble-mean
#' prediction from covariances between stored fractional-random training
#' weights and member predictions. This adapts the infinitesimal-jackknife
#' calculations for bootstrap-smoothed estimators in Efron (2014) and Wager,
#' Hastie, and Efron (2014) by replacing ordinary-bootstrap counts with the
#' mean-one FRW training weights used by SVEM.
#'
#' @param object A Gaussian \code{svem_model} fitted with
#'   \code{store_member_weights = TRUE}.
#' @param newdata A data frame of predictor values.
#' @param finite_B Logical; default \code{TRUE}. Subtract the ordinary-
#'   bootstrap-inspired finite-ensemble Monte Carlo correction
#'   \eqn{n B^{-2}\sum_b\{T_b(x)-\bar T(x)\}^2}. The untruncated raw estimate
#'   and the correction are always returned.
#'
#' @details
#' For member prediction \eqn{T_b(x)} and stored mean-one training weight
#' \eqn{W_{bi}}, the raw diagnostic is
#' \deqn{V_{IJ}(x)=\sum_i \operatorname{cov}_b\{W_{bi},T_b(x)\}^2.}
#' When \code{finite_B = TRUE}, the finite-ensemble correction is subtracted
#' and negative corrected values are truncated at zero for the reported
#' variance and standard error.
#'
#' This is an experimental diagnostic, not a calibrated confidence interval.
#' The cited IJ results are for ordinary bootstrap counts; the FRW substitution
#' used here is an adaptation. In particular, this function does not correct
#' prediction bias, selection nonregularity, or observation noise, and it must
#' not be used as a prediction interval for a new observation.
#'
#' @return A data frame with one row per row of \code{newdata} and columns
#'   \code{fit}, \code{variance}, \code{se}, \code{variance_raw},
#'   \code{finite_B_correction}, \code{truncated_at_zero}, and
#'   \code{calibrated} (always \code{FALSE}). Attributes record the number of
#'   observations, ensemble members, weighting scheme, method, and caution.
#'
#' @references
#' Efron, B. (2014). Estimation and Accuracy After Model Selection.
#' \emph{Journal of the American Statistical Association}, 109, 991--1007.
#'
#' Wager, S., Hastie, T., and Efron, B. (2014). Confidence Intervals for
#' Random Forests: The Jackknife and the Infinitesimal Jackknife.
#' \emph{Journal of Machine Learning Research}, 15, 1625--1651.
#'
#' @seealso \code{\link{SVEMnet}}, \code{\link{svem_forward}}
#' @export
svem_ij_variance <- function(object, newdata, finite_B = TRUE) {
  finite_B <- .svem_logical_scalar(finite_B, "finite_B")
  if (!inherits(object, "svem_model")) {
    stop("`object` must inherit from 'svem_model'.", call. = FALSE)
  }
  fam <- tolower(if (is.null(object$family)) "gaussian" else object$family)
  if (!identical(fam, "gaussian")) {
    stop("svem_ij_variance() is implemented for Gaussian SVEM fits only.",
         call. = FALSE)
  }
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.", call. = FALSE)
  if (is.null(object$coef_matrix) || !is.matrix(object$coef_matrix)) {
    stop("The fitted object does not contain per-member coefficients.",
         call. = FALSE)
  }
  mw <- object$member_weights
  if (is.null(mw) || is.null(mw$training)) {
    stop("Stored member weights are required; refit with ",
         "store_member_weights = TRUE.", call. = FALSE)
  }
  if (identical(object$weight_scheme, "Identity")) {
    stop("The IJ diagnostic requires random member weights; ",
         "weight_scheme = 'Identity' is not supported.", call. = FALSE)
  }

  W <- as.matrix(mw$training)
  B <- nrow(object$coef_matrix)
  n <- if (!is.null(object$nobs)) as.integer(object$nobs) else ncol(W)
  if (B < 2L) stop("At least two ensemble members are required.", call. = FALSE)
  if (nrow(W) != B || ncol(W) != n || any(!is.finite(W))) {
    stop("Stored training weights are not aligned with coef_matrix/nobs.",
         call. = FALSE)
  }

  newdata <- .svem_validate_newdata_classes(object, newdata)
  if (nrow(newdata) == 0L) {
    return(data.frame(
      fit = numeric(0L), variance = numeric(0L), se = numeric(0L),
      variance_raw = numeric(0L), finite_B_correction = numeric(0L),
      truncated_at_zero = logical(0L), calibrated = logical(0L)
    ))
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  terms_obj <- stats::delete.response(object$terms)
  if (identical(environment(terms_obj), baseenv())) {
    environment(terms_obj) <- asNamespace("stats")
  }
  xlev <- if (!is.null(object$xlevels) && is.list(object$xlevels)) {
    object$xlevels
  } else list()
  if (length(xlev)) {
    for (v in names(xlev)) {
      if (v %in% names(newdata)) {
        newdata[[v]] <- factor(as.character(newdata[[v]]), levels = xlev[[v]])
      }
    }
  }
  mf <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)
  ctr <- object$contrasts
  have_ctr <- !is.null(ctr)
  if (!have_ctr) {
    old_opts <- options("contrasts")
    on.exit(options(old_opts), add = TRUE)
    fit_opts <- object$schema$contrasts_options %||% old_opts$contrasts
    if (!is.null(fit_opts)) options(contrasts = fit_opts)
  } else {
    if (!is.list(ctr)) ctr <- as.list(ctr)
    ctr <- lapply(ctr, function(ci) {
      if (is.character(ci)) get(ci, mode = "function") else ci
    })
  }
  mm <- stats::model.matrix(
    terms_obj, data = mf,
    contrasts.arg = if (have_ctr) ctr else NULL
  )
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) mm <- mm[, -int_col, drop = FALSE]
  bad_rows <- rowSums(!is.finite(mm)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows),
            " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm[!is.finite(mm)] <- 0
  }
  train_cols <- if (!is.null(object$training_X)) {
    colnames(object$training_X)
  } else if (!is.null(object$schema$feature_names)) {
    object$schema$feature_names
  } else colnames(mm)
  mm_use <- matrix(0, nrow(mm), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  }
  storage.mode(mm_use) <- "double"

  coef_matrix <- object$coef_matrix
  intercepts <- coef_matrix[, 1L]
  betas <- coef_matrix[, -1L, drop = FALSE]
  if (!is.null(colnames(betas))) {
    missing_cols <- setdiff(train_cols, colnames(betas))
    if (length(missing_cols)) {
      stop("coef_matrix is missing training-design columns: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    betas <- betas[, train_cols, drop = FALSE]
  } else if (ncol(betas) != length(train_cols)) {
    stop("coef_matrix slope columns do not match the training design.",
         call. = FALSE)
  }
  Tmat <- mm_use %*% t(betas) +
    matrix(intercepts, nrow = nrow(mm_use), ncol = B, byrow = TRUE)

  fit <- rowMeans(Tmat)
  Tc <- Tmat - fit
  Wc <- sweep(W, 2L, colMeans(W), "-")
  cov_i <- Tc %*% Wc / B
  variance_raw <- rowSums(cov_i^2)
  correction <- n / B^2 * rowSums(Tc^2)
  corrected <- if (finite_B) variance_raw - correction else variance_raw
  truncated <- is.finite(corrected) & corrected < 0
  variance <- pmax(corrected, 0)
  out <- data.frame(
    fit = fit,
    variance = variance,
    se = sqrt(variance),
    variance_raw = variance_raw,
    finite_B_correction = if (finite_B) correction else rep(0, length(fit)),
    truncated_at_zero = truncated,
    calibrated = rep(FALSE, length(fit))
  )
  if (any(bad_rows)) out[bad_rows, 1:5] <- NA_real_
  attr(out, "nobs") <- n
  attr(out, "nBoot") <- B
  attr(out, "weight_scheme") <- object$weight_scheme
  attr(out, "method") <- "experimental FRW infinitesimal-jackknife analogue"
  attr(out, "caution") <- paste(
    "Not calibrated: the FRW substitution is an adaptation of ordinary",
    "bootstrap IJ calculations, and prediction bias and observation noise",
    "are not handled."
  )
  out
}

