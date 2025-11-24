#' Predict Method for SVEM Models (Gaussian and Binomial)
#'
#' Generate predictions from a fitted SVEM model (Gaussian or binomial), with
#' optional bootstrap uncertainty and family-appropriate output scales.
#'
#' This method dispatches on \code{object$family}:
#' \itemize{
#'   \item \strong{Gaussian}: returns predictions on the response (identity)
#'         scale. Optional linear calibration ("debias") learned at training
#'         may be applied.
#'   \item \strong{Binomial}: supports glmnet-style \code{type = "link"},
#'         \code{"response"}, or \code{"class"} predictions. No debiasing is
#'         applied; \code{type = "response"} returns probabilities in
#'         \eqn{[0,1]}.
#' }
#'
#' Uncertainty summaries (\code{se.fit}, \code{interval}) and all binomial
#' predictions are based on per-bootstrap member predictions obtained from
#' the coefficient matrix stored in \code{object$coef_matrix}. If
#' \code{coef_matrix} is \code{NULL}, these options are not available (and
#' binomial prediction will fail). For Gaussian models with
#' \code{se.fit = FALSE} and \code{interval = FALSE}, predictions are computed
#' directly from the aggregated coefficients.
#'
#' @section Design-matrix reconstruction:
#' The function rebuilds the design matrix for \code{newdata} to match the
#' training design:
#' \itemize{
#'   \item Uses the training \code{terms} (with environment set to
#'         \code{baseenv()}).
#'   \item Harmonizes factor and character predictors to the training
#'         \code{xlevels}.
#'   \item Reuses stored per-factor \code{contrasts} when available; otherwise
#'         falls back to saved global contrast options.
#'   \item Zero-fills any columns present at training but absent in
#'         \code{newdata}, and reorders columns to match the training order.
#' }
#' Rows containing unseen factor levels yield \code{NA} predictions (with a
#' warning).
#'
#' @section Aggregation and debiasing:
#' For Gaussian SVEM models:
#' \describe{
#'   \item{Point predictions}{When \code{se.fit = FALSE} and
#'         \code{interval = FALSE}, predictions are computed from the
#'         aggregated coefficients saved at fit time
#'         (\code{object$parms}; or \code{object$parms_debiased} when
#'         \code{debias = TRUE}). This is algebraically equivalent to
#'         averaging member predictions when the coefficients were formed as
#'         bootstrap means.}
#'   \item{Bootstrap-based summaries}{When \code{se.fit = TRUE} and/or
#'         \code{interval = TRUE}, predictions are computed from per-bootstrap
#'         member predictions using \code{object$coef_matrix}. For
#'         \code{debias = TRUE}, the linear calibration is applied to member
#'         predictions before summarizing.}
#' }
#'
#' For binomial SVEM models, predictions are always aggregated from member
#' predictions on the requested scale (probability or link) using
#' \code{coef_matrix}; the stored coefficient averages (\code{parms},
#' \code{parms_debiased}) are retained for diagnostics but are not used in
#' prediction. The \code{debias} argument is ignored and treated as
#' \code{FALSE} for binomial models.
#'
#' For Gaussian fits, if \code{debias = TRUE} and a calibration model
#' \code{lm(y ~ y_pred)} was learned at training, predictions (and, when
#' applicable, member predictions) are transformed by that calibration. This
#' debiasing is never applied for binomial fits.
#'
#' @section Uncertainty:
#' When \code{se.fit = TRUE}, standard errors are computed as the row-wise
#' standard deviations of member predictions on the requested scale. When
#' \code{interval = TRUE}, percentile intervals are computed from member
#' predictions on the requested scale, using the requested \code{level}. Both
#' require a non-null \code{coef_matrix}. For \code{type = "class"} (binomial),
#' uncertainty summaries are not available.
#'
#' @param object A fitted SVEM model (class \code{svem_model}; binomial models
#'   typically also inherit class \code{svem_binomial}). Created by
#'   \code{SVEMnet()}.
#' @param newdata A data frame of new predictor values.
#' @param type (Binomial only) One of:
#'   \itemize{
#'     \item \code{"response"} (default): predicted probabilities in
#'           \eqn{[0,1]}.
#'     \item \code{"link"}: linear predictor (log-odds).
#'     \item \code{"class"}: 0/1 class labels (threshold 0.5). Uncertainty
#'           summaries are not available for this type.
#'   }
#'   Ignored for Gaussian models.
#' @param debias (Gaussian only) Logical; default \code{FALSE}. If \code{TRUE},
#'   apply the linear calibration fit \code{lm(y ~ y_pred)} learned at
#'   training when available. Ignored (and internally set to \code{FALSE})
#'   for binomial models.
#' @param se.fit Logical; if \code{TRUE}, return bootstrap standard errors
#'   computed from member predictions (requires \code{coef_matrix}). Not
#'   available for \code{type = "class"}. For Gaussian models, this forces use
#'   of bootstrap member predictions instead of aggregate coefficients.
#' @param interval Logical; if \code{TRUE}, return percentile confidence
#'   limits from member predictions (requires \code{coef_matrix}). Not
#'   available for \code{type = "class"}. For Gaussian models, this forces use
#'   of bootstrap member predictions instead of aggregate coefficients.
#' @param level Confidence level for percentile intervals. Default
#'   \code{0.95}.
#' @param ... Currently unused.
#'
#' @return
#' If \code{se.fit = FALSE} and \code{interval = FALSE}:
#' \itemize{
#'   \item \strong{Gaussian}: a numeric vector of predictions on the response
#'         (identity) scale.
#'   \item \strong{Binomial}: a numeric vector for \code{type = "response"}
#'         (probabilities) or \code{type = "link"} (log-odds), or an integer
#'         vector of 0/1 labels for \code{type = "class"}.
#' }
#'
#' If \code{se.fit} and/or \code{interval} are \code{TRUE} (and
#' \code{type != "class"}), a list with components:
#' \itemize{
#'   \item \code{fit}: predictions on the requested scale.
#'   \item \code{se.fit}: bootstrap standard errors (when
#'         \code{se.fit = TRUE}).
#'   \item \code{lwr}, \code{upr}: percentile confidence limits (when
#'         \code{interval = TRUE}).
#' }
#'
#' Rows containing unseen or missing factor levels produce \code{NA}
#' predictions (and \code{NA} SEs/intervals), with a warning.
#'
#' @seealso \code{\link{SVEMnet}}
#'
#' @examples
#' ## ---- Gaussian example -------------------------------------------------
#' set.seed(1)
#' n  <- 60
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' y  <- 1 + 0.8 * X1 - 0.6 * X2 + 0.2 * X3 + rnorm(n, 0, 0.4)
#' dat <- data.frame(y, X1, X2, X3)
#'
#' fit_g <- SVEMnet(
#'   y ~ (X1 + X2 + X3)^2, dat,
#'   nBoot = 40, glmnet_alpha = c(1, 0.5),
#'   relaxed = TRUE, family = "gaussian"
#' )
#'
#' ## Aggregate-coefficient predictions (with and without debiasing)
#' p_g  <- predict(fit_g, dat)                # debias = FALSE (default)
#' p_gd <- predict(fit_g, dat, debias = TRUE) # apply calibration, if available
#'
#' ## Bootstrap-based uncertainty (requires coef_matrix)
#' out_g <- predict(
#'   fit_g, dat,
#'   debias   = TRUE,
#'   se.fit   = TRUE,
#'   interval = TRUE,
#'   level    = 0.90
#' )
#' str(out_g)
#'
#' \donttest{
#' ## ---- Binomial example ------------------------------------------------
#' set.seed(2)
#' n  <- 120
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' eta <- -0.3 + 1.1 * X1 - 0.8 * X2 + 0.5 * X1 * X3
#' p   <- plogis(eta)
#' yb  <- rbinom(n, 1, p)
#' db  <- data.frame(yb = yb, X1 = X1, X2 = X2, X3 = X3)
#'
#' fit_b <- SVEMnet(
#'   yb ~ (X1 + X2 + X3)^2, db,
#'   nBoot = 50, glmnet_alpha = c(1, 0.5),
#'   relaxed = TRUE, family = "binomial"
#' )
#'
#' ## Probabilities, link, and classes
#' p_resp <- predict(fit_b, db, type = "response")
#' p_link <- predict(fit_b, db, type = "link")
#' y_hat  <- predict(fit_b, db, type = "class")  # 0/1 labels (no SE or interval)
#'
#' ## Bootstrap-based uncertainty on the probability scale
#' out_b <- predict(
#'   fit_b, db,
#'   type     = "response",
#'   se.fit   = TRUE,
#'   interval = TRUE,
#'   level    = 0.90
#' )
#' str(out_b)
#' }
#'
#' @family SVEM methods
#' @importFrom stats plogis qlogis
#' @export
#' @method predict svem_model
predict.svem_model <- function(object, newdata,
                               type = c("response", "link", "class"),
                               debias = FALSE,
                               se.fit = FALSE,
                               interval = FALSE,
                               level = 0.95,
                               ...) {
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.")
  if (is.null(object$terms) || is.null(object$parms)) {
    stop("The fitted object is missing required components (terms/parms).")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0,1).")
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  fam <- tolower(object$family %||% "gaussian")

  # --- Family-specific handling ----------------------------------------------
  if (identical(fam, "binomial")) {
    type <- match.arg(type)
    if (type == "class" && (se.fit || interval)) {
      stop("type = 'class' does not support se.fit or interval. ",
           "Use type = 'response' or 'link' for uncertainty estimates.")
    }
    debias <- FALSE  # never debias binomial
  } else {
    # Gaussian: ignore type
    type <- NULL
  }

  # --- Build model frame/matrix exactly like at fit time ---------------------
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

  mf  <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)

  # Prefer exact per-factor contrasts saved at fit time
  ctr <- object$contrasts
  have_ctr <- !is.null(ctr)
  if (!have_ctr) {
    old_opts <- options("contrasts"); on.exit(options(old_opts), add = TRUE)
    fit_opts <- object$schema$contrasts_options %||% old_opts$contrasts
    if (!is.null(fit_opts)) options(contrasts = fit_opts)
  } else {
    if (!is.list(ctr)) ctr <- as.list(ctr)
    # Coerce any character contrast names back to functions
    ctr <- lapply(ctr, function(ci) {
      if (is.character(ci)) get(ci, mode = "function") else ci
    })
  }

  mm <- stats::model.matrix(
    terms_obj,
    data = mf,
    contrasts.arg = if (have_ctr) ctr else NULL
  )

  # Remove intercept col from mm; glmnet already had an intercept
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) mm <- mm[, -int_col, drop = FALSE]

  # Identify rows with non-finite (e.g., unseen levels -> NA)
  bad_rows <- rowSums(!is.finite(mm)) > 0L
  if (any(bad_rows)) {
    warning(sum(bad_rows),
            " row(s) in newdata contain unseen or missing levels; returning NA for those rows.")
    mm[!is.finite(mm)] <- 0
  }

  # Mirror training columns exactly; zero-fill missing
  train_cols <- if (!is.null(object$training_X)) {
    colnames(object$training_X)
  } else if (!is.null(object$schema$feature_names)) {
    object$schema$feature_names
  } else colnames(mm)

  extra_cols <- setdiff(colnames(mm), train_cols)
  if (length(extra_cols)) {
    nz <- extra_cols[colSums(abs(mm[, extra_cols, drop = FALSE])) > 0]
    if (length(nz)) {
      message("Ignoring ", length(nz), " column(s) not present at training: ",
              paste(utils::head(nz, 5L), collapse = ", "),
              if (length(nz) > 5L) " ..." else "")
    }
  }

  mm_use <- matrix(0, nrow(mm), length(train_cols))
  colnames(mm_use) <- train_cols
  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  storage.mode(mm_use) <- "double"

  # --- Aggregate-coefficient prediction via parms/parms_debiased ------------
  coefs <- if (!identical(fam, "binomial") && isTRUE(debias) && !is.null(object$parms_debiased))
    object$parms_debiased else object$parms

  beta <- stats::setNames(numeric(length(train_cols) + 1L),
                          c("(Intercept)", train_cols))
  common_beta <- intersect(names(coefs), names(beta))
  beta[common_beta] <- coefs[common_beta]
  eta_parms <- as.numeric(beta[1L] + mm_use %*% beta[-1L])

  # Transform helper (binomial only)
  transform_eta <- function(eta) {
    if (identical(fam, "binomial")) {
      if (type == "link") eta else plogis(eta)
    } else {
      # gaussian always returns on response (identity) scale
      eta
    }
  }

  # --- Do we need bootstrap members? ----------------------------------------
  have_members <- !is.null(object$coef_matrix) && is.matrix(object$coef_matrix)

  if (identical(fam, "binomial")) {
    # Binomial always uses member predictions (like old agg = "mean")
    if (!have_members) {
      stop("Binomial SVEM predictions require bootstrap member coefficients (coef_matrix). ",
           "Re-fit with nBoot >= 2 to populate coef_matrix.")
    }
  }

  need_members <- if (identical(fam, "binomial")) TRUE else (se.fit || interval)

  # --- Gaussian, no SE / interval: simple aggregate-only path ----------------
  if (identical(fam, "gaussian") && !need_members) {
    preds <- eta_parms
    if (any(bad_rows)) preds[bad_rows] <- NA_real_
    return(preds)
  }

  # --- Member-based predictions (binomial, or gaussian with SE/interval) -----
  if (!have_members) {
    stop("se.fit/interval require bootstrap member predictions (coef_matrix). ",
         "Re-fit with nBoot >= 2 to populate coef_matrix.")
  }

  coef_matrix <- object$coef_matrix  # nBoot x (p+1)

  if (!is.matrix(coef_matrix) || ncol(coef_matrix) < 1L) {
    stop("`coef_matrix` must be a non-empty matrix with an intercept column.")
  }

  m  <- nrow(mm_use)
  nb <- nrow(coef_matrix)

  intercepts <- coef_matrix[, 1]
  betas      <- coef_matrix[, -1, drop = FALSE]

  # Align slope coefficients to training design columns
  if (!is.null(colnames(betas))) {
    missing_cols <- setdiff(train_cols, colnames(betas))
    if (length(missing_cols)) {
      stop(
        "`coef_matrix` is missing slope coefficients for training design columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    betas <- betas[, train_cols, drop = FALSE]
  } else if (ncol(betas) != length(train_cols)) {
    stop(
      "Number of slope columns in `coef_matrix` (excluding intercept) ",
      "does not match number of training design columns."
    )
  }

  # Member predictions on eta scale
  eta_mat <- mm_use %*% t(betas) + matrix(intercepts, nrow = m, ncol = nb, byrow = TRUE)


  # Gaussian debias applies to members, too
  if (identical(fam, "gaussian") && isTRUE(debias) && !is.null(object$debias_fit)) {
    ab <- stats::coef(object$debias_fit)
    if (length(ab) >= 2 && all(is.finite(ab[1:2]))) {
      eta_mat <- ab[1] + ab[2] * eta_mat
    }
  }

  # Binomial class: aggregate on prob scale, then threshold (no se/interval)
  if (identical(fam, "binomial") && type == "class") {
    p_mean <- rowMeans(plogis(eta_mat))
    preds  <- ifelse(p_mean >= 0.5, 1L, 0L)
    if (any(bad_rows)) preds[bad_rows] <- NA_integer_
    return(preds)
  }

  # Transform to requested scale
  pred_mat <- transform_eta(eta_mat)

  # Mean prediction on requested scale
  preds_mean <- rowMeans(pred_mat)

  # SEs / percentile CIs on requested scale
  if (se.fit) {
    se <- apply(pred_mat, 1, stats::sd)
  }
  if (interval) {
    alpha <- 1 - level
    qs <- t(apply(pred_mat, 1, stats::quantile,
                  probs = c(alpha/2, 1 - alpha/2),
                  na.rm = TRUE, names = FALSE))
    lwr <- qs[, 1]; upr <- qs[, 2]
  }

  # Mask NA rows
  if (any(bad_rows)) {
    preds_mean[bad_rows] <- NA_real_
    if (se.fit)   se[bad_rows]   <- NA_real_
    if (interval) { lwr[bad_rows] <- NA_real_; upr[bad_rows] <- NA_real_ }
  }

  # Return based on requested outputs
  if (se.fit && interval) return(list(fit = preds_mean, se.fit = se, lwr = lwr, upr = upr))
  if (se.fit)             return(list(fit = preds_mean, se.fit = se))
  if (interval)           return(list(fit = preds_mean, lwr = lwr, upr = upr))
  return(preds_mean)
}
