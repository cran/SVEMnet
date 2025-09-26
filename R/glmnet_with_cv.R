#' Fit a glmnet Model with Cross-Validation
#'
#' Repeated K-fold CV over a per-alpha lambda path, with a proper
#' 1-SE rule across repeats. Preserves fields expected by predict_cv().
#' Optionally uses glmnet's built-in relaxed elastic net for both the warm-start
#' path and each CV fit. When relaxed=TRUE, the final coefficients are taken
#' from a cv.glmnet() object at the chosen lambda so that the returned model
#' reflects the relaxed solution (including its gamma).
#'
#' To avoid duplicate-argument errors, arguments like x, y, alpha, lambda,
#' foldid, and exclude are passed explicitly and removed from the dots
#' before calling glmnet::cv.glmnet().
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param glmnet_alpha Numeric vector of alphas, default c(0,0.25,0.5,0.75, 1).
#' @param standardize Logical passed to glmnet (default TRUE).
#' @param nfolds CV folds (default 10), internally constrained so >= ~3 obs/fold.
#' @param repeats Number of independent CV repeats (default 5).
#' @param choose_rule "1se" or "min" (default). In simulations "1se" lead to increased RMSE on holdout data when simulating small mixture designs.
#' @param seed Optional integer seed for reproducible fold IDs.
#' @param exclude Optional vector OR function for glmnet's exclude=.
#'   If a function, cv.glmnet applies it inside each training fold (glmnet >= 4.1-2).
#' @param relaxed Logical; if TRUE, call glmnet/cv.glmnet with relax=TRUE (default FALSE).
#' @param relax_gamma Optional numeric vector passed as gamma= to glmnet/cv.glmnet
#'   when relaxed=TRUE. If NULL, glmnet uses its internal default gamma path.
#' @param ... Args forwarded to both cv.glmnet() and glmnet(), e.g.
#'   family, weights, parallel, type.measure, intercept, maxit, lower.limits,
#'   upper.limits, penalty.factor, offset, standardize.response, keep, etc.
#'
#' @return A list with elements:
#' \itemize{
#'   \item parms Named numeric vector of coefficients (including "(Intercept)").
#'   \item glmnet_alpha Numeric vector of alphas searched.
#'   \item best_alpha Numeric; winning alpha.
#'   \item best_lambda Numeric; winning lambda.
#'   \item y_pred In-sample predictions from the returned coefficients.
#'   \item debias_fit Optional lm(y ~ y_pred) for Gaussian family.
#'   \item y_pred_debiased If debias_fit exists, its fitted values.
#'   \item cv_summary Per-alpha data frames with lambda, mean_cvm, sd_cvm, se_combined, n_repeats, idx_min, idx_1se.
#'   \item formula Original formula.
#'   \item terms Cleaned training terms (environment set to baseenv()).
#'   \item training_X Training design matrix without intercept.
#'   \item actual_y Training response vector.
#'   \item xlevels Factor levels seen during training (for safe predict).
#'   \item contrasts Contrasts used during training (for safe predict).
#'   \item schema list(feature_names, terms_str, xlevels, contrasts, terms_hash) for deterministic predict.
#'   \item note Character vector of notes (e.g., dropped rows, relaxed-coef source).
#'   \item meta List: nfolds, repeats, rule, family, relaxed, relax_cv_fallbacks, cv_object (if keep=TRUE for the final fit).
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1, -1, rep(0, p - 2))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' df_ex <- data.frame(y = y, X)
#' colnames(df_ex) <- c("y", paste0("x", 1:p))
#'
#' # Default: 1SE rule, repeats=1, non-relaxed
#' fit_1se <- glmnet_with_cv(y ~ ., df_ex, glmnet_alpha = c(0.5, 1),
#'                           nfolds = 5, repeats = 1, seed = 42)
#' str(fit_1se$parms)
#'
#' # v1-like behavior: choose_rule="min"
#' fit_min <- glmnet_with_cv(y ~ ., df_ex, glmnet_alpha = 1,
#'                           nfolds = 5, repeats = 1, choose_rule = "min", seed = 42)
#'
#' # Relaxed path with gamma search
#' fit_relax <- glmnet_with_cv(y ~ ., df_ex, glmnet_alpha = 1,
#'                             nfolds = 5, repeats = 1, relaxed = TRUE, seed = 42)
#'
#' @importFrom stats model.frame model.response model.matrix var coef predict lm delete.response
#' @importFrom glmnet glmnet cv.glmnet
#' @export
glmnet_with_cv <- function(formula, data,
                           glmnet_alpha = c(0,0.25,0.5,0.75, 1),
                           standardize = TRUE,
                           nfolds = 10,
                           repeats = 5,
                           choose_rule = c("min","1se"),
                           seed = NULL,
                           exclude = NULL,
                           relaxed = FALSE,
                           relax_gamma = NULL,
                           ...) {

  choose_rule <- match.arg(choose_rule)

  # --- Build model frame & matrix ---
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  drop_n <- nrow(as.data.frame(data)) - nrow(mf)
  drop_msg <- if (drop_n > 0) sprintf("Dropped %d row(s) due to missing values.", drop_n) else NULL

  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col)) X <- X[, -intercept_col, drop = FALSE]
  n <- nrow(X); p <- ncol(X)

  # --- Collect dots & align weights/offset to rows kept by na.omit ---
  dots <- list(...)
  omit <- attr(mf, "na.action")

  w_all_input <- if ("weights" %in% names(dots)) as.numeric(dots$weights) else NULL
  w <- NULL
  if (!is.null(w_all_input)) {
    w <- if (is.null(omit)) w_all_input else w_all_input[-omit]
    if (length(w) != n) w <- NULL
  }

  off_all_input <- if ("offset" %in% names(dots)) as.numeric(dots$offset) else NULL
  off <- NULL
  if (!is.null(off_all_input)) {
    off <- if (is.null(omit)) off_all_input else off_all_input[-omit]
    if (length(off) != n) off <- NULL
  }

  # --- Intercept-only path (early return) ---
  if (p == 0L) {
    intercept <- if (!is.null(w)) sum(w * y) / sum(w) else mean(y)
    # Clean terms and schema for consistency with normal path
    terms_clean <- attr(mf, "terms"); environment(terms_clean) <- baseenv()
    feature_names <- character(0)

    # Record factor AND character levels for completeness (even if unused)
    xlevels <- list()
    vars_in_terms <- base::all.vars(stats::delete.response(terms_clean))
    for (v in vars_in_terms) {
      if (v %in% colnames(mf)) {
        if (is.factor(mf[[v]])) {
          xlevels[[v]] <- levels(mf[[v]])
        } else if (is.character(mf[[v]])) {
          xlevels[[v]] <- sort(unique(as.character(mf[[v]])))
        }
      }
    }
    contrasts_used <- NULL

    terms_str <- tryCatch(paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
                          error = function(e) NA_character_)
    safe_hash <- function(s) {
      if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
      bytes <- charToRaw(paste0(s, collapse = ""))
      sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
    }
    schema <- list(
      feature_names = feature_names,
      terms_str     = terms_str,
      xlevels       = xlevels,
      contrasts     = contrasts_used,
      terms_hash    = safe_hash(terms_str)
    )
    return(list(
      parms = c("(Intercept)" = intercept),
      glmnet_alpha = glmnet_alpha,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      y_pred = rep(intercept, n),
      debias_fit = NULL,
      y_pred_debiased = NULL,
      cv_summary = list(),
      formula = formula,
      terms = terms_clean,
      training_X = X,
      actual_y = y,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = c("intercept_only", drop_msg),
      meta = list(nfolds = NA_integer_, repeats = NA_integer_,
                  rule = choose_rule, family = "gaussian", relaxed = relaxed,
                  relax_cv_fallbacks = 0L, cv_object = NULL)
    ))
  }

  # --- Safe folds & repeats ---
  max_folds  <- max(3L, floor(n / 3))
  nfolds_eff <- min(max(5L, min(nfolds, n)), max_folds)
  repeats    <- max(1L, as.integer(repeats))

  # --- Clean alphas ---
  glmnet_alpha <- unique(as.numeric(glmnet_alpha))
  glmnet_alpha <- glmnet_alpha[is.finite(glmnet_alpha)]
  glmnet_alpha <- glmnet_alpha[glmnet_alpha >= 0 & glmnet_alpha <= 1]

  if (isTRUE(relaxed) && any(glmnet_alpha == 0)) {
    warning("Dropping alpha = 0 (ridge) for relaxed fits; ridge + relaxed is not supported here.")
    glmnet_alpha <- glmnet_alpha[glmnet_alpha != 0]
  }
  if (!length(glmnet_alpha)) glmnet_alpha <- 1

  # --- Defaults; allow dots to override ---
  base_cv_args <- list(
    standardize  = standardize,
    family       = "gaussian",
    type.measure = "mse",
    maxit        = 1e6,
    relax        = isTRUE(relaxed)
  )
  if (isTRUE(relaxed) && !is.null(relax_gamma)) base_cv_args$gamma <- relax_gamma
  cv_base_args <- utils::modifyList(base_cv_args, dots, keep.null = TRUE)
  fam <- cv_base_args$family

  glmnet_formals <- names(formals(glmnet::glmnet))
  base_glmnet_args <- list(
    standardize = standardize,
    family      = fam,
    intercept   = TRUE,
    exclude     = exclude,
    maxit       = 1e6,
    relax       = isTRUE(relaxed)
  )
  if (isTRUE(relaxed) && !is.null(relax_gamma)) base_glmnet_args$gamma <- relax_gamma
  glmnet_fit_args_full <- utils::modifyList(base_glmnet_args, dots, keep.null = TRUE)
  glmnet_fit_args <- glmnet_fit_args_full[intersect(names(glmnet_fit_args_full), glmnet_formals)]
  glmnet_fit_args <- glmnet_fit_args[setdiff(names(glmnet_fit_args), c("x","y","alpha","lambda"))]

  if (!is.null(w)) {
    cv_base_args$weights <- w
    glmnet_fit_args$weights <- w
  } else {
    cv_base_args$weights <- NULL
    glmnet_fit_args$weights <- NULL
  }
  if (!is.null(off)) {
    cv_base_args$offset <- off
    glmnet_fit_args$offset <- off
  } else {
    cv_base_args$offset <- NULL
    glmnet_fit_args$offset <- NULL
  }

  is_gaussian <- (is.character(fam) && fam == "gaussian") ||
    (inherits(fam, "family") && isTRUE(fam$family == "gaussian"))

  # Helper to drop reserved names from arg lists before do.call
  drop_reserved <- function(lst, reserved) {
    if (!length(lst)) return(lst)
    lst[setdiff(names(lst), reserved)]
  }

  # Build one set of fold IDs reused across alphas and repeats (paired comparison)
  if (!is.null(seed)) set.seed(seed)
  foldids <- replicate(repeats, sample(rep(seq_len(nfolds_eff), length.out = n)),
                       simplify = FALSE)

  # --- One alpha -> repeated CV ---
  fit_alpha_repeated <- function(alpha_val) {
    reserved_cv <- c("x","y","alpha","lambda","foldid","exclude")

    cvms_list  <- list()
    cvSEs_list <- list()
    lam_ref    <- NULL
    fallback_count <- 0L

    for (r in seq_len(repeats)) {
      foldid <- foldids[[r]]

      fit_cv <- tryCatch(
        do.call(glmnet::cv.glmnet,
                c(list(x = X, y = y,
                       alpha  = alpha_val,
                       foldid = foldid,
                       exclude = exclude),
                  drop_reserved(cv_base_args, reserved_cv))),
        error = function(e) {
          if (isTRUE(relaxed)) {
            warning("cv.glmnet(relax=TRUE) failed; retrying without relax for this repeat/alpha.")
            args2 <- cv_base_args; args2$relax <- NULL; args2$gamma <- NULL
            fallback_count <<- fallback_count + 1L
            tryCatch(
              do.call(glmnet::cv.glmnet,
                      c(list(x = X, y = y,
                             alpha  = alpha_val,
                             foldid = foldid,
                             exclude = exclude),
                        drop_reserved(args2, reserved_cv))),
              error = function(e2) NULL
            )
          } else NULL
        }
      )
      if (is.null(fit_cv) || !length(fit_cv$cvm)) next

      lam <- as.numeric(fit_cv$lambda)
      cvm <- as.numeric(fit_cv$cvm)
      cvs <- as.numeric(fit_cv$cvsd)

      if (is.null(lam_ref)) {
        lam_ref <- lam
        cvms_list[[r]]  <- cvm
        cvSEs_list[[r]] <- cvs
      } else {
        # Align to the reference lambda vector, in case of minor path differences
        idx <- match(lam_ref, lam)
        cvm_aligned <- ifelse(is.na(idx), NA_real_, cvm[idx])
        cvs_aligned <- ifelse(is.na(idx), NA_real_, cvs[idx])
        cvms_list[[r]]  <- cvm_aligned
        cvSEs_list[[r]] <- cvs_aligned
      }
    }

    if (is.null(lam_ref)) return(NULL)
    cvms  <- do.call(cbind, cvms_list)
    cvSEs <- do.call(cbind, cvSEs_list)

    mean_cvm <- rowMeans(cvms, na.rm = TRUE)
    k_mean <- apply(cvms,  1L, function(z) sum(is.finite(z)))
    k_se   <- apply(cvSEs, 1L, function(z) sum(is.finite(z)))
    if (!any(is.finite(mean_cvm))) return(NULL)

    sd_cvm <- apply(cvms, 1L, function(z) {
      z <- z[is.finite(z)]
      if (length(z) <= 1L) return(NA_real_)
      stats::sd(z)
    })

    SE_within <- rowMeans(cvSEs^2, na.rm = TRUE)
    if (repeats == 1L) {
      SE_combined <- sqrt(SE_within / pmax(1L, k_se))
      SE_combined[!is.finite(SE_combined)] <- 0
    } else {
      SE_between  <- sd_cvm^2
      SE_combined <- sqrt( (SE_within / pmax(1L, k_se)) + (SE_between / pmax(1L, k_mean)) )
      SE_combined[!is.finite(SE_combined)] <- 0
    }

    idx_min <- which.min(mean_cvm)
    se_tol  <- if (is.finite(SE_combined[idx_min])) SE_combined[idx_min] else 0
    cand <- which(mean_cvm <= mean_cvm[idx_min] + se_tol)
    idx_1se <- if (length(cand)) max(cand) else idx_min

    list(lambda = lam_ref,
         mean_cvm = mean_cvm,
         sd_cvm = sd_cvm,
         se_combined = SE_combined,
         n_repeats = k_mean,
         idx_min = idx_min,
         idx_1se = idx_1se,
         fallback_count = fallback_count)
  }

  # --- Across alphas ---
  per_alpha <- list()
  for (a in glmnet_alpha) {
    agg <- fit_alpha_repeated(a)
    if (!is.null(agg)) per_alpha[[as.character(a)]] <- agg
  }

  # --- Defensive intercept-only check before ridge fallback ---
  if (ncol(X) == 0L) {
    intercept <- if (!is.null(w)) sum(w * y) / sum(w) else mean(y)
    terms_clean <- attr(mf, "terms"); environment(terms_clean) <- baseenv()

    # Record factor AND character levels
    xlevels <- list()
    vars_in_terms <- base::all.vars(stats::delete.response(terms_clean))
    for (v in vars_in_terms) {
      if (v %in% colnames(mf)) {
        if (is.factor(mf[[v]])) {
          xlevels[[v]] <- levels(mf[[v]])
        } else if (is.character(mf[[v]])) {
          xlevels[[v]] <- sort(unique(as.character(mf[[v]])))
        }
      }
    }
    contrasts_used <- NULL
    feature_names <- character(0)
    terms_str <- tryCatch(paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
                          error = function(e) NA_character_)
    safe_hash <- function(s) {
      if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
      bytes <- charToRaw(paste0(s, collapse = ""))
      sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
    }
    schema <- list(
      feature_names = feature_names,
      terms_str     = terms_str,
      xlevels       = xlevels,
      contrasts     = contrasts_used,
      terms_hash    = safe_hash(terms_str)
    )
    return(list(
      parms = c("(Intercept)" = intercept),
      glmnet_alpha = glmnet_alpha,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      y_pred = rep(intercept, nrow(X)),
      debias_fit = NULL,
      y_pred_debiased = NULL,
      cv_summary = list(),
      formula = formula,
      terms = terms_clean,
      training_X = X,
      actual_y = y,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = c("intercept_only", drop_msg, "ridge_fallback_bypassed"),
      meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                  family = fam, relaxed = relaxed, relax_cv_fallbacks = 0L, cv_object = NULL)
    ))
  }

  # --- Ridge fallback if everything failed (very rare) ---
  if (!length(per_alpha)) {
    warning("All cv.glmnet attempts failed; switching to ridge fallback with manual CV.")
    if (!is.null(seed)) set.seed(seed)
    foldid <- sample(rep(seq_len(nfolds_eff), length.out = n))
    lam_seq <- 10 ^ seq(3, -5, length.out = 120)

    mse_num <- mse_den <- rep(0, length(lam_seq))

    for (j in seq_along(lam_seq)) {
      lam <- lam_seq[j]
      for (fold in seq_len(nfolds_eff)) {
        tr_idx <- which(foldid != fold); te_idx <- which(foldid == fold)
        glmnet_fit_args_fold <- glmnet_fit_args
        if (!is.null(w)) glmnet_fit_args_fold$weights <- w[tr_idx]
        if (!is.null(off)) glmnet_fit_args_fold$offset <- off[tr_idx]
        fit_j <- tryCatch(
          do.call(glmnet::glmnet, c(
            list(x = X[tr_idx, , drop = FALSE], y = y[tr_idx],
                 alpha = 0, lambda = lam),
            glmnet_fit_args_fold
          )),
          error = function(e) NULL
        )
        if (is.null(fit_j)) next
        preds_te <- drop(stats::predict(fit_j, newx = X[te_idx, , drop = FALSE],
                                        s = lam, type = "response"))
        if (is.null(w)) {
          mse_num[j] <- mse_num[j] + sum((preds_te - y[te_idx])^2)
          mse_den[j] <- mse_den[j] + length(te_idx)
        } else {
          w_te <- w[te_idx]
          mse_num[j] <- mse_num[j] + sum(w_te * (preds_te - y[te_idx])^2)
          mse_den[j] <- mse_den[j] + sum(w_te)
        }
      }
    }
    mse_by_lambda <- mse_num / mse_den
    j_best <- which.min(mse_by_lambda)
    best_lambda <- lam_seq[j_best]

    final_ridge <- do.call(glmnet::glmnet, c(
      list(x = X, y = y, alpha = 0, lambda = best_lambda),
      glmnet_fit_args
    ))
    coef_mat   <- as.matrix(stats::coef(final_ridge, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_ridge, newx = X, s = best_lambda, type = "response"))

    debias_fit <- if (is_gaussian && length(y) >= 10 && stats::var(y_pred) > 0) stats::lm(y ~ y_pred) else NULL
    y_pred_debiased <- if (!is.null(debias_fit)) stats::predict(debias_fit) else NULL

    # Build schema
    terms_clean <- attr(mf, "terms"); environment(terms_clean) <- baseenv()
    feature_names <- colnames(X)

    xlevels <- list()
    vars_in_terms <- base::all.vars(stats::delete.response(terms_clean))
    for (v in vars_in_terms) {
      if (v %in% colnames(mf)) {
        if (is.factor(mf[[v]])) {
          xlevels[[v]] <- levels(mf[[v]])
        } else if (is.character(mf[[v]])) {
          xlevels[[v]] <- sort(unique(as.character(mf[[v]])))
        }
      }
    }
    contrasts_used <- attr(X, "contrasts")
    terms_str <- tryCatch(paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
                          error = function(e) NA_character_)
    safe_hash <- function(s) {
      if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
      bytes <- charToRaw(paste0(s, collapse = ""))
      sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
    }
    schema <- list(
      feature_names = feature_names,
      terms_str     = terms_str,
      xlevels       = xlevels,
      contrasts     = contrasts_used,
      terms_hash    = safe_hash(terms_str)
    )

    return(list(
      parms = best_coefs,
      glmnet_alpha = glmnet_alpha,
      best_alpha = 0,
      best_lambda = best_lambda,
      y_pred = y_pred,
      debias_fit = debias_fit,
      y_pred_debiased = y_pred_debiased,
      cv_summary = list(),
      formula = formula,
      terms = terms_clean,
      training_X = X,
      actual_y = y,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = c("ridge_fallback", drop_msg),
      meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                  family = fam, relaxed = relaxed, relax_cv_fallbacks = NA_integer_, cv_object = NULL)
    ))
  }

  # --- Pick alpha/lambda (1se or min) ---
  pick_for_alpha <- function(agg) {
    j <- if (identical(choose_rule, "min")) agg$idx_min else agg$idx_1se
    list(idx = j, score = agg$mean_cvm[j], lambda = agg$lambda[j])
  }
  alpha_scores <- lapply(per_alpha, pick_for_alpha)
  best_idx <- which.min(vapply(alpha_scores, function(x) x$score, numeric(1)))
  best_alpha <- as.numeric(names(alpha_scores)[best_idx])
  best_lambda <- alpha_scores[[best_idx]]$lambda

  # --- Final refit & outputs ---
  note_relax <- NULL
  cv_obj_to_keep <- NULL

  if (isTRUE(relaxed)) {
    reserved_cv <- c("x","y","alpha","lambda","foldid","exclude")
    cv_one <- do.call(glmnet::cv.glmnet,
                      c(list(x = X, y = y,
                             alpha  = best_alpha,
                             foldid = foldids[[1]],
                             exclude = exclude),
                        drop_reserved(cv_base_args, reserved_cv)))
    coef_mat   <- as.matrix(stats::coef(cv_one, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    mm <- cbind("(Intercept)" = 1, X)
    y_pred <- as.numeric(mm %*% best_coefs[match(colnames(mm), names(best_coefs))])
    note_relax <- "coefs from cv.glmnet (relaxed)"
    if (isTRUE(is.null(dots$keep)) || isTRUE(dots$keep)) cv_obj_to_keep <- cv_one
  } else {
    final_fit <- do.call(glmnet::glmnet, c(
      list(x = X, y = y, alpha = best_alpha, lambda = best_lambda),
      glmnet_fit_args
    ))
    coef_mat   <- as.matrix(stats::coef(final_fit, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_fit, newx = X, s = best_lambda, type = "response"))
  }

  debias_fit <- if (is_gaussian && length(y) >= 10 && stats::var(y_pred) > 0) stats::lm(y ~ y_pred) else NULL
  y_pred_debiased <- if (!is.null(debias_fit)) stats::predict(debias_fit) else NULL

  cv_summary <- lapply(per_alpha, function(agg) {
    data.frame(
      lambda      = agg$lambda,
      mean_cvm    = agg$mean_cvm,
      sd_cvm      = agg$sd_cvm,
      se_combined = agg$se_combined,
      n_repeats   = agg$n_repeats,
      idx_min     = agg$idx_min,
      idx_1se     = agg$idx_1se
    )
  })

  total_fallbacks <- sum(vapply(per_alpha, function(a) a$fallback_count, numeric(1)), na.rm = TRUE)

  # --- Clean terms and build schema for predict_cv() ---
  terms_clean <- attr(mf, "terms")
  environment(terms_clean) <- baseenv()

  feature_names <- colnames(X)
  xlevels <- list()
  vars_in_terms <- base::all.vars(stats::delete.response(terms_clean))
  for (v in vars_in_terms) {
    if (v %in% colnames(mf)) {
      if (is.factor(mf[[v]])) {
        xlevels[[v]] <- levels(mf[[v]])
      } else if (is.character(mf[[v]])) {
        xlevels[[v]] <- sort(unique(as.character(mf[[v]])))
      }
    }
  }
  contrasts_used <- attr(X, "contrasts")
  terms_str <- tryCatch(paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
                        error = function(e) NA_character_)
  safe_hash <- function(s) {
    if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
    bytes <- charToRaw(paste0(s, collapse = ""))
    sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
  }
  schema <- list(
    feature_names = feature_names,
    terms_str     = terms_str,
    xlevels       = xlevels,
    contrasts     = contrasts_used,
    terms_hash    = safe_hash(terms_str)
  )

  result<-list(
    parms = best_coefs,
    glmnet_alpha = glmnet_alpha,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    y_pred = y_pred,
    debias_fit = debias_fit,
    y_pred_debiased = y_pred_debiased,
    cv_summary = cv_summary,
    formula = formula,
    terms = terms_clean,
    training_X = X,
    actual_y = y,
    xlevels = xlevels,
    contrasts = contrasts_used,
    schema = schema,
    note = c(drop_msg, note_relax),
    meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                family = fam, relaxed = relaxed, relax_cv_fallbacks = total_fallbacks,
                cv_object = cv_obj_to_keep)
  )

  class(result) <- c("svem_cv", "svem_model")  # or just "svem_cv"

  return(result)
}
