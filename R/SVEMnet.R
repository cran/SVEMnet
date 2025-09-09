#' Fit an SVEMnet Model (with optional relaxed elastic net)
#'
#' Wrapper for glmnet (Friedman et al. 2010) to fit an ensemble of Elastic Net
#' models using the Self-Validated Ensemble Model method (SVEM; Lemkus et al. 2021),
#' with an option to use glmnet's built-in relaxed elastic net (Meinshausen 2007).
#' Supports searching over multiple alpha values in the Elastic Net penalty.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param nBoot Number of bootstrap iterations (default 200).
#' @param glmnet_alpha Elastic Net mixing parameter(s). May be a vector with entries in the range between 0 and 1, inclusive,
#'   where alpha = 1 is Lasso and alpha = 0 is Ridge. Defaults to c(0.25, 0.5, 0.75, 1).
#' @param weight_scheme Weighting scheme for SVEM (default "SVEM").
#'   One of "SVEM", "FWR", or "Identity".
#' @param objective Objective used to pick lambda on each bootstrap path
#'   (default "auto"). One of "auto", "wAIC", "wBIC", "wGIC", or "wSSE".
#' @param auto_ratio_cutoff Single cutoff for the automatic rule when
#'   objective = "auto" (default 1.3). Let r = n_X / p_X, where n_X is the
#'   number of training rows and p_X is the number of predictors in the model
#'   matrix after dropping the intercept column. If r >= auto_ratio_cutoff,
#'   SVEMnet uses wAIC; otherwise it uses wBIC.
#' @param gamma Penalty weight used only when objective = "wGIC" (numeric, default 2).
#'   Setting gamma = 2 reproduces wAIC. Larger values approach a BIC-like penalty.
#' @param relaxed Logical, TRUE or FALSE (default TRUE). When TRUE, use glmnet's
#'   relaxed elastic net path and select both lambda and relaxed gamma on each bootstrap.
#'   When FALSE, fit the standard glmnet path (equivalent to gamma = 1).
#'   Note: if relaxed = TRUE and glmnet_alpha includes 0 (ridge), alpha = 0 is dropped.
#' @param ... Additional args passed to glmnet() (e.g., penalty.factor,
#'   lower.limits, upper.limits, offset, standardize.response, etc.).
#'   Any user-supplied weights are ignored so SVEM can supply its own bootstrap weights.
#'   Any user-supplied standardize is ignored; SVEMnet always uses standardize = TRUE.
#'
#' @return An object of class svem_model with elements:
#' \itemize{
#'   \item parms: averaged coefficients (including intercept).
#'   \item parms_debiased: averaged coefficients adjusted by the calibration fit.
#'   \item debias_fit: lm(y ~ y_pred) calibration model used for debiasing (or NULL).
#'   \item coef_matrix: per-bootstrap coefficient matrix.
#'   \item nBoot, glmnet_alpha, best_alphas, best_lambdas, weight_scheme, relaxed.
#'   \item best_relax_gammas: per-bootstrap relaxed gamma chosen (NA if relaxed = FALSE).
#'   \item objective_input, objective_used, objective (same as objective_used),
#'         auto_used, auto_decision, auto_rule, gamma (when wGIC).
#'   \item dropped_alpha0_for_relaxed: whether alpha = 0 was removed because relaxed = TRUE.
#'   \item schema: list(feature_names, terms_str, xlevels, contrasts, terms_hash) for safe predict.
#'   \item sampling_schema: list(
#'         predictor_vars, var_classes,
#'         num_ranges = rbind(min=..., max=...) for numeric raw predictors,
#'         factor_levels = list(...) for factor/character raw predictors).
#'   \item diagnostics: list with k_summary (median and IQR of selected size),
#'         fallback_rate, n_eff_summary, alpha_freq, relax_gamma_freq.
#'   \item actual_y, training_X, y_pred, y_pred_debiased, nobs, nparm, formula, terms,
#'         xlevels, contrasts.
#' }
#'
#' @details
#' SVEM applies fractional bootstrap weights to training data and anti-correlated
#' weights for validation when weight_scheme = "SVEM". For each bootstrap, glmnet
#' paths are fit for each alpha in glmnet_alpha, and the lambda (and, if relaxed = TRUE,
#' relaxed gamma) minimizing a weighted validation criterion is selected.
#'
#' Predictors are always standardized internally via glmnet::glmnet(..., standardize = TRUE).
#'
#' When relaxed = TRUE, coef(fit, s = lambda, gamma = g) is used to obtain the
#' coefficient path at each relaxed gamma in the internal grid. Metrics are computed
#' from validation-weighted errors and model size is taken as the number of nonzero
#' coefficients including the intercept (support size), keeping selection consistent
#' between standard and relaxed paths. When relaxed = FALSE, gamma is fixed at 1.
#'
#' Automatic objective rule ("auto"): This function uses a single ratio cutoff,
#' auto_ratio_cutoff, applied to r = n_X / p_X, where p_X is computed from
#' the model matrix with the intercept column removed. If r >= auto_ratio_cutoff
#' the selection criterion is wAIC; otherwise it is wBIC.
#'
#' Implementation notes for safety:
#' - The training terms object is stored with environment set to baseenv() to avoid
#'   accidental lookups in the calling environment.
#' - A compact schema (feature names, xlevels, contrasts) is stored to let predict()
#'   reconstruct the design matrix deterministically.
#' - A lightweight sampling schema (numeric ranges and factor levels for raw predictors)
#'   is cached to enable random-table generation without needing the original data.
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview for structuring parts of the code
#' and documentation.
#'
#' @references
#' Gotwalt, C., & Ramsey, P. (2018). Model Validation Strategies for Designed Experiments Using
#' Bootstrapping Techniques With Applications to Biopharmaceuticals. JMP Discovery Conference.
#' https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849873/redirect_from_archived_page/true
#'
#' Karl, A. T. (2024). A randomized permutation whole-model test heuristic for Self-Validated
#' Ensemble Models (SVEM). Chemometrics and Intelligent Laboratory Systems, 249, 105122.
#' doi:10.1016/j.chemolab.2024.105122
#'
#' Karl, A., Wisnowski, J., & Rushing, H. (2022). JMP Pro 17 Remedies for Practical Struggles with
#' Mixture Experiments. JMP Discovery Conference. doi:10.13140/RG.2.2.34598.40003/1
#'
#' Lemkus, T., Gotwalt, C., Ramsey, P., & Weese, M. L. (2021). Self-Validated Ensemble Models for
#' Design of Experiments. Chemometrics and Intelligent Laboratory Systems, 219, 104439.
#' doi:10.1016/j.chemolab.2021.104439
#'
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020). Applications of the
#' Fractional-Random-Weight Bootstrap. The American Statistician, 74(4), 345-358.
#' doi:10.1080/00031305.2020.1731599
#'
#' Ramsey, P., Gaudard, M., & Levin, W. (2021). Accelerating Innovation with Space Filling Mixture
#' Designs, Neural Networks and SVEM. JMP Discovery Conference.
#' https://community.jmp.com/t5/Abstracts/Accelerating-Innovation-with-Space-Filling-Mixture-Designs/ev-p/756841
#'
#' Ramsey, P., & Gotwalt, C. (2018). Model Validation Strategies for Designed Experiments Using
#' Bootstrapping Techniques With Applications to Biopharmaceuticals. JMP Discovery Conference - Europe.
#' https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849647/redirect_from_archived_page/true
#'
#' Ramsey, P., Levin, W., Lemkus, T., & Gotwalt, C. (2021). SVEM: A Paradigm Shift in Design and
#' Analysis of Experiments. JMP Discovery Conference - Europe.
#' https://community.jmp.com/t5/Abstracts/SVEM-A-Paradigm-Shift-in-Design-and-Analysis-of-Experiments-2021/ev-p/756634
#'
#' Ramsey, P., & McNeill, P. (2023). CMC, SVEM, Neural Networks, DOE, and Complexity:
#' It's All About Prediction. JMP Discovery Conference.
#'
#' Friedman, J. H., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized
#' Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22.
#'
#' Meinshausen, N. (2007).
#' Relaxed Lasso. Computational Statistics & Data Analysis, 52(1), 374-393.
#'
#' @importFrom stats runif lm predict coef var model.frame model.matrix model.response delete.response IQR median
#' @importFrom glmnet glmnet
#'
#' @examples
#' set.seed(42)
#'
#' n  <- 30
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X3 <- rnorm(n)
#' eps <- rnorm(n, sd = 0.5)
#' y <- 1 + 2*X1 - 1.5*X2 + 0.5*X3 + 1.2*(X1*X2) + 0.8*(X1^2) + eps
#' dat <- data.frame(y, X1, X2, X3)
#'
#' mod_relax <- SVEMnet(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2),
#'   data          = dat,
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 75,
#'   objective     = "auto",
#'   weight_scheme = "SVEM",
#'   relaxed       = FALSE #to run faster to pass CRAN checks
#' )
#'
#' pred_in_raw <- predict(mod_relax, dat, debias = FALSE)
#' pred_in_db  <- predict(mod_relax, dat, debias = TRUE)
#'
#' @export
SVEMnet <- function(formula, data, nBoot = 200, glmnet_alpha = c(0.25, 0.5, 0.75, 1),
                    weight_scheme = c("SVEM", "FWR", "Identity"),
                    objective = c("auto", "wAIC", "wBIC", "wGIC", "wSSE"),
                    auto_ratio_cutoff = 1.3,
                    gamma = 2,
                    relaxed = TRUE,
                    ...) {

  # ---- argument checks ----
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)

  if (!is.logical(relaxed) || length(relaxed) != 1L || is.na(relaxed)) {
    stop("'relaxed' must be a single logical TRUE or FALSE.")
  }
  relax_flag <- isTRUE(relaxed)

  if (!is.numeric(nBoot) || length(nBoot) != 1L || !is.finite(nBoot) || nBoot < 1) stop("nBoot must be >= 1.")
  nBoot <- as.integer(nBoot)

  if (!is.numeric(glmnet_alpha) || any(!is.finite(glmnet_alpha))) stop("'glmnet_alpha' must be numeric and finite.")
  if (any(glmnet_alpha < 0 | glmnet_alpha > 1)) stop("'glmnet_alpha' values must be in [0, 1].")
  glmnet_alpha <- unique(as.numeric(glmnet_alpha)); if (!length(glmnet_alpha)) glmnet_alpha <- 1

  if (!is.numeric(auto_ratio_cutoff) || length(auto_ratio_cutoff) != 1L ||
      !is.finite(auto_ratio_cutoff) || auto_ratio_cutoff <= 0) {
    stop("'auto_ratio_cutoff' must be a single finite number > 0.")
  }
  auto_ratio_cutoff <- as.numeric(auto_ratio_cutoff)

  # ---- model frame / matrix ----
  data <- as.data.frame(data)
  mf   <- stats::model.frame(formula, data, na.action = stats::na.omit)
  if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")

  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)

  # drop intercept column (glmnet adds its own)
  int_idx <- which(colnames(X) %in% c("(Intercept)", "Intercept"))
  if (length(int_idx)) X <- X[, -int_idx, drop = FALSE]
  if (ncol(X) == 0L) stop("SVEMnet requires at least one predictor.")

  if (!is.numeric(y)) stop("SVEMnet currently supports numeric (Gaussian) y only.")
  y_numeric <- as.numeric(y)

  if (any(!is.finite(y_numeric)) || any(!is.finite(X))) stop("Non-finite values in response/predictors after NA handling.")
  storage.mode(X) <- "double"

  n <- nrow(X); p <- ncol(X); nobs <- n; nparm <- p + 1L

  # capture training xlevels and contrasts for prediction
  terms_full  <- attr(mf, "terms")
  terms_clean <- terms_full; environment(terms_clean) <- baseenv()

  # raw predictor symbols (no response, no expanded terms names)
  predictor_vars <- base::all.vars(stats::delete.response(terms_full))

  # classes and levels on raw columns as present in 'mf'
  var_classes <- setNames(vapply(predictor_vars, function(v) {
    if (v %in% colnames(mf)) class(mf[[v]])[1] else NA_character_
  }, character(1)), predictor_vars)

  xlevels <- list()
  factor_levels <- list()
  for (v in predictor_vars) {
    if (v %in% colnames(mf)) {
      if (is.factor(mf[[v]])) {
        xlevels[[v]] <- levels(mf[[v]])
        factor_levels[[v]] <- levels(mf[[v]])
      } else if (is.character(mf[[v]])) {
        lev <- sort(unique(as.character(mf[[v]])))
        factor_levels[[v]] <- lev
      }
    }
  }
  contrasts_used <- attr(X, "contrasts")

  # numeric ranges on raw predictors
  num_vars <- predictor_vars[vapply(predictor_vars, function(v) v %in% colnames(mf) && is.numeric(mf[[v]]), logical(1))]
  if (length(num_vars)) {
    rng_mat <- vapply(num_vars, function(v) {
      r <- range(mf[[v]], na.rm = TRUE)
      if (!all(is.finite(r)) || r[1] == r[2]) r <- c(min(mf[[v]], na.rm = TRUE), max(mf[[v]], na.rm = TRUE))
      r
    }, numeric(2))
    rownames(rng_mat) <- c("min","max")
    num_ranges <- as.matrix(rng_mat)
  } else {
    num_ranges <- matrix(numeric(0), nrow = 2, ncol = 0, dimnames = list(c("min","max"), NULL))
  }

  # ---- objective selection (auto) ----
  auto_rule      <- c(ratio_cutoff = auto_ratio_cutoff)
  auto_used      <- identical(objective, "auto")
  auto_decision  <- NA_character_
  objective_used <- if (auto_used) {
    r_tmp <- n / p
    if (is.finite(r_tmp) && r_tmp >= auto_ratio_cutoff) { auto_decision <- "wAIC"; "wAIC" } else { auto_decision <- "wBIC"; "wBIC" }
  } else objective

  # ---- wGIC sanity checks ----
  if (objective_used == "wGIC") {
    if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) stop("'gamma' must be >= 0.")
    if (gamma == 0) warning("gamma = 0 implies no complexity penalty (wGIC ~ wSSE).")
    if (gamma < 0.2) warning("gamma is very small (< 0.2); selection may behave close to wSSE.")
    bic_like <- log(max(2, n))
    if (gamma > 3 * bic_like) warning("gamma is much larger than a BIC-like penalty (3*log(n)); selection may underfit.")
  }

  # ---- drop ridge when relaxed ----
  dropped_alpha0_for_relaxed <- FALSE
  if (isTRUE(relax_flag) && any(glmnet_alpha == 0)) {
    warning("Dropping alpha = 0 (ridge) for relaxed fits; ridge + relaxed is not supported here.")
    glmnet_alpha <- glmnet_alpha[glmnet_alpha != 0]; if (!length(glmnet_alpha)) glmnet_alpha <- 1
    dropped_alpha0_for_relaxed <- TRUE
  }

  # ---- containers ----
  coef_matrix  <- matrix(NA_real_, nrow = nBoot, ncol = p + 1L)
  colnames(coef_matrix) <- c("(Intercept)", colnames(X))
  best_alphas       <- rep(NA_real_, nBoot)
  best_lambdas      <- rep(NA_real_, nBoot)
  best_relax_gammas <- rep(NA_real_, nBoot)
  k_sel_vec         <- rep(NA_integer_, nBoot)
  fallbacks         <- integer(nBoot)
  n_eff_keep        <- numeric(nBoot)

  # ---- capture and sanitize user '...' so SVEM controls weights ----
  dots <- list(...)
  if (length(dots)) {
    if ("weights" %in% names(dots)) { warning("Ignoring user 'weights'; SVEM uses bootstrap weights."); dots$weights <- NULL }
    if ("standardize" %in% names(dots)) { warning("Ignoring user 'standardize'; SVEMnet always standardizes."); dots$standardize <- NULL }
    protected <- c("x","y","alpha","intercept","relax","standardize","nlambda","maxit","lambda.min.ratio","gamma","lambda")
    dots <- dots[setdiff(names(dots), protected)]
  }

  .support_size_one <- function(beta_col, base_tol = 1e-7, count_intercept = TRUE) {
    if (!is.numeric(beta_col) || length(beta_col) < 1L) return(NA_integer_)
    slope <- beta_col[-1L]
    tol_j <- base_tol * max(1, max(abs(slope)))
    k_slope <- sum(abs(slope) > tol_j)
    (if (count_intercept) 1L else 0L) + k_slope
  }

  # ---- bootstrap loop ----
  for (i in seq_len(nBoot)) {
    eps <- .Machine$double.eps
    if (weight_scheme == "SVEM") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- -log1p(-U)
    } else if (weight_scheme == "FWR") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- w_train
    } else { w_train <- rep(1, n); w_valid <- rep(1, n) }
    w_train <- w_train * (n / sum(w_train)); w_valid <- w_valid * (n / sum(w_valid))

    best_val_score <- Inf
    best_alpha     <- NA_real_
    best_lambda    <- NA_real_
    best_beta_hat  <- rep(NA_real_, p + 1L)
    best_k         <- NA_integer_
    best_relax_g   <- if (isTRUE(relax_flag)) NA_real_ else 1

    sumw  <- sum(w_valid); sumw2 <- sum(w_valid^2)
    n_eff_raw <- (sumw^2) / (sumw2 + eps)
    n_eff_adm <- max(2, min(n, n_eff_raw))
    n_eff_keep[i] <- n_eff_raw

    relax_gamma_grid <- if (isTRUE(relax_flag)) c(0, 0.25, 0.5, 0.75, 1) else 1

    for (alpha in glmnet_alpha) {
      fit <- tryCatch({
        withCallingHandlers({
          do.call(glmnet::glmnet, c(list(x = X, y = y_numeric,
                                         alpha = alpha,
                                         weights = w_train,
                                         intercept = TRUE,
                                         standardize = TRUE,
                                         nlambda = 500,
                                         maxit = 1e6,
                                         relax = isTRUE(relax_flag)), dots))
        }, warning = function(w) base::invokeRestart("muffleWarning"))
      }, error = function(e) NULL)
      if (is.null(fit) || !length(fit$lambda)) next

      for (relax_g in relax_gamma_grid) {
        coef_path <- tryCatch({
          if (isTRUE(relax_flag)) as.matrix(stats::coef(fit, s = fit$lambda, gamma = relax_g))
          else                     as.matrix(stats::coef(fit, s = fit$lambda))
        }, error = function(e) NULL)
        if (is.null(coef_path) || nrow(coef_path) != (p + 1L)) next
        L <- ncol(coef_path); if (L == 0L) next

        pred_valid <- tryCatch({
          XB <- X %*% coef_path[-1L, , drop = FALSE]
          sweep(XB, 2, coef_path[1L, ], FUN = "+")
        }, error = function(e) NULL)
        if (is.null(pred_valid) || nrow(pred_valid) != n) next

        res <- pred_valid - y_numeric
        val_errors <- as.vector(crossprod(w_valid, res^2))
        val_errors[!is.finite(val_errors)] <- Inf
        adj_val_errors <- pmax(val_errors, .Machine$double.eps)

        k_raw <- integer(ncol(coef_path))
        for (jj in seq_len(ncol(coef_path))) k_raw[jj] <- .support_size_one(coef_path[, jj], 1e-7, TRUE)
        if (length(k_raw) != length(adj_val_errors)) next

        n_like <- sumw
        mse_w  <- adj_val_errors / n_like
        k_slope <- pmax(0L, k_raw - 1L)
        k_eff   <- 1L + k_slope
        logN_pen <- log(n_eff_adm)

        metric <- switch(objective_used,
                         "wSSE" = val_errors,
                         "wAIC" = { out <- rep(Inf, L); mask <- (k_slope < n_eff_adm); out[mask] <- n_like * log(mse_w[mask]) + 2 * k_eff[mask]; out },
                         "wBIC" = { out <- rep(Inf, L); mask <- (k_slope < n_eff_adm); out[mask] <- n_like * log(mse_w[mask]) + logN_pen * k_eff[mask]; out },
                         { out <- rep(Inf, L); mask <- (k_slope < n_eff_adm); out[mask] <- n_like * log(mse_w[mask]) + gamma * k_eff[mask]; out })
        metric[!is.finite(metric)] <- Inf
        if (!any(is.finite(metric))) next

        idx_min    <- which.min(metric)
        lambda_opt <- fit$lambda[idx_min]
        val_score  <- metric[idx_min]

        if (is.finite(val_score) && val_score < best_val_score) {
          best_val_score <- val_score
          best_alpha     <- alpha
          best_lambda    <- lambda_opt
          best_beta_hat  <- coef_path[, idx_min]
          best_k         <- k_raw[idx_min]
          best_relax_g   <- relax_g
        }
      }
    }

    if (anyNA(best_beta_hat) || !all(is.finite(best_beta_hat))) {
      fallbacks[i]  <- 1L
      mu_w          <- sum(w_train * y_numeric) / sum(w_train)
      best_beta_hat <- c(mu_w, rep(0, p))
      best_alpha    <- NA_real_
      best_lambda   <- NA_real_
      best_k        <- 1L
      best_relax_g  <- if (isTRUE(relax_flag)) NA_real_ else 1
    }

    coef_matrix[i, ]     <- best_beta_hat
    best_alphas[i]       <- best_alpha
    best_lambdas[i]      <- best_lambda
    best_relax_gammas[i] <- best_relax_g
    k_sel_vec[i]         <- best_k
  }

  # ---- finalize ----
  valid_rows <- rowSums(!is.finite(coef_matrix)) == 0
  if (!any(valid_rows)) stop("All bootstrap iterations failed to produce valid coefficients.")
  coef_matrix        <- coef_matrix[valid_rows, , drop = FALSE]
  best_alphas       <- best_alphas[valid_rows]
  best_lambdas      <- best_lambdas[valid_rows]
  best_relax_gammas <- best_relax_gammas[valid_rows]
  k_sel_vec         <- k_sel_vec[valid_rows]
  fallbacks         <- fallbacks[valid_rows]
  n_eff_keep        <- n_eff_keep[valid_rows]

  avg_coefficients <- colMeans(coef_matrix)
  y_pred <- as.vector(X %*% avg_coefficients[-1L] + avg_coefficients[1L])

  debias_fit <- NULL
  y_pred_debiased <- NULL
  if (nBoot >= 10 && stats::var(y_pred) > 0) {
    debias_fit <- stats::lm(y_numeric ~ y_pred)
    y_pred_debiased <- stats::predict(debias_fit)
  }

  parms_debiased <- avg_coefficients
  if (!is.null(debias_fit)) {
    ab <- try(stats::coef(debias_fit), silent = TRUE)
    if (!inherits(ab, "try-error") && length(ab) >= 2 && is.finite(ab[1]) && is.finite(ab[2])) {
      a <- unname(ab[1]); b <- unname(ab[2])
      int_name <- "(Intercept)"
      if (!is.null(names(parms_debiased)) && int_name %in% names(parms_debiased)) {
        parms_debiased[int_name] <- a + b * parms_debiased[int_name]
        if (length(parms_debiased) > 1L) {
          slope_names <- setdiff(names(parms_debiased), int_name)
          parms_debiased[slope_names] <- b * parms_debiased[slope_names]
        }
      } else {
        parms_debiased[1] <- a + b * parms_debiased[1]
        if (length(parms_debiased) > 1L) parms_debiased[-1] <- b * parms_debiased[-1]
      }
    }
  }

  diagnostics <- list(
    k_summary = c(k_median = stats::median(k_sel_vec), k_iqr = stats::IQR(k_sel_vec)),
    fallback_rate = mean(fallbacks),
    n_eff_summary = summary(n_eff_keep),
    relax_gamma_freq = if (isTRUE(relax_flag)) prop.table(table(round(best_relax_gammas, 2))) else NULL
  )
  tf <- table(best_alphas[is.finite(best_alphas)])
  diagnostics$alpha_freq <- if (length(tf)) as.numeric(tf) / sum(tf) else numeric()

  # ---- schema for safe predict ----
  feature_names <- colnames(X)
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

  # ---- sampling schema for random-table generation ----
  sampling_schema <- list(
    predictor_vars = predictor_vars,
    var_classes    = var_classes,
    num_ranges     = num_ranges,
    factor_levels  = factor_levels
  )

  result <- list(
    parms             = avg_coefficients,
    parms_debiased    = parms_debiased,
    debias_fit        = debias_fit,
    coef_matrix       = coef_matrix,
    nBoot             = nBoot,
    glmnet_alpha      = glmnet_alpha,
    best_alphas       = best_alphas,
    best_lambdas      = best_lambdas,
    best_relax_gammas = if (isTRUE(relax_flag)) best_relax_gammas else rep(NA_real_, length(best_alphas)),
    weight_scheme     = weight_scheme,
    relaxed           = isTRUE(relax_flag),
    dropped_alpha0_for_relaxed = dropped_alpha0_for_relaxed,
    objective_input   = objective,
    objective_used    = objective_used,
    objective         = objective_used,
    auto_used         = auto_used,
    auto_decision     = if (auto_used) auto_decision else NA_character_,
    auto_rule         = auto_rule,
    gamma             = if (objective_used == "wGIC") gamma else NA_real_,
    diagnostics       = diagnostics,
    actual_y          = y_numeric,
    training_X        = X,
    y_pred            = y_pred,
    y_pred_debiased   = y_pred_debiased,
    nobs              = nobs,
    nparm             = nparm,
    formula           = formula,
    terms             = terms_clean,
    xlevels           = xlevels,
    contrasts         = contrasts_used,
    schema            = schema,
    sampling_schema   = sampling_schema
  )
  class(result) <- c("svem_model", "SVEMnet")
  result
}
