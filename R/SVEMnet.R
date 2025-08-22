#' Fit an SVEMnet Model
#'
#' Wrapper for 'glmnet' (Friedman et al. 2010) to fit an ensemble of Elastic Net
#' models using the Self-Validated Ensemble Model method (SVEM, Lemkus et al. 2021).
#' Allows searching over multiple alpha values in the Elastic Net penalty.
#'
#' @param formula A formula specifying the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param nBoot Number of bootstrap iterations (default is 300).
#' @param glmnet_alpha Elastic Net mixing parameter(s) (default is \code{c(1)}).
#'   Can be a vector of alpha values, where \code{alpha = 1} is Lasso and \code{alpha = 0} is Ridge.
#' @param weight_scheme Weighting scheme for SVEM (default "SVEM").
#'   One of \code{"SVEM"}, \code{"FWR"}, or \code{"Identity"}.
#' @param objective Objective used to pick \code{lambda} on each bootstrap path (default "wAIC").
#'   One of \code{"wAIC"}, \code{"wBIC"}, \code{"wGIC"}, or \code{"wSSE"}.
#' @param gamma Penalty weight used only when \code{objective="wGIC"} (numeric, default \code{2}).
#'  Setting \code{gamma = 2} reproduces wAIC.
#' @param standardize Logical; passed to \code{glmnet} (default \code{TRUE}).
#' @param ... Additional args to \code{glmnet()}.
#'
#' @return An object of class \code{svem_model}.
#'
#' @details
#' The Self-Validated Ensemble Model (SVEM, Lemkus et al., 2021) framework provides a bootstrap
#' approach to improve predictions from base learners, including Elastic Net regression
#' as implemented in \code{glmnet}. In each of the \code{nBoot} iterations, SVEMnet applies
#' random exponentially distributed weights to the observations; anti-correlated weights are
#' used for validation when \code{weight_scheme="SVEM"}.
#'
#' SVEMnet allows \code{glmnet_alpha} to be a vector, enabling a search over multiple Elastic Net
#' mixing parameters within each bootstrap. The \code{objective} controls how the validation
#' criterion balances fit and complexity:
#'
#' \describe{
#'   \item{\code{"wSSE"}}{Weighted Sum of Squared Errors: uses the weighted validation SSE directly.}
#'   \item{\code{"wAIC"}}{Weighted AIC (Gaussian): \code{AIC = n * log(SSE_w / n) + 2 * k},
#'         where \code{n = sum(w_valid)} (after normalization) and \code{k} counts parameters
#'         including the intercept. Candidates require \code{k < n}.}
#'   \item{\code{"wBIC"}}{Weighted BIC-like criterion: \code{n * log(SSE_w / n) + log(n_eff) * k},
#'         with \code{n_eff = (sum(w_valid)^2) / sum(w_valid^2)} (Kish). For stability, \code{n_eff}
#'         is clipped to \code{[5, n]}. Candidates require \code{k < n} and \code{k < n_eff - 1}.}
#'   \item{\code{"wGIC"}}{Weighted Generalized Information Criterion:
#'         \code{n * log(SSE_w / n) + gamma * k}. Here \code{gamma} is a fixed nonnegative number.
#'         For robustness near the boundary, candidates require \code{k < n} and \code{k < n_eff - 1}.}
#' }
#'
#' \strong{Note on BIC:} In reweighted validation, information content varies with weight
#' heterogeneity; using \code{log(n_eff)} adapts the penalty to that effective size. With uniform
#' weights (Identity), \code{n_eff = n} and you recover standard BIC.
#'
#' A debiased fit is output (along with the standard fit). This is provided to allow the user to
#' match the output of JMP, which returns a debiased fit whenever \code{nBoot >= 10}.
#' The debiasing coefficients are always calculated by \code{SVEMnet()}, and the
#' \code{predict()} method determines whether the raw or debiased predictions are returned via its
#' \code{debias} argument. The default is \code{debias = FALSE}, based on performance on unpublished simulations.
#'
#' The returned object includes averaged coefficients (\code{parms}), debiased coefficients
#' (\code{parms_debiased}), the calibration fit (\code{debias_fit}), per-bootstrap coefficients,
#' chosen alphas and lambdas, the chosen \code{objective} (and \code{gamma} if applicable), and a
#' compact \code{diagnostics} list (median/IQR of selected model size and alpha frequencies).
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the
#' structure of some of the code and the roxygen documentation. The code for the significance test
#' is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @references
#' Gotwalt, C., & Ramsey, P. (2018). Model Validation Strategies for Designed Experiments Using
#' Bootstrapping Techniques With Applications to Biopharmaceuticals. \emph{JMP Discovery Conference}.
#' \url{https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849873/redirect_from_archived_page/true}
#'
#' Karl, A. T. (2024). A randomized permutation whole-model test heuristic for Self-Validated
#' Ensemble Models (SVEM). \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{249}, 105122.
#' \doi{10.1016/j.chemolab.2024.105122}
#'
#' Karl, A., Wisnowski, J., & Rushing, H. (2022). JMP Pro 17 Remedies for Practical Struggles with
#' Mixture Experiments. JMP Discovery Conference. \doi{10.13140/RG.2.2.34598.40003/1}
#'
#' Lemkus, T., Gotwalt, C., Ramsey, P., & Weese, M. L. (2021). Self-Validated Ensemble Models for
#' Design of Experiments. \emph{Chemometrics and Intelligent Laboratory Systems}, 219, 104439.
#' \doi{10.1016/j.chemolab.2021.104439}
#'
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020). Applications of the
#' Fractional-Random-Weight Bootstrap. \emph{The American Statistician}, 74(4), 345–358.
#' \doi{10.1080/00031305.2020.1731599}
#'
#' Ramsey, P., Gaudard, M., & Levin, W. (2021). Accelerating Innovation with Space Filling Mixture
#' Designs, Neural Networks and SVEM. \emph{JMP Discovery Conference}.
#' \url{ https://community.jmp.com/t5/Abstracts/Accelerating-Innovation-with-Space-Filling-Mixture-Designs/ev-p/756841}
#'
#' Ramsey, P., & Gotwalt, C. (2018). Model Validation Strategies for Designed Experiments Using
#' Bootstrapping Techniques With Applications to Biopharmaceuticals. \emph{JMP Discovery Conference - Europe}.
#' \url{https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849647/redirect_from_archived_page/true}
#'
#' Ramsey, P., Levin, W., Lemkus, T., & Gotwalt, C. (2021). SVEM: A Paradigm Shift in Design and
#' Analysis of Experiments. \emph{JMP Discovery Conference - Europe}.
#' \url{https://community.jmp.com/t5/Abstracts/SVEM-A-Paradigm-Shift-in-Design-and-Analysis-of-Experiments-2021/ev-p/756634}
#'
#' Ramsey, P., & McNeill, P. (2023). CMC, SVEM, Neural Networks, DOE, and Complexity:
#' It’s All About Prediction. \emph{JMP Discovery Conference}.
#'
#' Friedman, J. H., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized
#' Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 33(1), 1–22.
#'
#' @examples
#' # Simulate data
#' set.seed(0)
#' n <- 21
#' X1 <- runif(n)
#' X2 <- runif(n)
#' X3 <- runif(n)
#' y <- 1 + 2*X1 + 3*X2 + X1*X2 + X1^2  + rnorm(n)
#' data <- data.frame(y, X1, X2, X3)
#'
#' # Fit the SVEMnet model with a formula
#' model <- SVEMnet(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2) + I(X3^2),
#'   glmnet_alpha = c(1),
#'   data = data
#' )
#' coef(model)
#' plot(model)
#' predict(model, data)
#'
#' # Example: BIC-like penalty
#' # model_bic <- SVEMnet(y ~ X1 + X2 + X3, data = data, objective = "wBIC")
#' # Example: GIC with custom gamma
#' # model_gic <- SVEMnet(y ~ X1 + X2 + X3, data = data, objective = "wGIC", gamma = 4)
#'
#' @importFrom stats runif lm predict coef var model.frame model.matrix model.response
#' @importFrom glmnet glmnet
#' @export
SVEMnet <- function(formula, data, nBoot = 300, glmnet_alpha = c(1),
                    weight_scheme = c("SVEM", "FWR", "Identity"),
                    objective = c("wAIC", "wBIC", "wGIC", "wSSE"),
                    gamma = 2,
                    standardize = TRUE, ...) {

  ## ---- arg checks ----
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)

  if (!is.numeric(nBoot) || length(nBoot) != 1L || nBoot < 1 || !is.finite(nBoot)) {
    stop("nBoot must be a single finite number >= 1.")
  }
  nBoot <- as.integer(nBoot)

  if (length(glmnet_alpha) < 1L) stop("'glmnet_alpha' must contain at least one value.")
  if (!is.numeric(glmnet_alpha) || any(!is.finite(glmnet_alpha))) {
    stop("'glmnet_alpha' must be numeric and finite.")
  }
  if (any(glmnet_alpha < 0 | glmnet_alpha > 1)) {
    stop("'glmnet_alpha' values must be in [0, 1].")
  }
  glmnet_alpha <- unique(glmnet_alpha)

  ## ---- model frame / matrix (single NA policy) ----
  data <- as.data.frame(data)
  mf   <- model.frame(formula, data, na.action = stats::na.omit)
  if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")

  y <- model.response(mf)
  X <- model.matrix(formula, mf)

  ## Remove intercept column (glmnet adds its own)
  int_idx <- which(colnames(X) %in% c("(Intercept)", "Intercept"))
  if (length(int_idx)) X <- X[, -int_idx, drop = FALSE]
  if (ncol(X) == 0L) stop("SVEMnet requires at least one predictor.")

  y_numeric <- as.numeric(y)
  if (any(!is.finite(y_numeric)) || any(!is.finite(X))) {
    stop("Non-finite values in response/predictors after NA handling.")
  }
  storage.mode(X) <- "double"

  n <- nrow(X); p <- ncol(X)
  nobs  <- n
  nparm <- p + 1L

  ## ---- wGIC sanity checks now that n is known ----
  if (objective == "wGIC") {
    if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
      stop("'gamma' must be a single finite numeric value >= 0 for wGIC.")
    }
    if (gamma == 0) warning("gamma = 0 implies no complexity penalty (wGIC ~ wSSE); risk of overfitting.")
    if (gamma < 0.2) warning("gamma is very small (< 0.2); selection may behave close to wSSE.")
    bic_like <- log(max(2, n))
    if (gamma > 3 * bic_like) warning("gamma is much larger than a BIC-like penalty (3*log(n)); selection may underfit.")
  }

  ## ---- containers ----
  coef_matrix  <- matrix(NA_real_, nrow = nBoot, ncol = p + 1L)
  colnames(coef_matrix) <- c("(Intercept)", colnames(X))
  best_alphas  <- rep(NA_real_, nBoot)
  best_lambdas <- rep(NA_real_, nBoot)
  k_sel_vec    <- rep(NA_integer_, nBoot)

  ## ---- bootstrap loop ----
  for (i in seq_len(nBoot)) {
    ## fractional weights
    eps <- .Machine$double.eps
    if (weight_scheme == "SVEM") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U)
      w_valid <- -log1p(-U)
    } else if (weight_scheme == "FWR") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U)
      w_valid <- w_train
    } else { # "Identity"
      w_train <- rep(1, n)
      w_valid <- rep(1, n)
    }
    ## normalize to sum to n
    w_train <- w_train * (n / sum(w_train))
    w_valid <- w_valid * (n / sum(w_valid))

    ## trackers per bootstrap
    best_val_score <- Inf
    best_alpha     <- NA_real_
    best_lambda    <- NA_real_
    best_beta_hat  <- rep(NA_real_, p + 1L)
    best_k         <- NA_integer_

    ## ---- search over alpha ----
    for (alpha in glmnet_alpha) {
      fit <- tryCatch({
        withCallingHandlers({
          glmnet::glmnet(
            x = X, y = y_numeric,
            alpha = alpha,
            weights = w_train,
            intercept = TRUE,
            standardize = standardize,
            nlambda = 500,
            maxit = 1e6,
            ...
          )
        }, warning = function(w) invokeRestart("muffleWarning"))
      }, error = function(e) NULL)
      if (is.null(fit) || !length(fit$lambda) || !length(fit$df)) next

      ## predictions across the path
      pred_valid <- tryCatch({
        as.matrix(withCallingHandlers(stats::predict(fit, newx = X),
                                      warning = function(w) invokeRestart("muffleWarning")))
      }, error = function(e) NULL)
      if (is.null(pred_valid) || nrow(pred_valid) != n) next

      ## weighted SSE along the path
      res <- pred_valid - y_numeric
      val_errors <- as.vector(crossprod(w_valid, res^2))   # length = n_lambda
      val_errors[!is.finite(val_errors)] <- Inf
      adj_val_errors <- pmax(val_errors, .Machine$double.eps)

      ## degrees of freedom: df from glmnet excludes intercept -> add 1
      k_raw <- fit$df + 1L
      if (length(k_raw) != length(adj_val_errors)) next

      ## sums of validation weights and n_eff
      sumw  <- sum(w_valid)
      sumw2 <- sum(w_valid^2)
      n_eff <- (sumw^2) / (sumw2 + .Machine$double.eps)
      n_eff <- max(5, min(n, n_eff))  # clip to [5, n]

      ## weighted MSE (SSE_w / sum of weights)
      mse_w <- adj_val_errors / sumw

      ## conservative k to avoid near-boundary pathologies; keep >= 1 and <= n - 2
      k_eff <- pmax(1L, pmin(k_raw, n - 2L))

      ## objective metric
      if (objective == "wSSE") {

        metric <- val_errors

      } else if (objective == "wAIC") {

        metric <- rep(Inf, length(k_raw))
        mask <- (k_raw < n)
        metric[mask] <- sumw * log(mse_w[mask]) + 2 * k_eff[mask]

      } else if (objective == "wBIC") {

        metric <- rep(Inf, length(k_raw))
        mask <- (k_raw < n) & (k_eff < (n_eff - 1))
        metric[mask] <- sumw * log(mse_w[mask]) + log(n_eff) * k_eff[mask]

      } else {  # "wGIC"

        metric <- rep(Inf, length(k_raw))
        mask <- (k_raw < n) & (k_eff < (n_eff - 1))
        metric[mask] <- sumw * log(mse_w[mask]) + gamma * k_eff[mask]
      }

      metric[!is.finite(metric)] <- Inf
      if (all(!is.finite(metric))) next

      idx_min <- which.min(metric)
      lambda_opt <- fit$lambda[idx_min]
      val_score  <- metric[idx_min]

      if (is.finite(val_score) && val_score < best_val_score) {
        beta <- tryCatch({ as.matrix(stats::coef(fit, s = lambda_opt)) }, error = function(e) NULL)
        if (!is.null(beta)) {
          best_val_score <- val_score
          best_alpha     <- alpha
          best_lambda    <- lambda_opt
          best_beta_hat  <- drop(beta)  # (p+1), includes intercept
          best_k         <- k_raw[idx_min]
        }
      }
    } # alpha loop

    ## fallback this bootstrap if needed
    if (anyNA(best_beta_hat) || !all(is.finite(best_beta_hat))) {
      best_beta_hat <- c(mean(y_numeric), rep(0, p))
      best_alpha    <- NA_real_
      best_lambda   <- NA_real_
      best_k        <- 1L
    }

    coef_matrix[i, ] <- best_beta_hat
    best_alphas[i]   <- best_alpha
    best_lambdas[i]  <- best_lambda
    k_sel_vec[i]     <- best_k
  } # bootstrap loop

  ## ---- finalize ----
  valid_rows <- rowSums(!is.finite(coef_matrix)) == 0
  if (!any(valid_rows)) stop("All bootstrap iterations failed to produce valid coefficients.")
  coef_matrix  <- coef_matrix[valid_rows, , drop = FALSE]
  best_alphas  <- best_alphas[valid_rows]
  best_lambdas <- best_lambdas[valid_rows]
  k_sel_vec    <- k_sel_vec[valid_rows]

  avg_coefficients <- colMeans(coef_matrix)
  y_pred <- as.vector(X %*% avg_coefficients[-1L] + avg_coefficients[1L])

  debias_fit <- NULL
  y_pred_debiased <- NULL
  if (nBoot >= 10 && stats::var(y_pred) > 0) {
    debias_fit <- stats::lm(y_numeric ~ y_pred)
    y_pred_debiased <- stats::predict(debias_fit)
  }

  ## ---- debiased coefficients (beta' = a + b*beta0, b*beta) ----
  parms_debiased <- avg_coefficients
  if (!is.null(debias_fit)) {
    ab <- try(stats::coef(debias_fit), silent = TRUE)
    if (!inherits(ab, "try-error") && length(ab) >= 2 &&
        is.finite(ab[1]) && is.finite(ab[2])) {
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
        if (length(parms_debiased) > 1L) {
          parms_debiased[-1] <- b * parms_debiased[-1]
        }
      }
    }
  }

  diagnostics <- list(
    k_summary = c(k_median = stats::median(k_sel_vec, na.rm = TRUE),
                  k_iqr    = stats::IQR(k_sel_vec, na.rm = TRUE))
  )
  tf <- table(best_alphas[is.finite(best_alphas)])
  diagnostics$alpha_freq <- if (length(tf)) as.numeric(tf) / sum(tf) else numeric()

  result <- list(
    parms            = avg_coefficients,
    parms_debiased   = parms_debiased,
    debias_fit       = debias_fit,
    coef_matrix      = coef_matrix,
    nBoot            = nBoot,
    glmnet_alpha     = glmnet_alpha,
    best_alphas      = best_alphas,
    best_lambdas     = best_lambdas,
    weight_scheme    = weight_scheme,
    objective        = objective,
    gamma            = if (objective == "wGIC") gamma else NA_real_,
    diagnostics      = diagnostics,
    actual_y         = y_numeric,
    training_X       = X,
    y_pred           = y_pred,
    y_pred_debiased  = y_pred_debiased,
    nobs             = nobs,
    nparm            = nparm,
    formula          = formula,
    terms            = attr(mf, "terms")
  )
  class(result) <- c("svem_model", "SVEMnet")
  result
}
