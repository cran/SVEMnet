#' Fit an SVEMnet model (Self-Validated Ensemble Elastic Net)
#'
#' Fit a Self-Validated Ensemble Model (SVEM) with elastic net or relaxed
#' elastic net base learners using \code{glmnet}. Fractional random-weight
#' (FRW) train/validation weights are drawn on each bootstrap replicate,
#' a validation-weighted information criterion (wAIC, wBIC, or wSSE) is
#' minimized to select the penalty, and predictions are ensembled across
#' replicates. Gaussian and binomial responses are supported.
#'
#' You can pass either:
#' \itemize{
#'   \item a standard model formula, e.g.,
#'     \code{y ~ X1 + X2 + X3 + I(X1^2) + (X1 + X2 + X3)^2}, or
#'   \item a \code{bigexp_spec} created by \code{bigexp_terms()}, in which
#'     case SVEMnet will build the design matrix deterministically (locked
#'     types, levels, and contrasts) and, if requested, swap the response
#'     to fit multiple independent responses over the same expansion.
#' }
#'
#' @param formula A formula specifying the model to be fitted, or a
#'   \code{bigexp_spec} created by \code{bigexp_terms()}.
#' @param data A data frame containing the variables in the model.
#' @param nBoot Integer. Number of bootstrap replicates (default \code{200}).
#'   Each replicate draws FRW weights and fits a glmnet path.
#' @param glmnet_alpha Numeric vector of elastic net mixing parameters
#'   \code{alpha} in \code{[0, 1]}. \code{alpha = 1} is lasso, \code{alpha = 0}
#'   is ridge. Defaults to \code{c(0.5, 1)}. When \code{relaxed = TRUE},
#'   \code{alpha = 0} is automatically dropped (ridge + relaxed is not used).
#' @param weight_scheme Character. Weighting scheme for train/validation
#'   copies. One of:
#'   \itemize{
#'     \item \code{"SVEM"} (default): Self-Validated Ensemble Model weights.
#'       For each replicate and row, a shared uniform draw \code{U_i ~ Unif(0, 1)}
#'       is converted to anti-correlated FRW weights
#'       \code{w_train_i = -log(U_i)} and \code{w_valid_i = -log(1 - U_i)}.
#'       Each weight vector is then rescaled to have mean 1 (sum \code{n}).
#'     \item \code{"FRW_plain"}: Fractional random-weight regression without a
#'       separate validation copy. A single FRW vector
#'       \code{w_i = -log(U_i)} is used for both training and validation and
#'       rescaled to have mean 1 (sum \code{n}). This reproduces the FRW
#'       bootstrap regression of Xu et al. (2020) and related work.
#'     \item \code{"Identity"}: Uses unit weights for both training and
#'       validation (no resampling). In combination with \code{nBoot = 1} this
#'       wraps a single \code{glmnet} fit and selects the penalty by the chosen
#'       information criterion, while still using SVEMnet's expansion and
#'       diagnostics.
#'   }
#' @param objective Character. One of \code{"auto"}, \code{"wAIC"},
#'   \code{"wBIC"}, or \code{"wSSE"}.
#'   \itemize{
#'   \item \code{"wAIC"}: Gaussian AIC-like criterion based on the weighted SSE.
#'   \item \code{"wBIC"}: Gaussian BIC-like criterion based on the weighted SSE.
#'   }
#'   See Details for the exact criteria in the Gaussian and binomial cases.
#' @param relaxed Logical or character. Default \code{"auto"}.
#'   If \code{TRUE}, use glmnet's relaxed elastic-net path and select both
#'   the penalty \code{lambda} and the relaxed refit parameter \code{gamma}
#'   on each bootstrap. If \code{FALSE}, fit the standard glmnet path without
#'   the relaxed step. If \code{"auto"} (default), SVEMnet uses
#'   \code{relaxed = TRUE} for \code{family = "gaussian"} and
#'   \code{relaxed = FALSE} for \code{family = "binomial"}.
#' @param response Optional character. When \code{formula} is a \code{bigexp_spec},
#'   this names the response column to use on the left-hand side. Defaults to
#'   the response stored in the spec.
#' @param unseen How to treat factor levels not seen in the original
#'   \code{bigexp_spec} when \code{formula} is a \code{bigexp_spec}. One of
#'   \code{"warn_na"} (default; convert unseen levels to \code{NA} with a
#'   warning) or \code{"error"} (stop with an error).
#' @param family Character. One of \code{"gaussian"} (default) or
#'   \code{"binomial"}. For Gaussian models SVEMnet uses the identity link; for
#'   binomial it uses the canonical logit link. The binomial response must be
#'   numeric 0/1, logical, or a factor with exactly two levels (the second
#'   level is treated as 1).
#' @param ... Additional arguments passed to \code{glmnet()}, such as
#'   \code{penalty.factor}, \code{lower.limits}, \code{upper.limits},
#'   \code{offset}, or \code{standardize.response}. Any user-supplied
#'   \code{weights} are ignored (SVEMnet supplies its own bootstrap weights).
#'   Any user-supplied \code{standardize} is ignored; SVEMnet always calls
#'   \code{glmnet} with \code{standardize = TRUE}.
#'
#' @return An object of class \code{"svem_model"} (and \code{"svem_binomial"}
#'   when \code{family = "binomial"}) with components:
#' \itemize{
#'   \item \code{parms}: Vector of ensemble-averaged coefficients, including
#'     the intercept.
#'   \item \code{parms_debiased}: Vector of coefficients after optional
#'     debiasing (see Details; Gaussian only).
#'   \item \code{debias_fit}: If debiasing was performed, the calibration
#'     model \code{lm(y ~ y_pred)}; otherwise \code{NULL}.
#'   \item \code{coef_matrix}: Matrix of per-bootstrap coefficients
#'     (rows = bootstraps, columns = intercept and predictors).
#'   \item \code{nBoot}: Number of bootstrap replicates actually used.
#'   \item \code{glmnet_alpha}: Vector of alpha values considered.
#'   \item \code{best_alphas}: Per-bootstrap alpha selected by the criterion.
#'   \item \code{best_lambdas}: Per-bootstrap lambda selected by the criterion.
#'   \item \code{best_relax_gammas}: Per-bootstrap relaxed gamma selected when
#'     \code{relaxed = TRUE}; \code{NA} otherwise.
#'   \item \code{weight_scheme}: The weighting scheme that was used.
#'   \item \code{relaxed}: Logical flag indicating whether relaxed paths were
#'     used.
#'   \item \code{relaxed_input}: The user-supplied value for \code{relaxed}
#'     (one of \code{TRUE}, \code{FALSE}, or \code{"auto"}). The resolved flag
#'     actually used is reported in \code{relaxed}.
#'   \item \code{dropped_alpha0_for_relaxed}: Logical; \code{TRUE} if
#'     \code{alpha = 0} was dropped because \code{relaxed = TRUE}.
#'   \item \code{objective_input}: The objective requested by the user.
#'   \item \code{objective_used}: The objective actually used after applying
#'     the "auto" rule (for example \code{"wAIC"} or \code{"wBIC"}).
#'   \item \code{objective}: Same as \code{objective_used} (for convenience).
#'   \item \code{auto_used}: Logical; \code{TRUE} if \code{objective = "auto"}.
#'   \item \code{auto_decision}: The objective selected by the auto rule
#'     (wAIC or wBIC) when \code{auto_used = TRUE}.
#'   \item \code{diagnostics}: List with summary information, including:
#'     \itemize{
#'       \item \code{k_summary}: Median and IQR of selected model size
#'         (number of nonzero coefficients including intercept).
#'       \item \code{fallback_rate}: Proportion of bootstraps that fell back
#'         to an intercept-only fit.
#'        \item \code{n_eff_summary}: Summary of raw Kish effective validation
#'         sizes \eqn{n_eff = (\sum_i w_i^{valid})^2 / \sum_i (w_i^{valid})^2}
#'         across bootstraps (before truncation to form \eqn{n_eff_adm}).
#'       \item \code{alpha_freq}: Relative frequency of selected alpha values
#'         (if any).
#'       \item \code{relax_gamma_freq}: Relative frequency of selected relaxed
#'         gamma values (if \code{relaxed = TRUE} and any were selected).
#'     }
#'   \item \code{actual_y}: Numeric response vector used in fitting (0/1 for
#'     binomial).
#'   \item \code{training_X}: Numeric model matrix without the intercept
#'     column used for training.
#'   \item \code{y_pred}: Fitted values from the ensemble on the training
#'     data. For Gaussian this is on the response scale; for binomial it is
#'     on the probability scale.
#'   \item \code{y_pred_debiased}: Debiased fitted values on the training
#'     data (Gaussian only); \code{NULL} otherwise.
#'   \item \code{nobs}: Number of observations used in fitting.
#'   \item \code{nparm}: Number of parameters in the full expansion
#'     (intercept plus predictors).
#'   \item \code{formula}: The formula used for fitting (possibly derived
#'     from a \code{bigexp_spec}).
#'   \item \code{terms}: \code{terms} object used for building the design
#'     matrix, with environment set to \code{baseenv()} for safety.
#'   \item \code{xlevels}: Factor levels recorded at training time.
#'   \item \code{contrasts}: Contrasts used for building the design matrix.
#'   \item \code{schema}: Compact description for safe prediction, including
#'     \code{feature_names}, \code{terms_str}, \code{xlevels},
#'     \code{contrasts}, \code{contrasts_options}, and a simple hash.
#'   \item \code{sampling_schema}: Schema used to generate random candidate
#'     tables, including predictor names, variable classes, numeric ranges,
#'     and factor levels.
#'   \item \code{used_bigexp_spec}: Logical flag indicating whether a
#'     \code{bigexp_spec} was used.
#'   \item \code{family}: The fitted family ("gaussian" or "binomial").
#' }
#'
#' @details
#' SVEMnet implements Self-Validated Ensemble Models using elastic net and
#' relaxed elastic net base learners from glmnet. Each bootstrap replicate
#' draws fractional random weights, builds a train and validation copy, fits
#' a path over lambda (and optionally over alpha and relaxed gamma), and
#' selects a path point by minimizing a validation-weighted criterion. Final
#' predictions are obtained by averaging replicate predictions on the chosen
#' scale.
#'
#' By default, \code{relaxed = "auto"} resolves to \code{TRUE} for Gaussian
#' fits and \code{FALSE} for binomial fits.
#'
#' The function is typically used in small-n design-of-experiments (DOE)
#' workflows where classical train/validation splits and cross-validation
#' can be unstable. A common pattern is: (1) build a deterministic expansion
#' with \code{bigexp_terms()}, (2) fit SVEM models via \code{SVEMnet()},
#' (3) perform whole-model significance testing, and (4) call
#' \code{svem_score_random()} for constrained multi-response optimization.
#'
#' Weighting schemes:
#' \itemize{
#'   \item With \code{weight_scheme = "SVEM"}, SVEMnet uses a pair of
#'     anti-correlated FRW vectors for train and validation. All rows appear
#'     in every replicate, but train and validation contributions are
#'     separated through the shared uniform draws.
#'   \item With \code{weight_scheme = "FRW_plain"}, a single FRW vector is
#'     used for both train and validation, which reproduces FRW regression
#'     without a self-validation split. This is mainly provided for method
#'     comparison and teaching.
#'   \item With \code{weight_scheme = "Identity"}, both train and validation
#'     weights are 1. Setting \code{nBoot = 1} in this mode yields a single
#'     glmnet fit whose penalty is chosen by the selected information
#'     criterion, without any bootstrap variation.
#' }
#'
#' Selection criteria (Gaussian):
#' For \code{family = "gaussian"}, the validation loss is a weighted sum of
#' squared errors on the validation copy. Let
#' \eqn{\mathrm{SSE}_w = \sum_i w_i^{valid} r_i^2} denote the weighted SSE.
#' The criteria are:
#' \itemize{
#'   \item \code{"wSSE"}: loss-only selector that minimizes the weighted SSE
#'     \eqn{\mathrm{SSE}_w},
#'   \item \code{"wAIC"}: Gaussian AIC analog
#'     \eqn{C(\lambda) = n \log\{\mathrm{SSE}_w(\lambda)/n\} + 2 k},
#'   \item \code{"wBIC"}: Gaussian BIC analog
#'     \eqn{C(\lambda) = n \log\{\mathrm{SSE}_w(\lambda)/n\} +
#'       \log(n_eff_adm)\, k}.
#' }
#' The FRW validation weights are rescaled to have mean one, so that their sum
#' is always \eqn{\sum_i w_i^{valid} = n}. The AIC/BIC analogs therefore use
#' \eqn{n \log(\mathrm{SSE}_w/n)} as the Gaussian loss term, while the
#' \code{"wSSE"} selector uses \eqn{\mathrm{SSE}_w} directly.
#'
#' The effective validation size is computed from the FRW weights using Kish's
#' effective sample size \eqn{n_eff = (\sum_i w_i^{valid})^2 /
#'   \sum_i (w_i^{valid})^2} and then truncated to lie between 2 and \code{n}
#' to form \eqn{n_eff_adm}. The AIC-style selector uses a \eqn{2 k} penalty;
#' the BIC-style selector uses a \eqn{\log(n_eff_adm) k} penalty, so that the
#' loss term is scaled by total validation weight while the complexity penalty
#' reflects the effective amount of information under unequal weights. For
#' \code{"wAIC"} and \code{"wBIC"}, path points with more than
#' \eqn{n_eff_adm} non-intercept coefficients are treated as inadmissible when
#' evaluating the criterion.
#'
#' Because the FRW validation weights are random rather than fixed design
#' weights, these information-criterion scores are used heuristically for
#' relative model comparison within each FRW replicate, rather than as exact
#' AIC/BIC values.

#'
#' For diagnostics, SVEMnet reports the raw Kish effective sizes across
#' bootstraps (see \code{diagnostics$n_eff_summary}), while \eqn{n_eff_adm}
#' is used internally in the penalty and model-size guardrail. Near-interpolating
#' path points are screened out via a simple model size guardrail before
#' minimization. When \code{objective = "auto"}, SVEMnet uses \code{"wAIC"}.
#'
#' This structure (pseudo-likelihood using total weight and BIC penalty using a
#' Kish-type effective sample size) parallels survey-weighted information
#' criteria as in Lumley and Scott (2015) and Kish (1965).
#'
#' Selection criteria (binomial):
#' For \code{family = "binomial"}, the validation loss is the weighted
#' negative log-likelihood on the FRW validation copy (equivalently,
#' proportional to the binomial deviance up to a constant factor). Let
#' \eqn{\mathrm{NLL}} denote the weighted negative log-likelihood. The same
#' labels are used:
#' \itemize{
#'   \item \code{"wSSE"}: loss-only selector based on \eqn{\mathrm{NLL}}
#'     (the name is retained for backward compatibility),
#'   \item \code{"wAIC"}: deviance-style criterion
#'     \eqn{C(\lambda) = 2\,\mathrm{NLL}(\lambda) + 2 k},
#'   \item \code{"wBIC"}: deviance-style criterion
#'     \eqn{C(\lambda) = 2\,\mathrm{NLL}(\lambda) + \log(n_eff_adm)\, k}.
#' }
#' The effective validation size \eqn{n_eff_adm} and the model size guardrail
#' are handled as in the Gaussian case: for \code{"wAIC"} and \code{"wBIC"}
#' we compute a Kish effective size from the FRW validation weights, truncate
#' it to lie between 2 and \code{n}, and require the number of nonzero
#' coefficients (excluding the intercept) to be less than this effective size
#' when evaluating the criterion.
#'
#' Auto rule:
#' When \code{objective = "auto"}, SVEMnet selects the criterion by family:
#' \itemize{
#'   \item \code{family = "gaussian"} -> \code{"wAIC"}
#'   \item \code{family = "binomial"} -> \code{"wBIC"}
#' }
#'
#' Relaxed elastic net:
#' When \code{relaxed = TRUE}, SVEMnet calls glmnet with
#' \code{relax = TRUE} and traverses a small grid of relaxed refit values
#' (gamma). For each alpha and gamma, SVEMnet evaluates all lambda path
#' points on the validation copy and records the combination that minimizes
#' the selected criterion. Model size is always defined as the number of
#' nonzero coefficients including the intercept, so standard and relaxed
#' paths are scored on the same scale.
#'
#' Gaussian debiasing:
#' For Gaussian models, SVEMnet optionally performs a simple linear
#' calibration of ensemble predictions on the training data. When there is
#' sufficient variation in the fitted values and \code{nBoot} is at least
#' 10, the function fits \code{lm(y ~ y_pred)} and uses the coefficients to
#' construct debiased coefficients and debiased fitted values. Binomial
#' fits do not use debiasing; predictions are ensembled on the probability
#' or link scale directly.
#'
#' Implementation notes:
#' \itemize{
#'   \item Predictors are always standardized internally via
#'     \code{glmnet(..., standardize = TRUE)}.
#'   \item The terms object is stored with its environment set to
#'     \code{baseenv()} so that prediction does not accidentally capture
#'     objects from the calling environment.
#'   \item A compact schema (feature names, factor levels, contrasts, and
#'     a simple hash) is stored to allow \code{predict()} and companion
#'     functions to rebuild model matrices deterministically, even when the
#'     original data frame is not available.
#'   \item A separate sampling schema stores raw predictor ranges and
#'     factor levels for use in random candidate generation for optimization.
#' }
#'
#' @template ref-svem
#'
#' @importFrom stats runif lm predict coef var model.frame model.matrix
#'   model.response delete.response IQR median
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
#' y <- 1 + 2 * X1 - 1.5 * X2 + 0.5 * X3 + 1.2 * (X1 * X2) +
#'   0.8 * (X1^2) + eps
#' dat <- data.frame(y, X1, X2, X3)
#'
#' # Minimal hand-written expansion
#' mod_relax <- SVEMnet(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2),
#'   data          = dat,
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 75,
#'   objective     = "auto",
#'   weight_scheme = "SVEM",
#'   relaxed       = FALSE
#' )
#'
#' pred_in_raw <- predict(mod_relax, dat, debias = FALSE)
#' pred_in_db  <- predict(mod_relax, dat, debias = TRUE)
#'
#' \donttest{
#' # ---------------------------------------------------------------------------
#' # Big expansion (full factorial + polynomial surface + partial-cubic crosses)
#' # Build once, reuse for one or more responses
#' # ---------------------------------------------------------------------------
#' spec <- bigexp_terms(
#'   y ~ X1 + X2 + X3,
#'   data             = dat,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_3way  = FALSE
#' )
#'
#' # Fit using the spec (auto-prepares data)
#' fit_y <- SVEMnet(
#'   spec, dat,
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 50,
#'   objective     = "auto",
#'   weight_scheme = "SVEM"
#' )
#'
#' # A second, independent response over the same expansion
#' set.seed(99)
#' dat$y2 <- 0.5 + 1.4 * X1 - 0.6 * X2 + 0.2 * X3 + rnorm(n, 0, 0.4)
#' fit_y2 <- SVEMnet(
#'   spec, dat, response = "y2",
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 50,
#'   objective     = "auto",
#'   weight_scheme = "SVEM"
#' )
#'
#' svem_nonzero(fit_y2)
#'
#' p1 <- predict(fit_y,  dat)
#' p2 <- predict(fit_y2, dat, debias = TRUE)
#'
#' # Show that a new batch expands identically under the same spec
#' newdat <- data.frame(
#'   y  = y,
#'   X1 = X1 + rnorm(n, 0, 0.05),
#'   X2 = X2 + rnorm(n, 0, 0.05),
#'   X3 = X3 + rnorm(n, 0, 0.05)
#' )
#' prep_new <- bigexp_prepare(spec, newdat)
#' stopifnot(identical(
#'   colnames(model.matrix(spec$formula, bigexp_prepare(spec, dat)$data)),
#'   colnames(model.matrix(spec$formula, prep_new$data))
#' ))
#' preds_new <- predict(fit_y, prep_new$data)
#'
#' ## Binomial example
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
#'   nBoot        = 50,
#'   glmnet_alpha = c(1, 0.5),
#'   family       = "binomial"
#' )
#'
#' ## Probabilities, link, and classes
#' p_resp <- predict(fit_b, db, type = "response")
#' p_link <- predict(fit_b, db, type = "link")
#' y_hat  <- predict(fit_b, db, type = "class")  # 0/1 labels
#'
#' ## Mean aggregation with uncertainty on probability scale
#' out_b <- predict(
#'   fit_b, db,
#'   type     = "response",
#'   se.fit   = TRUE,
#'   interval = TRUE,
#'   level    = 0.9
#' )
#' str(out_b)
#'
#' #' ## Example with blocking (requires SVEMnet to store sampling_schema$blocking)
#' set.seed(2)
#' df_block <- data.frame(
#'   y1         = rnorm(40),
#'   y2         = rnorm(40),
#'   X1         = runif(40),
#'   X2         = runif(40),
#'   Operator   = factor(sample(paste0("Op", 1:3), 40, TRUE)),
#'   AmbientTmp = rnorm(40, mean = 22, sd = 2)
#' )
#'
#' spec_block <- bigexp_terms(
#'   y1 ~ X1 + X2,
#'   data             = df_block,
#'   factorial_order  = 2,
#'   polynomial_order = 2,
#'   blocking         = c("Operator", "AmbientTmp")
#' )
#'
#' fit_b1 <- SVEMnet(spec_block, df_block, response = "y1", nBoot = 30)
#' fit_b2 <- SVEMnet(spec_block, df_block, response = "y2", nBoot = 30)
#'
#' tab_block <- svem_random_table_multi(list(fit_b1, fit_b2), n = 500)
#' }
#' @export
SVEMnet <- function(formula, data,
                        nBoot = 200,
                        glmnet_alpha = c( 0.5, 1),
                        weight_scheme = c("SVEM", "FRW_plain", "Identity"),
                        objective = c("auto", "wAIC", "wBIC", "wSSE"),
                        relaxed = "auto",
                        response = NULL,
                        unseen = c("warn_na","error"),
                        family = c("gaussian", "binomial"),
                        ...) {

  ## ----------------- family handling -----------------
  if (inherits(family, "family")) {
    stop(
      "SVEMnet currently accepts family only as a character string: ",
      "family = 'gaussian' or family = 'binomial'."
    )
  }

  family   <- match.arg(family)
  fam_name <- family

  if (!fam_name %in% c("gaussian", "binomial")) {
    stop("SVEMnet currently supports family = 'gaussian' or 'binomial' only.")
  }


  ## ---- argument checks ----
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)
  unseen        <- match.arg(unseen)

  ## ---- relaxed handling (supports logical or "auto") ----
  relaxed_input <- relaxed
  relax_flag <- NA

  if (is.character(relaxed)) {
    relaxed <- match.arg(relaxed, c("auto"))
    # family has already been parsed into fam_name above
    relax_flag <- if (relaxed == "auto") (fam_name == "gaussian") else NA
  } else if (is.logical(relaxed) && length(relaxed) == 1L && !is.na(relaxed)) {
    relax_flag <- relaxed
  } else {
    stop("'relaxed' must be TRUE, FALSE, or the character string \"auto\" (default).")
  }


  if (!is.numeric(nBoot) || length(nBoot) != 1L || !is.finite(nBoot) || nBoot < 1) {
    stop("nBoot must be >= 1.")
  }
  nBoot <- as.integer(nBoot)

  if (!is.numeric(glmnet_alpha) || any(!is.finite(glmnet_alpha))) {
    stop("'glmnet_alpha' must be numeric and finite.")
  }
  if (any(glmnet_alpha < 0 | glmnet_alpha > 1)) {
    stop("'glmnet_alpha' values must be in [0, 1].")
  }
  glmnet_alpha <- unique(as.numeric(glmnet_alpha))
  if (!length(glmnet_alpha)) glmnet_alpha <- 1



  ## ---- model frame / matrix construction ----
  data <- as.data.frame(data)
  contrasts_opts_used <- NULL

  # Detect bigexp_spec either passed directly or attached to the formula
  using_spec     <- inherits(formula, "bigexp_spec")
  spec           <- NULL
  spec_from_attr <- FALSE

  if (using_spec) {
    # Case 1: SVEMnet(spec, data, ...)
    spec <- formula
  } else {
    # Case 2: SVEMnet(form_with_attr, data, ...)
    spec_attr <- attr(formula, "bigexp_spec", exact = TRUE)
    if (!is.null(spec_attr) && inherits(spec_attr, "bigexp_spec")) {
      using_spec     <- TRUE
      spec           <- spec_attr
      spec_from_attr <- TRUE
    }
  }

  if (using_spec) {
    # Decide which formula to use:
    # - If 'response' is supplied, always swap LHS to that name.
    # - Else, if the spec came via attribute, respect the formula's own LHS.
    # - Else (direct spec), use spec$formula as stored.
    if (!is.null(response)) {
      rhs_txt <- paste(deparse(spec$formula[[3L]]), collapse = " ")
      f_use   <- stats::as.formula(paste(response, "~", rhs_txt))
    } else if (spec_from_attr) {
      f_use <- formula
    } else {
      f_use <- spec$formula
    }

    prep <- bigexp_prepare(spec, data, unseen = unseen)

    spec_settings <- spec$settings
    if (is.null(spec_settings)) spec_settings <- list()

    # Blocking variables from the bigexp_spec (may be NULL / empty)
    blocking <- spec_settings$blocking
    if (is.null(blocking)) blocking <- character(0L)

    # Contrast options + per-factor contrasts from the spec
    if (!is.null(spec_settings$contrasts_options)) {
      spec_contrasts_opts <- spec_settings$contrasts_options
    } else {
      spec_contrasts_opts <- getOption("contrasts")
    }

    if (!is.null(spec_settings$contrasts)) {
      spec_contrasts <- spec_settings$contrasts
    } else {
      spec_contrasts <- NULL
    }

    contrasts_opts_used <- spec_contrasts_opts
    old_opts <- options(contrasts = spec_contrasts_opts)
    on.exit(options(old_opts), add = TRUE)

    mf <- stats::model.frame(f_use, prep$data, na.action = stats::na.omit)
    if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")
    X  <- stats::model.matrix(f_use, mf, contrasts.arg = spec_contrasts)

  } else {
    # Plain formula path, no bigexp_spec
    f_use <- formula
    mf    <- stats::model.frame(f_use, data, na.action = stats::na.omit)
    if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")
    contrasts_opts_used <- getOption("contrasts")
    X       <- stats::model.matrix(f_use, mf)
    blocking <- character(0L)   # no blocking when not using bigexp_spec
  }



  y <- stats::model.response(mf)

  ## ---- response handling by family ----
  .coerce_binomial_01 <- function(y) {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop("Binomial SVEMnet requires a factor with exactly 2 levels or 0/1 numeric y.")
      }
      return(as.integer(y == levels(y)[2L]))
    } else if (is.logical(y)) {
      return(as.integer(y))
    } else if (is.numeric(y)) {
      uy <- sort(unique(y))
      if (!all(uy %in% c(0,1))) {
        stop("Binomial SVEMnet requires numeric y coded as 0/1.")
      }
      return(as.integer(y))
    } else {
      stop("Unsupported y type for binomial; use 0/1 numeric, logical, or a 2-level factor.")
    }
  }

  if (fam_name == "gaussian") {
    if (!is.numeric(y)) stop("Gaussian SVEMnet requires numeric y.")
    y_vec <- as.numeric(y)
  } else if (fam_name == "binomial") {
    y_vec <- .coerce_binomial_01(y)
  }

  ## drop intercept column (glmnet adds its own)
  int_idx <- which(colnames(X) %in% c("(Intercept)", "Intercept"))
  if (length(int_idx)) X <- X[, -int_idx, drop = FALSE]
  if (ncol(X) == 0L) stop("SVEMnet requires at least one predictor.")

  if (any(!is.finite(y_vec)) || any(!is.finite(X))) {
    stop("Non-finite values in response/predictors after NA handling.")
  }
  storage.mode(X) <- "double"

  n <- nrow(X); p <- ncol(X); nobs <- n; nparm <- p + 1L

  ## capture training xlevels and contrasts for prediction
  terms_full  <- attr(mf, "terms")
  terms_clean <- terms_full; environment(terms_clean) <- baseenv()

  predictor_vars <- base::all.vars(stats::delete.response(terms_full))

  ## ---- guard: response must not also be a predictor ----
  resp_name <- tryCatch(
    as.character(f_use[[2L]]),
    error = function(e) NA_character_
  )
  if (!is.na(resp_name) && resp_name %in% predictor_vars) {
    stop(
      "SVEMnet does not allow a predictor with the same name as the response ('",
      resp_name, "'). Please rename your variables or adjust the formula."
    )
  }

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
        xlevels[[v]] <- lev
        factor_levels[[v]] <- lev
      }
    }
  }
  contrasts_used <- attr(X, "contrasts")
  # Helper for categorical mode (used for blocking factors)
  .mode_level <- function(x) {
    if (is.null(x)) return(NA_character_)
    x_chr <- as.character(x)
    x_chr <- x_chr[!is.na(x_chr)]
    if (!length(x_chr)) return(NA_character_)
    tab <- table(x_chr)
    if (!length(tab)) return(NA_character_)
    maxfreq <- max(tab)
    modes <- names(tab)[tab == maxfreq]
    modes[1L]  # break ties arbitrarily but deterministically
  }

  # Normalize blocking to be only predictors actually present
  blocking <- intersect(blocking, predictor_vars)

  # Modes for blocking categorical variables, if any
  block_cat_modes <- NULL
  if (length(blocking)) {
    # "Categorical" here means we stored factor_levels
    block_cat <- intersect(blocking, names(factor_levels))
    if (length(block_cat)) {
      block_cat_modes <- vector("list", length(block_cat))
      names(block_cat_modes) <- block_cat
      for (v in block_cat) {
        if (v %in% colnames(mf)) {
          mode_v <- .mode_level(mf[[v]])
          if (!is.na(mode_v)) {
            block_cat_modes[[v]] <- mode_v
          }
        }
      }
      # Drop if everything came back NA / empty
      if (!length(block_cat_modes) ||
          all(vapply(block_cat_modes, is.null, logical(1)))) {
        block_cat_modes <- NULL
      }
    }
  }


  num_vars <- predictor_vars[vapply(predictor_vars, function(v) {
    v %in% colnames(mf) && is.numeric(mf[[v]])
  }, logical(1))]
  if (length(num_vars)) {
    rng_mat <- vapply(num_vars, function(v) {
      r <- range(mf[[v]], na.rm = TRUE)
      if (!all(is.finite(r)) || r[1] == r[2]) {
        r <- c(min(mf[[v]], na.rm = TRUE), max(mf[[v]], na.rm = TRUE))
      }
      r
    }, numeric(2))
    rownames(rng_mat) <- c("min","max")
    num_ranges <- as.matrix(rng_mat)
  } else {
    num_ranges <- matrix(numeric(0), nrow = 2, ncol = 0,
                         dimnames = list(c("min","max"), NULL))
  }

  ## ---- objective selection (auto) ----
  auto_used     <- identical(objective, "auto")
  auto_decision <- NA_character_

  objective_used <- if (auto_used) {
    if (fam_name == "gaussian") "wAIC" else "wBIC"
  } else {
    objective
  }
  auto_decision <- if (auto_used) objective_used else NA_character_


  ## ---- drop ridge when relaxed ----
  dropped_alpha0_for_relaxed <- FALSE
  if (isTRUE(relax_flag) && any(glmnet_alpha == 0)) {
    warning("Dropping alpha = 0 (ridge) for relaxed fits; ridge + relaxed is not supported here.")
    glmnet_alpha <- glmnet_alpha[glmnet_alpha != 0]
    if (!length(glmnet_alpha)) glmnet_alpha <- 1
    dropped_alpha0_for_relaxed <- TRUE
  }

  ## ---- containers ----
  coef_matrix  <- matrix(NA_real_, nrow = nBoot, ncol = p + 1L)
  colnames(coef_matrix) <- c("(Intercept)", colnames(X))
  best_alphas       <- rep(NA_real_, nBoot)
  best_lambdas      <- rep(NA_real_, nBoot)
  best_relax_gammas <- rep(NA_real_, nBoot)
  k_sel_vec         <- rep(NA_integer_, nBoot)
  fallbacks         <- integer(nBoot)
  n_eff_keep        <- numeric(nBoot)

  ## ---- capture and sanitize user '...' ----
  dots <- list(...)
  if (length(dots)) {
    if ("weights" %in% names(dots)) {
      warning("Ignoring user 'weights'; SVEM uses bootstrap weights.")
      dots$weights <- NULL
    }
    if ("standardize" %in% names(dots)) {
      warning("Ignoring user 'standardize'; SVEMnet always standardizes.")
      dots$standardize <- NULL
    }
    protected <- c("x","y","alpha","intercept","relax","standardize",
                   "nlambda","maxit","lambda.min.ratio","lambda","family")
    dots <- dots[setdiff(names(dots), protected)]
  }

  .support_size_one <- function(beta_col,
                                base_tol        = 1e-7,
                                count_intercept = TRUE) {
    if (!is.numeric(beta_col) || length(beta_col) < 1L) {
      return(NA_integer_)
    }

    if (length(beta_col) == 1L) {
      return(if (count_intercept) 1L else 0L)
    }

    slope <- beta_col[-1L]
    slope <- slope[is.finite(slope)]

    if (!length(slope)) {
      return(if (count_intercept) 1L else 0L)
    }

    max_slope <- max(abs(slope))
    tol_rel   <- base_tol * max(1, max_slope)
    tol_abs   <- base_tol
    tol_j     <- max(tol_rel, tol_abs)

    k_slope <- sum(abs(slope) > tol_j)

    (if (count_intercept) 1L else 0L) + k_slope
  }


  ## ---- bootstrap loop ----
  for (i in seq_len(nBoot)) {
    eps <- .Machine$double.eps
    if (weight_scheme == "SVEM") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- -log1p(-U)
    } else if (weight_scheme == "FRW_plain") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- w_train
    } else {
      w_train <- rep(1, n); w_valid <- rep(1, n)
    }
    w_train <- w_train * (n / sum(w_train))
    w_valid <- w_valid * (n / sum(w_valid))

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

    relax_gamma_grid <- if (isTRUE(relax_flag)) c( 0.2, 0.6,  1) else 1

    for (alpha in glmnet_alpha) {
      fit <- tryCatch({
        withCallingHandlers({
          do.call(glmnet::glmnet,
                  c(list(x = X, y = y_vec,
                         alpha = alpha,
                         weights = w_train,
                         intercept = TRUE,
                         standardize = TRUE,
                         nlambda = 500,
                         maxit = 1e6,
                         relax = isTRUE(relax_flag),
                         family = fam_name),
                    dots))
        }, warning = function(w) base::invokeRestart("muffleWarning"))
      }, error = function(e) NULL)
      if (is.null(fit) || !length(fit$lambda)) next

      for (relax_g in relax_gamma_grid) {
        coef_path <- tryCatch({
          if (isTRUE(relax_flag)) {
            as.matrix(stats::coef(fit, s = fit$lambda, gamma = relax_g))
          } else {
            as.matrix(stats::coef(fit, s = fit$lambda))
          }
        }, error = function(e) NULL)
        if (is.null(coef_path) || nrow(coef_path) != (p + 1L)) next
        L <- ncol(coef_path); if (L == 0L) next

        pred_link <- tryCatch({
          XB <- X %*% coef_path[-1L, , drop = FALSE]
          sweep(XB, 2, coef_path[1L, ], FUN = "+")
        }, error = function(e) NULL)
        if (is.null(pred_link) || nrow(pred_link) != n) next

        ## family-specific validation loss
        if (fam_name == "gaussian") {
          res <- pred_link - y_vec
          val_core <- as.vector(crossprod(w_valid, res^2))  # weighted SSE
          val_core[!is.finite(val_core)] <- Inf
          adj_val_core <- pmax(val_core, .Machine$double.eps)
        } else {  # binomial: negative log-likelihood
          mu <- 1 / (1 + exp(-pred_link))
          mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
          yy <- matrix(y_vec, nrow = n, ncol = L)
          loglik <- colSums(w_valid * (yy * log(mu) + (1 - yy) * log(1 - mu)))
          val_core <- -loglik   # NLL
          val_core[!is.finite(val_core)] <- Inf
          adj_val_core <- pmax(val_core, .Machine$double.eps)
        }

        k_raw <- integer(L)
        for (jj in seq_len(L)) {
          k_raw[jj] <- .support_size_one(coef_path[, jj], 1e-7, TRUE)
        }
        if (length(k_raw) != length(adj_val_core)) next

        k_slope <- pmax(0L, k_raw - 1L)
        k_eff   <- 1L + k_slope
        logN_pen <- log(n_eff_adm)

        if (fam_name == "gaussian") {
          n_like <- sumw
          mse_w  <- adj_val_core / n_like
          metric <- switch(objective_used,
                           "wSSE" = adj_val_core,
                           "wAIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- n_like * log(mse_w[mask]) + 2 * k_eff[mask]
                             out
                           },
                           "wBIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- n_like * log(mse_w[mask]) +
                               logN_pen * k_eff[mask]
                             out
                           },
                           stop("Unknown objective: ", objective_used))
        } else {
          ## binomial: objective based on NLL
          metric <- switch(objective_used,
                           "wSSE" = adj_val_core,
                           "wAIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- 2 * adj_val_core[mask] + 2 * k_eff[mask]
                             out
                           },
                           "wBIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- 2 * adj_val_core[mask] +
                               logN_pen * k_eff[mask]
                             out
                           },
                           stop("Unknown objective: ", objective_used))
        }

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
      mu_w          <- sum(w_train * y_vec) / sum(w_train)
      if (fam_name == "gaussian") {
        best_beta_hat <- c(mu_w, rep(0, p))
      } else {
        mu_w <- pmin(pmax(mu_w, 1e-6), 1 - 1e-6)
        best_beta_hat <- c(qlogis(mu_w), rep(0, p))
      }
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

  ## ---- finalize ----
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
  lin_pred <- as.vector(X %*% avg_coefficients[-1L] + avg_coefficients[1L])

  if (fam_name == "gaussian") {
    y_pred <- lin_pred
  } else {
    y_pred <- 1 / (1 + exp(-lin_pred))  # in-model scale = probability
  }

  ## debias: only for gaussian (leave binomial as-is for now)
  debias_fit <- NULL
  y_pred_debiased <- NULL
  parms_debiased  <- avg_coefficients

  if (fam_name == "gaussian" && nBoot >= 10 && stats::var(y_pred) > 0) {
    debias_fit <- stats::lm(y_vec ~ y_pred)
    y_pred_debiased <- stats::predict(debias_fit)
    parms_debiased <- avg_coefficients
    ab <- try(stats::coef(debias_fit), silent = TRUE)
    if (!inherits(ab, "try-error") &&
        length(ab) >= 2 &&
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
        if (length(parms_debiased) > 1L) parms_debiased[-1] <- b * parms_debiased[-1]
      }
    }
  }

  diagnostics <- list(
    k_summary = c(k_median = stats::median(k_sel_vec),
                  k_iqr    = stats::IQR(k_sel_vec)),
    fallback_rate = mean(fallbacks),
    n_eff_summary = summary(n_eff_keep),
    relax_gamma_freq =
      if (isTRUE(relax_flag) && any(is.finite(best_relax_gammas))) {
        prop.table(table(round(best_relax_gammas, 2)))
      } else NULL
  )
  tf <- table(best_alphas[is.finite(best_alphas)])
  diagnostics$alpha_freq <- if (length(tf)) {
    out <- as.numeric(tf) / sum(tf)
    names(out) <- names(tf)
    out
  } else numeric()

  feature_names <- colnames(X)
  terms_str <- tryCatch(
    paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
    error = function(e) NA_character_
  )
  safe_hash <- function(s) {
    if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
    bytes <- charToRaw(paste0(s, collapse = ""))
    sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
  }
  schema <- list(
    feature_names     = feature_names,
    terms_str         = terms_str,
    xlevels           = xlevels,
    contrasts         = contrasts_used,
    contrasts_options = contrasts_opts_used,
    terms_hash        = safe_hash(terms_str)
  )

  sampling_schema <- list(
    predictor_vars  = predictor_vars,
    var_classes     = var_classes,
    num_ranges      = num_ranges,
    factor_levels   = factor_levels,
    blocking        = blocking,        # character(0) when no blocking
    block_cat_modes = block_cat_modes  # NULL when none / not available
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
    best_relax_gammas = if (isTRUE(relax_flag)) best_relax_gammas else
      rep(NA_real_, length(best_alphas)),
    weight_scheme     = weight_scheme,
    relaxed           = isTRUE(relax_flag),
    relaxed_input     = relaxed_input,
    dropped_alpha0_for_relaxed = dropped_alpha0_for_relaxed,
    objective_input   = objective,
    objective_used    = objective_used,
    objective         = objective_used,
    auto_used         = auto_used,
    auto_decision     = if (auto_used) auto_decision else NA_character_,
    diagnostics       = diagnostics,
    actual_y          = y_vec,
    training_X        = X,
    y_pred            = y_pred,
    y_pred_debiased   = y_pred_debiased,
    nobs              = nobs,
    nparm             = nparm,
    formula           = f_use,
    terms             = terms_clean,
    xlevels           = xlevels,
    contrasts         = contrasts_used,
    schema            = schema,
    sampling_schema   = sampling_schema,
    used_bigexp_spec  = using_spec,
    family            = fam_name
  )

  cls <- c("svem_model", "SVEMnet")
  if (fam_name == "binomial") {
    cls <- c("svem_binomial", cls)
  }
  class(result) <- cls
  result
}
