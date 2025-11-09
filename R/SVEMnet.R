#' Fit an SVEMnet Model (with optional relaxed elastic net)
#'
#' Wrapper for glmnet (Friedman et al. 2010) to fit an ensemble of Elastic Net
#' models using the Self-Validated Ensemble Model method (SVEM; Lemkus et al. 2021),
#' with an option to use glmnet's built-in relaxed elastic net (Meinshausen 2007).
#' Supports searching over multiple alpha values in the Elastic Net penalty.
#'
#' You can pass either:
#' - a standard model formula, e.g. y ~ X1 + X2 + X3 + I(X1^2) + (X1 + X2 + X3)^2
#' - a bigexp_spec created by bigexp_terms(), in which case SVEMnet will prepare
#'   the data deterministically (locked types/levels) and, if requested, swap
#'   the response to fit multiple independent responses over the same expansion.
#'
#' @param formula A formula specifying the model to be fitted, OR a bigexp_spec
#'   created by \code{bigexp_terms()}.
#' @param data A data frame containing the variables in the model.
#' @param nBoot Number of bootstrap iterations (default 200).
#' @param glmnet_alpha Elastic Net mixing parameter(s). May be a vector with
#'   entries in the range between 0 and 1, inclusive, where alpha = 1 is Lasso
#'   and alpha = 0 is Ridge. Defaults to \code{c(0.5, 1)}.
#' @param weight_scheme Character; weighting scheme used to generate bootstrap
#'   training and validation weights (default \code{"SVEM"}).
#'   One of:
#'   \itemize{
#'     \item \code{"SVEM"}: Self-Validated Ensemble Model weights. For each
#'       bootstrap, independent \eqn{U_i \sim \mathrm{Uniform}(0,1)} are drawn
#'       and converted to anti-correlated FRW copies
#'       \eqn{w_i^{\mathrm{train}} = -\log U_i} and
#'       \eqn{w_i^{\mathrm{valid}} = -\log(1-U_i)}, each rescaled to have mean 1.
#'       This is the default and implements the SVEM scheme of Lemkus et al.
#'     \item \code{"FRW_plain"}: Fractional–random–weight (FRW) regression
#'       without self-validation. A single FRW weight vector
#'       \eqn{w_i = -\log U_i} (rescaled to mean 1) is used for both training
#'       and validation. This reproduces the fractional-random-weight bootstrap
#'       regression of Xu et al. (2020) and related work, with one weighted fit
#'       and no self-validation split.
#'     \item \code{"Identity"}: Uses unit weights for both training and
#'       validation (no resampling). This is primarily useful with
#'       \code{nBoot = 1} when you want a single glmnet fit whose penalty
#'       is chosen via the selected information criterion (wAIC/wBIC/wSSE),
#'       while still using SVEMnet’s formula expansion and diagnostics.
#'   }
#' @param objective Objective used to pick lambda on each bootstrap path
#'   (default \code{"auto"}). One of \code{"auto"}, \code{"wAIC"}, \code{"wBIC"}, or \code{"wSSE"}.
#' @param auto_ratio_cutoff Single cutoff for the automatic rule when
#'   \code{objective = "auto"} (default 1.3). Let \code{r = n_X / p_X}, where \code{n_X} is the
#'   number of training rows and \code{p_X} is the number of predictors in the model
#'   matrix after dropping the intercept column. If \code{r >= auto_ratio_cutoff},
#'   SVEMnet uses wAIC; otherwise it uses wBIC.
#' @param relaxed Logical, TRUE or FALSE (default TRUE). When TRUE (for
#'   \code{family = "gaussian"}), use glmnet's relaxed elastic net path and
#'   select both lambda and relaxed gamma on each bootstrap. When FALSE,
#'   fit the standard glmnet path. Note: if \code{relaxed = TRUE} and
#'   \code{glmnet_alpha} includes 0 (ridge), alpha = 0 is dropped. For
#'   \code{family = "binomial"}, relaxed fits are currently disabled for
#'   stability: SVEMnet behaves as if \code{relaxed = FALSE}, and a warning
#'   is issued if \code{relaxed = TRUE} is requested.
#' @param family Character scalar specifying the model family.
#'   One of \code{"gaussian"} (default) or \code{"binomial"}.
#'   SVEMnet currently assumes the canonical identity link for Gaussian
#'   and the canonical logit link for binomial. For \code{"binomial"},
#'   the response must be numeric 0/1, logical, or a factor with exactly
#'   two levels (the second level is treated as 1).
#' @param response Optional character. When \code{formula} is a \code{bigexp_spec}, this names
#'   the response column to use on the LHS; defaults to the response stored in the spec.
#' @param unseen How to treat unseen factor levels when \code{formula} is a \code{bigexp_spec}:
#'   \code{"warn_na"} (default; convert to NA with a warning) or \code{"error"} (stop).
#' @param ... Additional args passed to \code{glmnet()} (e.g., \code{penalty.factor},
#'   \code{lower.limits}, \code{upper.limits}, \code{offset}, \code{standardize.response}, etc.).
#'   Any user-supplied \code{weights} are ignored so SVEM can supply its own bootstrap weights.
#'   Any user-supplied \code{standardize} is ignored; SVEMnet always uses \code{standardize = TRUE}.
#'
#' @return An object of class \code{svem_model} with elements:
#' \itemize{
#'   \item \code{parms}: averaged coefficients (including intercept).
#'   \item \code{parms_debiased}: averaged coefficients adjusted by the calibration fit.
#'   \item \code{debias_fit}: \code{lm(y ~ y_pred)} calibration model used for debiasing (or \code{NULL}).
#'   \item \code{coef_matrix}: per-bootstrap coefficient matrix.
#'   \item \code{nBoot}, \code{glmnet_alpha}, \code{best_alphas}, \code{best_lambdas}, \code{weight_scheme}, \code{relaxed}.
#'   \item \code{best_relax_gammas}: per-bootstrap relaxed gamma chosen (NA if \code{relaxed = FALSE}).
#'   \item \code{objective_input}, \code{objective_used}, \code{objective} (same as \code{objective_used}),
#'         \code{auto_used}, \code{auto_decision}, \code{auto_rule}.
#'   \item \code{dropped_alpha0_for_relaxed}: whether alpha = 0 was removed because \code{relaxed = TRUE}.
#'   \item \code{schema}: list(\code{feature_names}, \code{terms_str}, \code{xlevels}, \code{contrasts}, \code{terms_hash}) for safe predict.
#'   \item \code{sampling_schema}: list(
#'         \code{predictor_vars}, \code{var_classes},
#'         \code{num_ranges} = rbind(min=..., max=...) for numeric raw predictors,
#'         \code{factor_levels} = list(...) for factor/character raw predictors).
#'   \item \code{diagnostics}: list with \code{k_summary} (median and IQR of selected size),
#'         \code{fallback_rate}, \code{n_eff_summary}, \code{alpha_freq}, \code{relax_gamma_freq}.
#'   \item \code{actual_y}, \code{training_X}, \code{y_pred}, \code{y_pred_debiased}, \code{nobs}, \code{nparm},
#'         \code{formula}, \code{terms}, \code{xlevels}, \code{contrasts}.
#'   \item \code{used_bigexp_spec}: logical flag indicating whether a \code{bigexp_spec} was used.
#' }
#'
#' @details
#' In many applications, \code{SVEMnet()} is used as part of a closed-loop
#' optimization workflow:
#' models are fit on current experimental data, whole-model tests (WMT) are
#' optionally used to reweight response goals, and then
#' \code{svem_optimize_random()} proposes both optimal and exploration
#' candidates for the next experimental round.
#' See the \code{lipid_screen} help page for a worked example.
#'
#' SVEM applies fractional bootstrap weights to training data and anti-correlated
#' weights for validation when \code{weight_scheme = "SVEM"}. For each bootstrap, glmnet
#' paths are fit for each alpha in \code{glmnet_alpha}, and the lambda (and, if \code{relaxed = TRUE},
#' relaxed gamma) minimizing a weighted validation criterion is selected.
#'
#' \strong{Weighting schemes.}
#' With \code{weight_scheme = "SVEM"} (the default), SVEMnet uses the
#' fractional-random-weight (FRW) construction with an explicit
#' self-validation split: for each bootstrap replicate and observation
#' \eqn{i}, a shared \eqn{U_i \sim \mathrm{Uniform}(0,1)} is drawn and
#' converted into anti-correlated train/validation copies
#' \eqn{w_i^{\mathrm{train}} = -\log U_i} and
#' \eqn{w_i^{\mathrm{valid}} = -\log(1-U_i)}, each rescaled so that their
#' mean is 1. This keeps all rows in every fit while inducing a stable
#' out-of-bag–style validation set for selecting \eqn{\lambda} (and, if used,
#' the relaxed \eqn{\gamma}).
#'
#' With \code{weight_scheme = "FRW_plain"}, SVEMnet instead uses a single FRW
#' weight vector \eqn{w_i = -\log U_i} for \emph{both} training and validation
#' (one weighted fit, no self-validation split). It is included mainly for
#' historical comparison and method-teaching; for small-sample prediction we
#' recommend the default \code{"SVEM"} scheme.
#'
#' Finally, \code{weight_scheme = "Identity"} sets both training and
#' validation weights to 1. In combination with \code{nBoot = 1}, this
#' effectively wraps a single \code{glmnet} fit and chooses \eqn{\lambda}
#' (and, for Gaussian models, the relaxed \eqn{\gamma}) by the chosen
#' information criterion (wAIC or wBIC), without any bootstrap variation.
#' This can be useful when you want classical AIC/BIC selection on top of
#' a deterministic expansion, but do not want or need ensembling.
#'
#' Predictors are always standardized internally via
#' \code{glmnet::glmnet(..., standardize = TRUE)}.
#'
#' When \code{relaxed = TRUE} and \code{family = "gaussian"},
#' \code{coef(fit, s = lambda, gamma = g)} is used to obtain the coefficient
#' path at each relaxed gamma in the internal grid (by default
#' \code{c(0.2, 0.6, 1)}). Metrics are computed from validation-weighted
#' errors and model size is taken as the number of nonzero coefficients
#' including the intercept (support size), keeping selection consistent
#' between standard and relaxed paths. For \code{family = "binomial"},
#' relaxed fits are currently disabled for numerical stability, so only
#' the standard glmnet path is used even if \code{relaxed = TRUE}.
#'
#' Automatic objective rule (\code{"auto"}): This function uses a single ratio
#' cutoff, \code{auto_ratio_cutoff}, applied to \code{r = n_X / p_X}, where
#' \code{p_X} is computed from the model matrix with the intercept column
#' removed. If \code{r >= auto_ratio_cutoff} the selection criterion is wAIC;
#' otherwise it is wBIC.
#'
#' Implementation notes for safety:
#' - The training terms object is stored with environment set to \code{baseenv()} to avoid
#'   accidental lookups in the calling environment.
#' - A compact schema (feature names, xlevels, contrasts) is stored to let \code{predict()}
#'   reconstruct the design matrix deterministically.
#' - A lightweight sampling schema (numeric ranges and factor levels for raw predictors)
#'   is cached to enable random-table generation without needing the original data.
#'
#' For \code{family = "gaussian"}, the loss used in validation is a weighted SSE,
#' and wAIC/wBIC are computed from a Gaussian log-likelihood proxy.
#'
#' For \code{family = "binomial"}, the validation loss is the weighted binomial
#' deviance, and wAIC/wBIC are computed as \code{deviance + 2 * k} or
#' \code{deviance + log(n_eff) * k}, where \code{k} is the number of active
#' parameters (1 for the intercept plus the number of nonzero parameters) and
#' \code{n_eff} is the effective validation size. The response must be numeric
#' 0/1 or a two-level factor; internally it is converted to 0/1.
#'
#' @template ref-svem
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
#'   data            = dat,
#'   factorial_order = 3,  # up to 3-way interactions among X1, X2, X3
#'   polynomial_order = 3, # include I(X^3) for each continuous predictor
#'   include_pc_3way = FALSE
#' )
#'
#' # Fit using the spec (auto-prepares data)
#' fit_y <- SVEMnet(
#'   spec, dat,
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 50,
#'   objective     = "auto",
#'   weight_scheme = "SVEM",
#'   relaxed       = FALSE
#' )
#'
#' # A second, independent response over the same expansion
#' set.seed(99)
#' dat$y2 <- 0.5 + 1.4*X1 - 0.6*X2 + 0.2*X3 + rnorm(n, 0, 0.4)
#' fit_y2 <- SVEMnet(
#'   spec, dat, response = "y2",
#'   glmnet_alpha  = c(1, 0.5),
#'   nBoot         = 50,
#'   objective     = "auto",
#'   weight_scheme = "SVEM",
#'   relaxed       = FALSE
#' )
#'
#' p1  <- predict(fit_y,  dat)
#' p2  <- predict(fit_y2, dat, debias = TRUE)
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
#'
#' #' ## ---- Binomial example ------------------------------------------------
#' set.seed(2)
#' n  <- 120
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' eta <- -0.3 + 1.1*X1 - 0.8*X2 + 0.5*X1*X3
#' p   <- plogis(eta)
#' yb  <- rbinom(n, 1, p)
#' db  <- data.frame(yb = yb, X1 = X1, X2 = X2, X3 = X3)
#'
#' fit_b <- SVEMnet(
#'   yb ~ (X1 + X2 + X3)^2, db,
#'   nBoot = 50, glmnet_alpha = c(1, 0.5), relaxed = FALSE, family = "binomial"
#' )
#'
#' ## Probabilities, link, and classes
#' p_resp <- predict(fit_b, db, type = "response")
#' p_link <- predict(fit_b, db, type = "link")
#' y_hat  <- predict(fit_b, db, type = "class")      # 0/1 labels (no SE/interval)
#'
#' ## Mean-aggregation with uncertainty on probability scale
#' out_b <- predict(fit_b, db, type = "response",
#'                  se.fit = TRUE, interval = TRUE, level = 0.9)
#' str(out_b)
#' }
#'
#' @export
SVEMnet <- function(formula, data,
                        nBoot = 200,
                        glmnet_alpha = c( 0.5, 1),
                        weight_scheme = c("SVEM", "FRW_plain", "Identity"),
                        objective = c("auto", "wAIC", "wBIC", "wSSE"),
                        auto_ratio_cutoff = 1.3,
                        relaxed = TRUE,
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

  if (!is.logical(relaxed) || length(relaxed) != 1L || is.na(relaxed)) {
    stop("'relaxed' must be a single logical TRUE or FALSE.")
  }
  relax_flag <- isTRUE(relaxed)

  ## For stability, ignore relaxed fits for binomial
  if (fam_name == "binomial" && relax_flag) {
    warning(
      "SVEMnet: 'relaxed = TRUE' is ignored for family = 'binomial'; ",
      "using relaxed = FALSE for stability."
    )
    relax_flag <- FALSE
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

  if (!is.numeric(auto_ratio_cutoff) || length(auto_ratio_cutoff) != 1L ||
      !is.finite(auto_ratio_cutoff) || auto_ratio_cutoff <= 0) {
    stop("'auto_ratio_cutoff' must be a single finite number > 0.")
  }
  auto_ratio_cutoff <- as.numeric(auto_ratio_cutoff)

  ## ---- model frame / matrix construction ----
  data <- as.data.frame(data)
  using_spec <- inherits(formula, "bigexp_spec")
  contrasts_opts_used <- NULL

  if (using_spec) {
    spec <- formula
    if (is.null(response)) {
      f_use <- spec$formula
    } else {
      rhs_txt <- paste(deparse(spec$formula[[3]]), collapse = " ")
      f_use <- stats::as.formula(paste(response, "~", rhs_txt))
    }
    prep <- bigexp_prepare(spec, data, unseen = unseen)

    spec_settings <- spec$settings
    if (is.null(spec_settings)) spec_settings <- list()

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
    f_use <- formula
    mf    <- stats::model.frame(f_use, data, na.action = stats::na.omit)
    if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")
    contrasts_opts_used <- getOption("contrasts")
    X     <- stats::model.matrix(f_use, mf)
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
  auto_rule      <- c(ratio_cutoff = auto_ratio_cutoff)
  auto_used      <- identical(objective, "auto")
  auto_decision  <- NA_character_
  objective_used <- if (auto_used) {
    r_tmp <- n / p
    if (is.finite(r_tmp) && r_tmp >= auto_ratio_cutoff) {
      auto_decision <- "wAIC"; "wAIC"
    } else {
      auto_decision <- "wBIC"; "wBIC"
    }
  } else objective

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

  .support_size_one <- function(beta_col, base_tol = 1e-7, count_intercept = TRUE) {
    if (!is.numeric(beta_col) || length(beta_col) < 1L) return(NA_integer_)
    slope <- beta_col[-1L]
    tol_j <- base_tol * max(1, max(abs(slope)))
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
    best_relax_gammas = if (isTRUE(relax_flag)) best_relax_gammas else
      rep(NA_real_, length(best_alphas)),
    weight_scheme     = weight_scheme,
    relaxed           = isTRUE(relax_flag),
    dropped_alpha0_for_relaxed = dropped_alpha0_for_relaxed,
    objective_input   = objective,
    objective_used    = objective_used,
    objective         = objective_used,
    auto_used         = auto_used,
    auto_decision     = if (auto_used) auto_decision else NA_character_,
    auto_rule         = auto_rule,
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
