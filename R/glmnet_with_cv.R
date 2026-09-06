#' Fit a glmnet Model with Repeated Cross-Validation
#'
#' Repeated K-fold cross-validation over per-alpha lambda paths, with a
#' combined 1-SE rule across repeats. Preserves fields expected by
#' \code{predict.svem_model()} and internal prediction helpers. When
#' \code{relaxed = TRUE}, the complete relaxed cross-validation surfaces are
#' aggregated across repeats and alpha, lambda, and gamma are selected jointly.
#' This mirrors \code{cv.glmnet(relax = TRUE)} rather than choosing lambda from
#' the ordinary (gamma = 1) curve before choosing gamma.
#'
#' This function is a convenience wrapper around \code{glmnet()} and
#' \code{cv.glmnet()} that returns an object in the same structural format as
#' \code{SVEMnet()} (class \code{"svem_model"}). It is intended for:
#' \itemize{
#'   \item direct comparison of standard cross-validated glmnet fits to SVEMnet
#'         models using the same prediction and schema tools, or
#'   \item users who want a repeated-\code{cv.glmnet()} workflow without any
#'         SVEM weighting or bootstrap ensembling.
#' }
#' It is not called internally by the SVEM bootstrap routines.
#'
#' @param formula Model formula.
#' @param data Data frame containing the variables in the model.
#' @param glmnet_alpha Numeric vector of Elastic Net mixing parameters
#'   (alphas) in \code{[0,1]}; default \code{c(0.5, 1)}. When
#'   \code{relaxed = TRUE}, any \code{alpha = 0} (ridge) is dropped with a
#'   warning.
#' @param standardize Logical passed to \code{glmnet()} and \code{cv.glmnet()}
#'   (default \code{TRUE}).
#' @param nfolds Requested number of CV folds (default \code{10}). Internally
#'   constrained so that there are at least about 3 observations per fold and
#'   at least 5 folds when possible.
#' @param repeats Number of independent CV repeats (default \code{5}). Each
#'   repeat reuses the same folds across all alphas for paired comparisons.
#' @param choose_rule Character; how to choose the tuning point within each
#'   alpha:
#'   \itemize{
#'     \item \code{"min"}: tuning point minimizing the cross-validated
#'           criterion.
#'     \item \code{"1se"}: largest lambda (then largest gamma for relaxed
#'           fits) within 1 combined SE of the minimum, where the SE includes
#'           both within- and between-repeat variability.
#'   }
#'   Default is \code{"min"}. In small-mixture simulations, the 1-SE rule
#'   tended to increase RMSE on held-out data, so \code{"min"} is used as the
#'   default here.
#' @param seed Optional integer seed for reproducible fold IDs (and the
#'   ridge fallback, if used).
#' @param exclude Optional vector or function for \code{glmnet}'s
#'   \code{exclude=} argument. If a function, \code{cv.glmnet()} applies it
#'   inside each training fold (requires \code{glmnet} >= 4.1-2).
#' @param relaxed Logical; if \code{TRUE}, call \code{glmnet()} and
#'   \code{cv.glmnet()} with \code{relax = TRUE} and optionally a
#'   \code{gamma} path (default \code{FALSE}). The repeated CV criterion is
#'   aggregated over the full lambda-by-gamma surface and both parameters are
#'   selected jointly. If \code{cv.glmnet(relax = TRUE)} fails for a particular
#'   repeat/alpha, the function retries that fit without relaxation; incomplete
#'   non-relaxed fallback surfaces are not mixed into a successful relaxed
#'   search. The number of fallbacks is recorded in
#'   \code{meta$relax_cv_fallbacks}.
#' @param relax_gamma Optional numeric vector passed as \code{gamma=} to
#'   \code{glmnet()} and \code{cv.glmnet()} when \code{relaxed = TRUE}. If
#'   \code{NULL}, glmnet's internal default gamma grid is used.
#' @param family Model family: either \code{"gaussian"} or \code{"binomial"},
#'   or the corresponding \code{stats::gaussian()} or \code{stats::binomial()}
#'   family objects with canonical links. For Gaussian, \code{y} must be
#'   numeric. For binomial, \code{y} must be 0/1 numeric, logical, or a factor
#'   with exactly 2 levels (the second level is treated as 1). Non-canonical
#'   links are not supported.
#' @param ... Additional arguments forwarded to both \code{cv.glmnet()} and
#'   \code{glmnet()}, for example: \code{parallel}, \code{type.measure},
#'   \code{intercept}, \code{lower.limits}, \code{upper.limits},
#'   \code{penalty.factor}, \code{standardize.response}, \code{keep}, and so
#'   on. Any user-supplied \code{weights} or \code{offset} are ignored with a
#'   warning: \code{glmnet_with_cv()} returns objects meant for side-by-side
#'   comparison with \code{SVEMnet()} fits, whose weighting is controlled by
#'   the SVEM scheme, and downstream prediction helpers do not carry
#'   observation weights or offsets.
#'   Glmnet algorithm-control values such as \code{maxit}, \code{thresh},
#'   \code{dfmax}, and \code{pmax} may be supplied directly for compatibility
#'   or via \code{control = list(...)} with glmnet versions that support it;
#'   they are routed in a version-compatible way when possible.
#'   Argument names that the installed \code{glmnet}
#'   does not recognize are ignored with a warning (misspelling protection).
#'
#' @details
#' The basic workflow is:
#' \enumerate{
#'   \item Generate one set of CV fold IDs per repeat. Reuse each repeat's
#'         folds across all values in \code{glmnet_alpha}.
#'   \item For each alpha, run \code{repeats} \code{cv.glmnet()}
#'         fits, align the lambda paths and, for relaxed fits, the gamma paths,
#'         and aggregate the CV surfaces.
#'   \item At each lambda-gamma point, compute a combined SE that accounts for
#'         both within-repeat and between-repeat variability.
#'   \item Apply \code{choose_rule} (\code{"min"} or \code{"1se"}) to select
#'         lambda and gamma for that alpha, then choose the best alpha by
#'         comparing these per-alpha scores.
#' }
#'
#' Special cases and fallbacks:
#' \itemize{
#'   \item If there are no predictors after \code{model.matrix()} (an
#'         intercept-only model), the function returns an intercept-only fit
#'         without calling \code{glmnet()}, along with a minimal schema for
#'         safe prediction.
#'   \item If all \code{cv.glmnet()} attempts fail for every alpha (a rare
#'         edge case), the function falls back to a manual ridge
#'         (\code{alpha = 0}) CV search over a fixed lambda grid and returns
#'         the best ridge solution. For Gaussian models this search uses a
#'         mean-squared-error criterion; for binomial models it uses a
#'         negative log-likelihood (deviance-equivalent) criterion.
#' }
#'
#' Family-specific behavior:
#' \itemize{
#'   \item For the Gaussian family, an optional calibration \code{lm(y ~ y_pred)}
#'         is fit on the training data (when there is sufficient variation), and
#'         both \code{y_pred} and \code{y_pred_debiased} are stored.
#'   \item For the binomial family, \code{y_pred} is always on the probability
#'         (response) scale and debiasing is not applied. Both the primary
#'         cross-validation and any ridge fallback use deviance-style criteria
#'         (binomial negative log-likelihood) rather than squared error.
#' }
#'
#' Design-matrix schema and contrasts:
#' \itemize{
#'   \item The training \code{terms} are stored with a compact evaluation
#'         environment containing required formula transformations.
#'   \item Factor and character levels are recorded in \code{xlevels} for
#'         safe prediction.
#'   \item Per-factor contrasts are stored in \code{contrasts}, normalized
#'         so that any contrasts recorded as character names are converted
#'         back to contrast functions at prediction time.
#' }
#'
#' The returned object inherits classes \code{"svem_cv"} and \code{"svem_model"}
#' and is designed to be compatible with SVEMnet prediction and schema
#' utilities. It is a standalone, standard glmnet CV workflow that does not use
#' SVEM-style bootstrap weighting or ensembling.
#'
#' @return A list of class \code{c("svem_cv","svem_model")} with elements:
#' \itemize{
#'   \item \code{parms} Named numeric vector of coefficients (including
#'         \code{"(Intercept)"}).
#'   \item \code{glmnet_alpha} Numeric vector of alphas searched.
#'   \item \code{best_alpha} Numeric; winning alpha.
#'   \item \code{best_lambda} Numeric; winning lambda.
#'   \item \code{best_gamma} Numeric; winning relaxed mixing parameter, or
#'         1 for a non-relaxed fit. It is \code{NA} for structural
#'         intercept-only fits.
#'   \item \code{y_pred} In-sample predictions from the returned coefficients
#'         (fitted values for Gaussian; probabilities for binomial).
#'   \item \code{debias_fit} For Gaussian, an optional \code{lm(y ~ y_pred)}
#'         calibration model; \code{NULL} otherwise.
#'   \item \code{y_pred_debiased} If \code{debias_fit} exists, its fitted
#'         values; otherwise \code{NULL}.
#'   \item \code{cv_summary} Named list (one element per alpha) of data frames
#'         with columns \code{gamma}, \code{lambda}, \code{mean_cvm}, \code{sd_cvm},
#'         \code{se_combined}, \code{n_repeats}, \code{idx_min},
#'         \code{idx_1se}.
#'   \item \code{formula} Original modeling formula.
#'   \item \code{terms} Training \code{terms} object with a compact,
#'         serializable formula-evaluation environment.
#'   \item \code{training_X} Training design matrix (without intercept column).
#'   \item \code{actual_y} Training response vector used for glmnet:
#'         numeric \code{y} for Gaussian, or 0/1 numeric \code{y} for
#'         binomial.
#'   \item \code{xlevels} Factor and character levels seen during training
#'         (for safe prediction).
#'   \item \code{contrasts} Contrasts used for factor predictors during
#'         training.
#'   \item \code{schema} List
#'         \code{list(feature_names, terms_str, xlevels, contrasts, terms_hash)}
#'         for deterministic prediction.
#'   \item \code{note} Character vector of notes (for example, dropped rows,
#'         intercept-only path, ridge fallback, relaxed-coefficient source).
#'   \item \code{meta} List with fields such as \code{nfolds}, \code{repeats},
#'         \code{rule}, \code{family}, \code{relaxed},
#'         \code{relax_cv_fallbacks}, and \code{cv_object} (the final
#'         \code{cv.glmnet()} object when \code{relaxed = TRUE} and
#'         \code{keep = TRUE}, otherwise \code{NULL}).
#'   \item \code{diagnostics} List of simple diagnostics for the selected
#'         model, currently including:
#'         \itemize{
#'           \item \code{k_final}: number of coefficients estimated as
#'                 nonzero \emph{including} the intercept.
#'           \item \code{k_final_no_intercept}: number of nonzero
#'                 slope coefficients (excludes the intercept).
#'         }
#'   \item \code{family} Character scalar giving the resolved family
#'         (\code{"gaussian"} or \code{"binomial"}), mirroring
#'         \code{meta$family}.
#' }
#'
#' @template ref-svem
#'
#' @examples
#' set.seed(123)
#' n <- 60; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1, -1, rep(0, p - 2))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' df_ex <- data.frame(y = y, X)
#' colnames(df_ex) <- c("y", paste0("x", 1:p))
#'
#' # Gaussian example, v1-like behavior: choose_rule = "min"
#' fit_min <- glmnet_with_cv(
#'   y ~ ., df_ex,
#'   glmnet_alpha = 1,
#'   nfolds = 5,
#'   repeats = 1,
#'   choose_rule = "min",
#'   seed = 42,
#'   family = "gaussian"
#' )
#'
#' # Gaussian example, relaxed path with gamma search
#' fit_relax <- glmnet_with_cv(
#'   y ~ ., df_ex,
#'   glmnet_alpha = 1,
#'   nfolds = 5,
#'   repeats = 1,
#'   relaxed = TRUE,
#'   seed = 42,
#'   family = "gaussian"
#' )
#'
#' # Binomial example (numeric 0/1 response)
#' set.seed(456)
#' n2 <- 80; p2 <- 5
#' X2 <- matrix(rnorm(n2 * p2), n2, p2)
#' beta2 <- c(1.0, -1.5, rep(0, p2 - 2))
#' linpred <- as.numeric(X2 %*% beta2)
#' prob <- plogis(linpred)
#' y_bin <- rbinom(n2, size = 1, prob = prob)
#' df_bin <- data.frame(y = y_bin, X2)
#' colnames(df_bin) <- c("y", paste0("x", 1:p2))
#'
#' fit_bin <- glmnet_with_cv(
#'   y ~ ., df_bin,
#'   glmnet_alpha = 1,
#'   nfolds = 5,
#'   repeats = 1,
#'   seed = 99,
#'   family = "binomial"
#' )
#'
#' @importFrom stats model.frame model.response model.matrix var coef predict lm delete.response
#' @importFrom glmnet glmnet cv.glmnet
#' @export
glmnet_with_cv <- function(formula, data,
                           glmnet_alpha = c(0.5, 1),
                           standardize = TRUE,
                           nfolds = 10,
                           repeats = 5,
                           choose_rule = c("min", "1se"),
                           seed = NULL,
                           exclude = NULL,
                           relaxed = FALSE,
                           relax_gamma = NULL,
                           family = c("gaussian", "binomial"),
                           ...) {

  choose_rule <- match.arg(choose_rule)
  nfolds <- .svem_integer_scalar(nfolds, "nfolds", min = 3L)
  repeats <- .svem_integer_scalar(repeats, "repeats", min = 1L)
  standardize <- .svem_logical_scalar(standardize, "standardize")
  relaxed <- .svem_logical_scalar(relaxed, "relaxed")
  seed <- .svem_integer_scalar(seed, "seed", min = 0L, allow_null = TRUE)

  if (!is.numeric(glmnet_alpha) || any(!is.finite(glmnet_alpha))) {
    stop("`glmnet_alpha` must be numeric and finite.", call. = FALSE)
  }
  if (any(glmnet_alpha < 0 | glmnet_alpha > 1)) {
    stop("`glmnet_alpha` values must be in [0, 1].", call. = FALSE)
  }
  glmnet_alpha <- unique(as.numeric(glmnet_alpha))
  if (!length(glmnet_alpha)) glmnet_alpha <- 1

  if (!is.null(relax_gamma) &&
      (!is.numeric(relax_gamma) || !length(relax_gamma) ||
       any(!is.finite(relax_gamma)) ||
       any(relax_gamma < 0 | relax_gamma > 1))) {
    stop("`relax_gamma` must be NULL or contain finite values in [0, 1].",
         call. = FALSE)
  }

  # --- Helper: support size (number of active coefficients) ---
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


  # --- Build model frame & matrix ---
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  if (nrow(mf) < 2L) {
    stop("glmnet_with_cv() requires at least two complete observations.",
         call. = FALSE)
  }
  drop_n <- nrow(as.data.frame(data)) - nrow(mf)
  drop_msg <- if (drop_n > 0) sprintf("Dropped %d row(s) due to missing values.", drop_n) else NULL

  y <- stats::model.response(mf)
  terms_full <- attr(mf, "terms")
  if (!identical(attr(terms_full, "intercept"), 1L)) {
    stop(
      "glmnet_with_cv() follows the published intercept convention; formulas using `~ 0` or `- 1` are not supported.",
      call. = FALSE
    )
  }
  if (is.matrix(y) || is.array(y) || length(y) != nrow(mf)) {
    stop("glmnet_with_cv() requires a single univariate response; matrix/multivariate responses are not supported.",
         call. = FALSE)
  }

  compact_formula <- .svem_compact_formula_terms(
    formula, terms_full, data_names = names(as.data.frame(data))
  )
  formula_store <- compact_formula$formula
  terms_store <- compact_formula$terms

  X <- stats::model.matrix(formula, mf)
  # Capture the contrasts attribute BEFORE dropping the intercept column:
  # matrix subsetting keeps only dim/dimnames, so reading it afterwards
  # always returns NULL (same fix as "fix #1" in SVEMnet.R).
  contrasts_used_full <- attr(X, "contrasts")
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col)) X <- X[, -intercept_col, drop = FALSE]
  n <- nrow(X); p <- ncol(X)

  predictor_vars <- base::all.vars(stats::delete.response(terms_full))
  data_frame <- as.data.frame(data)
  data_frame_used <- data_frame
  omitted_rows <- as.integer(attr(mf, "na.action"))
  if (length(omitted_rows) && nrow(data_frame_used) >= max(omitted_rows)) {
    data_frame_used <- data_frame_used[-omitted_rows, , drop = FALSE]
  }
  var_classes <- stats::setNames(vapply(predictor_vars, function(v) {
    if (v %in% names(data_frame)) class(data_frame[[v]])[1L] else NA_character_
  }, character(1L)), predictor_vars)

  # --- Collect and sanitize dots ---
  dots <- list(...)

  # User weights and offsets are not supported: the returned object is meant
  # for side-by-side comparison with SVEMnet() fits (whose weighting is
  # controlled by the SVEM scheme), and the prediction helpers do not carry
  # observation weights or offsets.
  if ("weights" %in% names(dots)) {
    warning("Ignoring user 'weights'; glmnet_with_cv() does not support observation weights.")
    dots$weights <- NULL
  }
  if ("offset" %in% names(dots)) {
    warning("Ignoring user 'offset'; glmnet_with_cv() does not support offsets.")
    dots$offset <- NULL
  }

  # glmnet()/cv.glmnet() silently swallow unknown arguments; catch misspellings
  unknown <- .svem_unknown_glmnet_args(
    dots,
    funs = list(glmnet::glmnet, glmnet::cv.glmnet),
    context = "glmnet via glmnet_with_cv"
  )
  if (length(unknown)) dots <- dots[setdiff(names(dots), unknown)]

  # --- Family handling: restrict to gaussian / binomial with canonical links ---
  .coerce_binomial_01 <- function(y) {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop("Binomial glmnet_with_cv requires a factor with exactly 2 levels or 0/1 numeric/logical y.")
      }
      return(as.integer(y == levels(y)[2L]))
    } else if (is.logical(y)) {
      return(as.integer(y))
    } else if (is.numeric(y)) {
      uy <- sort(unique(y))
      if (!all(uy %in% c(0, 1))) {
        stop("Binomial glmnet_with_cv requires numeric y coded as 0/1.")
      }
      return(as.integer(y))
    } else {
      stop("Unsupported y type for binomial; use 0/1 numeric, logical, or a 2-level factor.")
    }
  }

  fam_raw <- family

  if (inherits(fam_raw, "family")) {
    fam_name <- fam_raw$family
    fam_link <- fam_raw$link

    if (identical(fam_name, "gaussian")) {
      if (!identical(fam_link, "identity")) {
        stop("glmnet_with_cv currently supports gaussian() only with the canonical identity link.")
      }
      fam <- "gaussian"
    } else if (identical(fam_name, "binomial")) {
      if (!identical(fam_link, "logit")) {
        stop("glmnet_with_cv currently supports binomial() only with the canonical logit link.")
      }
      fam <- "binomial"
    } else {
      stop("glmnet_with_cv currently supports family = 'gaussian' or 'binomial' only.")
    }
  } else {
    fam_char <- as.character(fam_raw)[1L]
    fam <- match.arg(fam_char, c("gaussian", "binomial"))
  }

  # Coerce y for modeling: numeric for gaussian; 0/1 numeric for binomial
  if (identical(fam, "gaussian")) {
    if (!is.numeric(y)) stop("Gaussian glmnet_with_cv requires numeric y.")
    y_glm <- as.numeric(y)
  } else {  # binomial
    y_glm <- .coerce_binomial_01(y)
  }

  is_gaussian <- identical(fam, "gaussian")

  if (any(!is.finite(y_glm)) || any(!is.finite(X))) {
    stop("Non-finite values in response/predictors after NA handling.",
         call. = FALSE)
  }
  storage.mode(X) <- "double"

  # --- Helper for intercept-only returns (used in both places) ---
  .cv_return_intercept_only <- function(
    fam, y_glm, X, mf, glmnet_alpha, choose_rule,
    drop_msg, nfolds_val, repeats_val, relaxed, note_extra = NULL
  ) {
    n_local <- nrow(X)

    if (identical(fam, "gaussian")) {
      intercept <- mean(y_glm)
      y_pred <- rep(intercept, n_local)
    } else {  # binomial: intercept on logit scale, y_pred as probabilities
      p_hat <- mean(y_glm)
      p_hat <- pmin(pmax(p_hat, 1e-6), 1 - 1e-6)
      intercept <- stats::qlogis(p_hat)
      y_pred <- rep(p_hat, n_local)
    }

    # Clean terms and schema for consistency with normal path
    terms_clean <- terms_store
    feature_names <- colnames(X)

    # Record factor AND character levels for completeness (even if unused)
    xlevels <- list()
    vars_in_terms <- predictor_vars
    for (v in vars_in_terms) {
      if (v %in% names(data_frame_used)) {
        if (is.factor(data_frame_used[[v]])) {
          xlevels[[v]] <- levels(data_frame_used[[v]])
        } else if (is.character(data_frame_used[[v]])) {
          xlevels[[v]] <- sort(unique(as.character(data_frame_used[[v]])))
        }
      }
    }
    contrasts_used <- contrasts_used_full
    if (!is.null(contrasts_used)) {
      contrasts_used <- lapply(contrasts_used, function(ci) {
        if (is.character(ci)) get(ci, mode = "function") else ci
      })
    }

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
      predictor_vars = predictor_vars,
      var_classes   = var_classes,
      terms_hash    = safe_hash(terms_str)
    )

    note_vec <- c("intercept_only", drop_msg, note_extra)
    note_vec <- note_vec[!is.na(note_vec)]

    # Diagnostics: model size (intercept only)
    k_final  <- 1L
    k_final0 <- 0L
    diagnostics <- list(
      k_final              = k_final,
      k_final_no_intercept = k_final0,
      selected_gamma       = NA_real_
    )

    res <- list(
      parms = c("(Intercept)" = intercept),
      glmnet_alpha = glmnet_alpha,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      best_gamma = NA_real_,
      y_pred = y_pred,
      debias_fit = NULL,
      y_pred_debiased = NULL,
      cv_summary = list(),
      formula = formula_store,
      terms = terms_clean,
      training_X = X,
      actual_y = y_glm,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = note_vec,
      meta = list(
        nfolds = nfolds_val,
        repeats = repeats_val,
        rule = choose_rule,
        family = fam,
        relaxed = relaxed,
        relax_cv_fallbacks = 0L,
        best_gamma = NA_real_,
        cv_object = NULL
      ),
      diagnostics = diagnostics,
      family = fam
    )
    class(res) <- c("svem_cv","svem_model")
    res
  }

  # --- Intercept-only path (early return) ---
  if (p == 0L) {
    return(.cv_return_intercept_only(
      fam         = fam,
      y_glm       = y_glm,
      X           = X,
      mf          = mf,
      glmnet_alpha = glmnet_alpha,
      choose_rule = choose_rule,
      drop_msg    = drop_msg,
      nfolds_val  = NA_integer_,
      repeats_val = NA_integer_,
      relaxed     = relaxed,
      note_extra  = NULL
    ))
  }

  if (identical(fam, "gaussian") && length(unique(y_glm)) < 2L) {
    return(.cv_return_intercept_only(
      fam          = fam,
      y_glm        = y_glm,
      X            = X,
      mf           = mf,
      glmnet_alpha = glmnet_alpha,
      choose_rule  = choose_rule,
      drop_msg     = drop_msg,
      nfolds_val   = NA_integer_,
      repeats_val  = NA_integer_,
      relaxed      = relaxed,
      note_extra   = "gaussian_constant_intercept_only"
    ))
  }

  # --- Binomial guardrails (small-sample stability) ---
  #
  # In tiny Bernoulli samples it is possible (especially under extreme
  # probabilities) to observe only one outcome class. In that case, logistic
  # glmnet / cv.glmnet can fail. We return an intercept-only model on the
  # probability scale rather than erroring.
  #
  # Additionally, some glmnet versions are more stable for binomial fits when
  # y is provided as a 2-column success/failure matrix. We keep y_glm as a 0/1
  # vector for storage, but use y_fit for model fitting.
  y_fit <- y_glm
  if (identical(fam, "binomial")) {
    uy <- unique(y_glm)
    uy <- uy[is.finite(uy)]
    if (length(uy) < 2L) {
      # Compute the effective folds/repeats for meta fields, even though we
      # are returning early.
      max_folds_tmp  <- max(3L, floor(n / 3))
      nfolds_eff_tmp <- min(max(5L, min(nfolds, n)), max_folds_tmp)
      repeats_tmp    <- repeats

      return(.cv_return_intercept_only(
        fam          = fam,
        y_glm        = y_glm,
        X            = X,
        mf           = mf,
        glmnet_alpha = glmnet_alpha,
        choose_rule  = choose_rule,
        drop_msg     = drop_msg,
        nfolds_val   = nfolds_eff_tmp,
        repeats_val  = repeats_tmp,
        relaxed      = relaxed,
        note_extra   = "binomial_one_class_intercept_only"
      ))
    }
    y_fit <- cbind(1 - y_glm, y_glm)
  }

  # --- Safe folds & repeats ---
  max_folds  <- max(3L, floor(n / 3))
  nfolds_eff <- min(max(5L, min(nfolds, n)), max_folds)

  # Binomial CV needs at least 2 observations in each class; otherwise
  # at least one fold's training set will contain a single class and
  # glmnet's binomial CV (especially relaxed CV) can fail.
  if (identical(fam, "binomial")) {
    n1 <- sum(y_glm == 1L, na.rm = TRUE)
    n0 <- sum(y_glm == 0L, na.rm = TRUE)
    if (n1 < 2L || n0 < 2L) {
      return(.cv_return_intercept_only(
        fam          = fam,
        y_glm        = y_glm,
        X            = X,
        mf           = mf,
        glmnet_alpha = glmnet_alpha,
        choose_rule  = choose_rule,
        drop_msg     = drop_msg,
        nfolds_val   = nfolds_eff,
        repeats_val  = repeats,
        relaxed      = relaxed,
        note_extra   = "binomial_min_class_lt2_intercept_only"
      ))
    }
  }

  # --- Remove the unsupported ridge member from a relaxed search ---
  if (isTRUE(relaxed) && any(glmnet_alpha == 0)) {
    warning("Dropping alpha = 0 (ridge) for relaxed fits; ridge + relaxed is not supported here.")
    glmnet_alpha <- glmnet_alpha[glmnet_alpha != 0]
    if (!length(glmnet_alpha)) glmnet_alpha <- 1
  }

  glmnet_design <- .svem_prepare_glmnet_design(
    X, args = dots, exclude = exclude
  )
  X_fit <- glmnet_design$x
  dots <- glmnet_design$args
  exclude_fit <- glmnet_design$exclude

  # --- Defaults; allow dots to override ---
  base_cv_args <- list(
    standardize  = standardize,
    family       = fam,
    type.measure = if (identical(fam, "binomial")) "deviance" else "mse",
    relax        = isTRUE(relaxed)
  )
  if (isTRUE(relaxed) && !is.null(relax_gamma)) base_cv_args$gamma <- relax_gamma
  cv_base_args <- utils::modifyList(base_cv_args, dots, keep.null = TRUE)
  cv_base_args <- .glmnet_prepare_call_args(
    cv_base_args,
    glmnet::cv.glmnet,
    default_control = list(maxit = 1e6),
    direct_control = c("thresh", "maxit", "dfmax", "pmax"),
    old_control = .glmnet_direct_control_args()
  )

  glmnet_formals <- names(formals(glmnet::glmnet))
  base_glmnet_args <- list(
    standardize = standardize,
    family      = fam,
    intercept   = TRUE,
    exclude     = exclude_fit,
    relax       = isTRUE(relaxed)
  )
  if (isTRUE(relaxed) && !is.null(relax_gamma)) base_glmnet_args$gamma <- relax_gamma
  glmnet_fit_args_full <- utils::modifyList(base_glmnet_args, dots, keep.null = TRUE)
  glmnet_fit_args_full <- .glmnet_prepare_call_args(
    glmnet_fit_args_full,
    glmnet::glmnet,
    default_control = list(maxit = 1e6)
  )
  glmnet_fit_args <- glmnet_fit_args_full[intersect(names(glmnet_fit_args_full), glmnet_formals)]
  glmnet_fit_args <- glmnet_fit_args[setdiff(names(glmnet_fit_args), c("x","y","alpha","lambda"))]

  # Helper to drop reserved names from arg lists before do.call
  drop_reserved <- function(lst, reserved) {
    if (!length(lst)) return(lst)
    lst[setdiff(names(lst), reserved)]
  }

  # Build one set of fold IDs per repeat, reused across alphas (paired comparison)
  # For binomial models, use stratified folds to avoid placing all
  # observations of a class into a single fold (which can make the
  # corresponding training set one-class and cause glmnet to fail).
  if (!is.null(seed)) set.seed(seed)
  make_foldid <- function() {
    if (identical(fam, "binomial")) {
      y01 <- as.integer(y_glm)
      idx1 <- which(y01 == 1L)
      idx0 <- which(y01 == 0L)
      foldid <- integer(length(y01))
      if (length(idx1)) foldid[idx1] <- sample(rep(seq_len(nfolds_eff), length.out = length(idx1)))
      if (length(idx0)) foldid[idx0] <- sample(rep(seq_len(nfolds_eff), length.out = length(idx0)))
      foldid
    } else {
      sample(rep(seq_len(nfolds_eff), length.out = n))
    }
  }
  foldids <- replicate(repeats, make_foldid(), simplify = FALSE)

  # --- One alpha -> repeated CV ---
  fit_alpha_repeated <- function(alpha_val) {
    reserved_cv <- c("x","y","alpha","lambda","foldid","exclude")

    repeat_surfaces <- vector("list", repeats)
    lam_ref <- NULL
    fallback_count <- 0L
    representative_cv <- NULL
    cv_name <- NULL

    extract_surface <- function(fit_cv) {
      if (inherits(fit_cv, "cv.relaxed") &&
          is.list(fit_cv$relaxed$statlist) &&
          length(fit_cv$relaxed$statlist)) {
        gammas <- as.numeric(fit_cv$relaxed$gamma)
        stats <- fit_cv$relaxed$statlist
        if (length(gammas) != length(stats)) return(NULL)
        out <- Map(function(st, gam) {
          list(
            gamma = gam,
            lambda = as.numeric(st$lambda),
            cvm = as.numeric(st$cvm),
            cvsd = as.numeric(st$cvsd)
          )
        }, stats, gammas)
        return(list(curves = out, is_relaxed = TRUE))
      }

      list(
        curves = list(list(
          gamma = 1,
          lambda = as.numeric(fit_cv$lambda),
          cvm = as.numeric(fit_cv$cvm),
          cvsd = as.numeric(fit_cv$cvsd)
        )),
        is_relaxed = FALSE
      )
    }

    for (r in seq_len(repeats)) {
      foldid <- foldids[[r]]

      fit_cv <- tryCatch(
        do.call(glmnet::cv.glmnet,
                c(list(x = X_fit, y = y_fit,
                       alpha  = alpha_val,
                       foldid = foldid,
                       exclude = exclude_fit),
                  drop_reserved(cv_base_args, reserved_cv))),
        error = function(e) {
          if (isTRUE(relaxed)) {
            warning("cv.glmnet(relax=TRUE) failed; retrying without relax for this repeat/alpha.")
            args2 <- cv_base_args; args2$relax <- NULL; args2$gamma <- NULL
            fallback_count <<- fallback_count + 1L
            tryCatch(
              do.call(glmnet::cv.glmnet,
                      c(list(x = X_fit, y = y_fit,
                             alpha  = alpha_val,
                             foldid = foldid,
                             exclude = exclude_fit),
                        drop_reserved(args2, reserved_cv))),
              error = function(e2) NULL
            )
          } else NULL
        }
      )
      if (is.null(fit_cv) || !length(fit_cv$cvm)) next

      surface <- extract_surface(fit_cv)
      if (is.null(surface)) next
      repeat_surfaces[[r]] <- surface
      if (is.null(lam_ref)) lam_ref <- as.numeric(fit_cv$lambda)
      if (is.null(cv_name)) cv_name <- fit_cv$name
      if (is.null(representative_cv) ||
          (surface$is_relaxed &&
           !inherits(representative_cv, "cv.relaxed"))) {
        representative_cv <- fit_cv
      }
    }

    if (is.null(lam_ref)) return(NULL)

    valid_surfaces <- Filter(Negate(is.null), repeat_surfaces)
    if (!length(valid_surfaces)) return(NULL)

    # A non-relaxed fallback has only gamma = 1. Mixing those incomplete
    # surfaces with complete relaxed surfaces would give different repeat
    # counts to different gamma values. Use the complete relaxed surfaces when
    # any are available; use gamma = 1 fallbacks only if every relaxed fit
    # failed.
    has_relaxed <- vapply(valid_surfaces, function(s) s[["is_relaxed"]],
                          logical(1L))
    if (isTRUE(relaxed) && any(has_relaxed)) {
      valid_surfaces <- valid_surfaces[has_relaxed]
    }

    gamma_grid <- sort(unique(unlist(lapply(
      valid_surfaces,
      function(s) vapply(s$curves, function(z) z[["gamma"]], numeric(1L))
    ))))
    if (!length(gamma_grid)) return(NULL)

    aggregate_curve <- function(gamma_value) {
      cvms <- matrix(NA_real_, nrow = length(lam_ref),
                     ncol = length(valid_surfaces))
      cvSEs <- cvms

      for (rr in seq_along(valid_surfaces)) {
        curves <- valid_surfaces[[rr]]$curves
        gi <- which(vapply(curves, function(z) {
          isTRUE(all.equal(z$gamma, gamma_value, tolerance = 1e-12))
        }, logical(1L)))
        if (!length(gi)) next
        curve <- curves[[gi[1L]]]
        idx <- match(lam_ref, curve$lambda)
        ok <- !is.na(idx)
        cvms[ok, rr] <- curve$cvm[idx[ok]]
        cvSEs[ok, rr] <- curve$cvsd[idx[ok]]
      }

      mean_cvm <- rowMeans(cvms, na.rm = TRUE)
      mean_cvm[!is.finite(mean_cvm)] <- NA_real_
      k_mean <- apply(cvms, 1L, function(z) sum(is.finite(z)))
      k_se <- apply(cvSEs, 1L, function(z) sum(is.finite(z)))
      sd_cvm <- apply(cvms, 1L, function(z) {
        z <- z[is.finite(z)]
        if (length(z) <= 1L) return(NA_real_)
        stats::sd(z)
      })

      SE_within <- rowMeans(cvSEs^2, na.rm = TRUE)
      SE_within[!is.finite(SE_within)] <- NA_real_
      if (length(valid_surfaces) == 1L) {
        SE_combined <- sqrt(SE_within / pmax(1L, k_se))
      } else {
        var_total <- sd_cvm^2
        tau2 <- pmax(0, var_total - SE_within)
        SE_combined <- sqrt(
          (SE_within / pmax(1L, k_se)) +
            (tau2 / pmax(1L, k_mean))
        )
      }
      SE_combined[!is.finite(SE_combined)] <- 0

      data.frame(
        gamma = rep(gamma_value, length(lam_ref)),
        lambda = lam_ref,
        mean_cvm = mean_cvm,
        sd_cvm = sd_cvm,
        se_combined = SE_combined,
        n_repeats = k_mean
      )
    }

    surface <- do.call(rbind, lapply(gamma_grid, aggregate_curve))
    surface <- surface[is.finite(surface$mean_cvm) &
                         surface$n_repeats > 0L, , drop = FALSE]
    if (!nrow(surface)) return(NULL)

    higher_is_better <- is.character(cv_name) && length(cv_name) == 1L &&
      cv_name %in% c("AUC", "C-index")
    selection_metric <- if (higher_is_better) {
      -surface$mean_cvm
    } else {
      surface$mean_cvm
    }

    min_value <- min(selection_metric)
    min_candidates <- which(selection_metric <= min_value)
    min_order <- order(
      surface$lambda[min_candidates],
      surface$gamma[min_candidates],
      decreasing = TRUE
    )
    idx_min <- min_candidates[min_order[1L]]

    se_tol <- surface$se_combined[idx_min]
    if (!is.finite(se_tol)) se_tol <- 0
    one_se_candidates <- which(selection_metric <= min_value + se_tol)
    one_se_order <- order(
      surface$lambda[one_se_candidates],
      surface$gamma[one_se_candidates],
      decreasing = TRUE
    )
    idx_1se <- one_se_candidates[one_se_order[1L]]

    list(
      gamma = surface$gamma,
      lambda = surface$lambda,
      lambda_path = lam_ref,
      mean_cvm = surface$mean_cvm,
      sd_cvm = surface$sd_cvm,
      se_combined = surface$se_combined,
      n_repeats = surface$n_repeats,
      idx_min = idx_min,
      idx_1se = idx_1se,
      higher_is_better = higher_is_better,
      fallback_count = fallback_count,
      cv_object = representative_cv
    )
  }

  # --- Across alphas ---
  per_alpha <- list()
  for (a in glmnet_alpha) {
    agg <- fit_alpha_repeated(a)
    if (!is.null(agg)) per_alpha[[as.character(a)]] <- agg
  }

  # --- Defensive intercept-only check before ridge fallback ---
  if (ncol(X) == 0L) {
    return(.cv_return_intercept_only(
      fam         = fam,
      y_glm       = y_glm,
      X           = X,
      mf          = mf,
      glmnet_alpha = glmnet_alpha,
      choose_rule = choose_rule,
      drop_msg    = drop_msg,
      nfolds_val  = nfolds_eff,
      repeats_val = repeats,
      relaxed     = relaxed,
      note_extra  = "ridge_fallback_bypassed"
    ))
  }

  # --- Ridge fallback if everything failed (very rare) ---
  if (!length(per_alpha)) {
    warning("All cv.glmnet attempts failed; switching to ridge fallback with manual CV.")
    if (!is.null(seed)) set.seed(seed)
    foldid <- sample(rep(seq_len(nfolds_eff), length.out = n))
    lam_seq <- 10 ^ seq(3, -5, length.out = 120)

    crit_num <- crit_den <- rep(0, length(lam_seq))

    # Ridge fallback should *not* use relaxation even if requested; ridge +
    # relax is not supported and is not needed for this emergency path.
    glmnet_fit_args_ridge <- glmnet_fit_args
    glmnet_fit_args_ridge$relax <- FALSE
    glmnet_fit_args_ridge$gamma <- NULL

    for (j in seq_along(lam_seq)) {
      lam <- lam_seq[j]
      for (fold in seq_len(nfolds_eff)) {
        tr_idx <- which(foldid != fold); te_idx <- which(foldid == fold)
        glmnet_fit_args_fold <- glmnet_fit_args
        glmnet_fit_args_fold$relax <- FALSE
        glmnet_fit_args_fold$gamma <- NULL

        y_tr_fit <- if (is.matrix(y_fit)) y_fit[tr_idx, , drop = FALSE] else y_fit[tr_idx]
        fit_j <- tryCatch(
          do.call(glmnet::glmnet, c(
            list(x = X_fit[tr_idx, , drop = FALSE], y = y_tr_fit,
                 alpha = 0, lambda = lam),
            glmnet_fit_args_fold
          )),
          error = function(e) NULL
        )
        if (is.null(fit_j)) next
        preds_te <- tryCatch(
          drop(stats::predict(fit_j, newx = X_fit[te_idx, , drop = FALSE],
                              s = lam, type = "response")),
          error = function(e) NULL
        )
        if (is.null(preds_te)) next

        if (is_gaussian) {
          # Gaussian: MSE criterion
          crit_num[j] <- crit_num[j] + sum((preds_te - y_glm[te_idx])^2)
          crit_den[j] <- crit_den[j] + length(te_idx)
        } else {
          # Binomial: negative log-likelihood / deviance-like criterion
          p_hat <- pmin(pmax(preds_te, 1e-8), 1 - 1e-8)
          yy <- y_glm[te_idx]
          crit_num[j] <- crit_num[j] - sum(yy * log(p_hat) + (1 - yy) * log(1 - p_hat))
          crit_den[j] <- crit_den[j] + length(te_idx)
        }
      }
    }
    crit_by_lambda <- crit_num / crit_den
    if (!any(is.finite(crit_by_lambda))) {
      return(.cv_return_intercept_only(
        fam          = fam,
        y_glm        = y_glm,
        X            = X,
        mf           = mf,
        glmnet_alpha = glmnet_alpha,
        choose_rule  = choose_rule,
        drop_msg     = drop_msg,
        nfolds_val   = nfolds_eff,
        repeats_val  = repeats,
        relaxed      = relaxed,
        note_extra   = "ridge_fallback_failed"
      ))
    }

    j_best <- which.min(crit_by_lambda)
    best_lambda <- lam_seq[j_best]

    final_ridge <- tryCatch(
      do.call(glmnet::glmnet, c(
        list(x = X_fit, y = y_fit, alpha = 0, lambda = best_lambda),
        glmnet_fit_args_ridge
      )),
      error = function(e) NULL
    )
    if (is.null(final_ridge)) {
      return(.cv_return_intercept_only(
        fam          = fam,
        y_glm        = y_glm,
        X            = X,
        mf           = mf,
        glmnet_alpha = glmnet_alpha,
        choose_rule  = choose_rule,
        drop_msg     = drop_msg,
        nfolds_val   = nfolds_eff,
        repeats_val  = repeats,
        relaxed      = relaxed,
        note_extra   = "ridge_final_fit_failed"
      ))
    }
    coef_mat   <- as.matrix(stats::coef(final_ridge, s = best_lambda))
    coef_mat   <- .svem_strip_glmnet_sentinel(coef_mat, glmnet_design$sentinel)
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_ridge, newx = X_fit, s = best_lambda, type = "response"))

    debias_fit <- if (is_gaussian && length(y_glm) >= 10 && stats::var(y_pred) > 0) stats::lm(y_glm ~ y_pred) else NULL
    y_pred_debiased <- if (!is.null(debias_fit)) stats::predict(debias_fit) else NULL

    # Build schema
    terms_clean <- terms_store
    feature_names <- colnames(X)

    xlevels <- list()
    vars_in_terms <- predictor_vars
    for (v in vars_in_terms) {
      if (v %in% names(data_frame_used)) {
        if (is.factor(data_frame_used[[v]])) {
          xlevels[[v]] <- levels(data_frame_used[[v]])
        } else if (is.character(data_frame_used[[v]])) {
          xlevels[[v]] <- sort(unique(as.character(data_frame_used[[v]])))
        }
      }
    }
    contrasts_used <- contrasts_used_full
    if (!is.null(contrasts_used)) {
      contrasts_used <- lapply(contrasts_used, function(ci) {
        if (is.character(ci)) get(ci, mode = "function") else ci
      })
    }
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
      predictor_vars = predictor_vars,
      var_classes   = var_classes,
      terms_hash    = safe_hash(terms_str)
    )

    # Diagnostics: model size from ridge fallback
    k_final  <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = TRUE)
    k_final0 <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = FALSE)
    diagnostics <- list(
      k_final              = k_final,
      k_final_no_intercept = k_final0,
      selected_gamma       = 1
    )

    res <- list(
      parms = best_coefs,
      glmnet_alpha = glmnet_alpha,
      best_alpha = 0,
      best_lambda = best_lambda,
      best_gamma = 1,
      y_pred = y_pred,
      debias_fit = debias_fit,
      y_pred_debiased = y_pred_debiased,
      cv_summary = list(),
      formula = formula_store,
      terms = terms_clean,
      training_X = X,
      actual_y = y_glm,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = c("ridge_fallback", drop_msg),
      meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                  family = fam, relaxed = relaxed, relax_cv_fallbacks = NA_integer_,
                  best_gamma = 1, cv_object = NULL),
      diagnostics = diagnostics,
      family = fam
    )
    class(res) <- c("svem_cv","svem_model")
    return(res)
  }

  # --- Pick alpha/lambda/gamma (1se or min) ---
  pick_for_alpha <- function(agg) {
    j <- if (identical(choose_rule, "min")) agg$idx_min else agg$idx_1se
    raw_score <- agg$mean_cvm[j]
    list(
      idx = j,
      score = if (isTRUE(agg$higher_is_better)) -raw_score else raw_score,
      raw_score = raw_score,
      lambda = agg$lambda[j],
      gamma = agg$gamma[j]
    )
  }
  alpha_scores <- lapply(per_alpha, pick_for_alpha)
  best_idx <- which.min(vapply(alpha_scores, function(x) x$score, numeric(1)))
  best_alpha <- as.numeric(names(alpha_scores)[best_idx])
  best_lambda <- alpha_scores[[best_idx]]$lambda
  best_gamma <- alpha_scores[[best_idx]]$gamma

  # --- Final refit & outputs ---
  note_relax <- NULL
  cv_obj_to_keep <- NULL

  if (isTRUE(relaxed)) {
    selected_agg <- per_alpha[[best_idx]]
    # Every cv.glmnet object already contains the full-data fit on its master
    # lambda path. Reuse the representative fit from the selected alpha rather
    # than performing another, statistically unrelated CV run.
    final_fit <- selected_agg$cv_object$glmnet.fit
    if (is.null(final_fit)) {
      final_fit <- tryCatch(
        do.call(glmnet::glmnet, c(
          list(
            x = X_fit,
            y = y_fit,
            alpha = best_alpha,
            lambda = selected_agg$lambda_path
          ),
          glmnet_fit_args
        )),
        error = function(e) NULL
      )
    }

    if (is.null(final_fit)) {
      warning("Final relaxed glmnet refit failed; refitting without relaxation.")
      final_nonrelaxed_args <- glmnet_fit_args
      final_nonrelaxed_args$relax <- FALSE
      final_nonrelaxed_args$gamma <- NULL
      final_fit <- do.call(glmnet::glmnet, c(
        list(
          x = X_fit,
          y = y_fit,
          alpha = best_alpha,
          lambda = selected_agg$lambda_path
        ),
        final_nonrelaxed_args
      ))
      coef_mat   <- as.matrix(stats::coef(final_fit, s = best_lambda))
      coef_mat   <- .svem_strip_glmnet_sentinel(coef_mat, glmnet_design$sentinel)
      best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
      y_pred <- drop(stats::predict(final_fit, newx = X_fit, s = best_lambda, type = "response"))
      best_gamma <- 1
      note_relax <- "final relaxed glmnet refit failed; coefs from glmnet (non-relaxed)"
    } else {
      coef_mat <- if (inherits(final_fit, "relaxed")) {
        as.matrix(stats::coef(final_fit, s = best_lambda, gamma = best_gamma))
      } else {
        best_gamma <- 1
        as.matrix(stats::coef(final_fit, s = best_lambda))
      }
      coef_mat <- .svem_strip_glmnet_sentinel(coef_mat, glmnet_design$sentinel)
      best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
      mm <- cbind("(Intercept)" = 1, X)
      y_pred <- as.numeric(mm %*% best_coefs[match(colnames(mm), names(best_coefs))])
      if (!is_gaussian && identical(fam, "binomial")) {
        # transform to probability scale for consistency with type = "response"
        y_pred <- 1 / (1 + exp(-y_pred))
      }
      note_relax <- paste0(
        "coefs from jointly selected repeated relaxed CV (gamma = ",
        format(best_gamma, trim = TRUE), ")"
      )
    }
    # Keep a representative selected-alpha cv.glmnet object only on explicit
    # request. The aggregate surface remains available in cv_summary.
    if (isTRUE(dots$keep)) cv_obj_to_keep <- selected_agg$cv_object
  } else {
    best_gamma <- 1
    final_fit <- do.call(glmnet::glmnet, c(
      list(x = X_fit, y = y_fit, alpha = best_alpha, lambda = best_lambda),
      glmnet_fit_args
    ))
    coef_mat   <- as.matrix(stats::coef(final_fit, s = best_lambda))
    coef_mat   <- .svem_strip_glmnet_sentinel(coef_mat, glmnet_design$sentinel)
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_fit, newx = X_fit, s = best_lambda, type = "response"))
  }

  # Diagnostics: model size for final selected model
  k_final  <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = TRUE)
  k_final0 <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = FALSE)

  debias_fit <- if (is_gaussian && length(y_glm) >= 10 && stats::var(y_pred) > 0) stats::lm(y_glm ~ y_pred) else NULL
  y_pred_debiased <- if (!is.null(debias_fit)) stats::predict(debias_fit) else NULL

  cv_summary <- lapply(per_alpha, function(agg) {
    data.frame(
      gamma       = agg$gamma,
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
  terms_clean <- terms_store

  feature_names <- colnames(X)
  xlevels <- list()
  vars_in_terms <- predictor_vars
  for (v in vars_in_terms) {
    if (v %in% names(data_frame_used)) {
      if (is.factor(data_frame_used[[v]])) {
        xlevels[[v]] <- levels(data_frame_used[[v]])
      } else if (is.character(data_frame_used[[v]])) {
        xlevels[[v]] <- sort(unique(as.character(data_frame_used[[v]])))
      }
    }
  }
  contrasts_used <- contrasts_used_full
  if (!is.null(contrasts_used)) {
    contrasts_used <- lapply(contrasts_used, function(ci) {
      if (is.character(ci)) get(ci, mode = "function") else ci
    })
  }
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
    predictor_vars = predictor_vars,
    var_classes   = var_classes,
    terms_hash    = safe_hash(terms_str)
  )

  diagnostics <- list(
    k_final              = k_final,
    k_final_no_intercept = k_final0,
    selected_gamma       = best_gamma
  )

  result <- list(
    parms = best_coefs,
    glmnet_alpha = glmnet_alpha,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_gamma = best_gamma,
    y_pred = y_pred,
    debias_fit = debias_fit,
    y_pred_debiased = y_pred_debiased,
    cv_summary = cv_summary,
    formula = formula_store,
    terms = terms_clean,
    training_X = X,
    actual_y = y_glm,
    xlevels = xlevels,
    contrasts = contrasts_used,
    schema = schema,
    note = c(drop_msg, note_relax),
    meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                family = fam, relaxed = relaxed, relax_cv_fallbacks = total_fallbacks,
                best_gamma = best_gamma, cv_object = cv_obj_to_keep),
    diagnostics = diagnostics,
    family = fam
  )

  class(result) <- c("svem_cv", "svem_model")

  return(result)
}
