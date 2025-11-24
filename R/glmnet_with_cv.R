#' Fit a glmnet Model with Repeated Cross-Validation
#'
#' Repeated K-fold cross-validation over a per-alpha lambda path, with a
#' combined 1-SE rule across repeats. Preserves fields expected by
#' \code{predict.svem_model()} and internal prediction helpers. Optionally uses
#' \code{glmnet}'s built-in relaxed elastic net for both the warm-start path and
#' each CV fit. When \code{relaxed = TRUE}, the final coefficients are taken
#' from a \code{cv.glmnet()} object at the chosen lambda so that the returned
#' model reflects the relaxed solution (including its chosen gamma).
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
#' @param choose_rule Character; how to choose lambda within each alpha:
#'   \itemize{
#'     \item \code{"min"}: lambda minimizing the cross-validated criterion.
#'     \item \code{"1se"}: largest lambda within 1 combined SE of the minimum,
#'           where the SE includes both within- and between-repeat variability.
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
#'   \code{gamma} path (default \code{FALSE}). If
#'   \code{cv.glmnet(relax = TRUE)} fails for a particular repeat/alpha, the
#'   function retries that fit without relaxation; the number of such
#'   fallbacks is recorded in \code{meta$relax_cv_fallbacks}.
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
#'   \code{glmnet()}, for example: \code{weights}, \code{parallel},
#'   \code{type.measure}, \code{intercept}, \code{maxit},
#'   \code{lower.limits}, \code{upper.limits}, \code{penalty.factor},
#'   \code{offset}, \code{standardize.response}, \code{keep}, and so on. If
#'   \code{family} is supplied here, it is ignored in favor of the explicit
#'   \code{family} argument.
#'
#' @details
#' The basic workflow is:
#' \enumerate{
#'   \item For each \code{alpha} in \code{glmnet_alpha}, generate a set of CV
#'         fold IDs (shared across alphas and repeats).
#'   \item For that alpha, run \code{repeats} independent \code{cv.glmnet()}
#'         fits, align the lambda paths, and aggregate the CV curves.
#'   \item At each lambda, compute a combined SE that accounts for both
#'         within-repeat and between-repeat variability.
#'   \item Apply \code{choose_rule} (\code{"min"} or \code{"1se"}) to select
#'         lambda for that alpha, then choose the best alpha by comparing these
#'         per-alpha scores.
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
#'   \item The training \code{terms} are stored with environment set to
#'         \code{baseenv()}.
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
#'   \item \code{y_pred} In-sample predictions from the returned coefficients
#'         (fitted values for Gaussian; probabilities for binomial).
#'   \item \code{debias_fit} For Gaussian, an optional \code{lm(y ~ y_pred)}
#'         calibration model; \code{NULL} otherwise.
#'   \item \code{y_pred_debiased} If \code{debias_fit} exists, its fitted
#'         values; otherwise \code{NULL}.
#'   \item \code{cv_summary} Named list (one element per alpha) of data frames
#'         with columns \code{lambda}, \code{mean_cvm}, \code{sd_cvm},
#'         \code{se_combined}, \code{n_repeats}, \code{idx_min},
#'         \code{idx_1se}.
#'   \item \code{formula} Original modeling formula.
#'   \item \code{terms} Training \code{terms} object with environment set to
#'         \code{baseenv()}.
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
#' n <- 100; p <- 10
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
#' n2 <- 150; p2 <- 8
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
#'   glmnet_alpha = c(0.5, 1),
#'   nfolds = 5,
#'   repeats = 2,
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

  # drop any family accidentally passed via ...
  if ("family" %in% names(dots)) {
    dots$family <- NULL
  }

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

  # --- Helper for intercept-only returns (used in both places) ---
  .cv_return_intercept_only <- function(
    fam, y_glm, X, mf, glmnet_alpha, choose_rule,
    drop_msg, nfolds_val, repeats_val, relaxed, w, note_extra = NULL
  ) {
    n_local <- nrow(X)

    if (identical(fam, "gaussian")) {
      intercept <- if (!is.null(w)) sum(w * y_glm) / sum(w) else mean(y_glm)
      y_pred <- rep(intercept, n_local)
    } else {  # binomial: intercept on logit scale, y_pred as probabilities
      p_hat <- if (!is.null(w)) sum(w * y_glm) / sum(w) else mean(y_glm)
      p_hat <- pmin(pmax(p_hat, 1e-6), 1 - 1e-6)
      intercept <- stats::qlogis(p_hat)
      y_pred <- rep(p_hat, n_local)
    }

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

    note_vec <- c("intercept_only", drop_msg, note_extra)
    note_vec <- note_vec[!is.na(note_vec)]

    # Diagnostics: model size (intercept only)
    k_final  <- 1L
    k_final0 <- 0L
    diagnostics <- list(
      k_final              = k_final,
      k_final_no_intercept = k_final0
    )

    res <- list(
      parms = c("(Intercept)" = intercept),
      glmnet_alpha = glmnet_alpha,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      y_pred = y_pred,
      debias_fit = NULL,
      y_pred_debiased = NULL,
      cv_summary = list(),
      formula = formula,
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
      w           = w,
      note_extra  = NULL
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
    family       = fam,
    type.measure = if (identical(fam, "binomial")) "deviance" else "mse",
    maxit        = 1e6,
    relax        = isTRUE(relaxed)
  )
  if (isTRUE(relaxed) && !is.null(relax_gamma)) base_cv_args$gamma <- relax_gamma
  cv_base_args <- utils::modifyList(base_cv_args, dots, keep.null = TRUE)

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
                c(list(x = X, y = y_glm,
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
                      c(list(x = X, y = y_glm,
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
      w           = w,
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

    for (j in seq_along(lam_seq)) {
      lam <- lam_seq[j]
      for (fold in seq_len(nfolds_eff)) {
        tr_idx <- which(foldid != fold); te_idx <- which(foldid == fold)
        glmnet_fit_args_fold <- glmnet_fit_args
        if (!is.null(w)) glmnet_fit_args_fold$weights <- w[tr_idx]
        if (!is.null(off)) glmnet_fit_args_fold$offset <- off[tr_idx]
        fit_j <- tryCatch(
          do.call(glmnet::glmnet, c(
            list(x = X[tr_idx, , drop = FALSE], y = y_glm[tr_idx],
                 alpha = 0, lambda = lam),
            glmnet_fit_args_fold
          )),
          error = function(e) NULL
        )
        if (is.null(fit_j)) next
        preds_te <- drop(stats::predict(fit_j, newx = X[te_idx, , drop = FALSE],
                                        s = lam, type = "response"))

        if (is_gaussian) {
          # Gaussian: MSE criterion
          if (is.null(w)) {
            crit_num[j] <- crit_num[j] + sum((preds_te - y_glm[te_idx])^2)
            crit_den[j] <- crit_den[j] + length(te_idx)
          } else {
            w_te <- w[te_idx]
            crit_num[j] <- crit_num[j] + sum(w_te * (preds_te - y_glm[te_idx])^2)
            crit_den[j] <- crit_den[j] + sum(w_te)
          }
        } else {
          # Binomial: negative log-likelihood / deviance-like criterion
          p_hat <- pmin(pmax(preds_te, 1e-8), 1 - 1e-8)
          yy <- y_glm[te_idx]
          if (is.null(w)) {
            crit_num[j] <- crit_num[j] - sum(yy * log(p_hat) + (1 - yy) * log(1 - p_hat))
            crit_den[j] <- crit_den[j] + length(te_idx)
          } else {
            w_te <- w[te_idx]
            crit_num[j] <- crit_num[j] - sum(w_te * (yy * log(p_hat) + (1 - yy) * log(1 - p_hat)))
            crit_den[j] <- crit_den[j] + sum(w_te)
          }
        }
      }
    }
    crit_by_lambda <- crit_num / crit_den
    j_best <- which.min(crit_by_lambda)
    best_lambda <- lam_seq[j_best]

    final_ridge <- do.call(glmnet::glmnet, c(
      list(x = X, y = y_glm, alpha = 0, lambda = best_lambda),
      glmnet_fit_args
    ))
    coef_mat   <- as.matrix(stats::coef(final_ridge, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_ridge, newx = X, s = best_lambda, type = "response"))

    debias_fit <- if (is_gaussian && length(y_glm) >= 10 && stats::var(y_pred) > 0) stats::lm(y_glm ~ y_pred) else NULL
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
      terms_hash    = safe_hash(terms_str)
    )

    # Diagnostics: model size from ridge fallback
    k_final  <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = TRUE)
    k_final0 <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = FALSE)
    diagnostics <- list(
      k_final              = k_final,
      k_final_no_intercept = k_final0
    )

    res <- list(
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
      actual_y = y_glm,
      xlevels = xlevels,
      contrasts = contrasts_used,
      schema = schema,
      note = c("ridge_fallback", drop_msg),
      meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                  family = fam, relaxed = relaxed, relax_cv_fallbacks = NA_integer_, cv_object = NULL),
      diagnostics = diagnostics,
      family = fam
    )
    class(res) <- c("svem_cv","svem_model")
    return(res)
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
                      c(list(x = X, y = y_glm,
                             alpha  = best_alpha,
                             foldid = foldids[[1]],
                             exclude = exclude),
                        drop_reserved(cv_base_args, reserved_cv)))
    coef_mat   <- as.matrix(stats::coef(cv_one, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    mm <- cbind("(Intercept)" = 1, X)
    y_pred <- as.numeric(mm %*% best_coefs[match(colnames(mm), names(best_coefs))])
    if (!is_gaussian && identical(fam, "binomial")) {
      # transform to probability scale for consistency with type = "response"
      y_pred <- 1 / (1 + exp(-y_pred))
    }
    note_relax <- "coefs from cv.glmnet (relaxed)"
    if (isTRUE(is.null(dots$keep)) || isTRUE(dots$keep)) cv_obj_to_keep <- cv_one
  } else {
    final_fit <- do.call(glmnet::glmnet, c(
      list(x = X, y = y_glm, alpha = best_alpha, lambda = best_lambda),
      glmnet_fit_args
    ))
    coef_mat   <- as.matrix(stats::coef(final_fit, s = best_lambda))
    best_coefs <- drop(coef_mat); names(best_coefs) <- rownames(coef_mat)
    y_pred <- drop(stats::predict(final_fit, newx = X, s = best_lambda, type = "response"))
  }

  # Diagnostics: model size for final selected model
  k_final  <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = TRUE)
  k_final0 <- .support_size_one(best_coefs, base_tol = 1e-7, count_intercept = FALSE)

  debias_fit <- if (is_gaussian && length(y_glm) >= 10 && stats::var(y_pred) > 0) stats::lm(y_glm ~ y_pred) else NULL
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
    terms_hash    = safe_hash(terms_str)
  )

  diagnostics <- list(
    k_final              = k_final,
    k_final_no_intercept = k_final0
  )

  result <- list(
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
    actual_y = y_glm,
    xlevels = xlevels,
    contrasts = contrasts_used,
    schema = schema,
    note = c(drop_msg, note_relax),
    meta = list(nfolds = nfolds_eff, repeats = repeats, rule = choose_rule,
                family = fam, relaxed = relaxed, relax_cv_fallbacks = total_fallbacks,
                cv_object = cv_obj_to_keep),
    diagnostics = diagnostics,
    family = fam
  )

  class(result) <- c("svem_cv", "svem_model")

  return(result)
}
