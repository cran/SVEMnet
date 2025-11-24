#' Member-level predictions for an SVEM model
#'
#' @description
#' Internal utility to compute per-bootstrap predictions for an
#' \code{svem_model} on new data. It rebuilds the design matrix in the same
#' way as \code{predict.svem_model()}, aligns columns to \code{training_X},
#' and applies \code{coef_matrix} to obtain an \eqn{n x B} matrix of
#' predictions.
#'
#' @param object An \code{svem_model} object as returned by \code{SVEMnet()},
#'   containing at least \code{terms}, \code{training_X}, \code{xlevels},
#'   \code{contrasts} (optional), \code{coef_matrix}, and \code{family}.
#' @param newdata A \code{data.frame} with the predictors required by
#'   \code{object}.
#' @param type Character string, one of \code{"response"} (default) or
#'   \code{"link"}. For binomial models, \code{"response"} returns
#'   probabilities via \code{plogis()} and \code{"link"} returns log-odds.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{pred_mat}}{Numeric matrix of dimension \eqn{n x B} containing
#'     member-level predictions on the requested scale.}
#'   \item{\code{bad_rows}}{Logical vector of length \eqn{n} indicating rows
#'     with missing / unseen levels; those rows are set to \code{NA} in
#'     \code{pred_mat}.}
#' }
#'
#' @keywords internal
#' @seealso \code{\link{SVEMnet}}, \code{\link{predict.svem_model}}
#' @noRd
.svem_member_predictions <- function(object,
                                     newdata,
                                     type = c("response", "link")) {
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.")
  }

  type <- match.arg(type)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Basic checks
  if (is.null(object$terms) || is.null(object$training_X) ||
      is.null(object$coef_matrix)) {
    stop("svem_model object must contain 'terms', 'training_X', and 'coef_matrix'.")
  }

  fam <- tolower(object$family %||% "gaussian")

  # Rebuild design matrix (mirror predict.svem_model)
  # Always use baseenv() to avoid depending on any stale environments
  # captured when the model was fitted.
  terms_obj <- stats::delete.response(object$terms)
  environment(terms_obj) <- baseenv()

  # Harmonize factor / character predictors to training levels
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

  mf <- stats::model.frame(terms_obj, data = newdata, na.action = stats::na.pass)

  # Use stored contrasts or saved global contrast options
  ctr <- object$contrasts
  have_ctr <- !is.null(ctr)

  if (!have_ctr) {
    old_opts <- options("contrasts")
    on.exit(options(old_opts), add = TRUE)

    fit_opts <- tryCatch(object$schema$contrasts_options,
                         error = function(e) NULL)
    fit_opts <- fit_opts %||% old_opts$contrasts
    if (!is.null(fit_opts)) {
      options(contrasts = fit_opts)
    }
    mm <- stats::model.matrix(terms_obj, data = mf)
  } else {
    # More robust reconstruction of contrasts from stored spec
    if (!is.list(ctr)) {
      stop("Invalid contrasts specification in fitted svem_model object.",
           call. = FALSE)
    }
    if (length(ctr)) {
      for (nm in names(ctr)) {
        val <- ctr[[nm]]
        if (is.character(val)) {
          # skip empty / NA names: use default contrast for this variable
          if (!length(val) || is.na(val)) {
            ctr[[nm]] <- NULL
            next
          }
          if (!exists(val, mode = "function")) {
            warning(
              "Contrast function '", val,
              "' not found; using default contrasts for variable '", nm, "'."
            )
            ctr[[nm]] <- NULL
          } else {
            ctr[[nm]] <- get(val, mode = "function")
          }
        }
      }
    }
    mm <- stats::model.matrix(terms_obj, data = mf, contrasts.arg = ctr)
  }

  # Drop intercept column: glmnet handled the intercept separately
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col)) {
    mm <- mm[, -int_col, drop = FALSE]
  }

  # Identify "bad" rows (non-finite entries in the model matrix)
  bad_rows <- !stats::complete.cases(mm) |
    rowSums(!is.finite(mm), na.rm = TRUE) > 0L
  if (any(bad_rows)) {
    mm[!is.finite(mm)] <- 0
  }

  # Align columns with training_X; zero-fill missing columns
  train_cols <- colnames(object$training_X)
  if (is.null(train_cols) || !length(train_cols)) {
    stop("`object$training_X` must have column names.")
  }

  m <- nrow(mm)

  mm_use <- matrix(0, nrow = m, ncol = length(train_cols))
  colnames(mm_use) <- train_cols

  common_cols <- intersect(colnames(mm), train_cols)
  if (length(common_cols)) {
    mm_use[, common_cols] <- mm[, common_cols, drop = FALSE]
  }

  storage.mode(mm_use) <- "double"

  # Member predictions from coef_matrix
  coef_matrix <- object$coef_matrix  # nBoot x (p + 1)
  if (!is.matrix(coef_matrix) || ncol(coef_matrix) < 1L) {
    stop("`object$coef_matrix` must be a numeric matrix with intercept and slopes.")
  }

  n_boot <- nrow(coef_matrix)

  intercepts <- coef_matrix[, 1]                # length B
  betas      <- coef_matrix[, -1, drop = FALSE] # B x p

  # Ensure coefficient columns align with training_X order
  if (!is.null(colnames(betas))) {
    missing_cols <- setdiff(train_cols, colnames(betas))
    if (length(missing_cols)) {
      stop(
        "`object$coef_matrix` is missing coefficients for training design columns: ",
        paste(missing_cols, collapse = ", "),
        ". This usually means the stored object is out of sync with the ",
        "training design (e.g. formula/specification or contrasts changed).",
        call. = FALSE
      )
    }
    betas <- betas[, train_cols, drop = FALSE]
  } else if (ncol(betas) != length(train_cols)) {
    stop(
      "Number of slope columns in `object$coef_matrix` (excluding intercept) ",
      "does not match the number of columns in `object$training_X` (",
      ncol(betas), " vs ", length(train_cols), ").\n",
      "This usually means the formula, deterministic expansion, contrasts, ",
      "or factor levels changed between fit and prediction.",
      call. = FALSE
    )
  }

  # Compute eta matrix: n x B
  eta_mat <- mm_use %*% t(betas) +
    matrix(intercepts, nrow = m, ncol = n_boot, byrow = TRUE)

  # Transform to requested scale
  if (identical(fam, "binomial")) {
    if (type == "response") {
      pred_mat <- stats::plogis(eta_mat)
    } else {
      pred_mat <- eta_mat
    }
  } else {
    # Gaussian (and other identity-link families): response equals link
    pred_mat <- eta_mat
  }

  # Set bad rows to NA in the returned matrix
  if (any(bad_rows)) {
    pred_mat[bad_rows, ] <- NA_real_
  }

  list(
    pred_mat = pred_mat,
    bad_rows = bad_rows
  )
}


#' Append mean-level "in spec" probabilities to an SVEM score table
#'
#' @description
#' Given an optimization score table from \code{svem_score_random()},
#' this helper recomputes bootstrap member
#' predictions (with \code{debias = FALSE}) for each response and appends:
#' \itemize{
#'   \item Per-response probabilities that the model-based mean lies within
#'         user-specified limits (\code{<resp>_p_in_spec_mean}).
#'   \item Per-response indicators that the point prediction is inside spec
#'         (\code{<resp>_in_spec_point}).
#'   \item A joint mean-level probability (\code{p_joint_mean}) and a joint
#'         point indicator (\code{joint_in_spec_point}).
#' }
#'
#' These are probabilities for the fitted conditional means under the SVEM
#' bootstrap ensemble. They are \emph{not} predictive probabilities for
#' individual observations and should not be interpreted as a formal QbD
#' design space.
#'
#' @param score_table A \code{data.frame} produced by
#'   \code{svem_score_random()}, containing
#'   the sampled predictors and per-response point predictions in columns
#'   \code{<resp>_pred} for each response \code{resp} in \code{names(objects)}.
#' @param objects Named list of fitted \code{svem_model} objects used in the
#'   optimization. Each must have a non-null \code{coef_matrix}.
#' @param specs Named list of specification lists, one per response in
#'   \code{objects}. Each entry is either \code{NULL} (no spec) or a list
#'   with numeric components \code{lower} and/or \code{upper}. If only one
#'   bound is provided, the other defaults to \code{-Inf} or \code{Inf}.
#'   Responses with no active bounds are excluded from the joint
#'   calculations.
#'
#' @return
#' A copy of \code{score_table} with additional columns:
#' \itemize{
#'   \item \code{<resp>_p_in_spec_mean}: estimated probability that the
#'         model-based mean for response \code{resp} lies within its spec.
#'   \item \code{<resp>_in_spec_point}: indicator (0/1) that the point
#'         prediction (from \code{<resp>_pred}) lies within its spec.
#'   \item \code{p_joint_mean}: product of the per-response mean-level
#'         probabilities over responses with active specs (or \code{NA} if
#'         none).
#'   \item \code{joint_in_spec_point}: indicator (0/1/NA) that all point
#'         predictions are in spec across responses with active specs.
#' }
#'
#' @noRd
svem_append_design_space_cols <- function(score_table,
                                          objects,
                                          specs) {
  # Basic checks
  if (!is.data.frame(score_table)) {
    stop("`score_table` must be a data.frame.")
  }
  if (!is.list(objects) || is.null(names(objects)) ||
      any(names(objects) == "")) {
    stop("`objects` must be a named list of svem_model objects.")
  }
  if (!is.null(specs) && !is.list(specs)) {
    stop("`specs` must be NULL or a named list.")
  }
  if (is.list(specs)) {
    if (is.null(names(specs)) || any(names(specs) == "")) {
      stop("`specs` must be a named list when provided.")
    }
    bad_names <- setdiff(names(specs), names(objects))
    if (length(bad_names)) {
      stop("The following names in `specs` do not match any response in `objects`: ",
           paste(bad_names, collapse = ", "))
    }
  } else {
    specs <- list()
  }

  # Ensure coef_matrix exists for all models
  for (nm in names(objects)) {
    obj <- objects[[nm]]
    if (!inherits(obj, "svem_model")) {
      stop("Element `", nm, "` of `objects` is not an svem_model.")
    }
    if (is.null(obj$coef_matrix)) {
      stop("Element `", nm, "` of `objects` has no coef_matrix; ",
           "refit with nBoot >= 2 to enable mean-level probability calculations.")
    }
  }

  n_row <- nrow(score_table)
  if (n_row == 0L) return(score_table)

  # Extract predictor columns from score_table via sampling_schema
  first_obj <- objects[[1L]]
  if (is.null(first_obj$sampling_schema) ||
      is.null(first_obj$sampling_schema$predictor_vars)) {
    stop("`objects` must contain models fitted with a valid sampling_schema. ",
         "Refit with an updated SVEMnet() that records sampling_schema.")
  }

  pred_vars <- first_obj$sampling_schema$predictor_vars
  pred_vars <- pred_vars[!is.na(pred_vars) & nzchar(pred_vars)]

  missing_pred_cols <- setdiff(pred_vars, colnames(score_table))
  if (length(missing_pred_cols)) {
    stop("The following predictor columns are missing from `score_table`: ",
         paste(missing_pred_cols, collapse = ", "))
  }

  X <- score_table[, pred_vars, drop = FALSE]

  # Ensure point prediction columns <resp>_pred are present for each response
  pred_cols_expected <- paste0(names(objects), "_pred")
  missing_pred <- setdiff(pred_cols_expected, colnames(score_table))
  if (length(missing_pred)) {
    stop(
      "The following prediction columns are missing from `score_table`: ",
      paste(missing_pred, collapse = ", "),
      ". Expected columns named <resp>_pred for each response."
    )
  }

  # Per-response mean-level probabilities and indicators
  p_cols   <- list()
  ind_cols <- list()

  for (resp in names(objects)) {
    obj <- objects[[resp]]
    sp  <- specs[[resp]]

    # No specs provided: skip
    if (is.null(sp)) {
      next
    }
    if (!is.list(sp)) {
      stop("Specs for response '", resp, "' must be NULL or a list with components ",
           "`lower` and/or `upper`.")
    }

    # Treat NA the same as NULL for bounds
    lo_raw <- sp$lower
    hi_raw <- sp$upper

    lo_is_missing <- is.null(lo_raw) || (length(lo_raw) == 1L && is.na(lo_raw))
    hi_is_missing <- is.null(hi_raw) || (length(hi_raw) == 1L && is.na(hi_raw))

    # Both bounds missing: treat as no active spec
    if (lo_is_missing && hi_is_missing) {
      next
    }

    lo <- if (lo_is_missing) -Inf else as.numeric(lo_raw)
    hi <- if (hi_is_missing)  Inf else as.numeric(hi_raw)

    if (!is.numeric(lo) || length(lo) != 1L ||
        !is.numeric(hi) || length(hi) != 1L) {
      stop("Specs for response '", resp, "' must have numeric `lower` and `upper` ",
           "(or NULL/NA to indicate one-sided or no spec).")
    }

    # Member predictions on response scale (no debias), with try-safety
    memb <- try(
      .svem_member_predictions(
        object  = obj,
        newdata = X,
        type    = "response"
      ),
      silent = TRUE
    )

    if (inherits(memb, "try-error") || !is.list(memb) || is.null(memb$pred_mat)) {
      if (isTRUE(getOption("svem.verbose", TRUE))) {
        msg_txt <- if (inherits(memb, "try-error")) as.character(memb)[1] else "invalid return from .svem_member_predictions()"
        message("Mean-level probability append: member predictions failed for response '",
                resp, "'. ", msg_txt)
      }
      p_cols[[resp]]   <- rep(NA_real_,    n_row)
      ind_cols[[resp]] <- rep(NA_integer_, n_row)
      next
    }

    pred_mat <- memb$pred_mat  # n x B

    if (!is.matrix(pred_mat) || nrow(pred_mat) != n_row) {
      if (isTRUE(getOption("svem.verbose", TRUE))) {
        message("Mean-level probability append: pred_mat has unexpected dimension for response '",
                resp, "'. Filling NA.")
      }
      p_cols[[resp]]   <- rep(NA_real_,    n_row)
      ind_cols[[resp]] <- rep(NA_integer_, n_row)
      next
    }

    # Probability that the fitted mean is in spec
    inside_mat <- (pred_mat >= lo & pred_mat <= hi)
    p_in_spec  <- rowMeans(inside_mat, na.rm = TRUE)
    p_in_spec[!is.finite(p_in_spec)] <- NA_real_

    # Point prediction in spec? (use <resp>_pred from score_table)
    pred_col <- paste0(resp, "_pred")
    mu_hat   <- score_table[[pred_col]]
    in_spec_point <- as.integer(mu_hat >= lo & mu_hat <= hi)
    in_spec_point[!is.finite(mu_hat)] <- NA_integer_

    p_cols[[resp]]   <- p_in_spec
    ind_cols[[resp]] <- in_spec_point
  }

  spec_resps <- names(p_cols)

  # Joint probabilities and indicators
  if (length(spec_resps) == 0L) {
    # No active specs: joint quantities undefined
    p_joint_mean <- rep(NA_real_,    n_row)
    joint_ind    <- rep(NA_integer_, n_row)
  } else {
    # Product of mean-level probabilities over responses with active specs
    p_joint_mean <- rep(1, n_row)
    for (resp in spec_resps) {
      pj <- p_cols[[resp]]
      na_idx <- is.na(pj)
      p_joint_mean <- p_joint_mean * pj
      p_joint_mean[na_idx] <- NA_real_
    }

    # Point indicator: 0 dominates NA; NA dominates 1
    ind_mat <- do.call(cbind, ind_cols[spec_resps])  # n_row x R

    zero_any <- rowSums(ind_mat == 0L, na.rm = TRUE) > 0L
    na_any   <- apply(ind_mat, 1L, function(z) any(is.na(z)))

    joint_ind <- integer(n_row)
    joint_ind[zero_any]             <- 0L
    joint_ind[!zero_any & na_any]   <- NA_integer_
    joint_ind[!zero_any & !na_any]  <- 1L
  }

  # Assemble appended columns
  if (length(spec_resps)) {
    p_df   <- as.data.frame(p_cols,   check.names = FALSE, optional = FALSE)
    ind_df <- as.data.frame(ind_cols, check.names = FALSE, optional = FALSE)

    colnames(p_df)   <- paste0(spec_resps, "_p_in_spec_mean")
    colnames(ind_df) <- paste0(spec_resps, "_in_spec_point")
  } else {
    p_df   <- data.frame()
    ind_df <- data.frame()
  }

  joint_df <- data.frame(
    p_joint_mean        = p_joint_mean,
    joint_in_spec_point = joint_ind
  )

  cbind(
    score_table,
    p_df,
    ind_df,
    joint_df
  )
}
