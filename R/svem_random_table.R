#' Generate a Random Prediction Table from a Fitted SVEMnet Model (no refit)
#'
#' Samples the original predictor factor space using ranges and levels cached
#' inside a fitted \code{svem_model} (from \code{SVEMnet()}) and computes
#' predictions at those points. Optional mixture groups let you sample
#' composition variables on a (possibly truncated) simplex via a Dirichlet draw.
#' No refitting is performed.
#'
#' @param object A fitted \code{svem_model} returned by \code{SVEMnet()}.
#'               The object must contain \code{$sampling_schema}, which is
#'               created by the updated \code{SVEMnet()}.
#' @param n Number of random points to generate (default \code{1000}).
#' @param mixture_groups Optional list of mixture factor groups. Each element is
#'   a list with components:
#'   \itemize{
#'     \item \code{vars}: character vector of mixture variable names (must be in the model).
#'     \item \code{lower}: numeric vector of lower bounds (same length as \code{vars}).
#'     \item \code{upper}: numeric vector of upper bounds (same length as \code{vars}).
#'     \item \code{total}: scalar specifying the group total (e.g., 1.0).
#'   }
#'   If omitted, all variables are sampled independently using the cached schema.
#' @param debias Logical; if \code{TRUE}, applies the calibration from
#'   \code{object$debias_fit} (if available) to the mean prediction (default \code{FALSE}).
#'
#' @return A data frame with sampled predictors and a predicted response column.
#'   The response column name matches the left-hand side of \code{object$formula}.
#'
#' @details
#' This function uses:
#' \itemize{
#'   \item \code{object$sampling_schema$num_ranges} for uniform sampling of numeric variables.
#'   \item \code{object$sampling_schema$factor_levels} for categorical sampling.
#'   \item \code{object$terms}, \code{object$xlevels}, and \code{object$contrasts}
#'         (via \code{predict.svem_model}) to encode the model matrix consistently.
#' }
#' Mixture groups are handled by drawing Dirichlet weights and mapping them to the
#' truncated simplex defined by \code{lower}, \code{upper}, and \code{total}; proposals
#' violating upper bounds are rejected (with oversampling to keep it efficient).
#'
#' @seealso \code{\link{SVEMnet}}, \code{\link{predict.svem_model}},
#'   \code{\link{svem_significance_test}}, \code{\link{svem_significance_test_parallel}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 40
#' X1 <- runif(n); X2 <- runif(n)
#' A <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)  # simple mixture-ish
#' F <- factor(sample(c("lo","hi"), n, TRUE))
#' y <- 1 + 2*X1 - X2 + 3*A + 1.5*B + 0.5*C + (F=="hi") + rnorm(n, 0, 0.3)
#' d <- data.frame(y, X1, X2, A, B, C, F)
#'
#' fit <- SVEMnet(y ~ X1 + X2 + A + B + C + F, d, nBoot = 50, glmnet_alpha = 1)
#'
#' # No mixture constraints (independent sampling using cached ranges/levels)
#' tbl1 <- svem_random_table_from_model(fit, n = 50)
#' head(tbl1)
#'
#' # With mixture constraints for A,B,C that sum to 1 and have bounds
#' mix <- list(list(vars = c("A","B","C"),
#'                  lower = c(0.1, 0.1, 0.1),
#'                  upper = c(0.7, 0.7, 0.7),
#'                  total = 1.0))
#' tbl2 <- svem_random_table_from_model(fit, n = 50, mixture_groups = mix)
#' head(tbl2)
#' }
#'
#' @importFrom lhs maximinLHS
#' @importFrom stats rgamma
#' @export
svem_random_table_from_model <- function(object, n = 1000, mixture_groups = NULL, debias = FALSE) {
  if (!inherits(object, "svem_model")) stop("object must be a 'svem_model' from SVEMnet().")
  if (is.null(object$sampling_schema) || !is.list(object$sampling_schema)) {
    stop("object$sampling_schema is missing. Refit with the updated SVEMnet() that records sampling_schema.")
  }
  ss <- object$sampling_schema
  predictor_vars <- ss$predictor_vars
  var_classes    <- ss$var_classes
  num_ranges     <- ss$num_ranges           # 2 x K matrix, rownames c('min','max')
  factor_levels  <- ss$factor_levels        # named list

  if (!length(predictor_vars)) stop("No predictor variables found in sampling_schema.")

  # Partition variables
  is_num <- names(var_classes)[var_classes %in% c("numeric","integer")]
  is_cat <- setdiff(predictor_vars, is_num)

  # Helper: truncated-simplex Dirichlet sampler (with upper-bound rejection)
  .sample_trunc_dirichlet <- function(n, lower, upper, total,
                                      alpha = NULL, oversample = 4L, max_tries = 10000L) {
    k <- length(lower)
    if (length(upper) != k) stop("upper must have same length as lower.")
    if (is.null(alpha)) alpha <- rep(1, k)

    min_sum <- sum(lower); max_sum <- sum(upper)
    if (total < min_sum - 1e-12 || total > max_sum + 1e-12) {
      stop("Infeasible mixture constraints: need sum(lower) <= total <= sum(upper).")
    }

    avail <- total - min_sum
    if (avail <= 1e-12) return(matrix(rep(lower, each = n), nrow = n))

    res <- matrix(NA_real_, nrow = n, ncol = k)
    filled <- 0L; tries <- 0L
    while (filled < n && tries < max_tries) {
      m <- max(oversample * (n - filled), 1L)
      g <- matrix(stats::rgamma(m * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
      W <- g / rowSums(g)
      cand <- matrix(lower, nrow = m, ncol = k, byrow = TRUE) + avail * W
      ok <- cand <= matrix(upper, nrow = m, ncol = k, byrow = TRUE)
      ok <- rowSums(ok) == k
      if (any(ok)) {
        keep <- which(ok)
        take <- min(length(keep), n - filled)
        res[(filled + 1):(filled + take), ] <- cand[keep[seq_len(take)], , drop = FALSE]
        filled <- filled + take
      }
      tries <- tries + 1L
    }
    if (filled < n) {
      stop("Could not sample enough feasible mixture points within max_tries. ",
           "Relax bounds or increase oversample/max_tries.")
    }
    res
  }

  # Track any mixture variables (to avoid double-sampling them as 'numeric')
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    # Validate groups and collect names
    for (grp in mixture_groups) {
      if (is.null(grp$vars)) stop("Each mixture group must contain a 'vars' character vector.")
      if (!all(grp$vars %in% predictor_vars)) {
        missing <- setdiff(grp$vars, predictor_vars)
        stop("Mixture variables not in model predictors: ", paste(missing, collapse = ", "))
      }
      mixture_vars <- c(mixture_vars, grp$vars)
    }
    # No variable may appear in more than one group
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }

  # ----- Sample non-mixture numeric variables (uniform over cached ranges) -----
  nonmix_num <- setdiff(is_num, mixture_vars)
  T_num <- NULL
  if (length(nonmix_num)) {
    # If any numeric var has a degenerate or missing range, fall back to [0,1]
    rng <- vapply(nonmix_num, function(v) {
      if (!is.null(num_ranges) && ncol(num_ranges) && v %in% colnames(num_ranges)) {
        r <- num_ranges[, v]
        if (!all(is.finite(r)) || r[1] >= r[2]) c(0, 1) else r
      } else {
        c(0, 1)
      }
    }, numeric(2))
    # Latin hypercube in [0,1] then map to ranges
    U <- as.matrix(lhs::maximinLHS(n, length(nonmix_num)))
    T_num <- matrix(NA_real_, nrow = n, ncol = length(nonmix_num))
    colnames(T_num) <- nonmix_num
    for (j in seq_along(nonmix_num)) {
      lo <- rng[1, j]; hi <- rng[2, j]
      T_num[, j] <- lo + (hi - lo) * U[, j]
    }
    T_num <- as.data.frame(T_num)
  }

  # ----- Sample mixture groups (Dirichlet on truncated simplex) -----
  T_mix <- NULL
  if (!is.null(mixture_groups) && length(mixture_groups)) {
    all_mix_vars <- unlist(lapply(mixture_groups, `[[`, "vars"))
    T_mix <- matrix(NA_real_, nrow = n, ncol = length(all_mix_vars))
    colnames(T_mix) <- all_mix_vars

    for (grp in mixture_groups) {
      vars  <- grp$vars
      k     <- length(vars)
      lower <- if (!is.null(grp$lower)) grp$lower else rep(0, k)
      upper <- if (!is.null(grp$upper)) grp$upper else rep(1, k)
      total <- if (!is.null(grp$total)) grp$total else 1

      if (length(lower) != k || length(upper) != k) {
        stop("lower and upper must each have length equal to the number of mixture variables (",
             paste(vars, collapse = ", "), ").")
      }
      vals <- .sample_trunc_dirichlet(n, lower, upper, total)
      colnames(vals) <- vars
      T_mix[, vars] <- vals
    }
    T_mix <- as.data.frame(T_mix)
  }

  # ----- Sample categorical variables (use cached level sets) -----
  T_cat <- NULL
  if (length(is_cat)) {
    T_cat <- vector("list", length(is_cat))
    names(T_cat) <- is_cat
    for (v in is_cat) {
      lev <- factor_levels[[v]]
      # If no cached levels (shouldn't happen with updated SVEMnet), fallback to two dummy levels
      if (is.null(lev) || !length(lev)) lev <- c("L1", "L2")
      T_cat[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
    }
    T_cat <- as.data.frame(T_cat, stringsAsFactors = FALSE)
  }

  # ----- Combine into a single raw data frame aligned to predictor names -----
  parts <- list(T_num, T_mix, T_cat)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) stop("No predictors could be sampled from the schema.")
  T_data <- do.call(cbind, parts)

  # Ensure all predictors present (fill missing with sensible defaults)
  missing_pred <- setdiff(predictor_vars, colnames(T_data))
  if (length(missing_pred)) {
    for (v in missing_pred) {
      if (v %in% names(var_classes) && var_classes[[v]] %in% c("numeric","integer")) {
        # fallback numeric in [0,1]
        T_data[[v]] <- runif(n)
      } else {
        # fallback factor with a single dummy level
        T_data[[v]] <- factor(rep("L1", n))
      }
    }
  }

  # Order columns like predictor_vars (not strictly necessary, but tidy)
  T_data <- T_data[, predictor_vars, drop = FALSE]

  # ----- Predict using the fitted model -----
  preds <- predict(object, newdata = T_data, debias = debias)
  if (is.list(preds) && !is.null(preds$fit)) preds <- preds$fit

  # Name response from stored formula
  resp_name <- tryCatch(as.character(object$formula[[2]]), error = function(e) "yhat")
  out <- T_data
  out[[resp_name]] <- as.numeric(preds)
  rownames(out) <- NULL
  out
}
