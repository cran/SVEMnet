#' Generate a Random Prediction Table from Multiple SVEMnet Models (no refit)
#'
#' Samples the original predictor factor space cached in fitted \code{svem_model}
#' objects and computes predictions from each model at the same random points.
#' This is intended for multiple responses built over the same factor space and
#' a deterministic factor expansion (for example via a shared
#' \code{\link{bigexp_terms}}), so that a shared sampling schema is available.
#'
#' No refitting is performed. Predictions are obtained by
#' averaging per-bootstrap member predictions on the requested scale.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Build a deterministic expansion (for example with \code{\link{bigexp_terms}})
#'         and fit several \code{SVEMnet()} models for different responses on
#'         the same factor space, using the same expansion / sampling settings.
#'   \item Ensure that the fitted models were created with a version of
#'         \code{SVEMnet()} that populates \code{$sampling_schema}.
#'   \item Collect the fitted models in a list and pass them to
#'         \code{svem_random_table_multi()}.
#'   \item Use \code{$data} (predictors), \code{$pred} (response columns), or
#'         \code{$all} (\code{cbind(data, pred)}) for downstream plotting,
#'         summarization, or cross-response comparison.
#' }
#'
#' @section Blocking variables:
#' If the models were fit using a \code{bigexp_spec} that included blocking
#' variables (for example \code{blocking = c("Operator", "Plate_ID")}) and
#' \code{SVEMnet()} stored these in \code{$sampling_schema$blocking}, then
#' \code{svem_random_table_multi()} will:
#' \itemize{
#'   \item treat those variables as blocking factors; and
#'   \item hold them fixed at a single value across the sampled table.
#' }
#'
#' Specifically:
#' \itemize{
#'   \item For blocking numeric variables, the function uses the midpoint of the
#'         recorded numeric range, \code{(min + max) / 2}, for all rows.
#'   \item For blocking categorical variables, the function uses a single
#'   reference level equal to the most frequent observed level (mode) in the
#'   training data, with ties broken deterministically; if the mode is
#'   unavailable, it falls back to the first stored level.
#' }
#'
#' Blocking variables are not allowed to appear in \code{mixture_groups}. If
#' any mixture group tries to use a blocking variable, the function stops with
#' an error.
#'
#' When no blocking information is present in \code{$sampling_schema} (for
#' example for models fit without a \code{bigexp_spec} or without blocking),
#' the behavior is unchanged from earlier versions: all predictors are sampled
#' according to the rules described under "Sampling strategy".
#'
#' @param objects A list of fitted \code{svem_model} objects returned by
#'   \code{SVEMnet()}. Each object must contain a valid \code{$sampling_schema}
#'   produced by the updated \code{SVEMnet()} implementation. A single model is
#'   also accepted and treated as a length-one list.
#' @param n Number of random points to generate (rows in the output tables).
#'   Default is \code{1000}.
#' @param mixture_groups Optional list of mixture constraint groups. Each group
#'   is a list with elements \code{vars}, \code{lower}, \code{upper}, \code{total}
#'   (see \emph{Notes on mixtures}). Mixture variables must be numeric-like and
#'   must also appear in the models' \code{predictor_vars} (that is, they must
#'   be used as predictors in all models).
#' @param debias Logical; if \code{TRUE}, apply each model's calibration during
#'   prediction when available (for Gaussian fits). This is passed to
#'   \code{predict.svem_model()}. Default is \code{FALSE}.
#' @param range_tol Numeric tolerance for comparing numeric ranges across models
#'   (used when checking that all \code{$sampling_schema$num_ranges} agree).
#'   Default is \code{1e-8}.
#' @param numeric_sampler Sampler for non-mixture numeric predictors:
#'   \code{"random"} (default), or \code{"uniform"}.
#'   \itemize{
#'     \item \code{"random"}: random Latin hypercube when the \pkg{lhs} package
#'           is available; otherwise independent uniforms via \code{runif()}.
#'     \item \code{"uniform"}: independent uniform draws within numeric ranges
#'           (fastest; no \pkg{lhs} dependency).
#'   }
#'
#' @return A list with three data frames:
#' \itemize{
#'   \item \code{data}: the sampled predictor settings, one row per random point.
#'   \item \code{pred}: one column per response, aligned to \code{data} rows.
#'   \item \code{all}: \code{cbind(data, pred)} for convenience.
#' }
#' Each prediction column is named by the model's response (left-hand side)
#' with a "_pred" suffix (for example, "y1_pred"). If that name would collide
#' with a predictor name or with another prediction column, the function stops
#' with an error and asks the user to rename the response or predictor.
#'
#' @details
#' All models must share an identical predictor schema. Specifically, their
#' \code{$sampling_schema} entries must agree on:
#' \itemize{
#'   \item The same \code{predictor_vars} in the same order.
#'   \item The same \code{var_classes} for each predictor.
#'   \item Identical factor \code{levels} and level order for all categorical
#'         predictors.
#'   \item Numeric \code{num_ranges} that match within \code{range_tol} for all
#'         continuous predictors.
#'   \item When present, the same \code{blocking} set (up to order).
#' }
#' The function stops with an informative error message if any of these checks fail.
#'
#' Models may be Gaussian or binomial. For binomial fits, predictions are
#' returned on the probability scale (that is, on the response scale) by default,
#' consistent with the default behaviour of \code{predict.svem_model()}.
#'
#' @section Sampling strategy:
#' Non-mixture numeric variables are sampled using the chosen \code{numeric_sampler}
#' within the numeric ranges recorded in \code{$sampling_schema$num_ranges}:
#' \itemize{
#'   \item \code{"random"}: random Latin hypercube when \pkg{lhs} is available,
#'         else independent uniforms on each range.
#'   \item \code{"uniform"}: independent uniform draws within numeric ranges
#'         (fastest; no \pkg{lhs} dependency).
#' }
#'
#' Mixture variables (if any) are sampled jointly within each specified group using
#' a truncated Dirichlet so that elementwise bounds and the total sum are satisfied.
#' Categorical variables are sampled from cached factor levels. Blocking variables
#' (if present) are held fixed (single level or single numeric value) and are not
#' randomized.
#'
#' The same random predictor table is fed to each model so response columns are
#' directly comparable.
#'
#' @section Notes on mixtures:
#' Each mixture group should list only numeric-like variables. Bounds are interpreted
#' on the original scale of those variables. If \code{total} equals the sum of lower
#' bounds, the sampler returns the lower-bound corner for that group. Infeasible
#' constraints (that is, \code{sum(lower) > total} or \code{sum(upper) < total})
#' produce an error.
#'
#' Mixture variables are removed from the pool of "non-mixture" numeric variables
#' before numeric sampling, so they are controlled entirely by the mixture
#' constraints and not also sampled independently. Mixture variables are not
#' allowed to be blocking variables.
#'
#' @seealso \code{\link{SVEMnet}}, \code{\link{predict.svem_model}},
#'   \code{\link{bigexp_terms}}, \code{\link{bigexp_formula}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' X1 <- runif(n); X2 <- runif(n)
#' A <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)
#' F <- factor(sample(c("lo","hi"), n, TRUE))
#'
#' ## Gaussian responses
#' y1 <- 1 + 2*X1 - X2 + 3*A + 1.5*B + 0.5*C + (F=="hi") + rnorm(n, 0, 0.3)
#' y2 <- 0.5 + 0.8*X1 + 0.4*X2 + rnorm(n, 0, 0.2)
#'
#' ## Binomial response (probability via logistic link)
#' eta  <- -0.5 + 1.2*X1 - 0.7*X2 + 0.8*(F=="hi") + 0.6*A
#' p    <- 1 / (1 + exp(-eta))
#' yb   <- rbinom(n, size = 1, prob = p)
#'
#' d  <- data.frame(y1, y2, yb, X1, X2, A, B, C, F)
#'
#' fit1 <- SVEMnet(y1 ~ X1 + X2 + A + B + C + F, d, nBoot = 40, family = "gaussian")
#' fit2 <- SVEMnet(y2 ~ X1 + X2 + A + B + C + F, d, nBoot = 40, family = "gaussian")
#' fitb <- SVEMnet(yb ~ X1 + X2 + A + B + C + F, d, nBoot = 40, family = "binomial")
#'
#' # Mixture constraint for A, B, C that sum to 1
#' mix <- list(list(vars  = c("A","B","C"),
#'                  lower = c(0,0,0),
#'                  upper = c(1,1,1),
#'                  total = 1))
#'
#' # Fast random sampler (shared predictor table; predictions bound as columns)
#' tab_fast <- svem_random_table_multi(
#'   objects         = list(y1 = fit1, y2 = fit2, yb = fitb),
#'   n               = 2000,
#'   mixture_groups  = mix,
#'   debias          = FALSE,
#'   numeric_sampler = "random"
#' )
#' head(tab_fast$all)
#'
#' # Check that the binomial predictions are on [0,1]
#' range(tab_fast$pred$yb_pred)
#'
#' # Uniform sampler (fastest)
#' tab_uni <- svem_random_table_multi(
#'   objects         = list(y1 = fit1, y2 = fit2, yb = fitb),
#'   n               = 2000,
#'   debias          = FALSE,
#'   numeric_sampler = "uniform"
#' )
#' head(tab_uni$all)
#'
#' ## Example with blocking (requires SVEMnet to store sampling_schema$blocking)
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
#'
#' ## Operator and AmbientTmp are held fixed across rows:
#' length(unique(tab_block$data$Operator))
#' range(tab_block$data$AmbientTmp)
#' }
#'
#' @importFrom lhs randomLHS
#' @importFrom stats rgamma
#' @export
svem_random_table_multi <- function(objects, n = 1000, mixture_groups = NULL,
                                    debias = FALSE, range_tol = 1e-8,
                                    numeric_sampler = c("random", "uniform")) {
  numeric_sampler <- match.arg(numeric_sampler)

  # ---- validate inputs ----
  if (inherits(objects, "svem_model")) objects <- list(objects)
  if (!is.list(objects) || length(objects) < 1L)
    stop("'objects' must be a non-empty list of 'svem_model' fits.")
  if (!all(vapply(objects, function(o) inherits(o, "svem_model"), logical(1))))
    stop("All elements of 'objects' must be 'svem_model' objects.")
  if (!all(vapply(objects, function(o) is.list(o$sampling_schema), logical(1))))
    stop("Each 'svem_model' must contain a valid $sampling_schema. Refit with updated SVEMnet().")

  # Reference schema from first object
  ref <- objects[[1L]]$sampling_schema
  ref_vars       <- ref$predictor_vars
  ref_classes    <- ref$var_classes
  ref_ranges     <- ref$num_ranges   # 2 x K matrix, rownames c("min","max")
  ref_fac_levels <- ref$factor_levels
  ref_blocking   <- ref$blocking
  if (is.null(ref_blocking)) ref_blocking <- character(0L)
  ref_blocking   <- intersect(unique(as.character(ref_blocking)), ref_vars)

  if (!length(ref_vars)) stop("No predictor variables found in reference sampling_schema.")

  # Helper to compare numeric ranges with tolerance
  ranges_equal <- function(A, B, tol) {
    if (is.null(A) && is.null(B)) return(TRUE)
    if (is.null(A) || is.null(B)) return(FALSE)
    AA <- as.matrix(A); BB <- as.matrix(B)
    if (!identical(rownames(AA), rownames(BB))) return(FALSE)
    if (!identical(colnames(AA), colnames(BB))) return(FALSE)
    if (!all(dim(AA) == dim(BB))) return(FALSE)
    diff <- AA - BB
    diff[!is.finite(diff)] <- Inf
    all(abs(diff) <= tol)
  }

  # Check all schemas match reference (predictors, classes, levels, ranges, blocking)
  for (k in seq_along(objects)) {
    s   <- objects[[k]]$sampling_schema
    msg <- sprintf("Object %d: ", k)

    if (!identical(s$predictor_vars, ref_vars))
      stop(msg, "predictor_vars do not match the reference (set and/or order).")

    if (!identical(s$var_classes, ref_classes))
      stop(msg, "var_classes do not match the reference.")

    s_levels <- s$factor_levels
    if (!identical(names(s_levels), names(ref_fac_levels)))
      stop(msg, "factor_levels names differ from the reference.")
    for (nm in names(ref_fac_levels)) {
      if (!identical(s_levels[[nm]], ref_fac_levels[[nm]]))
        stop(msg, "factor_levels for '", nm, "' differ from the reference.")
    }

    s_ranges <- s$num_ranges
    if (!is.null(ref_ranges) && ncol(ref_ranges) &&
        !is.null(s_ranges)   && ncol(s_ranges)) {
      s_ranges <- s_ranges[, colnames(ref_ranges), drop = FALSE]
    }
    if (!ranges_equal(ref_ranges, s_ranges, range_tol))
      stop(msg, "numeric ranges differ from the reference (beyond tolerance).")

    # Blocking must match across models (up to order)
    s_block <- s$blocking
    if (is.null(s_block)) s_block <- character(0L)
    s_block <- intersect(unique(as.character(s_block)), ref_vars)
    if (!identical(sort(s_block), sort(ref_blocking))) {
      stop(msg, "blocking set in sampling_schema does not match the reference.")
    }

    ## ---- OPTIONAL: soft check for block_cat_modes consistency ----
    if (!is.null(ref$block_cat_modes) && length(ref$block_cat_modes)) {
      ref_modes <- ref$block_cat_modes
      s_modes   <- s$block_cat_modes

      # If the other model has no modes stored, just warn and continue
      if (is.null(s_modes) || !length(s_modes)) {
        warning(
          msg,
          "sampling_schema$block_cat_modes is missing; ",
          "using modes from the first object only."
        )
      } else {
        ref_names <- names(ref_modes)
        s_names   <- names(s_modes)

        if (!identical(sort(ref_names), sort(s_names))) {
          warning(
            msg,
            "block_cat_modes names differ from the reference; ",
            "using reference modes for blocking variables."
          )
        } else {
          for (nm in ref_names) {
            if (!identical(ref_modes[[nm]], s_modes[[nm]])) {
              warning(
                msg,
                "block_cat_modes for '", nm,
                "' differ from the reference; using reference modes."
              )
              break
            }
          }
        }
      }
    }
  }


  # Aliases
  predictor_vars <- ref_vars
  var_classes    <- ref_classes
  num_ranges     <- ref_ranges
  factor_levels  <- ref_fac_levels
  blocking       <- ref_blocking

  numeric_like <- c("numeric", "double", "integer", "integer64")

  all_num <- names(var_classes)[var_classes %in% numeric_like]

  # ---- REQUIRE num_ranges for all numeric predictors ----
  if (length(all_num)) {
    if (is.null(num_ranges) || !is.matrix(num_ranges) || !ncol(num_ranges)) {
      stop("sampling_schema$num_ranges must be a 2 x K matrix with ranges for all numeric predictors.")
    }
    if (!all(c("min", "max") %in% rownames(num_ranges))) {
      stop("sampling_schema$num_ranges must have rownames 'min' and 'max'.")
    }
    missing_rng <- setdiff(all_num, colnames(num_ranges))
    if (length(missing_rng)) {
      stop("Missing numeric ranges in sampling_schema$num_ranges for predictors: ",
           paste(missing_rng, collapse = ", "))
    }
  }

  # Partition predictors by type and blocking
  block_num <- intersect(blocking, all_num)
  block_cat <- setdiff(blocking, block_num)

  is_num <- setdiff(all_num, blocking)
  is_cat <- setdiff(setdiff(predictor_vars, all_num), blocking)

  # ---- mixture validation ----
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) {
      if (is.null(grp$vars))
        stop("Each mixture group must contain a 'vars' character vector.")
      if (!all(grp$vars %in% predictor_vars)) {
        missing <- setdiff(grp$vars, predictor_vars)
        stop("Mixture variables not in model predictors: ",
             paste(missing, collapse = ", "))
      }
      bad_mix <- setdiff(grp$vars, names(var_classes)[var_classes %in% numeric_like])
      if (length(bad_mix)) {
        stop("Mixture variables must be numeric-like. Non-numeric mixture vars: ",
             paste(bad_mix, collapse = ", "))
      }
      # Mixture variables cannot be blocking variables
      if (length(intersect(grp$vars, blocking))) {
        stop("Mixture variables cannot be blocking variables. Offending vars: ",
             paste(intersect(grp$vars, blocking), collapse = ", "))
      }
      mixture_vars <- c(mixture_vars, grp$vars)
    }
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ",
           paste(dups, collapse = ", "))
    }
  }

  # Truncated Dirichlet sampler
  .sample_trunc_dirichlet <- function(n, lower, upper, total,
                                      alpha = NULL, oversample = 4L,
                                      max_tries = 10000L) {
    k <- length(lower)
    if (length(upper) != k)
      stop("upper must have the same length as lower.")
    if (is.null(alpha)) alpha <- rep(1, k)
    min_sum <- sum(lower); max_sum <- sum(upper)
    if (total < min_sum - 1e-12 || total > max_sum + 1e-12)
      stop("Infeasible mixture constraints: need sum(lower) <= total <= sum(upper).")
    avail <- total - min_sum
    if (avail <= 1e-12)
      return(matrix(rep(lower, each = n), nrow = n))
    res    <- matrix(NA_real_, nrow = n, ncol = k)
    filled <- 0L
    tries  <- 0L
    while (filled < n && tries < max_tries) {
      m <- max(oversample * (n - filled), 1L)
      g <- matrix(stats::rgamma(m * k, shape = alpha, rate = 1),
                  ncol = k, byrow = TRUE)
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
    if (filled < n)
      stop("Could not sample enough feasible mixture points within max_tries.")
    res
  }

  # ---- sample non-mixture numerics (excluding blocking) ----
  nonmix_num <- setdiff(is_num, mixture_vars)
  T_num <- NULL
  if (length(nonmix_num)) {
    rng <- vapply(nonmix_num, function(v) {
      r <- num_ranges[, v]
      if (!all(is.finite(r))) {
        stop("Numeric range for predictor '", v, "' must be finite.")
      }
      if (r[1] > r[2]) {
        stop("Numeric range for predictor '", v,
             "' must have min <= max. Check sampling_schema$num_ranges.")
      }
      as.numeric(r)
    }, numeric(2))
    rownames(rng) <- c("min", "max")
    lo    <- rng["min", ]
    hi    <- rng["max", ]
    width <- hi - lo
    q     <- length(nonmix_num)

    use_lhs <- function() isTRUE(requireNamespace("lhs", quietly = TRUE))
    U <- switch(numeric_sampler,
                "random" = {
                  if (use_lhs()) lhs::randomLHS(n, q)
                  else matrix(stats::runif(n * q), nrow = n, ncol = q)
                },
                "uniform" = {
                  matrix(stats::runif(n * q), nrow = n, ncol = q)
                }
    )
    T_num <- sweep(U, 2, width, `*`)
    T_num <- sweep(T_num, 2, lo, `+`)
    colnames(T_num) <- nonmix_num
    T_num <- as.data.frame(T_num)
  }

  # ---- sample mixture groups ----
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
      if (length(lower) != k || length(upper) != k)
        stop("lower/upper must each have length equal to mixture vars: ",
             paste(vars, collapse = ", "))
      vals <- .sample_trunc_dirichlet(n, lower, upper, total)
      colnames(vals) <- vars
      T_mix[, vars] <- vals
    }
    T_mix <- as.data.frame(T_mix)
  }

  # ---- sample categoricals (excluding blocking) ----
  T_cat <- NULL
  if (length(is_cat)) {
    T_cat <- vector("list", length(is_cat))
    names(T_cat) <- is_cat
    for (v in is_cat) {
      lev <- factor_levels[[v]]
      if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
      if (is.null(lev) || !length(lev)) lev <- c("L1", "L2")
      T_cat[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
    }
    T_cat <- as.data.frame(T_cat, stringsAsFactors = FALSE)
  }

  # ---- assemble predictors for non-blocking vars ----
  parts <- list(T_num, T_mix, T_cat)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts) && !length(blocking))
    stop("No predictors could be sampled from the schema.")
  T_data <- if (length(parts)) do.call(cbind, parts) else data.frame()

  # ---- ensure all predictors present; handle blocking specially ----
  missing_pred <- setdiff(predictor_vars, colnames(T_data))
  if (length(missing_pred)) {
    for (v in missing_pred) {
      is_block <- v %in% blocking
      is_num_like <- v %in% all_num

      if (is_block && is_num_like) {
        # Blocking numeric: use mid-point of recorded range
        r <- num_ranges[, v]
        val <- (r["min"] + r["max"]) / 2
        T_data[[v]] <- rep(as.numeric(val), n)
      } else if (is_block && !is_num_like) {
        # Blocking categorical: use mode if available, else first level
        lev <- factor_levels[[v]]
        if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
        if (is.null(lev) || !length(lev)) lev <- c("L1", "L2")

        # Prefer stored mode from sampling_schema
        mode_val <- NULL
        if (!is.null(ref$block_cat_modes) &&
            length(ref$block_cat_modes) &&
            !is.null(ref$block_cat_modes[[v]])) {
          mode_val <- as.character(ref$block_cat_modes[[v]])
        }

        if (is.null(mode_val) || is.na(mode_val) || !(mode_val %in% lev)) {
          mode_val <- lev[1L]
        }

        T_data[[v]] <- factor(rep(mode_val, n), levels = lev)

      } else if (!is_block && is_num_like) {
        # Non-blocking numeric fallback: mid-point of range
        r <- num_ranges[, v]
        if (!all(is.finite(r))) {
          stop("Numeric range for predictor '", v, "' must be finite.")
        }
        val <- (r["min"] + r["max"]) / 2
        T_data[[v]] <- rep(as.numeric(val), n)
      } else {
        # Non-blocking categorical fallback: random sample
        lev <- factor_levels[[v]]
        if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
        if (is.null(lev) || !length(lev)) lev <- c("L1", "L2")
        T_data[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
      }
    }
  }

  # Order like predictor_vars and ensure blocking columns exist (constant)
  T_data <- T_data[, predictor_vars, drop = FALSE]

  # ---- predict for each model on the shared T_data ----
  data_df <- T_data
  pred_df <- as.data.frame(matrix(NA_real_, nrow = nrow(data_df), ncol = 0))
  rownames(pred_df) <- NULL

  for (i in seq_along(objects)) {
    obj   <- objects[[i]]
    preds <- predict(obj, newdata = data_df, debias = debias)
    if (is.list(preds) && !is.null(preds$fit)) preds <- preds$fit

    # Response name (LHS of formula), with fallback
    resp <- tryCatch(
      as.character(obj$formula[[2L]]),
      error = function(e) paste0("resp", i)
    )

    # Always use "<response>_pred" as the prediction column name
    base_colname <- paste0(resp, "_pred")

    # If this would collide with a predictor, stop and ask user to rename
    if (base_colname %in% colnames(data_df)) {
      stop(
        "Prediction column name '", base_colname,
        "' would collide with an existing predictor name. ",
        "Please rename the response or the predictor to avoid using the '_pred' suffix."
      )
    }

    # If this would collide with an existing prediction column, stop
    if (base_colname %in% colnames(pred_df)) {
      stop(
        "Prediction column name '", base_colname,
        "' is duplicated across models. ",
        "Please ensure each response name is unique."
      )
    }

    pred_df[[base_colname]] <- as.numeric(preds)
  }

  rownames(data_df) <- NULL
  all_df <- cbind(data_df, pred_df)

  list(data = data_df, pred = pred_df, all = all_df)
}
