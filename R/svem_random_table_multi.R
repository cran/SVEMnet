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
#'         recorded numeric range, \code{(min + max) / 2}, for all rows. If the
#'         variable also has stored discrete support, the midpoint is snapped
#'         deterministically to the nearest allowed discrete value.
#'   \item For blocking categorical variables, the function uses a single
#'         reference level equal to the most frequent observed level (mode) in the
#'         training data, with ties broken deterministically; if the mode is
#'         unavailable, it falls back to the first stored level.
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
#' \strong{Discrete numeric predictors (automatic).}
#' If any supplied model stores discrete-numeric sampling information in its
#' \code{$sampling_schema}, this function will automatically respect it (no
#' separate user argument).
#'
#' In the updated \code{SVEMnet()} implementation this information is stored as:
#' \itemize{
#'   \item \code{$sampling_schema$discrete_numeric}: a character vector of discrete
#'         numeric variable names; and
#'   \item \code{$sampling_schema$discrete_levels}: a named list mapping those
#'         names to allowed numeric values.
#' }
#' (Older objects may use \code{$sampling_schema$discrete_values} instead of
#' \code{discrete_levels}; this function accepts both for backward compatibility.)
#'
#' Discrete numeric variables are sampled independently (uniform over allowed
#' values) and are excluded from Latin hypercube sampling; LHS (when used) is
#' applied only to the remaining continuous numeric predictors. Discrete numeric
#' variables are not allowed to be mixture variables.
#'
#' Models may be Gaussian or binomial. For binomial fits, predictions are
#' returned on the probability scale (that is, on the response scale) by default,
#' consistent with the default behavior of \code{predict.svem_model()}.
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

  n <- .svem_integer_scalar(n, "n", min = 1L)
  range_tol <- .svem_numeric_scalar(range_tol, "range_tol", lower = 0)
  debias <- .svem_logical_scalar(debias, "debias")

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

  # ---- helper: extract + normalize discrete numeric map from a sampling_schema ----
  .norm_discrete_map <- function(x) {
    if (is.null(x) || !is.list(x) || !length(x)) return(list())
    nm <- names(x)
    if (is.null(nm) || any(!nzchar(nm))) return(list())
    out <- vector("list", length(x))
    names(out) <- nm
    for (v in nm) {
      vals <- sort(unique(as.numeric(x[[v]])))
      vals <- vals[is.finite(vals)]
      out[[v]] <- vals
    }
    out <- out[vapply(out, length, integer(1)) > 0]
    out
  }

  .extract_discrete_map <- function(ss) {
    if (is.null(ss) || !is.list(ss)) return(list())

    dn <- ss$discrete_numeric
    # New name (SVEMnet stores this)
    dl <- ss$discrete_levels
    # Backward-compatible alias (older objects)
    dv <- ss$discrete_values

    # If both exist, prefer discrete_levels
    if (is.null(dl) && !is.null(dv)) dl <- dv

    # Case 1: discrete_numeric is a named list of supports (accept)
    if (is.list(dn) && !is.null(names(dn)) && all(nzchar(names(dn)))) {
      return(.norm_discrete_map(dn))
    }

    # Case 2: discrete_levels/discrete_values is a named list -> treat as discrete
    if (is.list(dl) && length(dl) && !is.null(names(dl)) && all(nzchar(names(dl)))) {
      # If dn is a character vector, subset to those names (if provided)
      if (is.character(dn) && length(dn)) {
        keep <- intersect(unique(as.character(dn)), names(dl))
        return(.norm_discrete_map(dl[keep]))
      }
      return(.norm_discrete_map(dl))
    }

    # Case 3: discrete_numeric is a character vector but no levels list -> error
    if (is.character(dn) && length(dn)) {
      stop(
        "sampling_schema$discrete_numeric is present (character vector) but ",
        "sampling_schema$discrete_levels (or legacy $discrete_values) is missing; ",
        "cannot resolve discrete supports."
      )
    }

    list()
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
    # names order shouldn't matter; compare as sets then compare by name
    if (!setequal(names(s_levels), names(ref_fac_levels)))
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
            "using reference modes for blocking variables. This is a low-priority note indicating which reference level is used for the blocking factor and can be ignored."
          )
        } else {
          for (nm in ref_names) {
            if (!identical(ref_modes[[nm]], s_modes[[nm]])) {
              warning(
                msg,
                "block_cat_modes for '", nm,
                "' differ from the reference; using reference modes. This is a low-priority note indicating which reference level is used for the blocking factor and can be ignored."
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

  # ---- AUTO discrete numeric map (from ANY supplied object) ----
  disc_maps <- lapply(objects, function(o) .extract_discrete_map(o$sampling_schema))
  has_disc  <- vapply(disc_maps, function(m) is.list(m) && length(m) > 0, logical(1))

  discrete_map <- list()
  if (any(has_disc)) {
    idx_ref <- which(has_disc)[1L]
    discrete_map <- disc_maps[[idx_ref]]

    # Validate consistency across models that also specify discrete info
    for (k in which(has_disc)) {
      mk <- disc_maps[[k]]
      if (!identical(sort(names(mk)), sort(names(discrete_map)))) {
        stop(
          "Discrete numeric specifications differ across models (variable sets differ). ",
          "Please refit models with a consistent sampling_schema."
        )
      }
      for (nm in names(discrete_map)) {
        a <- sort(unique(as.numeric(discrete_map[[nm]]))); a <- a[is.finite(a)]
        b <- sort(unique(as.numeric(mk[[nm]])));           b <- b[is.finite(b)]
        if (!identical(a, b)) {
          stop(
            "Discrete numeric specifications differ across models for variable '", nm, "'. ",
            "Please refit models with a consistent sampling_schema."
          )
        }
      }
    }

    bad_not_pred <- setdiff(names(discrete_map), predictor_vars)
    if (length(bad_not_pred)) {
      stop(
        "sampling_schema discrete numeric variables not found in predictors: ",
        paste(bad_not_pred, collapse = ", ")
      )
    }
    bad_not_num <- setdiff(names(discrete_map), all_num)
    if (length(bad_not_num)) {
      stop(
        "sampling_schema discrete numeric variables are not numeric-like: ",
        paste(bad_not_num, collapse = ", ")
      )
    }

    for (v in names(discrete_map)) {
      vals <- sort(unique(as.numeric(discrete_map[[v]])))
      vals <- vals[is.finite(vals)]
      if (!length(vals)) {
        stop("Discrete support for '", v, "' is empty or non-finite in sampling_schema.")
      }
      r <- num_ranges[, v]
      lo <- as.numeric(r["min"]); hi <- as.numeric(r["max"])
      if (all(is.finite(c(lo, hi)))) {
        if (any(vals < lo - 1e-12 | vals > hi + 1e-12)) {
          stop(
            "Discrete values for '", v, "' fall outside the recorded numeric range [",
            lo, ", ", hi, "] in sampling_schema$num_ranges."
          )
        }
      }
      discrete_map[[v]] <- vals
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

  # Discrete numeric variables cannot be mixture variables
  if (length(discrete_map)) {
    overlap_disc_mix <- intersect(names(discrete_map), mixture_vars)
    if (length(overlap_disc_mix)) {
      stop(
        "Discrete numeric variables are not allowed to be mixture variables. Offending vars: ",
        paste(overlap_disc_mix, collapse = ", ")
      )
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

  disc_num <- character(0L)
  if (length(discrete_map)) {
    disc_num <- intersect(nonmix_num, names(discrete_map))
  }
  cont_num <- setdiff(nonmix_num, disc_num)

  # continuous numeric block (LHS/uniform)
  T_num_cont <- NULL
  if (length(cont_num)) {
    rng <- vapply(cont_num, function(v) {
      r <- num_ranges[, v]
      if (!all(is.finite(r))) stop("Numeric range for predictor '", v, "' must be finite.")
      if (r[1] > r[2]) stop("Numeric range for predictor '", v, "' must have min <= max.")
      as.numeric(r)
    }, numeric(2))
    rownames(rng) <- c("min", "max")
    lo    <- rng["min", ]
    hi    <- rng["max", ]
    width <- hi - lo
    q     <- length(cont_num)

    use_lhs <- function() isTRUE(requireNamespace("lhs", quietly = TRUE))
    U <- switch(numeric_sampler,
                "random" = {
                  if (q == 0) matrix(numeric(0), nrow = n, ncol = 0) else
                    if (use_lhs()) lhs::randomLHS(n, q) else matrix(stats::runif(n * q), nrow = n, ncol = q)
                },
                "uniform" = {
                  if (q == 0) matrix(numeric(0), nrow = n, ncol = 0) else
                    matrix(stats::runif(n * q), nrow = n, ncol = q)
                }
    )

    if (q > 0) {
      T_num_cont <- sweep(U, 2, width, `*`)
      T_num_cont <- sweep(T_num_cont, 2, lo, `+`)
      colnames(T_num_cont) <- cont_num
      T_num_cont <- as.data.frame(T_num_cont)
    }
  }

  # discrete numeric block
  T_num_disc <- NULL
  if (length(disc_num)) {
    T_num_disc <- vector("list", length(disc_num))
    names(T_num_disc) <- disc_num
    for (v in disc_num) {
      vals <- discrete_map[[v]]
      vals <- sort(unique(as.numeric(vals)))
      vals <- vals[is.finite(vals)]
      if (!length(vals)) stop("Discrete support for '", v, "' is empty or non-finite.")
      # index-based sampling avoids sample()'s 1:x behavior for length-one supports
      T_num_disc[[v]] <- vals[sample.int(length(vals), n, replace = TRUE)]
    }
    T_num_disc <- as.data.frame(T_num_disc, stringsAsFactors = FALSE)
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
        stop("lower/upper must each have length equal to mixture vars: ", paste(vars, collapse = ", "))
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
      if (is.null(lev) || !length(lev)) {
        stop("No recorded levels for categorical predictor '", v,
             "' in sampling_schema or xlevels; cannot sample it. ",
             "Refit the model with a version of SVEMnet() that records factor levels.")
      }
      # Reproduce the training class: predictors fitted as ordered factors
      # must be sampled as ordered or predict()'s class validator rejects
      # the table.
      T_cat[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev,
                           ordered = identical(unname(var_classes[v]), "ordered"))
    }
    T_cat <- as.data.frame(T_cat, stringsAsFactors = FALSE)
  }

  # ---- assemble predictors for non-blocking vars ----
  parts <- list(T_num_cont, T_num_disc, T_mix, T_cat)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts) && !length(blocking))
    stop("No predictors could be sampled from the schema.")

  # Robust init when only blocking (or only missing_pred fill) will populate columns
  T_data <- if (length(parts)) do.call(cbind, parts) else data.frame(.row = seq_len(n))
  if (".row" %in% names(T_data)) {
    # keep only as a scaffold; it will be dropped by final subsetting anyway
  }

  # ---- ensure all predictors present; handle blocking specially ----
  missing_pred <- setdiff(predictor_vars, colnames(T_data))
  if (length(missing_pred)) {
    for (v in missing_pred) {
      is_block    <- v %in% blocking
      is_num_like <- v %in% all_num

      if (is_block && is_num_like) {
        r <- num_ranges[, v]
        mid <- (r["min"] + r["max"]) / 2
        val <- as.numeric(mid)

        if (!is.null(discrete_map[[v]])) {
          vals <- sort(unique(as.numeric(discrete_map[[v]])))
          vals <- vals[is.finite(vals)]
          if (length(vals)) {
            val <- vals[which.min(abs(vals - val))]
          }
        }
        T_data[[v]] <- rep(as.numeric(val), n)

      } else if (is_block && !is_num_like) {
        lev <- factor_levels[[v]]
        if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
        if (is.null(lev) || !length(lev)) {
          stop("No recorded levels for blocking variable '", v,
               "' in sampling_schema or xlevels; cannot construct its reference level. ",
               "Refit the model with a version of SVEMnet() that records factor levels.")
        }

        mode_val <- NULL
        if (!is.null(ref$block_cat_modes) &&
            length(ref$block_cat_modes) &&
            !is.null(ref$block_cat_modes[[v]])) {
          mode_val <- as.character(ref$block_cat_modes[[v]])
        }
        if (is.null(mode_val) || is.na(mode_val) || !(mode_val %in% lev)) {
          mode_val <- lev[1L]
        }
        T_data[[v]] <- factor(rep(mode_val, n), levels = lev,
                              ordered = identical(unname(var_classes[v]), "ordered"))

      } else if (!is_block && is_num_like) {
        if (!is.null(discrete_map[[v]])) {
          vals <- sort(unique(as.numeric(discrete_map[[v]])))
          vals <- vals[is.finite(vals)]
          if (!length(vals)) stop("Discrete support for '", v, "' is empty or non-finite.")
          # index-based sampling avoids sample()'s 1:x behavior for length-one supports
          T_data[[v]] <- vals[sample.int(length(vals), n, replace = TRUE)]
        } else {
          r <- num_ranges[, v]
          if (!all(is.finite(r))) stop("Numeric range for predictor '", v, "' must be finite.")
          val <- (r["min"] + r["max"]) / 2
          T_data[[v]] <- rep(as.numeric(val), n)
        }

      } else {
        lev <- factor_levels[[v]]
        if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
        if (is.null(lev) || !length(lev)) {
          stop("No recorded levels for categorical predictor '", v,
               "' in sampling_schema or xlevels; cannot sample it. ",
               "Refit the model with a version of SVEMnet() that records factor levels.")
        }
        T_data[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev,
                              ordered = identical(unname(var_classes[v]), "ordered"))
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

    # deparse-based label stays a length-1 string for transformed responses
    # such as log(y), where as.character(formula[[2L]]) would be length 2 and
    # crash the scalar %in% checks below
    resp <- .svem_response_label(obj$formula)
    if (!is.character(resp) || length(resp) != 1L || is.na(resp) || !nzchar(resp)) {
      resp <- paste0("resp", i)
    }

    base_colname <- paste0(resp, "_pred")

    if (base_colname %in% colnames(data_df)) {
      stop(
        "Prediction column name '", base_colname,
        "' would collide with an existing predictor name. ",
        "Please rename the response or the predictor to avoid using the '_pred' suffix."
      )
    }
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
