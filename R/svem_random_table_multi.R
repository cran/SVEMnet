#' Generate a Random Prediction Table from Multiple SVEMnet Models (no refit)
#'
#' Samples the original predictor factor space cached in fitted \code{svem_model}
#' objects and computes predictions from each model at the same random points.
#' Intended for multiple responses built over the same factor space and a
#' deterministic factor expansion.
#'
#' @param objects A list of fitted \code{svem_model} objects returned by \code{SVEMnet()}.
#'   Each object must contain \code{$sampling_schema} produced by the updated \code{SVEMnet()}.
#'   A single model is also accepted and treated as a length-one list.
#' @param n Number of random points to generate. Default is \code{1000}.
#' @param mixture_groups Optional list of mixture constraint groups. Each group is a list
#'   with elements \code{vars}, \code{lower}, \code{upper}, \code{total}. The variables in
#'   \code{vars} must be numeric-like predictors present in all models. The sampler uses a
#'   truncated Dirichlet so that elementwise bounds are respected and \code{sum(vars) = total}.
#' @param debias Logical; if \code{TRUE}, apply each model's calibration during prediction
#'   when available. Default is \code{FALSE}.
#' @param range_tol Numeric tolerance for comparing numeric ranges across models.
#'   Default is \code{1e-8}.
#'
#' @return A list with three data frames:
#' \itemize{
#'   \item \code{data}: the sampled predictor settings, one row per random point.
#'   \item \code{pred}: one column per response, aligned to \code{data} rows.
#'   \item \code{all}: \code{cbind(data, pred)} for convenience.
#' }
#' Each prediction column is named by the model's response (left-hand side). If a response
#' name would collide with a predictor name, \code{".pred"} is appended.
#'
#' @details
#' All models must share an identical predictor schema:
#' \itemize{
#'   \item The same \code{predictor_vars} in the same order
#'   \item The same \code{var_classes} for each predictor
#'   \item Identical factor \code{levels} and level order
#'   \item Numeric \code{ranges} that match within \code{range_tol}
#' }
#' The function stops with an informative error message if any of these checks fail.
#'
#' @section Sampling strategy:
#' Non-mixture numeric variables are sampled with a maximin Latin hypercube over the
#' cached numeric ranges. Mixture variables are sampled jointly within each specified
#' group using a truncated Dirichlet so that elementwise bounds and the total sum
#' are satisfied. Categorical variables are sampled from the cached factor levels.
#' The same random predictor table is fed to each model so response columns are directly comparable.
#'
#' @section Notes on mixtures:
#' Each mixture group should list only numeric-like variables. Bounds are interpreted
#' on the original scale of those variables. If \code{total} equals the sum of lower
#' bounds, the sampler returns the lower-bound corner for that group.
#'
#' @seealso \code{SVEMnet}, \code{predict.svem_model}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' X1 <- runif(n); X2 <- runif(n)
#' A <- runif(n); B <- runif(n); C <- pmax(0, 1 - A - B)
#' F <- factor(sample(c("lo","hi"), n, TRUE))
#' y1 <- 1 + 2*X1 - X2 + 3*A + 1.5*B + 0.5*C + (F=="hi") + rnorm(n, 0, 0.3)
#' y2 <- 0.5 + 0.8*X1 + 0.4*X2 + rnorm(n, 0, 0.2)
#' d  <- data.frame(y1, y2, X1, X2, A, B, C, F)
#'
#' fit1 <- SVEMnet(y1 ~ X1 + X2 + A + B + C + F, d, nBoot = 40)
#' fit2 <- SVEMnet(y2 ~ X1 + X2 + A + B + C + F, d, nBoot = 40)
#'
#' # Mixture constraint for A, B, C that sum to 1
#' mix <- list(list(vars=c("A","B","C"),
#'                  lower=c(0,0,0), upper=c(1,1,1), total=1))
#'
#' res <- svem_random_table_multi(
#'   list(fit1, fit2), n = 200, mixture_groups = mix,
#'   debias = FALSE
#' )
#' head(res$all)
#'
#' # Lipid screening example with composition group
#' data(lipid_screen)
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data = lipid_screen, factorial_order = 3
#' )
#' fP <- bigexp_formula(spec, "Potency")
#' fS <- bigexp_formula(spec, "Size")
#' fD <- bigexp_formula(spec, "PDI")
#' mP <- SVEMnet(fP, lipid_screen, nBoot = 30)
#' mS <- SVEMnet(fS, lipid_screen, nBoot = 30)
#' mD <- SVEMnet(fD, lipid_screen, nBoot = 30)
#'
#' mixL <- list(list(
#'   vars  = c("Cholesterol","PEG","Ionizable","Helper"),
#'   lower = c(0.10, 0.01, 0.10, 0.10),
#'   upper = c(0.60, 0.05, 0.60, 0.60),
#'   total = 1
#' ))
#'
#' tab <- svem_random_table_multi(
#'   objects        = list(Potency = mP, Size = mS, PDI = mD),
#'   n              = 1000,
#'   mixture_groups = mixL,
#'   debias         = FALSE
#' )
#' head(tab$all)
#' }
#'
#' @importFrom lhs maximinLHS
#' @importFrom stats rgamma
#' @export
svem_random_table_multi <- function(objects, n = 1000, mixture_groups = NULL,
                                    debias = FALSE, range_tol = 1e-8) {
  # ---- validate inputs ----
  if (inherits(objects, "svem_model")) objects <- list(objects)
  if (!is.list(objects) || length(objects) < 1L)
    stop("'objects' must be a non-empty list of 'svem_model' fits.")
  if (!all(vapply(objects, function(o) inherits(o, "svem_model"), logical(1))))
    stop("All elements of 'objects' must be 'svem_model' objects.")

  # All must have sampling_schema
  if (!all(vapply(objects, function(o) is.list(o$sampling_schema), logical(1))))
    stop("Each 'svem_model' must contain a valid $sampling_schema. Refit with updated SVEMnet().")

  # Reference schema from first object
  ref <- objects[[1L]]$sampling_schema
  ref_vars       <- ref$predictor_vars
  ref_classes    <- ref$var_classes
  ref_ranges     <- ref$num_ranges   # 2 x K matrix, rownames c("min","max")
  ref_fac_levels <- ref$factor_levels
  if (!length(ref_vars)) stop("No predictor variables found in reference sampling_schema.")

  # Helper to compare numeric ranges with tolerance
  ranges_equal <- function(A, B, tol) {
    if (is.null(dim(A)) && is.null(dim(B))) return(TRUE)
    if (is.null(A) || is.null(B)) return(FALSE)
    if (!identical(rownames(A), rownames(B))) return(FALSE)
    if (!identical(colnames(A), colnames(B))) return(FALSE)
    AA <- as.matrix(A); BB <- as.matrix(B)
    if (!all(dim(AA) == dim(BB))) return(FALSE)
    diff <- AA - BB
    diff[!is.finite(diff)] <- Inf
    all(abs(diff) <= tol)
  }

  # Check all schemas match reference (predictors, classes, levels, ranges)
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
  }

  # Prevent duplicate responses
  resp_names <- vapply(objects, function(o) {
    nm <- tryCatch(as.character(o$formula[[2L]]), error = function(e) NA_character_)
    if (!is.character(nm) || !nzchar(nm)) NA_character_ else nm
  }, character(1))
  if (anyNA(resp_names))
    stop("At least one model has an invalid or missing response name.")
  if (any(duplicated(resp_names))) {
    dupes <- unique(resp_names[duplicated(resp_names)])
    stop("Duplicate response names across models: ", paste(dupes, collapse = ", "),
         ". Fit responses must be unique.")
  }

  # Aliases
  predictor_vars <- ref_vars
  var_classes    <- ref_classes
  num_ranges     <- ref_ranges
  factor_levels  <- ref_fac_levels
  numeric_like   <- c("numeric", "integer", "integer64")

  is_num <- names(var_classes)[var_classes %in% numeric_like]
  is_cat <- setdiff(predictor_vars, is_num)

  # ---- mixture validation ----
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) {
      if (is.null(grp$vars)) stop("Each mixture group must contain a 'vars' character vector.")
      if (!all(grp$vars %in% predictor_vars)) {
        missing <- setdiff(grp$vars, predictor_vars)
        stop("Mixture variables not in model predictors: ", paste(missing, collapse = ", "))
      }
      bad_mix <- setdiff(grp$vars, names(var_classes)[var_classes %in% numeric_like])
      if (length(bad_mix)) {
        stop("Mixture variables must be numeric-like. Non-numeric mixture vars: ",
             paste(bad_mix, collapse = ", "))
      }
      mixture_vars <- c(mixture_vars, grp$vars)
    }
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }

  # Truncated Dirichlet sampler
  .sample_trunc_dirichlet <- function(n, lower, upper, total,
                                      alpha = NULL, oversample = 4L, max_tries = 10000L) {
    k <- length(lower)
    if (length(upper) != k) stop("upper must have the same length as lower.")
    if (is.null(alpha)) alpha <- rep(1, k)
    min_sum <- sum(lower); max_sum <- sum(upper)
    if (total < min_sum - 1e-12 || total > max_sum + 1e-12)
      stop("Infeasible mixture constraints: need sum(lower) <= total <= sum(upper).")
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
    if (filled < n) stop("Could not sample enough feasible mixture points within max_tries.")
    res
  }

  # ---- sample non-mixture numerics ----
  nonmix_num <- setdiff(is_num, mixture_vars)
  T_num <- NULL
  if (length(nonmix_num)) {
    rng <- vapply(nonmix_num, function(v) {
      if (!is.null(num_ranges) && ncol(num_ranges) && v %in% colnames(num_ranges)) {
        r <- num_ranges[, v]
        if (!all(is.finite(r)) || r[1] >= r[2]) c(0, 1) else r
      } else c(0, 1)
    }, numeric(2))
    U <- as.matrix(lhs::maximinLHS(n, length(nonmix_num)))
    T_num <- matrix(NA_real_, nrow = n, ncol = length(nonmix_num))
    colnames(T_num) <- nonmix_num
    for (j in seq_along(nonmix_num)) {
      lo <- rng[1, j]; hi <- rng[2, j]
      T_num[, j] <- lo + (hi - lo) * U[, j]
    }
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
        stop("lower/upper must each have length equal to mixture vars: ", paste(vars, collapse = ", "))
      vals <- .sample_trunc_dirichlet(n, lower, upper, total)
      colnames(vals) <- vars
      T_mix[, vars] <- vals
    }
    T_mix <- as.data.frame(T_mix)
  }

  # ---- sample categoricals ----
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

  # ---- assemble predictors ----
  parts <- list(T_num, T_mix, T_cat)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) stop("No predictors could be sampled from the schema.")
  T_data <- do.call(cbind, parts)

  # Ensure all predictors present
  missing_pred <- setdiff(predictor_vars, colnames(T_data))
  if (length(missing_pred)) {
    for (v in missing_pred) {
      if (v %in% names(var_classes) && var_classes[[v]] %in% numeric_like) {
        T_data[[v]] <- runif(n)
      } else {
        lev <- factor_levels[[v]]
        if (is.null(lev) || !length(lev)) lev <- objects[[1]]$xlevels[[v]]
        if (is.null(lev) || !length(lev)) lev <- c("L1", "L2")
        T_data[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
      }
    }
  }

  # Order like predictor_vars
  T_data <- T_data[, predictor_vars, drop = FALSE]

  # ---- predict for each model on the shared T_data ----
  data_df <- T_data
  pred_df <- as.data.frame(matrix(NA_real_, nrow = nrow(data_df), ncol = 0))
  rownames(pred_df) <- NULL

  for (i in seq_along(objects)) {
    obj   <- objects[[i]]
    preds <- predict(obj, newdata = data_df, debias = debias)
    if (is.list(preds) && !is.null(preds$fit)) preds <- preds$fit
    resp  <- tryCatch(as.character(obj$formula[[2L]]), error = function(e) paste0("resp", i))
    colname <- resp
    if (colname %in% colnames(data_df)) colname <- paste0(colname, ".pred")
    pred_df[[colname]] <- as.numeric(preds)
  }

  rownames(data_df) <- NULL
  all_df <- cbind(data_df, pred_df)

  list(data = data_df, pred = pred_df, all = all_df)
}
