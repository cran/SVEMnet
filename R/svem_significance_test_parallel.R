.svem_wmt_set_rngkind <- function(sample_kind) {
  suppressWarnings(RNGkind(kind = "L'Ecuyer-CMRG", sample.kind = sample_kind))
  invisible(NULL)
}

.svem_wmt_seed_add <- function(seed, add) {
  if (is.null(seed)) return(NULL)
  out <- (as.double(seed) + as.double(add)) %% as.double(.Machine$integer.max)
  as.integer(ifelse(out <= 0, 1, out))
}

.svem_wmt_make_iter_seeds <- function(n, base_seed = NULL,
                                      sample_kind = "Rejection") {
  .svem_wmt_set_rngkind(sample_kind)
  if (!is.null(base_seed)) set.seed(as.integer(base_seed))
  sample.int(.Machine$integer.max, n)
}

.svem_wmt_worker_setup <- function(contrasts_opts, sample_kind, normal_kind) {
  suppressWarnings(
    RNGkind(
      kind = "L'Ecuyer-CMRG",
      normal.kind = normal_kind,
      sample.kind = sample_kind
    )
  )
  options(contrasts = contrasts_opts)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  } else {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
  }
  NULL
}

.svem_wmt_lapply <- function(tasks, fun, nCore, contrasts_opts, sample_kind,
                             normal_kind = "Inversion", ...) {
  if (nCore == 1L) {
    old_contrasts <- getOption("contrasts")
    on.exit(options(contrasts = old_contrasts), add = TRUE)
    options(contrasts = contrasts_opts)
    return(lapply(tasks, fun, ...))
  }

  cl <- parallel::makeCluster(nCore)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(
    cl, .svem_wmt_worker_setup,
    contrasts_opts = contrasts_opts,
    sample_kind = sample_kind,
    normal_kind = normal_kind
  )
  parallel::parLapplyLB(cl, tasks, fun, ...)
}

.svem_wmt_fit_task <- function(task, f_use, data_fit, resp_name, nBoot,
                               glmnet_alpha, weight_scheme, objective, relaxed,
                               dots, T_data, y_mean, rng_sample_kind,
                               rng_normal_kind = "Inversion",
                               fit_fun = NULL, predict_fun = NULL) {
  if (is.null(fit_fun)) fit_fun <- getExportedValue("SVEMnet", "SVEMnet")
  if (is.null(predict_fun)) predict_fun <- stats::predict

  suppressWarnings(
    RNGkind(
      kind = "L'Ecuyer-CMRG",
      normal.kind = rng_normal_kind,
      sample.kind = rng_sample_kind
    )
  )
  set.seed(task$seed)

  fit_data <- data_fit
  if (identical(task$type, "permutation")) {
    # Draw this once, before any fitting attempt, so retries hold the exact
    # permutation fixed and alter only the inner SVEM random weights.
    perm_idx <- sample.int(nrow(data_fit), replace = FALSE)
    fit_data[[resp_name]] <- data_fit[[resp_name]][perm_idx]
  }

  failures <- character(0L)
  for (attempt in seq_len(3L)) {
    # Attempt one retains the historical per-slot stream. Retry streams are
    # deterministic and schedule-independent.
    if (attempt > 1L) {
      retry_seed <- (
        as.double(task$seed) + 3000001 + (attempt - 2L) * 104729
      ) %% as.double(.Machine$integer.max)
      set.seed(as.integer(ifelse(retry_seed <= 0, 1, retry_seed)))
    }

    fit <- tryCatch(
      do.call(fit_fun, c(list(
        formula = f_use,
        data = fit_data,
        nBoot = nBoot,
        glmnet_alpha = glmnet_alpha,
        weight_scheme = weight_scheme,
        objective = objective,
        relaxed = relaxed
      ), dots)),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      failures <- c(failures, paste0("fit: ", conditionMessage(fit)))
      next
    }

    pred <- tryCatch(
      predict_fun(fit, newdata = T_data, debias = FALSE, se.fit = TRUE),
      error = function(e) e
    )
    if (inherits(pred, "error")) {
      failures <- c(failures, paste0("prediction: ", conditionMessage(pred)))
      next
    }

    fit_values <- as.numeric(pred$fit)
    member_sd <- as.numeric(pred$se.fit)
    if (length(fit_values) != nrow(T_data) ||
        length(member_sd) != nrow(T_data)) {
      failures <- c(failures, "prediction: unexpected output length")
      next
    }
    if (any(!is.finite(fit_values))) {
      failures <- c(failures, "prediction: non-finite fitted values")
      next
    }
    if (any(!is.finite(member_sd) | member_sd <= 0)) {
      failures <- c(failures, "prediction: nonpositive or non-finite member SD")
      next
    }

    surface <- (fit_values - y_mean) / member_sd
    if (any(!is.finite(surface))) {
      failures <- c(failures, "prediction: non-finite standardized surface")
      next
    }
    return(list(
      ok = TRUE,
      surface = surface,
      retries = attempt - 1L,
      retry_failures = failures
    ))
  }

  list(
    ok = FALSE,
    surface = NULL,
    retries = 2L,
    retry_failures = failures
  )
}

.svem_wmt_abort_failed_slots <- function(results, label) {
  ok <- vapply(results, `[[`, logical(1L), "ok")
  if (all(ok)) return(invisible(NULL))

  bad <- which(!ok)
  reasons <- unlist(lapply(results[bad], `[[`, "retry_failures"),
                    use.names = FALSE)
  reason_table <- sort(table(reasons), decreasing = TRUE)
  summary <- if (length(reason_table)) {
    paste0(names(reason_table), " (", as.integer(reason_table), "x)",
           collapse = "; ")
  } else {
    "unknown failure"
  }
  stop(
    length(bad), " of ", length(results), " requested ", label,
    " fit slot(s) failed after three attempts; no fits were dropped. ",
    "Failure summary: ", summary,
    call. = FALSE
  )
}

.svem_wmt_component_count <- function(evalues, percent) {
  cumulative <- cumsum(evalues) / sum(evalues) * 100
  selected <- which(cumulative > percent)[1L]
  if (is.na(selected)) length(evalues) else selected
}

.svem_wmt_reduce_surfaces <- function(M_Y, M_pi_Y, percent) {
  M_Y <- as.matrix(M_Y)
  M_pi_Y <- as.matrix(M_pi_Y)
  if (!nrow(M_Y) || !nrow(M_pi_Y) || ncol(M_Y) != ncol(M_pi_Y)) {
    stop("Original and permutation surface matrices have incompatible dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(M_Y)) || any(!is.finite(M_pi_Y))) {
    stop("WMT surface matrices must be finite.", call. = FALSE)
  }

  col_means <- colMeans(M_pi_Y)
  col_sds <- apply(M_pi_Y, 2L, stats::sd)
  keep <- is.finite(col_sds) & col_sds > 0
  removed <- which(!keep)
  if (length(removed)) {
    warning(
      "Removed ", length(removed), " zero-variance permutation-grid column(s): ",
      paste(removed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!any(keep)) {
    stop("All permutation-grid columns have zero or non-finite variance.",
         call. = FALSE)
  }

  M_Y <- M_Y[, keep, drop = FALSE]
  M_pi_Y <- M_pi_Y[, keep, drop = FALSE]
  col_means <- col_means[keep]
  col_sds <- col_sds[keep]
  tilde_M_pi_Y <- sweep(sweep(M_pi_Y, 2L, col_means, "-"),
                        2L, col_sds, "/")
  tilde_M_Y <- sweep(sweep(M_Y, 2L, col_means, "-"),
                     2L, col_sds, "/")

  svd_res <- svd(tilde_M_pi_Y)
  evalues_all <- (svd_res$d^2) / (nrow(tilde_M_pi_Y) - 1L)
  ev_sum <- sum(evalues_all)
  if (!length(evalues_all) || !is.finite(ev_sum) || ev_sum <= 0) {
    stop("Permutation surfaces have no positive SVD variance.", call. = FALSE)
  }

  # Karl (2024) rescales the correlation-PCA eigenvalues to sum to the
  # number of retained grid columns. This rescaling leaves the variance
  # fractions unchanged but is part of the published distance convention.
  p_cols <- ncol(tilde_M_pi_Y)
  evalues_all <- evalues_all / ev_sum * p_cols
  k_idx <- .svem_wmt_component_count(evalues_all, percent)
  k_idx <- max(1L, min(k_idx, ncol(svd_res$v)))
  evalues <- evalues_all[seq_len(k_idx)]
  evectors <- svd_res$v[, seq_len(k_idx), drop = FALSE]
  evalues[!is.finite(evalues) | evalues <= 0] <- 1e-12

  Z_pi <- tilde_M_pi_Y %*% evectors
  Z_Y <- tilde_M_Y %*% evectors
  d_pi_Y <- sqrt(rowSums((Z_pi^2) /
                           rep(evalues, each = nrow(Z_pi))))
  d_Y <- sqrt(rowSums((Z_Y^2) /
                        rep(evalues, each = nrow(Z_Y))))

  list(
    d_Y = d_Y,
    d_pi_Y = d_pi_Y,
    removed_grid_columns = removed,
    retained_grid_columns = which(keep),
    retained_components = k_idx,
    eigenvalues = evalues
  )
}

.svem_wmt_shasho_upper_tail <- function(q, mu, sigma, nu, tau) {
  if (any(!is.finite(c(q, mu, sigma, nu, tau))) ||
      any(sigma <= 0) || any(tau <= 0)) {
    stop("SHASHo tail inputs must be finite with positive sigma and tau.",
         call. = FALSE)
  }
  z <- (q - mu) / sigma
  transformed <- suppressWarnings(sinh(tau * asinh(z) - nu))
  p <- stats::pnorm(transformed, lower.tail = FALSE)
  pmin(1, pmax(0, p))
}

#' SVEM whole-model significance test with mixture support (parallel)
#'
#' Perform a permutation-based whole-model significance test for a continuous
#' (Gaussian) SVEM fit, with optional mixture-factor groups and parallel SVEM
#' refits.
#'
#' The procedure follows Karl (2024): it generates a space-filling evaluation
#' grid in the factor space, fits multiple SVEM models on the original data and
#' on permuted responses, standardizes grid predictions, reduces them via an
#' SVD-based low-rank representation, and summarizes each fit by a
#' Mahalanobis-type distance in the reduced space. A flexible SHASHo
#' distribution is then fit to the permutation distances and used to obtain a
#' whole-model \eqn{p}-value for the observed surface.
#'
#' Because the test is based on a finite number of permutations and a fitted
#' null distribution, the reported \eqn{p}-values are approximate and are
#' intended as a diagnostic measure of global factor signal, not as exact
#' hypothesis tests.
#'
#' With \code{nCore = 1}, SVEM refits run in-process. With more than one core,
#' the function uses a private PSOCK cluster and \code{parallel::parLapplyLB()}.
#' It does not register or change a \code{foreach} backend.
#'
#' Reproducible parallel RNG (Windows/macOS/Linux): when \code{seed} is supplied,
#' the function sets its private computation stream to
#' \code{RNGkind("L'Ecuyer-CMRG", sample.kind = rng_sample_kind)} and generates
#' a deterministic, per-iteration seed schedule on the master. Each fit slot
#' calls \code{set.seed()} with its assigned seed before performing random draws
#' (including permutations and the bootstrap randomness inside \code{SVEMnet()}).
#' This makes results reproducible regardless of worker scheduling and
#' independent of the number of cores. The caller's RNG kind and
#' \code{.Random.seed} are restored on exit, including error exits.
#' A failed fit slot is retried twice with deterministic sub-seeds while its
#' response permutation remains fixed. If all three attempts fail, the test
#' aborts with a failure summary rather than silently reducing \code{nSVEM} or
#' \code{nPerm}.
#' Likewise, a non-converged SHASHo null-distribution fit is reported as an
#' error rather than being used to calculate a whole-model \eqn{p}-value.
#'
#' The function can optionally reuse a deterministic, locked expansion built
#' with \code{bigexp_terms()}. Supply \code{spec} (and optionally
#' \code{response}) to ensure that categorical levels, contrasts, and the
#' polynomial/interaction structure are identical across repeated calls and
#' across multiple responses sharing the same factor space.
#'
#' Although the implementation calls \code{SVEMnet()} internally and will
#' technically run for any supported \code{family}, the significance test is
#' \emph{designed} for continuous (Gaussian) responses and should be interpreted
#' in that setting.
#'
#' @param formula A model formula. If \code{spec} is provided, the right-hand
#'   side is ignored and replaced by the locked expansion in \code{spec}.
#' @param data A data frame containing the variables in the model.
#' @param mixture_groups Optional list describing one or more mixture-factor
#'   groups. Each element should be a list with components:
#'   \itemize{
#'     \item \code{vars}: character vector of column names;
#'     \item \code{lower}: numeric vector of lower bounds (same length as \code{vars});
#'     \item \code{upper}: numeric vector of upper bounds (same length as \code{vars});
#'     \item \code{total}: scalar specifying the sum of the mixture variables.
#'   }
#'   All mixture variables must appear in exactly one group. Defaults to \code{NULL}.
#' @param nPoint Number of random evaluation points in the factor space
#'   (default \code{2000}); must be a positive integer.
#' @param nSVEM Number of SVEM fits on the original (unpermuted) data used to
#'   summarize the observed surface (default \code{10}); must be at least one.
#' @param nPerm Number of SVEM fits on permuted responses used to build the null
#'   reference distribution (default \code{150}); must be at least five.
#'   Values below 20 are permitted for tests but trigger a warning because they
#'   are too small for production inference.
#' @param percent Percentage of variance to capture in the SVD of the permutation
#'   surfaces (default \code{90}); must lie strictly between zero and 100.
#' @param nBoot Number of bootstrap iterations within each inner SVEM fit
#'   (default \code{100}); must be at least two because member standard
#'   deviations are required.
#' @param glmnet_alpha Numeric vector of \code{glmnet} alpha values
#'   (default \code{c(1)}).
#' @param weight_scheme Weighting scheme for the inner SVEM fits. Only
#'   \code{"SVEM"} is accepted: the whole-model test in Karl (2024) is defined
#'   for the anti-correlated fractional-random-weight scheme, so the other
#'   \code{SVEMnet()} weighting schemes are intentionally not offered here.
#' @param objective Objective used inside \code{SVEMnet()} to pick the bootstrap
#'   path solution. One of \code{"wAIC"}, \code{"wBIC"}, or
#'   \code{"wSSE"} (default \code{"wAIC"}).
#' @param relaxed Logical; default \code{FALSE}. When \code{TRUE}, inner
#'   \code{SVEMnet()} fits use \code{glmnet}'s relaxed elastic-net path and
#'   select both \code{lambda} and relaxed \code{gamma} on each bootstrap.
#'   When \code{FALSE}, the standard \code{glmnet} path is used. If
#'   \code{relaxed = TRUE} and \code{glmnet_alpha} includes \code{0}, ridge
#'   (\code{alpha = 0}) is dropped by \code{SVEMnet()} for relaxed fits.
#' @param verbose Logical; if \code{TRUE}, display progress messages
#'   (default \code{TRUE}).
#' @param nCore Number of CPU cores for processing. Default is \code{1}, which
#'   runs in-process and is suitable for CRAN and other shared systems. Values
#'   greater than one use a private cluster local to this call.
#' @param seed Optional integer seed for reproducible RNG (default \code{NULL}).
#'   Deterministic per-fit seeds yield reproducibility regardless of parallel
#'   scheduling and core count.
#' @param rng_sample_kind Sampling algorithm used with the
#'   \code{"L'Ecuyer-CMRG"} RNG. The default \code{"Rejection"} uses R's modern
#'   unbiased discrete sampler. Use \code{"Rounding"} to reproduce legacy WMT
#'   streams from SVEMnet 3.3.1 and earlier.
#' @param spec Optional \code{bigexp_spec} created by \code{bigexp_terms()}.
#'   If provided, the test reuses its locked expansion. The working formula
#'   becomes \code{bigexp_formula(spec, response_name)}, where
#'   \code{response_name} is taken from \code{response} if supplied, otherwise
#'   from the left-hand side of \code{formula}. Categorical sampling uses
#'   \code{spec$levels}, and numeric sampling prefers \code{spec$num_range}
#'   when available.
#'   Discrete numeric predictors recorded by \code{bigexp_terms()}
#'   (\code{spec$settings$discrete_numeric} + \code{spec$settings$discrete_levels})
#'   are sampled only from their recorded allowed levels when building the
#'   evaluation grid.
#'   Blocking variables recorded in \code{spec$settings$blocking} are treated
#'   as nuisance factors and held fixed at a single reference value across the
#'   evaluation grid (numeric blocks at the range midpoint, snapped to the
#'   discrete support when recorded; categorical blocks at the most frequent
#'   training level), matching \code{svem_random_table_multi()}. They are
#'   never swept, so between-block variation does not contribute to the
#'   whole-model distances.
#' @param response Optional character name for the response variable to use when
#'   \code{spec} is supplied. If omitted, the response is taken from the
#'   left-hand side of \code{formula}.
#' @param use_spec_contrasts Logical; default \code{TRUE}. When \code{spec} is
#'   supplied and \code{use_spec_contrasts = TRUE}, the function replays
#'   \code{spec$settings$contrasts_options} on the parallel workers for
#'   deterministic factor coding.
#' @param ... Additional arguments passed to \code{SVEMnet()} and then to
#'   \code{glmnet()} (for example: \code{penalty.factor},
#'   \code{lower.limits}, \code{upper.limits}, \code{standardize.response}, etc.).
#'   The \code{relaxed} setting is controlled by the \code{relaxed} argument of
#'   this function and any \code{relaxed} value passed via \code{...} is ignored
#'   with a warning. User \code{weights} and \code{offset} are not supported
#'   (SVEM controls its own weighting). Argument names that neither
#'   \code{SVEMnet()} nor the installed \code{glmnet} recognizes are ignored
#'   with a warning (misspelling protection).
#'
#' @return An object of class \code{"svem_significance_test"}, a list with
#'   components:
#'   \itemize{
#'     \item \code{p_value}: median whole-model \eqn{p}-value across the
#'       \code{nSVEM} original SVEM fits.
#'     \item \code{p_values}: numeric vector of length \code{nSVEM} with the
#'       per-fit \eqn{p}-values.
#'     \item \code{d_Y}: numeric vector of distances for the original SVEM fits.
#'     \item \code{d_pi_Y}: numeric vector of distances for the permutation fits.
#'     \item \code{distribution_fit}: fitted SHASHo distribution object.
#'     \item \code{data_d}: data frame of distances and source labels
#'       (original vs permutation), suitable for plotting.
#'     \item \code{diagnostics}: RNG mode, deterministic retry counts and
#'       reasons, SHASHo convergence information, removed/retained
#'       evaluation-grid columns, and the retained SVD component count.
#'   }
#'
#' @seealso \code{\link{SVEMnet}}, \code{\link{bigexp_terms}},
#'   \code{\link{bigexp_formula}}
#' @template ref-svem
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist SHASHo
#' @importFrom stats model.frame model.response model.matrix delete.response terms
#' @importFrom stats median complete.cases rgamma coef predict sd
#' @importFrom parallel makeCluster stopCluster clusterCall parLapplyLB
#'
#' @examples
#' \dontrun{
#'   ## Production WMT inference is intentionally not run in examples because
#'   ## it requires repeated ensemble fits for many response permutations.
#'   set.seed(1)
#'
#'   # Small toy data with a 3-component mixture A, B, C
#'   n <- 40
#'   sample_trunc_dirichlet <- function(n, lower, upper, total) {
#'     k <- length(lower)
#'     stopifnot(length(upper) == k, total >= sum(lower), total <= sum(upper))
#'     avail <- total - sum(lower)
#'     if (avail <= 0) return(matrix(rep(lower, each = n), nrow = n))
#'     out <- matrix(NA_real_, n, k)
#'     i <- 1L
#'     while (i <= n) {
#'       g <- rgamma(k, 1, 1)
#'       w <- g / sum(g)
#'       x <- lower + avail * w
#'       if (all(x <= upper + 1e-12)) { out[i, ] <- x; i <- i + 1L }
#'     }
#'     out
#'   }
#'
#'   lower <- c(0.10, 0.20, 0.05)
#'   upper <- c(0.60, 0.70, 0.50)
#'   total <- 1.0
#'   ABC   <- sample_trunc_dirichlet(n, lower, upper, total)
#'   A <- ABC[, 1]; B <- ABC[, 2]; C <- ABC[, 3]
#'   X <- runif(n)
#'   F <- factor(sample(c("red", "blue"), n, replace = TRUE))
#'   y <- 2 + 3*A + 1.5*B + 1.2*C + 0.5*X + 1*(F == "red") + rnorm(n, sd = 0.3)
#'   dat <- data.frame(y = y, A = A, B = B, C = C, X = X, F = F)
#'
#'   mix_spec <- list(list(
#'     vars  = c("A", "B", "C"),
#'     lower = lower,
#'     upper = upper,
#'     total = total
#'   ))
#'
#'   ## Example 1: direct formula interface (no locked expansion spec)
#'   res1 <- svem_significance_test_parallel(
#'     y ~ A + B + C + X + F,
#'     data           = dat,
#'     mixture_groups = mix_spec,
#'     glmnet_alpha   = 1,
#'     weight_scheme  = "SVEM",
#'     objective      = "auto",
#'     relaxed        = FALSE,   # default, shown for clarity
#'     nCore          = 1,
#'     seed           = 123,
#'     verbose        = FALSE
#'   )
#'   res1$p_value
#'
#'   ## Example 2: using a deterministic bigexp expansion spec
#'   ## Build a wide expansion once and reuse it via `spec`
#'   spec <- bigexp_terms(
#'     y ~ A + B + C + X + F,
#'     data             = dat,
#'     factorial_order  = 2,  # up to 2-way interactions
#'     polynomial_order = 2   # up to quadratic terms in continuous vars
#'   )
#'
#'   ## Run the same significance test, but with the locked expansion:
#'   ## - `formula` is still required, but its RHS is ignored when `spec` is given
#'   ## - `response` tells the helper which LHS to use with `spec`
#'   res2 <- svem_significance_test_parallel(
#'     y ~ A + B + C + X + F,
#'     data               = dat,
#'     mixture_groups     = mix_spec,
#'     glmnet_alpha       = 1,
#'     weight_scheme      = "SVEM",
#'     objective          = "auto",
#'     relaxed            = FALSE,
#'     nCore              = 1,
#'     seed               = 123,
#'     spec               = spec,
#'     response           = "y",
#'     use_spec_contrasts = TRUE,
#'     verbose            = FALSE
#'   )
#'   res2$p_value
#' }
#' @export
svem_significance_test_parallel <- function(
    formula, data, mixture_groups = NULL,
    nPoint = 2000, nSVEM = 10, nPerm = 150,
    percent = 90, nBoot = 100,
    glmnet_alpha = c(1),
    weight_scheme = c("SVEM"),
    objective = c("wAIC", "wBIC", "wSSE","auto"),
    relaxed = FALSE,
    verbose = TRUE,
    nCore = 1L,
    seed = NULL,
    rng_sample_kind = c("Rejection", "Rounding"),
    spec = NULL,
    response = NULL,
    use_spec_contrasts = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- basic choices ---
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)
  rng_sample_kind <- match.arg(rng_sample_kind)
  nPoint <- .svem_integer_scalar(nPoint, "nPoint", min = 1L)
  nSVEM <- .svem_integer_scalar(nSVEM, "nSVEM", min = 1L)
  nPerm <- .svem_integer_scalar(nPerm, "nPerm", min = 5L)
  nBoot <- .svem_integer_scalar(nBoot, "nBoot", min = 2L)
  nCore <- .svem_integer_scalar(nCore, "nCore", min = 1L)
  seed <- .svem_integer_scalar(seed, "seed", min = 0L, allow_null = TRUE)
  percent <- .svem_numeric_scalar(
    percent, "percent", lower = 0, upper = 100,
    lower_open = TRUE, upper_open = TRUE
  )
  relaxed <- .svem_logical_scalar(relaxed, "relaxed")
  verbose <- .svem_logical_scalar(verbose, "verbose")
  use_spec_contrasts <- .svem_logical_scalar(
    use_spec_contrasts, "use_spec_contrasts"
  )
  if (nPerm < 20L) {
    warning(
      "`nPerm < 20` is intended only for tests; use the default of 150 ",
      "for production WMT inference.",
      call. = FALSE
    )
  }
  data <- as.data.frame(data)

  # FIX: "auto" objective must assign a real objective
  if (identical(objective, "auto")) {
    objective <- "wAIC"
    if (isTRUE(verbose)) message("objective='auto' -> using 'wAIC'.")
  }

  # Robust single-variable response name; as.character(formula[[2]]) is
  # length > 1 for transformed responses such as log(y) and would crash
  # scalar contexts below on R >= 4.2. The permutation scheme permutes the
  # raw response column, so exactly one response variable is required.
  .resp_var_from_formula <- function(f) {
    rv <- .svem_response_vars(f)
    if (length(rv) != 1L) {
      stop(
        "svem_significance_test_parallel() requires a formula whose ",
        "left-hand side contains exactly one response variable (e.g. y ~ ... ",
        "or log(y) ~ ...); found: ",
        if (length(rv)) paste(rv, collapse = ", ") else "none", "."
      )
    }
    rv
  }

  # Determine response and working formula (optionally from spec)
  if (!is.null(spec)) {
    resp_name <- if (!is.null(response)) {
      if (!is.character(response) || length(response) != 1L || !nzchar(response))
        stop("response must be a non-empty character scalar when provided.")
      response
    } else {
      .resp_var_from_formula(formula)
    }
    if (!resp_name %in% names(data)) stop("Response '", resp_name, "' not found in 'data'.")
    f_use <- bigexp_formula(spec, resp_name)
  } else {
    f_use <- formula
    resp_name <- .resp_var_from_formula(formula)
  }

  wmt_terms <- stats::terms(f_use, data = data)
  if (!identical(attr(wmt_terms, "intercept"), 1L) ||
      (!is.null(spec) && identical(spec$settings$intercept, FALSE))) {
    stop(
      "svem_significance_test_parallel() follows the published intercept convention; ",
      "formulas using `~ 0`, `- 1`, or `bigexp_terms(intercept = FALSE)` are not supported.",
      call. = FALSE
    )
  }

  # Choose contrasts options (spec's, if requested)
  contrasts_opts <- if (!is.null(spec) && isTRUE(use_spec_contrasts)) {
    spec$settings$contrasts_options %||% getOption("contrasts")
  } else {
    getOption("contrasts")
  }

  # --- private master RNG for serial prep (design, mixtures, factors) ---
  oldKinds <- RNGkind()
  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_random_seed <- if (had_random_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    suppressWarnings(RNGkind(
      kind = oldKinds[1L],
      normal.kind = oldKinds[2L],
      sample.kind = oldKinds[3L]
    ))
    if (had_random_seed) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  .svem_wmt_set_rngkind(rng_sample_kind)
  rng_normal_kind <- RNGkind()[2L]
  if (!is.null(seed)) set.seed(seed)

  # --- enforce single-threaded BLAS/OpenMP BEFORE spawning workers -------------
  thread_vars <- c(
    "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  old_thread_env <- Sys.getenv(thread_vars, unset = NA_character_)
  on.exit({
    for (v in thread_vars) {
      oldv <- old_thread_env[[v]]
      if (is.na(oldv) || !nzchar(oldv)) {
        Sys.unsetenv(v)
      } else {
        do.call(Sys.setenv, setNames(list(oldv), v))
      }
    }
  }, add = TRUE)

  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))

  # Sanitize ... so explicit 'relaxed' here cannot be overridden
  dots <- list(...)

  ## --- guard against binomial family: only Gaussian is supported here ----
  if ("family" %in% names(dots)) {
    fam <- dots$family
    fam_name <- tryCatch({
      if (inherits(fam, "family")) {
        fam$family
      } else if (is.function(fam)) {
        fam()$family
      } else {
        as.character(fam)[1L]
      }
    }, error = function(e) NA_character_)

    if (!is.na(fam_name) && tolower(fam_name) == "binomial") {
      stop("svem_significance_test_parallel() is designed for continuous (Gaussian) responses; ",
           "'binomial' family is not supported.")
    }
  }

  if ("relaxed" %in% names(dots)) {
    warning("Ignoring 'relaxed' in '...'; use the 'relaxed' argument of svem_significance_test_parallel().")
    dots$relaxed <- NULL
  }

  if ("complexity" %in% names(dots)) {
    warning("Ignoring 'complexity' in '...'; the whole-model test always uses ",
            "the support-count complexity measure of Karl (2024).")
    dots$complexity <- NULL
  }

  # Warn here (on the master) about arguments SVEMnet() would drop on the
  # workers, where warnings are invisible.
  if ("weights" %in% names(dots)) {
    warning("Ignoring user 'weights'; SVEM uses its own bootstrap weights.")
    dots$weights <- NULL
  }
  if ("offset" %in% names(dots)) {
    warning("Ignoring user 'offset'; SVEMnet does not support offsets.")
    dots$offset <- NULL
  }

  # Validate remaining '...' names on the master before spawning workers:
  # SVEMnet()/glmnet() warnings raised in PSOCK workers are not visible on the
  # master, so a misspelled argument would otherwise vanish silently.
  unknown <- .svem_unknown_glmnet_args(
    dots,
    funs = list(SVEMnet, glmnet::glmnet),
    context = "SVEMnet/glmnet via svem_significance_test_parallel"
  )
  if (length(unknown)) dots <- dots[setdiff(names(dots), unknown)]

  # ---------------------------------------------------------------------------
  # Training design pieces (ROBUST TO NA IN RESPONSE/PREDICTORS):
  # - build model.frame with na.pass so we can explicitly filter
  # - fit/permutation data are restricted to complete cases for the model
  # ---------------------------------------------------------------------------
  mf0 <- stats::model.frame(f_use, data, na.action = stats::na.pass)

  keep <- stats::complete.cases(mf0)
  n_drop <- sum(!keep)

  if (n_drop > 0L && isTRUE(verbose)) {
    message("Dropping ", n_drop, " row(s) with NA in model variables before fitting/permuting.")
  }
  if (all(!keep)) stop("No complete cases available after removing NA rows in model variables.")

  data_fit <- data[keep, , drop = FALSE]
  mf <- mf0[keep, , drop = FALSE]
  y  <- stats::model.response(mf)

  if (length(y) < 2L) stop("Not enough complete cases to run the significance test.")
  if (is.matrix(y) || is.array(y) || !is.numeric(y)) {
    stop("The whole-model test requires one numeric continuous response.",
         call. = FALSE)
  }
  if (any(!is.finite(y))) {
    stop("The whole-model-test response must contain only finite values after NA removal.",
         call. = FALSE)
  }

  # For sampling, decide which raw columns are categorical vs continuous
  if (!is.null(spec)) {
    predictor_vars    <- spec$vars
    is_cat            <- spec$is_cat
    categorical_vars  <- names(is_cat)[is_cat]
    continuous_vars   <- names(is_cat)[!is_cat]
    disc_vars         <- spec$settings$discrete_numeric %||% character(0L)
    disc_levels       <- spec$settings$discrete_levels  %||% list()
    # Blocking variables are nuisance factors: hold them fixed at a single
    # reference value in the evaluation grid (matching svem_random_table_multi)
    # instead of sweeping them, so between-block variation does not inflate
    # the observed surface distances relative to the permutation surfaces.
    blocking_vars     <- intersect(as.character(spec$settings$blocking %||% character(0L)),
                                   predictor_vars)
    categorical_vars  <- setdiff(categorical_vars, blocking_vars)
    continuous_vars   <- setdiff(continuous_vars, blocking_vars)
  } else {
    predictor_vars   <- base::all.vars(stats::delete.response(stats::terms(f_use, data = data_fit)))
    # is.factor() (not class(z)[1L]) so ordered factors are categorical too;
    # class(ordered)[1L] is "ordered", which previously fell through to the
    # continuous sampler and produced all-NA predictions.
    is_cat_col       <- vapply(data_fit[predictor_vars], function(z) {
      is.factor(z) || is.character(z) || is.logical(z)
    }, logical(1L))
    categorical_vars <- predictor_vars[is_cat_col]
    continuous_vars  <- setdiff(predictor_vars, categorical_vars)
    disc_vars        <- character(0L)
    disc_levels      <- list()
    blocking_vars    <- character(0L)
  }

  # Mixture bookkeeping and validation (mirrors svem_random_table_multi):
  # every mixture variable must be a continuous numeric predictor of the
  # model, otherwise the evaluation grid silently leaves the simplex.
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) {
      if (is.null(grp$vars) || !is.character(grp$vars) || !length(grp$vars)) {
        stop("Each mixture group must contain a nonempty 'vars' character vector.")
      }
      missing_mix <- setdiff(grp$vars, predictor_vars)
      if (length(missing_mix)) {
        stop("Mixture variables not in model predictors: ",
             paste(missing_mix, collapse = ", "))
      }
      if (length(intersect(grp$vars, blocking_vars))) {
        stop("Mixture variables cannot be blocking variables. Offending vars: ",
             paste(intersect(grp$vars, blocking_vars), collapse = ", "))
      }
      bad_mix <- setdiff(grp$vars, continuous_vars)
      if (length(bad_mix)) {
        stop("Mixture variables must be continuous numeric predictors. ",
             "Offending vars: ", paste(bad_mix, collapse = ", "))
      }
      mixture_vars <- c(mixture_vars, grp$vars)
    }
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }

  # Disallow discrete numeric predictors as mixture variables (matches scoring behavior)
  if (length(disc_vars) && length(intersect(mixture_vars, disc_vars))) {
    stop(
      "Discrete numeric predictors cannot be used as mixture variables. Offenders: ",
      paste(intersect(mixture_vars, disc_vars), collapse = ", ")
    )
  }

  nonmix_continuous_vars <- setdiff(continuous_vars, mixture_vars)

  # Non-mixture continuous via maximin LHS over ranges; but respect discrete numeric support
  T_continuous <- NULL
  if (length(nonmix_continuous_vars) > 0) {
    disc_nonmix <- intersect(nonmix_continuous_vars, disc_vars)
    cont_nonmix <- setdiff(nonmix_continuous_vars, disc_nonmix)

    parts_cont <- list()

    # 1) Truly-continuous vars: LHS over ranges (prefer spec$num_range)
    if (length(cont_nonmix) > 0) {
      if (!is.null(spec) && !is.null(spec$num_range) && ncol(spec$num_range) > 0) {
        rng_mat <- spec$num_range[, colnames(spec$num_range) %in% cont_nonmix, drop = FALSE]
        missing <- setdiff(cont_nonmix, colnames(rng_mat))
        if (length(missing)) {
          add <- sapply(data_fit[missing], function(col) range(col, na.rm = TRUE))
          rownames(add) <- c("min","max")
          rng_mat <- cbind(rng_mat, add)
        }
        rng_mat <- rng_mat[, cont_nonmix, drop = FALSE]
      } else {
        rng_mat <- sapply(data_fit[cont_nonmix], function(col) range(col, na.rm = TRUE))
        rownames(rng_mat) <- c("min","max")
      }

      T_continuous_raw <- as.matrix(lhs::maximinLHS(nPoint, length(cont_nonmix)))
      T_cont <- matrix(NA_real_, nrow = nPoint, ncol = length(cont_nonmix))
      colnames(T_cont) <- cont_nonmix
      for (i in seq_along(cont_nonmix)) {
        lo <- rng_mat["min", i]; hi <- rng_mat["max", i]
        T_cont[, i] <- T_continuous_raw[, i] * (hi - lo) + lo
      }
      parts_cont[[length(parts_cont) + 1L]] <- as.data.frame(T_cont)
    }

    # 2) Discrete numeric vars: sample from allowed levels
    if (length(disc_nonmix) > 0) {
      Td <- lapply(disc_nonmix, function(v) {
        lv <- disc_levels[[v]]
        if (is.null(lv) || !length(lv)) {
          # fallback: infer from data_fit if spec didn't record levels for some reason
          x <- as.numeric(data_fit[[v]])
          lv <- sort(unique(x[is.finite(x)]))
        }
        lv <- as.numeric(lv)
        lv <- lv[is.finite(lv)]
        if (length(lv) < 1L) stop("No finite discrete levels available for '", v, "'.")
        # index-based sampling avoids sample()'s 1:x behavior for length-one supports
        lv[sample.int(length(lv), nPoint, replace = TRUE)]
      })
      Td <- as.data.frame(Td, check.names = FALSE, optional = FALSE)
      names(Td) <- disc_nonmix
      parts_cont[[length(parts_cont) + 1L]] <- Td
    }

    T_continuous <- do.call(cbind, parts_cont)
    T_continuous <- as.data.frame(T_continuous)
    # enforce the original column order
    T_continuous <- T_continuous[, nonmix_continuous_vars, drop = FALSE]
  }

  # Truncated Dirichlet sampler for mixture groups
  .sample_trunc_dirichlet <- function(n, lower, upper, total,
                                      alpha = NULL, oversample = 4L, max_tries = 10000L) {
    k <- length(lower)
    if (length(upper) != k) stop("upper must have the same length as lower.")
    if (is.null(alpha)) alpha <- rep(1, k)

    min_sum <- sum(lower); max_sum <- sum(upper)
    if (total < min_sum - 1e-12 || total > max_sum + 1e-12) {
      stop("Infeasible mixture constraints: need sum(lower) <= total <= sum(upper).")
    }

    avail <- total - min_sum
    if (avail <= 1e-12) {
      return(matrix(rep(lower, each = n), nrow = n))
    }

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
           "Try relaxing upper bounds or increasing 'oversample'/'max_tries'.")
    }
    res
  }

  # Mixture evaluation points
  T_mixture <- NULL
  if (!is.null(mixture_groups)) {
    mix_all_vars <- unlist(lapply(mixture_groups, `[[`, "vars"))
    T_mixture <- matrix(NA_real_, nrow = nPoint, ncol = length(mix_all_vars))
    colnames(T_mixture) <- mix_all_vars

    for (grp in mixture_groups) {
      vars  <- grp$vars
      k     <- length(vars)
      lower <- if (!is.null(grp$lower)) grp$lower else rep(0, k)
      upper <- if (!is.null(grp$upper)) grp$upper else rep(1, k)
      total <- if (!is.null(grp$total)) grp$total else 1

      if (length(lower) != k || length(upper) != k) {
        stop("lower and upper must each have length equal to the number of mixture variables (",
             paste(vars, collapse = ","), ").")
      }

      vals <- .sample_trunc_dirichlet(nPoint, lower, upper, total)
      colnames(vals) <- vars
      T_mixture[, vars] <- vals
    }
    T_mixture <- as.data.frame(T_mixture)
  }

  # Categorical sampling
  T_categorical <- NULL
  if (length(categorical_vars) > 0) {
    T_categorical <- vector("list", length(categorical_vars))
    names(T_categorical) <- categorical_vars
    for (v in categorical_vars) {
      # The inner fits see data_fit's raw classes, so the grid must reproduce
      # them: an ordered training factor requires an ordered grid column or
      # the predict-time class validator rejects the evaluation points.
      ord <- is.ordered(data_fit[[v]])
      if (!is.null(spec)) {
        lv <- spec$levels[[v]]
        if (is.null(lv)) lv <- sort(unique(as.character(data_fit[[v]])))
        T_categorical[[v]] <- factor(sample(lv, nPoint, replace = TRUE),
                                     levels = lv, ordered = ord)
      } else {
        x <- data_fit[[v]]
        if (is.factor(x)) {
          obs_lev <- levels(base::droplevels(x))
          T_categorical[[v]] <- factor(
            sample(obs_lev, nPoint, replace = TRUE),
            levels = levels(x), ordered = ord
          )
        } else {
          obs_lev <- sort(unique(as.character(x)))
          T_categorical[[v]] <- factor(
            sample(obs_lev, nPoint, replace = TRUE),
            levels = obs_lev
          )
        }
      }
      if (is.logical(data_fit[[v]])) {
        T_categorical[[v]] <- as.logical(as.character(T_categorical[[v]]))
      }
    }
    T_categorical <- as.data.frame(T_categorical, stringsAsFactors = FALSE)
  }

  # Assemble evaluation grid
  parts <- list(T_continuous, T_mixture, T_categorical)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0 && !length(blocking_vars)) stop("No predictors provided.")
  T_data <- if (length(parts)) {
    do.call(cbind, parts)
  } else {
    as.data.frame(matrix(nrow = nPoint, ncol = 0))
  }

  # Pin blocking variables to a single reference value across the grid,
  # mirroring svem_random_table_multi(): numeric blocks at the range midpoint
  # (snapped to the discrete support when recorded), categorical blocks at the
  # most frequent training level (ties broken by level order).
  if (length(blocking_vars)) {
    for (v in blocking_vars) {
      if (isTRUE(spec$is_cat[[v]])) {
        lv <- spec$levels[[v]]
        if (is.null(lv) || !length(lv)) lv <- sort(unique(as.character(data_fit[[v]])))
        tab <- table(factor(as.character(data_fit[[v]]), levels = lv))
        mode_val <- if (length(tab) && any(tab > 0)) names(tab)[which.max(tab)] else lv[1L]
        T_data[[v]] <- if (is.logical(data_fit[[v]])) {
          rep(as.logical(mode_val), nPoint)
        } else {
          factor(rep(mode_val, nPoint), levels = lv,
                 ordered = is.ordered(data_fit[[v]]))
        }
      } else {
        r <- if (!is.null(spec$num_range) && v %in% colnames(spec$num_range)) {
          spec$num_range[, v]
        } else {
          range(as.numeric(data_fit[[v]]), na.rm = TRUE)
        }
        val <- (as.numeric(r[1L]) + as.numeric(r[2L])) / 2
        if (v %in% disc_vars) {
          lv <- as.numeric(disc_levels[[v]])
          lv <- lv[is.finite(lv)]
          if (length(lv)) val <- lv[which.min(abs(lv - val))]
        }
        T_data[[v]] <- rep(val, nPoint)
      }
    }
  }

  y_mean <- mean(y, na.rm = TRUE)

  # --- Per-iteration seeds (schedule-independent reproducibility) --------------
  # Use separate namespaces so permutations don't change when nSVEM changes.
  seed_Y_base <- .svem_wmt_seed_add(seed, 1000001L)
  seed_piY_base <- .svem_wmt_seed_add(seed, 2000001L)
  iter_seeds_Y <- .svem_wmt_make_iter_seeds(
    nSVEM, base_seed = seed_Y_base, sample_kind = rng_sample_kind
  )
  iter_seeds_piY <- .svem_wmt_make_iter_seeds(
    nPerm, base_seed = seed_piY_base, sample_kind = rng_sample_kind
  )

  tasks <- c(
    lapply(seq_len(nSVEM), function(i) {
      list(type = "original", index = i, seed = iter_seeds_Y[i])
    }),
    lapply(seq_len(nPerm), function(i) {
      list(type = "permutation", index = i, seed = iter_seeds_piY[i])
    })
  )
  if (isTRUE(verbose)) {
    mode <- if (nCore == 1L) "in-process" else paste0("on ", nCore, " workers")
    message(
      "Fitting ", nSVEM, " original and ", nPerm,
      " permutation SVEM model(s) ", mode, "."
    )
  }
  fit_results <- .svem_wmt_lapply(
    tasks = tasks,
    fun = .svem_wmt_fit_task,
    nCore = nCore,
    contrasts_opts = contrasts_opts,
    sample_kind = rng_sample_kind,
    normal_kind = rng_normal_kind,
    f_use = f_use,
    data_fit = data_fit,
    resp_name = resp_name,
    nBoot = nBoot,
    glmnet_alpha = glmnet_alpha,
    weight_scheme = weight_scheme,
    objective = objective,
    relaxed = relaxed,
    dots = dots,
    T_data = T_data,
    y_mean = y_mean,
    rng_sample_kind = rng_sample_kind,
    rng_normal_kind = rng_normal_kind
  )
  original_results <- fit_results[seq_len(nSVEM)]
  permutation_results <- fit_results[nSVEM + seq_len(nPerm)]
  .svem_wmt_abort_failed_slots(original_results, "original-data")
  .svem_wmt_abort_failed_slots(permutation_results, "permutation")

  M_Y <- do.call(rbind, lapply(original_results, `[[`, "surface"))
  M_pi_Y <- do.call(rbind, lapply(permutation_results, `[[`, "surface"))
  reduction <- .svem_wmt_reduce_surfaces(M_Y, M_pi_Y, percent)
  d_Y <- reduction$d_Y
  d_pi_Y <- reduction$d_pi_Y

  # SHASHo fit
  shasho_warnings <- character(0L)
  shasho_error <- NULL
  distribution_fit <- suppressMessages(
    withCallingHandlers(
      tryCatch(
        gamlss::gamlss(
          d_pi_Y ~ 1,
          family = gamlss.dist::SHASHo(
            mu.link = "identity", sigma.link = "log",
            nu.link = "identity", tau.link = "log"
          ),
          control = gamlss::gamlss.control(n.cyc = 1000, trace = FALSE)
        ),
        error = function(e) {
          shasho_error <<- conditionMessage(e)
          NULL
        }
      ),
      warning = function(w) {
        shasho_warnings <<- c(shasho_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  shasho_warnings <- unique(shasho_warnings)
  warning_text <- if (length(shasho_warnings)) {
    paste0(" Warning(s): ", paste(shasho_warnings, collapse = "; "))
  } else {
    ""
  }
  if (is.null(distribution_fit)) {
    stop(
      "Failed to fit the SHASHo permutation-distance distribution",
      if (!is.null(shasho_error)) paste0(": ", shasho_error) else ".",
      warning_text,
      call. = FALSE
    )
  }
  if (!isTRUE(distribution_fit$converged)) {
    stop(
      "The SHASHo permutation-distance fit did not converge; no WMT p-value ",
      "was returned.", warning_text,
      call. = FALSE
    )
  }
  if (length(shasho_warnings)) {
    warning(
      "The converged SHASHo permutation-distance fit reported warning(s): ",
      paste(shasho_warnings, collapse = "; "),
      call. = FALSE
    )
  }

  mu    <- as.numeric(stats::coef(distribution_fit, what = "mu"))
  sigma <- exp(as.numeric(stats::coef(distribution_fit, what = "sigma")))
  nu    <- as.numeric(stats::coef(distribution_fit, what = "nu"))
  tau   <- exp(as.numeric(stats::coef(distribution_fit, what = "tau")))

  # pSHASHo(lower.tail = FALSE) currently computes 1 - CDF internally, which
  # loses extreme right-tail probabilities to cancellation. The transformed
  # normal survival probability is algebraically identical to SHASHo's CDF.
  p_values <- .svem_wmt_shasho_upper_tail(
    d_Y, mu = mu, sigma = sigma, nu = nu, tau = tau
  )
  p_value  <- stats::median(p_values)

  data_d <- data.frame(
    D = c(d_Y, d_pi_Y),
    Source_Type = c(rep("Original", length(d_Y)), rep("Permutation", length(d_pi_Y))),
    Response = resp_name
  )

  retry_failures <- unlist(
    lapply(fit_results, `[[`, "retry_failures"),
    use.names = FALSE
  )
  retry_failure_counts <- if (length(retry_failures)) {
    sort(table(retry_failures), decreasing = TRUE)
  } else {
    structure(integer(0L), names = character(0L))
  }
  retry_counts_Y <- vapply(original_results, `[[`, integer(1L), "retries")
  retry_counts_piY <- vapply(
    permutation_results, `[[`, integer(1L), "retries"
  )

  results_list <- list(
    p_value = p_value,
    p_values = p_values,
    d_Y = d_Y,
    d_pi_Y = d_pi_Y,
    distribution_fit = distribution_fit,
    data_d = data_d,
    diagnostics = list(
      rng_kind = "L'Ecuyer-CMRG",
      rng_normal_kind = rng_normal_kind,
      rng_sample_kind = rng_sample_kind,
      retry_counts = list(
        original = retry_counts_Y,
        permutation = retry_counts_piY,
        total = sum(retry_counts_Y) + sum(retry_counts_piY)
      ),
      retry_failure_reasons = retry_failure_counts,
      shasho_converged = TRUE,
      shasho_warnings = shasho_warnings,
      removed_grid_columns = reduction$removed_grid_columns,
      retained_grid_columns = reduction$retained_grid_columns,
      retained_components = reduction$retained_components
    )
  )
  class(results_list) <- "svem_significance_test"
  results_list
}
