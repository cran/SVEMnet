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
#' All SVEM refits (for the original and permuted responses) are run in
#' parallel using \code{foreach} + \code{doParallel}. Random draws (including
#' permutations and evaluation-grid sampling) are made reproducible across
#' workers using \code{doRNG} together with
#' \code{RNGkind("L'Ecuyer-CMRG", sample.kind = "Rounding")} when a
#' \code{seed} is supplied.
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
#'   (default \code{2000}).
#' @param nSVEM Number of SVEM fits on the original (unpermuted) data used to
#'   summarize the observed surface (default \code{10}).
#' @param nPerm Number of SVEM fits on permuted responses used to build the null
#'   reference distribution (default \code{150}).
#' @param percent Percentage of variance to capture in the SVD of the permutation
#'   surfaces (default \code{90}).
#' @param nBoot Number of bootstrap iterations within each inner SVEM fit
#'   (default \code{100}).
#' @param glmnet_alpha Numeric vector of \code{glmnet} alpha values
#'   (default \code{c(1)}).
#' @param weight_scheme Weighting scheme for SVEM (default \code{"SVEM"}).
#'   Passed to \code{SVEMnet()}.
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
#' @param nCore Number of CPU cores for parallel processing. Default is
#'   \code{parallel::detectCores() - 2}, with a floor of \code{1}.
#' @param seed Optional integer seed for reproducible parallel RNG (default
#'   \code{NULL}). When supplied, the master RNG kind is set to
#'   \code{"L'Ecuyer-CMRG"} with \code{sample.kind = "Rounding"}, and
#'   \code{doRNG::registerDoRNG()} is used so that the \code{\%dorng\%} loops are
#'   reproducible regardless of scheduling.
#' @param spec Optional \code{bigexp_spec} created by \code{bigexp_terms()}.
#'   If provided, the test reuses its locked expansion. The working formula
#'   becomes \code{bigexp_formula(spec, response_name)}, where
#'   \code{response_name} is taken from \code{response} if supplied, otherwise
#'   from the left-hand side of \code{formula}. Categorical sampling uses
#'   \code{spec$levels}, and numeric sampling prefers \code{spec$num_range}
#'   when available.
#' @param response Optional character name for the response variable to use when
#'   \code{spec} is supplied. If omitted, the response is taken from the
#'   left-hand side of \code{formula}.
#' @param use_spec_contrasts Logical; default \code{TRUE}. When \code{spec} is
#'   supplied and \code{use_spec_contrasts = TRUE}, the function replays
#'   \code{spec$settings$contrasts_options} on the parallel workers for
#'   deterministic factor coding.
#' @param ... Additional arguments passed to \code{SVEMnet()} and then to
#'   \code{glmnet()} (for example: \code{penalty.factor}, \code{offset},
#'   \code{lower.limits}, \code{upper.limits}, \code{standardize.response}, etc.).
#'   The \code{relaxed} setting is controlled by the \code{relaxed} argument of
#'   this function and any \code{relaxed} value passed via \code{...} is ignored
#'   with a warning.
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
#'   }
#'
#' @seealso \code{\link{SVEMnet}}, \code{\link{bigexp_terms}},
#'   \code{\link{bigexp_formula}}
#' @template ref-svem
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist SHASHo pSHASHo
#' @importFrom stats model.frame model.response model.matrix delete.response terms
#' @importFrom stats median complete.cases rgamma coef predict sd
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng% registerDoRNG
#' @importFrom parallel makeCluster stopCluster detectCores clusterCall
#'
#' @examples
#' \donttest{
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
#'     nCore          = 2,
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
#'     nCore              = 2,
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
    nCore = parallel::detectCores()-2,
    seed = NULL,
    spec = NULL,
    response = NULL,
    use_spec_contrasts = TRUE,
    ...
) {
  # --- basic choices ---
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)
  data <- as.data.frame(data)
  if (objective == "auto") (objective == "AIC")
  # Determine response and working formula (optionally from spec)
  if (!is.null(spec)) {
    resp_name <- if (!is.null(response)) {
      if (!is.character(response) || length(response) != 1L || !nzchar(response))
        stop("response must be a non-empty character scalar when provided.")
      response
    } else {
      as.character(formula[[2]])
    }
    if (!resp_name %in% names(data)) stop("Response '", resp_name, "' not found in 'data'.")
    f_use <- bigexp_formula(spec, resp_name)
  } else {
    f_use <- formula
    resp_name <- as.character(formula[[2]])
  }

  # Choose contrasts options (spec's, if requested)
  contrasts_opts <- if (!is.null(spec) && isTRUE(use_spec_contrasts)) {
    spec$settings$contrasts_options %||% getOption("contrasts")
  } else {
    getOption("contrasts")
  }

  # --- master RNG for serial prep (design, mixtures, factors) ---
  if (!is.null(seed)) {
    oldKinds <- RNGkind()
    on.exit({
      RNGkind(kind = oldKinds[1L], normal.kind = oldKinds[2L], sample.kind = oldKinds[3L])
    }, add = TRUE)

    # For R >= 3.6, sample.kind="Rounding" avoids deprecation warnings
    ok <- try(suppressWarnings(RNGkind("L'Ecuyer-CMRG", sample.kind = "Rounding")), silent = TRUE)
    if (inherits(ok, "try-error")) RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  # --- enforce single-threaded BLAS/OpenMP BEFORE spawning workers -------------
  Sys.setenv(
    OMP_NUM_THREADS        = "1",
    MKL_NUM_THREADS        = "1",
    OPENBLAS_NUM_THREADS   = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS    = "1"
  )

  # --- cluster setup -----------------------------------------------------------
  nCore <- max(1L, as.integer(`%||%`(nCore, parallel::detectCores())))
  cl <- parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  if (!is.null(seed)) {
    parallel::clusterCall(cl, function() {
      ok <- try(suppressWarnings(RNGkind("L'Ecuyer-CMRG", sample.kind = "Rounding")), silent = TRUE)
      if (inherits(ok, "try-error")) RNGkind("L'Ecuyer-CMRG")
      NULL
    })
  }

  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Enforce single-threaded BLAS/OpenMP on workers too (belt & suspenders)
  parallel::clusterCall(cl, function() {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    } else {
      Sys.setenv(
        OMP_NUM_THREADS        = "1",
        MKL_NUM_THREADS        = "1",
        OPENBLAS_NUM_THREADS   = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS    = "1"
      )
    }
    NULL
  })

  # Set contrasts options on workers
  parallel::clusterCall(cl, function(opts) { options(contrasts = opts); NULL }, contrasts_opts)

  # Register doRNG to make foreach iterations reproducible, independent of scheduling
  if (!is.null(seed)) doRNG::registerDoRNG(seed)

  # Sanitize ... so explicit 'relaxed' here cannot be overridden
  dots <- list(...)

  ## --- guard against binomial family: only Gaussian is supported here ----
  if ("family" %in% names(dots)) {
    fam <- dots$family
    fam_name <- tryCatch({
      if (inherits(fam, "family")) {
        fam$family                 # e.g. binomial() object
      } else if (is.function(fam)) {
        fam()$family               # e.g. family = binomial
      } else {
        as.character(fam)[1L]      # e.g. "binomial"
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

  # Training design pieces
  mf <- stats::model.frame(f_use, data)
  y  <- stats::model.response(mf)

  # For sampling, decide which raw columns are categorical vs continuous
  if (!is.null(spec)) {
    predictor_vars   <- spec$vars
    is_cat           <- spec$is_cat
    categorical_vars <- names(is_cat)[is_cat]
    continuous_vars  <- names(is_cat)[!is_cat]
  } else {
    predictor_vars  <- base::all.vars(stats::delete.response(stats::terms(f_use, data = data)))
    predictor_types <- sapply(data[predictor_vars], function(z) class(z)[1L])
    categorical_vars <- predictor_vars[predictor_types %in% c("factor", "character", "logical")]
    continuous_vars  <- setdiff(predictor_vars, categorical_vars)
  }

  # Mixture bookkeeping
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) mixture_vars <- c(mixture_vars, grp$vars)
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }
  nonmix_continuous_vars <- setdiff(continuous_vars, mixture_vars)

  # Non-mixture continuous via maximin LHS over ranges (prefer spec$num_range)
  if (length(nonmix_continuous_vars) > 0) {
    if (!is.null(spec) && !is.null(spec$num_range) && ncol(spec$num_range) > 0) {
      rng_mat <- spec$num_range[, colnames(spec$num_range) %in% nonmix_continuous_vars, drop = FALSE]
      missing <- setdiff(nonmix_continuous_vars, colnames(rng_mat))
      if (length(missing)) {
        add <- sapply(data[missing], function(col) range(col, na.rm = TRUE))
        rownames(add) <- c("min","max")
        rng_mat <- cbind(rng_mat, add)
      }
      rng_mat <- rng_mat[, nonmix_continuous_vars, drop = FALSE]
    } else {
      rng_mat <- sapply(data[nonmix_continuous_vars], function(col) range(col, na.rm = TRUE))
      rownames(rng_mat) <- c("min","max")
    }

    T_continuous_raw <- as.matrix(lhs::maximinLHS(nPoint, length(nonmix_continuous_vars)))
    T_continuous <- matrix(NA_real_, nrow = nPoint, ncol = length(nonmix_continuous_vars))
    colnames(T_continuous) <- nonmix_continuous_vars
    for (i in seq_along(nonmix_continuous_vars)) {
      lo <- rng_mat["min", i]; hi <- rng_mat["max", i]
      T_continuous[, i] <- T_continuous_raw[, i] * (hi - lo) + lo
    }
    T_continuous <- as.data.frame(T_continuous)
  } else {
    T_continuous <- NULL
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
      if (!is.null(spec)) {
        lv <- spec$levels[[v]]
        if (is.null(lv)) lv <- sort(unique(as.character(data[[v]])))
        T_categorical[[v]] <- factor(sample(lv, nPoint, replace = TRUE), levels = lv)
      } else {
        x <- data[[v]]
        if (is.factor(x)) {
          obs_lev <- levels(base::droplevels(x))
          T_categorical[[v]] <- factor(
            sample(obs_lev, nPoint, replace = TRUE),
            levels = levels(x)
          )
        } else {
          obs_lev <- sort(unique(as.character(x)))
          T_categorical[[v]] <- factor(
            sample(obs_lev, nPoint, replace = TRUE),
            levels = obs_lev
          )
        }
      }
    }
    T_categorical <- as.data.frame(T_categorical, stringsAsFactors = FALSE)
  }

  # Assemble evaluation grid
  parts <- list(T_continuous, T_mixture, T_categorical)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) stop("No predictors provided.")
  T_data <- do.call(cbind, parts)

  y_mean <- mean(y)

  # --- Originals: parallel SVEM fits (reproducible with %dorng%) ---------------
  if (isTRUE(verbose)) message("Fitting SVEM models to original data (parallel)...")
  M_Y <- foreach::foreach(
    i = 1:nSVEM,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats"),
    .options.snow = list(preschedule = FALSE)
  ) %dorng% {
    svem_model <- tryCatch({
      do.call(SVEMnet::SVEMnet, c(list(
        formula = f_use, data = data, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
        weight_scheme = weight_scheme, objective = objective,
        relaxed = relaxed
      ), dots))
    }, error = function(e) {
      message("Error in SVEMnet during SVEM fitting: ", e$message)
      return(NULL)
    })
    if (is.null(svem_model)) return(rep(NA_real_, nPoint))

    pred_res <- predict(svem_model, newdata = T_data, debias = FALSE, se.fit = TRUE)
    f_hat_Y_T <- pred_res$fit
    s_hat_Y_T <- pred_res$se.fit
    s_hat_Y_T[s_hat_Y_T == 0] <- 1e-6
    (f_hat_Y_T - y_mean) / s_hat_Y_T
  }

  # --- Permutations: parallel SVEM fits (reproducible with %dorng%) ------------
  if (isTRUE(verbose)) message("Starting permutation testing (parallel)...")
  start_time_perm <- Sys.time()
  M_pi_Y <- foreach::foreach(
    jloop = 1:nPerm,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats"),
    .options.snow = list(preschedule = FALSE)
  ) %dorng% {
    y_perm <- sample(y, replace = FALSE)
    data_perm <- data
    data_perm[[resp_name]] <- y_perm

    svem_model_perm <- tryCatch({
      do.call(SVEMnet::SVEMnet, c(list(
        formula = f_use, data = data_perm, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
        weight_scheme = weight_scheme, objective = objective,
        relaxed = relaxed
      ), dots))
    }, error = function(e) {
      message("Error in SVEMnet during permutation fitting: ", e$message)
      return(NULL)
    })
    if (is.null(svem_model_perm)) return(rep(NA_real_, nPoint))

    pred_res <- predict(svem_model_perm, newdata = T_data, debias = FALSE, se.fit = TRUE)
    f_hat_piY_T <- pred_res$fit
    s_hat_piY_T <- pred_res$se.fit
    s_hat_piY_T[s_hat_piY_T == 0] <- 1e-6
    h_piY <- (f_hat_piY_T - y_mean) / s_hat_piY_T

    if (isTRUE(verbose) && (jloop %% 10 == 0 || jloop == nPerm)) {
      elapsed_time <- Sys.time() - start_time_perm
      elapsed_secs <- as.numeric(elapsed_time, units = "secs")
      estimated_total_secs <- (elapsed_secs / jloop) * nPerm
      remaining_secs <- pmax(0, estimated_total_secs - elapsed_secs)
      remaining_time_formatted <- sprintf(
        "%02d:%02d:%02d",
        floor(remaining_secs / 3600),
        floor((remaining_secs %% 3600) / 60),
        floor(remaining_secs %% 60)
      )
      message(sprintf("Permutation %d/%d completed. Estimated time remaining: %s",
                      jloop, nPerm, remaining_time_formatted))
    }

    h_piY
  }

  # Gather and check
  M_Y    <- M_Y[stats::complete.cases(M_Y), , drop = FALSE]
  M_pi_Y <- M_pi_Y[stats::complete.cases(M_pi_Y), , drop = FALSE]
  if (nrow(M_Y) == 0)    stop("All SVEM fits on the original data failed.")
  if (nrow(M_pi_Y) == 0) stop("All SVEM fits on permuted data failed.")

  # Normalize by permutation mean/sd
  col_means_M_pi_Y <- colMeans(M_pi_Y, na.rm = TRUE)
  col_sds_M_pi_Y   <- apply(M_pi_Y, 2, sd, na.rm = TRUE)
  col_sds_M_pi_Y[col_sds_M_pi_Y == 0] <- 1e-6
  tilde_M_pi_Y <- scale(M_pi_Y, center = col_means_M_pi_Y, scale = col_sds_M_pi_Y)

  M_Y_centered <- sweep(M_Y, 2, col_means_M_pi_Y, "-")
  tilde_M_Y    <- sweep(M_Y_centered, 2, col_sds_M_pi_Y, "/")

  # SVD and distances (paper-consistent scaling)
  svd_res <- svd(tilde_M_pi_Y)
  s <- svd_res$d
  V <- svd_res$v

  n_perm <- nrow(tilde_M_pi_Y)
  p_cols <- ncol(tilde_M_pi_Y)

  # 1) correlation-PCA eigenvalues: (s^2)/(n_perm - 1)
  evalues_all <- (s^2) / max(n_perm - 1, 1L)

  # 2) Rescale to make eigenvalues sum to p (matches paper)
  ev_sum <- sum(evalues_all)
  if (!is.finite(ev_sum) || ev_sum <= 0) ev_sum <- 1e-12
  evalues_all <- evalues_all / ev_sum * p_cols

  # Explained-variance percent under this convention
  cumsum_evalues <- cumsum(evalues_all) / p_cols * 100
  percent <- max(min(percent, 99.999), 0.1)
  k_idx <- which(cumsum_evalues >= percent)[1L]
  if (!length(k_idx) || is.na(k_idx)) k_idx <- length(evalues_all)
  k_idx <- min(k_idx, ncol(V))
  k_idx <- max(k_idx, 1L)

  evalues  <- evalues_all[seq_len(k_idx)]
  evectors <- V[, seq_len(k_idx), drop = FALSE]

  # Force matrices to avoid accidental data.frames
  tilde_M_pi_Y <- as.matrix(tilde_M_pi_Y)
  tilde_M_Y    <- as.matrix(tilde_M_Y)

  # Project
  Z_pi <- tilde_M_pi_Y %*% evectors   # nPerm x k
  Z_Y  <- tilde_M_Y    %*% evectors   # nSVEM x k

  # Numerical guard: zero/NA eigenvalues
  evalues[!is.finite(evalues) | evalues <= 0] <- 1e-12

  T2_perm <- rowSums((Z_pi^2) / rep(evalues, each = nrow(Z_pi)))
  d_pi_Y  <- sqrt(T2_perm)

  T2_Y <- rowSums((Z_Y^2) / rep(evalues, each = nrow(Z_Y)))
  d_Y  <- sqrt(T2_Y)

  if (length(d_pi_Y) == 0) stop("No valid permutation distances to fit a distribution.")

  # SHASHo fit
  suppressMessages(
    distribution_fit <- tryCatch({
      gamlss::gamlss(
        d_pi_Y ~ 1,
        family = gamlss.dist::SHASHo(mu.link = "identity", sigma.link = "log",
                                     nu.link = "identity", tau.link = "log"),
        control = gamlss::gamlss.control(n.cyc = 1000, trace = FALSE)
      )
    }, error = function(e) {
      message("Error in fitting SHASHo distribution: ", e$message)
      NULL
    })
  )
  if (is.null(distribution_fit)) stop("Failed to fit SHASHo distribution.")

  mu    <- as.numeric(stats::coef(distribution_fit, what = "mu"))
  sigma <- exp(as.numeric(stats::coef(distribution_fit, what = "sigma")))
  nu    <- as.numeric(stats::coef(distribution_fit, what = "nu"))
  tau   <- exp(as.numeric(stats::coef(distribution_fit, what = "tau")))

  p_values <- 1 - gamlss.dist::pSHASHo(d_Y, mu = mu, sigma = sigma, nu = nu, tau = tau)
  p_value  <- stats::median(p_values)

  data_d <- data.frame(
    D = c(d_Y, d_pi_Y),
    Source_Type = c(rep("Original", length(d_Y)), rep("Permutation", length(d_pi_Y))),
    Response = resp_name
  )

  results_list <- list(
    p_value = p_value,
    p_values = p_values,
    d_Y = d_Y,
    d_pi_Y = d_pi_Y,
    distribution_fit = distribution_fit,
    data_d = data_d
  )
  class(results_list) <- "svem_significance_test"
  results_list
}

# internal helper
`%||%` <- function(a, b) if (!is.null(a)) a else b
