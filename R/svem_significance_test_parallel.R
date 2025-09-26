#' SVEM Significance Test with Mixture Support (Parallel Version)
#'
#' Whole-model significance test using SVEM with support for mixture factor groups,
#' parallelizing the SVEM fits for originals and permutations.
#'
#' Identical in logic to \code{svem_significance_test()} but runs the expensive
#' SVEM refits in parallel using \code{foreach} + \code{doParallel}. Random draws
#' (including permutations) use \code{RNGkind("L'Ecuyer-CMRG")} for parallel-suitable streams.
#'
#' This parallel version can optionally reuse a deterministic, locked expansion built
#' with \code{bigexp_terms()}. Provide \code{spec} (and optionally \code{response}) to
#' ensure that categorical levels, contrasts, and the polynomial/interaction structure are
#' identical across repeated calls and across multiple responses of the same factor space.
#'
#' @param formula A formula specifying the model to be tested. If \code{spec} is provided,
#'   the right-hand side is ignored and replaced by the locked expansion in \code{spec}.
#' @param data A data frame containing the variables in the model.
#' @param mixture_groups Optional list describing one or more mixture factor
#'   groups. Each element of the list should be a list with components
#'   \code{vars} (character vector of column names), \code{lower} (numeric vector of
#'   lower bounds of the same length as \code{vars}), \code{upper} (numeric vector
#'   of upper bounds of the same length), and \code{total} (scalar specifying the
#'   sum of the mixture variables). All mixture variables must be included in \code{vars},
#'   and no variable can appear in more than one mixture group. Defaults to \code{NULL}.
#' @param nPoint Number of random points in the factor space (default: 2000).
#' @param nSVEM Number of SVEM fits on the original data (default: 10).
#' @param nPerm Number of SVEM fits on permuted responses for the reference
#'   distribution (default: 150).
#' @param percent Percentage of variance to capture in the SVD (default: 90).
#' @param nBoot Number of bootstrap iterations within each SVEM fit (default: 100).
#' @param glmnet_alpha The alpha parameter(s) for glmnet (default: \code{c(1)}).
#' @param weight_scheme Weighting scheme for SVEM (default: "SVEM").
#' @param objective Objective used inside \code{SVEMnet()} to pick the bootstrap
#'   path solution. One of \code{"auto"}, \code{"wAIC"}, \code{"wBIC"}, \code{"wSSE"}
#'   (default: \code{"auto"}).  (Note: \code{"wGIC"} is no longer supported.)
#' @param auto_ratio_cutoff Single cutoff for the automatic rule when
#'   \code{objective = "auto"} (default 1.3). With \code{r = n_X/p_X}, if
#'   \code{r >= auto_ratio_cutoff} use wAIC; else wBIC. Passed to \code{SVEMnet()}.
#' @param relaxed Logical; default \code{FALSE}. When \code{TRUE}, inner
#'   \code{SVEMnet()} fits use glmnet's relaxed elastic net path and select both
#'   lambda and relaxed gamma on each bootstrap. When \code{FALSE}, the standard
#'   glmnet path is used. This value is passed through to \code{SVEMnet()}.
#'   Note: if \code{relaxed = TRUE} and \code{glmnet_alpha} includes \code{0}, ridge
#'   (\code{alpha = 0}) is dropped by \code{SVEMnet()} for relaxed fits.
#' @param verbose Logical; if \code{TRUE}, displays progress messages (default: \code{TRUE}).
#' @param nCore Number of CPU cores for parallel processing (default: all available cores).
#' @param seed Optional integer seed for reproducible parallel RNG (default: \code{NULL}).
#' @param spec Optional \code{bigexp_spec} created by \code{bigexp_terms()}. If provided,
#'   the test reuses its locked expansion. The working formula becomes
#'   \code{bigexp_formula(spec, response_name)}, where \code{response_name} is taken from
#'   \code{response} if supplied, otherwise from the left-hand side of \code{formula}.
#'   Categorical sampling uses \code{spec$levels} and numeric sampling prefers \code{spec$num_range}
#'   when available.
#' @param response Optional character name for the response variable to use when \code{spec}
#'   is supplied. If omitted, the response is taken from the left-hand side of \code{formula}.
#' @param use_spec_contrasts Logical; default \code{TRUE}. When \code{spec} is supplied and
#'   \code{use_spec_contrasts=TRUE}, the function replays \code{spec$settings$contrasts_options}
#'   on the parallel workers for deterministic coding.
#' @param ... Additional arguments passed to \code{SVEMnet()} and then to \code{glmnet()}
#'   (for example: \code{penalty.factor}, \code{offset}, \code{lower.limits},
#'   \code{upper.limits}, \code{standardize.response}, etc.). The \code{relaxed}
#'   setting is controlled by the \code{relaxed} argument of this function and
#'   any \code{relaxed} value passed via \code{...} is ignored with a warning.
#'
#' @return A list of class \code{svem_significance_test} containing the test results.
#'
#' @seealso \code{svem_significance_test}, \code{bigexp_terms}, \code{bigexp_formula}
#'
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist SHASHo pSHASHo
#' @importFrom stats model.frame model.response model.matrix delete.response terms
#' @importFrom stats median complete.cases rgamma coef predict
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterSetRNGStream detectCores clusterCall
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
#'   # Parallel significance test (default relaxed = FALSE)
#'   res <- svem_significance_test_parallel(
#'     y ~ A + B + C + X + F,
#'     data           = dat,
#'     mixture_groups = mix_spec,
#'     glmnet_alpha   = c(1),
#'     weight_scheme  = "SVEM",
#'     objective      = "auto",
#'     auto_ratio_cutoff = 1.3,
#'     relaxed        = FALSE,   # default, shown for clarity
#'     nCore          = 2,
#'     seed           = 123,
#'     verbose        = FALSE
#'   )
#'   print(res$p_value)
#' }
#' @export
svem_significance_test_parallel <- function(
    formula, data, mixture_groups = NULL,
    nPoint = 2000, nSVEM = 10, nPerm = 150,
    percent = 90, nBoot = 100,
    glmnet_alpha = c(1),
    weight_scheme = c("SVEM"),
    objective = c("auto", "wAIC", "wBIC", "wSSE"),
    auto_ratio_cutoff = 1.3,
    relaxed = FALSE,
    verbose = TRUE,
    nCore = parallel::detectCores(),
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

  # --- parallel RNG & cluster setup ---
  RNGkind("L'Ecuyer-CMRG")
  if (!is.null(seed)) set.seed(seed)

  nCore <- max(1L, as.integer(`%||%`(nCore, parallel::detectCores())))
  cl <- parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)

  # Set contrasts options on workers without referencing a missing symbol
  parallel::clusterCall(cl, function(opts) { options(contrasts = opts); NULL }, contrasts_opts)

  # Sanitize ... so explicit 'relaxed' here cannot be overridden
  dots <- list(...)
  if ("relaxed" %in% names(dots)) {
    warning("Ignoring 'relaxed' in '...'; use the 'relaxed' argument of svem_significance_test_parallel().")
    dots$relaxed <- NULL
  }

  # Training design pieces
  mf <- stats::model.frame(f_use, data)
  y  <- stats::model.response(mf)

  # For sampling, decide which raw columns are categorical vs continuous
  if (!is.null(spec)) {
    predictor_vars  <- spec$vars
    is_cat          <- spec$is_cat
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
    if (!is.null(spec) && ncol(spec$num_range) > 0) {
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

  # --- Originals: parallel SVEM fits ---
  if (isTRUE(verbose)) message("Fitting SVEM models to original data (parallel)...")
  M_Y <- foreach::foreach(
    i = 1:nSVEM,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats")
  ) %dopar% {
    svem_model <- tryCatch({
      do.call(SVEMnet::SVEMnet, c(list(
        formula = f_use, data = data, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
        weight_scheme = weight_scheme, objective = objective,
        auto_ratio_cutoff = auto_ratio_cutoff,
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

  # --- Permutations: parallel SVEM fits ---
  if (isTRUE(verbose)) message("Starting permutation testing (parallel)...")
  start_time_perm <- Sys.time()
  M_pi_Y <- foreach::foreach(
    jloop = 1:nPerm,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats")
  ) %dopar% {
    y_perm <- sample(y, replace = FALSE)
    data_perm <- data
    data_perm[[resp_name]] <- y_perm

    svem_model_perm <- tryCatch({
      do.call(SVEMnet::SVEMnet, c(list(
        formula = f_use, data = data_perm, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
        weight_scheme = weight_scheme, objective = objective,
        auto_ratio_cutoff = auto_ratio_cutoff,
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

  # SVD and distances (hardened, same logic)
  svd_res <- svd(tilde_M_pi_Y)
  s <- svd_res$d
  V <- svd_res$v

  evalues_temp <- s^2
  # explained-variance percent (same as before; the ncol() factor cancels)
  cumsum_evalues <- cumsum(evalues_temp) / sum(evalues_temp) * 100

  k_idx <- which(cumsum_evalues >= percent)[1L]
  if (!length(k_idx) || is.na(k_idx)) k_idx <- length(evalues_temp)
  # guard against asking for more PCs than exist
  k_idx <- min(k_idx, ncol(V))

  evalues  <- evalues_temp[seq_len(k_idx)]
  evectors <- V[, seq_len(k_idx), drop = FALSE]

  # force matrices to avoid accidental data.frames
  tilde_M_pi_Y <- as.matrix(tilde_M_pi_Y)
  tilde_M_Y    <- as.matrix(tilde_M_Y)

  # Project, then divide component-wise by eigenvalues (avoid diag() shape traps)
  Z_pi <- tilde_M_pi_Y %*% evectors   # nPerm x k
  Z_Y  <- tilde_M_Y    %*% evectors   # nSVEM x k

  # numerical guard: zero/NA eigenvalues
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
