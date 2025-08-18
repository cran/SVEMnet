#' SVEM Significance Test with Mixture Support (Parallel Version)
#'
#' Whole-model significance test using SVEM with support for mixture factor groups,
#' parallelizing the SVEM fits for originals and permutations.
#'
#' @inheritParams svem_significance_test
#' @param nCore Number of CPU cores for parallel processing (default: all available cores).
#' @param seed Optional integer seed for reproducible parallel RNG (default: NULL).
#' @return A list of class `svem_significance_test` containing the test results.
#' @details
#' Identical to \code{svem_significance_test()} but runs the expensive
#' SVEM refits in parallel using \code{foreach} + \code{doParallel}. Random draws
#' (including permutations) use \code{RNGkind("L'Ecuyer-CMRG")} for parallel-suitable streams.
#'
#' @seealso svem_significance_test svem_significance_test_parallel
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist SHASHo pSHASHo
#' @importFrom stats model.response delete.response terms reformulate median complete.cases rgamma
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterSetRNGStream detectCores
#' @export
svem_significance_test_parallel <- function(
    formula, data, mixture_groups = NULL,
    nPoint = 2000, nSVEM = 5, nPerm = 125,
    percent = 85, nBoot = 200,
    glmnet_alpha = c(1),
    weight_scheme = c("SVEM"),
    objective = c("wAIC", "wSSE"),
    verbose = TRUE,
    nCore = parallel::detectCores(),
    seed = NULL, ...
) {
  objective <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)

  RNGkind("L'Ecuyer-CMRG")
  if (!is.null(seed)) set.seed(seed)

  cl <- parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)

  data <- as.data.frame(data)

  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) X <- X[, -intercept_col, drop = FALSE]

  predictor_vars  <- all.vars(delete.response(terms(formula, data = data)))
  predictor_types <- sapply(data[predictor_vars], class)
  continuous_vars <- predictor_vars[!predictor_types %in% c("factor", "character")]
  categorical_vars <- predictor_vars[predictor_types %in% c("factor", "character")]

  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) mixture_vars <- c(mixture_vars, grp$vars)
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }
  nonmix_continuous_vars <- setdiff(continuous_vars, mixture_vars)

  if (length(nonmix_continuous_vars) > 0) {
    ranges <- sapply(data[nonmix_continuous_vars], function(col) range(col, na.rm = TRUE))
    T_continuous_raw <- as.matrix(lhs::maximinLHS(nPoint, length(nonmix_continuous_vars)))
    T_continuous <- matrix(NA_real_, nrow = nPoint, ncol = length(nonmix_continuous_vars))
    colnames(T_continuous) <- nonmix_continuous_vars
    for (i in seq_along(nonmix_continuous_vars)) {
      T_continuous[, i] <- T_continuous_raw[, i] * (ranges[2, i] - ranges[1, i]) + ranges[1, i]
    }
    T_continuous <- as.data.frame(T_continuous)
  } else {
    T_continuous <- NULL
  }

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

  # Categorical sampling (use observed levels; keep training levels attribute for factors)
  T_categorical <- NULL
  if (length(categorical_vars) > 0) {
    T_categorical <- vector("list", length(categorical_vars))
    names(T_categorical) <- categorical_vars
    for (v in categorical_vars) {
      x <- data[[v]]
      if (is.factor(x)) {
        obs_lev <- levels(base::droplevels(x))  # <-- FIX: use base::droplevels
        T_categorical[[v]] <- factor(
          sample(obs_lev, nPoint, replace = TRUE),
          levels = levels(x)                     # keep original full level set
        )
      } else {
        obs_lev <- sort(unique(as.character(x)))
        T_categorical[[v]] <- factor(
          sample(obs_lev, nPoint, replace = TRUE),
          levels = obs_lev
        )
      }
    }
    T_categorical <- as.data.frame(T_categorical, stringsAsFactors = FALSE)
  }



  parts <- list(T_continuous, T_mixture, T_categorical)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) stop("No predictors provided.")
  T_data <- do.call(cbind, parts)

  y_mean <- mean(y)
  if (verbose) message("Fitting SVEM models to original data with mixture handling (parallel)...")

  M_Y <- foreach::foreach(
    i = 1:nSVEM,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats"),
    .export = c("SVEMnet")
  ) %dopar% {
    svem_model <- tryCatch({
      SVEMnet::SVEMnet(formula, data = data, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
                       weight_scheme = weight_scheme, objective = objective, ...)
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

  if (verbose) message("Starting permutation testing (parallel)...")
  start_time_perm <- Sys.time()

  M_pi_Y <- foreach::foreach(
    jloop = 1:nPerm,
    .combine = rbind,
    .packages = c("SVEMnet", "glmnet", "stats"),
    .export = c("SVEMnet")
  ) %dopar% {
    y_perm <- sample(y, replace = FALSE)
    data_perm <- data
    data_perm[[as.character(formula[[2]])]] <- y_perm

    svem_model_perm <- tryCatch({
      SVEMnet::SVEMnet(formula, data = data_perm, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
                       weight_scheme = weight_scheme, objective = objective, ...)
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

    if (verbose && (jloop %% 10 == 0 || jloop == nPerm)) {
      elapsed_time <- Sys.time() - start_time_perm
      elapsed_secs <- as.numeric(elapsed_time, units = "secs")
      estimated_total_secs <- (elapsed_secs / jloop) * nPerm
      remaining_secs <- estimated_total_secs - elapsed_secs
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

  M_Y    <- M_Y[complete.cases(M_Y), , drop = FALSE]
  M_pi_Y <- M_pi_Y[complete.cases(M_pi_Y), , drop = FALSE]
  if (nrow(M_Y) == 0)    stop("All SVEM fits on the original data failed.")
  if (nrow(M_pi_Y) == 0) stop("All SVEM fits on permuted data failed.")

  col_means_M_pi_Y <- colMeans(M_pi_Y, na.rm = TRUE)
  col_sds_M_pi_Y   <- apply(M_pi_Y, 2, sd, na.rm = TRUE)
  col_sds_M_pi_Y[col_sds_M_pi_Y == 0] <- 1e-6
  tilde_M_pi_Y <- scale(M_pi_Y, center = col_means_M_pi_Y, scale = col_sds_M_pi_Y)

  M_Y_centered <- sweep(M_Y, 2, col_means_M_pi_Y, "-")
  tilde_M_Y    <- sweep(M_Y_centered, 2, col_sds_M_pi_Y, "/")

  svd_res <- svd(tilde_M_pi_Y)
  s <- svd_res$d
  V <- svd_res$v

  evalues_temp <- s^2
  evalues_temp <- evalues_temp / sum(evalues_temp) * ncol(tilde_M_pi_Y)
  cumsum_evalues <- cumsum(evalues_temp) / sum(evalues_temp) * 100
  k_idx <- which(cumsum_evalues >= percent)[1]
  if (is.na(k_idx)) k_idx <- length(evalues_temp)
  evalues  <- evalues_temp[1:k_idx]
  evectors <- V[, 1:k_idx]

  T2_perm <- rowSums((tilde_M_pi_Y %*% evectors %*% diag(1 / evalues)) * (tilde_M_pi_Y %*% evectors))
  d_pi_Y  <- sqrt(T2_perm)

  T2_Y <- rowSums((tilde_M_Y %*% evectors %*% diag(1 / evalues)) * (tilde_M_Y %*% evectors))
  d_Y  <- sqrt(T2_Y)

  if (length(d_pi_Y) == 0) stop("No valid permutation distances to fit a distribution.")

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

  mu    <- as.numeric(coef(distribution_fit, what = "mu"))
  sigma <- exp(as.numeric(coef(distribution_fit, what = "sigma")))
  nu    <- as.numeric(coef(distribution_fit, what = "nu"))
  tau   <- exp(as.numeric(coef(distribution_fit, what = "tau")))

  p_values <- 1 - gamlss.dist::pSHASHo(d_Y, mu = mu, sigma = sigma, nu = nu, tau = tau)
  p_value  <- median(p_values)

  response_name <- as.character(formula[[2]])
  data_d <- data.frame(
    D = c(d_Y, d_pi_Y),
    Source_Type = c(rep("Original", length(d_Y)), rep("Permutation", length(d_pi_Y))),
    Response = response_name
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
  return(results_list)
}

