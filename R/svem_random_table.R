#
# NOTE: This utility function
# allows users to explore the fitted SVEMnet model by generating a
# random grid of predictor values and computing the corresponding model
# predictions.  Mixture factor groups can be specified so that
# composition variables are sampled from a Dirichlet distribution.
#' Generate a Random Prediction Table for a Fitted SVEMnet Model
#'
#' This utility function generates a random sample of points from the
#' predictor space and computes the corresponding predicted responses from
#' a fitted SVEMnet model.  It can be used to explore the fitted
#' response surface in a way analogous to JMP's "Output Random Table"
#' feature.  The function recognizes mixture factor groups and draws
#' Dirichlet-distributed compositions within the specified bounds so
#' that mixture variables sum to a user-supplied total.  Continuous
#' non-mixture variables are sampled uniformly across their observed
#' ranges using a maximin Latin hypercube design, and categorical
#' variables are sampled from their observed levels.  No random noise is
#' added to the predicted responses.
#'
#' @param formula A formula specifying the fitted model.  This should be the
#'   same formula used when fitting the SVEMnet model.
#' @param data A data frame containing the variables in the model.
#' @param n Number of random points to generate (default: 1000).
#' @param mixture_groups Optional list describing mixture factor groups.  Each
#'   element should be a list with components `vars` (character vector of
#'   mixture variable names), `lower` (numeric vector of lower bounds),
#'   `upper` (numeric vector of upper bounds) and `total` (scalar sum).
#'   See `svem_significance_test_with_mixture()` for details.  Defaults
#'   to `NULL` (no mixtures).
#' @param nBoot Number of bootstrap iterations to use when fitting the
#'   SVEMnet model (default: 200).
#' @param glmnet_alpha Elastic net mixing parameter(s) passed to `SVEMnet`
#'   (default: `c(1)`).
#' @param weight_scheme Weighting scheme for SVEM (default: "SVEM").
#' @param objective Objective function for SVEM ("wAIC" or "wSSE";
#'   default: "wAIC").
#' @param debias Logical; if `TRUE`, the debiasing coefficients of the
#'   fitted model are applied when predicting (default: `FALSE`).
#' @param ... Additional arguments passed to `SVEMnet()` and then to
#'   `glmnet()`.
#' @return A data frame containing the sampled predictor values and the
#'   corresponding predicted responses.  The response column is named
#'   according to the left-hand side of `formula`.
#' @details
#' This function first fits an SVEMnet model using the supplied
#' parameters.  It then generates a random grid of points in the
#' predictor space, honouring mixture constraints if `mixture_groups`
#' is provided.  Predictions are computed from the fitted model on these
#' points.  No random noise is added; the predictions come directly from
#' the model.  If you wish to explore the uncertainty of predictions,
#' consider adding noise separately or using the standard error output
#' from `predict.svem_model()`.
#'
#' @seealso `SVEMnet`, `predict.svem_model`, `svem_significance_test_with_mixture`.
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 40
#'
#' # Helper to generate training data mixtures with bounds
#' sample_trunc_dirichlet <- function(n, lower, upper, total) {
#'   k <- length(lower)
#'   min_sum <- sum(lower); max_sum <- sum(upper)
#'   stopifnot(total >= min_sum, total <= max_sum)
#'   avail <- total - min_sum
#'   out <- matrix(NA_real_, n, k)
#'   i <- 1
#'   while (i <= n) {
#'     g <- rgamma(k, 1, 1)
#'     w <- g / sum(g)
#'     x <- lower + avail * w
#'     if (all(x <= upper + 1e-12)) {
#'       out[i, ] <- x; i <- i + 1
#'     }
#'   }
#'   out
#' }
#'
#' # Three mixture factors (A, B, C) with distinct bounds; sum to total = 1
#' lower <- c(0.10, 0.20, 0.05)
#' upper <- c(0.60, 0.70, 0.50)
#' total <- 1.0
#' ABC   <- sample_trunc_dirichlet(n, lower, upper, total)
#' A <- ABC[, 1]; B <- ABC[, 2]; C <- ABC[, 3]
#'
#' # Additional predictors
#' X <- runif(n)
#' F <- factor(sample(letters[1:3], n, replace = TRUE))
#'
#' # Response
#' y <- 1 + 2*A + 3*B + 1.5*C + 0.5*X +
#'      ifelse(F == "a", 0, ifelse(F == "b", 1, -1)) +
#'      rnorm(n, sd = 0.3)
#'
#' dat <- data.frame(y = y, A = A, B = B, C = C, X = X, F = F)
#'
#' # Mixture specification for the random table generator
#' mix_spec <- list(
#'   list(
#'     vars  = c("A", "B", "C"),
#'     lower = c(0.10, 0.20, 0.05),
#'     upper = c(0.60, 0.70, 0.50),
#'     total = 1.0
#'   )
#' )
#'
#' # Fit SVEMnet and generate 50 random points
#' rand_tab <- svem_random_table(
#'   y ~ A + B + C + X + F,
#'   data           = dat,
#'   n              = 50,
#'   mixture_groups = mix_spec,
#'   nBoot          = 50,
#'   glmnet_alpha   = c(1),
#'   weight_scheme  = "SVEM",
#'   objective      = "wAIC",
#'   debias         = FALSE
#' )
#'
#' # Check mixture validity in the generated table
#' stopifnot(all(abs((rand_tab$A + rand_tab$B + rand_tab$C) - 1) < 1e-8))
#' summary(rand_tab[c("A","B","C")])
#' head(rand_tab)
#' }
#' @export
svem_random_table <- function(formula, data, n = 1000, mixture_groups = NULL,
                              nBoot = 200, glmnet_alpha = c(1),
                              weight_scheme = c("SVEM"),
                              objective = c("wAIC", "wSSE"),
                              debias = FALSE, ...) {
  objective <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)
  data <- as.data.frame(data)

  # Fit the model
  svem_model <- SVEMnet(formula, data = data, nBoot = nBoot,
                        glmnet_alpha = glmnet_alpha,
                        weight_scheme = weight_scheme,
                        objective = objective, ...)

  # Parse model terms
  mf <- model.frame(formula, data)
  y_var <- as.character(formula[[2]])
  predictor_vars <- all.vars(delete.response(terms(formula, data = data)))

  # Identify types
  is_cat <- vapply(data[predictor_vars], function(x) is.factor(x) || is.character(x), logical(1))
  is_num <- vapply(data[predictor_vars], is.numeric, logical(1))
  categorical_vars <- predictor_vars[is_cat]
  continuous_vars  <- predictor_vars[is_num]

  # Collect mixture vars (and guard against overlaps)
  mixture_vars <- character(0)
  if (!is.null(mixture_groups)) {
    for (grp in mixture_groups) mixture_vars <- c(mixture_vars, grp$vars)
    if (any(duplicated(mixture_vars))) {
      dups <- unique(mixture_vars[duplicated(mixture_vars)])
      stop("Mixture variables appear in multiple groups: ", paste(dups, collapse = ", "))
    }
  }
  nonmix_continuous <- setdiff(continuous_vars, mixture_vars)

  # Non-mixture continuous via maximin LHS
  if (length(nonmix_continuous) > 0) {
    ranges <- vapply(data[nonmix_continuous],
                     function(col) {
                       r <- range(col, na.rm = TRUE)
                       if (!is.finite(r[1]) || !is.finite(r[2])) r <- c(0, 1)
                       r
                     },
                     numeric(2))
    T_continuous_raw <- as.matrix(lhs::maximinLHS(n, length(nonmix_continuous)))
    T_continuous <- matrix(NA_real_, nrow = n, ncol = length(nonmix_continuous))
    colnames(T_continuous) <- nonmix_continuous
    for (i in seq_along(nonmix_continuous)) {
      T_continuous[, i] <- T_continuous_raw[, i] * (ranges[2, i] - ranges[1, i]) + ranges[1, i]
    }
    T_continuous <- as.data.frame(T_continuous)
  } else {
    T_continuous <- NULL
  }

  # Mixture groups: truncated-simplex Dirichlet with rejection on upper
  .sample_trunc_dirichlet <- function(n, lower, upper, total,
                                      alpha = NULL, oversample = 4, max_tries = 10000L) {
    k <- length(lower)
    if (length(upper) != k) stop("upper must have length equal to lower.")
    if (is.null(alpha)) alpha <- rep(1, k)

    min_sum <- sum(lower)
    max_sum <- sum(upper)
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
        keep_idx <- which(ok)
        take <- min(length(keep_idx), n - filled)
        res[(filled + 1):(filled + take), ] <- cand[keep_idx[seq_len(take)], , drop = FALSE]
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
    all_mix_vars <- unlist(lapply(mixture_groups, `[[`, "vars"))
    T_mixture <- matrix(NA_real_, nrow = n, ncol = length(all_mix_vars))
    colnames(T_mixture) <- all_mix_vars

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
      vals <- .sample_trunc_dirichlet(n, lower, upper, total)
      colnames(vals) <- vars
      T_mixture[, vars] <- vals
    }
    T_mixture <- as.data.frame(T_mixture)
  }

  # Categorical sampling (force factors even if source was character)
  T_categorical <- NULL
  if (length(categorical_vars) > 0) {
    T_categorical <- vector("list", length(categorical_vars))
    names(T_categorical) <- categorical_vars
    for (v in categorical_vars) {
      x <- data[[v]]
      if (is.factor(x)) {
        lev <- levels(x)
        T_categorical[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
      } else {
        lev <- sort(unique(as.character(x)))
        T_categorical[[v]] <- factor(sample(lev, n, replace = TRUE), levels = lev)
      }
    }
    T_categorical <- as.data.frame(T_categorical, stringsAsFactors = FALSE)
  }

  # Combine
  parts <- list(T_continuous, T_mixture, T_categorical)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) stop("No predictors provided.")
  T_data <- do.call(cbind, parts)

  # Align column order to model predictors if possible
  common <- intersect(predictor_vars, colnames(T_data))
  T_data <- T_data[, common, drop = FALSE]

  # Predict
  preds <- predict(svem_model, newdata = T_data, debias = debias)
  if (is.list(preds) && !is.null(preds$fit)) preds <- preds$fit

  out <- T_data
  out[[y_var]] <- as.numeric(preds)
  rownames(out) <- NULL
  out
}
