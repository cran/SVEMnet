#' SVEM Significance Test
#'
#' Performs a whole-model significance test using the SVEM framework, handling both continuous and categorical predictors.
#'
#' @param formula A formula specifying the model to be tested.
#' @param data A data frame containing the variables in the model.
#' @param nPoint The number of random points to generate in the factor space (default: 2000).
#' @param nSVEM The number of SVEM models to fit to the original data (default: 5).
#' @param nPerm The number of SVEM models to fit to permuted data for reference distribution (default: 125).
#' @param percent The percentage of variance to capture in the SVD (default: 90).
#' @param nBoot The number of bootstrap iterations within SVEM (default: 200).
#' @param glmnet_alpha The alpha parameter(s) for glmnet (default: \code{c(1)}).
#' @param weight_scheme The weight scheme to use in SVEM (default: "SVEM"). Valid options are "SVEM" and "FWR".
#' @param debias Logical; debiasing option passed to \code{\link{SVEMnet}} (default: \code{FALSE}).
#' @param objective Character; the objective function to use in \code{\link{SVEMnet}}. Options are "wAIC" or "wSSE" (default: "wSSE").
#' @param verbose Logical; if \code{TRUE}, displays progress messages (default: \code{FALSE}).
#' @param ... Additional arguments passed to the underlying \code{SVEMnet()} and then \code{glmnet()} functions.
#' @return A list containing the test results.
#' @details
#' The `svem_significance_test` function implements a whole-model test designed to gauge the significance of a fitted SVEM model compared to the null hypothesis of a constant response surface. This method helps identify responses that have relatively stronger or weaker relationships with study factors.
#'
#' The test constructs standardized predictions by centering the SVEM predictions by the response mean and scaling by the ensemble standard deviation. A reference distribution is created by fitting the SVEM model to multiple randomized permutations of the response vector. The Mahalanobis distances of the original and permuted models are calculated using a reduced-rank singular value decomposition.
#'
#' The R code to perform this test (using matrices of nSVEM and nPerm predictions) is taken from the supplementary material of Karl (2024).
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @references
#' Karl, A. T. (2024). A randomized permutation whole-model test heuristic for Self-Validated Ensemble Models (SVEM). \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{249}, 105122. \doi{10.1016/j.chemolab.2024.105122}
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(0)
#' n <- 21
#' X1 <- runif(n)
#' X2 <- runif(n)
#' X3 <- runif(n)
#' y <- 1 + X1 +  X2 + X1 * X2 + X1^2 + rnorm(n)
#' data <- data.frame(y, X1, X2, X3)
#'
#' # Perform the SVEM significance test
#' test_result <- svem_significance_test(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2) + I(X3^2),
#'   data = data,
#'   nPoint = 2000,
#'   nSVEM = 5,
#'   nPerm = 125,
#'   nBoot = 200
#'
#' )
#'
#' # View the p-value
#' test_result$p_value
#'
#' test_result2 <- svem_significance_test(
#'   y ~ (X1 + X2 )^2 + I(X1^2) + I(X2^2),
#'   data = data,
#'   nPoint = 2000,
#'   nSVEM = 5,
#'   nPerm = 125,
#'   nBoot = 200
#' )
#'
#' # View the p-value
#' test_result2$p_value
#'
#' # Plot the Mahalanobis distances
#' plot(test_result,test_result2)
#' }
#'
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist GG SHASHo pGG pSHASHo
#' @importFrom stats model.response delete.response terms reformulate median

#' @export
svem_significance_test <- function(formula, data, nPoint = 2000, nSVEM = 5, nPerm = 125,
                                   percent = 90, nBoot = 200, glmnet_alpha = c(1),
                                   weight_scheme = c("SVEM", "FWR"), debias = FALSE,
                                   objective = c("wSSE", "wAIC"),verbose = FALSE,...) {
  # Match the arguments
  objective <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)

  # Ensure data is a data frame
  data <- as.data.frame(data)

  # Create model frame and extract response and predictors
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)

  # Remove intercept column if present
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X <- X[, -intercept_col, drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  # Get predictor variable names
  predictor_vars <- all.vars(delete.response(terms(formula, data = data)))

  # Identify categorical (factor) and continuous predictors
  predictor_types <- sapply(data[predictor_vars], class)
  continuous_vars <- predictor_vars[!predictor_types %in% c("factor", "character")]
  categorical_vars <- predictor_vars[predictor_types %in% c("factor", "character")]

  # Get ranges of continuous predictors from the data
  if (length(continuous_vars) > 0) {
    ranges <- sapply(data[continuous_vars], function(col) range(col, na.rm = TRUE))
  }

  # Step 1: Generate T_factors using LHS over the ranges of continuous predictors
  if (length(continuous_vars) > 0) {
    T_continuous_raw <- as.matrix(maximinLHS(nPoint, length(continuous_vars)))
    T_continuous <- matrix(NA, nrow = nPoint, ncol = length(continuous_vars))
    colnames(T_continuous) <- continuous_vars
    for (i in seq_along(continuous_vars)) {
      T_continuous[, i] <- T_continuous_raw[, i] * (ranges[2, i] - ranges[1, i]) + ranges[1, i]
    }
  } else {
    T_continuous <- NULL
  }

  # For categorical variables, sample levels randomly
  if (length(categorical_vars) > 0) {
    T_categorical <- as.data.frame(matrix(NA, nrow = nPoint, ncol = length(categorical_vars)))
    colnames(T_categorical) <- categorical_vars
    for (i in seq_along(categorical_vars)) {
      levels_i <- unique(data[[categorical_vars[i]]])
      T_categorical[[categorical_vars[i]]] <- sample(levels_i, nPoint, replace = TRUE)
    }
  } else {
    T_categorical <- NULL
  }

  # Combine continuous and categorical predictors
  if (!is.null(T_continuous) && !is.null(T_categorical)) {
    T_data <- cbind(as.data.frame(T_continuous), T_categorical)
  } else if (!is.null(T_continuous)) {
    T_data <- as.data.frame(T_continuous)
  } else if (!is.null(T_categorical)) {
    T_data <- T_categorical
  } else {
    stop("No predictors provided.")
  }

  # Create a formula without the response variable
  terms_obj <- terms(formula)
  rhs_formula <- reformulate(attr(terms_obj, "term.labels"))

  # Transform T_data into model matrix using the RHS-only formula
  T_mf <- model.frame(rhs_formula, data = T_data)
  T <- model.matrix(rhs_formula, data = T_mf)

  # Remove intercept column if present
  intercept_col_T <- which(colnames(T) == "(Intercept)")
  if (length(intercept_col_T) > 0) {
    T <- T[, -intercept_col_T, drop = FALSE]
  }

  y_mean <- mean(y)

  # Steps 2 and 3: Fit SVEMnet nSVEM times to get M_Y
  M_Y <- matrix(0, nrow = nSVEM, ncol = nPoint)
  if (verbose) {
    message("Fitting SVEM models to original data...")
  }
  for (i in 1:nSVEM) {
    svem_model <- tryCatch({
      SVEMnet(formula, data = data, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
              weight_scheme = weight_scheme, objective = objective,...)
    }, error = function(e) {
      message("Error in SVEMnet during SVEM fitting: ", e$message)
      NULL
    })

    if (is.null(svem_model)) next

    # Use predict() function
    pred_res <- predict(svem_model, newdata = T_data, debias = debias, se.fit = TRUE)
    f_hat_Y_T <- pred_res$fit
    s_hat_Y_T <- pred_res$se.fit
    # Handle zero standard errors to avoid division by zero
    s_hat_Y_T[s_hat_Y_T == 0] <- 1e-6
    h_Y <- (f_hat_Y_T - y_mean) / s_hat_Y_T
    M_Y[i, ] <- h_Y
  }

  # Steps 4 to 6: Fit SVEMnet to permuted data to get M_pi_Y
  M_pi_Y <- matrix(0, nrow = nPerm, ncol = nPoint)
  if (verbose) {
    message("Starting permutation testing...")
  }
  # Record the start time
  start_time_perm <- Sys.time()

  for (j in 1:nPerm) {
    y_perm <- sample(y, replace = FALSE)
    data_perm <- data
    data_perm[[as.character(formula[[2]])]] <- y_perm

    svem_model_perm <- tryCatch({
      SVEMnet(formula, data = data_perm, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
              weight_scheme = weight_scheme, objective = objective,...)
    }, error = function(e) {
      message("Error in SVEMnet during permutation fitting: ", e$message)
      NULL
    })

    if (is.null(svem_model_perm)) next

    # Use predict() function
    pred_res <- predict(svem_model_perm, newdata = T_data, debias = debias, se.fit = TRUE)
    f_hat_piY_T <- pred_res$fit
    s_hat_piY_T <- pred_res$se.fit
    # Handle zero standard errors
    s_hat_piY_T[s_hat_piY_T == 0] <- 1e-6
      h_piY <- (f_hat_piY_T - y_mean) / s_hat_piY_T
    M_pi_Y[j, ] <- h_piY

    # Progress Indicator
    if (verbose && (j %% 10 == 0 || j == nPerm)) {
      elapsed_time <- Sys.time() - start_time_perm
      elapsed_secs <- as.numeric(elapsed_time, units = "secs")
      estimated_total_secs <- (elapsed_secs / j) * nPerm
      remaining_secs <- estimated_total_secs - elapsed_secs

      # Format the remaining time
      remaining_time_formatted <- sprintf(
        "%02d:%02d:%02d",
        floor(remaining_secs / 3600),
        floor((remaining_secs %% 3600) / 60),
        floor(remaining_secs %% 60)
      )

      # Print the progress message
      message(sprintf("Permutation %d/%d completed. Estimated time remaining: %s",
                      j, nPerm, remaining_time_formatted))
    }
  }

  # Step 7: Standardize M_pi_Y
  col_means_M_pi_Y <- colMeans(M_pi_Y)
  col_sds_M_pi_Y <- apply(M_pi_Y, 2, sd)
  col_sds_M_pi_Y[col_sds_M_pi_Y == 0] <- 1e-6  # Handle zero standard deviations

    tilde_M_pi_Y <- scale(M_pi_Y, center = TRUE, scale = TRUE)

  # Step 8: Standardize M_Y using M_pi_Y stats
  M_Y_centered <- sweep(M_Y, 2, col_means_M_pi_Y, "-")

    tilde_M_Y <- sweep(M_Y_centered, 2, col_sds_M_pi_Y, "/")

  # Step 9: Perform SVD on tilde_M_pi_Y
  svd_res <- svd(tilde_M_pi_Y)
  U <- svd_res$u
  s <- svd_res$d
  V <- svd_res$v

  # Step 10: Compute scaled eigenvalues to sum to number of columns
  evalues_temp <- s^2
  evalues_temp <- evalues_temp / sum(evalues_temp) * ncol(tilde_M_pi_Y)

  # Step 11: Determine number of components to capture desired variance
  cumsum_evalues <- cumsum(evalues_temp) / sum(evalues_temp) * 100
  k <- which(cumsum_evalues > percent)[1]
  evalues <- evalues_temp[1:k]
  evectors <- V[, 1:k]

  # Step 12: Calculate Mahalanobis distances for permutations
  T2_perm <- diag(tilde_M_pi_Y %*% evectors %*% diag(evalues^(-1)) %*% t(evectors) %*% t(tilde_M_pi_Y))
  d_pi_Y <- sqrt(T2_perm)

  # Step 13: Calculate Mahalanobis distances for original data
  T2_Y <- diag(tilde_M_Y %*% evectors %*% diag(evalues^(-1)) %*% t(evectors) %*% t(tilde_M_Y))
  d_Y <- sqrt(T2_Y)

  # Step 14: Fit the specified null distribution to d_pi_Y
  if (length(d_pi_Y) == 0) {
    stop("No valid d_pi_Y values to fit a distribution.")
  }


    suppressMessages(
      distribution_fit <- tryCatch({
        gamlss::gamlss(
          d_pi_Y ~ 1,
          family = gamlss.dist::SHASHo(mu.link = "identity", sigma.link = "log", nu.link = "identity",
                                       tau.link = "log"),
          control = gamlss::gamlss.control(n.cyc = 1000, trace = FALSE)
        )
      }, error = function(e) {
        message("Error in fitting SHASHo distribution: ", e$message)
        NULL
      })
    )
    if (is.null(distribution_fit)) stop("Failed to fit SHASHo distribution.")

    # Extract parameters
    mu <- as.numeric(coef(distribution_fit, what = "mu"))
    sigma <- exp(as.numeric(coef(distribution_fit, what = "sigma")))
    nu <- as.numeric(coef(distribution_fit, what = "nu"))
    tau <- exp(as.numeric(coef(distribution_fit, what = "tau")))

    # Compute p-values using the pSHASHo function
    p_values <- 1 - gamlss.dist::pSHASHo(d_Y, mu = mu, sigma = sigma, nu = nu, tau = tau)
    p_value <- median(p_values)



  # **Include the response variable name**
  response_name <- as.character(formula[[2]])

  # **Combine distances into a data frame for plotting**
  data_d <- data.frame(
    D = c(d_Y, d_pi_Y),
    Source_Type = c(rep("Original", length(d_Y)), rep("Permutation", length(d_pi_Y))),
    Response = response_name
  )

  # Return results
  results_list <- list(
    p_value = p_value,
    p_values = p_values,
    d_Y = d_Y,
    d_pi_Y = d_pi_Y,
    distribution_fit = distribution_fit,
    data_d = data_d  # Include data for plotting
  )

  # Assign class to the results_list
  class(results_list) <- "svem_significance_test"

  return(results_list)
}
