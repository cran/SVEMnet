#' SVEM Significance Test (Parallel Version)
#'
#' Performs a whole-model significance test using the SVEM framework, handling both continuous and categorical predictors, with parallel computation.
#'
#' @param formula A formula specifying the model to be tested.
#' @param data A data frame containing the variables in the model.
#' @param nPoint The number of random points to generate in the factor space (default: 2000).
#' @param nSVEM The number of SVEM models to fit to the original data (default: 7).
#' @param nPerm The number of SVEM models to fit to permuted data for reference distribution (default: 200).
#' @param percent The percentage of variance to capture in the SVD (default: 85).
#' @param nBoot The number of bootstrap iterations within SVEM (default: 200).
#' @param glmnet_alpha The alpha parameter(s) for glmnet (default: \code{c(1)}).
#' @param weight_scheme The weight scheme to use in SVEM (default: "SVEM").
#' @param objective Character; the objective function to use in \code{\link{SVEMnet}}. Options are "wAIC" or "wSSE" (default: "wAIC").
#' @param verbose Logical; if \code{TRUE}, displays progress messages (default: \code{TRUE}).
#' @param nCore The number of CPU cores to use for parallel processing (default: all available cores).
#' @param seed An integer seed for random number generation (default: NULL).
#' @param ... Additional arguments passed to the underlying \code{SVEMnet()} and then \code{glmnet()} functions.
#' @return A list containing the test results.
#' @details
#' The `svem_significance_test_parallel` function implements a whole-model test designed to gauge the significance of a fitted SVEM model compared to the null hypothesis of a constant response surface, with parallel computation. This method helps identify responses that have relatively stronger or weaker relationships with study factors.
#'
#' The test constructs standardized predictions by centering the SVEM predictions (obtained from \code{SVEMnet()}) by the response mean and scaling by the ensemble standard deviation. A reference distribution is created by fitting the SVEM model to multiple randomized permutations of the response vector. The Mahalanobis distances of the original and permuted models are calculated using a reduced-rank singular value decomposition.
#'
#' The R code to perform this test (using matrices of nSVEM and nPerm predictions) is taken from the supplementary material of Karl (2024).
#'
#' This function assumes that there are no restrictions among the factors (e.g. mixture factors). The method will work with restrictions, but the code would need to be changed to ensure the \code{nPoint} points respect the factor restriction(s). For example, \code{rdirichlet()} could be used for the mixture factors.
#'
#' The \code{SVEMnet} parameter \code{debias} is hard coded to \code{FALSE} for this test. Unpublished simulation work suggests that setting \code{debias=TRUE} reduces the power of the test (without affecting the Type I error rate).
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024).
#'
#' @references
#' Karl, A. T. (2024). A randomized permutation whole-model test heuristic for Self-Validated Ensemble Models (SVEM). \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{249}, 105122. \doi{10.1016/j.chemolab.2024.105122}
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(0)
#' n <- 30
#' X1 <- runif(n)
#' X2 <- runif(n)
#' X3 <- runif(n)
#' y <- 1 + X1 +  X2 + X1 * X2 + X1^2 + rnorm(n)
#' data <- data.frame(y, X1, X2, X3)
#'
#' #CRAN requires a max of nCore=2 for example. Recommend using default nCore to use entire CPU.
#'
#' # Perform the SVEM significance test
#' test_result <- svem_significance_test_parallel(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2) + I(X3^2),
#'   data = data,
#'   nPoint = 2000,
#'   nSVEM = 9,
#'   nPerm = 250,
#'   nBoot = 200,
#'   nCore = 2
#' )
#'
#' # View the p-value
#' print(test_result)
#'
#'
#' test_result2 <- svem_significance_test_parallel(
#'   y ~ (X1 + X2 )^2 + I(X1^2) + I(X2^2),
#'   data = data,
#'   nPoint = 2000,
#'   nSVEM = 9,
#'   nPerm = 250,
#'   nBoot = 200,
#'   nCore = 2
#' )
#'
#' # View the p-value
#' print(test_result2)
#'
#'
#' # Plot the Mahalanobis distances
#' plot(test_result,test_result2)
#' }
#'
#' @importFrom lhs maximinLHS
#' @importFrom gamlss gamlss gamlss.control
#' @importFrom gamlss.dist SHASHo pSHASHo
#' @importFrom stats model.response delete.response terms reformulate median complete.cases
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterSetRNGStream detectCores
#' @export
svem_significance_test_parallel <- function(formula, data, nPoint = 2000, nSVEM = 7, nPerm = 200,
                                            percent = 85, nBoot = 200, glmnet_alpha = c(1),
                                            weight_scheme = c("SVEM"),
                                            objective = c("wAIC", "wSSE"), verbose = TRUE,
                                            nCore = parallel::detectCores(), seed = NULL, ...) {
  # Match the arguments
  objective <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)

  # Set RNGkind
  RNGkind("L'Ecuyer-CMRG")

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Set up the cluster
  cl <- parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))  # Ensure cluster is stopped when function exits

  # Distribute RNG streams to each worker
  parallel::clusterSetRNGStream(cl, iseed = seed)

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
    T_continuous_raw <- as.matrix(lhs::maximinLHS(nPoint, length(continuous_vars)))
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

  # Create a formula without the response variable, including data argument
  terms_obj <- terms(formula, data = data)
  environment(terms_obj) <- baseenv()  # Set environment to baseenv()
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
  if (verbose) {
    message("Fitting SVEM models to original data...")
  }

  M_Y <- foreach::foreach(i = 1:nSVEM, .combine = rbind,
                          .packages = c("SVEMnet", "stats"), .export = c("SVEMnet")) %dopar% {
                            svem_model <- tryCatch({
                              SVEMnet::SVEMnet(formula, data = data, nBoot = nBoot, glmnet_alpha = glmnet_alpha,
                                               weight_scheme = weight_scheme, objective = objective, ...)
                            }, error = function(e) {
                              message("Error in SVEMnet during SVEM fitting: ", e$message)
                              return(NULL)
                            })

                            if (is.null(svem_model)) {
                              # Return a vector of NAs to maintain the structure
                              return(rep(NA, nPoint))
                            }

                            pred_res <- predict(svem_model, newdata = T_data, debias = FALSE, se.fit = TRUE)
                            f_hat_Y_T <- pred_res$fit
                            s_hat_Y_T <- pred_res$se.fit
                            s_hat_Y_T[s_hat_Y_T == 0] <- 1e-6
                            h_Y <- (f_hat_Y_T - y_mean) / s_hat_Y_T
                            return(h_Y)
                          }

  # Steps 4 to 6: Fit SVEMnet to permuted data to get M_pi_Y
  if (verbose) {
    message("Starting permutation testing...")
  }
  # Record the start time
  start_time_perm <- Sys.time()

  M_pi_Y <- foreach::foreach(jloop = 1:nPerm, .combine = rbind,
                             .packages = c("SVEMnet", "stats"), .export = c("SVEMnet")) %dopar% {
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

                               if (is.null(svem_model_perm)) {
                                 # Return a vector of NAs to maintain the structure
                                 return(rep(NA, nPoint))
                               }

                               pred_res <- predict(svem_model_perm, newdata = T_data, debias = FALSE, se.fit = TRUE)
                               f_hat_piY_T <- pred_res$fit
                               s_hat_piY_T <- pred_res$se.fit
                               s_hat_piY_T[s_hat_piY_T == 0] <- 1e-6
                               h_piY <- (f_hat_piY_T - y_mean) / s_hat_piY_T

                               # Progress Indicator
                               if (verbose && (jloop %% 10 == 0 || jloop == nPerm)) {
                                 elapsed_time <- Sys.time() - start_time_perm
                                 elapsed_secs <- as.numeric(elapsed_time, units = "secs")
                                 estimated_total_secs <- (elapsed_secs / jloop) * nPerm
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
                                                 jloop, nPerm, remaining_time_formatted))
                               }

                               return(h_piY)
                             }

  # Replace any rows with NAs (from failed iterations)
  M_Y <- M_Y[complete.cases(M_Y), , drop = FALSE]
  M_pi_Y <- M_pi_Y[complete.cases(M_pi_Y), , drop = FALSE]

  if (nrow(M_Y) == 0) {
    stop("All SVEM model fittings on the original data failed.")
  }
  if (nrow(M_pi_Y) == 0) {
    stop("All SVEM model fittings on the permuted data failed.")
  }

  # Step 7: Standardize M_pi_Y
  col_means_M_pi_Y <- colMeans(M_pi_Y, na.rm = TRUE)
  col_sds_M_pi_Y <- apply(M_pi_Y, 2, sd, na.rm = TRUE)
  col_sds_M_pi_Y[col_sds_M_pi_Y == 0] <- 1e-6  # Handle zero standard deviations

  tilde_M_pi_Y <- scale(M_pi_Y, center = TRUE, scale = TRUE)

  # Step 8: Standardize M_Y using M_pi_Y stats
  M_Y_centered <- sweep(M_Y, 2, col_means_M_pi_Y, "-")
  tilde_M_Y <- sweep(M_Y_centered, 2, col_sds_M_pi_Y, "/")

  # Step 9: Perform SVD on tilde_M_pi_Y
  svd_res <- svd(tilde_M_pi_Y)
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
  T2_perm <- rowSums((tilde_M_pi_Y %*% evectors %*% diag(1 / evalues)) * (tilde_M_pi_Y %*% evectors))
  d_pi_Y <- sqrt(T2_perm)

  # Step 13: Calculate Mahalanobis distances for original data
  T2_Y <- rowSums((tilde_M_Y %*% evectors %*% diag(1 / evalues)) * (tilde_M_Y %*% evectors))
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

  # Include the response variable name
  response_name <- as.character(formula[[2]])

  # Combine distances into a data frame for plotting
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
