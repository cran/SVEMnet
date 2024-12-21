# Load necessary libraries
library(SVEMnet)
library(lhs)
library(ggplot2)
library(dplyr)
library(pbapply)
library(parallel)
library(tidyr)  # For pivot_wider

# Define vectors for p, d, n, sd
p_values <- seq(3, 6, 1)        # Number of parameters
d_values <- seq(0.1, 0.9, 0.1)  # Density (proportion of active parameters)
n_values <- seq(15, 50, 5)      # Number of design points
sd_values <- c(0.25, 0.5, 1, 1.5)  # Standard deviations of noise

nSim <- 20                  # Number of simulations per setting
n_holdout <- 1000           # Number of holdout points

# Create a grid of all combinations of p, d, n, sd
param_grid <- expand.grid(p = p_values, d = d_values, n = n_values, sd = sd_values)

# Prepare a list of simulation parameters
sim_params_list <- list()
sim_id <- 1  # Initialize simulation ID

for (i in 1:nrow(param_grid)) {
  p <- param_grid$p[i]
  d <- param_grid$d[i]
  n <- param_grid$n[i]
  sd <- param_grid$sd[i]

  for (sim in 1:nSim) {
    sim_params_list[[length(sim_params_list) + 1]] <- list(
      sim_id = sim_id,
      p = p,
      d = d,
      n = n,
      sd = sd
    )
    sim_id <- sim_id + 1
  }
}

total_iterations <- length(sim_params_list)
cat("Total simulations to run:", total_iterations, "\n")

# Set up parallel backend using 'parallel' package
num_cores <- 8
RNGkind("L'Ecuyer-CMRG")
set.seed(0)
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(SVEMnet)
  library(lhs)
  library(dplyr)
})

# Export necessary variables to the cluster
clusterExport(cl, varlist = c("n_holdout"))

# Set up RNG for parallel processing

clusterSetRNGStream(cl, 0)

# Function to perform one simulation for a given setting
simulate_one <- function(sim_params) {
  sim_id <- sim_params$sim_id
  p <- sim_params$p
  d <- sim_params$d
  n <- sim_params$n
  sd <- sim_params$sd

  # Generate design data X using LHS
  X <- randomLHS(n, p)
  colnames(X) <- paste0("V", 1:p)

  # Generate holdout data X_T using LHS
  X_T <- randomLHS(n_holdout, p)
  colnames(X_T) <- paste0("V", 1:p)

  # Select active parameters
  n_active <- max(1, floor(p * d))
  active_params <- sample(1:p, n_active)

  # Generate coefficients
  beta <- numeric(p)
  beta[active_params] <- rexp(n_active) - rexp(n_active)

  # Generate response variable y
  y <- X %*% beta + rnorm(n, mean = 0, sd = sd)
  y <- as.numeric(y)  # Convert to vector

  # Compute true y for holdout data X_T
  y_T_true <- X_T %*% beta
  y_T_true <- as.numeric(y_T_true)

  # Define formula for SVEMnet
  data_train <- data.frame(y = y, X)
  data_holdout <- data.frame(X_T)
  colnames(data_holdout) <- colnames(X)

  # Include interactions and quadratic terms
  formula <- as.formula(paste(
    "y ~ (", paste(colnames(X), collapse = " + "), ")^2 + ",
    paste("I(", colnames(X), "^2)", collapse = " + ")
  ))

  # Initialize data frame to store results for this simulation
  sim_results <- data.frame()

  # Set debias = FALSE and objective = "wAIC"
  debias <- FALSE
  objective <- "wAIC"

  # Loop over model types
  for (model_type in c("SVEM", "FWR", "Identity")) {
    # Set weight_scheme and nBoot
    if (model_type == "SVEM") {
      weight_scheme <- "SVEM"
      nBoot <- 200
    } else if (model_type == "FWR") {
      weight_scheme <- "FWR"
      nBoot <- 200
    } else if (model_type == "Identity") {
      weight_scheme <- "Identity"
      nBoot <- 1
    }

    # Fit SVEMnet model
    model <- SVEMnet(
      formula = formula,
      data = data_train,
      objective = objective,
      nBoot = nBoot,
      glmnet_alpha = c(1),  # Lasso
      weight_scheme = weight_scheme
    )

    # Predict on holdout data
    y_pred_T <- predict(
      object = model,
      newdata = data_holdout,
      debias = debias
    )

    # Compute RMSE over X_T
    RMSE <- sqrt(mean((y_T_true - y_pred_T)^2))
    logRMSE <- log(RMSE)

    # Record results
    sim_results <- rbind(
      sim_results,
      data.frame(
        sim_id = sim_id,
        p = p,
        d = d,
        n = n,
        sd = sd,
        model_type = model_type,
        logRMSE = logRMSE
      )
    )
  }

  return(sim_results)
}

# Run simulations using pblapply with progress bar
results_list <- pblapply(sim_params_list, simulate_one, cl = cl)

# Stop the cluster
stopCluster(cl)

# Combine all results into a data frame
results <- bind_rows(results_list)

# Convert sim_id to integer
results$sim_id <- as.integer(as.character(results$sim_id))

# Compute the mean logRMSE for each simulation ID
results <- results %>%
  group_by(sim_id) %>%
  mutate(sim_mean_logRMSE = mean(logRMSE)) %>%
  ungroup()

# Compute residuals
results <- results %>%
  mutate(residual_logRMSE = logRMSE - sim_mean_logRMSE)

# Compute average residuals for each model_type
average_residuals <- results %>%
  group_by(model_type) %>%
  summarise(mean_residual_logRMSE = mean(residual_logRMSE),
            .groups = 'drop')

# Print average residuals
print(average_residuals)

# Plot residuals side by side for the three models
results$model_type <- factor(results$model_type, levels = c("SVEM", "FWR", "Identity"))

ggplot(results, aes(x = model_type, y = logRMSE)) +
  geom_boxplot() +
  xlab("Model Type") +
  ylab("LogRMSE") +
  ggtitle("LogRMSE by Model") +
  theme_minimal()

ggplot(results, aes(x = model_type, y = residual_logRMSE)) +
  geom_boxplot() +
  xlab("Model Type") +
  ylab("Residual LogRMSE") +
  ggtitle("Residuals (after removing simulation mean logRMSE)") +
  theme_minimal()

# Compute paired differences between models for each simulation
results_wide <- results %>%
  select(sim_id, model_type, logRMSE) %>%
  pivot_wider(names_from = model_type, values_from = logRMSE)

# Compute differences
results_wide <- results_wide %>%
  mutate(
    diff_SVEM_FWR = SVEM - FWR,
    diff_SVEM_Identity = SVEM - Identity,
    diff_FWR_Identity = FWR - Identity
  )

# Summarize differences
differences_summary <- results_wide %>%
  summarise(
    mean_diff_SVEM_FWR = mean(diff_SVEM_FWR, na.rm = TRUE),
    mean_diff_SVEM_Identity = mean(diff_SVEM_Identity, na.rm = TRUE),
    mean_diff_FWR_Identity = mean(diff_FWR_Identity, na.rm = TRUE),
    sd_diff_SVEM_FWR = sd(diff_SVEM_FWR, na.rm = TRUE),
    sd_diff_SVEM_Identity = sd(diff_SVEM_Identity, na.rm = TRUE),
    sd_diff_FWR_Identity = sd(diff_FWR_Identity, na.rm = TRUE)
  )

print(differences_summary)

# Optionally, perform paired t-tests
t_test_SVEM_FWR <- t.test(results_wide$SVEM, results_wide$FWR, paired = TRUE)
t_test_SVEM_Identity <- t.test(results_wide$SVEM, results_wide$Identity, paired = TRUE)
t_test_FWR_Identity <- t.test(results_wide$FWR, results_wide$Identity, paired = TRUE)

# Print t-test results
print(t_test_SVEM_FWR)
print(t_test_SVEM_Identity)
print(t_test_FWR_Identity)
