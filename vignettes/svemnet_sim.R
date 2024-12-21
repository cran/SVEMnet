# Load necessary libraries
library(SVEMnet)
library(lhs)
library(ggplot2)
library(dplyr)
library(pbapply)
library(parallel)

# Define vectors for p, d, n, sd
p_values <- seq(3, 6, 1)        # Number of parameters
d_values <- seq(0.1, 0.9, 0.1)  # Density (proportion of active parameters)
n_values <- seq(15, 50, 5)      # Number of design points
sd_values <- c(.25,0.5, 1, 1.5)       # Standard deviations of noise

nSim <- 20                  # Number of simulations per setting
n_holdout <- 1000               # Number of holdout points

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
set.seed(1234)
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(SVEMnet)
  library(lhs)
  library(dplyr)
})

# Export necessary variables to the cluster
clusterExport(cl, varlist = c("n_holdout"))

# Set up RNG for parallel processing

clusterSetRNGStream(cl, 1234)

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

  # Loop over objectives and debias options
  for (objective in c("wAIC", "wSSE")) {
    for (debias in c(TRUE, FALSE)) {
      # Fit SVEMnet model
      model <- SVEMnet(
        formula = formula,
        data = data_train,
        objective = objective,
        nBoot = 200,
        glmnet_alpha = c(1),  # Lasso
        weight_scheme = "SVEM"
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
          objective = objective,
          debias = debias,
          logRMSE = logRMSE
        )
      )
    }
  }

  return(sim_results)
}

# Run simulations using pblapply with progress bar
results_list <- pblapply(sim_params_list, simulate_one, cl = cl)

# Stop the cluster
stopCluster(cl)

# Combine all results into a data frame
results <- bind_rows(results_list)

# Convert sim_id to integer (since it may come back as a factor)
results$sim_id <- as.integer(as.character(results$sim_id))

# Compute the mean logRMSE for each simulation ID
results <- results %>%
  group_by(sim_id) %>%
  mutate(sim_mean_logRMSE = mean(logRMSE)) %>%
  ungroup()

# Compute residuals
results <- results %>%
  mutate(residual_logRMSE = logRMSE - sim_mean_logRMSE)

# Create a new variable that concatenates 'objective' and 'debias'
results <- results %>%
  mutate(obj_debias = paste0(objective, "_", debias))

# Compute average residuals for each combination of 'obj_debias'
average_residuals <- results %>%
  group_by(obj_debias) %>%
  summarise(mean_residual_logRMSE = mean(residual_logRMSE),
            .groups = 'drop')

# Print average residuals
print(average_residuals)

# Plot residuals side by side for the four combinations
# Set factor levels for desired order
results$obj_debias <- factor(results$obj_debias,
                             levels = c("wAIC_TRUE", "wAIC_FALSE", "wSSE_TRUE", "wSSE_FALSE"))

ggplot(results, aes(x = obj_debias, y = logRMSE)) +
  geom_boxplot() +
  xlab("Objective_Debias Combination") +
  ylab("LogRMSE") +
  ggtitle("logRMSE over all simulations") +
  theme_minimal()


ggplot(results, aes(x = obj_debias, y = residual_logRMSE)) +
  geom_boxplot() +
  xlab("Objective_Debias Combination") +
  ylab("Residual LogRMSE") +
  ggtitle("Residuals (after removing simulation mean logRMSE)") +
  theme_minimal()
