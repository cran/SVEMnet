# Install if not already installed
# install.packages("gtools")

library(SVEMnet)
library(ggplot2)
library(dplyr)
library(parallel)
library(tidyr)
library(gtools)  # for rdirichlet

# Define vectors for p, d, n, sd, and outlier_prop
p_values <- seq(4, 7, 1)                # Number of parameters (mixture components)
d_values <- seq(.1, .9, .1)             # Density (proportion of terms to be active)
n_values <- seq(20, 70, 10)            # Number of design points
sd_values <- c(0.25, 0.5, 1)       # Standard deviations of noise
outlier_prop_values <- c(0,0.1)    # Proportions of outliers

nSim <- 5           # Number of simulations per setting
n_holdout <- 2000         # Number of holdout points
mult <- 200               # Oversampling multiple for candidate set generation

# Create a grid of all combinations of p, d, n, sd, outlier_prop
param_grid <- expand.grid(
  p = p_values,
  d = d_values,
  n = n_values,
  sd = sd_values,
  outlier_prop = outlier_prop_values
)

# Prepare a list of simulation parameters
sim_params_list <- list()
sim_id <- 1

for (i in 1:nrow(param_grid)) {
  p <- param_grid$p[i]
  d <- param_grid$d[i]
  n <- param_grid$n[i]
  sd <- param_grid$sd[i]
  outlier_prop <- param_grid$outlier_prop[i]

  for (sim in 1:nSim) {
    sim_params_list[[length(sim_params_list) + 1]] <- list(
      sim_id = sim_id,
      p = p,
      d = d,
      n = n,
      sd = sd,
      outlier_prop = outlier_prop
    )
    sim_id <- sim_id + 1
  }
}

total_iterations <- length(sim_params_list)
cat("Total simulations to run:", total_iterations, "\n")

# Set up parallel backend
num_cores <- 10
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 1)
clusterEvalQ(cl, {
  library(SVEMnet)
  library(dplyr)
  library(gtools)
})

clusterExport(cl, varlist = c("n_holdout","mult"))

simulate_one <- function(sim_params) {
  sim_id <- sim_params$sim_id
  p <- sim_params$p
  d <- sim_params$d
  n <- sim_params$n
  sd <- sim_params$sd
  outlier_prop <- sim_params$outlier_prop

  # Step 1: Generate a large candidate set of points in the simplex
  N <- n * mult
  candidates <- rdirichlet(N, rep(1, p))

  # Step 2: k-means to get n cluster centroids as final design
  km_result <- kmeans(candidates, centers = n, nstart = 10)

  # Final mixture design
  X <- km_result$centers
  colnames(X) <- paste0("V", 1:p)

  # Use a fourth-order model: (X1 + X2 + ... + Xp)^4
  formula_str <- paste("y ~ (", paste(colnames(X), collapse = " + "), ")^",as.character(3))
  formula <- as.formula(formula_str)

  #model formula
  formula_str_model <- paste("y ~ (", paste(colnames(X), collapse = " + "), ")^3")
  formula_model <- as.formula(formula_str)

  # Temporary data frame with y=0 just to build model matrix
  data_train <- data.frame(y = 0, X)
  mf <- model.frame(formula, data_train)
  X_full <- model.matrix(formula, mf)

  # Count how many terms (excluding intercept)
  M <- ncol(X_full) - 1
  n_active <- max(1, floor(M * d))

  # Select active terms from all non-intercept terms
  active_terms <- sample(2:(M+1), n_active, replace = FALSE)

  # Create coefficients
  beta <- numeric(M+1)
  beta[active_terms] <- rexp(n_active) - rexp(n_active)

  # Generate response y
  data_train$y <- drop(X_full %*% beta + rnorm(n, mean = 0, sd = sd))

  # Introduce outliers
  num_outliers <- ceiling(outlier_prop * n)
  if (num_outliers > 0) {
    outlier_indices <- sample(seq_len(n), num_outliers, replace = FALSE)
    data_train$y[outlier_indices] <- data_train$y[outlier_indices] + rnorm(num_outliers, mean = 0, sd = 3 * sd)
  }

  # Holdout data
  X_T <- rdirichlet(n_holdout, rep(1, p))
  colnames(X_T) <- paste0("V", 1:p)
  data_holdout <- data.frame(X_T)

  # Remove response for holdout
  terms_obj <- terms(formula, data = data_train)
  terms_obj_no_resp <- delete.response(terms_obj)
  mf_T <- model.frame(terms_obj_no_resp, data_holdout)
  X_T_full <- model.matrix(terms_obj_no_resp, mf_T)
  y_T_true <- drop(X_T_full %*% beta)

  # Fit models
  model_svem <- SVEMnet(
    formula = formula_model,
    data = data_train,
    standardize = TRUE,
    objective = "wAIC",
    nBoot = 200,
    glmnet_alpha = c(0,0.5,1),
    weight_scheme = "SVEM"
  )

  model_glmnet_cv <- glmnet_with_cv(
    formula = formula_model,
    data = data_train,
    glmnet_alpha = c(0,0.5,1),
    standardize = TRUE
  )

  model_svem_identity <- SVEMnet(
    formula = formula_model,
    data = data_train,
    objective = "wAIC",
    standardize=TRUE,
    nBoot = 1,
    glmnet_alpha = c(0,0.5,1),
    weight_scheme = "Identity"
  )

  pred_list <- list()

  # SVEMnet debias=FALSE
  y_pred_svem_false <- predict(model_svem, newdata = data_holdout, debias = FALSE)
  RMSE_svem_false <- sqrt(mean((y_T_true - y_pred_svem_false)^2))
  logRMSE_svem_false <- log(RMSE_svem_false)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "SVEMnet", debias = FALSE, logRMSE = logRMSE_svem_false
  )

  # SVEMnet debias=TRUE
  y_pred_svem_true <- predict(model_svem, newdata = data_holdout, debias = TRUE)
  RMSE_svem_true <- sqrt(mean((y_T_true - y_pred_svem_true)^2))
  logRMSE_svem_true <- log(RMSE_svem_true)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "SVEMnet", debias = TRUE, logRMSE = logRMSE_svem_true
  )

  # glmnet_cv debias=FALSE
  y_pred_glmnet_false <- predict_cv(model_glmnet_cv, newdata = data_holdout, debias = FALSE)
  RMSE_glmnet_false <- sqrt(mean((y_T_true - y_pred_glmnet_false)^2))
  logRMSE_glmnet_false <- log(RMSE_glmnet_false)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "glmnet_cv", debias = FALSE, logRMSE = logRMSE_glmnet_false
  )

  # glmnet_cv debias=TRUE
  y_pred_glmnet_true <- predict_cv(model_glmnet_cv, newdata = data_holdout, debias = TRUE)
  RMSE_glmnet_true <- sqrt(mean((y_T_true - y_pred_glmnet_true)^2))
  logRMSE_glmnet_true <- log(RMSE_glmnet_true)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "glmnet_cv", debias = TRUE, logRMSE = logRMSE_glmnet_true
  )

  # SVEMnet_Identity debias=FALSE
  y_pred_svem_identity_false <- predict(model_svem_identity, newdata = data_holdout, debias = FALSE)
  RMSE_svem_identity_false <- sqrt(mean((y_T_true - y_pred_svem_identity_false)^2))
  logRMSE_svem_identity_false <- log(RMSE_svem_identity_false)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "SVEMnet_Identity", debias = FALSE, logRMSE = logRMSE_svem_identity_false
  )

  # SVEMnet_Identity debias=TRUE
  y_pred_svem_identity_true <- predict(model_svem_identity, newdata = data_holdout, debias = TRUE)
  RMSE_svem_identity_true <- sqrt(mean((y_T_true - y_pred_svem_identity_true)^2))
  logRMSE_svem_identity_true <- log(RMSE_svem_identity_true)
  pred_list[[length(pred_list)+1]] <- data.frame(
    sim_id = sim_id, p = p, d = d, n = n, sd = sd, outlier_prop = outlier_prop,
    model = "SVEMnet_Identity", debias = TRUE, logRMSE = logRMSE_svem_identity_true
  )

  sim_results <- bind_rows(pred_list)
  return(sim_results)
}

# Run load-balanced parallel computations
results_list <- parLapplyLB(cl, sim_params_list, simulate_one)
stopCluster(cl)

results <- bind_rows(results_list)

results <- results %>%
  group_by(sim_id) %>%
  mutate(sim_mean_logRMSE = mean(logRMSE)) %>%
  ungroup() %>%
  mutate(residual_logRMSE = logRMSE - sim_mean_logRMSE)

average_residuals <- results %>%
  group_by(model, debias) %>%
  summarise(mean_residual_logRMSE = mean(residual_logRMSE), .groups = 'drop')

print(average_residuals)

results$model_debias <- paste0(results$model, "_debias=", results$debias)
results$model_debias <- factor(results$model_debias,
                               levels = c("SVEMnet_debias=FALSE", "SVEMnet_debias=TRUE",
                                          "glmnet_cv_debias=FALSE", "glmnet_cv_debias=TRUE",
                                          "SVEMnet_Identity_debias=FALSE", "SVEMnet_Identity_debias=TRUE"))

# Your ggplot calls:
ggplot(results, aes(y = model_debias, x = logRMSE, fill = as.factor(debias))) +
  geom_boxplot() +
  ylab("Model & Debias Setting") +
  xlab("LogRMSE") +
  ggtitle("LogRMSE by Model and Debias Setting") +
  theme_minimal() +
  coord_flip() +
  scale_fill_discrete(name = "Debias", labels = c("FALSE", "TRUE"))

ggplot(results, aes(y = model_debias, x = residual_logRMSE, fill = as.factor(debias))) +
  geom_boxplot() +
  ylab("Model & Debias Setting") +
  xlab("Residual LogRMSE") +
  ggtitle("Residuals (after removing simulation mean)") +
  theme_minimal() +
  coord_flip() +
  scale_fill_discrete(name = "Debias", labels = c("FALSE", "TRUE"))

results_wide <- results %>%
  dplyr::select(sim_id, model_debias, logRMSE) %>%
  pivot_wider(names_from = model_debias, values_from = logRMSE)

results_wide <- results_wide %>%
  mutate(
    diff_SVEM_wAIC_debiasTRUE_FALSE = `SVEMnet_debias=TRUE` - `SVEMnet_debias=FALSE`,
    diff_glmnet_cv_debiasTRUE_FALSE = `glmnet_cv_debias=TRUE` - `glmnet_cv_debias=FALSE`,
    diff_SVEM_glmnet_cv_FALSE = `SVEMnet_debias=FALSE` - `glmnet_cv_debias=FALSE`,
    diff_SVEM_glmnet_cv_TRUE  = `SVEMnet_debias=TRUE` - `glmnet_cv_debias=TRUE`
  )

differences_summary <- results_wide %>%
  summarise(
    mean_diff_SVEM_debias = mean(diff_SVEM_wAIC_debiasTRUE_FALSE, na.rm = TRUE),
    sd_diff_SVEM_debias = sd(diff_SVEM_wAIC_debiasTRUE_FALSE, na.rm = TRUE),
    mean_diff_glmnet_debias = mean(diff_glmnet_cv_debiasTRUE_FALSE, na.rm = TRUE),
    sd_diff_glmnet_debias = sd(diff_glmnet_cv_debiasTRUE_FALSE, na.rm = TRUE),
    mean_diff_SVEM_glmnet_F = mean(diff_SVEM_glmnet_cv_FALSE, na.rm = TRUE),
    sd_diff_SVEM_glmnet_F = sd(diff_SVEM_glmnet_cv_FALSE, na.rm = TRUE),
    mean_diff_SVEM_glmnet_T = mean(diff_SVEM_glmnet_cv_TRUE, na.rm = TRUE),
    sd_diff_SVEM_glmnet_T = sd(diff_SVEM_glmnet_cv_TRUE, na.rm = TRUE)
  )

print(differences_summary)

t_test_SVEM_debias <- t.test(results_wide$`SVEMnet_debias=TRUE`, results_wide$`SVEMnet_debias=FALSE`, paired = TRUE)
t_test_glmnet_debias <- t.test(results_wide$`glmnet_cv_debias=TRUE`, results_wide$`glmnet_cv_debias=FALSE`, paired = TRUE)
t_test_SVEM_vs_glmnet_FALSE <- t.test(results_wide$`SVEMnet_debias=FALSE`, results_wide$`glmnet_cv_debias=FALSE`, paired = TRUE)
t_test_SVEM_vs_glmnet_TRUE  <- t.test(results_wide$`SVEMnet_debias=TRUE`, results_wide$`glmnet_cv_debias=TRUE`, paired = TRUE)

print(t_test_SVEM_debias)
print(t_test_glmnet_debias)
print(t_test_SVEM_vs_glmnet_FALSE)
print(t_test_SVEM_vs_glmnet_TRUE)
