# ===== svemnet_sim_mixture_v2.R  =====
library(SVEMnet)
library(ggplot2)
library(dplyr)
library(parallel)
library(tidyr)
library(gtools)  # rdirichlet

p_values <- seq(4, 7, 1)
d_values <- seq(0.1, 0.9, 0.2)
n_values <- seq(20, 50, 10)
sd_values <- c(0.25, 0.5, 1)
outlier_prop_values <- c(0, 0.1)

nSim      <- 5
n_holdout <- 1000
mult      <- 200

param_grid <- expand.grid(p = p_values, d = d_values, n = n_values, sd = sd_values, outlier_prop = outlier_prop_values)

sim_params_list <- list(); sim_id <- 1L
for (i in seq_len(nrow(param_grid))) {
  for (sim in seq_len(nSim)) {
    sim_params_list[[length(sim_params_list) + 1L]] <- cbind(param_grid[i, ], sim_id = sim_id)
    sim_id <- sim_id + 1L
  }
}
cat("Total simulations to run:", length(sim_params_list), "\n")

num_cores <- max(1, parallel::detectCores(logical = FALSE) - 3)
RNGkind("L'Ecuyer-CMRG"); set.seed(1)
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 1)
clusterEvalQ(cl, { library(SVEMnet); library(dplyr); library(gtools) })
clusterExport(cl, varlist = c("n_holdout","mult"))

simulate_one <- function(sim_params) {
  p <- sim_params$p; d <- sim_params$d; n <- sim_params$n; sd <- sim_params$sd; outlier_prop <- sim_params$outlier_prop; sim_id <- sim_params$sim_id

  N <- n * mult
  candidates <- gtools::rdirichlet(N, rep(1, p))
  km <- kmeans(candidates, centers = n, nstart = 10)
  X <- km$centers; colnames(X) <- paste0("V", seq_len(p))

  # Truth: 4th-order; Model: 3rd-order (intentional misspecification)
  formula_truth_str <- paste("y ~ (", paste(colnames(X), collapse = " + "), ")^4")
  formula_model_str <- paste("y ~ (", paste(colnames(X), collapse = " + "), ")^3")
  formula_truth <- as.formula(formula_truth_str)
  formula_model <- as.formula(formula_model_str)

  df_train <- data.frame(y = 0, X)
  X_full <- model.matrix(formula_truth, model.frame(formula_truth, df_train))
  M <- ncol(X_full) - 1L
  n_active <- max(1, floor(M * d))

  active_terms <- sample(2:(M+1), n_active, replace = FALSE)
  beta <- numeric(M + 1L); beta[active_terms] <- rexp(n_active) - rexp(n_active)

  df_train$y <- drop(X_full %*% beta + rnorm(n, 0, sd))
  n_out <- ceiling(outlier_prop * n)
  if (n_out > 0) {
    idx <- sample.int(n, n_out)
    df_train$y[idx] <- df_train$y[idx] + rnorm(n_out, 0, 3 * sd)
  }

  XT <- gtools::rdirichlet(n_holdout, rep(1, p)); colnames(XT) <- colnames(X)
  df_hold <- data.frame(XT)
  terms_no_y <- delete.response(terms(formula_truth, data = df_train))
  X_T_full <- model.matrix(terms_no_y, model.frame(terms_no_y, df_hold))
  yT_true <- drop(X_T_full %*% beta)

  # Fits
  fit_svem <- SVEMnet(formula = formula_model, data = df_train,
                      standardize = TRUE, objective = "wAIC",
                      nBoot = 200, glmnet_alpha = c(0,0.5,1), weight_scheme = "SVEM")
  fit_id   <- SVEMnet(formula = formula_model, data = df_train,
                      standardize = TRUE, objective = "wAIC",
                      nBoot = 1, glmnet_alpha = c(0,0.5,1), weight_scheme = "Identity")

  # cv.glmnet wrapper (assumed available in the package)
  fit_cv <- glmnet_with_cv(formula = formula_model, data = df_train,
                           glmnet_alpha = c(0,0.5,1), standardize = TRUE)

  rows <- list(); k <- 1L
  # SVEM debias=FALSE/TRUE
  for (deb in c(FALSE, TRUE)) {
    yhat <- predict(fit_svem, newdata = df_hold, debias = deb)
    rows[[k]] <- data.frame(sim_id, p, d, n, sd, outlier_prop, model = "SVEMnet", debias = deb,
                            logRMSE = log(sqrt(mean((yT_true - yhat)^2)))); k <- k + 1L
  }
  # cv.glmnet debias=FALSE/TRUE
  for (deb in c(FALSE, TRUE)) {
    yhat <- predict_cv(fit_cv, newdata = df_hold, debias = deb)
    rows[[k]] <- data.frame(sim_id, p, d, n, sd, outlier_prop, model = "glmnet_cv", debias = deb,
                            logRMSE = log(sqrt(mean((yT_true - yhat)^2)))); k <- k + 1L
  }
  # Identity debias=FALSE/TRUE
  for (deb in c(FALSE, TRUE)) {
    yhat <- predict(fit_id, newdata = df_hold, debias = deb)
    rows[[k]] <- data.frame(sim_id, p, d, n, sd, outlier_prop, model = "SVEMnet_Identity", debias = deb,
                            logRMSE = log(sqrt(mean((yT_true - yhat)^2)))); k <- k + 1L
  }
  dplyr::bind_rows(rows)
}

results_list <- parLapplyLB(cl, sim_params_list, simulate_one)
stopCluster(cl)
results <- dplyr::bind_rows(results_list)

results <- results %>% group_by(sim_id) %>%
  mutate(sim_mean_logRMSE = mean(logRMSE), residual_logRMSE = logRMSE - sim_mean_logRMSE) %>%
  ungroup()

results$model_debias <- paste0(results$model, "_debias=", results$debias)
results$model_debias <- factor(results$model_debias,
                               levels = c("SVEMnet_debias=FALSE","SVEMnet_debias=TRUE",
                                          "glmnet_cv_debias=FALSE","glmnet_cv_debias=TRUE",
                                          "SVEMnet_Identity_debias=FALSE","SVEMnet_Identity_debias=TRUE"))

print(results %>% group_by(model, debias) %>%
        summarise(mean_residual_logRMSE = mean(residual_logRMSE), .groups = "drop"))

ggplot(results, aes(y = model_debias, x = logRMSE, fill = as.factor(debias))) +
  geom_boxplot() + ylab("Model & Debias Setting") + xlab("LogRMSE") +
  ggtitle("LogRMSE by Model and Debias Setting") + theme_minimal() +
  coord_flip() + scale_fill_discrete(name = "Debias", labels = c("FALSE","TRUE"))
