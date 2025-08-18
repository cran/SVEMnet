# ===== svemnet_sim_identity.R  =====
library(SVEMnet)
library(lhs)
library(ggplot2)
library(dplyr)
library(pbapply)
library(parallel)
library(tidyr)

p_values  <- seq(3, 6, 1)
d_values  <- seq(0.1, 0.9, 0.1)
n_values  <- seq(15, 50, 5)
sd_values <- c(0.25, 0.5, 1, 1.5)

nSim      <- 20
n_holdout <- 1000

param_grid <- expand.grid(p = p_values, d = d_values, n = n_values, sd = sd_values)

sim_params_list <- list(); sim_id <- 1L
for (i in seq_len(nrow(param_grid))) {
  for (sim in seq_len(nSim)) {
    sim_params_list[[length(sim_params_list) + 1L]] <- cbind(param_grid[i, ], sim_id = sim_id)
    sim_id <- sim_id + 1L
  }
}
cat("Total simulations to run:", length(sim_params_list), "\n")

num_cores <- max(1, parallel::detectCores(logical = FALSE) - 2)
RNGkind("L'Ecuyer-CMRG"); set.seed(0)
cl <- makeCluster(num_cores)
clusterEvalQ(cl, { library(SVEMnet); library(lhs); library(dplyr) })
clusterExport(cl, varlist = c("n_holdout"))
clusterSetRNGStream(cl, 0)

simulate_one <- function(sim_params) {
  p  <- sim_params$p; d <- sim_params$d; n <- sim_params$n; sd <- sim_params$sd; sim_id <- sim_params$sim_id
  X  <- lhs::randomLHS(n, p); colnames(X) <- paste0("V", seq_len(p))
  XT <- lhs::randomLHS(n_holdout, p); colnames(XT) <- colnames(X)

  n_active <- max(1, floor(p * d))
  active_params <- sample.int(p, n_active)
  beta <- numeric(p); beta[active_params] <- rexp(n_active) - rexp(n_active)

  y <- as.numeric(X %*% beta + rnorm(n, 0, sd))
  yT_true <- as.numeric(XT %*% beta)

  data_train   <- data.frame(y = y, X)
  data_holdout <- data.frame(XT)

  formula <- as.formula(paste0(
    "y ~ (", paste(colnames(X), collapse = " + "), ")^2 + ",
    paste("I(", colnames(X), "^2)", collapse = " + ")
  ))

  out <- list(); k <- 1L
  for (model_type in c("SVEM","FWR","Identity")) {
    if (model_type == "SVEM")     { weight_scheme <- "SVEM";     nBoot <- 200 }
    if (model_type == "FWR")      { weight_scheme <- "FWR";      nBoot <- 200 }
    if (model_type == "Identity") { weight_scheme <- "Identity"; nBoot <- 1   }

    fit <- SVEMnet(
      formula = formula, data = data_train,
      objective = "wAIC", nBoot = nBoot,
      glmnet_alpha = 1, weight_scheme = weight_scheme,
      standardize = TRUE
    )
    yhat <- predict(fit, newdata = data_holdout, debias = FALSE)
    rmse <- sqrt(mean((yT_true - yhat)^2))

    out[[k]] <- data.frame(sim_id = sim_id, p, d, n, sd, model_type, logRMSE = log(rmse))
    k <- k + 1L
  }
  dplyr::bind_rows(out)
}

results_list <- pblapply(sim_params_list, simulate_one, cl = cl)
stopCluster(cl)
results <- dplyr::bind_rows(results_list)
results$sim_id <- as.integer(results$sim_id)

results <- results %>% group_by(sim_id) %>%
  mutate(sim_mean_logRMSE = mean(logRMSE),
         residual_logRMSE = logRMSE - sim_mean_logRMSE) %>% ungroup()

print(results %>% group_by(model_type) %>%
        summarise(mean_residual_logRMSE = mean(residual_logRMSE), .groups = "drop"))

results$model_type <- factor(results$model_type, levels = c("SVEM","FWR","Identity"))

ggplot(results, aes(x = model_type, y = logRMSE)) + geom_boxplot() +
  xlab("Model Type") + ylab("LogRMSE") + ggtitle("LogRMSE by Model") + theme_minimal()

ggplot(results, aes(x = model_type, y = residual_logRMSE)) + geom_boxplot() +
  xlab("Model Type") + ylab("Residual LogRMSE") +
  ggtitle("Residuals (after removing simulation mean logRMSE)") + theme_minimal()

# Paired diffs
wide <- results %>% select(sim_id, model_type, logRMSE) %>% tidyr::pivot_wider(names_from = model_type, values_from = logRMSE)
wide <- wide %>% mutate(
  diff_SVEM_FWR      = SVEM - FWR,
  diff_SVEM_Identity = SVEM - Identity,
  diff_FWR_Identity  = FWR - Identity
)
print(wide %>% summarise(
  mean_diff_SVEM_FWR = mean(diff_SVEM_FWR, na.rm=TRUE),
  mean_diff_SVEM_Identity = mean(diff_SVEM_Identity, na.rm=TRUE),
  mean_diff_FWR_Identity  = mean(diff_FWR_Identity, na.rm=TRUE)
))
