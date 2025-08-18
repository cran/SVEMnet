## ============================================================
## SVEM (wAICc & wAIC): debias=FALSE vs debias=TRUE
## Focus: small training sizes; head-to-head test RMSE
## Parallelized with future.apply (N_WORKERS)
## Includes paired analyses + interaction (objective × debias)
## ============================================================

## ----- Packages -----
if (!requireNamespace("SVEMnet", quietly = TRUE)) install.packages("SVEMnet")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
# Optional progress bar:
# if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr")

library(SVEMnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(future.apply)
# library(progressr)

## ----- Controls -----
set.seed(42)
REPS         <- 10
NBOOT        <- 200
R2_GRID      <- seq(.1,.9,.2)
RHO_GRID     <- c(0, -0.9, 0.9)
P_GRID       <- c(4,5,6,7)
MODELS       <- c( "full_quadratic", "full_cubic")  # added cubic
NTEST        <- 1000
DENSITY_GRID <- seq(.1,.9,.2)
N_WORKERS    <- 7

## ----- Helpers -----
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

safe_predict_svem <- function(fit, newdata, debias = FALSE, agg = "mean") {
  out <- try(predict(fit, newdata = newdata, debias = debias, agg = agg), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, nrow(newdata)))
  as.numeric(out)
}

## Count p_full (columns in model matrix – intercept)
count_p_full <- function(fm, p) {
  rhs <- stats::reformulate(attr(stats::delete.response(stats::terms(fm)), "term.labels"))
  tmp <- as.data.frame(matrix(rnorm(4 * p), nrow = 4))
  names(tmp) <- paste0("X", seq_len(p))
  mm <- model.matrix(rhs, data = tmp)
  sum(colnames(mm) != "(Intercept)")
}

build_formulas <- function(p) {
  vars <- paste0("X", seq_len(p))
  main_terms <- paste(vars, collapse = " + ")
  int2_terms <- paste0("(", main_terms, ")^2")  # up to 2-way interactions
  int3_terms <- paste0("(", main_terms, ")^3")  # up to 3-way interactions
  sq_terms   <- paste(sprintf("I(%s^2)", vars), collapse = " + ")
  cube_terms <- paste(sprintf("I(%s^3)", vars), collapse = " + ")
  list(
    main_effects   = as.formula(paste0("y ~ ", main_terms)),
    main_plus_int  = as.formula(paste0("y ~ ", int2_terms)),
    full_quadratic = as.formula(paste0("y ~ ", int2_terms, " + ", sq_terms)),
    full_cubic     = as.formula(paste0("y ~ ", int3_terms, " + ", sq_terms, " + ", cube_terms))
  )
}

gen_X <- function(n, p, rho = 0.5) {
  Z <- matrix(rnorm(n * p), n, p)
  if (p >= 2 && abs(rho) > 0) {
    for (j in 2:p) {
      Z[, j] <- rho * scale(Z[, j - 1], TRUE, TRUE) +
        sqrt(1 - rho^2) * scale(Z[, j], TRUE, TRUE)
    }
  }
  X <- Z
  colnames(X) <- paste0("X", seq_len(p))
  as.data.frame(X)
}

## Density-driven truth generator using the selected formula fm
make_sparse_data_R2 <- function(n, p, rho, target_R2, fm, density = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xdf <- gen_X(n, p, rho); names(Xdf) <- paste0("X", seq_len(p))
  terms_obj <- stats::delete.response(stats::terms(fm, data = transform(Xdf, y = 0)))
  mf  <- model.frame(terms_obj, data = Xdf)
  MM  <- model.matrix(terms_obj, mf)   # n x (M+1) with intercept at column 1
  
  M <- ncol(MM) - 1L
  if (M < 1L) stop("Formula produces no predictors beyond intercept.")
  n_active   <- max(1L, floor(density * M))
  active_idx <- sample.int(M, n_active, replace = FALSE) + 1L  # skip intercept
  
  beta <- numeric(M + 1L)
  beta[active_idx] <- rexp(n_active) - rexp(n_active)
  y_signal <- drop(MM %*% beta)
  
  sd_sig <- sqrt(max(var(y_signal), .Machine$double.eps))
  sd_eps <- sd_sig * sqrt((1 - target_R2) / target_R2)
  y <- y_signal + rnorm(n, 0, sd_eps)
  
  out <- cbind.data.frame(y = y, Xdf)
  out$y_signal <- y_signal
  rownames(out) <- sprintf("row%06d", seq_len(n))
  out
}

## Generic: fit SVEM once for a given objective, then predict for debias=FALSE/TRUE
fit_and_predict_obj <- function(fm, train_used, test_df, objective, nBoot = NBOOT, agg = "mean") {
  t0 <- proc.time()[3]
  mod <- try(SVEMnet::SVEMnet(
    formula = fm, data = train_used,
    nBoot = nBoot,
    glmnet_alpha = c(0, .5, 1),
    weight_scheme = "SVEM",
    objective = objective,  # "wAICc" or "wAIC"
    standardize = TRUE
  ), silent = TRUE)
  elapsed <- round(proc.time()[3] - t0, 2)
  if (inherits(mod, "try-error")) {
    pred_false <- rep(NA_real_, nrow(test_df))
    pred_true  <- rep(NA_real_, nrow(test_df))
  } else {
    pred_false <- safe_predict_svem(mod, test_df, debias = FALSE, agg = agg)
    pred_true  <- safe_predict_svem(mod, test_df, debias = TRUE,  agg = agg)
  }
  list(pred_false = pred_false, pred_true = pred_true, time = elapsed)
}

## ----- Build scenario grid -----
grid <- expand.grid(
  p       = P_GRID,
  model   = MODELS,
  rho     = RHO_GRID,
  R2_tgt  = R2_GRID,
  dens    = DENSITY_GRID,
  rep     = seq_len(REPS),
  stringsAsFactors = FALSE
)

## ----- Parallel plan -----
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = N_WORKERS)

## ----- One scenario run -----
run_one <- function(row) {
  p      <- row$p
  model  <- row$model
  rho    <- row$rho
  R2_tgt <- row$R2_tgt
  dens   <- row$dens
  repi   <- row$rep
  
  fms    <- build_formulas(p)
  fm     <- fms[[model]]
  p_full <- count_p_full(fm, p)
  
  # n around the boundary; clamp to [12, 40]
  n_grid <- sort(unique(pmin(pmax(p_full + c(-12,-8,-6,-4,-2,0,2), 12), 40)))
  
  # per-row deterministic seed
  seed_i <- 3000 +
    19 * repi + 100 * p +
    1000 * match(model, MODELS) +
    round(1000 * R2_tgt) +
    round(1000 * rho) +
    round(1000 * dens)
  set.seed(seed_i)
  
  n_tr <- sample(n_grid, 1L)
  n_te <- NTEST
  
  df <- make_sparse_data_R2(
    n = n_tr + n_te, p = p, rho = rho, target_R2 = R2_tgt,
    fm = fm, density = dens, seed = seed_i + 1L
  )
  
  set.seed(seed_i + 2L)
  idx      <- sample(seq_len(nrow(df)), size = n_tr)
  train_df <- df[idx, ]
  test_df  <- df[-idx, ]
  
  keep <- complete.cases(model.frame(fm, data = train_df, na.action = stats::na.pass))
  if (sum(keep) < 2) return(NULL)
  train_used <- train_df[keep, , drop = FALSE]
  n_used     <- nrow(train_used)
  
  r2_true <- function(d) var(d$y_signal) / var(d$y)
  R2_train_true <- r2_true(train_used)
  R2_test_true  <- r2_true(test_df)
  
  # Fit both objectives once each; predict debias FALSE & TRUE
  sv_wAICc <- fit_and_predict_obj(fm, train_used, test_df, "wAICc", nBoot = NBOOT)
  sv_wAIC  <- fit_and_predict_obj(fm, train_used, test_df, "wAIC",  nBoot = NBOOT)
  
  data.frame(
    p = p, model = model, rho = rho, R2_target = R2_tgt, density = dens,
    p_full = p_full, n_train = n_tr, n_train_used = n_used, n_test = n_te,
    n_train_minus_p_full = n_tr - p_full,
    ratio_n_over_p = n_used / p_full,
    above_boundary = as.integer(n_used > p_full),
    R2_true_train = R2_train_true, R2_true_test = R2_test_true,
    rmse_test_wAICc_debiasFalse = rmse(test_df$y, sv_wAICc$pred_false),
    rmse_test_wAICc_debiasTrue  = rmse(test_df$y, sv_wAICc$pred_true),
    rmse_test_wAIC_debiasFalse  = rmse(test_df$y, sv_wAIC$pred_false),
    rmse_test_wAIC_debiasTrue   = rmse(test_df$y, sv_wAIC$pred_true),
    time_sec_fit_wAICc          = sv_wAICc$time,
    time_sec_fit_wAIC           = sv_wAIC$time,
    stringsAsFactors = FALSE
  )
}

## ----- Run in parallel -----
res_list <- future.apply::future_lapply(
  seq_len(nrow(grid)),
  function(i) run_one(grid[i, , drop = FALSE]),
  future.seed = TRUE
)

res <- dplyr::bind_rows(Filter(Negate(is.null), res_list))
stopifnot(nrow(res) > 0)

## ---------- Summaries (marginal + paired + interaction) ----------
need_cols <- c("rmse_test_wAICc_debiasFalse","rmse_test_wAICc_debiasTrue",
               "rmse_test_wAIC_debiasFalse","rmse_test_wAIC_debiasTrue")
stopifnot(all(need_cols %in% names(res)))

## Long format (for marginal summaries & ANOM)
res_long <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row, model, p, rho, R2_target, density,
    above_boundary, ratio_n_over_p,
    wAICc_FALSE = log(pmax(rmse_test_wAICc_debiasFalse, .Machine$double.eps)),
    wAICc_TRUE  = log(pmax(rmse_test_wAICc_debiasTrue,  .Machine$double.eps)),
    wAIC_FALSE  = log(pmax(rmse_test_wAIC_debiasFalse,  .Machine$double.eps)),
    wAIC_TRUE   = log(pmax(rmse_test_wAIC_debiasTrue,   .Machine$double.eps))
  ) %>%
  tidyr::pivot_longer(
    cols = c(wAICc_FALSE, wAICc_TRUE, wAIC_FALSE, wAIC_TRUE),
    names_to = "setting", values_to = "lmrse"
  ) %>%
  mutate(
    objective = factor(ifelse(grepl("^wAICc", setting), "wAICc", "wAIC"),
                       levels = c("wAICc","wAIC")),
    debias    = factor(ifelse(grepl("TRUE$",  setting), "TRUE", "FALSE"),
                       levels = c("FALSE","TRUE")),
    setting   = factor(paste0(objective, " (debias=", debias, ")"),
                       levels = c("wAICc (debias=FALSE)","wAICc (debias=TRUE)",
                                  "wAIC (debias=FALSE)", "wAIC (debias=TRUE)"))
  )

## ---------- Which combo is best overall? (per-row winner among 4)
winners <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row,
    `wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue
  ) %>%
  tidyr::pivot_longer(-Row, names_to = "setting", values_to = "rmse") %>%
  group_by(Row) %>%
  filter(is.finite(rmse), rmse == min(rmse, na.rm = TRUE)) %>%
  ungroup()

best_counts <- winners %>% count(setting, name = "wins") %>%
  mutate(share = wins / sum(wins))
cat("\n==== Which (objective, debias) wins most often? ====\n")
print(best_counts, row.names = FALSE, digits = 4)

## ---------- Marginal summary (log RMSE) — (not paired; info only)
summ_marg <- res_long %>%
  group_by(setting) %>%
  summarise(
    mean_lmrse = mean(lmrse, na.rm = TRUE),
    sd_lmrse   = sd(lmrse,   na.rm = TRUE),
    n          = dplyr::n(),
    se         = sd_lmrse / sqrt(pmax(n, 1)),
    ci_lo      = mean_lmrse - 1.96 * se,
    ci_hi      = mean_lmrse + 1.96 * se,
    .groups    = "drop"
  ) %>%
  arrange(mean_lmrse)
cat("\n================ SUMMARY (log RMSE; marginal) ================\n")
print(summ_marg, row.names = FALSE, digits = 4)

## ---------- Debias effect within each objective (paired)
paired_effect <- function(num, den, label) {
  eps <- .Machine$double.eps
  logdiff <- log(pmax(num, eps)) - log(pmax(den, eps))
  logdiff <- logdiff[is.finite(logdiff)]
  n <- length(logdiff)
  mu <- mean(logdiff); se <- sd(logdiff)/sqrt(n); ci <- mu + c(-1,1) * qt(0.975, n-1) * se
  gm <- exp(mu); gci <- exp(ci)
  wr <- mean(exp(logdiff) < 1) # <1 favors numerator
  tibble::tibble(contrast = label, n = n,
                 gm_ratio = gm, gm_lo = gci[1], gm_hi = gci[2], winrate = wr)
}

debias_wAICc <- paired_effect(res$rmse_test_wAICc_debiasTrue,  res$rmse_test_wAICc_debiasFalse,
                              "wAICc: TRUE / FALSE")
debias_wAIC  <- paired_effect(res$rmse_test_wAIC_debiasTrue,   res$rmse_test_wAIC_debiasFalse,
                              "wAIC: TRUE / FALSE")

cat("\n==== Debias effect within each objective (GM ratios; <1 favors debias=TRUE) ====\n")
print(bind_rows(debias_wAICc, debias_wAIC), row.names = FALSE, digits = 4)

## ---------- Objective effect at fixed debias (paired)
obj_FALSE <- paired_effect(res$rmse_test_wAICc_debiasFalse, res$rmse_test_wAIC_debiasFalse,
                           "Objective: wAICc / wAIC  (debias=FALSE)")
obj_TRUE  <- paired_effect(res$rmse_test_wAICc_debiasTrue,  res$rmse_test_wAIC_debiasTrue,
                           "Objective: wAICc / wAIC  (debias=TRUE)")

cat("\n==== Objective effect at each debias (GM ratios; <1 favors wAICc) ====\n")
print(bind_rows(obj_FALSE, obj_TRUE), row.names = FALSE, digits = 4)

## ---------- Interaction: (debias TRUE vs FALSE) differs by objective?
## Difference-in-differences on log scale:
##   (log TRUE - log FALSE)_wAICc  -  (log TRUE - log FALSE)_wAIC
logdiff_debias_wAICc <- with(res, log(pmax(rmse_test_wAICc_debiasTrue,  .Machine$double.eps)) -
                               log(pmax(rmse_test_wAICc_debiasFalse, .Machine$double.eps)))
logdiff_debias_wAIC  <- with(res, log(pmax(rmse_test_wAIC_debiasTrue,   .Machine$double.eps)) -
                               log(pmax(rmse_test_wAIC_debiasFalse,  .Machine$double.eps)))
did <- logdiff_debias_wAICc - logdiff_debias_wAIC
did <- did[is.finite(did)]
n_did <- length(did)
mu_did <- mean(did); se_did <- sd(did)/sqrt(n_did)
ci_did <- mu_did + c(-1,1) * qt(0.975, n_did - 1) * se_did
cat(sprintf("\n==== Interaction (difference-in-differences on log RMSE): mean=%.4f, 95%% CI [%.4f, %.4f]\n",
            mu_did, ci_did[1], ci_did[2]))
cat(sprintf("On ratio scale, exp(mean)=%.3f  [%.3f, %.3f]  (values <1 mean stronger debias benefit for wAICc)\n",
            exp(mu_did), exp(ci_did[1]), exp(ci_did[2])))

## ---------- Row-blocked model for interaction (ANOVA style)
res_long_anom <- res_long %>%
  filter(is.finite(lmrse)) %>%
  mutate(Row = as.factor(Row))
fit_interact <- lm(lmrse ~ objective * debias + Row, data = res_long_anom)
cat("\n==== Two-way (objective × debias) with Row blocking ====\n")
print(anova(fit_interact))

## ---------- By-boundary (marginal means)
ab_key <- res %>% mutate(Row = dplyr::row_number()) %>% select(Row, above_boundary)
res_long_ab <- res_long %>% select(-dplyr::any_of(c("above_boundary"))) %>%
  left_join(ab_key, by = "Row")

summ_by_boundary <- res_long_ab %>%
  group_by(above_boundary, setting) %>%
  summarise(
    mean_lmrse = mean(lmrse, na.rm = TRUE),
    sd_lmrse   = sd(lmrse,   na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )
cat("\n---- Means by boundary (n_train_used > p_full; marginal) ----\n")
print(summ_by_boundary, row.names = FALSE, digits = 4)

## ---------- Win-rate by boundary (which setting wins among 4; paired within row)
win_by_boundary4 <- res %>%
  mutate(Row = dplyr::row_number(),
         above_boundary = as.integer(n_train_used > p_full)) %>%
  transmute(
    Row, above_boundary,
    `wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue
  ) %>%
  tidyr::pivot_longer(-c(Row, above_boundary), names_to = "setting", values_to = "rmse") %>%
  group_by(Row) %>%
  filter(is.finite(rmse), rmse == min(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  count(above_boundary, setting, name = "wins") %>%
  group_by(above_boundary) %>%
  mutate(share = wins / sum(wins)) %>%
  ungroup()

cat("\n---- Win shares by boundary (among 4 settings) ----\n")
print(win_by_boundary4, row.names = FALSE, digits = 4)

## ---------- ANOM-style, row-blocked plot (now 4 settings)
means_df <- res_long_anom %>%
  group_by(setting) %>%
  summarise(mean_lmrse = mean(lmrse, na.rm = TRUE), .groups = "drop")

aov_tbl  <- anova(lm(lmrse ~ setting + Row, data = res_long_anom))
MSE      <- aov_tbl[["Mean Sq"]][nrow(aov_tbl)]
df_res   <- aov_tbl[["Df"]][nrow(aov_tbl)]
t_methods<- nlevels(res_long_anom$setting)
b_blocks <- nlevels(res_long_anom$Row)
grand_mu <- mean(res_long_anom$lmrse, na.rm = TRUE)
se_group <- sqrt(MSE * (t_methods - 1) / (t_methods * b_blocks))
crit     <- qt(1 - 0.05 / (2 * t_methods), df = df_res)
UCL      <- grand_mu + crit * se_group
LCL      <- grand_mu - crit * se_group

p_anom <- ggplot(means_df, aes(x = setting, y = mean_lmrse)) +
  geom_hline(yintercept = grand_mu, linetype = 2) +
  geom_hline(yintercept = UCL,      linetype = 3) +
  geom_hline(yintercept = LCL,      linetype = 3) +
  geom_segment(aes(xend = setting, y = grand_mu, yend = mean_lmrse), linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Blocked ANOM-style plot (4 settings)", x = NULL, y = "Mean log RMSE") +
  theme_bw()

p_pairs <- ggplot(res_long_anom, aes(x = setting, y = lmrse, group = Row)) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.5, position = position_dodge(width = 0.05)) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  labs(title = "Paired runs by Row", x = NULL, y = "log RMSE") +
  theme_bw()

(p_anom / p_pairs)

## ---------- Optional visuals you already use for stratification can be reused
## e.g., rebuild 'pair' objects per contrast if you want heatmaps by n/p × R2, etc.
