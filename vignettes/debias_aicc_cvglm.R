## ============================================================
## SVEM (wAICc & wAIC) vs glmnet:
## debias = FALSE vs TRUE for all methods
## Focus: small training sizes; head-to-head test RMSE
## Parallelized with future.apply (N_WORKERS)
## Includes paired analyses + interactions + extra slices
## ============================================================

## ----- Packages -----
if (!requireNamespace("SVEMnet", quietly = TRUE)) install.packages("SVEMnet")
if (!requireNamespace("glmnet",  quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("tibble",  quietly = TRUE)) install.packages("tibble")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
# Optional progress bar:
# if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr")

library(SVEMnet)
library(glmnet)     # glmnet_with_cv likely depends on this
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(future.apply)
# library(progressr)

## ----- Controls -----
set.seed(67)
REPS         <- 10
NBOOT        <- 200
R2_GRID      <- seq(.1, .9, .2)
RHO_GRID     <- c(0, -0.9, 0.9)
P_GRID       <- c(4, 5, 6, 7)
MODELS       <- c("full_quadratic", "full_cubic")  # quadratic + cubic truths
NTEST        <- 1000
DENSITY_GRID <- seq(.1, .9, .2)
N_WORKERS    <- 7

## ----- Helpers -----
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

safe_predict_svem <- function(fit, newdata, debias = FALSE, agg = "mean") {
  out <- try(predict(fit, newdata = newdata, debias = debias, agg = agg), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, nrow(newdata)))
  as.numeric(out)
}

## Robustly call your glmnet_with_cv() + predict_cv() (they must exist in the environment)
safe_predict_cv <- function(object, newdata, debias = FALSE) {
  out <- try(predict_cv(object, newdata = newdata, debias = debias), silent = TRUE)
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
  y <- y_signal + stats::rnorm(n, 0, sd_eps)

  out <- cbind.data.frame(y = y, Xdf)
  out$y_signal <- y_signal
  rownames(out) <- sprintf("row%06d", seq_len(n))
  out
}

## SVEM: fit once for a given objective, then predict for debias=FALSE/TRUE
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

## glmnet comparator using your glmnet_with_cv() + predict_cv()
fit_and_predict_glmnet_cv <- function(fm, train_used, test_df, alphas = c(0, 0.5, 1)) {
  if (!exists("glmnet_with_cv")) stop("glmnet_with_cv() not found in environment.")
  if (!exists("predict_cv"))    stop("predict_cv() not found in environment.")
  t0 <- proc.time()[3]
  mod <- try(glmnet_with_cv(formula = fm, data = train_used, glmnet_alpha = alphas), silent = TRUE)
  elapsed <- round(proc.time()[3] - t0, 2)
  if (inherits(mod, "try-error")) {
    return(list(pred_false = rep(NA_real_, nrow(test_df)),
                pred_true  = rep(NA_real_, nrow(test_df)),
                time = elapsed))
  }
  pred_false <- safe_predict_cv(mod, test_df, debias = FALSE)
  pred_true  <- safe_predict_cv(mod, test_df, debias = TRUE)
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
  n_grid <- sort(unique(pmin(pmax(p_full + c(-12, -8, -6, -4, -2, 0, 2), 12), 40)))

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

  # SVEM (wAICc & wAIC)
  sv_wAICc <- fit_and_predict_obj(fm, train_used, test_df, "wAICc", nBoot = NBOOT)
  sv_wAIC  <- fit_and_predict_obj(fm, train_used, test_df, "wAIC",  nBoot = NBOOT)

  # glmnet comparator (cv + built-in debias)
  glmn <- fit_and_predict_glmnet_cv(fm, train_used, test_df)

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
    rmse_test_glmnet_debiasFalse= rmse(test_df$y, glmn$pred_false),
    rmse_test_glmnet_debiasTrue = rmse(test_df$y, glmn$pred_true),

    time_sec_fit_wAICc          = sv_wAICc$time,
    time_sec_fit_wAIC           = sv_wAIC$time,
    time_sec_fit_glmnet         = glmn$time,
    stringsAsFactors = FALSE
  )
}

## ----- Run in parallel -----
# progressr::with_progress({
#   pbar <- progressr::progressor(along = seq_len(nrow(grid)))
#   res_list <- future.apply::future_lapply(
#     seq_len(nrow(grid)),
#     function(i) { on.exit(pbar(), add = TRUE); run_one(grid[i, , drop = FALSE]) },
#     future.seed = TRUE
#   )
# })
res_list <- future.apply::future_lapply(
  seq_len(nrow(grid)),
  function(i) run_one(grid[i, , drop = FALSE]),
  future.seed = TRUE
)

res <- dplyr::bind_rows(Filter(Negate(is.null), res_list))
stopifnot(nrow(res) > 0)

## ---------- Summaries (marginal + paired + interactions) ----------
need_cols <- c("rmse_test_wAICc_debiasFalse","rmse_test_wAICc_debiasTrue",
               "rmse_test_wAIC_debiasFalse","rmse_test_wAIC_debiasTrue",
               "rmse_test_glmnet_debiasFalse","rmse_test_glmnet_debiasTrue")
stopifnot(all(need_cols %in% names(res)))

## Long format (for marginal summaries & ANOM)
eps <- .Machine$double.eps
res_long <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row, model, p, rho, R2_target, density,
    above_boundary, ratio_n_over_p,
    wAICc_FALSE = log(pmax(rmse_test_wAICc_debiasFalse, eps)),
    wAICc_TRUE  = log(pmax(rmse_test_wAICc_debiasTrue,  eps)),
    wAIC_FALSE  = log(pmax(rmse_test_wAIC_debiasFalse,  eps)),
    wAIC_TRUE   = log(pmax(rmse_test_wAIC_debiasTrue,   eps)),
    glmnet_FALSE= log(pmax(rmse_test_glmnet_debiasFalse,eps)),
    glmnet_TRUE = log(pmax(rmse_test_glmnet_debiasTrue, eps))
  ) %>%
  tidyr::pivot_longer(
    cols = c(wAICc_FALSE, wAICc_TRUE, wAIC_FALSE, wAIC_TRUE, glmnet_FALSE, glmnet_TRUE),
    names_to = "setting", values_to = "lmrse"
  ) %>%
  mutate(
    method    = case_when(grepl("^glmnet", setting) ~ "glmnet",
                          grepl("^wAICc", setting) ~ "SVEM-wAICc",
                          TRUE                      ~ "SVEM-wAIC"),
    debias    = ifelse(grepl("TRUE$",  setting), "TRUE", "FALSE"),
    method    = factor(method, levels = c("SVEM-wAICc","SVEM-wAIC","glmnet")),
    debias    = factor(debias, levels = c("FALSE","TRUE")),
    setting   = factor(paste0(method, " (debias=", debias, ")"),
                       levels = c("SVEM-wAICc (debias=FALSE)","SVEM-wAICc (debias=TRUE)",
                                  "SVEM-wAIC (debias=FALSE)","SVEM-wAIC (debias=TRUE)",
                                  "glmnet (debias=FALSE)","glmnet (debias=TRUE)"))
  )

## ---------- Which combo is best overall? (per-row winner among 6)
winners <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row,
    `SVEM-wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `SVEM-wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `SVEM-wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `SVEM-wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue,
    `glmnet (debias=FALSE)`     = rmse_test_glmnet_debiasFalse,
    `glmnet (debias=TRUE)`      = rmse_test_glmnet_debiasTrue
  ) %>%
  tidyr::pivot_longer(-Row, names_to = "setting", values_to = "rmse") %>%
  group_by(Row) %>%
  filter(is.finite(rmse), rmse == min(rmse, na.rm = TRUE)) %>%
  ungroup()

best_counts <- winners %>% count(setting, name = "wins") %>%
  mutate(share = wins / sum(wins)) %>%
  arrange(desc(share))
cat("\n==== Which (method, debias) wins most often? ====\n")
print(best_counts, row.names = FALSE, digits = 4)

## ---------- Marginal summary (log RMSE) — info only
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

## ---------- Debias effect within each method (paired)
paired_effect <- function(num, den, label) {
  logdiff <- log(pmax(num, eps)) - log(pmax(den, eps))
  logdiff <- logdiff[is.finite(logdiff)]
  n <- length(logdiff)
  mu <- mean(logdiff); se <- sd(logdiff)/sqrt(n); ci <- mu + c(-1,1) * qt(0.975, n-1) * se
  tibble(contrast = label, n = n,
         gm_ratio = exp(mu), gm_lo = exp(ci[1]), gm_hi = exp(ci[2]),
         winrate = mean(exp(logdiff) < 1))
}

debias_wAICc <- paired_effect(res$rmse_test_wAICc_debiasTrue,  res$rmse_test_wAICc_debiasFalse, "SVEM-wAICc: TRUE / FALSE")
debias_wAIC  <- paired_effect(res$rmse_test_wAIC_debiasTrue,   res$rmse_test_wAIC_debiasFalse,  "SVEM-wAIC: TRUE / FALSE")
debias_glm   <- paired_effect(res$rmse_test_glmnet_debiasTrue, res$rmse_test_glmnet_debiasFalse,"glmnet: TRUE / FALSE")

cat("\n==== Debias effect within each method (GM ratios; <1 favors debias=TRUE) ====\n")
print(bind_rows(debias_wAICc, debias_wAIC, debias_glm), row.names = FALSE, digits = 4)

## ---------- Objective effect at fixed debias (SVEM only; paired)
obj_FALSE <- paired_effect(res$rmse_test_wAICc_debiasFalse, res$rmse_test_wAIC_debiasFalse,
                           "SVEM: wAICc / wAIC (debias=FALSE)")
obj_TRUE  <- paired_effect(res$rmse_test_wAICc_debiasTrue,  res$rmse_test_wAIC_debiasTrue,
                           "SVEM: wAICc / wAIC (debias=TRUE)")
cat("\n==== Objective effect (SVEM) at each debias (GM ratios; <1 favors wAICc) ====\n")
print(bind_rows(obj_FALSE, obj_TRUE), row.names = FALSE, digits = 4)

## ---------- SVEM vs glmnet at fixed debias (paired; optional but handy)
svem_wAICc_vs_glm_FALSE <- paired_effect(res$rmse_test_wAICc_debiasFalse, res$rmse_test_glmnet_debiasFalse,
                                         "SVEM-wAICc / glmnet  (debias=FALSE)")
svem_wAICc_vs_glm_TRUE  <- paired_effect(res$rmse_test_wAICc_debiasTrue,  res$rmse_test_glmnet_debiasTrue,
                                         "SVEM-wAICc / glmnet  (debias=TRUE)")
svem_wAIC_vs_glm_FALSE  <- paired_effect(res$rmse_test_wAIC_debiasFalse,  res$rmse_test_glmnet_debiasFalse,
                                         "SVEM-wAIC / glmnet   (debias=FALSE)")
svem_wAIC_vs_glm_TRUE   <- paired_effect(res$rmse_test_wAIC_debiasTrue,   res$rmse_test_glmnet_debiasTrue,
                                         "SVEM-wAIC / glmnet   (debias=TRUE)")
cat("\n==== SVEM vs glmnet (GM ratios; <1 favors SVEM) ====\n")
print(bind_rows(svem_wAICc_vs_glm_FALSE, svem_wAICc_vs_glm_TRUE,
                svem_wAIC_vs_glm_FALSE, svem_wAIC_vs_glm_TRUE),
      row.names = FALSE, digits = 4)

## ---------- Interaction: method × debias (Row-blocked + DiD style)
# assumes res_long is already built
res_long_dm <- res_long %>%
  dplyr::group_by(Row) %>%
  dplyr::mutate(lmrse_dm = lmrse - mean(lmrse)) %>%
  dplyr::ungroup()

fit_interact_fast <- lm(lmrse_dm ~ method * debias, data = res_long_dm)
anova(fit_interact_fast)

# If you want the “setting” one-way with Row FE:
fit_setting_fast <- lm(lmrse_dm ~ setting, data = res_long_dm)
anova(fit_setting_fast)

## ---------- By-boundary (marginal means over 6 settings)
ab_key <- res %>% mutate(Row = dplyr::row_number()) %>% select(Row, above_boundary)
res_long_ab <- res_long %>%
  select(-dplyr::any_of(c("above_boundary"))) %>%
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

## ---------- Win-rate by boundary (which setting wins among 6; paired within row)
win_by_boundary6 <- res %>%
  mutate(Row = dplyr::row_number(),
         above_boundary = as.integer(n_train_used > p_full)) %>%
  transmute(
    Row, above_boundary,
    `SVEM-wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `SVEM-wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `SVEM-wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `SVEM-wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue,
    `glmnet (debias=FALSE)`     = rmse_test_glmnet_debiasFalse,
    `glmnet (debias=TRUE)`      = rmse_test_glmnet_debiasTrue
  ) %>%
  tidyr::pivot_longer(-c(Row, above_boundary), names_to = "setting", values_to = "rmse") %>%
  group_by(Row) %>%
  filter(is.finite(rmse), rmse == min(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  count(above_boundary, setting, name = "wins") %>%
  group_by(above_boundary) %>%
  mutate(share = wins / sum(wins)) %>%
  ungroup()

cat("\n---- Win shares by boundary (among 6 settings) ----\n")
print(win_by_boundary6, row.names = FALSE, digits = 4)

## ---------- Extra slices: wins by density & by R2_target
win_by_density <- winners %>%
  left_join(res %>% mutate(Row = dplyr::row_number(), density_fac = factor(density)), by = "Row") %>%
  count(density_fac, setting, name = "wins") %>%
  group_by(density_fac) %>% mutate(share = wins / sum(wins)) %>% ungroup()
cat("\n---- Win shares by density ----\n")
print(win_by_density, row.names = FALSE, digits = 4)

win_by_R2 <- winners %>%
  left_join(res %>% mutate(Row = dplyr::row_number(), R2_fac = factor(R2_target)), by = "Row") %>%
  count(R2_fac, setting, name = "wins") %>%
  group_by(R2_fac) %>% mutate(share = wins / sum(wins)) %>% ungroup()
cat("\n---- Win shares by R2_target ----\n")
print(win_by_R2, row.names = FALSE, digits = 4)

## ---------- ANOM-style, row-blocked plot (now 6 settings)
means_df <- res_long_anom %>%
  group_by(setting) %>%
  summarise(mean_lmrse = mean(lmrse, na.rm = TRUE), .groups = "drop")

## ---- FAST ANOM pieces (avoid lm(... + Row)) ----
# 1) Row-demean ( already ran this earlier)
res_long_dm <- res_long %>%
  dplyr::group_by(Row) %>%
  dplyr::mutate(lmrse_dm = lmrse - mean(lmrse)) %>%
  dplyr::ungroup()

# 2) One-way on setting using the demeaned response
fit_setting_fast <- lm(lmrse_dm ~ setting, data = res_long_dm)
aov_tbl <- anova(fit_setting_fast)

# 3) Extract the same quantities you were using before
MSE      <- aov_tbl[["Mean Sq"]][nrow(aov_tbl)]
df_res   <- aov_tbl[["Df"]][nrow(aov_tbl)]
t_methods<- nlevels(res_long_anom$setting)
b_blocks <- nlevels(res_long_anom$Row)
grand_mu <- mean(res_long_anom$lmrse, na.rm = TRUE)

# ANOM group SE and limits (same formulas as before)
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
  labs(title = "Blocked ANOM-style plot (6 settings)", x = NULL, y = "Mean log RMSE") +
  theme_bw()

p_pairs <- ggplot(res_long_anom, aes(x = setting, y = lmrse, group = Row)) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.5, position = position_dodge(width = 0.05)) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  labs(title = "Paired runs by Row", x = NULL, y = "log RMSE") +
  theme_bw()

p_anom

## ---------- (Optional) Timing summary
time_summary <- res %>%
  summarise(
    median_wAICc = median(time_sec_fit_wAICc, na.rm = TRUE),
    median_wAIC  = median(time_sec_fit_wAIC,  na.rm = TRUE),
    median_glm   = median(time_sec_fit_glmnet,na.rm = TRUE)
  )
cat("\n---- Median fit times (sec) ----\n")
print(time_summary, row.names = FALSE, digits = 4)



## Win shares by n/p deciles (how “over/under” the boundary affects winners)
win_by_np <- winners %>%
  left_join(res %>% mutate(Row = dplyr::row_number(),
                           np_decile = ntile(ratio_n_over_p, 10)), by = "Row") %>%
  count(np_decile, setting, name = "wins") %>%
  group_by(np_decile) %>% mutate(share = wins / sum(wins)) %>% ungroup()

print(win_by_np, n = 60)

## Winners by model family (quadratic vs. cubic)
win_by_model <- winners %>%
  left_join(res %>% mutate(Row = dplyr::row_number(),
                           model = factor(model)), by = "Row") %>%
  count(model, setting, name = "wins") %>%
  group_by(model) %>% mutate(share = wins / sum(wins)) %>% ungroup()

print(win_by_model, n = 20)








paired_by <- function(res, num_col, den_col, by) {
  eps <- .Machine$double.eps
  res %>%
    mutate(Row = dplyr::row_number()) %>%
    select(Row, !!by, num = all_of(num_col), den = all_of(den_col)) %>%
    filter(is.finite(num), is.finite(den)) %>%
    group_by(!!sym(by)) %>%
    summarise(
      n = n(),
      gm_ratio = exp(mean(log(pmax(num, eps)) - log(pmax(den, eps)))),
      gm_lo = exp(t.test(log(pmax(num, eps)) - log(pmax(den, eps)))$conf.int[1]),
      gm_hi = exp(t.test(log(pmax(num, eps)) - log(pmax(den, eps)))$conf.int[2]),
      winrate = mean(num < den),
      .groups = "drop"
    )
}

# Prepare stratifiers
res_np <- res %>% mutate(np_decile = ntile(ratio_n_over_p, 10))
res_den <- res %>% mutate(density_fac = factor(density))
res_r2  <- res %>% mutate(R2_fac = factor(R2_target))

# SVEM-wAICc debias effect
db_wAICc_np <- paired_by(res_np,  "rmse_test_wAICc_debiasTrue",  "rmse_test_wAICc_debiasFalse", "np_decile")
db_wAICc_den<- paired_by(res_den, "rmse_test_wAICc_debiasTrue",  "rmse_test_wAICc_debiasFalse", "density_fac")
db_wAICc_r2 <- paired_by(res_r2,  "rmse_test_wAICc_debiasTrue",  "rmse_test_wAICc_debiasFalse", "R2_fac")

# glmnet debias effect
db_glm_np  <- paired_by(res_np,  "rmse_test_glmnet_debiasTrue", "rmse_test_glmnet_debiasFalse", "np_decile")
db_glm_den <- paired_by(res_den, "rmse_test_glmnet_debiasTrue", "rmse_test_glmnet_debiasFalse", "density_fac")
db_glm_r2  <- paired_by(res_r2,  "rmse_test_glmnet_debiasTrue", "rmse_test_glmnet_debiasFalse", "R2_fac")

# SVEM-wAIC debias effect
db_wAIC_np  <- paired_by(res_np,  "rmse_test_wAIC_debiasTrue",  "rmse_test_wAIC_debiasFalse",  "np_decile")
db_wAIC_den <- paired_by(res_den, "rmse_test_wAIC_debiasTrue",  "rmse_test_wAIC_debiasFalse",  "density_fac")
db_wAIC_r2  <- paired_by(res_r2,  "rmse_test_wAIC_debiasTrue",  "rmse_test_wAIC_debiasFalse",  "R2_fac")

db_wAICc_np; db_glm_np; db_wAIC_np
db_wAICc_den; db_glm_den; db_wAIC_den
db_wAICc_r2; db_glm_r2; db_wAIC_r2



svem_vs_glm_by <- function(res, by, debias_flag = c("False","True")) {
  eps <- .Machine$double.eps
  col_svem <- paste0("rmse_test_wAICc_debias", debias_flag)
  col_glm  <- paste0("rmse_test_glmnet_debias", debias_flag)
  res %>%
    mutate(Row = dplyr::row_number()) %>%
    select(Row, !!by, svem = all_of(col_svem), glm = all_of(col_glm)) %>%
    filter(is.finite(svem), is.finite(glm)) %>%
    group_by(!!sym(by)) %>%
    summarise(
      n = n(),
      gm_ratio = exp(mean(log(pmax(svem, eps)) - log(pmax(glm, eps)))),
      gm_lo = exp(t.test(log(pmax(svem, eps)) - log(pmax(glm, eps)))$conf.int[1]),
      gm_hi = exp(t.test(log(pmax(svem, eps)) - log(pmax(glm, eps)))$conf.int[2]),
      winrate = mean(svem < glm),
      .groups = "drop"
    ) %>%
    mutate(debias = debias_flag)
}

sv_glm_np_F <- svem_vs_glm_by(res_np,  "np_decile",  "False")
sv_glm_np_T <- svem_vs_glm_by(res_np,  "np_decile",  "True")
sv_glm_den_F<- svem_vs_glm_by(res_den, "density_fac","False")
sv_glm_den_T<- svem_vs_glm_by(res_den, "density_fac","True")
sv_glm_r2_F <- svem_vs_glm_by(res_r2,  "R2_fac",     "False")
sv_glm_r2_T <- svem_vs_glm_by(res_r2,  "R2_fac",     "True")

rbind(sv_glm_np_F, sv_glm_np_T) %>% arrange(np_decile, debias)
na_rates <- res %>%
  transmute(
    glmnet_FALSE = is.na(rmse_test_glmnet_debiasFalse),
    glmnet_TRUE  = is.na(rmse_test_glmnet_debiasTrue),
    wAICc_FALSE  = is.na(rmse_test_wAICc_debiasFalse),
    wAICc_TRUE   = is.na(rmse_test_wAICc_debiasTrue),
    wAIC_FALSE   = is.na(rmse_test_wAIC_debiasFalse),
    wAIC_TRUE    = is.na(rmse_test_wAIC_debiasTrue)
  ) %>%
  summarise(across(everything(), ~mean(.)))
print(na_rates, digits = 3)
frontier <- res_long %>%
  group_by(setting) %>%
  summarise(mean_lmrse = mean(lmrse, na.rm = TRUE), .groups="drop") %>%
  mutate(median_time = case_when(
    grepl("^SVEM-wAICc", setting) ~ median(res$time_sec_fit_wAICc, na.rm = TRUE),
    grepl("^SVEM-wAIC",  setting) ~ median(res$time_sec_fit_wAIC,  na.rm = TRUE),
    grepl("^glmnet",     setting) ~ median(res$time_sec_fit_glmnet,na.rm = TRUE)
  ))
frontier
# ggplot(frontier, aes(median_time, mean_lmrse, label = setting)) + geom_point() + ggrepel::geom_text_repel()





winner_margins <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row,
    `SVEM-wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `SVEM-wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `SVEM-wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `SVEM-wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue,
    `glmnet (debias=FALSE)`     = rmse_test_glmnet_debiasFalse,
    `glmnet (debias=TRUE)`      = rmse_test_glmnet_debiasTrue
  ) %>%
  tidyr::pivot_longer(-Row, names_to="setting", values_to="rmse") %>%
  group_by(Row) %>%
  summarise(
    best_setting = setting[which.min(rmse)],
    best_rmse = min(rmse, na.rm=TRUE),
    second_rmse = sort(rmse, partial=2)[2],
    margin_pct = (second_rmse / best_rmse) - 1,
    .groups="drop"
  )

winner_margins %>%
  summarise(
    median_margin_pct = median(margin_pct, na.rm=TRUE),
    p90_margin_pct    = quantile(margin_pct, 0.9, na.rm=TRUE)
  )


winner_margins <- res %>%
  mutate(Row = dplyr::row_number()) %>%
  transmute(
    Row,
    `SVEM-wAICc (debias=FALSE)` = rmse_test_wAICc_debiasFalse,
    `SVEM-wAICc (debias=TRUE)`  = rmse_test_wAICc_debiasTrue,
    `SVEM-wAIC (debias=FALSE)`  = rmse_test_wAIC_debiasFalse,
    `SVEM-wAIC (debias=TRUE)`   = rmse_test_wAIC_debiasTrue,
    `glmnet (debias=FALSE)`     = rmse_test_glmnet_debiasFalse,
    `glmnet (debias=TRUE)`      = rmse_test_glmnet_debiasTrue
  ) %>%
  tidyr::pivot_longer(-Row, names_to="setting", values_to="rmse") %>%
  group_by(Row) %>%
  summarise(
    best_setting = setting[which.min(rmse)],
    best_rmse = min(rmse, na.rm=TRUE),
    second_rmse = { vals <- rmse[is.finite(rmse)] ; sort(vals, partial=2)[2] },
    margin_pct = (second_rmse / best_rmse) - 1,
    .groups="drop"
  )

# Overall margins (adds “how often is it basically a tie?”)
winner_margins %>%
  summarise(
    median_margin_pct = median(margin_pct, na.rm=TRUE),
    p90_margin_pct    = quantile(margin_pct, 0.9, na.rm=TRUE),
    share_lt_2pct     = mean(margin_pct < 0.02, na.rm=TRUE),
    share_lt_5pct     = mean(margin_pct < 0.05, na.rm=TRUE)
  )

# Margins by winner (who wins “by a lot”?)
winner_margins %>%
  group_by(best_setting) %>%
  summarise(
    n = n(),
    med = median(margin_pct, na.rm=TRUE),
    p90 = quantile(margin_pct, 0.9, na.rm=TRUE),
    share_gt_5pct  = mean(margin_pct > 0.05, na.rm=TRUE),
    share_gt_10pct = mean(margin_pct > 0.10, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(med))

margins_np_r2 <- winner_margins %>%
  left_join(res %>% mutate(
    Row = dplyr::row_number(),
    np_decile = dplyr::ntile(ratio_n_over_p, 10),
    R2_fac = factor(R2_target)
  ), by="Row") %>%
  group_by(R2_fac, np_decile) %>%
  summarise(
    med_margin = median(margin_pct, na.rm=TRUE),
    p90_margin = quantile(margin_pct, 0.9, na.rm=TRUE),
    .groups="drop"
  )
margins_np_r2

tau <- 0.02  # try 0.02 (2%) or 0.05 (5%)

policy_eval <- winner_margins %>%
  mutate(choice = ifelse(best_setting == "SVEM-wAICc (debias=FALSE)" & margin_pct > tau,
                         "SVEM-wAICc (debias=FALSE)",
                         "glmnet (debias=FALSE)")) %>%
  left_join(
    res %>% mutate(Row = dplyr::row_number()) %>%
      transmute(Row,
                `SVEM-wAICc (debias=FALSE)`=rmse_test_wAICc_debiasFalse,
                `SVEM-wAICc (debias=TRUE)` =rmse_test_wAICc_debiasTrue,
                `SVEM-wAIC (debias=FALSE)` =rmse_test_wAIC_debiasFalse,
                `SVEM-wAIC (debias=TRUE)`  =rmse_test_wAIC_debiasTrue,
                `glmnet (debias=FALSE)`    =rmse_test_glmnet_debiasFalse,
                `glmnet (debias=TRUE)`     =rmse_test_glmnet_debiasTrue
      ) %>% tidyr::pivot_longer(-Row, names_to="setting", values_to="rmse"),
    by=c("Row","choice"="setting")
  ) %>%
  rename(rmse_choice = rmse) %>%
  mutate(regret_pct = (rmse_choice / best_rmse) - 1) %>%
  left_join(
    res %>% mutate(Row = dplyr::row_number()) %>%
      select(Row, time_sec_fit_wAICc, time_sec_fit_wAIC, time_sec_fit_glmnet),
    by="Row"
  ) %>%
  mutate(time_choice = dplyr::case_when(
    grepl("^glmnet", choice)      ~ time_sec_fit_glmnet,
    grepl("^SVEM-wAICc", choice)  ~ time_sec_fit_wAICc,
    grepl("^SVEM-wAIC", choice)   ~ time_sec_fit_wAIC
  ))

policy_eval %>%
  summarise(
    median_regret_pct = median(regret_pct, na.rm=TRUE),
    p90_regret_pct    = quantile(regret_pct, 0.9, na.rm=TRUE),
    share_no_regret   = mean(regret_pct <= 0, na.rm=TRUE),
    median_time_sec   = median(time_choice, na.rm=TRUE)
  )


