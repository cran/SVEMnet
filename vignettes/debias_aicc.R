## ============================================================
## SVEM (wAICc only): debias=FALSE vs debias=TRUE
## Focus: small training sizes; compare test RMSE head-to-head
## Parallelized with future.apply (N_WORKERS)
## Includes paired analyses (log-diff CI, Wilcoxon) while keeping marginal summary
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
# library(progressr)  # uncomment if you want progress bars

## ----- Controls -----
set.seed(42)
REPS         <- 10
NBOOT        <- 200
R2_GRID      <- c(0.30, 0.50, 0.70, 0.90)
RHO_GRID     <- c(0, -0.5, 0.5, -0.9, 0.9)
P_GRID       <- c(4, 5,6, 7)
MODELS       <- c("main_plus_int", "full_quadratic")
NTEST        <- 1000
DENSITY_GRID <- c(.05, .10, .20, .30, .50, .75)
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
  int_terms  <- paste0("(", main_terms, ")^2")
  sq_terms   <- paste(sprintf("I(%s^2)", vars), collapse = " + ")
  list(
    main_effects   = as.formula(paste0("y ~ ", main_terms)),
    main_plus_int  = as.formula(paste0("y ~ ", int_terms)),
    full_quadratic = as.formula(paste0("y ~ ", int_terms, " + ", sq_terms))
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

## Fit SVEM once (objective = wAICc), then predict with debias = FALSE/TRUE
fit_and_predict_wAICc <- function(fm, train_used, test_df, nBoot = NBOOT, agg = "mean") {
  t0 <- proc.time()[3]
  mod <- try(SVEMnet::SVEMnet(
    formula = fm, data = train_used,
    nBoot = nBoot,
    glmnet_alpha = c(0, .5, 1),
    weight_scheme = "SVEM",
    objective = "wAICc",
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

## Optional progress bar:
# progressr::handlers(global = TRUE)
# progressr::with_progress({
#   pbar <- progressr::progressor(along = seq_len(nrow(grid)))
#   res_list <- future.apply::future_lapply(
#     seq_len(nrow(grid)),
#     function(i) { on.exit(pbar(), add = TRUE); run_one(grid[i, , drop = FALSE]) },
#     future.seed = TRUE
#   )
# })

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
  
  # Fit once (wAICc), predict debias FALSE & TRUE
  sv_wAICc <- fit_and_predict_wAICc(fm, train_used, test_df, nBoot = NBOOT)
  
  data.frame(
    p = p, model = model, rho = rho, R2_target = R2_tgt, density = dens,
    p_full = p_full, n_train = n_tr, n_train_used = n_used, n_test = n_te,
    n_train_minus_p_full = n_tr - p_full,
    ratio_n_over_p = n_used / p_full,
    above_boundary = as.integer(n_used > p_full),
    R2_true_train = R2_train_true, R2_true_test = R2_test_true,
    rmse_test_wAICc_debiasFalse = rmse(test_df$y, sv_wAICc$pred_false),
    rmse_test_wAICc_debiasTrue  = rmse(test_df$y, sv_wAICc$pred_true),
    time_sec_fit_wAICc          = sv_wAICc$time,
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

## ---------- Summary & win-rates (marginal + paired where noted) ----------
stopifnot(all(c("rmse_test_wAICc_debiasFalse","rmse_test_wAICc_debiasTrue") %in% names(res)))

# Long format for plotting and *marginal* summaries (kept for information)
res_long <- res %>%
  mutate(
    Row         = dplyr::row_number(),
    lmrse_false = log(pmax(rmse_test_wAICc_debiasFalse, .Machine$double.eps)),
    lmrse_true  = log(pmax(rmse_test_wAICc_debiasTrue,  .Machine$double.eps))
  ) %>%
  tidyr::pivot_longer(
    cols = c(lmrse_false, lmrse_true),
    names_to = "metric", values_to = "lmrse"
  ) %>%
  mutate(
    setting = dplyr::recode(metric,
                            lmrse_false = "wAICc (debias=FALSE)",
                            lmrse_true  = "wAICc (debias=TRUE)"),
    setting = factor(setting, levels = c("wAICc (debias=FALSE)", "wAICc (debias=TRUE)"))
  )

# Head-to-head win-rate: debias TRUE better than FALSE? (paired)
winrate_true_vs_false <- with(res,
                              mean(rmse_test_wAICc_debiasTrue < rmse_test_wAICc_debiasFalse, na.rm = TRUE)
)

# Overall marginal summary (log RMSE) — not paired, kept for info
summ <- res_long %>%
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
  mutate(winrate_vs_other = ifelse(setting == "wAICc (debias=TRUE)", winrate_true_vs_false, 1 - winrate_true_vs_false))

cat("\n================ SUMMARY (log RMSE; marginal) ================\n")
print(summ, row.names = FALSE, digits = 4)

## ---------- By-boundary summaries (marginal; robust) ----------
# Build a single clean above_boundary key from res
if ("above_boundary" %in% names(res)) {
  ab_key <- res %>%
    mutate(Row = dplyr::row_number()) %>%
    select(Row, above_boundary)
} else if (all(c("n_train_used","p_full") %in% names(res))) {
  ab_key <- res %>%
    mutate(Row = dplyr::row_number(),
           above_boundary = as.integer(n_train_used > p_full)) %>%
    select(Row, above_boundary)
} else if ("ratio_n_over_p" %in% names(res)) {
  ab_key <- res %>%
    mutate(Row = dplyr::row_number(),
           above_boundary = as.integer(ratio_n_over_p > 1)) %>%
    select(Row, above_boundary)
} else {
  ab_key <- tibble::tibble(Row = seq_len(nrow(res)), above_boundary = NA_integer_)
}

# Ensure res_long has exactly one 'above_boundary' column after join
res_long <- res_long %>%
  select(-dplyr::any_of(c("above_boundary",
                          "above_boundary.x",
                          "above_boundary.y"))) %>%
  left_join(ab_key, by = "Row")

summ_by_boundary <- res_long %>%
  group_by(above_boundary, setting) %>%
  summarise(
    mean_lmrse = mean(lmrse, na.rm = TRUE),
    sd_lmrse   = sd(lmrse,   na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )

cat("\n---- Means by boundary (n_train_used > p_full; marginal) ----\n")
print(summ_by_boundary, row.names = FALSE, digits = 4)


# Boundary-specific win-rate (TRUE beats FALSE?) — paired
win_by_boundary <- res %>%
  mutate(
    above_boundary = as.integer(n_train_used > p_full),
    win_true = rmse_test_wAICc_debiasTrue < rmse_test_wAICc_debiasFalse
  ) %>%
  group_by(above_boundary) %>%
  summarise(
    winrate_true = mean(win_true, na.rm = TRUE),
    n            = sum(!is.na(win_true)),
    .groups      = "drop"
  )

cat("\n---- Win-rate (debias=TRUE better than FALSE) by boundary (paired) ----\n")
print(win_by_boundary, row.names = FALSE, digits = 4)

## ========== ANOM-style, row-blocked plot (paired linear model) ==========
res_long_anom <- res_long %>%
  dplyr::mutate(Row = as.factor(Row)) %>%
  dplyr::filter(is.finite(lmrse))

fit_anom <- lm(lmrse ~ setting + Row, data = res_long_anom)
aov_tbl  <- anova(fit_anom)

MSE    <- aov_tbl[["Mean Sq"]][nrow(aov_tbl)]
df_res <- aov_tbl[["Df"]][nrow(aov_tbl)]

t_methods <- nlevels(res_long_anom$setting)
b_blocks  <- nlevels(res_long_anom$Row)
grand_mu  <- mean(res_long_anom$lmrse, na.rm = TRUE)

# Var(ȳ_i. − ȳ..) = σ^2 * (t − 1) / (t * b) under RBD
se_group <- sqrt(MSE * (t_methods - 1) / (t_methods * b_blocks))

alpha <- 0.05
crit  <- qt(1 - alpha / (2 * t_methods), df = df_res)
UCL   <- grand_mu + crit * se_group
LCL   <- grand_mu - crit * se_group

means_df <- res_long_anom %>%
  dplyr::group_by(setting) %>%
  dplyr::summarise(mean_lmrse = mean(lmrse, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(flag = dplyr::case_when(
    mean_lmrse > UCL ~ "Above UCL",
    mean_lmrse < LCL ~ "Below UCL",
    TRUE             ~ "Within Limits"
  ))

p_anom <- ggplot(means_df, aes(x = setting, y = mean_lmrse)) +
  geom_hline(yintercept = grand_mu, linetype = 2) +
  geom_hline(yintercept = UCL,      linetype = 3) +
  geom_hline(yintercept = LCL,      linetype = 3) +
  geom_segment(aes(xend = setting, y = grand_mu, yend = mean_lmrse), linewidth = 1) +
  geom_point(aes(color = flag), size = 3) +
  scale_color_manual(values = c("Within Limits" = "black",
                                "Above UCL"     = "red",
                                "Below UCL"     = "red")) +
  labs(
    title = "Blocked ANOM-style plot for lmrse (Row = block)",
    x = NULL, y = "Mean lmrse", color = NULL
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p_pairs <- ggplot(res_long_anom, aes(x = setting, y = lmrse, group = Row)) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.6, position = position_dodge(width = 0.05)) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  labs(title = "Paired runs by Row", x = NULL, y = "lmrse") +
  theme_bw()

(p_anom / p_pairs)

## ---------- Geometric-mean RMSE ratio (TRUE / FALSE) with bootstrap (paired) ----------
eps <- .Machine$double.eps
rmse_ratio <- with(res, pmax(rmse_test_wAICc_debiasTrue,  eps) /
                     pmax(rmse_test_wAICc_debiasFalse, eps))
logdiff    <- log(rmse_ratio)

gm_ratio <- exp(mean(logdiff, na.rm = TRUE))           # primary paired effect size
q75      <- quantile(rmse_ratio, 0.75, na.rm = TRUE)   # tail
q90      <- quantile(rmse_ratio, 0.90, na.rm = TRUE)
winrate  <- mean(rmse_ratio < 1, na.rm = TRUE)         # TRUE better than FALSE

# simple bootstrap CI for GM ratio (paired)
set.seed(123)
B <- 2000
idx <- replicate(B, sample.int(length(logdiff), replace = TRUE))
gm_boot <- apply(idx, 2, function(ii) exp(mean(logdiff[ii], na.rm = TRUE)))
gm_ci   <- quantile(gm_boot, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf(
  "\nPaired geometric-mean RMSE ratio (debias=TRUE / debias=FALSE): %.3f  (95%% CI %.3f – %.3f)\n",
  gm_ratio, gm_ci[1], gm_ci[2]
))
cat(sprintf("Win rate (debias=TRUE better): %.1f%%\n", 100*winrate))
cat(sprintf("Tail ratios: 75th=%.3f, 90th=%.3f\n", q75, q90))

# Trim top/bottom 5% of logdiff before averaging
trim <- function(x, p=0.05) x[x >= quantile(x, p, na.rm=TRUE) & x <= quantile(x, 1-p, na.rm=TRUE)]
gm_ratio_trim <- exp(mean(trim(logdiff), na.rm = TRUE))
cat(sprintf("Trimmed GM ratio (5%%/5%%): %.3f\n", gm_ratio_trim))

## ---------- NEW: Paired analytic CI + Wilcoxon + paired summary table ----------
logdiff <- with(res, log(pmax(rmse_test_wAICc_debiasTrue,  .Machine$double.eps)) -
                  log(pmax(rmse_test_wAICc_debiasFalse, .Machine$double.eps)))
logdiff <- logdiff[is.finite(logdiff)]

n <- length(logdiff)
if (n > 1) {
  mu <- mean(logdiff)
  se <- sd(logdiff)/sqrt(n)
  ci <- mu + c(-1,1) * qt(0.975, df = n-1) * se
  cat(sprintf("\nPaired log-RMSE diff (TRUE − FALSE): mean=%.4f, 95%% CI [%.4f, %.4f]\n",
              mu, ci[1], ci[2]))
  cat(sprintf("On the RMSE ratio scale: %.3f  [%.3f, %.3f]\n",
              exp(mu), exp(ci[1]), exp(ci[2])))
  
  wil <- suppressWarnings(wilcox.test(logdiff, mu = 0, alternative = "two.sided", exact = FALSE))
  cat(sprintf("Wilcoxon signed-rank p=%.4g\n", wil$p.value))
  
  paired_summ <- data.frame(
    n = n,
    gm_ratio_TRUE_over_FALSE = exp(mu),
    ci_lo = exp(ci[1]),
    ci_hi = exp(ci[2]),
    winrate_TRUE_better = mean(exp(logdiff) < 1)
  )
  cat("\n==== PAIRED SUMMARY (primary) ====\n")
  print(paired_summ, row.names = FALSE, digits = 4)
} else {
  cat("\nNot enough paired observations for analytic CI / Wilcoxon.\n")
}

## ---------- (Optional) Visualize paired effects directly ----------
# df_diff <- res %>%
#   transmute(
#     Row = dplyr::row_number(),
#     log_ratio = log(pmax(rmse_test_wAICc_debiasTrue,  .Machine$double.eps)) -
#                 log(pmax(rmse_test_wAICc_debiasFalse, .Machine$double.eps)),
#     above_boundary = as.integer(n_train_used > p_full)
#   ) %>% filter(is.finite(log_ratio))
#
# ggplot(df_diff, aes(x = log_ratio)) +
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_histogram(bins = 40) +
#   labs(title = "Distribution of paired log RMSE differences (TRUE − FALSE)",
#        x = "log(RMSE_TRUE) − log(RMSE_FALSE)  (negative favors debias=TRUE)",
#        y = "Count") +
#   theme_bw()







eps <- .Machine$double.eps

# Row-level paired objects
pair <- res %>%
  transmute(
    Row = dplyr::row_number(),
    p, model, rho, R2_target, density,
    p_full, n_train, n_train_used, n_train_minus_p_full, ratio_n_over_p,
    rmse_false = pmax(rmse_test_wAICc_debiasFalse, eps),
    rmse_true  = pmax(rmse_test_wAICc_debiasTrue,  eps)
  ) %>%
  mutate(
    ratio    = rmse_true / rmse_false,
    logdiff  = log(ratio),
    win      = ratio < 1,
    rho_mag  = factor(dplyr::case_when(
      rho == 0       ~ "0",
      abs(rho) == .5 ~ "0.5",
      abs(rho) >= .9 ~ "0.9"
    ), levels = c("0","0.5","0.9")),
    rho_sign = factor(ifelse(rho < 0, "neg", ifelse(rho > 0, "pos", "zero")),
                      levels = c("neg","zero","pos")),
    R2_fac     = factor(R2_target, levels = sort(unique(R2_target))),
    p_fac      = factor(p,        levels = sort(unique(p))),
    model_fac  = factor(model,    levels = c("main_effects","main_plus_int","full_quadratic")),
    density_fac= factor(density,  levels = sort(unique(density))),
    # Bins for distance to boundary
    boundary_bin = cut(ratio_n_over_p,
                       breaks = c(0, 0.75, 1, 1.25, 1.5, 2, Inf),
                       labels = c("<0.75","0.75–1","1–1.25","1.25–1.5","1.5–2",">2"),
                       include.lowest = TRUE),
    boundary_decile = dplyr::ntile(ratio_n_over_p, 10L)
  ) %>%
  dplyr::filter(is.finite(logdiff))

# Helper: paired summary within a group
paired_summary <- function(df) {
  n   <- nrow(df)
  mu  <- mean(df$logdiff)
  se  <- sd(df$logdiff)/sqrt(pmax(n, 1))
  ci  <- mu + c(-1,1) * qt(0.975, df = pmax(n - 1, 1)) * se
  wins <- sum(df$win, na.rm = TRUE)
  p_sign <- 2 * pbinom(q = min(wins, n - wins), size = n, prob = 0.5)
  tibble::tibble(
    n = n,
    gm_ratio = exp(mu),
    gm_lo = exp(ci[1]),
    gm_hi = exp(ci[2]),
    winrate = mean(df$win),
    sign_test_p = p_sign
  )
}

## 1) Stratify by distance to boundary × true R^2 (paired GM ratio & winrate)
strat_boundary_R2 <- pair %>%
  group_by(boundary_bin, R2_fac) %>%
  group_modify(~paired_summary(.x)) %>%
  ungroup()

print(strat_boundary_R2, n = 100)

# Visualize: heatmaps for GM ratio and winrate
p_gm <- ggplot(strat_boundary_R2,
               aes(x = R2_fac, y = boundary_bin, fill = gm_ratio)) +
  geom_tile() + geom_text(aes(label = sprintf("%.3f", gm_ratio)), size = 3) +
  labs(title = "GM RMSE ratio (TRUE/FALSE) by boundary & R2_target",
       x = "R2_target", y = "n/p bin", fill = "GM ratio (<1 favors TRUE)") +
  theme_bw()

p_wr <- ggplot(strat_boundary_R2,
               aes(x = R2_fac, y = boundary_bin, fill = winrate)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", winrate)), size = 3) +
  labs(title = "Win rate (TRUE better) by boundary & R2_target",
       x = "R2_target", y = "n/p bin", fill = "Win rate") +
  theme_bw()

p_gm / p_wr

## 2) Stratify by model × p × correlation magnitude (paired)
strat_model_p_rho <- pair %>%
  group_by(model_fac, p_fac, rho_mag) %>%
  group_modify(~paired_summary(.x)) %>%
  ungroup() %>%
  arrange(model_fac, p_fac, rho_mag)

print(strat_model_p_rho, n = 100)

ggplot(strat_model_p_rho,
       aes(x = p_fac, y = gm_ratio, group = rho_mag, shape = rho_mag)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line() +
  facet_wrap(~ model_fac) +
  labs(title = "GM RMSE ratio by model × p × |rho|",
       x = "p", y = "GM ratio (TRUE/FALSE)") +
  theme_bw()

## 3) Stratify by sparsity (density) × correlation sign (paired)
strat_density_rhosign <- pair %>%
  group_by(density_fac, rho_sign) %>%
  group_modify(~paired_summary(.x)) %>%
  ungroup()

print(strat_density_rhosign, n = 100)

ggplot(strat_density_rhosign,
       aes(x = density_fac, y = gm_ratio, group = rho_sign, linetype = rho_sign)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line() +
  labs(title = "GM RMSE ratio by sparsity × correlation sign",
       x = "Density", y = "GM ratio (TRUE/FALSE)", linetype = "rho sign") +
  theme_bw()

## 4) Smooth effect of n/p on paired log-diff, faceted by model (paired)
ggplot(pair, aes(x = ratio_n_over_p, y = logdiff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  facet_wrap(~ model_fac, scales = "free_x") +
  labs(title = "Paired log RMSE diff vs n/p (negative favors debias=TRUE)",
       x = "n_train_used / p_full", y = "log(RMSE_TRUE/RMSE_FALSE)") +
  theme_bw()

## 5) Where does TRUE *lose*? Top/bottom 5% tails, describe them (paired)
tails <- pair %>%
  mutate(rank = dplyr::percent_rank(ratio)) %>%
  mutate(tail = dplyr::case_when(
    rank <= 0.05 ~ "best_TRUE",
    rank >= 0.95 ~ "worst_TRUE",
    TRUE ~ "middle"
  ))

tail_summary <- tails %>%
  filter(tail != "middle") %>%
  group_by(tail) %>%
  summarise(
    n = dplyr::n(),
    gm_ratio = exp(mean(logdiff)),
    winrate = mean(win),
    med_np = median(ratio_n_over_p),
    med_R2t = median(R2_target),
    .groups = "drop"
  )
print(tail_summary)

# Categorical skew inside tails
tails %>%
  filter(tail != "middle") %>%
  count(tail, model_fac) %>%
  tidyr::pivot_wider(names_from = tail, values_from = n, values_fill = 0) %>%
  print()

tails %>%
  filter(tail != "middle") %>%
  count(tail, density_fac) %>%
  tidyr::pivot_wider(names_from = tail, values_from = n, values_fill = 0) %>%
  print()

## 6) Light-weight models to quantify drivers (paired)
# (a) Probability debias=TRUE wins
m_win <- glm(win ~ scale(ratio_n_over_p) + model_fac + p_fac +
               R2_fac + density_fac + rho_sign + rho_mag,
             data = pair, family = binomial())
summary(m_win)

# ---- FIXED: build newdata with matching factor levels for predict() ----
L <- list(
  ratio_n_over_p = seq(0.3, 3.0, length.out = 50),
  
  # vary model across facets (keep factor levels identical to training data)
  model_fac  = factor(levels(pair$model_fac), levels = levels(pair$model_fac)),
  
  # hold other factors at typical/middle levels (still factors!)
  p_fac      = factor(levels(pair$p_fac)[ceiling(length(levels(pair$p_fac))/2)],
                      levels = levels(pair$p_fac)),
  R2_fac     = factor(levels(pair$R2_fac)[ceiling(length(levels(pair$R2_fac))/2)],
                      levels = levels(pair$R2_fac)),
  density_fac= factor(levels(pair$density_fac)[ceiling(length(levels(pair$density_fac))/2)],
                      levels = levels(pair$density_fac)),
  rho_sign   = factor("zero", levels = levels(pair$rho_sign)),
  rho_mag    = factor("0",    levels = levels(pair$rho_mag))
)

newdat <- tidyr::expand_grid(!!!L)
newdat$win_hat <- predict(m_win, newdata = newdat, type = "response")

ggplot(newdat, aes(x = ratio_n_over_p, y = win_hat)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_line() +
  facet_wrap(~ model_fac) +
  labs(title = "Estimated P(debias=TRUE beats FALSE) vs n/p (from GLM)",
       x = "n_train_used / p_full", y = "Win probability") +
  theme_bw()

# (b) Effect size model on log ratio
m_log <- lm(logdiff ~ scale(ratio_n_over_p) + model_fac + p_fac +
              R2_fac + density_fac + rho + I(rho^2),
            data = pair)
summary(m_log)

## 7) Fine-grained boundary deciles × R2 × model (paired; larger table)
strat_decile <- pair %>%
  group_by(boundary_decile, R2_fac, model_fac) %>%
  group_modify(~paired_summary(.x)) %>%
  ungroup()
print(strat_decile, n = 100)

ggplot(strat_decile,
       aes(x = boundary_decile, y = gm_ratio, group = R2_fac)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line() +
  facet_wrap(~ model_fac) +
  labs(title = "GM ratio across n/p deciles × R2 × model",
       x = "n/p decile (1=lowest)", y = "GM ratio (TRUE/FALSE)") +
  theme_bw()

## 8) Calibration vs train/test signal strength (paired)
# Bin by train/test "true" R^2 quantiles to see stability
calib <- pair %>%
  mutate(
    qR2_tr = cut(R2_target, breaks = quantile(R2_target, probs = seq(0,1,0.25)),
                 include.lowest = TRUE),
    qR2_te = cut(R2_target, breaks = quantile(R2_target, probs = seq(0,1,0.25)),
                 include.lowest = TRUE)  # same here; choose R2_true_test if preferred
  ) %>%
  group_by(qR2_tr, model_fac) %>%
  group_modify(~paired_summary(.x)) %>%
  ungroup()

print(calib, n = 50)

ggplot(calib, aes(x = qR2_tr, y = gm_ratio)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line(aes(group = model_fac)) +
  facet_wrap(~ model_fac) +
  labs(title = "GM ratio by train R2 quartile × model",
       x = "Train R2 quartile (by target_R2)", y = "GM ratio") +
  theme_bw()
