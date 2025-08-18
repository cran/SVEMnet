## ============================================================
## SVEM (wAIC / wAICc / wSSE, no debias) vs glmnet_with_cv (defaults)
## Focus: small training sizes; compare test RMSE head-to-head
## Parallelized with future.apply (5 workers)
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
 library(progressr)  # uncomment if you want progress bars

## ----- Controls -----
set.seed(425)
REPS         <- 20
NBOOT        <- 200
R2_GRID      <- c(0.30, 0.50, 0.70, 0.90)
RHO_GRID     <- c(0, -0.5, 0.5, -0.9, 0.9)  # correlation among factors
P_GRID       <- c(3, 5, 7)
MODELS       <- c("main_effects","main_plus_int", "full_quadratic")
NTEST        <- 1000
DENSITY_GRID <- c(.05, .10, .20, .30, .50,.75)  # fraction of non-intercept terms active
N_WORKERS    <- 7

## ----- Helpers -----
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

safe_predict_svem <- function(fit, newdata, agg = "mean") {
  out <- try(predict(fit, newdata = newdata, debias = FALSE, agg = agg), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, nrow(newdata)))
  as.numeric(out)
}
safe_predict_cv <- function(fitcv, newdata) {
  out <- try(predict_cv(fitcv, newdata = newdata, debias = FALSE), silent = TRUE)
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

## Convenience wrapper to fit SVEM for a given objective (namespaced; robust)
fit_and_predict_svem <- function(fm, train_used, test_df, objective, nBoot = NBOOT, agg = "mean") {
  t0 <- proc.time()[3]
  mod <- try(SVEMnet::SVEMnet(
    formula = fm, data = train_used,
    nBoot = nBoot,
    glmnet_alpha = c(0,.5,1),
    weight_scheme = "SVEM",
    objective = objective,
    standardize = TRUE
  ), silent = TRUE)
  elapsed <- round(proc.time()[3] - t0, 2)
  if (inherits(mod, "try-error")) {
    preds <- rep(NA_real_, nrow(test_df))
  } else {
    preds <- safe_predict_svem(mod, test_df, agg = agg)
  }
  list(pred = preds, time = elapsed)
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

## ----- Parallel plan (5 workers) -----
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = N_WORKERS)

## Optional progress bar:
# handlers(global = TRUE)
# with_progress({
#   pbar <- progressor(along = seq_len(nrow(grid)))
#   ... call future_lapply with on.exit(pbar()) ...
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
  n_grid <- sort(unique(pmin(pmax(p_full + c(-12,-8,-6,-4,-2,0,2,4,6,8,12), 12), 40)))

  # per-row deterministic seed (so sampling is reproducible too)
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

  # SVEM variants (alpha=1; no debias)
  sv_waic  <- fit_and_predict_svem(fm, train_used, test_df, "wAIC",  nBoot = NBOOT)
  sv_waicc <- fit_and_predict_svem(fm, train_used, test_df, "wAICc", nBoot = NBOOT)
  sv_wsse  <- fit_and_predict_svem(fm, train_used, test_df, "wSSE",  nBoot = NBOOT)

  # glmnet_with_cv (alpha=1; no debias)
  t1 <- proc.time()[3]
  fit_cv <- glmnet_with_cv(
    formula = fm, data = train_used,
    glmnet_alpha = c(0,.5,1),
    standardize = TRUE
  )
  pr_te_cv <- safe_predict_cv(fit_cv, test_df)
  time_cv  <- round(proc.time()[3] - t1, 2)

  data.frame(
    p = p, model = model, rho = rho, R2_target = R2_tgt, density = dens,
    p_full = p_full, n_train = n_tr, n_train_used = n_used, n_test = n_te,
    n_train_minus_p_full = n_tr - p_full,
    ratio_n_over_p = n_used / p_full,
    above_boundary = as.integer(n_used > p_full),
    R2_true_train = R2_train_true, R2_true_test = R2_test_true,
    rmse_test_svem_waic  = rmse(test_df$y, sv_waic$pred),
    rmse_test_svem_waicc = rmse(test_df$y, sv_waicc$pred),
    rmse_test_svem_wsse  = rmse(test_df$y, sv_wsse$pred),
    rmse_test_glmnet_cv  = rmse(test_df$y, pr_te_cv),
    time_sec_svem_waic   = sv_waic$time,
    time_sec_svem_waicc  = sv_waicc$time,
    time_sec_svem_wsse   = sv_wsse$time,
    time_sec_cv          = time_cv,
    stringsAsFactors = FALSE
  )
}

## ----- Run in parallel -----
res_list <- future.apply::future_lapply(
  seq_len(nrow(grid)),
  function(i) {
    # if using progressr:
    # on.exit(pbar(), add = TRUE)
    run_one(grid[i, , drop = FALSE])
  },
  future.seed = TRUE
)

res <- dplyr::bind_rows(Filter(Negate(is.null), res_list))
stopifnot(nrow(res) > 0)

## ---------- Summary & win-rates (robust, single block) ----------
stopifnot(all(c("rmse_test_svem_waic","rmse_test_svem_waicc","rmse_test_svem_wsse","rmse_test_glmnet_cv") %in% names(res)))

# Long format for plotting and summaries
res_long <- res %>%
  mutate(
    Row           = dplyr::row_number(),
    lmrse_waic    = log(pmax(rmse_test_svem_waic,  .Machine$double.eps)),
    lmrse_waicc   = log(pmax(rmse_test_svem_waicc, .Machine$double.eps)),
    lmrse_wsse    = log(pmax(rmse_test_svem_wsse,  .Machine$double.eps)),
    lmrse_cv      = log(pmax(rmse_test_glmnet_cv,  .Machine$double.eps))
  ) %>%
  tidyr::pivot_longer(
    cols = c(lmrse_waic, lmrse_waicc, lmrse_wsse, lmrse_cv),
    names_to = "metric", values_to = "lmrse"
  ) %>%
  mutate(
    settings = dplyr::recode(metric,
                             lmrse_waic  = "SVEM_wAIC",
                             lmrse_waicc = "SVEM_wAICc",
                             lmrse_wsse  = "SVEM_wSSE",
                             lmrse_cv    = "glmnet_cv"),
    settings = factor(settings, levels = c("SVEM_wAIC","SVEM_wAICc","SVEM_wSSE","glmnet_cv"))
  )

# Head-to-head win-rates vs glmnet_cv (paired by row)
win_rates_vs_cv <- res %>%
  summarise(
    SVEM_wAIC  = mean(rmse_test_svem_waic  < rmse_test_glmnet_cv, na.rm = TRUE),
    SVEM_wAICc = mean(rmse_test_svem_waicc < rmse_test_glmnet_cv, na.rm = TRUE),
    SVEM_wSSE  = mean(rmse_test_svem_wsse  < rmse_test_glmnet_cv, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "settings", values_to = "winrate_vs_cv")

# Overall summary (log RMSE)
summ <- res_long %>%
  group_by(settings) %>%
  summarise(
    mean_lmrse = mean(lmrse, na.rm = TRUE),
    sd_lmrse   = sd(lmrse,   na.rm = TRUE),
    n          = dplyr::n(),
    se         = sd_lmrse / sqrt(pmax(n, 1)),
    ci_lo      = mean_lmrse - 1.96 * se,
    ci_hi      = mean_lmrse + 1.96 * se,
    .groups    = "drop"
  ) %>%
  left_join(win_rates_vs_cv, by = "settings")

cat("\n================ SUMMARY (log RMSE) ================\n")
print(summ, row.names = FALSE, digits = 4)

## ---------- By-boundary summaries (attach boundary to long) ----------
if ("above_boundary" %in% names(res)) {
  bound_key <- res %>%
    mutate(Row = dplyr::row_number()) %>%
    select(Row, above_boundary)
} else if (all(c("n_train_used","p_full") %in% names(res))) {
  bound_key <- res %>%
    mutate(Row = dplyr::row_number(),
           above_boundary = as.integer(n_train_used > p_full)) %>%
    select(Row, above_boundary)
} else if ("ratio_n_over_p" %in% names(res)) {
  bound_key <- res %>%
    mutate(Row = dplyr::row_number(),
           above_boundary = as.integer(ratio_n_over_p > 1)) %>%
    select(Row, above_boundary)
} else {
  bound_key <- res %>%
    mutate(Row = dplyr::row_number(), above_boundary = NA_integer_) %>%
    select(Row, above_boundary)
}

res_long <- res_long %>% left_join(bound_key, by = "Row")

if ("above_boundary.x" %in% names(res_long) || "above_boundary.y" %in% names(res_long)) {
  res_long <- res_long %>%
    dplyr::mutate(
      above_boundary = dplyr::coalesce(
        if ("above_boundary.y" %in% names(.)) .data$above_boundary.y else NA_integer_,
        if ("above_boundary.x" %in% names(.)) .data$above_boundary.x else NA_integer_
      )
    ) %>%
    dplyr::select(-dplyr::any_of(c("above_boundary.x","above_boundary.y")))
}

# Means of log RMSE by boundary & setting
summ_by_boundary <- res_long %>%
  group_by(above_boundary, settings) %>%
  summarise(
    mean_lmrse = mean(lmrse, na.rm = TRUE),
    sd_lmrse   = sd(lmrse,   na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )

cat("\n---- Means by boundary (n_train_used > p_full) ----\n")
print(summ_by_boundary, row.names = FALSE, digits = 4)

# Boundary-specific win-rate for each SVEM variant vs glmnet
win_by_boundary <- res %>%
  mutate(
    above_boundary = if ("above_boundary" %in% names(.)) above_boundary
    else if (all(c("n_train_used","p_full") %in% names(.))) as.integer(n_train_used > p_full)
    else if ("ratio_n_over_p" %in% names(.)) as.integer(ratio_n_over_p > 1)
    else NA_integer_,
    win_waic  = rmse_test_svem_waic  < rmse_test_glmnet_cv,
    win_waicc = rmse_test_svem_waicc < rmse_test_glmnet_cv,
    win_wsse  = rmse_test_svem_wsse  < rmse_test_glmnet_cv
  ) %>%
  group_by(above_boundary) %>%
  summarise(
    winrate_waic  = mean(win_waic,  na.rm = TRUE),
    winrate_waicc = mean(win_waicc, na.rm = TRUE),
    winrate_wsse  = mean(win_wsse,  na.rm = TRUE),
    n             = sum(!is.na(win_waic) | !is.na(win_waicc) | !is.na(win_wsse)),
    .groups       = "drop"
  )

cat("\n---- Win-rate vs glmnet (by boundary) ----\n")
print(win_by_boundary, row.names = FALSE, digits = 4)

## ========== ANOM-style, row-blocked plot ==========
res_long_anom <- res_long %>%
  dplyr::mutate(Row = as.factor(Row)) %>%
  dplyr::filter(is.finite(lmrse))

fit_anom <- lm(lmrse ~ settings + Row, data = res_long_anom)
aov_tbl  <- anova(fit_anom)

MSE    <- aov_tbl[["Mean Sq"]][nrow(aov_tbl)]
df_res <- aov_tbl[["Df"]][nrow(aov_tbl)]

t_methods <- nlevels(res_long_anom$settings)
b_blocks  <- nlevels(res_long_anom$Row)
grand_mu  <- mean(res_long_anom$lmrse, na.rm = TRUE)

# Var(ȳ_i. − ȳ..) = σ^2 * (t − 1) / (t * b) under RBD
se_group <- sqrt(MSE * (t_methods - 1) / (t_methods * b_blocks))

alpha <- 0.05
crit  <- qt(1 - alpha / (2 * t_methods), df = df_res)
UCL   <- grand_mu + crit * se_group
LCL   <- grand_mu - crit * se_group

means_df <- res_long_anom %>%
  dplyr::group_by(settings) %>%
  dplyr::summarise(mean_lmrse = mean(lmrse, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(flag = dplyr::case_when(
    mean_lmrse > UCL ~ "Above UCL",
    mean_lmrse < LCL ~ "Below LCL",
    TRUE             ~ "Within Limits"
  ))

p_anom <- ggplot(means_df, aes(x = settings, y = mean_lmrse)) +
  geom_hline(yintercept = grand_mu, linetype = 2) +
  geom_hline(yintercept = UCL,      linetype = 3) +
  geom_hline(yintercept = LCL,      linetype = 3) +
  geom_segment(aes(xend = settings, y = grand_mu, yend = mean_lmrse), linewidth = 1) +
  geom_point(aes(color = flag), size = 3) +
  scale_color_manual(values = c("Within Limits" = "black",
                                "Above UCL"     = "red",
                                "Below LCL"     = "red")) +
  labs(
    title = "Blocked ANOM-style plot for lmrse (Row = block)",
    subtitle = sprintf("Grand mean = %.3f | Limits = [%.3f, %.3f] | df = %d",
                       grand_mu, LCL, UCL, df_res),
    x = NULL, y = "Mean lmrse", color = NULL
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p_pairs <- ggplot(res_long_anom, aes(x = settings, y = lmrse, group = Row)) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.6, position = position_dodge(width = 0.05)) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  labs(title = "Paired runs by Row", x = NULL, y = "lmrse") +
  theme_bw()

(p_anom / p_pairs)

## ---------- Geometric-mean RMSE ratio (SVEM wAIC vs CV) ----------
eps <- .Machine$double.eps
rmse_ratio <- with(res, pmax(rmse_test_svem_waic, eps) / pmax(rmse_test_glmnet_cv, eps))
logdiff    <- log(rmse_ratio)

gm_ratio <- exp(mean(logdiff, na.rm = TRUE))           # primary
q75      <- quantile(rmse_ratio, 0.75, na.rm = TRUE)   # tail
q90      <- quantile(rmse_ratio, 0.90, na.rm = TRUE)
winrate  <- mean(rmse_ratio < 1, na.rm = TRUE)

# simple bootstrap CI for GM ratio
set.seed(123)
B <- 2000
idx <- replicate(B, sample.int(length(logdiff), replace = TRUE))
gm_boot <- apply(idx, 2, function(ii) exp(mean(logdiff[ii], na.rm = TRUE)))
gm_ci   <- quantile(gm_boot, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf(
  "\nPaired geometric-mean RMSE ratio (SVEM/CV): %.3f  (95%% CI %.3f – %.3f)\n",
  gm_ratio, gm_ci[1], gm_ci[2]
))
cat(sprintf("Win rate (SVEM better): %.1f%%\n", 100*winrate))
cat(sprintf("Tail ratios: 75th=%.3f, 90th=%.3f\n", q75, q90))

# Trim top/bottom 5% of logdiff before averaging
trim <- function(x, p=0.05) x[x >= quantile(x, p, na.rm=TRUE) & x <= quantile(x, 1-p, na.rm=TRUE)]
gm_ratio_trim <- exp(mean(trim(logdiff), na.rm = TRUE))
cat(sprintf("Trimmed GM ratio (5%%/5%%): %.3f\n", gm_ratio_trim))
