## ======================================================================
## SVEMnet boundary sims: wAIC vs wSSE + debias(TRUE/FALSE)
## - standardize = TRUE
## - glmnet_alpha = c(1)
## - No leakage; test is disjoint from train; debias refit check
## - n_train ranges chosen to cross n_train_used vs p_full boundary
## - Outputs:
##   * res ........ per-fit records (you can model this)
##   * A tidy summary with 4 win-rates:
##       winrate_debias_TRUE, winrate_debias_FALSE,
##       winrate_wAIC, winrate_wSSE
##   * Two simple plots (comment out prints if running headless)
## ======================================================================

if (!requireNamespace("SVEMnet", quietly = TRUE)) install.packages("SVEMnet")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(SVEMnet)
library(ggplot2)

## ---------- Logging / Progress ----------
VERBOSE     <- TRUE     # per-fit messages
ECHO_EVERY  <- 1        # print every k-th fit
PROGRESS    <- TRUE     # progress bar

stamp   <- function() format(Sys.time(), "%H:%M:%S")
log_msg <- function(...) if (VERBOSE) {
  cat(sprintf("[%s] %s\n", stamp(), sprintf(...))); flush.console()
}
options(warn = 1)       # surface warnings immediately

set.seed(20250814)

## ------------------ Controls ------------------
REPS      <- 10           # increase to 20+ for more power
NBOOT     <- 200
R2_GRID   <- c(0.30, 0.70, 0.90)
RHO_GRID  <- c(0.2, 0.6, 0.9)
P_GRID    <- c(4, 8, 12)
NTEST     <- 1000
MODELS    <- c( "main_plus_int", "full_quadratic")

## ------------------ Helpers -------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

safe_predict <- function(fit, newdata, debias=FALSE) {
  out <- try(predict(fit, newdata=newdata, debias=debias), silent=TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, nrow(newdata)))
  as.numeric(out)
}

## Count p_full (columns in model.matrix minus intercept) WITHOUT response 'y'
count_p_full <- function(fm, p) {
  rhs_fml <- stats::reformulate(
    attr(stats::delete.response(stats::terms(fm)), "term.labels")
  )
  tmp <- as.data.frame(matrix(rnorm(4 * p), nrow = 4))
  names(tmp) <- paste0("X", seq_len(p))
  mm <- model.matrix(rhs_fml, data = tmp)
  sum(colnames(mm) != "(Intercept)")
}

verify_no_leakage <- function(fit, train_df_used, test_df) {
  ok_disjoint <- !any(rownames(train_df_used) %in% rownames(test_df))
  raw_pred_train <- safe_predict(fit, train_df_used, debias=FALSE)
  y_train <- train_df_used$y
  ref <- lm(y_train ~ raw_pred_train)
  co_ref <- coef(ref)
  if (is.null(fit$debias_fit)) {
    return(list(ok=FALSE, ok_disjoint=ok_disjoint, coef_close=FALSE))
  }
  co_pkg <- coef(fit$debias_fit)
  same_len <- length(co_pkg) == length(co_ref)
  coef_close <- same_len && all(abs(co_pkg - co_ref) < 1e-8)
  list(ok = (ok_disjoint && coef_close), ok_disjoint=ok_disjoint, coef_close=coef_close)
}

build_formulas <- function(p) {
  vars <- paste0("X", seq_len(p))
  main_terms <- paste(vars, collapse=" + ")
  int_terms  <- paste0("(", main_terms, ")^2")
  sq_terms   <- paste(sprintf("I(%s^2)", vars), collapse=" + ")
  list(
    main_effects   = as.formula(paste0("y ~ ", main_terms)),
    main_plus_int  = as.formula(paste0("y ~ ", int_terms)),
    full_quadratic = as.formula(paste0("y ~ ", int_terms, " + ", sq_terms))
  )
}

## Sparse truth
truth_for_p <- function(p) {
  list(
    active_main  = seq_len(min(3, p)),
    active_pairs = if (p >= 3) rbind(c(1,2), c(2,3)) else matrix(numeric(), ncol=2),
    active_sq    = if (p >= 2) 2L else integer(0)
  )
}

## Equal scales; AR-ish correlation
gen_X <- function(n, p, rho=0.5) {
  Z <- matrix(rnorm(n*p), n, p)
  if (p >= 2 && abs(rho) > 0) {
    for (j in 2:p) {
      Z[, j] <- rho * scale(Z[, j-1], TRUE, TRUE) + sqrt(1-rho^2) * scale(Z[, j], TRUE, TRUE)
    }
  }
  X <- Z
  colnames(X) <- paste0("X", seq_len(p))
  as.data.frame(X)
}

## Data with targeted in-sample R^2
make_sparse_data_R2 <- function(n, p, rho, target_R2,
                                active_main, active_pairs, active_sq,
                                coef_main=2, coef_int=1.5, coef_sq=1.0, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- gen_X(n, p, rho)
  mu <- 1
  y_sig <- rep(mu, n)
  for (j in active_main) y_sig <- y_sig + coef_main * X[[paste0("X", j)]]
  if (length(active_sq)) for (j in active_sq) y_sig <- y_sig + coef_sq * (X[[paste0("X", j)]]^2)
  if (nrow(active_pairs) > 0) {
    for (k in seq_len(nrow(active_pairs))) {
      i <- active_pairs[k,1]; j <- active_pairs[k,2]
      y_sig <- y_sig + coef_int * (X[[paste0("X", i)]] * X[[paste0("X", j)]])
    }
  }
  vs <- var(y_sig); sd_sig <- sqrt(max(vs, .Machine$double.eps))
  sd_eps <- sd_sig * sqrt((1 - target_R2) / target_R2)
  y <- y_sig + rnorm(n, 0, sd_eps)
  df <- cbind.data.frame(y=y, X)
  df$y_signal <- y_sig
  rownames(df) <- sprintf("row%06d", seq_len(n))
  df
}

## ------------------ RUN ------------------
rows <- list(); rid <- 1
leaks <- list()
TOTAL_FITS <- sum(sapply(P_GRID, function(p) length(MODELS)*length(RHO_GRID)*length(R2_GRID)*REPS))
fit_counter <- 0L
if (PROGRESS) pb <- txtProgressBar(min = 0, max = TOTAL_FITS, style = 3)

for (p in P_GRID) {
  fms <- build_formulas(p)
  thr <- truth_for_p(p)

  for (model in MODELS) {
    fm <- fms[[model]]
    p_full <- count_p_full(fm, p)

    ## choose n_train to cross the boundary; clamp [10, 80]
    n_grid_raw <- round(c(0.5, 0.75, 1.0, 1.25, 1.75, 2.5) * p_full)
    n_grid <- sort(unique(pmin(pmax(n_grid_raw, 10), 80)))

    for (rho in RHO_GRID) {
      for (R2_tgt in R2_GRID) {
        for (rep in seq_len(REPS)) {
          n_tr <- sample(n_grid, 1L)
          n_te <- NTEST

          df <- make_sparse_data_R2(
            n = n_tr + n_te, p = p, rho = rho, target_R2 = R2_tgt,
            active_main  = thr$active_main,
            active_pairs = thr$active_pairs,
            active_sq    = thr$active_sq,
            seed = 1000 + 7*rep + p + round(100*R2_tgt) + round(100*rho)
          )

          idx <- sample(seq_len(nrow(df)), size = n_tr)
          train_df <- df[idx, ]
          test_df  <- df[-idx, ]

          ## PRE-FILTER complete cases for THIS formula (avoid NA crash)
          mf_all <- model.frame(fm, data = train_df, na.action = na.pass)
          keep <- complete.cases(mf_all)
          if (sum(keep) < 2) {
            warning(sprintf("Skipping: n_train=%d but <2 complete cases for %s (p=%d)", n_tr, model, p))
            next
          }
          train_used <- train_df[keep, , drop = FALSE]
          n_used <- nrow(train_used)

          ## realized true R^2
          r2_true <- function(d) var(d$y_signal)/var(d$y)
          R2_train_true <- r2_true(train_used)
          R2_test_true  <- r2_true(test_df)

          fit_counter <- fit_counter + 1L
          if (PROGRESS) setTxtProgressBar(pb, fit_counter)
          if (fit_counter %% max(1L, ECHO_EVERY) == 0L) {
            log_msg("Start #%d | p=%d (p_full=%d) | model=%s | rho=%.2f | R2=%.2f | n_train=%d (used=%d) | ratio=%.2f",
                    fit_counter, p, p_full, model, rho, R2_tgt, n_tr, n_used, n_used/p_full)
          }

          ## ---- Fit wAIC ----
          t0 <- proc.time()[3]
          fit_waic <- SVEMnet(
            formula = fm, data = train_used,
            nBoot = NBOOT,
            glmnet_alpha = c(1),
            weight_scheme = "SVEM",
            objective = "wAIC",
            standardize = TRUE
          )
          t1 <- proc.time()[3]

          ## ---- Fit wSSE ----
          t2 <- proc.time()[3]
          fit_wsse <- SVEMnet(
            formula = fm, data = train_used,
            nBoot = NBOOT,
            glmnet_alpha = c(1),
            weight_scheme = "SVEM",
            objective = "wSSE",
            standardize = TRUE
          )
          t3 <- proc.time()[3]

          ## Predictions
          pr_tr_raw_waic <- safe_predict(fit_waic, train_used, debias = FALSE)
          pr_te_raw_waic <- safe_predict(fit_waic, test_df,   debias = FALSE)
          pr_tr_deb_waic <- safe_predict(fit_waic, train_used, debias = TRUE)
          pr_te_deb_waic <- safe_predict(fit_waic, test_df,    debias = TRUE)

          pr_tr_raw_wsse <- safe_predict(fit_wsse, train_used, debias = FALSE)
          pr_te_raw_wsse <- safe_predict(fit_wsse, test_df,   debias = FALSE)
          pr_tr_deb_wsse <- safe_predict(fit_wsse, train_used, debias = TRUE)
          pr_te_deb_wsse <- safe_predict(fit_wsse, test_df,    debias = TRUE)

          ## RMSEs
          rmse_tr_raw_waic <- rmse(train_used$y, pr_tr_raw_waic)
          rmse_te_raw_waic <- rmse(test_df$y,    pr_te_raw_waic)
          rmse_tr_deb_waic <- rmse(train_used$y, pr_tr_deb_waic)
          rmse_te_deb_waic <- rmse(test_df$y,    pr_te_deb_waic)

          rmse_tr_raw_wsse <- rmse(train_used$y, pr_tr_raw_wsse)
          rmse_te_raw_wsse <- rmse(test_df$y,    pr_te_raw_wsse)
          rmse_tr_deb_wsse <- rmse(train_used$y, pr_tr_deb_wsse)
          rmse_te_deb_wsse <- rmse(test_df$y,    pr_te_deb_wsse)

          log_msg("Done  #%d | AIC raw=%.3f deb=%.3f | SSE raw=%.3f deb=%.3f | Δ_deb(raw): AIC=%.3f SSE=%.3f | %.1fs+%.1fs",
                  fit_counter, rmse_te_raw_waic, rmse_te_deb_waic,
                  rmse_te_raw_wsse, rmse_te_deb_wsse,
                  rmse_te_deb_waic - rmse_te_raw_waic,
                  rmse_te_deb_wsse - rmse_te_raw_wsse,
                  (t1 - t0), (t3 - t2))

          ## Save row
          rows[[rid]] <- data.frame(
            p = p, p_full = p_full, model = model, rho = rho, R2_target = R2_tgt,
            n_train = n_tr, n_train_used = n_used, n_test = n_te,
            ratio_n_over_p = n_used / p_full,
            above_boundary = as.integer(n_used > p_full),
            R2_true_train = R2_train_true, R2_true_test = R2_test_true,
            rmse_test_raw_waic = rmse_te_raw_waic,
            rmse_test_deb_waic = rmse_te_deb_waic,
            rmse_test_raw_wsse = rmse_te_raw_wsse,
            rmse_test_deb_wsse = rmse_te_deb_wsse,
            rmse_train_raw_waic = rmse_tr_raw_waic,
            rmse_train_deb_waic = rmse_tr_deb_waic,
            rmse_train_raw_wsse = rmse_tr_raw_wsse,
            rmse_train_deb_wsse = rmse_tr_deb_wsse,
            time_sec_waic = round(t1 - t0, 2),
            time_sec_wsse = round(t3 - t2, 2),
            stringsAsFactors = FALSE
          )
          rid <- rid + 1

          ## light leakage check
          if (rep <= 2) {
            chk <- verify_no_leakage(fit_waic, train_used, test_df)
            leaks[[length(leaks) + 1]] <- data.frame(
              p = p, model = model, rho = rho, R2_target = R2_tgt, n_train = n_tr, n_train_used = n_used,
              ok = chk$ok, ok_disjoint = chk$ok_disjoint, coef_close = chk$coef_close
            )
          }
        }
      }
    }
  }
}

if (!length(rows)) stop("No successful fits were recorded. Check logs/warnings above.")
res <- do.call(rbind, rows)
leak_df <- if (length(leaks)) do.call(rbind, leaks) else NULL

## Ratio bins for readability on plots/tables
res$ratio_bin <- cut(
  res$ratio_n_over_p,
  breaks = c(-Inf, 0.75, 1.0, 1.25, 1.75, Inf),
  labels = c("<0.75", "0.75–1.0", "1.0–1.25", "1.25–1.75", ">1.75"),
  right = TRUE
)

## ------------------ SUMMARY (with 4 win-rates) ------------------

## CI helper
ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x); if (!n) return(c(mean=NA, lo=NA, hi=NA, n=0))
  m <- mean(x); se <- sd(x)/sqrt(n)
  c(mean=m, lo=m-1.96*se, hi=m+1.96*se, n=n)
}

## Grouping key
grp_vars <- c("model","R2_target","ratio_bin","above_boundary")

## Build summary per group using by() to avoid unpack glitches
summ <- do.call(rbind, by(res, INDICES = res[grp_vars], FUN = function(d){
  ## win indicators (row-wise)
  deb_win_waic <- as.numeric(d$rmse_test_deb_waic < d$rmse_test_raw_waic)
  deb_win_wsse <- as.numeric(d$rmse_test_deb_wsse < d$rmse_test_raw_wsse)
  deb_lose_waic <- as.numeric(d$rmse_test_deb_waic > d$rmse_test_raw_waic)
  deb_lose_wsse <- as.numeric(d$rmse_test_deb_wsse > d$rmse_test_raw_wsse)

  waic_win_raw <- as.numeric(d$rmse_test_raw_waic < d$rmse_test_raw_wsse)
  waic_win_deb <- as.numeric(d$rmse_test_deb_waic < d$rmse_test_deb_wsse)
  wsse_win_raw <- as.numeric(d$rmse_test_raw_wsse < d$rmse_test_raw_waic)
  wsse_win_deb <- as.numeric(d$rmse_test_deb_wsse < d$rmse_test_deb_waic)

  ## pooled winrates (redundant on purpose)
  winrate_debias_TRUE  <- mean(c(deb_win_waic,  deb_win_wsse),  na.rm = TRUE)
  winrate_debias_FALSE <- mean(c(deb_lose_waic, deb_lose_wsse), na.rm = TRUE)
  winrate_wAIC         <- mean(c(waic_win_raw,  waic_win_deb),  na.rm = TRUE)
  winrate_wSSE         <- mean(c(wsse_win_raw,  wsse_win_deb),  na.rm = TRUE)

  ## RMSE means + CI
  a_raw <- ci95(d$rmse_test_raw_waic)
  a_deb <- ci95(d$rmse_test_deb_waic)
  s_raw <- ci95(d$rmse_test_raw_wsse)
  s_deb <- ci95(d$rmse_test_deb_wsse)

  ## debias deltas (for context)
  dA <- ci95(d$rmse_test_deb_waic - d$rmse_test_raw_waic)
  dS <- ci95(d$rmse_test_deb_wsse - d$rmse_test_raw_wsse)

  ## realized true R2
  rtr <- ci95(d$R2_true_train); rte <- ci95(d$R2_true_test)

  ## timings
  tt  <- ci95(rowMeans(cbind(d$time_sec_waic, d$time_sec_wsse), na.rm = TRUE))

  data.frame(
    model = d$model[1],
    R2_target = d$R2_target[1],
    ratio_bin = d$ratio_bin[1],
    above_boundary = d$above_boundary[1],

    R2_true_train_mean = rtr["mean"], R2_true_test_mean = rte["mean"],

    rmse_raw_waic_mean = a_raw["mean"],  rmse_raw_waic_ci95_lo = a_raw["lo"],  rmse_raw_waic_ci95_hi = a_raw["hi"],
    rmse_deb_waic_mean = a_deb["mean"],  rmse_deb_waic_ci95_lo = a_deb["lo"],  rmse_deb_waic_ci95_hi = a_deb["hi"],
    rmse_raw_wsse_mean = s_raw["mean"],  rmse_raw_wsse_ci95_lo = s_raw["lo"],  rmse_raw_wsse_ci95_hi = s_raw["hi"],
    rmse_deb_wsse_mean = s_deb["mean"],  rmse_deb_wsse_ci95_lo = s_deb["lo"],  rmse_deb_wsse_ci95_hi = s_deb["hi"],

    ## the four requested winrates
    winrate_debias_TRUE  = winrate_debias_TRUE,
    winrate_debias_FALSE = winrate_debias_FALSE,
    winrate_wAIC         = winrate_wAIC,
    winrate_wSSE         = winrate_wSSE,

    ## debias deltas by objective (context)
    deb_raw_diff_AIC_mean = dA["mean"], deb_raw_diff_AIC_ci95_lo = dA["lo"], deb_raw_diff_AIC_ci95_hi = dA["hi"],
    deb_raw_diff_SSE_mean = dS["mean"], deb_raw_diff_SSE_ci95_lo = dS["lo"], deb_raw_diff_SSE_ci95_hi = dS["hi"],

    mean_time = tt["mean"],
    row.names = NULL,
    check.names = FALSE
  )
}))
rownames(summ) <- NULL
summ <- summ[order(summ$model, summ$R2_target, summ$ratio_bin, summ$above_boundary), ]

## Leakage summary
cat("\n---- Leakage summary (should be ~100%) ----\n")
if (!is.null(leak_df)) {
  print(aggregate(cbind(ok, ok_disjoint, coef_close) ~ model, data = leak_df,
                  FUN = function(x) mean(x, na.rm=TRUE)), row.names = FALSE)
} else {
  cat("No leakage samples collected.\n")
}

## Print compact summary columns
ord <- c("model","R2_target","ratio_bin","above_boundary",
         "R2_true_train_mean","R2_true_test_mean",
         "rmse_raw_waic_mean","rmse_deb_waic_mean",
         "rmse_raw_wsse_mean","rmse_deb_wsse_mean",
         "winrate_debias_TRUE","winrate_debias_FALSE",
         "winrate_wAIC","winrate_wSSE",
         "deb_raw_diff_AIC_mean","deb_raw_diff_AIC_ci95_lo","deb_raw_diff_AIC_ci95_hi",
         "deb_raw_diff_SSE_mean","deb_raw_diff_SSE_ci95_lo","deb_raw_diff_SSE_ci95_hi",
         "mean_time")
cat("\n================ SUMMARY (with 4 win-rates) ================\n")
print(summ[, ord], row.names = FALSE, digits = 4)

## ------------------ Plots (optional) ------------------

## Debias − Raw (negative helps) for each objective
plot_df <- within(res, {
  deb_raw_AIC  <- rmse_test_deb_waic - rmse_test_raw_waic
  deb_raw_SSE  <- rmse_test_deb_wsse - rmse_test_raw_wsse
})
ggA <- ggplot(plot_df, aes(ratio_n_over_p, deb_raw_AIC)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = TRUE, method = "loess", span = 0.7) +
  facet_grid(R2_target ~ model) +
  labs(x = expression(paste(n[train[used]], " / ", p[full])),
       y = "Debias − Raw (TEST RMSE) under wAIC",
       title = "When does debiasing help? (wAIC)") +
  theme_bw()

ggS <- ggplot(plot_df, aes(ratio_n_over_p, deb_raw_SSE)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = TRUE, method = "loess", span = 0.7) +
  facet_grid(R2_target ~ model) +
  labs(x = expression(paste(n[train[used]], " / ", p[full])),
       y = "Debias − Raw (TEST RMSE) under wSSE",
       title = "When does debiasing help? (wSSE)") +
  theme_bw()

print(ggA)
print(ggS)

## ------------------ Export (for your own modeling) ------------------
keep_cols <- c("p","p_full","model","rho","R2_target",
               "n_train","n_train_used","n_test",
               "ratio_n_over_p","above_boundary",
               "R2_true_train","R2_true_test",
               "rmse_test_raw_waic","rmse_test_deb_waic",
               "rmse_test_raw_wsse","rmse_test_deb_wsse",
               "time_sec_waic","time_sec_wsse")
res_export <- res[, keep_cols]
cat("\nFirst few rows of res_export:\n")
print(utils::head(res_export), row.names = FALSE)

## Write to CSV if you like:
## write.csv(res_export, "res_export.csv", row.names = FALSE)

## ---- Oneway-style plots: block on Row and show ANOM-like needles --------
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# If your results frame is called res_export, just do:
# res <- res_export

# 1) Long format with a Row (block) ID and log(RMSE)
df <- res_export %>%
  mutate(Row = row_number()) %>%                           # block label = run id
  pivot_longer(
    cols = c(rmse_test_deb_waic, rmse_test_deb_wsse,
             rmse_test_raw_waic, rmse_test_raw_wsse),
    names_to  = "settings",
    values_to = "rmse"
  ) %>%
  mutate(
    settings = factor(settings,
                      levels = c("rmse_test_deb_waic","rmse_test_deb_wsse",
                                 "rmse_test_raw_waic","rmse_test_raw_wsse")),
    lmrse    = log(pmax(rmse, .Machine$double.eps))
  )

# 2) Block-centering (remove Row mean, then add back overall mean so the scale is familiar)
overall_mean <- mean(df$lmrse, na.rm = TRUE)
df <- df %>%
  group_by(Row) %>%
  mutate(lmrse_block = lmrse - mean(lmrse, na.rm = TRUE) + overall_mean) %>%
  ungroup()

# 3) Fit two-way fixed effects model to get a pooled MSE (settings + block)
fit <- lm(lmrse_block ~ settings + factor(Row), data = df)
MSE <- summary(fit)$sigma^2
df_resid <- df.residual(fit)

# Effective sample size per setting (usually each setting appears once per Row)
n_i <- as.numeric(table(df$settings))
SE_i <- sqrt(MSE / n_i)            # SE for each setting mean (RCBD-style)
names(SE_i) <- names(table(df$settings))

# 4) Means per setting and an ANOM-like CI
alpha <- 0.05
tcrit <- qt(1 - alpha/2, df = df_resid)

means <- df %>%
  group_by(settings) %>%
  summarise(mean_lmrse = mean(lmrse_block, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    SE   = SE_i[as.character(settings)],
    lo   = mean_lmrse - tcrit * SE,
    hi   = mean_lmrse + tcrit * SE
  )

# Grand mean and a constant decision band
grand <- mean(df$lmrse_block, na.rm = TRUE)
SE_bar <- mean(SE_i)                      # one band for the panel
UDL <- grand + tcrit * SE_bar
LDL <- grand - tcrit * SE_bar
band_df <- data.frame(LDL = LDL, UDL = UDL)

# 5) Panel 1: “Block Centered lmrse” (violin + dots)
p1 <- ggplot(df, aes(settings, lmrse_block)) +
  geom_violin(fill = "grey90", color = "grey40", scale = "width") +
  geom_jitter(width = 0.10, height = 0, alpha = 0.25, size = 1) +
  geom_hline(yintercept = grand, linetype = 2, color = "grey35") +
  labs(title = "Oneway Analysis of lmrse by settings",
       subtitle = "Block: Row (run ID) — block-centered values",
       y = "Block-Centered lmrse", x = "settings") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 6) Panel 2: “Analysis of Means” (lollipops + CI + decision band)
p2 <- ggplot(means, aes(settings, mean_lmrse)) +
  # decision band
  geom_rect(data = band_df,
            aes(xmin = -Inf, xmax = Inf, ymin = LDL, ymax = UDL),
            inherit.aes = FALSE, alpha = 0.12, fill = "skyblue") +
  geom_hline(yintercept = grand, color = "grey25") +
  # “needles”
  geom_segment(aes(x = settings, xend = settings, y = grand, yend = mean_lmrse),
               linewidth = 0.6, color = "black") +
  # point + CI
  geom_point(size = 3, color = "firebrick") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  annotate("text", x = Inf, y = UDL, hjust = 1.05, vjust = -0.4,
           label = sprintf("UDL=%.4f", UDL), size = 3.3) +
  annotate("text", x = Inf, y = grand, hjust = 1.05, vjust = -0.4,
           label = sprintf("Avg = %.4f", grand), size = 3.3) +
  annotate("text", x = Inf, y = LDL, hjust = 1.05, vjust = 1.2,
           label = sprintf("LDL=%.4f", LDL), size = 3.3) +
  labs(title = "Analysis of Means (ANOM-style)",
       subtitle = expression(paste("95% t-band; pooled MSE from ",
                                   "lm(lmrse_block ~ settings + block)")),
       y = "lmrse", x = "settings") +
  coord_cartesian(ylim = range(c(means$lo, means$hi, LDL, UDL))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 7)
(p1 / p2) + plot_annotation(
  title = "Blocked Oneway (Row as Block) for lmrse vs. settings",
  theme = theme(plot.title = element_text(face = "bold"))
)

