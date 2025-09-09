###############################################################################
# Parallel SVEMnet simulation (Windows-safe, 8 cores)
# - Run sizes: every integer 14..26
# - Models: SVEMnet wAIC (lasso), wBIC (lasso), wSSE (lasso), glmnet CV lasso
# - Outputs: clear console summaries + transition table (AIC vs BIC by n_total)
###############################################################################

# ----------- Packages -----------
suppressPackageStartupMessages({
  library(SVEMnet)     # your package
  library(glmnet)      # benchmark CV
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(furrr)
})

# ----------- User-configurable knobs -----------
CORES        <- 15                 # Windows-safe parallel workers
OUT_ITERS    <- 300                # repeats per n_total (adjust)
N_TOTAL_SEQ  <- 14:26             # every integer from 14 to 26
N_BOOT       <- 300               # SVEMnet bootstrap reps
HOLDOUT_N    <- 10000               # holdout size per run
DIST_N       <- 10000             # size to estimate sd(TrueY) for R^2
ALPHAS_LASSO <- 1                 # lasso only

set.seed(20251201)                # top-level seed for reproducibility

# ----------- Parallel plan (Windows-friendly) -----------
future::plan(future::multisession, workers = CORES)

# ----------- Utilities -----------
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(stats::rgamma(n * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
  x / rowSums(x)
}

# sample mixture points A..D on simplex with bounds:
# A in [0.1, 0.4], B,C,D in [0, 0.8]
sample_mixture <- function(n) {
  out <- matrix(NA_real_, nrow = n, ncol = 4)
  i <- 1
  while (i <= n) {
    x <- as.numeric(rdirichlet(1, c(1, 1, 1, 1)))
    if (x[1] >= 0.1 && x[1] <= 0.4 && all(x[2:4] <= 0.8)) {
      out[i, ] <- x
      i <- i + 1
    }
  }
  out <- as.data.frame(out)
  names(out) <- c("A", "B", "C", "D")
  out
}

# Generate pv vector like your JSL (Laplace-ish with sparsity patterns)
gen_pv <- function() {
  rexp_signed <- function() stats::rexp(1) - stats::rexp(1)
  p1_4 <- replicate(4, rexp_signed())
  p5_8 <- replicate(4, rexp_signed() * rbinom(1, 1, 0.8))
  p9_25 <- replicate(17, rexp_signed() * rbinom(1, 1, 0.5))
  c(p1_4, p5_8, p9_25)
}

# True response surface (mirrors your JSL structure)
true_response <- function(df, pv) {
  # df has A,B,C,D and E (numeric {0, 0.002})
  s  <- 1 - 0.1
  zA <- (df$A - 0.1) / s
  zB <- df$B / s
  zC <- df$C / s
  zD <- df$D / s
  Esign <- ifelse(df$E == 0, 1, -1)

  # linear + E interactions
  part1 <- pv[1]*zA + pv[2]*zB + pv[3]*zC + pv[4]*zD +
    (pv[5]*zA + pv[6]*zB + pv[7]*zC + pv[8]*zD) * Esign

  # pair terms (scaled by 4)
  part2 <- 4 * ( pv[9]*zA*zB + pv[10]*zA*zC + pv[11]*zA*zD +
                   pv[12]*zB*zC + pv[13]*zB*zD + pv[14]*zC*zD )

  # triple terms (scaled by 27)
  part3 <- 27 * ( pv[15]*zA*zB*zC + pv[16]*zA*zB*zD +
                    pv[17]*zA*zC*zD + pv[18]*zB*zC*zD )

  # untargeted Scheffé cubic terms (scaled by 27)
  part4 <- 27 * ( pv[19]*zB*zA*(zA - zB) +
                    pv[20]*zC*zA*(zA - zC) +
                    pv[21]*zC*zB*(zB - zC) +
                    pv[22]*zD*zA*(zA - zD) +
                    pv[23]*zD*zB*(zB - zD) +
                    pv[24]*zD*zC*(zC - zD) )

  # quartic product (scaled by 256)
  part5 <- 256 * pv[25] * zA*zB*zC*zD

  part1 + part2 + part3 + part4 + part5
}

# Build one simulated dataset (training + holdout + noise), return per-setting rows
run_one <- function(run_id, n_total) {
  # 1) Random coefficients and theoretical R^2
  pv  <- gen_pv()
  R2_choices <- c(0.3, 0.5, 0.7, 0.9)
  r2 <- sample(R2_choices, size = 1)

  # 2) Estimate sd(TrueY) over the domain
  dist_pts <- sample_mixture(DIST_N)
  dist_pts$E <- sample(c(0, 0.002), size = DIST_N, replace = TRUE)
  dist_true <- true_response(dist_pts, pv)
  y_sd_global <- stats::sd(dist_true)

  # 3) Training
  tr <- sample_mixture(n_total)
  tr$E <- sample(c(0, 0.002), size = n_total, replace = TRUE)
  tr_true <- true_response(tr, pv)

  # Noise to match R^2
  err_sd <- y_sd_global * sqrt((1 - r2) / r2)
  Y <- tr_true + stats::rnorm(n_total, sd = err_sd)

  train <- tr %>%
    mutate(Y = Y) %>%
    mutate(E = factor(E))  # treat E categorically

  # 4) Holdout (evaluate on noiseless TrueY)
  hold <- sample_mixture(HOLDOUT_N)
  hold$E <- sample(c(0, 0.002), size = HOLDOUT_N, replace = TRUE)
  hold_true <- true_response(hold, pv)
  hold <- hold %>% mutate(E = factor(E))

  # Scaling (per-run)
  sd_hold_true <- stats::sd(hold_true)

  # 5) Common formula (rich enough but not too huge)
  #    Includes main effects, all 2-way interactions among A..D, interactions with E,
  #    plus 3-way among A..D and 4-way ABCD (lasso will regularize strongly at small n)
  form <- Y ~ (A + B + C + D + E)^2 + A:B:C + A:B:D + A:C:D + B:C:D + A:B:C:D

  # 6) Fit models
  # SVEMnet (lasso) with different objectives
  fit_svem <- function(objective) {
    SVEMnet::SVEMnet(
      formula     = form,
      data        = train,
      glmnet_alpha= ALPHAS_LASSO,
      nBoot       = N_BOOT,
      objective   = objective,
      standardize = TRUE
    )
  }

  m_aic <- fit_svem("wAIC")
  m_bic <- fit_svem("wBIC")
  m_sse <- fit_svem("wSSE")

  # glmnet CV lasso benchmark
  mm_train <- model.matrix(form, data = train)
  y_train  <- train$Y
  # drop intercept column (glmnet adds it internally)
  if ("(Intercept)" %in% colnames(mm_train)) {
    mm_train <- mm_train[, colnames(mm_train) != "(Intercept)", drop = FALSE]
  }
  cvfit <- glmnet::cv.glmnet(
    x = mm_train, y = y_train, alpha = 1, nfolds = 5,
    standardize = TRUE, intercept = TRUE, type.measure = "mse"
  )

  # 7) Predictions on holdout (compare to TrueY)
  preds <- list()

  preds$SVEM_wAIC_lasso <- as.numeric(predict(m_aic, newdata = hold))
  preds$SVEM_wBIC_lasso <- as.numeric(predict(m_bic, newdata = hold))
  preds$SVEM_wSSE_lasso <- as.numeric(predict(m_sse, newdata = hold))

  mm_hold <- model.matrix(update(form, NULL ~ .), data = hold) # same RHS as form
  if ("(Intercept)" %in% colnames(mm_hold)) {
    mm_hold <- mm_hold[, colnames(mm_hold) != "(Intercept)", drop = FALSE]
  }
  preds$glmnet_cv_lasso <- as.numeric(predict(cvfit, newx = mm_hold, s = "lambda.min"))

  # 8) Metrics (Holdout; normalized by sd(TrueY_holdout))
  metric_one <- function(yhat, y_true, sd_scale) {
    err  <- yhat - y_true
    rase <- sqrt(mean(err^2))
    aae  <- mean(abs(err))
    c(
      NRASE_Holdout = rase / sd_scale,
      NAAE_Holdout  = aae  / sd_scale
    )
  }

  rows <- purrr::imap_dfr(preds, function(p, name) {
    m <- metric_one(p, hold_true, sd_hold_true)
    tibble(
      RunID           = factor(run_id),
      n_total         = n_total,
      TheoreticalR2   = factor(r2),
      Holdout_SDTrueY = sd_hold_true,
      Setting         = factor(name, levels = c("SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "SVEM_wSSE_lasso", "glmnet_cv_lasso")),
      NRASE_Holdout   = unname(m["NRASE_Holdout"]),
      NAAE_Holdout    = unname(m["NAAE_Holdout"])
    )
  })

  rows
}

# ----------- Build job grid and run in parallel -----------
jobs <- tibble(
  n_total = rep(N_TOTAL_SEQ, each = OUT_ITERS),
  iter    = rep(seq_len(OUT_ITERS), times = length(N_TOTAL_SEQ))
) %>%
  mutate(RunID = sprintf("run%05d", row_number()))

cat(sprintf("Running %d jobs across %d workers... (n_total in %s)\n",
            nrow(jobs), CORES, paste(range(N_TOTAL_SEQ), collapse = ":")))
flush.console()

df_list <- furrr::future_pmap(
  .l = list(jobs$RunID, jobs$n_total),
  .f = ~ run_one(..1, ..2),
  .options = furrr::furrr_options(seed = TRUE)
)

df <- dplyr::bind_rows(df_list)
df$Setting <- droplevels(df$Setting)

# ----------- Helper summaries -----------
se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

summary_by <- function(d) {
  d %>%
    group_by(Setting) %>%
    summarise(
      runs      = n(),
      mean_NRASE= mean(NRASE_Holdout, na.rm = TRUE),
      se_NRASE  = se(NRASE_Holdout),
      mean_NAAE = mean(NAAE_Holdout, na.rm = TRUE),
      se_NAAE   = se(NAAE_Holdout),
      .groups   = "drop"
    ) %>%
    arrange(mean_NRASE)
}

win_rate_tbl <- function(d) {
  winners <- d %>%
    group_by(RunID, n_total) %>%
    slice_min(NRASE_Holdout, with_ties = TRUE) %>%
    ungroup() %>%
    mutate(flag = 1L)

  d %>%
    select(RunID, n_total, Setting) %>%
    left_join(winners %>% select(RunID, n_total, Setting, flag),
              by = c("RunID", "n_total", "Setting")) %>%
    group_by(Setting) %>%
    summarise(
      wins     = sum(replace(flag, is.na(flag), 0L)),
      runs     = n_distinct(paste(RunID, n_total)),
      win_rate = wins / runs,
      .groups  = "drop"
    ) %>%
    arrange(desc(win_rate))
}

avg_rank_tbl <- function(d) {
  d %>%
    group_by(RunID, n_total) %>%
    mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
    ungroup() %>%
    group_by(Setting) %>%
    summarise(
      mean_rank = mean(rk, na.rm = TRUE),
      se_rank   = se(rk),
      .groups   = "drop"
    ) %>%
    arrange(mean_rank)
}

paired_compare <- function(d, a, b, label = "") {
  wide <- d %>%
    filter(Setting %in% c(a, b)) %>%
    select(RunID, n_total, Setting, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = Setting, values_from = NRASE_Holdout)
  if (!all(c(a, b) %in% names(wide))) {
    cat("Skipping paired compare (missing columns)\n")
    return(invisible(NULL))
  }
  delta <- wide[[a]] - wide[[b]]
  delta <- delta[is.finite(delta)]
  if (length(delta) < 3) {
    cat("Not enough pairs for", label, "\n"); return(invisible(NULL))
  }
  tt <- t.test(delta)
  wt <- suppressWarnings(wilcox.test(delta, exact = FALSE))
  cat(sprintf("\n-- %s: %s - %s (N=%d) --\n  meanΔ = %+0.4f   t(%d) = %0.2f  p = %.3g   |   medianΔ = %+0.4f  Wilcoxon p = %.3g\n",
              label, a, b, length(delta),
              mean(delta), as.integer(tt$parameter), unname(tt$statistic), tt$p.value,
              stats::median(delta), wt$p.value))
}

pct_vs_baseline <- function(d, baseline = "glmnet_cv_lasso") {
  b <- baseline
  wide <- d %>%
    select(RunID, n_total, Setting, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = Setting, values_from = NRASE_Holdout)

  others <- intersect(c("SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "SVEM_wSSE_lasso"), names(wide))
  out <- map_dfr(others, function(s) {
    if (!all(c(b, s) %in% names(wide))) return(NULL)
    rel <- 100 * (wide[[s]] - wide[[b]]) / wide[[b]]
    tibble(Setting = s, mean_pct = mean(rel, na.rm = TRUE), se_pct = se(rel))
  })
  arrange(out, mean_pct)
}

# ----------- Printouts -----------
cat("\n==================== OVERALL ====================\n")
overall <- df
print(summary_by(overall))
cat("\nWin rate (tie-aware):\n"); print(win_rate_tbl(overall))
cat("\nAverage rank:\n");      print(avg_rank_tbl(overall))
paired_compare(overall, "SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "Overall (AIC vs BIC)")
cat("\nPercent Δ vs benchmark (glmnet_cv_lasso):\n"); print(pct_vs_baseline(overall))

# ---- Transition block: AIC vs BIC for every n_total in 14..26 ----
cat("\n================ TRANSITION: n_total 14..26 (AIC vs BIC) ================\n")

paired_deltas_by_n <- df %>%
  filter(Setting %in% c("SVEM_wAIC_lasso", "SVEM_wBIC_lasso")) %>%
  select(RunID, n_total, Setting, NRASE_Holdout) %>%
  distinct() %>%
  pivot_wider(names_from = Setting, values_from = NRASE_Holdout) %>%
  arrange(n_total, RunID)

transition_tbl <- paired_deltas_by_n %>%
  group_by(n_total) %>%
  group_modify(~{
    d   <- .x$`SVEM_wAIC_lasso` - .x$`SVEM_wBIC_lasso`
    rel <- 100 * d / .x$`SVEM_wBIC_lasso`

    t_p <- if (sum(is.finite(d)) > 2) t.test(d)$p.value else NA_real_
    w_p <- if (any(is.finite(d))) suppressWarnings(wilcox.test(d, exact = FALSE)$p.value) else NA_real_
    N   <- sum(is.finite(d))

    tibble(
      N            = N,
      mean_delta   = mean(d,   na.rm = TRUE),
      median_delta = median(d, na.rm = TRUE),
      t_p          = t_p,
      wilcox_p     = w_p,
      p_AIC_worse  = mean(d > 0, na.rm = TRUE),
      mean_pct     = mean(rel, na.rm = TRUE),
      se_pct       = se(rel)
    )
  }) %>%
  ungroup() %>%
  arrange(n_total)

print(transition_tbl)

# Tie-aware wins per n_total (AIC/BIC/benchmark context)
winners_by_n <- df %>%
  filter(Setting %in% c("SVEM_wAIC_lasso","SVEM_wBIC_lasso","glmnet_cv_lasso")) %>%
  group_by(n_total, RunID) %>%
  slice_min(NRASE_Holdout, with_ties = TRUE) %>%
  ungroup() %>%
  count(n_total, Setting, name = "wins") %>%
  group_by(n_total) %>%
  mutate(total_runs = sum(wins), win_rate = wins / total_runs) %>%
  arrange(n_total, desc(win_rate))

cat("\n--- Tie-aware win rates by n_total (AIC/BIC + benchmark) ---\n")
print(winners_by_n)

# Suggested threshold: earliest n_total with meanΔ<=0 and Wilcoxon p<=0.05
boundary_rows <- transition_tbl %>%
  filter(is.finite(mean_delta), mean_delta <= 0,
         is.finite(wilcox_p),  wilcox_p <= 0.05) %>%
  arrange(n_total)

cat("\n--- Suggested transition threshold(s) ---\n")
if (nrow(boundary_rows)) {
  cat("Earliest n_total where AIC ≤ BIC (mean) AND Wilcoxon p ≤ 0.05:\n")
  print(boundary_rows %>% select(n_total, mean_delta, median_delta, wilcox_p, p_AIC_worse, mean_pct, se_pct))
} else {
  cat("No n_total in 14..26 met both conditions (mean_delta≤0 AND Wilcoxon p≤0.05).\n")
  cat("Use the table above to choose a cutoff consistent with your regret tolerance.\n")
}

cat("\n================= RECOMMENDATION =================\n")
cat("* Use the TRANSITION table (14..26) to pick the AIC/BIC cutoff.\n")
cat("* Benchmark (glmnet CV lasso) is context only.\n")
cat("\nDone.\n")
