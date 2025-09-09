###############################################################################
# Parallel SVEMnet vs CV (with relax/standard and alpha scenarios)
# - Train sizes: 15, 25, 35, 45, 55, 65, 75
# - Methods: SVEMnet objective in {auto, wAIC, wBIC} x relaxed {TRUE, FALSE}
#            glmnet_with_cv with relaxed {TRUE, FALSE}
# - Alpha scenarios:
#       "lasso" => glmnet_alpha = 1
#       "mix"   => glmnet_alpha = c(0.5, 1)
# - Holdout metrics versus the true (noiseless) surface
###############################################################################

suppressPackageStartupMessages({
  library(SVEMnet)     # SVEMnet(), glmnet_with_cv(), predict_cv()
  library(glmnet)      # backend
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(furrr)
  library(ggplot2)
  library(scales)
})

# ------------------- knobs -------------------
CORES        <- 12                    # parallel workers
OUT_ITERS    <- 50                    # repeats per n_total
N_TOTAL_SEQ  <- seq(15, 75, by = 10)  # train sizes to sweep
N_BOOT       <- 300                   # SVEMnet bootstrap reps
HOLDOUT_N    <- 10000                 # holdout size per run
DIST_N       <- 10000                 # for estimating sd(TrueY) to set noise

ALPHAS_LASSO <- 1
ALPHAS_MIX   <- c(0.5, 1)

# glmnet_with_cv wrapper knobs
CV_NFOLDS    <- 10
CV_REPEATS   <- 3
CV_CHOICE    <- "1se"

set.seed(20251202)
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

# Generate pv vector
gen_pv <- function() {
  rexp_signed <- function() stats::rexp(1) - stats::rexp(1)
  p1_4  <- replicate(4, rexp_signed())
  p5_8  <- replicate(4, rexp_signed() * rbinom(1, 1, 0.8))
  p9_25 <- replicate(17, rexp_signed() * rbinom(1, 1, 0.5))
  c(p1_4, p5_8, p9_25)
}

# True response surface
true_response <- function(df, pv) {
  s  <- 1 - 0.1
  zA <- (df$A - 0.1) / s
  zB <- df$B / s
  zC <- df$C / s
  zD <- df$D / s
  Esign <- ifelse(df$E == 0, 1, -1)

  part1 <- pv[1]*zA + pv[2]*zB + pv[3]*zC + pv[4]*zD +
    (pv[5]*zA + pv[6]*zB + pv[7]*zC + pv[8]*zD) * Esign

  part2 <- 4 * ( pv[9]*zA*zB + pv[10]*zA*zC + pv[11]*zA*zD +
                   pv[12]*zB*zC + pv[13]*zB*zD + pv[14]*zC*zD )

  part3 <- 27 * ( pv[15]*zA*zB*zC + pv[16]*zA*zB*zD +
                    pv[17]*zA*zC*zD + pv[18]*zB*zC*zD )

  part4 <- 27 * ( pv[19]*zB*zA*(zA - zB) +
                    pv[20]*zC*zA*(zA - zC) +
                    pv[21]*zC*zB*(zB - zC) +
                    pv[22]*zD*zA*(zA - zD) +
                    pv[23]*zD*zB*(zB - zD) +
                    pv[24]*zD*zC*(zC - zD) )

  part5 <- 256 * pv[25] * zA*zB*zC*zD

  part1 + part2 + part3 + part4 + part5
}

# ----------- One simulated dataset (does both alpha scenarios) -----------
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

  # 3) Training with noise to hit target R^2
  tr <- sample_mixture(n_total)
  tr$E <- sample(c(0, 0.002), size = n_total, replace = TRUE)
  tr_true <- true_response(tr, pv)

  err_sd <- y_sd_global * sqrt((1 - r2) / r2)
  Y <- tr_true + stats::rnorm(n_total, sd = err_sd)

  train <- tr %>%
    mutate(Y = Y) %>%
    mutate(E = factor(E))

  # 4) Holdout (noiseless truth)
  hold <- sample_mixture(HOLDOUT_N)
  hold$E <- sample(c(0, 0.002), size = HOLDOUT_N, replace = TRUE)
  hold_true <- true_response(hold, pv)
  hold <- hold %>% mutate(E = factor(E))
  sd_hold_true <- stats::sd(hold_true)

  # 5) Common formula
  form <- Y ~ (A + B + C + D + E)^2 + A:B:C + A:B:D + A:C:D + B:C:D + A:B:C:D

  # metric helper
  metric_one <- function(yhat, y_true, sd_scale) {
    err  <- yhat - y_true
    c(NRASE_Holdout = sqrt(mean(err^2)) / sd_scale,
      NAAE_Holdout  = mean(abs(err))   / sd_scale)
  }

  # ---- inner runner for one alpha scenario ----
  fit_scenario <- function(alpha_vec, alpha_tag) {

    # SVEM fits: objective in {auto, wAIC, wBIC} x relaxed {FALSE, TRUE}
    fit_svem <- function(obj, relax_flag, label) {
      withCallingHandlers({
        SVEMnet::SVEMnet(
          formula        = form,
          data           = train,
          glmnet_alpha   = alpha_vec,
          nBoot          = N_BOOT,
          objective      = obj,
          weight_scheme  = "SVEM",
          standardize    = TRUE,
          relaxed        = relax_flag
        )
      }, warning = function(w) invokeRestart("muffleWarning"))
    }

    m_auto_std   <- fit_svem("auto", FALSE, "SVEM_auto_std")
    m_auto_relax <- fit_svem("auto", TRUE,  "SVEM_auto_relax")
    m_AIC_std    <- fit_svem("wAIC", FALSE, "SVEM_wAIC_std")
    m_AIC_relax  <- fit_svem("wAIC", TRUE,  "SVEM_wAIC_relax")
    m_BIC_std    <- fit_svem("wBIC", FALSE, "SVEM_wBIC_std")
    m_BIC_relax  <- fit_svem("wBIC", TRUE,  "SVEM_wBIC_relax")

    # glmnet_with_cv fits: corrected argument name 'relaxed'
    fit_cv <- function(relax_flag) {
      suppressWarnings(
        glmnet_with_cv(
          formula       = form,
          data          = train,
          glmnet_alpha  = alpha_vec,
          nfolds        = CV_NFOLDS,
          repeats       = CV_REPEATS,
          choose_rule   = CV_CHOICE,
          standardize   = TRUE,
          relaxed       = relax_flag
        )
      )
    }
    fit_cv_std   <- fit_cv(FALSE)
    fit_cv_relax <- fit_cv(TRUE)

    # Predictions on holdout (debiased = FALSE for clean truth comparison)
    preds <- list(
      SVEM_auto_std   = as.numeric(predict(m_auto_std,   hold, debias = FALSE)),
      SVEM_auto_relax = as.numeric(predict(m_auto_relax, hold, debias = FALSE)),
      SVEM_wAIC_std   = as.numeric(predict(m_AIC_std,    hold, debias = FALSE)),
      SVEM_wAIC_relax = as.numeric(predict(m_AIC_relax,  hold, debias = FALSE)),
      SVEM_wBIC_std   = as.numeric(predict(m_BIC_std,    hold, debias = FALSE)),
      SVEM_wBIC_relax = as.numeric(predict(m_BIC_relax,  hold, debias = FALSE)),
      CV_std          = as.numeric(predict_cv(fit_cv_std,   hold, debias = FALSE)),
      CV_relax        = as.numeric(predict_cv(fit_cv_relax, hold, debias = FALSE))
    )

    # objective used (SVEM only). For auto this will be wAIC or wBIC.
    svem_obj_used <- c(
      SVEM_auto_std   = m_auto_std$objective_used,
      SVEM_auto_relax = m_auto_relax$objective_used,
      SVEM_wAIC_std   = m_AIC_std$objective_used,
      SVEM_wAIC_relax = m_AIC_relax$objective_used,
      SVEM_wBIC_std   = m_BIC_std$objective_used,
      SVEM_wBIC_relax = m_BIC_relax$objective_used,
      CV_std          = NA_character_,
      CV_relax        = NA_character_
    )
    relax_flag_map <- c(
      SVEM_auto_std   = FALSE, SVEM_auto_relax = TRUE,
      SVEM_wAIC_std   = FALSE, SVEM_wAIC_relax = TRUE,
      SVEM_wBIC_std   = FALSE, SVEM_wBIC_relax = TRUE,
      CV_std          = FALSE, CV_relax        = TRUE
    )
    method_family_map <- c(
      SVEM_auto_std   = "SVEM", SVEM_auto_relax = "SVEM",
      SVEM_wAIC_std   = "SVEM", SVEM_wAIC_relax = "SVEM",
      SVEM_wBIC_std   = "SVEM", SVEM_wBIC_relax = "SVEM",
      CV_std          = "CV",   CV_relax        = "CV"
    )

    purrr::imap_dfr(preds, function(p, name) {
      m <- metric_one(p, hold_true, sd_hold_true)
      tibble(
        RunID           = factor(run_id),
        n_total         = n_total,
        TheoreticalR2   = factor(r2),
        Holdout_SDTrueY = sd_hold_true,
        Setting         = factor(name, levels = c(
          "SVEM_auto_std","SVEM_auto_relax",
          "SVEM_wAIC_std","SVEM_wAIC_relax",
          "SVEM_wBIC_std","SVEM_wBIC_relax",
          "CV_std","CV_relax"
        )),
        AlphaScenario   = factor(alpha_tag, levels = c("lasso", "mix")),
        NRASE_Holdout   = unname(m["NRASE_Holdout"]),
        NAAE_Holdout    = unname(m["NAAE_Holdout"]),
        ObjectiveUsed   = factor(svem_obj_used[[name]], levels = c("wAIC","wBIC","wGIC","wSSE")),
        RelaxUsed       = relax_flag_map[[name]],
        MethodFamily    = method_family_map[[name]],
        CV_Choice       = CV_CHOICE
      )
    })
  }

  # run both alpha scenarios and bind
  bind_rows(
    fit_scenario(ALPHAS_LASSO, "lasso"),
    fit_scenario(ALPHAS_MIX,   "mix")
  )
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
    group_by(AlphaScenario, Setting) %>%
    summarise(
      runs       = n(),
      mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
      se_NRASE   = se(NRASE_Holdout),
      mean_NAAE  = mean(NAAE_Holdout, na.rm = TRUE),
      se_NAAE    = se(NAAE_Holdout),
      .groups    = "drop"
    ) %>%
    arrange(AlphaScenario, mean_NRASE)
}

win_rate_tbl <- function(d) {
  winners <- d %>%
    group_by(AlphaScenario, RunID, n_total) %>%
    slice_min(NRASE_Holdout, with_ties = TRUE) %>%
    ungroup() %>%
    mutate(flag = 1L)

  d %>%
    select(AlphaScenario, RunID, n_total, Setting) %>%
    left_join(winners %>% select(AlphaScenario, RunID, n_total, Setting, flag),
              by = c("AlphaScenario","RunID","n_total","Setting")) %>%
    group_by(AlphaScenario, Setting) %>%
    summarise(
      wins     = sum(replace(flag, is.na(flag), 0L)),
      runs     = n_distinct(paste(AlphaScenario, RunID, n_total)),
      win_rate = wins / runs,
      .groups  = "drop"
    ) %>%
    arrange(AlphaScenario, desc(win_rate))
}

avg_rank_tbl <- function(d) {
  d %>%
    group_by(AlphaScenario, RunID, n_total) %>%
    mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
    ungroup() %>%
    group_by(AlphaScenario, Setting) %>%
    summarise(
      mean_rank = mean(rk, na.rm = TRUE),
      se_rank   = se(rk),
      .groups   = "drop"
    ) %>%
    arrange(AlphaScenario, mean_rank)
}

paired_compare <- function(d, a, b, label = "", by_alpha = TRUE) {
  base <- d %>%
    filter(Setting %in% c(a, b)) %>%
    select(AlphaScenario, RunID, n_total, Setting, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = Setting, values_from = NRASE_Holdout)

  do_one <- function(tb, subtitle) {
    if (!all(c(a, b) %in% names(tb))) {
      cat("Skipping paired compare (missing columns)\n"); return(invisible(NULL))
    }
    delta <- tb[[a]] - tb[[b]]
    delta <- delta[is.finite(delta)]
    if (length(delta) < 3) {
      cat("Not enough pairs for", subtitle, "\n"); return(invisible(NULL))
    }
    tt <- t.test(delta)
    wt <- suppressWarnings(wilcox.test(delta, exact = FALSE))
    cat(sprintf("\n-- %s %s: %s - %s (N=%d) --\n  meanDelta = %+0.4f   t(%d) = %0.2f  p = %.3g   |   medianDelta = %+0.4f  Wilcoxon p = %.3g\n",
                label, subtitle, a, b, length(delta),
                mean(delta), as.integer(tt$parameter), unname(tt$statistic), tt$p.value,
                stats::median(delta), wt$p.value))
  }

  if (by_alpha) {
    base %>% group_by(AlphaScenario) %>% group_walk(~ do_one(.x, paste0("[", .y$AlphaScenario, "]")))
  } else {
    do_one(base, "")
  }
}

# ----------- PRINT: Overall and by alpha -----------
cat("\n==================== OVERALL (by AlphaScenario) ====================\n")
overall <- df
print(summary_by(overall))

cat("\nWin rate (tie-aware):\n")
print(win_rate_tbl(overall))

cat("\nAverage rank:\n")
print(avg_rank_tbl(overall))

# Key paired comparisons answering the main questions
paired_compare(overall, "SVEM_auto_relax", "SVEM_auto_std", "SVEM: relax vs std")    # does relaxed help SVEM?
paired_compare(overall, "CV_relax",        "CV_std",        "CV: relax vs std")      # does relaxed help CV?
paired_compare(overall, "SVEM_auto_relax", "CV_relax",      "Relax: SVEM(auto) vs CV")
paired_compare(overall, "SVEM_auto_std",   "CV_std",        "Std: SVEM(auto) vs CV")

# Compare AUTO vs fixed objectives (SVEM)
paired_compare(overall, "SVEM_auto_std",   "SVEM_wAIC_std",  "SVEM std: auto vs wAIC")
paired_compare(overall, "SVEM_auto_std",   "SVEM_wBIC_std",  "SVEM std: auto vs wBIC")
paired_compare(overall, "SVEM_auto_relax", "SVEM_wAIC_relax","SVEM relax: auto vs wAIC")
paired_compare(overall, "SVEM_auto_relax", "SVEM_wBIC_relax","SVEM relax: auto vs wBIC")

cat("\n==================== MEAN NRASE vs n_total (faceted by AlphaScenario) ====================\n")
mean_by_n <- df %>%
  group_by(AlphaScenario, n_total, Setting) %>%
  summarise(mean_NRASE = mean(NRASE_Holdout), .groups = "drop")

print(
  ggplot(mean_by_n, aes(n_total, mean_NRASE, color = Setting, group = Setting)) +
    geom_line() + geom_point() +
    facet_wrap(~ AlphaScenario, ncol = 1, scales = "free_y") +
    labs(x = "n_total", y = "Mean NRASE (holdout)", color = "Method")
)

# Cutover: SVEM_auto_std vs CV_std per n_total and alpha scenario
cat("\n==================== CUTOVER: SVEM_auto_std vs CV_std by n_total & alpha ====================\n")
cutover_tbl <- df %>%
  select(AlphaScenario, RunID, n_total, Setting, NRASE_Holdout) %>%
  pivot_wider(names_from = Setting, values_from = NRASE_Holdout) %>%
  mutate(Delta_SVEMstd_vs_CVstd = SVEM_auto_std - CV_std) %>%
  group_by(AlphaScenario, n_total) %>%
  summarise(mean_delta = mean(Delta_SVEMstd_vs_CVstd, na.rm = TRUE), .groups = "drop") %>%
  arrange(AlphaScenario, n_total)
print(cutover_tbl)
cat("\nRows with mean_delta > 0 (CV_std better on average):\n")
print(cutover_tbl %>% filter(mean_delta > 0))

# Relax penalties by n_total and alpha scenario
cat("\n==================== RELAX PENALTY by n_total & alpha ====================\n")
deltas_by_n <- df %>%
  select(AlphaScenario, RunID, n_total, Setting, NRASE_Holdout) %>%
  distinct() %>%
  pivot_wider(names_from = Setting, values_from = NRASE_Holdout) %>%
  mutate(
    Delta_SVEM_relax_vs_std = SVEM_auto_relax - SVEM_auto_std,
    Delta_CV_relax_vs_std   = CV_relax        - CV_std
  )

relax_tbl <- deltas_by_n %>%
  group_by(AlphaScenario, n_total) %>%
  summarise(
    SVEM_relax_penalty = mean(Delta_SVEM_relax_vs_std, na.rm = TRUE),
    CV_relax_penalty   = mean(Delta_CV_relax_vs_std,   na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(AlphaScenario, n_total)
print(relax_tbl, n = Inf)

# Objective mix that "auto" actually chose (SVEM rows only)
cat("\n==================== AUTO OBJECTIVE MIX (SVEM rows only) ====================\n")
obj_mix <- df %>%
  filter(grepl("^SVEM", Setting)) %>%
  count(AlphaScenario, n_total, ObjectiveUsed, name = "n") %>%
  group_by(AlphaScenario, n_total) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  complete(AlphaScenario, n_total, ObjectiveUsed, fill = list(n = 0, pct = 0)) %>%
  arrange(AlphaScenario, n_total, ObjectiveUsed)
print(obj_mix, n = Inf)

print(
  ggplot(obj_mix %>% filter(!is.na(ObjectiveUsed)),
         aes(x = n_total, y = pct, fill = ObjectiveUsed)) +
    geom_col(position = "fill") +
    facet_wrap(~ AlphaScenario, ncol = 1) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "n_total", y = "Objective share", fill = "Auto objective")
)

cat("\n==================== AUTO CHOICES (auto settings only) ====================\n")
auto_only_mix <- df %>%
  filter(Setting %in% c("SVEM_auto_std","SVEM_auto_relax")) %>%
  count(AlphaScenario, n_total, ObjectiveUsed, name = "n") %>%
  group_by(AlphaScenario, n_total) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  complete(AlphaScenario, n_total, ObjectiveUsed, fill = list(n = 0, pct = 0)) %>%
  arrange(AlphaScenario, n_total, ObjectiveUsed)
print(auto_only_mix, n = Inf)

# ==================== NEW: Alpha scenario comparisons (lasso vs mix) ====================
# Paired tests within each method Setting across RunID and n_total.

alpha_compare_table <- function(d) {
  d %>%
    select(Setting, RunID, n_total, AlphaScenario, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = AlphaScenario, values_from = NRASE_Holdout) %>%
    mutate(delta_lasso_minus_mix = lasso - mix) %>%
    group_by(Setting) %>%
    summarise(
      N_pairs    = sum(is.finite(delta_lasso_minus_mix)),
      mean_delta = mean(delta_lasso_minus_mix, na.rm = TRUE),
      median_delta = median(delta_lasso_minus_mix, na.rm = TRUE),
      t_stat     = if (N_pairs > 2) unname(t.test(delta_lasso_minus_mix)$statistic) else NA_real_,
      t_df       = if (N_pairs > 2) as.integer(unname(t.test(delta_lasso_minus_mix)$parameter)) else NA_integer_,
      t_p        = if (N_pairs > 2) unname(t.test(delta_lasso_minus_mix)$p.value) else NA_real_,
      wilcox_p   = if (N_pairs > 2) suppressWarnings(wilcox.test(delta_lasso_minus_mix, exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(mean_delta)
}

alpha_compare_by_n <- function(d) {
  d %>%
    select(Setting, RunID, n_total, AlphaScenario, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = AlphaScenario, values_from = NRASE_Holdout) %>%
    mutate(delta_lasso_minus_mix = lasso - mix) %>%
    group_by(Setting, n_total) %>%
    summarise(
      N_pairs    = sum(is.finite(delta_lasso_minus_mix)),
      mean_delta = mean(delta_lasso_minus_mix, na.rm = TRUE),
      median_delta = median(delta_lasso_minus_mix, na.rm = TRUE),
      t_p        = if (N_pairs > 2) unname(t.test(delta_lasso_minus_mix)$p.value) else NA_real_,
      wilcox_p   = if (N_pairs > 2) suppressWarnings(wilcox.test(delta_lasso_minus_mix, exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(Setting, n_total)
}

alpha_win_rate_tbl <- function(d) {
  d %>%
    select(Setting, RunID, n_total, AlphaScenario, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = AlphaScenario, values_from = NRASE_Holdout) %>%
    mutate(win_mix = as.integer(mix < lasso)) %>%
    group_by(Setting) %>%
    summarise(
      N_pairs = sum(is.finite(win_mix)),
      mix_win_rate = mean(win_mix, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mix_win_rate))
}

cat("\n==================== ALPHA EFFECT (lasso - mix) per Setting ====================\n")
alpha_overall <- alpha_compare_table(df)
print(alpha_overall, n = Inf)

cat("\nWin rate of mix over lasso (per Setting):\n")
print(alpha_win_rate_tbl(df), n = Inf)

cat("\n==================== ALPHA EFFECT BY n_total (lasso - mix) ====================\n")
alpha_by_n <- alpha_compare_by_n(df)
print(alpha_by_n, n = Inf)

# Optional plot: mean delta (lasso - mix) vs n_total by Setting
# Negative means mix is better (lower NRASE).
print(
  ggplot(alpha_by_n, aes(n_total, mean_delta, color = Setting, group = Setting)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() + geom_point() +
    labs(x = "n_total", y = "Mean (lasso - mix) NRASE", color = "Method")
)

cat("\nDone.\n")
