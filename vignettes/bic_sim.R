# --- simulate_svemnet_lasso_with_wSSE.R --------------------------------------
# Adds SVEMnet objective "wSSE" (LASSO) to the previous script.
# Compares: SVEM_wAIC_lasso, SVEM_wBIC_lasso, SVEM_wSSE_lasso vs glmnet_cv_lasso.
# -----------------------------------------------------------------------------

# Packages
pkgs <- c("SVEMnet","glmnet","data.table","dplyr","tibble","tidyr","purrr","stringr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  install.packages(p, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

set.seed(1234)

# ------------------------------- Controls -------------------------------------
OUT_ITERS   <- 50                 # bump up when ready
N_TOTAL_SEQ <- seq(15, 95, by=10) # 15,25,35,45,...,95
HOLDOUT_N   <- 800
R2_LEVELS   <- c(0.3, 0.5, 0.7, 0.9)
SMALL_N_MAX <- 40                 # key threshold for the decision

# --------------------------- Helper generators --------------------------------
rmixture_bounded <- function(n) {
  out <- matrix(NA_real_, nrow = n, ncol = 4); colnames(out) <- c("A","B","C","D")
  i <- 1
  while (i <= n) {
    A <- runif(1, 0.1, 0.4)
    rest <- 1 - A
    u <- rexp(3); u <- u / sum(u) * rest
    B <- u[1]; C <- u[2]; D <- u[3]
    if (B <= 0.8 && C <= 0.8 && D <= 0.8) { out[i,] <- c(A,B,C,D); i <- i+1 }
  }
  as.data.frame(out)
}
space_fill_mixture <- function(n) rmixture_bounded(n)
sample_E <- function(n) sample(c("0","0.002"), n, replace = TRUE)

draw_pv <- function(nparm = 25) {
  lap1 <- function() rexp(1) - rexp(1)
  pv <- numeric(nparm)
  pv[1] <- lap1()
  for (j in 2:5)  pv[j] <- lap1() * rbinom(1,1,0.8)
  for (j in 6:nparm) pv[j] <- lap1() * rbinom(1,1,0.5)
  pv
}

true_y <- function(A,B,C,D,E, pv) {
  s  <- 0.9
  zA <- (A - 0.1)/s; zB <- B/s; zC <- C/s; zD <- D/s
  E_sign <- ifelse(E == "0",  1, -1)
  part1 <- pv[1]*zA + pv[2]*zB + pv[3]*zC + pv[4]*zD +
    E_sign*pv[5]*zA + E_sign*pv[6]*zB + E_sign*pv[7]*zC + E_sign*pv[8]*zD
  part2 <- 4 * ( pv[9]*zA*zB + pv[10]*zA*zC + pv[11]*zA*zD +
                   pv[12]*zB*zC + pv[13]*zB*zD + pv[14]*zC*zD )
  part3 <- 27 * ( pv[15]*zA*zB*zC + pv[16]*zA*zB*zD +
                    pv[17]*zA*zC*zD + pv[18]*zB*zC*zD )
  part4 <- 27 * ( pv[19]*zB*zA*(zA - zB) + pv[20]*zC*zA*(zA - zC) +
                    pv[21]*zC*zB*(zB - zC) + pv[22]*zD*zA*(zA - zD) +
                    pv[23]*zD*zB*(zB - zD) + pv[24]*zD*zC*(zC - zD) )
  part5 <- 256 * pv[25]*zA*zB*zC*zD
  part1 + part2 + part3 + part4 + part5
}

# Feature construction (consistent across train/holdout)
build_X <- function(df) {
  En <- ifelse(df$E == "0.002", 1, 0)
  with(df, {
    s  <- 0.9
    zA <- (A - 0.1)/s; zB <- B/s; zC <- C/s; zD <- D/s
    X <- cbind(
      A,B,C,D, En,
      A_En = A*En, B_En = B*En, C_En = C*En, D_En = D*En,
      A_B = A*B, A_C = A*C, A_D = A*D, B_C = B*C, B_D = B*D, C_D = C*D,
      A_B_C = A*B*C, A_B_D = A*B*D, A_C_D = A*C*D, B_C_D = B*C*D,
      A_B_C_D = A*B*C*D,
      uSC1 = zB*zA*(zA - zB), uSC2 = zC*zA*(zA - zC), uSC3 = zC*zB*(zB - zC),
      uSC4 = zD*zA*(zA - zD), uSC5 = zD*zB*(zB - zD), uSC6 = zD*zC*(zC - zD)
    )
    as.data.frame(X)
  })
}

# Metrics
compute_metrics <- function(y, yhat) {
  sse <- sum((y - yhat)^2); sst <- sum( (y - mean(y))^2 )
  rsq <- if (sst > 0) 1 - sse/sst else NA_real_
  rmse <- sqrt(mean((y - yhat)^2)); mae <- mean(abs(y - yhat))
  list(RSquare = rsq, RASE = rmse, AAE = mae)
}

# Prediction from named coefficient vector (handles (Intercept)/Intercept)
predict_from_parms <- function(parms, Xdf) {
  b <- parms
  b0 <- 0
  int_name <- intersect(names(b), c("(Intercept)","Intercept"))
  if (length(int_name)) { b0 <- as.numeric(b[int_name[1]]); b <- b[setdiff(names(b), int_name[1])] }
  cols_present <- intersect(names(b), colnames(Xdf))
  xb <- as.matrix(Xdf[, cols_present, drop = FALSE]) %*% as.numeric(b[cols_present])
  as.numeric(b0 + xb)
}

# Fitters (LASSO only)
fit_svemnet <- function(df, objective = c("wAIC","wBIC","wSSE")) {
  objective <- match.arg(objective)
  m <- SVEMnet::SVEMnet(Response ~ ., data = df,
                        glmnet_alpha = c(1), nBoot = 300, objective = objective)
  list(parms = m$parms,
       label = paste0("SVEM_", objective, "_lasso"))
}

fit_glmnet_cv <- function(df) {
  g <- SVEMnet::glmnet_with_cv(Response ~ ., data = df, glmnet_alpha = c(1))
  list(parms = g$parms, label = "glmnet_cv_lasso")
}

# ------------------------------- One run --------------------------------------
run_one <- function(run_id, n_total) {
  pv <- draw_pv(25)
  r2 <- sample(R2_LEVELS, 1)

  # Reference SD(TrueY) for noise calibration (analogous to your 15k sample)
  ref <- space_fill_mixture(15000); ref$E <- sample_E(nrow(ref))
  ref$TrueY <- with(ref, true_y(A,B,C,D,E, pv))
  y_sd_ref  <- sd(ref$TrueY)
  error_sd  <- y_sd_ref * sqrt((1 - r2) / r2)

  # Training design
  tr <- space_fill_mixture(n_total); tr$E <- sample_E(nrow(tr))
  tr$TrueY <- with(tr, true_y(A,B,C,D,E, pv))
  tr$Y     <- tr$TrueY + rnorm(nrow(tr), 0, error_sd)

  # Holdout design
  ho <- space_fill_mixture(HOLDOUT_N); ho$E <- sample_E(nrow(ho))
  ho$TrueY <- with(ho, true_y(A,B,C,D,E, pv))

  # Features
  X_tr <- build_X(tr); X_ho <- build_X(ho)
  df_tr <- data.frame(X_tr, Response = tr$Y)

  # Models to fit  (now includes wSSE)
  fits <- list(
    fit_svemnet(df_tr, "wAIC"),
    fit_svemnet(df_tr, "wBIC"),
    fit_svemnet(df_tr, "wSSE"),
    fit_glmnet_cv(df_tr)      # benchmark only
  )

  hold_sd_true <- sd(ho$TrueY)

  rows <- purrr::map_dfr(fits, function(f) {
    yhat_ho <- tryCatch(predict_from_parms(f$parms, X_ho), error = function(e) rep(NA_real_, nrow(X_ho)))
    yhat_tr <- tryCatch(predict_from_parms(f$parms, X_tr), error = function(e) rep(NA_real_, nrow(X_tr)))

    m_ho <- compute_metrics(ho$TrueY, yhat_ho)
    m_tr <- compute_metrics(tr$TrueY, yhat_tr)

    tibble::tibble(
      RunID = run_id,
      n_total = n_total,
      TheoreticalR2 = r2,
      Setting = f$label,
      # Holdout metrics (+ standardized)
      RSq_Holdout   = m_ho$RSquare,
      RASE_Holdout  = m_ho$RASE,
      AAE_Holdout   = m_ho$AAE,
      NRASE_Holdout = m_ho$RASE / hold_sd_true,
      NAAE_Holdout  = m_ho$AAE  / hold_sd_true,
      # Training (optional)
      RSq_Train  = m_tr$RSquare,
      RASE_Train = m_tr$RASE,
      AAE_Train  = m_tr$AAE
    )
  })
  rows
}

# ----------------------------- Run simulation ---------------------------------
all_rows <- list(); k <- 1
for (it in seq_len(OUT_ITERS)) {
  for (n_total in N_TOTAL_SEQ) {
    rid <- sprintf("run%05d", k)
    message("Simulating ", rid, " (n_total=", n_total, ")")
    all_rows[[length(all_rows)+1]] <- run_one(rid, n_total)
    k <- k + 1
  }
}
df <- data.table::rbindlist(all_rows)
df$Setting       <- factor(df$Setting,
                           levels = c("SVEM_wAIC_lasso","SVEM_wBIC_lasso","SVEM_wSSE_lasso","glmnet_cv_lasso"))
df$RunID         <- factor(df$RunID)
df$TheoreticalR2 <- factor(df$TheoreticalR2, levels = sort(unique(df$TheoreticalR2)))
df <- df %>% mutate(n_bucket = ifelse(n_total <= SMALL_N_MAX, "n<=40", "n>40"))

# ------------------------------ Analysis utils --------------------------------
summary_by <- function(dd) {
  dd %>%
    group_by(Setting) %>%
    summarise(
      runs       = n_distinct(RunID),
      mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
      se_NRASE   = sd(NRASE_Holdout,   na.rm = TRUE) / sqrt(n_distinct(RunID)),
      mean_NAAE  = mean(NAAE_Holdout,  na.rm = TRUE),
      se_NAAE    = sd(NAAE_Holdout,    na.rm = TRUE) / sqrt(n_distinct(RunID)),
      .groups = "drop"
    ) %>% arrange(mean_NRASE)
}

win_rate_tbl <- function(dd) {
  winners <- dd %>%
    group_by(RunID) %>%
    slice_min(NRASE_Holdout, with_ties = TRUE) %>%
    ungroup() %>% mutate(flag = 1L)

  dd %>%
    select(RunID, Setting) %>%
    left_join(winners %>% select(RunID, Setting, flag), by = c("RunID","Setting")) %>%
    group_by(Setting) %>%
    summarise(
      wins     = sum(replace(flag, is.na(flag), 0L)),
      runs     = n_distinct(RunID),
      win_rate = wins / runs,
      .groups  = "drop"
    ) %>% arrange(desc(win_rate))
}

avg_rank_tbl <- function(dd) {
  dd %>%
    group_by(RunID) %>%
    mutate(rank = rank(NRASE_Holdout, ties.method = "average")) %>%
    ungroup() %>%
    group_by(Setting) %>%
    summarise(
      mean_rank = mean(rank, na.rm = TRUE),
      se_rank   = sd(rank,  na.rm = TRUE)/sqrt(n_distinct(RunID)),
      .groups = "drop"
    ) %>% arrange(mean_rank)
}

paired_compare <- function(dd, s1, s2, label="Overall") {
  wide <- dd %>%
    select(RunID, Setting, NRASE_Holdout) %>%
    tidyr::pivot_wider(names_from = Setting, values_from = NRASE_Holdout)
  if (!all(c(s1,s2) %in% colnames(wide))) return(invisible(NULL))
  d <- wide[[s1]] - wide[[s2]]
  d <- d[is.finite(d)]
  if (length(d) <= 2) return(invisible(NULL))
  tt <- t.test(d)
  ww <- suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE))
  cat(sprintf("\n-- %s: %s - %s (N=%d) --\n", label, s1, s2, length(d)))
  cat(sprintf("  meanΔ = %+0.4f   t(%d) = %0.2f  p = %.3g   |   medianΔ = %+0.4f  Wilcoxon p = %.3g\n",
              mean(d), tt$parameter, tt$statistic, tt$p.value, median(d), ww$p.value))
}

pct_vs_baseline <- function(dd, baseline="glmnet_cv_lasso") {
  wide <- dd %>%
    select(RunID, Setting, NRASE_Holdout) %>%
    tidyr::pivot_wider(names_from = Setting, values_from = NRASE_Holdout) %>%
    tibble::as_tibble()
  others <- setdiff(colnames(wide), c("RunID", baseline))
  if (!baseline %in% colnames(wide)) return(tibble())
  out <- lapply(others, function(s) {
    d <- 100 * (wide[[s]] - wide[[baseline]]) / wide[[baseline]]
    tibble(Setting = s,
           mean_pct = mean(d, na.rm = TRUE),
           se_pct   = sd(d,   na.rm = TRUE)/sqrt(sum(is.finite(d))))
  }) %>% bind_rows() %>% arrange(mean_pct)
  out
}

# ------------------------------- Printouts ------------------------------------
cat("\n==================== OVERALL ====================\n")
overall <- df
print(summary_by(overall))
cat("\nWin rate (tie-aware):\n"); print(win_rate_tbl(overall))
cat("\nAverage rank:\n"); print(avg_rank_tbl(overall))

# Primary decision remains AIC vs BIC (wSSE will appear in tables above)
paired_compare(overall, "SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "Overall (AIC vs BIC)")
cat("\nPercent Δ vs benchmark (glmnet_cv_lasso):\n"); print(pct_vs_baseline(overall))

cat("\n================== STRATIFIED ===================\n")
small <- df %>% filter(n_bucket == "n<=40")
large <- df %>% filter(n_bucket == "n>40")

cat("\n--- n_total <= 40 (primary interest) ---\n")
print(summary_by(small))
cat("\nWin rate:\n"); print(win_rate_tbl(small))
cat("\nAverage rank:\n"); print(avg_rank_tbl(small))
paired_compare(small, "SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "n<=40 (AIC vs BIC)")
cat("\nPercent Δ vs benchmark:\n"); print(pct_vs_baseline(small))

cat("\n--- n_total > 40 ---\n")
print(summary_by(large))
cat("\nWin rate:\n"); print(win_rate_tbl(large))
cat("\nAverage rank:\n"); print(avg_rank_tbl(large))
paired_compare(large, "SVEM_wAIC_lasso", "SVEM_wBIC_lasso", "n>40 (AIC vs BIC)")
cat("\nPercent Δ vs benchmark:\n"); print(pct_vs_baseline(large))

cat("\n================= RECOMMENDATION =================\n")
cat("* Use results in 'n_total <= 40' block to pick SVEM default (focus: AIC vs BIC).\n")
cat("* Benchmark (glmnet_cv_lasso) and wSSE are for additional context.\n")
cat("\nDone.\n")
