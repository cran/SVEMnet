###############################################################################
# Focused simulation: AUTO only, compare five target combos
# Combos:
#   1) SVEM(auto, alpha=c(1,.5),  relaxed=TRUE)   -> SVEM_auto_relax_mix
#   2) SVEM(auto, alpha=c(1,.5),  relaxed=FALSE)  -> SVEM_auto_std_mix
#   3) SVEM(auto, alpha=1,        relaxed=FALSE)  -> SVEM_auto_std_lasso
#   4) CV   (alpha=c(1,.5),       relaxed=FALSE)  -> CV_std_mix
#   5) CV   (alpha=1,             relaxed=FALSE)  -> CV_std_lasso
#
# We run *both* alpha scenarios (lasso=1, mix=c(0.5,1)) for each run, then
# carve out the 5 target combos via a FocusSetting column.
###############################################################################

suppressPackageStartupMessages({
  library(SVEMnet)     # SVEMnet(), glmnet_with_cv(), predict_cv()
  library(glmnet)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(furrr)
  library(ggplot2)
  library(scales)
})

# ------------------- knobs -------------------
CORES        <- 10                     # parallel workers
OUT_ITERS    <- 100                    # repeats per n_total
N_TOTAL_SEQ  <- seq(15, 75, by = 10)   # train sizes to sweep
N_BOOT       <- 300                    # SVEM bootstrap reps
HOLDOUT_N    <- 10000                  # holdout size per run
DIST_N       <- 10000                  # for estimating sd(TrueY) to set noise

ALPHAS_LASSO <- 1
ALPHAS_MIX   <- c(0.5, 1)              # IMPORTANT: no alpha=0

# glmnet_with_cv knobs
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

# Randomized sparsity with Beta priors (tunable concentration)
gen_pv <- function(beta5_8 = c(8, 2),   # mean ~0.80
                   beta9_25 = c(5, 5))  # mean ~0.50
{
  stopifnot(length(beta5_8)==2, length(beta9_25)==2)
  rexp_signed <- function(n) stats::rexp(n) - stats::rexp(n)

  p5_8  <- stats::rbeta(1, beta5_8[1], beta5_8[2])
  p9_25 <- stats::rbeta(1, beta9_25[1], beta9_25[2])

  p1_4   <- rexp_signed(4)
  p5_8v  <- rexp_signed(4)  * rbinom(4, 1, p5_8)
  p9_25v <- rexp_signed(17) * rbinom(17, 1, p9_25)

  list(pv = c(p1_4, p5_8v, p9_25v), p5_8 = p5_8, p9_25 = p9_25, scheme = "beta")
}

# True response surface
true_response <- function(df, pv) {
  if (is.list(pv)) pv <- pv$pv
  stopifnot(is.numeric(pv), length(pv) == 25)
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

# ----------- One simulated dataset (runs both alpha scenarios) -----------
run_one <- function(run_id, n_total) {

  # 1) Random coefficients and theoretical R^2
  pv_info <- gen_pv()
  pv      <- pv_info$pv
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

    # SVEM: AUTO only, relaxed FALSE and TRUE
    fit_svem <- function(relax_flag) {
      withCallingHandlers({
        SVEMnet::SVEMnet(
          formula        = form,
          data           = train,
          glmnet_alpha   = alpha_vec,      # {1} or {0.5,1}
          nBoot          = N_BOOT,
          objective      = "auto",         # AUTO only
          weight_scheme  = "SVEM",
          standardize    = TRUE,
          relaxed        = relax_flag
        )
      }, warning = function(w) invokeRestart("muffleWarning"))
    }

    # CV: relaxed = FALSE only
    fit_cv <- function() {
      suppressWarnings(
        glmnet_with_cv(
          formula       = form,
          data          = train,
          glmnet_alpha  = alpha_vec,
          nfolds        = CV_NFOLDS,
          repeats       = CV_REPEATS,
          choose_rule   = CV_CHOICE,
          standardize   = TRUE,
          relaxed       = FALSE
        )
      )
    }

    # Fits (names DO NOT include alpha; alpha lives in AlphaScenario)
    m_svem_std   <- fit_svem(FALSE)
    m_svem_relax <- fit_svem(TRUE)
    m_cv_std     <- fit_cv()

    # Predictions on holdout (debiased = FALSE for clean truth comparison)
    preds <- list(
      SVEM_auto_std   = as.numeric(predict(m_svem_std,   hold, debias = FALSE)),
      SVEM_auto_relax = as.numeric(predict(m_svem_relax, hold, debias = FALSE)),
      CV_std          = as.numeric(predict_cv(m_cv_std,  hold, debias = FALSE))
    )

    svem_obj_used <- c(
      SVEM_auto_std   = m_svem_std$objective_used,
      SVEM_auto_relax = m_svem_relax$objective_used,
      CV_std          = NA_character_
    )
    relax_flag_map    <- c(SVEM_auto_std = FALSE, SVEM_auto_relax = TRUE, CV_std = FALSE)
    method_family_map <- c(SVEM_auto_std = "SVEM", SVEM_auto_relax = "SVEM", CV_std = "CV")

    purrr::imap_dfr(preds, function(p, name) {
      m <- metric_one(p, hold_true, sd_hold_true)
      tibble(
        RunID           = factor(run_id),
        n_total         = n_total,
        TheoreticalR2   = factor(r2),
        Holdout_SDTrueY = sd_hold_true,
        Setting         = factor(name, levels = c("SVEM_auto_std","SVEM_auto_relax","CV_std")),
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

# ---------- robust worker setup + safe runner ----------
furrr_pkgs <- c("SVEMnet","glmnet","dplyr","tidyr","tibble","purrr","stats")

safe_run_one <- function(run_id, n_total) {
  tryCatch(
    {
      suppressPackageStartupMessages({
        library(SVEMnet); library(glmnet); library(dplyr)
        library(tidyr);   library(tibble); library(purrr)
      })
      res <- run_one(run_id, n_total)
      stopifnot(is.data.frame(res), all(c("Setting","NRASE_Holdout") %in% names(res)))
      res
    },
    error = function(e) {
      tibble(
        RunID           = factor(run_id),
        n_total         = n_total,
        TheoreticalR2   = NA,
        Holdout_SDTrueY = NA_real_,
        Setting         = factor("__ERROR__", levels = "__ERROR__"),
        AlphaScenario   = factor(NA_character_, levels = c("lasso","mix")),
        NRASE_Holdout   = NA_real_,
        NAAE_Holdout    = NA_real_,
        ObjectiveUsed   = factor(NA_character_, levels = c("wAIC","wBIC","wGIC","wSSE")),
        RelaxUsed       = NA,
        MethodFamily    = NA_character_,
        CV_Choice       = CV_CHOICE,
        ErrorMsg        = conditionMessage(e)
      )
    },
    warning = function(w) invokeRestart("muffleWarning")
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
  .f = ~ safe_run_one(..1, ..2),
  .options = furrr::furrr_options(seed = TRUE, packages = furrr_pkgs)
)

df <- dplyr::bind_rows(df_list)

# If anything failed, print a quick summary:
if (any(df$Setting == "__ERROR__")) {
  message("Some jobs returned errors; top messages:")
  print(
    df %>%
      dplyr::filter(Setting == "__ERROR__") %>%
      dplyr::count(ErrorMsg, sort = TRUE) %>%
      dplyr::slice_head(n = 10)
  )
}

df$Setting <- droplevels(df$Setting)

# ------------------- Focus on the 5 requested combos -------------------
# Build FocusSetting using Setting + AlphaScenario
df <- df %>%
  mutate(
    FocusSetting = dplyr::case_when(
      Setting == "SVEM_auto_relax" & AlphaScenario == "mix"   ~ "SVEM_auto_relax_mix",
      Setting == "SVEM_auto_std"   & AlphaScenario == "mix"   ~ "SVEM_auto_std_mix",
      Setting == "SVEM_auto_std"   & AlphaScenario == "lasso" ~ "SVEM_auto_std_lasso",
      Setting == "CV_std"          & AlphaScenario == "mix"   ~ "CV_std_mix",
      Setting == "CV_std"          & AlphaScenario == "lasso" ~ "CV_std_lasso",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(FocusSetting), Setting != "__ERROR__") %>%
  mutate(FocusSetting = factor(FocusSetting,
                               levels = c("SVEM_auto_relax_mix",
                                          "SVEM_auto_std_mix",
                                          "SVEM_auto_std_lasso",
                                          "CV_std_mix",
                                          "CV_std_lasso")))

# ----------- Helper summaries -----------
se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

summary_by <- function(d) {
  d %>%
    group_by(n_total, FocusSetting) %>%
    summarise(
      runs       = n(),
      mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
      se_NRASE   = se(NRASE_Holdout),
      mean_NAAE  = mean(NAAE_Holdout, na.rm = TRUE),
      se_NAAE    = se(NAAE_Holdout),
      .groups    = "drop"
    ) %>%
    arrange(n_total, mean_NRASE)
}

win_rate_tbl <- function(d) {
  winners <- d %>%
    group_by(RunID, n_total) %>%
    slice_min(NRASE_Holdout, with_ties = TRUE) %>%
    ungroup() %>%
    mutate(flag = 1L)

  d %>%
    select(RunID, n_total, FocusSetting) %>%
    left_join(winners %>% select(RunID, n_total, FocusSetting, flag),
              by = c("RunID","n_total","FocusSetting")) %>%
    group_by(FocusSetting) %>%
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
    group_by(FocusSetting) %>%
    summarise(
      mean_rank = mean(rk, na.rm = TRUE),
      se_rank   = se(rk),
      .groups   = "drop"
    ) %>%
    arrange(mean_rank)
}

paired_compare <- function(d, a, b, label = "") {
  base <- d %>%
    filter(FocusSetting %in% c(a, b)) %>%
    select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
    distinct() %>%
    pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout)

  if (!all(c(a, b) %in% names(base))) {
    cat("Skipping paired compare (missing columns)\n"); return(invisible(NULL))
  }
  delta <- base[[a]] - base[[b]]
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

# ----------- PRINT: Rankings & key head-to-heads -----------
cat("\n==================== SUMMARY by n_total ====================\n")
print(summary_by(df), n = Inf)

cat("\n==================== OVERALL WIN RATE ====================\n")
print(win_rate_tbl(df))

cat("\n==================== AVERAGE RANK ====================\n")
print(avg_rank_tbl(df))

# Head-to-heads matching your questions
paired_compare(df, "SVEM_auto_relax_mix", "SVEM_auto_std_mix", "SVEM: relax (mix) vs std (mix)")
paired_compare(df, "SVEM_auto_std_lasso", "SVEM_auto_std_mix", "SVEM std: lasso vs mix")
paired_compare(df, "CV_std_lasso",       "CV_std_mix",        "CV std: lasso vs mix")
paired_compare(df, "SVEM_auto_std_mix",  "CV_std_mix",        "SVEM std mix vs CV std mix")
paired_compare(df, "SVEM_auto_relax_mix","CV_std_mix",        "SVEM relax mix vs CV std mix")

# Auto objective mix (SVEM rows only)
cat("\n==================== AUTO OBJECTIVE MIX (SVEM only) ====================\n")
obj_mix <- df %>%
  filter(grepl("^SVEM", FocusSetting)) %>%
  count(n_total, FocusSetting, ObjectiveUsed, name = "n") %>%
  group_by(n_total, FocusSetting) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  complete(n_total, FocusSetting, ObjectiveUsed, fill = list(n = 0, pct = 0)) %>%
  arrange(n_total, FocusSetting, ObjectiveUsed)
print(obj_mix, n = Inf)

print(
  ggplot(obj_mix %>% filter(!is.na(ObjectiveUsed)),
         aes(x = n_total, y = pct, fill = ObjectiveUsed)) +
    geom_col(position = "fill") +
    facet_wrap(~ FocusSetting, ncol = 1) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "n_total", y = "Objective share", fill = "Auto objective")
)

cat("\nDone.\n")



















analyze_sim <- function(df, seed = 123, bt_ref = NULL, B_boot = 2000, B_break = 400) {
  stopifnot(all(c("RunID","n_total","FocusSetting","TheoreticalR2","NRASE_Holdout") %in% names(df)))
  set.seed(seed)

  # ---------- helpers ----------
  se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  eps_win_tbl <- function(d, eps = 0.01) {
    d %>%
      dplyr::group_by(RunID, n_total) %>%
      dplyr::mutate(best = min(NRASE_Holdout, na.rm = TRUE),
                    regret = NRASE_Holdout - best,
                    within = regret <= eps) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(FocusSetting) %>%
      dplyr::summarise(N = dplyr::n_distinct(paste(RunID, n_total)),
                       pct_within = mean(within, na.rm = TRUE),
                       mean_regret = mean(regret, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(pct_within))
  }
  boot_ci <- function(x, stat = mean, B = B_boot) {
    if (length(x) == 0 || all(!is.finite(x))) return(setNames(rep(NA_real_,3), c("2.5%","50%","97.5%")))
    s <- replicate(B, stat(sample(x[is.finite(x)], length(x[is.finite(x)]), replace = TRUE)))
    stats::quantile(s, c(`2.5%`=.025, `50%`=.5, `97.5%`=.975), na.rm = TRUE)
  }
  boot_delta_ci <- function(d, a, b, B = B_boot) {
    wide <- d %>%
      dplyr::filter(FocusSetting %in% c(a,b)) %>%
      dplyr::select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
      dplyr::mutate(delta = .data[[a]] - .data[[b]]) %>%
      dplyr::pull(delta) %>% stats::na.omit()
    qs <- boot_ci(wide, stat = mean, B = B)
    tibble::tibble(pair = paste(a,"-",b),
                   N = length(wide),
                   mean = mean(wide),
                   lo = unname(qs["2.5%"]),
                   med = unname(qs["50%"]),
                   hi = unname(qs["97.5%"]))
  }
  win_rate_tbl <- function(d) {
    winners <- d %>%
      dplyr::group_by(RunID, n_total) %>%
      dplyr::slice_min(NRASE_Holdout, with_ties = TRUE) %>%
      dplyr::ungroup() %>% dplyr::mutate(flag = 1L)
    d %>%
      dplyr::select(RunID, n_total, FocusSetting) %>%
      dplyr::left_join(winners %>% dplyr::select(RunID, n_total, FocusSetting, flag),
                       by = c("RunID","n_total","FocusSetting")) %>%
      dplyr::group_by(FocusSetting) %>%
      dplyr::summarise(wins     = sum(replace(flag, is.na(flag), 0L)),
                       runs     = dplyr::n_distinct(paste(RunID, n_total)),
                       win_rate = wins / runs,
                       .groups  = "drop") %>%
      dplyr::arrange(dplyr::desc(win_rate))
  }
  avg_rank_tbl <- function(d) {
    d %>%
      dplyr::group_by(RunID, n_total) %>%
      dplyr::mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(FocusSetting) %>%
      dplyr::summarise(mean_rank = mean(rk, na.rm = TRUE),
                       se_rank   = se(rk),
                       .groups   = "drop") %>%
      dplyr::arrange(mean_rank)
  }
  pairwise_theta_matrix <- function(d) {
    meths <- levels(d$FocusSetting) %||% unique(d$FocusSetting)
    wide <- d %>%
      dplyr::select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout)
    out <- matrix(NA_real_, nrow = length(meths), ncol = length(meths),
                  dimnames = list(meths, meths))
    for (i in meths) for (j in meths) {
      if (i == j || !all(c(i,j) %in% names(wide))) next
      dij <- wide[[i]] - wide[[j]]
      out[i,j] <- mean(dij < 0, na.rm = TRUE) + 0.5 * mean(dij == 0, na.rm = TRUE)
    }
    out
  }
  bt_skill_glm <- function(d, bt_ref = NULL) {
    # Winners per block (break ties deterministically by method name)
    winners <- d %>%
      dplyr::group_by(RunID, n_total) %>%
      dplyr::arrange(NRASE_Holdout, FocusSetting, .by_group = TRUE) %>%
      dplyr::slice(1) %>% dplyr::ungroup() %>%
      dplyr::transmute(RunID, n_total, winner = FocusSetting)
    rest <- d %>%
      dplyr::select(RunID, n_total, loser = FocusSetting) %>%
      dplyr::anti_join(winners, by = c("RunID","n_total","loser"="winner"))
    pairs <- dplyr::left_join(winners, rest, by = c("RunID","n_total"))
    meths <- levels(d$FocusSetting) %||% unique(d$FocusSetting)
    if (is.null(bt_ref)) {
      wr <- win_rate_tbl(d); bt_ref <- wr$FocusSetting[1]
    }
    winner_mm <- stats::model.matrix(~ winner - 1, data = pairs)
    loser_mm  <- stats::model.matrix(~ loser  - 1, data = pairs)
    colnames(winner_mm) <- sub("^winner", "", colnames(winner_mm))
    colnames(loser_mm)  <- sub("^loser",  "", colnames(loser_mm))
    Xdiff <- winner_mm[, meths, drop = FALSE] - loser_mm[, meths, drop = FALSE]
    keep_cols <- setdiff(colnames(Xdiff), bt_ref) # drop ref for identifiability
    y <- rep(1, nrow(Xdiff))
    fit <- stats::glm(y ~ Xdiff[, keep_cols, drop = FALSE] - 1, family = stats::binomial())
    coefs <- rep(NA_real_, length(meths)); names(coefs) <- meths
    coefs[bt_ref] <- 0
    coefs[keep_cols] <- stats::coef(fit)
    tibble::tibble(method = names(coefs),
                   logit_skill = unname(coefs)) %>%
      dplyr::mutate(skill_odds = exp(logit_skill - max(logit_skill))) %>%
      dplyr::arrange(dplyr::desc(skill_odds))
  }
  break_even_boot <- function(d, a = "SVEM_auto_relax_mix", b = "SVEM_auto_std_mix", B = B_break) {
    base <- d %>%
      dplyr::filter(FocusSetting %in% c(a,b)) %>%
      dplyr::select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
      dplyr::mutate(delta = .data[[a]] - .data[[b]])
    # spline fit
    m <- stats::lm(delta ~ splines::ns(n_total, df = 3), data = base)
    grid <- tibble::tibble(n_total = seq(min(base$n_total), max(base$n_total), by = 1),
                           delta_hat = stats::predict(m, newdata = tibble::tibble(n_total = seq(min(base$n_total), max(base$n_total), 1))))
    est <- grid$n_total[which.min(abs(grid$delta_hat))]
    boots <- replicate(B, {
      idx <- sample(seq_len(nrow(base)), replace = TRUE)
      mb  <- try(stats::lm(delta ~ splines::ns(n_total, df = 3), data = base[idx, ]), silent = TRUE)
      if (inherits(mb, "try-error")) return(NA_real_)
      gh  <- try(stats::predict(mb, newdata = grid), silent = TRUE)
      if (inherits(gh, "try-error")) return(NA_real_)
      grid$n_total[which.min(abs(gh))]
    })
    list(est = est,
         ci  = stats::quantile(boots, c(.025,.5,.975), na.rm = TRUE))
  }

  # ---------- PRINT: summaries ----------
  cat("\n==================== SUMMARY by n_total ====================\n")
  sum_by_n <- df %>%
    dplyr::group_by(n_total, FocusSetting) %>%
    dplyr::summarise(runs       = dplyr::n(),
                     mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
                     se_NRASE   = se(NRASE_Holdout),
                     mean_NAAE  = mean(NAAE_Holdout, na.rm = TRUE),
                     se_NAAE    = se(NAAE_Holdout),
                     .groups    = "drop") %>%
    dplyr::arrange(n_total, mean_NRASE)
  print(sum_by_n, n = Inf)

  cat("\n==================== OVERALL WIN RATE ====================\n")
  win_tbl <- win_rate_tbl(df)
  print(win_tbl)

  cat("\n==================== AVERAGE RANK ====================\n")
  avg_rank <- avg_rank_tbl(df)
  print(avg_rank)

  # ---------- Oracle gap / regret ----------
  cat("\n==================== ORACLE GAP (regret to blockwise best) ====================\n")
  reg_df <- df %>%
    dplyr::group_by(RunID, n_total) %>%
    dplyr::mutate(best = min(NRASE_Holdout, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(regret = NRASE_Holdout - best)
  cat("\n-- Per n_total --\n")
  regret_by_n <- reg_df %>%
    dplyr::group_by(n_total, FocusSetting) %>%
    dplyr::summarise(mean_regret = mean(regret, na.rm = TRUE),
                     se = se(regret), .groups="drop")
  print(regret_by_n, n = Inf)

  cat("\n-- Per TheoreticalR2 --\n")
  regret_by_r2 <- reg_df %>%
    dplyr::group_by(TheoreticalR2, FocusSetting) %>%
    dplyr::summarise(mean_regret = mean(regret, na.rm = TRUE),
                     se = se(regret), .groups="drop")
  print(regret_by_r2, n = Inf)

  # ---------- Pairwise P(i<j) theta ----------
  cat("\n==================== PAIRWISE P(i beats j) MATRIX (Mann–Whitney θ) ====================\n")
  theta_mat <- pairwise_theta_matrix(df)
  print(theta_mat)

  # ---------- Bradley–Terry skill ----------
  cat("\n==================== Bradley–Terry skill (GLM, block winners vs rest) ====================\n")
  bt_tbl <- bt_skill_glm(df, bt_ref = bt_ref)
  print(bt_tbl)

  # ---------- AUTO objective choice: counts & proportions ----------
  if ("ObjectiveUsed" %in% names(df)) {
    cat("\n==================== AUTO objective by n_total & R^2 ====================\n")
    auto_tab <- df %>%
      dplyr::filter(grepl("^SVEM", FocusSetting)) %>%
      dplyr::count(n_total, TheoreticalR2, ObjectiveUsed) %>%
      dplyr::group_by(n_total, TheoreticalR2) %>%
      dplyr::mutate(p = n/sum(n)) %>% dplyr::ungroup() %>%
      dplyr::arrange(n_total, TheoreticalR2, ObjectiveUsed)
    print(auto_tab, n = Inf)
  } else {
    auto_tab <- tibble::tibble()
  }

  # ---------- SVEM Δ(relax - std) analyses ----------
  if (all(c("SVEM_auto_relax_mix","SVEM_auto_std_mix") %in% df$FocusSetting)) {
    cat("\n==================== SVEM Δ(relax - std) on mix: meta-regression ====================\n")
    delta_relax <- df %>%
      dplyr::filter(FocusSetting %in% c("SVEM_auto_relax_mix","SVEM_auto_std_mix")) %>%
      dplyr::select(TheoreticalR2, n_total, RunID, FocusSetting, NRASE_Holdout) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
      dplyr::mutate(delta = SVEM_auto_relax_mix - SVEM_auto_std_mix)
    print(summary(stats::lm(delta ~ n_total * TheoreticalR2, data = delta_relax)))

    # per-n head-to-head
    cat("\n-- Per n_total head-to-head (relax - std) --\n")
    svem_relax_vs_std_by_n <- delta_relax %>%
      dplyr::group_by(n_total) %>%
      dplyr::summarise(
        N = sum(is.finite(delta)),
        mean_delta   = mean(delta, na.rm = TRUE),
        median_delta = stats::median(delta, na.rm = TRUE),
        t_p          = if (N > 2) unname(stats::t.test(delta)$p.value) else NA_real_,
        wilcox_p     = if (N > 2) suppressWarnings(stats::wilcox.test(delta, exact = FALSE)$p.value) else NA_real_,
        .groups = "drop"
      )
    print(svem_relax_vs_std_by_n, n = Inf)

    # per-R2 head-to-head
    cat("\n-- Per R^2 head-to-head (relax - std) --\n")
    svem_relax_vs_std_by_r2 <- delta_relax %>%
      dplyr::group_by(TheoreticalR2) %>%
      dplyr::summarise(
        N = sum(is.finite(delta)),
        mean_delta   = mean(delta, na.rm = TRUE),
        median_delta = stats::median(delta, na.rm = TRUE),
        t_p          = if (N > 2) unname(stats::t.test(delta)$p.value) else NA_real_,
        wilcox_p     = if (N > 2) suppressWarnings(stats::wilcox.test(delta, exact = FALSE)$p.value) else NA_real_,
        .groups = "drop"
      )
    print(svem_relax_vs_std_by_r2, n = Inf)
  } else {
    delta_relax <- tibble::tibble()
    svem_relax_vs_std_by_n <- svem_relax_vs_std_by_r2 <- tibble::tibble()
  }

  # ---------- Learning-curve fit: a/sqrt(n) + b ----------
  cat("\n==================== Learning-curve fit (a/sqrt(n) + b) ====================\n")
  lc_tbl <- df %>%
    dplyr::group_by(FocusSetting, n_total) %>%
    dplyr::summarise(m = mean(NRASE_Holdout), .groups="drop") %>%
    dplyr::group_by(FocusSetting) %>%
    dplyr::summarise(
      a  = stats::coef(stats::lm(m ~ I(1/sqrt(n_total))))[2],
      b  = stats::coef(stats::lm(m ~ I(1/sqrt(n_total))))[1],
      r2 = summary(stats::lm(m ~ I(1/sqrt(n_total))))$adj.r.squared,
      .groups="drop"
    )
  print(lc_tbl)

  # ---------- K-means on profiles (mean NRASE per method at each grid) ----------
  cat("\n==================== K-means clusters of performance profiles ====================\n")
  prof_grid <- df %>%
    dplyr::group_by(n_total, TheoreticalR2, FocusSetting) %>%
    dplyr::summarise(m = mean(NRASE_Holdout), .groups = "drop")
  prof_mat <- prof_grid %>% tidyr::pivot_wider(names_from = FocusSetting, values_from = m)
  km <- stats::kmeans(stats::na.omit(prof_mat %>% dplyr::select(-n_total,-TheoreticalR2)), centers = 3, nstart = 10)
  cl_tbl <- dplyr::bind_cols(prof_mat %>% dplyr::select(n_total, TheoreticalR2), cluster = km$cluster) %>%
    dplyr::count(cluster, n_total, TheoreticalR2)
  print(cl_tbl, n = Inf)

  # ---------- Distribution quantiles ----------
  cat("\n==================== Quantiles (NRASE) by n_total ====================\n")
  quant_tbl <- df %>%
    dplyr::group_by(n_total, FocusSetting) %>%
    dplyr::summarise(q25 = stats::quantile(NRASE_Holdout, .25, na.rm = TRUE),
                     med = stats::median(NRASE_Holdout, na.rm = TRUE),
                     q75 = stats::quantile(NRASE_Holdout, .75, na.rm = TRUE),
                     .groups="drop")
  print(quant_tbl, n = Inf)

  # ---------- Break-even n (relax vs std on mix) ----------
  cat("\n==================== Break-even n for SVEM (relax mix vs std mix) ====================\n")
  be <- break_even_boot(df)
  cat(sprintf("Point estimate (nearest n_total): %s \n", be$est))
  print(be$ci)

  # ---------- AUTO vs Oracle(AIC,BIC) gap placeholder ----------
  cat("\n==================== AUTO vs Oracle{AIC,BIC} gap (if available) ====================\n")
  if ("ObjectiveUsed" %in% names(df)) {
    # With AUTO-only sims, AIC/BIC oracles aren't present; keep as placeholder zeroes:
    print(tibble::tibble(mean_gap = 0, q95 = 0))
  } else {
    print(tibble::tibble())
  }

  # ---------- ε-within-best & regret ----------
  cat("\n==================== ε-within-best & regret (per method) ====================\n")
  eps_tbl <- dplyr::bind_rows(
    eps_win_tbl(df, 0.005) %>% dplyr::mutate(eps = 0.005),
    eps_win_tbl(df, 0.01)  %>% dplyr::mutate(eps = 0.01),
    eps_win_tbl(df, 0.02)  %>% dplyr::mutate(eps = 0.02)
  )
  print(eps_tbl, n = Inf)

  # Win-rate bootstrap CIs
  cat("\n-- Win rate bootstrap CIs --\n")
  winners <- df %>% dplyr::group_by(RunID, n_total) %>% dplyr::slice_min(NRASE_Holdout, with_ties = TRUE) %>% dplyr::ungroup()
  win_ci <- df %>%
    dplyr::select(RunID, n_total, FocusSetting) %>%
    dplyr::left_join(winners %>% dplyr::mutate(win=1L), by=c("RunID","n_total","FocusSetting")) %>%
    dplyr::mutate(win = replace(win, is.na(win), 0L)) %>%
    dplyr::group_by(FocusSetting) %>%
    dplyr::summarise(
      win_rate = mean(win),
      N        = dplyr::n_distinct(paste(RunID,n_total)),
      { qs <- boot_ci(win); tibble::tibble(lo=unname(qs["2.5%"]),
                                           med=unname(qs["50%"]),
                                           hi=unname(qs["97.5%"])) },
      .groups = "drop"
    ) %>% dplyr::arrange(dplyr::desc(win_rate))
  print(win_ci)

  # Key delta CI: SVEM_auto_relax_mix - CV_std_mix
  cat("\n-- Key delta CI: SVEM_auto_relax_mix - CV_std_mix --\n")
  key_ci <- boot_delta_ci(df, "SVEM_auto_relax_mix", "CV_std_mix")
  print(key_ci)

  # ---------- Plots ----------
  cat("\n==================== Plot: Rank heatmap (avg rank by n and R^2) ====================\n")
  rank_grid <- df %>%
    dplyr::group_by(TheoreticalR2, n_total, RunID) %>%
    dplyr::mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(TheoreticalR2, n_total, FocusSetting) %>%
    dplyr::summarise(mean_rank = mean(rk), .groups = "drop")
  print(
    ggplot2::ggplot(rank_grid, ggplot2::aes(n_total, TheoreticalR2, fill = mean_rank)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ FocusSetting) +
      ggplot2::scale_fill_viridis_c(direction = -1) +
      ggplot2::labs(x="n_total", y="Theoretical R²", fill="Avg rank")
  )

  cat("\n==================== Plot: ECDF of Regret ====================\n")
  print(
    ggplot2::ggplot(reg_df, ggplot2::aes(NRASE_Holdout - best, colour = FocusSetting)) +
      ggplot2::stat_ecdf(na.rm = TRUE) +
      ggplot2::coord_cartesian(xlim = c(0, stats::quantile(reg_df$NRASE_Holdout - reg_df$best, .95, na.rm = TRUE))) +
      ggplot2::labs(x = "Regret (NRASE − best in block)", y = "ECDF")
  )

  # ---------- Effect sizes & rank model ----------
  cat("\n==================== Cohen's d by n_total (SVEM relax mix vs std mix) ====================\n")
  cohen_d <- function(x,y) {
    nx<-sum(is.finite(x)); ny<-sum(is.finite(y))
    mx<-mean(x,na.rm=TRUE); my<-mean(y,na.rm=TRUE)
    sx<-stats::sd(x,na.rm=TRUE); sy<-stats::sd(y,na.rm=TRUE)
    sp <- sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx+ny-2))
    (mx - my) / sp
  }
  effect_by_n <- if (all(c("SVEM_auto_relax_mix","SVEM_auto_std_mix") %in% df$FocusSetting)) {
    df %>%
      dplyr::filter(FocusSetting %in% c("SVEM_auto_relax_mix","SVEM_auto_std_mix")) %>%
      dplyr::group_by(n_total) %>%
      dplyr::summarise(d = cohen_d(NRASE_Holdout[FocusSetting=="SVEM_auto_relax_mix"],
                                   NRASE_Holdout[FocusSetting=="SVEM_auto_std_mix"]))
  } else tibble::tibble()
  print(effect_by_n, n = Inf)

  cat("\n==================== Linear model on rank ====================\n")
  rank_mod_df <- df %>%
    dplyr::group_by(RunID, n_total) %>%
    dplyr::mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
    dplyr::ungroup()
  print(summary(stats::lm(rk ~ FocusSetting + factor(n_total) + TheoreticalR2, data = rank_mod_df)))

  # ---------- Stability ----------
  cat("\n==================== Stability across runs (SD, IQR) ====================\n")
  stab_tbl <- df %>%
    dplyr::group_by(FocusSetting, n_total) %>%
    dplyr::summarise(sd_NRASE = stats::sd(NRASE_Holdout),
                     iqr = stats::IQR(NRASE_Holdout), .groups="drop")
  print(stab_tbl, n = Inf)

  # ---------- Relax-minus-std by AUTO objective ----------
  if (all(c("SVEM_auto_relax_mix","SVEM_auto_std_mix") %in% df$FocusSetting) && "ObjectiveUsed" %in% names(df)) {
    cat("\n==================== SVEM relax-minus-std by selected AUTO objective ====================\n")
    relax_delta_by_obj <- df %>%
      dplyr::filter(FocusSetting %in% c("SVEM_auto_relax_mix","SVEM_auto_std_mix")) %>%
      dplyr::select(RunID, n_total, ObjectiveUsed, FocusSetting, NRASE_Holdout) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
      dplyr::mutate(delta = SVEM_auto_relax_mix - SVEM_auto_std_mix) %>%
      dplyr::group_by(n_total, ObjectiveUsed) %>%
      dplyr::summarise(N = dplyr::n(),
                       mean_delta = mean(delta, na.rm=TRUE),
                       se = stats::sd(delta)/sqrt(dplyr::n()),
                       .groups = "drop")
    print(relax_delta_by_obj, n = Inf)
  } else {
    relax_delta_by_obj <- tibble::tibble()
  }

  # ---------- Margin of victory ----------
  cat("\n==================== Margin of victory (2nd minus 1st) ====================\n")
  mov_tbl <- df %>%
    dplyr::group_by(RunID, n_total) %>%
    dplyr::arrange(NRASE_Holdout, .by_group = TRUE) %>%
    dplyr::summarise(winner = dplyr::first(FocusSetting),
                     mov = dplyr::nth(NRASE_Holdout, 2) - dplyr::first(NRASE_Holdout),
                     .groups = "drop") %>%
    dplyr::group_by(winner) %>%
    dplyr::summarise(mean_mov = mean(mov), median_mov = stats::median(mov), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_mov))
  print(mov_tbl)

  # ---------- Tie rate ----------
  cat("\n==================== Tie rate by n_total ====================\n")
  tie_rate <- df %>%
    dplyr::group_by(RunID, n_total) %>%
    dplyr::mutate(best = min(NRASE_Holdout, na.rm = TRUE)) %>%
    dplyr::summarise(tie = mean(abs(NRASE_Holdout - best) <= 0.005), .groups="drop") %>%
    dplyr::group_by(n_total) %>%
    dplyr::summarise(tie_rate = mean(tie), .groups="drop")
  print(tie_rate, n = Inf)

  # ---------- Boxplots ----------
  cat("\n==================== Plot: Boxplots by n_total ====================\n")
  print(
    ggplot2::ggplot(df, ggplot2::aes(FocusSetting, NRASE_Holdout, fill = FocusSetting)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.2) +
      ggplot2::facet_wrap(~ n_total, nrow = 2, scales = "free_y") +
      ggplot2::guides(fill = "none") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1))
  )

  # ---------- Learning curve params & pairwise means by n ----------
  cat("\n==================== Pairwise deltas by n_total (means) ====================\n")
  pair_means_by_n <- (function(d){
    meths <- levels(d$FocusSetting) %||% unique(d$FocusSetting)
    purrr::map_dfr(combn(meths, 2, simplify = FALSE), function(ab){
      wide <- d %>%
        dplyr::filter(FocusSetting %in% ab) %>%
        dplyr::select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
        dplyr::distinct() %>% tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
        dplyr::group_by(n_total) %>%
        dplyr::summarise(mean_delta = mean(.data[[ab[1]]] - .data[[ab[2]]], na.rm=TRUE),
                         .groups="drop")
      wide$pair <- paste(ab[1], "-", ab[2]); wide
    })
  })(df)
  print(pair_means_by_n, n = Inf)

  # ---------- Optionally save CSVs (commented) ----------
  # readr::write_csv(sum_by_n, "sim_summary_by_n.csv")
  # readr::write_csv(win_tbl,  "sim_win_tbl.csv")
  # readr::write_csv(win_ci,   "sim_win_ci.csv")

  invisible(list(
    summary_by_n = sum_by_n,
    win_tbl = win_tbl,
    avg_rank = avg_rank,
    regret_by_n = regret_by_n,
    regret_by_r2 = regret_by_r2,
    theta_mat = theta_mat,
    bt_tbl = bt_tbl,
    auto_tab = auto_tab,
    delta_relax = delta_relax,
    svem_relax_vs_std_by_n = if (exists("svem_relax_vs_std_by_n")) svem_relax_vs_std_by_n else NULL,
    svem_relax_vs_std_by_r2 = if (exists("svem_relax_vs_std_by_r2")) svem_relax_vs_std_by_r2 else NULL,
    learning_curve = lc_tbl,
    kmeans_clusters = cl_tbl,
    quantiles = quant_tbl,
    break_even = be,
    eps_tbl = eps_tbl,
    win_ci = win_ci,
    key_ci = key_ci,
    effect_by_n = effect_by_n,
    rank_lm = summary(stats::lm(rk ~ FocusSetting + factor(n_total) + TheoreticalR2, data = rank_mod_df)),
    stability = stab_tbl,
    relax_delta_by_obj = if (exists("relax_delta_by_obj")) relax_delta_by_obj else NULL,
    mov_tbl = mov_tbl,
    tie_rate = tie_rate,
    pair_means_by_n = pair_means_by_n
  ))
}

# Run once after your simulation produced df
out <- analyze_sim(df)


analyze_profiles_pca <- function(df) {
  stopifnot(all(c("n_total","TheoreticalR2","FocusSetting","NRASE_Holdout") %in% names(df)))
  prof_grid <- df %>%
    dplyr::group_by(n_total, TheoreticalR2, FocusSetting) %>%
    dplyr::summarise(m = mean(NRASE_Holdout), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = FocusSetting, values_from = m)
  X <- prof_grid %>% dplyr::select(-n_total, -TheoreticalR2) %>% as.data.frame()
  X <- X[, sapply(X, function(v) sum(is.finite(v)) == length(v)), drop = FALSE]
  pr <- prcomp(scale(X), center = TRUE, scale. = TRUE)

  print(summary(pr))
  loadings <- as.data.frame(pr$rotation[,1:2])
  loadings$Method <- rownames(loadings)
  scores <- as.data.frame(pr$x[,1:2])
  scores$n_total <- prof_grid$n_total
  scores$TheoreticalR2 <- prof_grid$TheoreticalR2

  print(
    ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, label = paste0("n=", n_total, ", R2=", TheoreticalR2))) +
      ggplot2::geom_point() + ggplot2::geom_text(nudge_y = 0.1, size = 3, alpha = 0.7) +
      ggplot2::labs(title = "Grid profiles: PCA scores")
  )
  print(
    ggplot2::ggplot(loadings, ggplot2::aes(PC1, PC2, label = Method)) +
      ggplot2::geom_point() + ggplot2::geom_text(nudge_y = 0.02, size = 4) +
      ggplot2::labs(title = "Method loadings on PC1/PC2")
  )
  invisible(list(pr = pr, scores = scores, loadings = loadings))
}




















small_sim_summary <- function(df, save_prefix = NULL) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble)
    library(ggplot2); library(scales); library(splines); library(purrr)
  })

  # --- helpers ---
  se   <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  bci  <- function(x, stat = mean, B = 2000) {
    s <- replicate(B, stat(sample(x, length(x), replace = TRUE)))
    c(lo = unname(quantile(s, .025, na.rm = TRUE)),
      med= unname(quantile(s, .5,   na.rm = TRUE)),
      hi = unname(quantile(s, .975, na.rm = TRUE)))
  }

  # Ensure FocusSetting exists (maps your 5 combos); otherwise, create it
  if (!"FocusSetting" %in% names(df)) {
    if (!all(c("Setting","AlphaScenario") %in% names(df))) {
      stop("df needs FocusSetting (or Setting + AlphaScenario).")
    }
    df <- df %>%
      mutate(
        FocusSetting = dplyr::case_when(
          Setting == "SVEM_auto_relax" & AlphaScenario == "mix"   ~ "SVEM_auto_relax_mix",
          Setting == "SVEM_auto_std"   & AlphaScenario == "mix"   ~ "SVEM_auto_std_mix",
          Setting == "SVEM_auto_std"   & AlphaScenario == "lasso" ~ "SVEM_auto_std_lasso",
          Setting == "CV_std"          & AlphaScenario == "mix"   ~ "CV_std_mix",
          Setting == "CV_std"          & AlphaScenario == "lasso" ~ "CV_std_lasso",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(FocusSetting)) %>%
      mutate(FocusSetting = factor(
        FocusSetting,
        levels = c("SVEM_auto_relax_mix",
                   "SVEM_auto_std_mix",
                   "SVEM_auto_std_lasso",
                   "CV_std_mix",
                   "CV_std_lasso")
      ))
  } else {
    df <- df %>% mutate(FocusSetting = droplevels(FocusSetting))
  }

  # ---------- TABLES ----------
  # 1) Mean error by n_total (core learning-curve table)
  tbl_by_n <- df %>%
    group_by(n_total, FocusSetting) %>%
    summarise(
      runs       = n(),
      mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
      se_NRASE   = se(NRASE_Holdout),
      mean_NAAE  = mean(NAAE_Holdout, na.rm = TRUE),
      se_NAAE    = se(NAAE_Holdout),
      .groups    = "drop"
    ) %>%
    arrange(n_total, mean_NRASE)

  # 2) Overall win rate + bootstrap CI; and average rank
  winners <- df %>%
    group_by(RunID, n_total) %>%
    slice_min(NRASE_Holdout, with_ties = TRUE) %>% ungroup()

  tbl_win <- df %>%
    select(RunID, n_total, FocusSetting) %>%
    left_join(winners %>% mutate(win = 1L),
              by = c("RunID","n_total","FocusSetting")) %>%
    mutate(win = replace(win, is.na(win), 0L)) %>%
    group_by(FocusSetting) %>%
    summarise(
      runs     = n_distinct(paste(RunID, n_total)),
      win_rate = mean(win),
      lo = bci(win)["lo"], med = bci(win)["med"], hi = bci(win)["hi"],
      .groups = "drop"
    ) %>% arrange(desc(win_rate))

  tbl_rank <- df %>%
    group_by(RunID, n_total) %>%
    mutate(rk = rank(NRASE_Holdout, ties.method = "average")) %>%
    ungroup() %>%
    group_by(FocusSetting) %>%
    summarise(mean_rank = mean(rk), se_rank = se(rk), .groups = "drop") %>%
    arrange(mean_rank)

  # 3) One key pairwise delta for the paper: SVEM relax mix vs CV mix (overall)
  pair_delta <- (function(d, a, b, lab = NULL) {
    wide <- d %>%
      filter(FocusSetting %in% c(a, b)) %>%
      select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
      distinct() %>% tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
      mutate(delta = .data[[a]] - .data[[b]]) %>% pull(delta) %>% na.omit()
    ci <- bci(wide, stat = mean)
    tibble(
      pair = if (is.null(lab)) paste(a,"-",b) else lab,
      N = length(wide),
      mean = mean(wide),
      lo = ci["lo"], med = ci["med"], hi = ci["hi"]
    )
  })(df, "SVEM_auto_relax_mix", "CV_std_mix", "SVEM(auto,relax,mix) − CV(mix)")

  # 4) Auto objective usage by n_total (SVEM only, %)
  tbl_objective <- df %>%
    filter(grepl("^SVEM", FocusSetting)) %>%
    count(n_total, ObjectiveUsed, name = "n") %>%
    group_by(n_total) %>%
    mutate(p = n / sum(n)) %>% ungroup() %>%
    arrange(n_total, ObjectiveUsed)

  # 5) Head-to-head: SVEM relax (mix) vs std (mix) by n_total (mean deltas)
  tbl_relax_vs_std_by_n <- df %>%
    filter(FocusSetting %in% c("SVEM_auto_relax_mix","SVEM_auto_std_mix")) %>%
    select(RunID, n_total, FocusSetting, NRASE_Holdout) %>%
    distinct() %>% tidyr::pivot_wider(names_from = FocusSetting, values_from = NRASE_Holdout) %>%
    mutate(delta = SVEM_auto_relax_mix - SVEM_auto_std_mix) %>%
    group_by(n_total) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              median_delta = median(delta, na.rm = TRUE),
              .groups = "drop")

  # 6) Oracle gap (regret) by n_total (smaller is better)
  tbl_regret_by_n <- df %>%
    group_by(RunID, n_total) %>%
    mutate(best = min(NRASE_Holdout, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(regret = NRASE_Holdout - best) %>%
    group_by(n_total, FocusSetting) %>%
    summarise(mean_regret = mean(regret, na.rm = TRUE),
              se = se(regret), .groups = "drop") %>%
    arrange(n_total, mean_regret)

  # ---------- PLOTS (clean, journal-friendly) ----------
  # A) Learning curves: mean NRASE ± 1 SE
  p_learning <- ggplot(tbl_by_n,
                       aes(x = n_total, y = mean_NRASE, color = FocusSetting, group = FocusSetting)) +
    geom_line(linewidth = 0.9) +        # <- was size = 0.9
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_NRASE - se_NRASE, ymax = mean_NRASE + se_NRASE),
                  width = 0.0, alpha = 0.7) +
    labs(x = "Training runs (n_total)", y = "Mean NRASE (±1 SE)", color = "Method") +
    theme_minimal(base_size = 12)


  # B) Win rates with bootstrap CIs
  p_win <- ggplot(tbl_win, aes(x = reorder(FocusSetting, -win_rate), y = win_rate, fill = FocusSetting)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "Win rate (95% bootstrap CI)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  # C) Objective mix vs n_total (SVEM only)
  p_obj <- df %>%
    filter(grepl("^SVEM", FocusSetting), !is.na(ObjectiveUsed)) %>%
    count(n_total, ObjectiveUsed, name = "n") %>%
    group_by(n_total) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
    ggplot(aes(x = n_total, y = p, fill = ObjectiveUsed)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "n_total", y = "Objective share", fill = "Auto objective") +
    theme_minimal(base_size = 12)

  # D) Compact boxplots by n_total (optional but handy)
  p_box <- ggplot(df, aes(FocusSetting, NRASE_Holdout, fill = FocusSetting)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~ n_total, nrow = 2, scales = "free_y") +
    guides(fill = "none") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = "NRASE")

  # ---------- Optional disk save ----------
  if (!is.null(save_prefix)) {
    ggsave(paste0(save_prefix, "_learning_curves.png"), p_learning, width = 7, height = 4.2, dpi = 300)
    ggsave(paste0(save_prefix, "_winrates.png"),        p_win,      width = 6.2, height = 3.8, dpi = 300)
    ggsave(paste0(save_prefix, "_objective_mix.png"),   p_obj,      width = 6.2, height = 3.8, dpi = 300)
    ggsave(paste0(save_prefix, "_boxplots.png"),        p_box,      width = 9.0, height = 5.8, dpi = 300)
  }

  # ---------- Print the most salient tables ----------
  cat("\n== Mean error by n_total (NRASE, NAAE) ==\n"); print(tbl_by_n, n = Inf)
  cat("\n== Win rates (with bootstrap CIs) ==\n");     print(tbl_win)
  cat("\n== Average rank (lower is better) ==\n");     print(tbl_rank)
  cat("\n== Key delta: SVEM(auto,relax,mix) − CV(mix) (NRASE) ==\n"); print(pair_delta)
  cat("\n== SVEM objective usage by n_total ==\n");    print(tbl_objective, n = Inf)
  cat("\n== SVEM relax (mix) − std (mix) by n_total ==\n"); print(tbl_relax_vs_std_by_n, n = Inf)
  cat("\n== Oracle regret (mean) by n_total ==\n");    print(tbl_regret_by_n, n = Inf)

  # ---------- Return a tidy bundle ----------
  invisible(list(
    tables = list(
      by_n        = tbl_by_n,
      win_rates   = tbl_win,
      avg_rank    = tbl_rank,
      key_delta   = pair_delta,
      objective   = tbl_objective,
      relax_vs_std_by_n = tbl_relax_vs_std_by_n,
      regret_by_n = tbl_regret_by_n
    ),
    plots = list(
      learning_curves = p_learning,
      win_rates       = p_win,
      objective_mix   = p_obj,
      boxplots        = p_box
    )
  ))
}


out_small <- small_sim_summary(df)                 # just print + return
# or save figures:
out_small <- small_sim_summary(df, save_prefix = "svem_paper_fig")



