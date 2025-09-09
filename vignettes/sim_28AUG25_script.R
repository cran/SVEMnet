# Focused small-n RELAX vs STD study (no CV; richer alpha test; N_BOOT=500)
run_gap_simulation <- function(
    design = c("smalln","light","medium","full"),
    CORES = 12,
    N_BOOT = 500,                 # <- bumped to 500
    OUT_ITERS = NULL,
    seed = 20260101
) {
  design <- match.arg(design)

  suppressPackageStartupMessages({
    library(SVEMnet)
    library(glmnet)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(furrr)
    library(ggplot2)
    library(scales)
    library(stringr)
  })

  # ---------- Design presets (new: smalln) ----------
  if (is.null(OUT_ITERS)) {
    OUT_ITERS <- switch(design,
                        smalln = 50,            # higher reps for precision at tiny n
                        light  = 20, medium = 30, full = 30
    )
  }

  grid <- switch(
    design,
    smalln = list(
      N_TOTAL_SEQ   = c(20,25,30,35,40,45),
      R2_LEVELS     = c(0.5,0.7, 0.9),
      NOISE_TYPES   = c("gaussian","t3","contam"),
      OUTLIER_RATES = c(0, 0.05),
      MISSPEC       = c("none"),       # misspec dropped to save budget
      WEIGHTS       = c("SVEM"),
      OBJECTIVES    = c("auto"),
      ALPHA_SCEN    = c("lasso","rich")  # lasso=1 vs rich=c(1,.75,.5,.25)
    ),
    light = list(
      N_TOTAL_SEQ   = c(15,25,35,45,55,65,75,100,150),
      R2_LEVELS     = c(0.5, 0.9),
      NOISE_TYPES   = c("gaussian","t3"),
      OUTLIER_RATES = c(0, 0.05),
      MISSPEC       = c("none","smooth"),
      WEIGHTS       = c("SVEM"),
      OBJECTIVES    = c("auto"),
      ALPHA_SCEN    = c("lasso","rich")
    ),
    medium = list(
      N_TOTAL_SEQ   = c(15,25,35,45,55,65,75,100,150,200),
      R2_LEVELS     = c(0.3, 0.5, 0.9),
      NOISE_TYPES   = c("gaussian","t3","contam"),
      OUTLIER_RATES = c(0, 0.05),
      MISSPEC       = c("none"),
      WEIGHTS       = c("SVEM"),
      OBJECTIVES    = c("auto"),
      ALPHA_SCEN    = c("lasso","rich")
    ),
    full = list(
      N_TOTAL_SEQ   = c(25,35,65,100,200),
      R2_LEVELS     = c(0.3,0.5,0.9),
      NOISE_TYPES   = c("gaussian","t3","contam","hetero"),
      OUTLIER_RATES = c(0, 0.05),
      MISSPEC       = c("none","smooth"),
      WEIGHTS       = c("SVEM"),
      OBJECTIVES    = c("auto"),
      ALPHA_SCEN    = c("lasso","rich")
    )
  )

  # ---------- Utilities ----------
  rdirichlet <- function(n, alpha) {
    k <- length(alpha)
    x <- matrix(stats::rgamma(n * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
    x / rowSums(x)
  }

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
    names(out) <- c("A","B","C","D")
    out
  }

  gen_pv <- function(beta5_8 = c(8, 2), beta9_25 = c(5, 5)) {
    stopifnot(length(beta5_8)==2, length(beta9_25)==2)
    rexp_signed <- function(n) stats::rexp(n) - stats::rexp(n)
    p5_8  <- stats::rbeta(1, beta5_8[1], beta5_8[2])
    p9_25 <- stats::rbeta(1, beta9_25[1], beta9_25[2])
    p1_4   <- rexp_signed(4)
    p5_8v  <- rexp_signed(4)  * rbinom(4, 1, p5_8)
    p9_25v <- rexp_signed(17) * rbinom(17, 1, p9_25)
    list(pv = c(p1_4, p5_8v, p9_25v), p5_8 = p5_8, p9_25 = p9_25, scheme = "beta")
  }

  true_response <- function(df, pv, misspec = c("none","smooth")) {
    if (is.list(pv)) pv <- pv$pv
    stopifnot(is.numeric(pv) && length(pv)==25)
    misspec <- match.arg(misspec)
    s  <- 1 - 0.1
    zA <- (df$A - 0.1)/s; zB <- df$B/s; zC <- df$C/s; zD <- df$D/s
    Esign <- ifelse(df$E==0, 1, -1)

    part1 <- pv[1]*zA + pv[2]*zB + pv[3]*zC + pv[4]*zD +
      (pv[5]*zA + pv[6]*zB + pv[7]*zC + pv[8]*zD) * Esign
    part2 <- 4*(pv[9]*zA*zB + pv[10]*zA*zC + pv[11]*zA*zD +
                  pv[12]*zB*zC + pv[13]*zB*zD + pv[14]*zC*zD)
    part3 <- 27*(pv[15]*zA*zB*zC + pv[16]*zA*zB*zD +
                   pv[17]*zA*zC*zD + pv[18]*zB*zC*zD)
    part4 <- 27*(pv[19]*zB*zA*(zA - zB) +
                   pv[20]*zC*zA*(zA - zC) +
                   pv[21]*zC*zB*(zB - zC) +
                   pv[22]*zD*zA*(zA - zD) +
                   pv[23]*zD*zB*(zB - zD) +
                   pv[24]*zD*zC*(zC - zD))
    part5 <- 256 * pv[25] * zA*zB*zC*zD
    base <- part1 + part2 + part3 + part4 + part5

    if (misspec == "smooth") {
      base <- base + 0.20*sin(6*pi*df$A) + 0.20*cos(6*pi*df$B) + 0.10*(df$C^2 - df$D^2)
    }
    base
  }

  add_noise <- function(mu, target_sd, type = c("gaussian","t3","contam","hetero")) {
    type <- match.arg(type); n <- length(mu)
    if (type == "gaussian") {
      eps <- rnorm(n, sd = target_sd)
    } else if (type == "t3") {
      z <- rt(n, df = 3); eps <- z / sd(z) * target_sd
    } else if (type == "contam") {
      base <- rnorm(n, sd = target_sd)
      big  <- rnorm(n, sd = 3*target_sd)
      mix  <- rbinom(n, 1, 0.05)
      eps <- ifelse(mix==1, big, base)
    } else {
      scale <- 0.5*target_sd + 0.8*target_sd * scales::rescale(abs(mu), to = c(0,1))
      eps <- rnorm(n, sd = scale)
    }
    eps
  }

  inject_y_outliers <- function(y, rate = 0, scale = 10) {
    if (rate <= 0) return(y)
    n <- length(y)
    k <- rbinom(n, 1, rate)
    y + k * (rt(n, df = 3) * sd(y, na.rm = TRUE) * scale / 3)
  }

  # ---------- Fitting helpers (SVEM only) ----------
  fit_svem <- function(form, train, alpha_vec, objective, weight_scheme, relaxed, N_BOOT) {
    withCallingHandlers({
      SVEMnet::SVEMnet(
        formula        = form,
        data           = train,
        glmnet_alpha   = alpha_vec,
        nBoot          = N_BOOT,
        objective      = objective,     # "auto"
        weight_scheme  = weight_scheme, # "SVEM"
        standardize    = TRUE,
        relaxed        = relaxed
      )
    }, warning = function(w) invokeRestart("muffleWarning"))
  }

  # ---------- Metrics ----------
  metric_one <- function(yhat, y_true, sd_scale) {
    err  <- yhat - y_true
    tibble(
      NRASE_Holdout = sqrt(mean(err^2)) / sd_scale,
      NAAE_Holdout  = mean(abs(err))    / sd_scale
    )
  }

  # ---------- One scenario/run ----------
  run_one <- function(run_id, n_total, r2, noise_type, out_rate, misspec, weight_scheme, objective, alpha_mode) {
    pv <- gen_pv()

    HOLDOUT_N <- 10000
    hold <- sample_mixture(HOLDOUT_N); hold$E <- sample(c(0, 0.002), size = HOLDOUT_N, replace = TRUE)
    hold_true <- true_response(hold, pv, misspec)
    hold <- hold %>% mutate(E = factor(E))
    sd_hold_true <- sd(hold_true)

    DIST_N <- 10000
    dist_pts <- sample_mixture(DIST_N); dist_pts$E <- sample(c(0, 0.002), DIST_N, replace = TRUE)
    dist_true <- true_response(dist_pts, pv, misspec)
    y_sd_global <- sd(dist_true)

    tr <- sample_mixture(n_total); tr$E <- sample(c(0, 0.002), n_total, replace = TRUE)
    tr_true <- true_response(tr, pv, misspec)
    err_sd <- y_sd_global * sqrt((1 - r2) / r2)

    eps <- add_noise(mu = tr_true, target_sd = err_sd, type = noise_type)
    Y <- tr_true + eps
    Y <- inject_y_outliers(Y, rate = out_rate, scale = 10)

    train <- tr %>% mutate(Y = Y, E = factor(E))

    form <- Y ~ (A + B + C + D + E)^2 + A:B:C + A:B:D + A:C:D + B:C:D + A:B:C:D

    alpha_vec <- if (alpha_mode == "lasso") 1 else c(1, 0.75, 0.5, 0.25)

    get_safe <- function(expr) tryCatch(expr, error = function(e) e)

    # SVEM (std/relax) only
    m_svem_std   <- get_safe(fit_svem(form, train, alpha_vec, objective, weight_scheme, relaxed = FALSE, N_BOOT))
    m_svem_relax <- get_safe(fit_svem(form, train, alpha_vec, objective, weight_scheme, relaxed = TRUE,  N_BOOT))

    # predictions
    ps <- function(m, hold) if (inherits(m, "error")) NA_real_ else as.numeric(predict(m, hold, debias = FALSE))

    preds <- list(
      SVEM_std   = ps(m_svem_std,   hold),
      SVEM_relax = ps(m_svem_relax, hold)
    )

    svem_obj_used <- c(
      SVEM_std   = if (inherits(m_svem_std, "error")) NA_character_ else m_svem_std$objective_used,
      SVEM_relax = if (inherits(m_svem_relax,"error")) NA_character_ else m_svem_relax$objective_used
    )

    purrr::imap_dfr(preds, function(yhat, name) {
      if (all(!is.finite(yhat))) {
        tibble(
          RunID = run_id, n_total = n_total, TheoreticalR2 = r2,
          Noise = noise_type, OutlierRate = out_rate, Misspec = misspec,
          Weights = weight_scheme, Objective = objective, AlphaScenario = alpha_mode,
          Setting = name, NRASE_Holdout = NA_real_, NAAE_Holdout = NA_real_,
          ObjectiveUsed = svem_obj_used[[name]]
        )
      } else {
        m <- metric_one(yhat, hold_true, sd_hold_true)
        tibble(
          RunID = run_id, n_total = n_total, TheoreticalR2 = r2,
          Noise = noise_type, OutlierRate = out_rate, Misspec = misspec,
          Weights = weight_scheme, Objective = objective, AlphaScenario = alpha_mode,
          Setting = name, NRASE_Holdout = m$NRASE_Holdout, NAAE_Holdout = m$NAAE_Holdout,
          ObjectiveUsed = svem_obj_used[[name]]
        )
      }
    })
  }

  # ---------- Scenario grid ----------
  SC <- expand.grid(
    n_total       = grid$N_TOTAL_SEQ,
    TheoreticalR2 = grid$R2_LEVELS,
    Noise         = grid$NOISE_TYPES,
    OutlierRate   = grid$OUTLIER_RATES,
    Misspec       = grid$MISSPEC,
    Weights       = grid$WEIGHTS,
    Objective     = grid$OBJECTIVES,
    AlphaScenario = grid$ALPHA_SCEN,
    stringsAsFactors = FALSE
  )
  if (design != "full") SC <- SC %>% filter(!(Noise=="hetero")) %>% arrange(dplyr::across(everything()))

  jobs <- SC %>%
    mutate(iter = list(seq_len(OUT_ITERS))) %>% tidyr::unnest(iter) %>%
    mutate(RunID = sprintf("run%06d", row_number()))

  set.seed(seed)
  future::plan(future::multisession, workers = CORES)
  cat(sprintf("\nRunning %d jobs across %d workers...\n", nrow(jobs), CORES)); flush.console()

  df <- furrr::future_pmap_dfr(
    .l = list(jobs$RunID, jobs$n_total, jobs$TheoreticalR2, jobs$Noise,
              jobs$OutlierRate, jobs$Misspec, jobs$Weights, jobs$Objective, jobs$AlphaScenario),
    .f = run_one,
    .options = furrr::furrr_options(seed = TRUE)
  ) %>%
    mutate(across(c(Noise, Misspec, Weights, Objective, AlphaScenario, Setting), ~factor(.)))

  # ---------- Helpers ----------
  se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  boot_ci <- function(x, stat = mean, B = 2000) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(`2.5%`=NA, `50%`=NA, `97.5%`=NA))
    s <- replicate(B, stat(sample(x, length(x), replace = TRUE)))
    stats::quantile(s, c(.025,.5,.975), na.rm = TRUE)
  }
  pair_level <- function(d, a, b) {
    base <- d %>%
      filter(Setting %in% c(a,b)) %>%
      select(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec, Weights, Objective, AlphaScenario,
             Setting, NRASE_Holdout) %>%
      distinct() %>%
      tidyr::pivot_wider(names_from = Setting, values_from = NRASE_Holdout)
    if (!all(c(a,b) %in% names(base))) return(tibble())
    base %>% mutate(delta = .data[[a]] - .data[[b]])
  }

  # ---------- Summaries (SVEM only) ----------
  cat("\n==================== OVERALL MEANS ====================\n")
  overall_tbl <- df %>%
    group_by(Setting) %>%
    summarise(N = sum(is.finite(NRASE_Holdout)),
              mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
              se_NRASE = se(NRASE_Holdout),
              mean_NAAE = mean(NAAE_Holdout, na.rm = TRUE),
              se_NAAE = se(NAAE_Holdout),
              .groups = "drop") %>% arrange(mean_NRASE)
  print(overall_tbl)

  cat("\n==================== MEAN by n_total (collapsed) ====================\n")
  by_n_tbl <- df %>%
    group_by(n_total, Setting) %>%
    summarise(N = sum(is.finite(NRASE_Holdout)),
              mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
              se_NRASE = se(NRASE_Holdout),
              .groups = "drop") %>% arrange(n_total, mean_NRASE)
  print(by_n_tbl, n = Inf)

  cat("\n==================== SVEM: AUTO (objective used) ====================\n")
  obj_mix <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax"), Objective == "auto") %>%
    count(n_total, ObjectiveUsed, name = "n") %>%
    group_by(n_total) %>% mutate(p = n/sum(n)) %>% ungroup() %>%
    arrange(n_total, desc(p))
  print(obj_mix, n = Inf)

  cat("\n==================== RELAX − STD (AUTO) ====================\n")
  relax_deltas <- pair_level(df %>% filter(Objective=="auto"), "SVEM_relax", "SVEM_std")
  relax_tbl <- relax_deltas %>%
    group_by(n_total, TheoreticalR2, Noise, OutlierRate, AlphaScenario) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              med_delta  = median(delta, na.rm = TRUE),
              harm_02 = mean(delta > 0.02, na.rm = TRUE),
              harm_05 = mean(delta > 0.05, na.rm = TRUE),
              harm_10 = mean(delta > 0.10, na.rm = TRUE),
              tie_01  = mean(abs(delta) <= 0.01, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(n_total, TheoreticalR2, Noise, OutlierRate, AlphaScenario)
  print(relax_tbl, n = 50)

  # Small-n focus panel (n in [20..50]) aggregated to decide default
  cat("\n==================== SMALL-n DECISION PANEL (n=20–50) ====================\n")
  smalln_set <- unique(df$n_total)
  smalln_set <- smalln_set[smalln_set >= 20 & smalln_set <= 50]
  decision_tbl <- relax_deltas %>%
    filter(n_total %in% smalln_set) %>%
    group_by(TheoreticalR2, Noise, OutlierRate, AlphaScenario) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              med_delta  = median(delta, na.rm = TRUE),
              lo_mean = boot_ci(delta)[1], hi_mean = boot_ci(delta)[3],
              P_relax_le_std = mean(delta <= 0, na.rm = TRUE),
              P_harm_gt_0.05 = mean(delta > 0.05, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(Noise, TheoreticalR2, OutlierRate, AlphaScenario)
  print(decision_tbl, n = Inf)

  # Mann–Whitney θ = P(RELAX < STD) + 0.5 P(=)
  mann_theta <- function(x,y){ m <- outer(x,y,"-"); mean((m<0) + 0.5*(m==0), na.rm=TRUE) }
  theta_tbl <- relax_deltas %>%
    group_by(n_total, TheoreticalR2, Noise, OutlierRate, AlphaScenario) %>%
    summarise(theta = mann_theta(SVEM_relax, SVEM_std), .groups="drop") %>%
    group_by(TheoreticalR2, Noise, OutlierRate, AlphaScenario) %>%
    summarise(theta_smalln = mean(theta[n_total %in% smalln_set], na.rm=TRUE), .groups="drop")
  cat("\n==================== MANN–WHITNEY θ (small-n mean across n) ====================\n")
  print(theta_tbl, n = Inf)

  # Alpha competition summary in small-n range
  cat("\n==================== ALPHA: lasso vs rich (small-n) ====================\n")
  alpha_effect <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax"),
           Objective=="auto",
           n_total %in% smalln_set,
           AlphaScenario %in% c("lasso","rich")) %>%
    group_by(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Setting, AlphaScenario) %>%
    summarise(NR = mean(NRASE_Holdout, na.rm = TRUE), .groups="drop") %>%
    tidyr::pivot_wider(names_from = AlphaScenario, values_from = NR) %>%
    mutate(delta_rich_minus_lasso = rich - lasso) %>%
    group_by(Setting, TheoreticalR2, Noise, OutlierRate) %>%
    summarise(N = sum(is.finite(delta_rich_minus_lasso)),
              mean_delta = mean(delta_rich_minus_lasso, na.rm=TRUE),
              med_delta  = median(delta_rich_minus_lasso, na.rm=TRUE),
              win_rich_le_lasso = mean(delta_rich_minus_lasso <= 0, na.rm=TRUE),
              .groups="drop") %>%
    arrange(Setting, Noise, TheoreticalR2, OutlierRate)
  print(alpha_effect, n = Inf)

  # Minimal plot objects (silent)
  p_by_n <- ggplot(by_n_tbl, aes(x = n_total, y = mean_NRASE, colour = Setting)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    labs(x = "n_total", y = "Mean NRASE", colour = "Method")

  invisible(list(
    df = df,
    overall_tbl = overall_tbl,
    by_n_tbl = by_n_tbl,
    obj_mix = obj_mix,
    relax_tbl = relax_tbl,
    decision_tbl = decision_tbl,
    theta_smalln = theta_tbl,
    alpha_effect_smalln = alpha_effect,
    plots = list(by_n = p_by_n),
    design = list(grid = grid, scenarios = SC, OUT_ITERS = OUT_ITERS, CORES = CORES, N_BOOT = N_BOOT)
  ))
}

# ----------------- HOW TO RUN (small-n focus) -----------------
 out <- run_gap_simulation(design = "smalln", CORES = 15)




 library(dplyr)

 dec <- out$decision_tbl
 # Aggregate across Noise & OutlierRate (keep R2 split) to make the call
 call_tbl <- dec %>%
   group_by(TheoreticalR2, AlphaScenario) %>%
   summarise(
     mean_delta = mean(mean_delta, na.rm=TRUE),
     hi_mean    = mean(hi_mean, na.rm=TRUE),
     P_relax_le_std = mean(P_relax_le_std, na.rm=TRUE),
     P_harm_gt_0.05 = mean(P_harm_gt_0.05, na.rm=TRUE),
     .groups="drop"
   ) %>%
   arrange(AlphaScenario, TheoreticalR2)

 print(call_tbl, n = Inf)
