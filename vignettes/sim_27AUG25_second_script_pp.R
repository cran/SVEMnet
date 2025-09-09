# ===================== GAP-FILLING SIMULATION (w/ CV_relax EXTRAS) =====================
# Adds: "Is CV_relax worse?" summary, win rates incl CV_relax, oracle gap (4 methods),
# and a compact 4x4 Mann–Whitney θ matrix. Returns all new tables in the output list.

run_gap_simulation <- function(
    design = c("light","medium","full","A1"),
    CORES = 12,
    N_BOOT = 300,
    OUT_ITERS = NULL,
    seed = 202601
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

  # ---------- Design presets ----------
  if (is.null(OUT_ITERS)) OUT_ITERS <- switch(design, light = 20, medium = 30, full = 40,A1 = 25)
  grid <- switch(
    design,
    light = list(
      N_TOTAL_SEQ  = c(15,25,35,45,55,65,75,100,150),
      R2_LEVELS    = c(0.5, 0.9),
      NOISE_TYPES  = c("gaussian","t3"),
      OUTLIER_RATES= c(0, 0.05),
      MISSPEC      = c("none","smooth"),
      WEIGHTS      = c("SVEM"),
      OBJECTIVES   = c("auto","wSSE"),
      ALPHA_SCEN   = c("mix")  # {0.5,1}
    ),
    medium = list(
      N_TOTAL_SEQ  = c(15,25,35,45,55,65,75,100,150,200),
      R2_LEVELS    = c(0.3, 0.5, 0.9),
      NOISE_TYPES  = c("gaussian","t3","contam"),
      OUTLIER_RATES= c(0, 0.05),
      MISSPEC      = c("none","smooth"),
      WEIGHTS      = c("SVEM"),
      OBJECTIVES   = c("auto","wSSE"),
      ALPHA_SCEN   = c("mix","lasso")
    ),
    full = list(
      N_TOTAL_SEQ  = c(15,25,35,45,55,65,75,100,150,200,300),
      R2_LEVELS    = c(0.3,0.5,0.7,0.9),
      NOISE_TYPES  = c("gaussian","t3","contam","hetero"),
      OUTLIER_RATES= c(0, 0.05),
      MISSPEC      = c("none","smooth"),
      WEIGHTS      = c("SVEM"),
      OBJECTIVES   = c("auto","wSSE"),
      ALPHA_SCEN   = c("mix","lasso")
    ),
    A1 = list(
      N_TOTAL_SEQ  = c(25,45,100),
      R2_LEVELS    = c(0.4,0.7,0.9),
      NOISE_TYPES  = c("gaussian","t3","contam","hetero"),
      OUTLIER_RATES= c(0, 0.07),
      MISSPEC      = c("none","smooth"),
      WEIGHTS      = c("SVEM"),
      OBJECTIVES   = c("auto"),
      ALPHA_SCEN   = c("mix","lasso")
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

  # ---------- Fitting helpers ----------
  fit_svem <- function(form, train, alpha_vec, objective, weight_scheme, relaxed, N_BOOT) {
    withCallingHandlers({
      SVEMnet::SVEMnet(
        formula        = form,
        data           = train,
        glmnet_alpha   = alpha_vec,
        nBoot          = N_BOOT,
        objective      = objective,     # "auto" or "wSSE"
        weight_scheme  = weight_scheme, # "SVEM" or "equal"
        standardize    = TRUE,
        relaxed        = relaxed
      )
    }, warning = function(w) invokeRestart("muffleWarning"))
  }

  fit_cv <- function(form, train, alpha_vec, relaxed = FALSE) {
    suppressWarnings(
      glmnet_with_cv(
        formula       = form,
        data          = train,
        glmnet_alpha  = alpha_vec,
        nfolds        = 10,
        repeats       = 3,
        choose_rule   = "1se",
        standardize   = TRUE,
        relaxed       = relaxed
      )
    )
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
    alpha_vec <- if (alpha_mode == "mix") c(0.25,0.5, 1) else 1

    get_safe <- function(expr) tryCatch(expr, error = function(e) e)

    # SVEM (std/relax)
    m_svem_std   <- get_safe(fit_svem(form, train, alpha_vec, objective, weight_scheme, relaxed = FALSE, N_BOOT))
    m_svem_relax <- get_safe(fit_svem(form, train, alpha_vec, objective, weight_scheme, relaxed = TRUE,  N_BOOT))

    # CV (std/relax)
    m_cv_std     <- get_safe(fit_cv(form, train, alpha_vec, relaxed = FALSE))
    m_cv_relax   <- get_safe(fit_cv(form, train, alpha_vec, relaxed = TRUE))

    # predictions
    p <- function(m, hold) if (inherits(m, "error")) NA_real_ else as.numeric(predict_cv(m, hold, debias = FALSE))
    ps <- function(m, hold) if (inherits(m, "error")) NA_real_ else as.numeric(predict(m, hold, debias = FALSE))

    preds <- list(
      SVEM_std   = ps(m_svem_std,   hold),
      SVEM_relax = ps(m_svem_relax, hold),
      CV_std     = p(m_cv_std,      hold),
      CV_relax   = p(m_cv_relax,    hold)
    )

    svem_obj_used <- c(
      SVEM_std   = if (inherits(m_svem_std, "error")) NA_character_ else m_svem_std$objective_used,
      SVEM_relax = if (inherits(m_svem_relax,"error")) NA_character_ else m_svem_relax$objective_used,
      CV_std     = NA_character_,
      CV_relax   = NA_character_
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

  # ---------- Summaries ----------
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

  cat("\n==================== SVEM: AUTO vs SSE (collapsed) ====================\n")
  svem_only <- df %>% filter(Setting %in% c("SVEM_std","SVEM_relax"))
  auto_vs_sse <- svem_only %>%
    group_by(Setting, Objective) %>%
    summarise(mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
              se_NRASE = se(NRASE_Holdout), N = sum(is.finite(NRASE_Holdout)),
              .groups = "drop") %>%
    arrange(Setting, mean_NRASE)
  print(auto_vs_sse)

  cat("\n==================== AUTO objective mix (AIC vs BIC) ====================\n")
  obj_mix <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax"), Objective == "auto") %>%
    count(n_total, ObjectiveUsed, name = "n") %>%
    group_by(n_total) %>% mutate(p = n/sum(n)) %>% ungroup() %>%
    arrange(n_total, desc(p))
  print(obj_mix, n = Inf)

  cat("\n==================== SVEM(relax - std), AUTO only ====================\n")
  relax_deltas <- pair_level(df %>% filter(Objective=="auto"), "SVEM_relax", "SVEM_std")
  relax_tbl <- relax_deltas %>%
    group_by(n_total) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              med_delta  = median(delta, na.rm = TRUE),
              .groups = "drop")
  print(relax_tbl, n = Inf)

  # ============== NEW: Is CV_relax worse than CV_std? (Overall + by n) ==============
  cat("\n==================== Is CV_relax worse than CV_std? ====================\n")
  cv_deltas_all <- pair_level(df, "CV_relax", "CV_std")
  if (nrow(cv_deltas_all)) {
    cvd <- cv_deltas_all$delta
    tt <- t.test(cvd)
    wt <- suppressWarnings(wilcox.test(cvd, exact = FALSE))
    ci <- boot_ci(cvd, stat = mean, B = 2000)
    overall_cv_relax_summary <- tibble(
      N = sum(is.finite(cvd)),
      mean_delta = mean(cvd, na.rm = TRUE),
      median_delta = median(cvd, na.rm = TRUE),
      t_p = unname(tt$p.value),
      wilcox_p = unname(wt$p.value),
      lo = unname(ci[1]), med = unname(ci[2]), hi = unname(ci[3]),
      win_rate_relax_le_std = mean(cvd <= 0, na.rm = TRUE)
    )
    print(overall_cv_relax_summary)

    cat("\n-- By n_total --\n")
    cv_relax_tbl <- cv_deltas_all %>%
      group_by(n_total) %>%
      summarise(
        N = sum(is.finite(delta)),
        mean_delta = mean(delta, na.rm = TRUE),
        median_delta = median(delta, na.rm = TRUE),
        win_rate_relax_le_std = mean(delta <= 0, na.rm = TRUE),
        .groups = "drop"
      )
    print(cv_relax_tbl, n = Inf)
  } else {
    message("No CV_relax vs CV_std pairs found.")
    overall_cv_relax_summary <- tibble()
    cv_relax_tbl <- tibble()
  }

  cat("\n==================== Head-to-head: SVEM_relax (AUTO, mix) vs CV_relax (mix) ====================\n")
  h2h_relax <- pair_level(df %>% filter(Objective=="auto", AlphaScenario=="mix"),
                          "SVEM_relax", "CV_relax")
  h2h_relax_tbl <- h2h_relax %>%
    group_by(n_total, Noise, OutlierRate, Misspec, TheoreticalR2) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              med_delta  = median(delta, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(n_total, Noise, OutlierRate, Misspec, TheoreticalR2)
  print(h2h_relax_tbl, n = 50)

  cat("\n==================== Head-to-head: SVEM_std (AUTO, mix) vs CV_std (mix) ====================\n")
  h2h_std <- pair_level(df %>% filter(Objective=="auto", AlphaScenario=="mix"),
                        "SVEM_std", "CV_std")
  h2h_std_tbl <- h2h_std %>%
    group_by(n_total, Noise, OutlierRate, Misspec, TheoreticalR2) %>%
    summarise(N = sum(is.finite(delta)),
              mean_delta = mean(delta, na.rm = TRUE),
              med_delta  = median(delta, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(n_total, Noise, OutlierRate, Misspec, TheoreticalR2)
  print(h2h_std_tbl, n = 50)

  # ============== NEW: Win-rates across ALL four methods ==============
  cat("\n==================== Win rate (all four methods) ====================\n")
  win_rate_all <- df %>%
    group_by(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec, Weights, Objective, AlphaScenario) %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax","CV_std","CV_relax")) %>%
    mutate(best = min(NRASE_Holdout, na.rm = TRUE),
           is_best = NRASE_Holdout <= best + 1e-12) %>%
    ungroup() %>%
    group_by(Setting) %>%
    summarise(wins = sum(is_best, na.rm = TRUE),
              blocks = n_distinct(paste(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec, Weights, Objective, AlphaScenario)),
              win_rate = wins / blocks,
              .groups = "drop") %>%
    arrange(desc(win_rate))
  print(win_rate_all)

  # ============== NEW: Oracle gap (regret) by n_total, all four ==============
  cat("\n==================== ORACLE GAP by n_total (all four) ====================\n")
  regret_by_n <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax","CV_std","CV_relax")) %>%
    group_by(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec, Weights, Objective, AlphaScenario) %>%
    mutate(best = min(NRASE_Holdout, na.rm = TRUE),
           regret = NRASE_Holdout - best) %>%
    ungroup() %>%
    group_by(n_total, Setting) %>%
    summarise(mean_regret = mean(regret, na.rm = TRUE),
              se = se(regret),
              .groups = "drop") %>%
    arrange(n_total, mean_regret)
  print(regret_by_n, n = Inf)

  # ============== NEW: Pairwise P(i beats j) θ matrix (four methods) ==============
  cat("\n==================== Pairwise P(i beats j) (θ) among four ====================\n")
  # Build wide blocks
  four <- df %>%
    filter(Setting %in% c("SVEM_relax","SVEM_std","CV_relax","CV_std")) %>%
    select(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec, Weights, Objective, AlphaScenario, Setting, NRASE_Holdout) %>%
    distinct()
  meths <- c("SVEM_relax","SVEM_std","CV_relax","CV_std")

  mann_theta <- function(x, y) {
    # P(X<Y) + 0.5 P(X=Y)
    # robust to ties
    m <- outer(x, y, FUN = "-")
    mean((m < 0) + 0.5*(m == 0), na.rm = TRUE)
  }

  theta_4x4 <- matrix(NA_real_, nrow = 4, ncol = 4, dimnames = list(meths, meths))
  # pool across all blocks (same way as earlier, but for four)
  wide <- four %>%
    tidyr::pivot_wider(names_from = Setting, values_from = NRASE_Holdout)

  for (i in seq_along(meths)) for (j in seq_along(meths)) {
    if (i == j) next
    xi <- wide[[meths[i]]]; yj <- wide[[meths[j]]]
    theta_4x4[i, j] <- mann_theta(xi, yj)
  }
  print(round(theta_4x4, 3))

  # ---------- Robustness slices & learning curves ----------
  cat("\n==================== Robustness slices (AUTO, SVEM only) ====================\n")
  robust_tbl <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax"), Objective=="auto") %>%
    group_by(Noise, OutlierRate, Misspec, Setting) %>%
    summarise(mean_NRASE = mean(NRASE_Holdout, na.rm = TRUE),
              se_NRASE = se(NRASE_Holdout),
              .groups = "drop") %>%
    arrange(Noise, OutlierRate, Misspec, mean_NRASE)
  print(robust_tbl, n = Inf)

  cat("\n==================== Learning-curve fit (a/sqrt(n)+b) ====================\n")
  lc_tbl <- df %>%
    group_by(n_total, Setting) %>%
    summarise(m = mean(NRASE_Holdout, na.rm = TRUE), .groups = "drop") %>%
    group_by(Setting) %>%
    group_modify(~{
      dd <- .
      if (n_distinct(dd$n_total) < 4 || any(!is.finite(dd$m))) {
        return(tibble(a = NA_real_, b = NA_real_, r2 = NA_real_))
      }
      fit <- try(lm(m ~ I(1/sqrt(n_total)), data = dd), silent = TRUE)
      if (inherits(fit, "try-error")) return(tibble(a = NA,b = NA,r2=NA))
      co <- coef(fit); a <- unname(co[2]); b <- unname(co[1])
      tibble(a = a, b = b, r2 = summary(fit)$r.squared)
    }) %>% ungroup() %>% arrange(Setting)
  print(lc_tbl, n = Inf)

  # Minimal plots (silent; print if desired)
  p_by_n <- ggplot(by_n_tbl, aes(x = n_total, y = mean_NRASE, colour = Setting)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    labs(x = "n_total", y = "Mean NRASE (collapsed)", colour = "Method")

  p_obj <- ggplot(obj_mix %>% filter(!is.na(ObjectiveUsed)),
                  aes(x = n_total, y = p, fill = ObjectiveUsed)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "n_total", y = "AUTO objective share", fill = "ObjectiveUsed")

  invisible(list(
    df = df,
    overall_tbl = overall_tbl,
    by_n_tbl = by_n_tbl,
    auto_vs_sse = auto_vs_sse,
    obj_mix = obj_mix,
    relax_tbl = relax_tbl,
    # NEW: CV_relax summaries
    cv_relax_summary = overall_cv_relax_summary,
    cv_relax_by_n = cv_relax_tbl,
    # NEW: win-rates and regret including CV_relax
    win_rate_all = win_rate_all,
    regret_by_n = regret_by_n,
    theta_4x4 = theta_4x4,
    # Existing head-to-heads
    h2h_relax_tbl = h2h_relax_tbl,
    h2h_std_tbl = h2h_std_tbl,
    robust_tbl = robust_tbl,
    lc_tbl = lc_tbl,
    plots = list(by_n = p_by_n, obj_share = p_obj),
    design = list(grid = grid, scenarios = SC, OUT_ITERS = OUT_ITERS, CORES = CORES, N_BOOT = N_BOOT)
  ))
}

# ----------------- HOW TO RUN -----------------
# out <- run_gap_simulation()                    # light preset (fastest)
# out <- run_gap_simulation(design = "medium")   # bigger grid
# out <- run_gap_simulation(design = "full")     # largest grid
out <- run_gap_simulation(design = "A1")














# ===================== SVEM-ONLY: EXTRA ANALYSIS =====================
# Works directly on the output of run_gap_simulation(...).
# Prints compact tables and returns them in a list you can include in the paper.

analyze_gap_extras <- function(out) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble); library(purrr); library(stringr)
  })
  df <- out$df

  `%||%` <- function(a,b) if (!is.null(a)) a else b
  se <- function(x) sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x)))
  mann_theta <- function(x, y) {
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x)==0 || length(y)==0) return(NA_real_)
    m <- outer(x, y, `-`)
    mean((m < 0) + 0.5*(m == 0), na.rm = TRUE)  # P(X<Y)+0.5P(=)
  }

  # SVEM-only, AUTO only
  svem_auto <- df %>%
    filter(Setting %in% c("SVEM_std","SVEM_relax"), Objective=="auto") %>%
    select(RunID, n_total, TheoreticalR2, Noise, OutlierRate, Misspec,
           Setting, NRASE_Holdout, NAAE_Holdout, ObjectiveUsed) %>%
    distinct()

  # Wide-by-block for direct head-to-heads (SVEM_relax - SVEM_std)
  wide <- svem_auto %>%
    pivot_wider(names_from = Setting, values_from = c(NRASE_Holdout, NAAE_Holdout, ObjectiveUsed)) %>%
    mutate(
      d_NRASE = NRASE_Holdout_SVEM_relax - NRASE_Holdout_SVEM_std,
      d_NAAE  = NAAE_Holdout_SVEM_relax  - NAAE_Holdout_SVEM_std
    )

  cat("\n==================== SVEM RELAX - STD (AUTO): NRASE & NAAE by n_total ====================\n")
  relax_by_n <- wide %>%
    group_by(n_total) %>%
    summarise(
      N         = sum(is.finite(d_NRASE)),
      mean_dNR  = mean(d_NRASE, na.rm=TRUE),
      med_dNR   = median(d_NRASE, na.rm=TRUE),
      wr_relax_le_std_NR = mean(d_NRASE <= 0, na.rm=TRUE),
      mean_dNA  = mean(d_NAAE,  na.rm=TRUE),
      med_dNA   = median(d_NAAE, na.rm=TRUE),
      wr_relax_le_std_NA = mean(d_NAAE <= 0, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(n_total)
  print(relax_by_n, n = Inf)

  cat("\n==================== SVEM RELAX - STD (AUTO): NRASE by R2 and n_total ====================\n")
  relax_by_n_r2 <- wide %>%
    group_by(n_total, TheoreticalR2) %>%
    summarise(
      N       = sum(is.finite(d_NRASE)),
      mean_d  = mean(d_NRASE, na.rm=TRUE),
      med_d   = median(d_NRASE, na.rm=TRUE),
      winrate = mean(d_NRASE <= 0, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(n_total, TheoreticalR2)
  print(relax_by_n_r2, n = Inf)

  cat("\n==================== Conditional on AUTO choice: BIC-only vs AIC-only blocks ====================\n")
  # Only keep blocks where BOTH SVEM variants report the same ObjectiveUsed
  same_obj <- wide %>%
    filter(!is.na(ObjectiveUsed_SVEM_relax), !is.na(ObjectiveUsed_SVEM_std),
           ObjectiveUsed_SVEM_relax == ObjectiveUsed_SVEM_std) %>%
    mutate(Obj = ObjectiveUsed_SVEM_relax)

  cond_tbl <- same_obj %>%
    group_by(Obj, n_total) %>%
    summarise(
      N       = sum(is.finite(d_NRASE)),
      mean_d  = mean(d_NRASE, na.rm=TRUE),
      med_d   = median(d_NRASE, na.rm=TRUE),
      winrate = mean(d_NRASE <= 0, na.rm=TRUE),  # relax better or equal
      .groups="drop"
    ) %>% arrange(Obj, n_total)
  print(cond_tbl, n = Inf)

  cat("\n==================== Mann–Whitney θ = P(NRASE_relax < NRASE_std) by n_total ====================\n")
  theta_by_n <- svem_auto %>%
    select(RunID, n_total, NRASE_Holdout, Setting) %>%
    pivot_wider(names_from = Setting, values_from = NRASE_Holdout) %>%
    group_by(n_total) %>%
    summarise(theta = mann_theta(SVEM_relax, SVEM_std), .groups="drop") %>%
    arrange(n_total)
  print(theta_by_n, n = Inf)

  cat("\n==================== Oracle regret (SVEM only, AUTO) by n_total ====================\n")
  regret_svem <- wide %>%
    transmute(
      RunID, n_total,
      regret_relax = NRASE_Holdout_SVEM_relax -
        pmin(NRASE_Holdout_SVEM_relax, NRASE_Holdout_SVEM_std, na.rm = TRUE),
      regret_std   = NRASE_Holdout_SVEM_std -
        pmin(NRASE_Holdout_SVEM_relax, NRASE_Holdout_SVEM_std, na.rm = TRUE)
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("regret_"),
      names_to = "Setting", values_to = "regret"
    ) %>%
    mutate(Setting = dplyr::recode(Setting,
                                   regret_relax = "SVEM_relax",
                                   regret_std   = "SVEM_std")) %>%
    group_by(n_total, Setting) %>%
    summarise(mean_regret = mean(regret, na.rm = TRUE),
              se = se(regret),
              .groups = "drop") %>%
    arrange(n_total, mean_regret)

  print(regret_svem, n = Inf)


  cat("\n==================== Robustness slices: RELAX advantage (NRASE) by Noise × Outliers × Misspec ====================\n")
  relax_robust <- wide %>%
    group_by(Noise, OutlierRate, Misspec) %>%
    summarise(
      N      = sum(is.finite(d_NRASE)),
      mean_d = mean(d_NRASE, na.rm=TRUE),
      med_d  = median(d_NRASE, na.rm=TRUE),
      winrate= mean(d_NRASE <= 0, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(Noise, OutlierRate, Misspec)
  print(relax_robust, n = Inf)

  cat("\n==================== Objective gating check: P(AIC) ~ n_total + R2 (SVEM only, AUTO) ====================\n")
  obj_df <- svem_auto %>% filter(!is.na(ObjectiveUsed)) %>%
    mutate(AIC = as.integer(ObjectiveUsed == "wAIC"))
  if (nrow(obj_df) > 0) {
    fit <- try(stats::glm(AIC ~ factor(n_total) + TheoreticalR2, family=binomial(), data=obj_df), silent=TRUE)
    if (!inherits(fit,"try-error")) print(summary(fit))
  } else {
    cat("No ObjectiveUsed recorded; cannot fit gating model.\n")
  }

  cat("\n==================== Correlation NRASE vs NAAE (SVEM only, AUTO) ====================\n")
  cor_tbl <- svem_auto %>%
    group_by(Setting) %>%
    summarise(
      N = sum(is.finite(NRASE_Holdout) & is.finite(NAAE_Holdout)),
      cor_NR_NA = cor(NRASE_Holdout, NAAE_Holdout, use="complete.obs"),
      .groups="drop"
    )
  print(cor_tbl)

  cat("\n==================== Margin of victory |d(NRASE)| by n_total (SVEM only) ====================\n")
  mov_tbl <- wide %>%
    mutate(mov = abs(d_NRASE)) %>%
    group_by(n_total) %>%
    summarise(N = sum(is.finite(mov)),
              mean_mov = mean(mov, na.rm=TRUE),
              q50 = median(mov, na.rm=TRUE),
              q90 = quantile(mov, 0.90, na.rm=TRUE),
              .groups="drop")
  print(mov_tbl, n = Inf)

  invisible(list(
    relax_by_n = relax_by_n,
    relax_by_n_r2 = relax_by_n_r2,
    cond_obj_tbl = cond_tbl,
    theta_by_n = theta_by_n,
    regret_svem = regret_svem,
    relax_robust = relax_robust,
    obj_gate_model = if (exists("fit")) fit else NULL,
    cor_tbl = cor_tbl,
    mov_tbl = mov_tbl
  ))
}

