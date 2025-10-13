#' Random-Search Optimizer with Goals, Weights, Optional CIs, and Diverse Candidates
#'
#' Draws random points via \code{svem_random_table_multi}, scores them using user
#' goals and weights, returns the best design point, and optionally proposes
#' \code{k_candidates} diverse high-score candidates by clustering the top fraction
#' of rows with Gower distance and PAM medoids. Medoids are representative
#' existing rows, so proposed candidates are guaranteed to be feasible under all
#' sampling and mixture constraints.
#'
#' @param objects Named list of \code{svem_model} objects (from \code{SVEMnet()}).
#'   List names must match the response names (left-hand sides) of the models.
#' @param goals Named list per response of the form
#'   \code{list(goal = "max"|"min"|"target", weight = nonnegative number, target = number when goal = "target")}.
#'   Weights are normalized to sum to one internally.
#' @param n Number of random samples to draw.
#' @param mixture_groups Optional mixture constraints forwarded to \code{svem_random_table_multi()}.
#'   Each group may include \code{vars}, \code{lower}, \code{upper}, and \code{total}.
#' @param debias Logical; if \code{TRUE}, use debiased predictions for scoring where available.
#' @param agg Aggregation for point predictions, one of \code{"parms"} or \code{"mean"}.
#'   This is passed to \code{predict.svem_model} when applicable.
#' @param ci Logical; if \code{TRUE}, compute percentile confidence intervals when available
#'   via \code{predict(..., interval = TRUE)} or \code{predict_with_ci(...)}.
#' @param level Confidence level for percentile intervals (default \code{0.95}).
#' @param top_frac Fraction in \code{(0, 1]} of highest-score rows to cluster (default \code{0.01}).
#' @param k_candidates Integer number of diverse candidates (medoids) to return (default \code{0}).
#'   If \code{0}, no clustering is performed.
#' @param verbose Logical; print a compact summary of the run and results.
#' @param combine How to combine normalized per-response scores. Use \code{"geom"} for
#'   weighted geometric mean (default) or \code{"mean"} for weighted arithmetic mean.
#' @param geom_floor Small positive floor in \code{(0,1)} used only for \code{combine = "geom"}.
#'   Each normalized score z in \code{[0,1]} is mapped to \code{z_adj = geom_floor + (1 - geom_floor) * z}
#'   before taking logs to avoid \code{log(0)}. Default \code{1e-6}.
#' @param numeric_sampler Sampler for non-mixture numeric predictors passed to
#'   \code{svem_random_table_multi}. One of \code{"random"} (default), \code{"maximin"}, or \code{"uniform"}.
#'   "random" uses \code{lhs::randomLHS} when available, otherwise plain \code{runif}.
#'
#' @return A list with:
#' \describe{
#'   \item{best_idx}{Row index of the selected best design in the sampled table.}
#'   \item{best_x}{Predictors at the best design.}
#'   \item{best_pred}{Named numeric vector of predicted responses at \code{best_x}.}
#'   \item{best_ci}{Data frame of percentile limits when \code{ci = TRUE}; otherwise \code{NULL}.}
#'   \item{candidates}{Data frame of \code{k_candidates} diverse candidates (medoids; existing rows) with
#'         predictors, predictions, optional CIs, and score; \code{NULL} if \code{k_candidates = 0}.}
#'   \item{score_table}{Sampled table with response columns, normalized contributions,
#'         weighted contributions, optional log-weighted terms (when \code{combine = "geom"}), and final score.}
#'   \item{weights}{Normalized weights used in scoring.}
#'   \item{goals}{Tidy data frame describing each response goal, weight, and target.}
#' }
#'
#' @section Scoring:
#' Each response is normalized on the sampled range with \code{normalize01} mapping to \code{[0,1]}:
#' \itemize{
#'   \item \code{goal = "max"} uses \code{normalize01(y)}.
#'   \item \code{goal = "min"} uses \code{normalize01(-y)}.
#'   \item \code{goal = "target"} uses \code{normalize01(-abs(y - target))}.
#' }
#' Then the normalized values are combined across responses:
#' \itemize{
#'   \item \code{combine = "mean"}: weighted arithmetic mean, computed as \code{sum_r w_r * z_r}.
#'   \item \code{combine = "geom"}: weighted geometric mean, computed as
#'         \code{exp( sum_r w_r * log(z_adj_r) )} with \code{z_adj_r = geom_floor + (1 - geom_floor) * z_r}.
#' }
#'
#' @section Diverse candidates:
#' We take the top \code{top_frac} fraction by score, compute Gower distances on predictors,
#' and run PAM to get medoids. Returning medoids rather than centroids ensures each
#' candidate corresponds to an actual sampled setting that satisfies constraints.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 120
#' X1 <- runif(n); X2 <- runif(n)
#' F  <- factor(sample(c("lo","hi"), n, TRUE))
#' y1 <- 1 + 1.5*X1 - 0.8*X2 + 0.4*(F=="hi") + rnorm(n, 0, 0.2)
#' y2 <- 0.7 + 0.4*X1 + 0.4*X2 - 0.2*(F=="hi") + rnorm(n, 0, 0.2)
#' dat <- data.frame(y1, y2, X1, X2, F)
#'
#' m1 <- SVEMnet(y1 ~ X1 + X2 + F, dat, nBoot = 30)
#' m2 <- SVEMnet(y2 ~ X1 + X2 + F, dat, nBoot = 30)
#' objs <- list(y1 = m1, y2 = m2)
#'
#' goals <- list(
#'   y1 = list(goal = "max",    weight = 0.6),
#'   y2 = list(goal = "target", weight = 0.4, target = 0.9)
#' )
#'
#' # Default combine = "geom" and numeric_sampler = "random"
#' out <- svem_optimize_random(
#'   objects      = objs,
#'   goals        = goals,
#'   n            = 5000,
#'   agg          = "mean",
#'   debias       = FALSE,
#'   ci           = TRUE,
#'   level        = 0.95,
#'   k_candidates = 5,
#'   top_frac     = 0.02,
#'   numeric_sampler = "random",
#'   verbose      = TRUE
#' )
#' out$best_x; head(out$candidates)
#'
#' # Use arithmetic mean instead, with fastest uniform sampler
#' out_mean <- svem_optimize_random(
#'   objects      = objs,
#'   goals        = goals,
#'   n            = 5000,
#'   combine      = "mean",
#'   numeric_sampler = "uniform",
#'   verbose      = FALSE
#' )
#'
#' # Mixture-constrained lipid example (composition sums to 1)
#' data(lipid_screen)
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data = lipid_screen, factorial_order = 3
#' )
#' fP <- bigexp_formula(spec, "Potency")
#' fS <- bigexp_formula(spec, "Size")
#' fD <- bigexp_formula(spec, "PDI")
#' mP <- SVEMnet(fP, lipid_screen, nBoot = 40)
#' mS <- SVEMnet(fS, lipid_screen, nBoot = 40)
#' mD <- SVEMnet(fD, lipid_screen, nBoot = 40)
#' objs2 <- list(Potency = mP, Size = mS, PDI = mD)
#'
#' goals2 <- list(
#'   Potency = list(goal = "max", weight = 0.7),
#'   Size    = list(goal = "min", weight = 0.2),
#'   PDI     = list(goal = "min", weight = 0.1)
#' )
#'
#' mixL <- list(list(
#'   vars  = c("Cholesterol","PEG","Ionizable","Helper"),
#'   lower = c(0.10, 0.01, 0.10, 0.10),
#'   upper = c(0.60, 0.05, 0.60, 0.60),
#'   total = 1
#' ))
#'
#' opt <- svem_optimize_random(
#'   objects        = objs2,
#'   goals          = goals2,
#'   n              = 8000,
#'   mixture_groups = mixL,
#'   agg            = "mean",
#'   debias         = FALSE,
#'   ci             = TRUE,
#'   level          = 0.95,
#'   k_candidates   = 5,
#'   top_frac       = 0.01,
#'   numeric_sampler = "random",
#'   verbose        = TRUE
#' )
#' opt$best_x; head(opt$candidates)
#' }
#'
#' @importFrom stats coef model.frame model.matrix na.pass sd
#' @importFrom utils head
#' @importFrom cluster daisy pam
#' @export
svem_optimize_random <- function(objects,
                                 goals,
                                 n = 50000,
                                 mixture_groups = NULL,
                                 debias = FALSE,
                                 agg = c("parms", "mean"),
                                 ci = TRUE,
                                 level = 0.95,
                                 top_frac = 0.01,
                                 k_candidates = 0,
                                 verbose = TRUE,
                                 combine = c("geom", "mean"),
                                 geom_floor = 1e-6,
                                 numeric_sampler = c("random","maximin","uniform")) {
  agg <- match.arg(agg)
  combine <- match.arg(combine)
  numeric_sampler <- match.arg(numeric_sampler)

  # --- validate inputs ---
  if (!is.list(objects) || !length(objects))
    stop("objects must be a nonempty named list of svem_model objects.")
  if (is.null(names(objects)) || any(!nzchar(names(objects))))
    stop("objects must be a named list; names should be the response names.")
  if (!is.list(goals) || !length(goals))
    stop("goals must be a named list.")
  if (is.null(names(goals)) || any(!nzchar(names(goals))))
    stop("goals must be a named list with names matching response names.")
  if (!all(vapply(objects, inherits, logical(1), what = "svem_model")))
    stop("All elements of objects must be svem_model objects.")
  if (any(duplicated(names(objects))))
    stop("Duplicate response names in objects.")

  resp_names <- names(objects)
  if (length(setdiff(resp_names, names(goals)))) {
    stop("Missing goals for responses: ",
         paste(setdiff(resp_names, names(goals)), collapse = ", "))
  }

  # tidy goal table
  goal_df <- data.frame(response = resp_names,
                        goal     = NA_character_,
                        weight   = NA_real_,
                        target   = NA_real_,
                        stringsAsFactors = FALSE)
  for (i in seq_along(resp_names)) {
    r  <- resp_names[i]
    gi <- goals[[r]]
    if (!is.list(gi) || is.null(gi$goal) || is.null(gi$weight))
      stop("Each goals[[response]] must have goal and weight. Offender: ", r)
    g  <- tolower(as.character(gi$goal))
    if (!g %in% c("max", "min", "target"))
      stop("goal for ", r, " must be 'max', 'min', or 'target'.")
    w  <- as.numeric(gi$weight)
    if (!is.finite(w) || w < 0)
      stop("weight for ", r, " must be nonnegative and finite.")
    tval <- if (g == "target") {
      if (is.null(gi$target))
        stop("target must be provided for ", r, " when goal = 'target'.")
      as.numeric(gi$target)
    } else NA_real_
    goal_df$goal[i]   <- g
    goal_df$weight[i] <- w
    goal_df$target[i] <- tval
  }

  # normalize weights
  sw <- sum(goal_df$weight)
  goal_df$weight <- if (sw > 0) goal_df$weight / sw else rep(1 / nrow(goal_df), nrow(goal_df))
  wts <- stats::setNames(goal_df$weight, goal_df$response)

  # --- sample and predict (sampler selectable) ---
  sampled_raw <- svem_random_table_multi(
    objects = objects, n = n, mixture_groups = mixture_groups,
    debias = debias, numeric_sampler = numeric_sampler
  )

  # support both structured and plain data.frame returns
  if (is.list(sampled_raw) && all(c("data", "pred", "all") %in% names(sampled_raw))) {
    sampled <- sampled_raw$all
  } else if (is.data.frame(sampled_raw)) {
    sampled <- sampled_raw
  } else {
    stop("svem_random_table_multi returned an unexpected type.")
  }

  resp_cols <- goal_df$response
  miss <- setdiff(resp_cols, colnames(sampled))
  if (length(miss))
    stop("svem_random_table_multi did not return response columns for: ",
         paste(miss, collapse = ", "))

  # --- scoring ---
  normalize01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] <= rng[1]) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }

  contrib_cols  <- list()
  weighted_cols <- list()
  for (r in resp_cols) {
    g <- goal_df$goal[goal_df$response == r]
    z <- switch(g,
                max    = normalize01(sampled[[r]]),
                min    = normalize01(-sampled[[r]]),
                target = {
                  tval <- goal_df$target[goal_df$response == r]
                  normalize01(-abs(sampled[[r]] - tval))
                }
    )
    contrib_cols[[paste0(r, "_norm")]] <- z
    weighted_cols[[paste0(r, "_w")]]   <- wts[r] * z
  }

  scored <- sampled
  for (nm in names(contrib_cols))  scored[[nm]] <- contrib_cols[[nm]]
  for (nm in names(weighted_cols)) scored[[nm]] <- weighted_cols[[nm]]

  if (combine == "mean") {
    scored$score <- Reduce(`+`, weighted_cols)
  } else {
    # Geometric mean: exp( sum_r w_r * log(z_adj) ), where z_adj = geom_floor + (1 - geom_floor) * z
    log_terms <- list()
    for (r in resp_cols) {
      z  <- contrib_cols[[paste0(r, "_norm")]]
      z_adj <- (1 - geom_floor) * z + geom_floor
      lw <- as.numeric(wts[r]) * log(z_adj)
      log_terms[[r]] <- lw
      scored[[paste0(r, "_logw")]] <- lw
    }
    scored$score <- exp(Reduce(`+`, log_terms))
  }

  # --- best row and predictions ---
  best_idx <- which.max(scored$score)
  predictor_cols <- setdiff(colnames(scored),
                            c(resp_cols, names(contrib_cols), names(weighted_cols), "score",
                              paste0(resp_cols, "_logw")))
  best_x <- scored[best_idx, predictor_cols, drop = FALSE]

  best_pred <- vapply(resp_cols, function(r) {
    pr <- predict(objects[[r]], newdata = best_x, debias = debias, agg = agg)
    if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
  }, numeric(1))
  names(best_pred) <- resp_cols

  # CIs if available
  .get_ci <- function(obj, newx) {
    if (!ci) return(c(lwr = NA_real_, upr = NA_real_))
    res1 <- try(predict(obj, newdata = newx, debias = debias,
                        interval = TRUE, level = level, agg = "mean"), silent = TRUE)
    if (!inherits(res1, "try-error") && is.list(res1) &&
        !is.null(res1$lwr) && !is.null(res1$upr)) {
      return(c(lwr = as.numeric(res1$lwr[1]), upr = as.numeric(res1$upr[1])))
    }
    c(lwr = NA_real_, upr = NA_real_)
  }

  best_ci <- NULL
  if (ci) {
    rows <- lapply(resp_cols, function(r) {
      ci_vals <- .get_ci(objects[[r]], best_x)
      data.frame(response = r, lwr = ci_vals["lwr"], upr = ci_vals["upr"],
                 stringsAsFactors = FALSE)
    })
    best_ci <- do.call(rbind, rows)
  }

  # --- diverse candidates (medoids of top fraction) ---
  candidates <- NULL
  if (k_candidates > 0) {
    if (!(is.numeric(top_frac) && length(top_frac) == 1 && top_frac > 0 && top_frac <= 1))
      stop("top_frac must be in (0,1].")
    m_top <- max(1L, min(nrow(scored), ceiling(top_frac * nrow(scored))))
    ord   <- order(scored$score, decreasing = TRUE)
    top_idx <- ord[seq_len(m_top)]
    top_X   <- scored[top_idx, predictor_cols, drop = FALSE]

    d <- cluster::daisy(top_X, metric = "gower")
    k <- min(k_candidates, m_top)
    pam_fit <- cluster::pam(d, k = k, diss = TRUE)

    med_top_pos    <- pam_fit$id.med
    med_global_idx <- top_idx[med_top_pos]

    build_row <- function(i_idx) {
      xrow <- scored[i_idx, predictor_cols, drop = FALSE]
      preds <- vapply(resp_cols, function(r) {
        pr <- predict(objects[[r]], newdata = xrow, debias = debias, agg = agg)
        if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
      }, numeric(1))
      out <- cbind(xrow, as.data.frame(as.list(preds), optional = TRUE))
      if (ci) {
        for (r in resp_cols) {
          ci_vals <- .get_ci(objects[[r]], xrow)
          out[[paste0(r, "_lwr")]] <- ci_vals["lwr"]
          out[[paste0(r, "_upr")]] <- ci_vals["upr"]
        }
      }
      out$score <- scored$score[i_idx]
      out
    }

    cand_list <- lapply(med_global_idx, build_row)
    candidates <- do.call(rbind, cand_list)
    rownames(candidates) <- NULL
  }

  if (verbose) {
    cat("SVEM random-search optimization\n")
    cat("Points sampled:", n, "\n")
    cat("Aggregation:", agg, "  Debias:", debias, "  CI:", ci,
        if (ci) paste0(" (level=", level, ")"), "\n")
    cat("Combine:", combine,
        if (combine == "geom") paste0("  geom_floor=", format(geom_floor, digits = 6)), "\n")
    cat("Numeric sampler:", numeric_sampler, "\n")
    if (k_candidates > 0)
      cat("Diverse candidates:", k_candidates, "  Top fraction:", top_frac, "\n")
  }

  list(
    best_idx    = best_idx,
    best_x      = best_x,
    best_pred   = best_pred,
    best_ci     = best_ci,
    candidates  = candidates,
    score_table = scored,
    weights     = wts,
    goals       = goal_df
  )
}
