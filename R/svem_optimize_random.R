#' Random-search optimizer with desirabilities, WMT reweighting, CIs,
#' optimal + exploration candidates, and scored originals
#'
#' Draw random points via \code{svem_random_table_multi}, map each response to a
#' desirability in \code{[0,1]} using Derringer–Suich (DS) curves, combine them into
#' a single score, optionally reweight response importances by whole-model test
#' (WMT) p-values, and return: the best design, diverse high-score optimal candidates
#' (PAM medoids of the top fraction), and a second set of exploration candidates that
#' target high predicted uncertainty. Medoids are rows of the sampled table, so all
#' candidates are feasible under sampling and mixture constraints. If \code{data} is
#' provided, the function also returns that same table augmented with per-row predictions,
#' DS terms, the combined score, and an uncertainty measure.
#'
#' @section Binomial handling:
#' For responses fit with \code{family = "binomial"}, this function expects
#' predictions on the probability scale. Predicted fits and percentile CI bounds
#' (when available) are internally clamped to \code{[0,1]} before desirability and
#' uncertainty calculations. To protect current Gaussian behavior, no link-scale
#' transforms are applied. \strong{Reweighting via WMT is not supported when any
#' responses are binomial}; if \code{reweight_by_wmt = TRUE} and at least one response
#' is binomial, the function stops with an informative error.
#'
#' @param objects Named list of \code{svem_model} objects (from \code{SVEMnet()}).
#'   Names are treated as response identifiers (typically matching the left-hand
#'   sides of the model formulas).
#' @param goals Named list per response. Each entry must include: \code{goal} (one of
#'   \code{"max"}, \code{"min"}, \code{"target"}) and a nonnegative \code{weight}. For
#'   \code{goal = "target"}, also provide \code{target}. Optional per-response DS controls:
#'   for \code{"max"}/\code{"min"}: \code{lower_acceptable} (L), \code{upper_acceptable} (U),
#'   \code{shape} (>= 0); for \code{"target"}: \code{tol} (symmetric) or
#'   \code{tol_left}/\code{tol_right}, and \code{shape_left}/\code{shape_right}. If anchors
#'   or tolerances are not provided, robust defaults are inferred from the sampled responses
#'   using the 2nd-98th percentile range.
#' @param data Optional data frame used when \code{reweight_by_wmt = TRUE} and to
#'   produce \code{original_data_scored}. If \code{reweight_by_wmt = TRUE} and
#'   \code{data} is not supplied (or is not a data frame), the function stops.
#'   Each model's stored formula is evaluated on \code{data} for the WMT via
#'   \code{svem_significance_test_parallel()}. Any \code{mixture_groups} are forwarded.
#' @param n Number of random samples to draw.
#' @param mixture_groups Optional mixture constraints forwarded to
#'   \code{svem_random_table_multi()}. Each group may include \code{vars},
#'   \code{lower}, \code{upper}, and \code{total}.
#' @param level Confidence level for percentile intervals. Default: \code{0.95}.
#' @param top_frac Fraction in \code{(0, 1]} of highest-score rows to cluster for
#'   optimal candidates. Default: \code{0.02}.
#' @param k_candidates Number of diverse optimal candidates (medoids) to return.
#'   If \code{0}, no optimality clustering is performed. Default: \code{5}.
#' @param verbose Logical; print a compact summary of the run and results.
#' @param combine How to combine per-response desirabilities. Use \code{"geom"} for
#'   weighted geometric mean (default) or \code{"mean"} for weighted arithmetic mean.
#'   For \code{combine = "geom"}, a small floor is applied before logging to avoid \code{log(0)}.
#' @param numeric_sampler Sampler for non-mixture numeric predictors passed to
#'   \code{svem_random_table_multi()}. One of \code{"random"} (default), \code{"maximin"},
#'   or \code{"uniform"}. \code{"random"} uses \code{lhs::randomLHS()} when available,
#'   otherwise plain \code{runif()}.
#' @param reweight_by_wmt Logical; if \code{TRUE}, compute whole-model p-values (WMT)
#'   for each response on \code{data} and downweight responses with weak factor relationships
#'   before scoring. Requires \code{data}. \strong{Not allowed if any responses are binomial.}
#' @param wmt_transform Transformation used to turn p-values into multipliers when
#'   \code{reweight_by_wmt = TRUE}. One of \code{"neglog10"}, \code{"one_minus_p"}.
#'   Multipliers are floored internally to avoid zeroing weights and then renormalized
#'   to sum to one with the user weights.
#' @param wmt_control Optional list to override WMT defaults passed to
#'   \code{svem_significance_test_parallel()}. Recognized entries:
#'   \code{nPoint}, \code{nSVEM}, \code{nPerm}, \code{percent}, \code{nBoot},
#'   \code{glmnet_alpha}, \code{weight_scheme}, \code{objective},
#'   \code{auto_ratio_cutoff}, \code{relaxed}, \code{verbose}, \code{nCore},
#'   \code{seed}, \code{spec}, \code{response}, \code{use_spec_contrasts}.
#'   By default, \code{svem_optimize_random()} uses \code{wmt_control = list(seed = 123)},
#'   so WMT calls are reproducible; you can override this by passing your own
#'   \code{wmt_control} (with or without a \code{seed}).
#' @param k_exploration_candidates Number of diverse exploration candidates (medoids)
#'   to return. If \code{0}, no exploration clustering is performed. Default: \code{5}.
#' @param exploration_top_frac Fraction in \code{(0, 1]} of rows with the largest
#'   uncertainty measure to cluster for exploration candidates. Default: \code{0.05}.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{best}{One-row data frame at the winning design with predictors,
#'     per-response predictions, per-response percentile CIs (if available),
#'     the combined \code{score}, and \code{uncertainty_measure}.}
#'   \item{best_idx}{Row index of the selected best design in the sampled table.}
#'   \item{best_x}{Predictors at the best design.}
#'   \item{best_pred}{Named numeric vector of predicted responses at \code{best_x}.}
#'   \item{best_ci}{Data frame of percentile limits at \code{best_x}.}
#'   \item{candidates}{Data frame of \code{k_candidates} diverse optimal candidates
#'     (medoids; existing rows) with predictors, predictions, percentile CIs, the combined
#'     \code{score}, and \code{uncertainty_measure}; \code{NULL} if \code{k_candidates = 0}.}
#'   \item{exploration_best}{One-row data frame at the exploration target, with
#'     predictors, per-response predictions, percentile CIs, \code{score}, and
#'     \code{uncertainty_measure}.}
#'   \item{exploration_best_idx}{Row index with the largest \code{uncertainty_measure}.}
#'   \item{exploration_best_x}{Predictors at the exploration target.}
#'   \item{exploration_best_pred}{Predicted responses at \code{exploration_best_x}.}
#'   \item{exploration_best_ci}{Percentile CIs at \code{exploration_best_x}.}
#'   \item{exploration_candidates}{Data frame of \code{k_exploration_candidates} diverse
#'     high-uncertainty candidates (medoids; existing rows) with predictors, predictions,
#'     percentile CIs, \code{uncertainty_measure}, and \code{score}; \code{NULL} if
#'     \code{k_exploration_candidates = 0}.}
#'   \item{score_table}{Sampled table with responses, per-response desirabilities,
#'     weighted terms, optional log-weighted terms (when \code{combine = "geom"}), CI widths,
#'     normalized CI widths, weighted CI widths, the \code{uncertainty_measure}, and final
#'     \code{score}.}
#'   \item{original_data_scored}{If \code{data} is provided: \code{data} augmented with
#'     predicted responses, per-response desirabilities, combined \code{score}, and
#'     \code{uncertainty_measure}. Otherwise \code{NULL}.}
#'   \item{weights_original}{User-provided weights normalized to sum to one before WMT reweighting.}
#'   \item{weights_final}{Final weights after WMT multipliers and renormalization.}
#'   \item{wmt_p_values}{Named vector of per-response whole-model p-values when
#'     \code{reweight_by_wmt = TRUE}; otherwise \code{NULL}.}
#'   \item{wmt_multipliers}{Named vector of per-response WMT multipliers when
#'     \code{reweight_by_wmt = TRUE}; otherwise \code{NULL}.}
#'   \item{goals}{Data frame describing each response goal, weight, target, and
#'     echoing original and final weights and (when applicable) WMT information.}
#' }
#'
#' @details
#' A typical closed-loop workflow for formulation or process optimization is:
#' \enumerate{
#'   \item Fit one or more \code{SVEMnet()} models for responses of interest.
#'   \item Optionally run whole-model testing (WMT) to reweight response goals.
#'   \item Call \code{svem_optimize_random()} to generate:
#'     \itemize{
#'       \item high-scoring "optimal" candidates for immediate testing, and
#'       \item high-uncertainty exploration candidates to improve the models.
#'     }
#'   \item Run these candidates in the lab, append the new data, refit the models,
#'         and repeat as needed.
#' }
#'
#' See the package vignette for a full worked example of this closed-loop workflow.
#'
#' \strong{Multi-response scoring.} Each response is mapped to a DS desirability
#' \eqn{d_r \in [0,1]}. Anchors L and U (and target-band tolerances) default to robust
#' values derived from the sampled 2nd-98th percentile range when not supplied. Desirabilities are
#' combined across responses using either a weighted arithmetic mean (\code{combine = "mean"})
#' or a weighted geometric mean (\code{combine = "geom"}), with a small fixed floor applied
#' inside the log to avoid \code{log(0)}.
#'
#' \strong{Whole-model reweighting (WMT).} When \code{reweight_by_wmt = TRUE}, each response
#' receives a multiplier from its whole-model p-value computed by
#' \code{svem_significance_test_parallel()} on \code{data}. Final weights are proportional to
#' the product of the user weight and the multiplier, then renormalized to sum to one.
#' Supported transforms: \code{"neglog10"} (aggressive) and \code{"one_minus_p"} (conservative).
#' Multipliers are floored internally.
#'
#' \strong{Uncertainty and exploration.} The \code{uncertainty_measure} is the weighted sum of
#' robustly normalized CI widths across responses (each width normalized using the sampled
#' 2nd-98th percentile range, then weighted by the final weights). The largest row is the exploration
#' target; PAM medoids over the top \code{exploration_top_frac} by this measure are returned
#' as exploration candidates. Both optimal and exploration candidate tables include
#' \code{score} and \code{uncertainty_measure}.
#'
#' \emph{Implementation notes.} Point predictions use ensemble-mean aggregation
#' (\code{agg = "mean"}) with \code{debias = FALSE}, both inside
#' \code{svem_random_table_multi()} and in the candidate summaries. Percentile CIs use
#' \code{agg = "mean"}. The geometric combiner uses a fixed floor of \code{1e-6}; the WMT
#' multiplier floor is \code{1e-3}. For binomial responses, fits and CI bounds are
#' clamped to \code{[0,1]}.
#'
#' @template ref-svem
#'
#' @examples
#' \donttest{
#' ## --- Larger Gaussian-only example ---
#' set.seed(1)
#' n  <- 120
#' X1 <- runif(n); X2 <- runif(n)
#' F  <- factor(sample(c("lo","hi"), n, TRUE))
#' y1 <- 1 + 1.5*X1 - 0.8*X2 + 0.4*(F=="hi") + rnorm(n, 0, 0.2)
#' y2 <- 0.7 + 0.4*X1 + 0.4*X2 - 0.2*(F=="hi") + rnorm(n, 0, 0.2)
#' dat <- data.frame(y1, y2, X1, X2, F)
#'
#' m1 <- SVEMnet(y1 ~ X1 + X2 + F, dat, nBoot = 30, family = "gaussian")
#' m2 <- SVEMnet(y2 ~ X1 + X2 + F, dat, nBoot = 30, family = "gaussian")
#' objs <- list(y1 = m1, y2 = m2)
#'
#' goals <- list(
#'   y1 = list(goal = "max",    weight = 0.6),
#'   y2 = list(goal = "target", weight = 0.4, target = 0.9)
#' )
#'
#' out <- svem_optimize_random(
#'   objects      = objs,
#'   goals        = goals,
#'   n            = 5000,
#'   level        = 0.95,
#'   k_candidates = 5,
#'   top_frac     = 0.02,
#'   k_exploration_candidates = 5,
#'   exploration_top_frac     = 0.01,
#'   numeric_sampler = "random",
#'   verbose      = TRUE
#' )
#' out$best
#' head(out$candidates)
#' out$exploration_best
#' head(out$exploration_candidates)
#'
#' ## Optional: reweight goals by whole-model p-values (Gaussian-only).
#' out_wmt <- svem_optimize_random(
#'   objects         = objs,
#'   goals           = goals,
#'   data            = dat,
#'   n               = 5000,
#'   level           = 0.95,
#'   k_candidates    = 5,
#'   top_frac        = 0.02,
#'   k_exploration_candidates = 5,
#'   exploration_top_frac     = 0.01,
#'   numeric_sampler = "random",
#'   reweight_by_wmt = TRUE,
#'   wmt_transform   = "neglog10",
#'   verbose         = TRUE
#' )
#' out_wmt$weights_original
#' out_wmt$weights_final
#' out_wmt$wmt_p_values
#' head(out_wmt$candidates)
#' head(out_wmt$exploration_candidates)
#' head(out_wmt$original_data_scored)
#'
#' ## --- Smaller mixed example: one Gaussian + one Binomial (probability scale) ---
#' set.seed(42)
#' n  <- 80
#' X1 <- runif(n); X2 <- runif(n); G <- factor(sample(c("lo","hi"), n, TRUE))
#'
#' # Gaussian response
#' yg <- 2 + 1.1*X1 - 0.7*X2 + 0.5*(G=="hi") + rnorm(n, 0, 0.25)
#'
#' # Binomial response (probability via logistic link)
#' eta <- -0.3 + 1.2*X1 - 0.4*X2 + 0.6*(G=="hi")
#' p   <- 1/(1 + exp(-eta))
#' yb  <- rbinom(n, 1, p)
#'
#' dmix <- data.frame(yg, yb, X1, X2, G)
#'
#' mg <- SVEMnet(yg ~ X1 + X2 + G, dmix, nBoot = 30, family = "gaussian")
#' mb <- SVEMnet(yb ~ X1 + X2 + G, dmix, nBoot = 30, family = "binomial", relaxed = FALSE)
#'
#' objs_mix <- list(yg = mg, yb = mb)
#' goals_mix <- list(
#'   yg = list(goal = "max", weight = 0.5),
#'   yb = list(goal = "max", weight = 0.5)  # maximize event probability
#' )
#'
#' out_mix <- svem_optimize_random(
#'   objects      = objs_mix,
#'   goals        = goals_mix,
#'   n            = 3000,
#'   level        = 0.95,
#'   k_candidates = 3,
#'   top_frac     = 0.03,
#'   numeric_sampler = "random",
#'   reweight_by_wmt = FALSE,  # required when any response is binomial
#'   verbose      = TRUE
#' )
#' out_mix$best
#' head(out_mix$candidates)
#' }
#'
#' @importFrom stats coef model.frame model.matrix na.pass sd quantile setNames
#' @importFrom utils head
#' @importFrom cluster daisy pam
#' @seealso [SVEMnet()], [svem_random_table_multi()], [predict.svem_model()]
#' @export
svem_optimize_random <- function(objects,
                                 goals,
                                 data = NULL,                 # optional training data for WMT and for appending scored columns
                                 n = 50000,
                                 mixture_groups = NULL,
                                 level = 0.95,
                                 top_frac = 0.02,
                                 k_candidates = 5,
                                 verbose = TRUE,
                                 combine = c("geom", "mean"),
                                 numeric_sampler = c("random","maximin","uniform"),
                                 reweight_by_wmt = FALSE,
                                 wmt_transform = c("neglog10","one_minus_p"),
                                 wmt_control = list(seed = 123),
                                 k_exploration_candidates = 5,
                                 exploration_top_frac = 0.05) {

  # ---------- fixed defaults ----------
  agg_mode      <- "mean"  # hardcoded aggregation for point preds
  debias_flag   <- FALSE    # hardcoded
  geom_floor    <- 1e-6     # hardcoded floor for geometric combine
  wmt_strength  <- 1        # hardcoded transform strength
  wmt_floor     <- 1e-3     # hardcoded multiplier floor

  # robust desirability defaults
  .q_lo         <- 0.02   # robust lower quantile anchor
  .q_hi         <- 0.98   # robust upper quantile anchor
  .tol_frac     <- 0.10   # default target tolerance as fraction of robust span
  .exp_left     <- 1.0    # linear shapes
  .exp_right    <- 1.0

  combine         <- match.arg(combine)
  numeric_sampler <- match.arg(numeric_sampler)
  wmt_transform   <- match.arg(wmt_transform)

  # ---------- validate inputs ----------
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
  miss_goals <- setdiff(resp_names, names(goals))
  if (length(miss_goals))
    stop("Missing goals for responses: ", paste(miss_goals, collapse = ", "))

  # ---- detect families (protect Gaussian behavior; treat binomial on probability scale) ----
  resp_family <- vapply(objects, function(o) {
    fam <- tryCatch(o$family, error = function(e) NA_character_)
    if (!is.character(fam) || !length(fam)) NA_character_ else fam[1]
  }, character(1))
  is_binomial_resp <- tolower(resp_family) %in% "binomial"

  # If any binomial present and WMT requested, do NOT proceed (explicit requirement)
  if (isTRUE(reweight_by_wmt) && any(is_binomial_resp, na.rm = TRUE)) {
    offending <- paste(resp_names[which(is_binomial_resp)], collapse = ", ")
    stop("reweight_by_wmt = TRUE is not supported when any responses are binomial.\n",
         "Binomial responses detected: ", offending,
         "\nPlease disable reweight_by_wmt or run WMT on Gaussian-only models.")
  }

  # ---------- goals & user weights ----------
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
    g <- tolower(as.character(gi$goal))
    if (!g %in% c("max","min","target"))
      stop("goal for ", r, " must be 'max', 'min', or 'target'.")
    w <- as.numeric(gi$weight)
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
  # normalize user weights
  sw <- sum(goal_df$weight)
  weights_user <- if (sw > 0) goal_df$weight / sw else rep(1 / nrow(goal_df), nrow(goal_df))
  names(weights_user) <- goal_df$response

  # ---------- optional whole-model reweighting (Gaussian-only path) ----------
  wmt_p_vals <- setNames(rep(NA_real_, length(resp_names)), resp_names)
  wmt_mult   <- setNames(rep(1.0,      length(resp_names)), resp_names)

  if (isTRUE(reweight_by_wmt)) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("reweight_by_wmt=TRUE requires 'data' (the training data frame) to be supplied to svem_optimize_random().")
    }
    if (isTRUE(verbose)) cat("Whole-model reweighting (WMT): computing p-values using provided data (n=",
                             nrow(data), ")...\n", sep = "")

    # defaults for the test; user can override via wmt_control
    ctrl_def <- list(
      nPoint = 2000, nSVEM = 10, nPerm = 150, percent = 90,
      nBoot = 100, glmnet_alpha = c(1), weight_scheme = "SVEM",
      objective = "auto", auto_ratio_cutoff = 1.3,
      relaxed = FALSE, verbose = verbose,
      nCore = parallel::detectCores(), seed = NULL,
      spec = NULL, response = NULL, use_spec_contrasts = TRUE
    )
    if (!is.null(wmt_control)) {
      nm <- intersect(names(ctrl_def), names(wmt_control))
      for (k in nm) ctrl_def[[k]] <- wmt_control[[k]]
    }

    for (r in resp_names) {
      obj <- objects[[r]]
      # get a formula for this response
      fml <- NULL
      if (!is.null(obj$formula) && inherits(obj$formula, "formula")) {
        fml <- obj$formula
      } else if (!is.null(obj$call) && !is.null(obj$call$formula)) {
        fml <- tryCatch(if (inherits(obj$call$formula, "formula")) obj$call$formula
                        else eval(obj$call$formula, parent.frame()),
                        error = function(e) NULL)
      }
      if (is.null(fml))
        stop("Could not obtain a formula from svem_model '", r, "' for whole-model test. Ensure SVEMnet stored a formula.")

      p <- NA_real_
      res_wmt <- tryCatch({
        do.call(svem_significance_test_parallel,
                c(list(formula = fml,
                       data = data,
                       mixture_groups = mixture_groups),
                  ctrl_def))
      }, error = function(e) {
        if (isTRUE(verbose)) cat(" -", r, ": WMT error ->", conditionMessage(e), "\n")
        NULL
      })
      if (!is.null(res_wmt) && is.finite(res_wmt$p_value)) {
        p <- as.numeric(res_wmt$p_value)
        p <- min(max(p, 1e-16), 1 - 1e-12)  # clamp
      }

      wmt_p_vals[r] <- p

      mult_raw <- if (is.na(p)) 1.0 else switch(
        wmt_transform,
        neglog10    = (-log10(p))^wmt_strength,
        one_minus_p = (1 - p)^wmt_strength
      )
      mult <- max(mult_raw, wmt_floor)
      wmt_mult[r] <- mult

      if (isTRUE(verbose)) {
        cat(sprintf(" - %s: p=%s, multiplier(%s)=%.4g\n",
                    r,
                    if (is.na(p)) "NA" else format(p, digits = 4),
                    wmt_transform, mult))
      }
    }
  }

  # combine user weights with WMT multipliers, then renormalize
  weights_final_raw <- weights_user * wmt_mult
  sw2 <- sum(weights_final_raw)
  weights_final <- if (sw2 > 0 && is.finite(sw2)) weights_final_raw / sw2 else weights_user
  wts <- stats::setNames(weights_final, goal_df$response)

  # ---------- sample & predict ----------
  sampled_raw <- svem_random_table_multi(
    objects = objects, n = n, mixture_groups = mixture_groups,
    debias = debias_flag, numeric_sampler = numeric_sampler
  )
  sampled <- if (is.list(sampled_raw) && all(c("data","pred","all") %in% names(sampled_raw))) {
    sampled_raw$all
  } else if (is.data.frame(sampled_raw)) {
    sampled_raw
  } else stop("svem_random_table_multi returned an unexpected type.")

  resp_cols <- goal_df$response
  if (length(setdiff(resp_cols, colnames(sampled))))
    stop("svem_random_table_multi did not return all response columns.")

  # predictor columns are from the original sampled table
  predictor_cols <- setdiff(colnames(sampled), resp_cols)

  # ---------- DESIRABILITY (Derringer-Suich) ----------
  .q <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
  .clip01 <- function(z) pmin(pmax(z, 0), 1)

  .ds_max <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    z <- ifelse(y <= L, 0, ifelse(y >= U, 1, ((y - L) / (U - L))^s))
    .clip01(z)
  }
  .ds_min <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    z <- ifelse(y <= L, 1, ifelse(y >= U, 0, ((U - y) / (U - L))^s))
    .clip01(z)
  }
  .ds_target <- function(y, T0, L, U, sL = 1, sR = 1) {
    if (!is.finite(T0) || !is.finite(L) || !is.finite(U) || !(L < T0 && T0 < U))
      return(rep(0, length(y)))
    left  <- ifelse(y <  L, 0, ifelse(y <= T0, ((y - L) / (T0 - L))^sL, NA_real_))
    right <- ifelse(y >  U, 0, ifelse(y >= T0, ((U - y) / (U - T0))^sR, NA_real_))
    z <- ifelse(is.na(left), right, ifelse(is.na(right), left, 1))
    .clip01(z)
  }

  # cache desirability parameters so we can reuse on original data
  ds_params <- setNames(vector("list", length(resp_cols)), resp_cols)

  # ---------- scoring (optimality) ----------
  scored <- sampled
  contrib_cols  <- list()
  weighted_cols <- list()

  for (r in resp_cols) {
    gi <- goals[[r]]
    g  <- goal_df$goal[goal_df$response == r]
    y  <- sampled[[r]]

    # For binomial responses, clamp predictions to [0,1] safeguard
    if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
      y <- pmin(pmax(as.numeric(y), 0), 1)
    }

    L_def <- .q(y, .q_lo); U_def <- .q(y, .q_hi)
    if (!is.finite(L_def)) L_def <- min(y, na.rm = TRUE)
    if (!is.finite(U_def)) U_def <- max(y, na.rm = TRUE)
    if (U_def <= L_def) U_def <- L_def + .Machine$double.eps

    if (g == "max") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      z <- .ds_max(y, L, U, s)
      ds_params[[r]] <- list(type="max", L=L, U=U, s=s)
    } else if (g == "min") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      z <- .ds_min(y, L, U, s)
      ds_params[[r]] <- list(type="min", L=L, U=U, s=s)
    } else { # target
      T0 <- as.numeric(goal_df$target[goal_df$response == r])
      tol <- if (!is.null(gi$tol)) as.numeric(gi$tol) else NA_real_
      tol_left  <- if (!is.null(gi$tol_left))  as.numeric(gi$tol_left)  else tol
      tol_right <- if (!is.null(gi$tol_right)) as.numeric(gi$tol_right) else tol
      if (!is.finite(tol_left)  || tol_left  <= 0) tol_left  <- .tol_frac * (U_def - L_def)
      if (!is.finite(tol_right) || tol_right <= 0) tol_right <- .tol_frac * (U_def - L_def)
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else (T0 - tol_left)
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else (T0 + tol_right)
      sL <- if (!is.null(gi$shape_left))  as.numeric(gi$shape_left)  else .exp_left
      sR <- if (!is.null(gi$shape_right)) as.numeric(gi$shape_right) else .exp_right
      if (!(is.finite(L) && is.finite(U) && L < T0 && T0 < U)) {
        band <- .tol_frac * (U_def - L_def)
        L <- T0 - band; U <- T0 + band
      }
      z <- .ds_target(y, T0, L, U, sL, sR)
      ds_params[[r]] <- list(type="target", T0=T0, L=L, U=U, sL=sL, sR=sR)
    }

    contrib_cols[[paste0(r, "_des")]] <- z
    weighted_cols[[paste0(r, "_w")]]  <- as.numeric(wts[r]) * z
  }

  for (nm in names(contrib_cols))  scored[[nm]] <- contrib_cols[[nm]]
  for (nm in names(weighted_cols)) scored[[nm]] <- weighted_cols[[nm]]

  if (combine == "mean") {
    scored$score <- Reduce(`+`, weighted_cols)
  } else {
    for (r in resp_cols) {
      z <- contrib_cols[[paste0(r, "_des")]]
      z_adj <- (1 - geom_floor) * z + geom_floor
      scored[[paste0(r, "_logw")]] <- as.numeric(wts[r]) * log(z_adj)
    }
    scored$score <- exp(Reduce(`+`, lapply(resp_cols, function(r) scored[[paste0(r, "_logw")]])))
  }

  # ---------- uncertainty measure (weighted, robust-normalized CI widths) ----------
  .normalize01_robust <- function(x) {
    a <- .q(x, .q_lo); b <- .q(x, .q_hi)
    if (!is.finite(a) || !is.finite(b) || b <= a) return(rep(0.5, length(x)))
    (pmin(pmax(x, a), b) - a) / (b - a)
  }

  ciw_cols <- list(); ciw_norm_cols <- list(); ciw_w_cols <- list()
  newx <- scored[, predictor_cols, drop = FALSE]
  for (r in resp_cols) {
    obj <- objects[[r]]
    pr <- try(predict(obj, newdata = newx, debias = debias_flag,
                      interval = TRUE, level = level, agg = "mean"),
              silent = TRUE)
    lwr <- upr <- rep(NA_real_, nrow(newx))
    if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
      lwr <- as.numeric(pr$lwr); upr <- as.numeric(pr$upr)
    }
    # For binomial, clamp bounds to [0,1] before computing widths
    if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
      lwr <- pmin(pmax(lwr, 0), 1)
      upr <- pmin(pmax(upr, 0), 1)
    }
    width <- upr - lwr
    ciw_cols[[paste0(r, "_ciw")]] <- width
    zc <- .normalize01_robust(width)
    ciw_norm_cols[[paste0(r, "_ciw_norm")]] <- zc
    ciw_w_cols[[paste0(r, "_ciw_w")]] <- as.numeric(wts[r]) * zc
  }
  for (nm in names(ciw_cols))      scored[[nm]] <- ciw_cols[[nm]]
  for (nm in names(ciw_norm_cols)) scored[[nm]] <- ciw_norm_cols[[nm]]
  for (nm in names(ciw_w_cols))    scored[[nm]] <- ciw_w_cols[[nm]]
  scored$uncertainty_measure <- if (length(ciw_w_cols)) Reduce(`+`, ciw_w_cols) else NA_real_

  # ---------- optimal best + CIs (plus a one-row summary 'best') ----------
  best_idx <- which.max(scored$score)
  best_x   <- scored[best_idx, predictor_cols, drop = FALSE]

  # point predictions at best_x
  best_pred <- vapply(resp_cols, function(r) {
    pr <- predict(objects[[r]], newdata = best_x, debias = debias_flag, agg = agg_mode)
    if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
  }, numeric(1))
  names(best_pred) <- resp_cols

  # percentile CIs at best_x (if available)
  best_ci <- do.call(rbind, lapply(resp_cols, function(r) {
    pr <- try(predict(objects[[r]], newdata = best_x, debias = debias_flag,
                      interval = TRUE, level = level, agg = "mean"),
              silent = TRUE)
    lwr <- upr <- NA_real_
    if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
      lwr <- as.numeric(pr$lwr[1]); upr <- as.numeric(pr$upr[1])
    }
    data.frame(response = r, lwr = lwr, upr = upr, stringsAsFactors = FALSE)
  }))

  # assemble a convenient one-row data.frame like 'candidates'
  best <- best_x
  for (r in resp_cols) {
    best[[r]] <- as.numeric(best_pred[[r]])
  }
  for (i in seq_len(nrow(best_ci))) {
    r   <- best_ci$response[i]
    lwr <- best_ci$lwr[i]
    upr <- best_ci$upr[i]
    best[[paste0(r, "_lwr")]] <- lwr
    best[[paste0(r, "_upr")]] <- upr
  }
  best$score <- scored$score[best_idx]
  best$uncertainty_measure <- scored$uncertainty_measure[best_idx]
  row.names(best) <- NULL


  # ---------- diverse optimal candidates (include uncertainty_measure) ----------
  candidates <- NULL
  if (k_candidates > 0) {
    if (!(is.numeric(top_frac) && length(top_frac) == 1 && top_frac > 0 && top_frac <= 1))
      stop("top_frac must be in (0,1].")
    m_top  <- max(1L, min(nrow(scored), ceiling(top_frac * nrow(scored))))
    ord    <- order(scored$score, decreasing = TRUE)
    top_idx <- ord[seq_len(m_top)]
    top_X   <- scored[top_idx, predictor_cols, drop = FALSE]

    d <- cluster::daisy(top_X, metric = "gower")
    k <- min(k_candidates, m_top)
    pam_fit <- cluster::pam(d, k = k, diss = TRUE)

    med_top_pos    <- pam_fit$id.med
    med_global_idx <- top_idx[med_top_pos]

    candidates <- do.call(rbind, lapply(med_global_idx, function(i_idx) {
      xrow <- scored[i_idx, predictor_cols, drop = FALSE]
      preds <- vapply(resp_cols, function(r) {
        pr <- predict(objects[[r]], newdata = xrow, debias = debias_flag, agg = agg_mode)
        if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
      }, numeric(1))
      out <- cbind(xrow, as.data.frame(as.list(preds), optional = TRUE))
      for (r in resp_cols) {
        pr <- try(predict(objects[[r]], newdata = xrow, debias = debias_flag,
                          interval = TRUE, level = level, agg = "mean"),
                  silent = TRUE)
        lwr <- upr <- NA_real_
        if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
          lwr <- as.numeric(pr$lwr[1]); upr <- as.numeric(pr$upr[1])
          if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
            lwr <- pmin(pmax(lwr, 0), 1)
            upr <- pmin(pmax(upr, 0), 1)
          }
        }
        out[[paste0(r, "_lwr")]] <- lwr
        out[[paste0(r, "_upr")]] <- upr
      }
      out$score <- scored$score[i_idx]
      out$uncertainty_measure <- scored$uncertainty_measure[i_idx]
      out
    }))
    rownames(candidates) <- NULL
  }

  # ---------- exploration target + diverse exploration candidates ----------
  exploration_best_idx <- if (any(is.finite(scored$uncertainty_measure))) which.max(scored$uncertainty_measure) else NA_integer_
  exploration_best_x   <- if (is.finite(exploration_best_idx)) scored[exploration_best_idx, predictor_cols, drop = FALSE] else NULL

  exploration_best_pred <- NULL
  exploration_best_ci   <- NULL
  if (!is.null(exploration_best_x)) {
    exploration_best_pred <- vapply(resp_cols, function(r) {
      pr <- predict(objects[[r]], newdata = exploration_best_x, debias = debias_flag, agg = agg_mode)
      if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
    }, numeric(1))
    names(exploration_best_pred) <- resp_cols
    exploration_best_ci <- do.call(rbind, lapply(resp_cols, function(r) {
      pr <- try(predict(objects[[r]], newdata = exploration_best_x, debias = debias_flag,
                        interval = TRUE, level = level, agg = "mean"),
                silent = TRUE)
      lwr <- upr <- NA_real_
      if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
        lwr <- as.numeric(pr$lwr[1]); upr <- as.numeric(pr$upr[1])
        if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
          lwr <- pmin(pmax(lwr, 0), 1)
          upr <- pmin(pmax(upr, 0), 1)
        }
      }
      data.frame(response = r, lwr = lwr, upr = upr, stringsAsFactors = FALSE)
    }))
  }

  # assemble a convenient one-row data.frame like 'candidates' for exploration
  exploration_best <- NULL
  if (!is.null(exploration_best_x)) {
    exploration_best <- exploration_best_x

    # add point predictions
    for (r in resp_cols) {
      exploration_best[[r]] <- as.numeric(exploration_best_pred[[r]])
    }

    # add CIs (respect binomial clamping if applicable)
    for (i in seq_len(nrow(exploration_best_ci))) {
      r   <- exploration_best_ci$response[i]
      lwr <- exploration_best_ci$lwr[i]
      upr <- exploration_best_ci$upr[i]
      # clamp to [0,1] for binomial responses
      if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
        lwr <- pmin(pmax(lwr, 0), 1)
        upr <- pmin(pmax(upr, 0), 1)
      }
      exploration_best[[paste0(r, "_lwr")]] <- lwr
      exploration_best[[paste0(r, "_upr")]] <- upr
    }

    # add score and uncertainty for that row
    exploration_best$score <- scored$score[exploration_best_idx]
    exploration_best$uncertainty_measure <- scored$uncertainty_measure[exploration_best_idx]
    row.names(exploration_best) <- NULL
  }


  exploration_candidates <- NULL
  if (!is.null(exploration_best_x) && k_exploration_candidates > 0) {
    if (!(is.numeric(exploration_top_frac) && length(exploration_top_frac) == 1 &&
          exploration_top_frac > 0 && exploration_top_frac <= 1))
      stop("exploration_top_frac must be in (0,1].")
    m_top <- max(1L, min(nrow(scored), ceiling(exploration_top_frac * nrow(scored))))
    ord   <- order(scored$uncertainty_measure, decreasing = TRUE)
    top_idx <- ord[seq_len(m_top)]
    top_X   <- scored[top_idx, predictor_cols, drop = FALSE]

    d <- cluster::daisy(top_X, metric = "gower")
    k <- min(k_exploration_candidates, m_top)
    pam_fit <- cluster::pam(d, k = k, diss = TRUE)

    med_top_pos    <- pam_fit$id.med
    med_global_idx <- top_idx[med_top_pos]

    exploration_candidates <- do.call(rbind, lapply(med_global_idx, function(i_idx) {
      xrow <- scored[i_idx, predictor_cols, drop = FALSE]
      preds <- vapply(resp_cols, function(r) {
        pr <- predict(objects[[r]], newdata = xrow, debias = debias_flag, agg = agg_mode)
        if (is.list(pr) && !is.null(pr$fit)) pr$fit[1] else as.numeric(pr[1])
      }, numeric(1))
      out <- cbind(xrow, as.data.frame(as.list(preds), optional = TRUE))
      for (r in resp_cols) {
        pr <- try(predict(objects[[r]], newdata = xrow, debias = debias_flag,
                          interval = TRUE, level = level, agg = "mean"),
                  silent = TRUE)
        lwr <- upr <- NA_real_
        if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
          lwr <- as.numeric(pr$lwr[1]); upr <- as.numeric(pr$upr[1])
          if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
            lwr <- pmin(pmax(lwr, 0), 1)
            upr <- pmin(pmax(upr, 0), 1)
          }
        }
        out[[paste0(r, "_lwr")]] <- lwr
        out[[paste0(r, "_upr")]] <- upr
      }
      out$uncertainty_measure <- scored$uncertainty_measure[i_idx]
      out$score               <- scored$score[i_idx]
      out
    }))
    rownames(exploration_candidates) <- NULL
  }

  # ---------- optional: score + uncertainty on the ORIGINAL data ----------
  original_data_scored <- NULL
  if (!is.null(data) && is.data.frame(data)) {
    # (1) predictions per response on 'data'
    pred_df <- data.frame(row_id = seq_len(nrow(data)))
    ci_lwr  <- list(); ci_upr <- list()

    for (r in resp_names) {
      pr <- try(predict(objects[[r]], newdata = data, debias = debias_flag,
                        interval = TRUE, level = level, agg = "mean"),
                silent = TRUE)
      fit_col <- if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$fit)) {
        as.numeric(pr$fit)
      } else {
        as.numeric(try(predict(objects[[r]], newdata = data, debias = debias_flag, agg = agg_mode), silent = TRUE))
      }
      # clamp probability fits for binomial responses
      if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
        fit_col <- pmin(pmax(fit_col, 0), 1)
      }
      pred_df[[r]] <- fit_col

      if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
        lwr <- as.numeric(pr$lwr); upr <- as.numeric(pr$upr)
        if (isTRUE(is_binomial_resp[which(resp_names == r)])) {
          lwr <- pmin(pmax(lwr, 0), 1)
          upr <- pmin(pmax(upr, 0), 1)
        }
        ci_lwr[[r]] <- lwr
        ci_upr[[r]] <- upr
      } else {
        ci_lwr[[r]] <- rep(NA_real_, nrow(data))
        ci_upr[[r]] <- rep(NA_real_, nrow(data))
      }
    }

    # (2) desirability using the SAME ds_params learned on the sampled table
    des_cols <- list()
    for (r in resp_names) {
      y <- pred_df[[r]]
      p <- ds_params[[r]]
      if (p$type == "max") {
        des_cols[[paste0(r, "_des")]] <- .ds_max(y, p$L, p$U, p$s)
      } else if (p$type == "min") {
        des_cols[[paste0(r, "_des")]] <- .ds_min(y, p$L, p$U, p$s)
      } else { # target
        des_cols[[paste0(r, "_des")]] <- .ds_target(y, p$T0, p$L, p$U, p$sL, p$sR)
      }
    }

    # (3) combine into score
    if (combine == "mean") {
      score_vec <- Reduce(`+`, lapply(resp_names, function(r) as.numeric(wts[r]) * des_cols[[paste0(r, "_des")]]))
    } else {
      logsum <- Reduce(`+`, lapply(resp_names, function(r) {
        z  <- des_cols[[paste0(r, "_des")]]
        zf <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts[r]) * log(zf)
      }))
      score_vec <- exp(logsum)
    }

    # (4) uncertainty_measure for original data (robust-normalized CI widths)
    unc_w <- list()
    for (r in resp_names) {
      width <- ci_upr[[r]] - ci_lwr[[r]]
      a <- .q(width, .q_lo); b <- .q(width, .q_hi)
      if (!is.finite(a) || !is.finite(b) || b <= a) {
        zc <- rep(0.5, length(width))
      } else {
        zc <- (pmin(pmax(width, a), b) - a) / (b - a)
      }
      unc_w[[r]] <- as.numeric(wts[r]) * zc
    }
    unc_vec <- Reduce(`+`, unc_w)

    # (5) assemble: original data + predictions + desirabilities + score + uncertainty

    # Rename prediction columns with _pred suffix to avoid name clashes
    pred_mat <- pred_df[resp_names]
    colnames(pred_mat) <- paste0(resp_names, "_pred")

    original_data_scored <- cbind(
      data,
      pred_mat,
      setNames(as.data.frame(des_cols, stringsAsFactors = FALSE), names(des_cols))
    )
    original_data_scored$score <- score_vec
    original_data_scored$uncertainty_measure <- unc_vec

  }

  # ---------- logging ----------
  if (isTRUE(verbose)) {
    cat("SVEM random-search optimization\n")
    cat("Points sampled:", n, "\n")
    cat("Percentile CI level:", level, "\n")
    cat("Combine:", combine,
        if (combine == "geom") paste0("  geom_floor=", format(geom_floor, digits = 6)), "\n")
    cat("Numeric sampler:", numeric_sampler, "\n")
    cat("Desirability mapping: Derringer-Suich with robust anchors (2%-98%)",
        " unless per-response anchors/tolerances are provided in goals.\n", sep = "")
    if (any(is_binomial_resp)) {
      cat("Detected binomial responses: using probability scale with clamped predictions/CIs in [0,1].\n")
    }
    if (isTRUE(reweight_by_wmt)) {
      cat("WMT reweighting: transform=", wmt_transform,
          " strength=1 floor=1e-03\n", sep = "")
    }
    if (k_candidates > 0)
      cat("Diverse optimal candidates:", k_candidates, "  Top fraction:", top_frac, "\n")
    if (k_exploration_candidates > 0)
      cat("Diverse exploration candidates:", k_exploration_candidates,
          "  Top fraction:", exploration_top_frac, "\n")
  }

  # ---------- output ----------
  list(
    # optimality track
    best_idx              = best_idx,
    best_x                = best_x,
    best_pred             = best_pred,
    best_ci               = best_ci,
    best                  = best,
    candidates            = candidates,                 # includes score + uncertainty_measure

    # exploration track
    exploration_best      = exploration_best,
    exploration_best_idx  = exploration_best_idx,
    exploration_best_x    = exploration_best_x,
    exploration_best_pred = exploration_best_pred,
    exploration_best_ci   = exploration_best_ci,
    exploration_candidates= exploration_candidates,     # includes uncertainty_measure + score

    # full sampled table with all diagnostics
    score_table           = scored,                     # includes per-response desirabilities and uncertainty_measure

    # (optional) scored original data
    original_data_scored  = original_data_scored,

    # weights
    weights_original      = weights_user,               # user-normalized
    weights_final         = wts,                        # after WMT multipliers + renormalization
    wmt_p_values          = if (isTRUE(reweight_by_wmt)) wmt_p_vals else NULL,
    wmt_multipliers       = if (isTRUE(reweight_by_wmt)) wmt_mult   else NULL,

    # goals echo
    goals                 = within(goal_df, {
      weight_original <- weights_user[response]
      weight_final    <- wts[response]
      wmt_p           <- if (isTRUE(reweight_by_wmt)) wmt_p_vals[response] else NA_real_
      wmt_mult_col    <- if (isTRUE(reweight_by_wmt)) wmt_mult[response]   else 1.0
    })
  )
}
