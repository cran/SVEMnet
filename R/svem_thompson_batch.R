#' Propose a batch of runs by Thompson sampling from SVEM bootstrap ensembles
#'
#' @description
#' Experimental sequential-design proposer for multi-response SVEM workflows.
#' This extends the core workflow described in Karl (2026),
#' \doi{10.1016/j.chemolab.2026.105660}; its performance is not established
#' by the results in that package paper. Given fitted
#' \code{\link{SVEMnet}} models and per-response desirability goals,
#' \code{svem_thompson_batch()} proposes \code{batch_size} new candidate runs
#' by parallel Thompson sampling (Thompson, 1933; Kandasamy et al., 2018)
#' from the SVEM bootstrap ensembles: for each batch slot, one bootstrap
#' member is drawn at random for each response, the drawn members' predictions
#' are mapped to Derringer--Suich desirabilities and combined into a single
#' score, and the candidate maximizing that score is added to the batch
#' (previously selected candidates are excluded).
#'
#' Each slot re-draws the bootstrap members, allowing candidates preferred
#' by different members to enter the batch without an explicit exploration
#' parameter. Bootstrap members are used as heuristic function draws; they
#' are not samples from a fitted Bayesian posterior. This complements the
#' score/uncertainty workflow of \code{\link{svem_score_random}} and
#' \code{\link{svem_select_from_score_table}}: CI-width exploration targets
#' regions of large global uncertainty, while Thompson sampling refines the
#' current contenders for the optimum.
#'
#' @details
#' \subsection{Candidate table}{
#' Candidates may be supplied directly via \code{candidates} (a data frame
#' containing at least the predictor columns used by the models). When
#' \code{candidates = NULL}, \code{n} random candidates are drawn from the
#' models' shared sampling schema with
#' \code{\link{svem_random_table_multi}}, honoring \code{mixture_groups} and
#' \code{numeric_sampler} exactly as in \code{\link{svem_score_random}}.
#' }
#'
#' \subsection{Goals, desirabilities, and combination}{
#' \code{goals} uses the same specification as
#' \code{\link{svem_score_random}}: each response gets \code{goal}
#' (\code{"max"}, \code{"min"}, or \code{"target"}), a nonnegative
#' \code{weight}, and optional Derringer--Suich anchors
#' (\code{lower_acceptable}/\code{upper_acceptable}/\code{shape}, or
#' \code{target}/\code{tol}/\code{tol_left}/\code{tol_right}/
#' \code{shape_left}/\code{shape_right}). When anchors are not supplied,
#' robust defaults are inferred from the ensemble-mean predictions over the
#' candidate table using the q0.02--q0.98 span, mirroring
#' \code{\link{svem_score_random}}. Fixing anchors explicitly is recommended
#' for sequential use so that the desirability scale does not drift between
#' iterations. Desirabilities are combined by a weighted geometric mean
#' (\code{combine = "geom"}, with a small floor inside the log) or a weighted
#' arithmetic mean (\code{combine = "mean"}); weights are normalized to sum
#' to one.
#' }
#'
#' \subsection{Thompson draws}{
#' For each response the full member-by-candidate prediction matrix is formed
#' once from \code{object$coef_matrix} (all models must be fitted with
#' \code{nBoot >= 2} so that \code{coef_matrix} is available). For slot
#' \eqn{b = 1, \dots,} \code{batch_size} and response \eqn{r}, a member index
#' \eqn{m_{b,r}} is drawn uniformly at random and the slot score at candidate
#' \eqn{x} is
#' \deqn{ s_b(x) = \prod_r d_r\{\hat y^{(m_{b,r})}_r(x)\}^{w_r} }
#' (or the weighted arithmetic mean when \code{combine = "mean"}). The
#' proposal for slot \eqn{b} is the eligible candidate maximizing
#' \eqn{s_b(x)}; ties keep the earliest candidate row. Draws are controlled by
#' the R random number generator, so use \code{set.seed()} for
#' reproducibility.
#'
#' Predictions are always computed with debiasing disabled
#' (\code{debias = FALSE}), matching \code{\link{svem_score_random}}.
#' Binomial responses are handled on the probability scale.
#' }
#'
#' \subsection{Typical sequential workflow}{
#' \enumerate{
#'   \item Fit one \code{SVEMnet()} model per response on the current data.
#'   \item Call \code{svem_thompson_batch()} to propose the next batch of
#'         runs (optionally reserving some slots for medoids from
#'         \code{\link{svem_select_from_score_table}}, e.g. top-\code{score}
#'         medoids for exploitation or widest-CI medoids for global
#'         gap-filling).
#'   \item Run the proposed settings, append the new data, refit, repeat.
#' }
#' }
#'
#' @param objects List of \code{svem_model} objects (from
#'   \code{\link{SVEMnet}}), one per response, as in
#'   \code{\link{svem_score_random}}. When unnamed, response names are
#'   inferred from the model formulas. All models must have bootstrap
#'   coefficients (\code{coef_matrix}; fit with \code{nBoot >= 2}) and, when
#'   \code{candidates = NULL}, a shared sampling schema.
#' @param goals List of per-response goal specifications in the same format
#'   as \code{\link{svem_score_random}} (see Details).
#' @param batch_size Integer number of runs to propose. Default \code{4}.
#' @param candidates Optional data frame of candidate settings containing the
#'   predictor columns required by the models. When \code{NULL} (default),
#'   \code{n} candidates are drawn with \code{\link{svem_random_table_multi}}.
#' @param n Number of random candidates to draw when \code{candidates} is
#'   \code{NULL}. Default \code{5000}.
#' @param mixture_groups Optional mixture and simplex constraints passed to
#'   \code{\link{svem_random_table_multi}} (ignored when \code{candidates}
#'   is supplied).
#' @param numeric_sampler Character; how numeric predictors are sampled when
#'   \code{candidates = NULL}: \code{"random"} (Latin hypercube when
#'   \pkg{lhs} is available) or \code{"uniform"}. As in
#'   \code{\link{svem_score_random}}.
#' @param combine How to combine per-response desirabilities into the slot
#'   score: \code{"geom"} (weighted geometric mean; default) or
#'   \code{"mean"} (weighted arithmetic mean).
#' @param verbose Logical; if \code{TRUE} (default), print a compact summary
#'   of the proposed batch.
#'
#' @return
#' An object of class \code{"svem_thompson_batch"}: a list with components
#' \describe{
#'   \item{\code{proposals}}{Data frame with \code{batch_size} rows: the
#'     selected candidate settings (all columns of the candidate table), plus
#'     \code{slot} (batch position), \code{ts_score} (the slot's Thompson
#'     score at the selected candidate), \code{<resp>_draw} (the drawn
#'     member's prediction at the selected candidate),
#'     \code{<resp>_des} (its desirability), and \code{<resp>_pred}
#'     (the ensemble-mean prediction, for reference).}
#'   \item{\code{members}}{Integer matrix (\code{batch_size} by number of
#'     responses) of the bootstrap member indices drawn for each slot.}
#'   \item{\code{candidate_index}}{Integer vector: row indices of the
#'     proposals in the candidate table.}
#'   \item{\code{candidates}}{The candidate table used (predictor columns).}
#'   \item{\code{weights}}{Normalized response weights.}
#'   \item{\code{ds_params}}{Per-response Derringer--Suich parameters
#'     actually used (after anchor inference), useful for keeping anchors
#'     fixed across sequential iterations.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @references
#' Thompson, W. R. (1933). On the likelihood that one unknown probability
#' exceeds another in view of the evidence of two samples.
#' \emph{Biometrika}, 25(3--4), 285--294. \doi{10.1093/biomet/25.3-4.285}
#'
#' Kandasamy, K., Krishnamurthy, A., Schneider, J., and Poczos, B. (2018).
#' Parallelised Bayesian optimisation via Thompson sampling.
#' \emph{Proceedings of the 21st International Conference on Artificial
#' Intelligence and Statistics (AISTATS)}, PMLR 84, 133--142.
#'
#' Russo, D. J., Van Roy, B., Kazerouni, A., Osband, I., and Wen, Z. (2018).
#' A tutorial on Thompson sampling. \emph{Foundations and Trends in Machine
#' Learning}, 11(1), 1--96. \doi{10.1561/2200000070}
#'
#' Derringer, G., and Suich, R. (1980). Simultaneous optimization of several
#' response variables. \emph{Journal of Quality Technology}, 12(4), 214--219.
#' \doi{10.1080/00224065.1980.11980968}
#'
#' @seealso
#' \code{\link{SVEMnet}}, \code{\link{svem_score_random}},
#' \code{\link{svem_select_from_score_table}},
#' \code{\link{svem_random_table_multi}}
#' @template ref-svem
#'
#' @examples
#' ## Small runnable two-response proposal
#' set.seed(21)
#' n_toy <- 36
#' toy <- data.frame(x1 = runif(n_toy), x2 = runif(n_toy))
#' toy$yield <- 1 + 1.5 * toy$x1 - toy$x2 + rnorm(n_toy, sd = 0.15)
#' toy$impurity <- 1.2 - toy$x1 + 0.4 * toy$x2 + rnorm(n_toy, sd = 0.10)
#'
#' fit_yield <- SVEMnet(yield ~ x1 + x2, toy, nBoot = 5,
#'                      glmnet_alpha = 1, relaxed = FALSE)
#' fit_impurity <- SVEMnet(impurity ~ x1 + x2, toy, nBoot = 5,
#'                         glmnet_alpha = 1, relaxed = FALSE)
#'
#' toy_prop <- svem_thompson_batch(
#'   objects = list(yield = fit_yield, impurity = fit_impurity),
#'   goals = list(
#'     yield = list(goal = "max", weight = 0.7),
#'     impurity = list(goal = "min", weight = 0.3)
#'   ),
#'   batch_size = 2,
#'   n = 100,
#'   verbose = FALSE
#' )
#' toy_prop$proposals
#'
#' \dontrun{
#' ## Production-scale lipid workflow. This is not run automatically because
#' ## it fits wide ensembles and searches thousands of constrained candidates.
#' ## ------------------------------------------------------------------------
#' ## Thompson-sampling batch proposals for the lipid screen
#' ## ------------------------------------------------------------------------
#'
#' data(lipid_screen)
#'
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 2,
#'   polynomial_order = 2
#' )
#'
#' set.seed(1)
#' fit_pot <- SVEMnet(bigexp_formula(spec, "Potency"), lipid_screen)
#' fit_siz <- SVEMnet(bigexp_formula(spec, "Size"),    lipid_screen)
#'
#' objs <- list(Potency = fit_pot, Size = fit_siz)
#'
#' goals <- list(
#'   Potency = list(goal = "max", weight = 0.7),
#'   Size    = list(goal = "min", weight = 0.3)
#' )
#'
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' set.seed(2)
#' prop <- svem_thompson_batch(
#'   objects        = objs,
#'   goals          = goals,
#'   batch_size     = 4,
#'   n              = 3000,
#'   mixture_groups = mix
#' )
#'
#' prop$proposals[, c("PEG", "Helper", "Ionizable", "Cholesterol",
#'                    "ts_score", "Potency_draw", "Size_draw")]
#'
#' ## Keep the inferred desirability anchors fixed at later iterations by
#' ## passing explicit anchors extracted from prop$ds_params, e.g.:
#' p <- prop$ds_params$Potency
#' goals$Potency <- list(goal = "max", weight = 0.7,
#'                       lower_acceptable = p$L, upper_acceptable = p$U)
#'
#' ## ------------------------------------------------------------------------
#' ## Hybrid batch: exploit medoids + Thompson-sampling slots
#' ## ------------------------------------------------------------------------
#'
#' set.seed(3)
#' scored <- svem_score_random(objs, goals, n = 3000,
#'                             mixture_groups = mix, verbose = FALSE)
#' sel <- svem_select_from_score_table(scored$score_table,
#'                                     target = "score", k = 2)
#'
#' set.seed(4)
#' ts <- svem_thompson_batch(objs, goals, batch_size = 2,
#'                           n = 3000, mixture_groups = mix,
#'                           verbose = FALSE)
#'
#' pred_cols <- attr(scored$score_table, "svem_predictor_cols")
#' batch <- rbind(sel$candidates[, pred_cols],
#'                ts$proposals[, pred_cols])
#' batch
#' }
#'
#' @export
svem_thompson_batch <- function(objects,
                                goals,
                                batch_size = 4,
                                candidates = NULL,
                                n = 5000,
                                mixture_groups = NULL,
                                numeric_sampler = c("random", "uniform"),
                                combine = c("geom", "mean"),
                                verbose = TRUE) {

  mc <- match.call()

  # ---- constants shared with svem_score_random ----
  geom_floor         <- 1e-6
  .q_lo              <- 0.02
  .q_hi              <- 0.98
  .tol_frac          <- 0.10
  .ds_span_abs_floor <- 1e-6
  .ds_span_rel_floor <- 0.05

  combine         <- match.arg(combine)
  numeric_sampler <- match.arg(numeric_sampler)

  # ---- validate objects ----
  if (inherits(objects, "svem_model")) objects <- list(objects)
  if (!is.list(objects) || !length(objects))
    stop("objects must be a nonempty list of svem_model objects.")
  if (!all(vapply(objects, inherits, logical(1), what = "svem_model")))
    stop("All elements of 'objects' must be svem_model objects.")

  batch_size <- .svem_integer_scalar(batch_size, "batch_size", min = 1L)
  verbose <- .svem_logical_scalar(verbose, "verbose")

  # ---- infer response names (same rules as svem_score_random) ----
  lhs_names <- vapply(objects, function(o) {
    if (!is.null(o$formula) && inherits(o$formula, "formula")) {
      .svem_response_label(o$formula)
    } else NA_character_
  }, character(1))
  valid_lhs <- !is.na(lhs_names) & nzchar(lhs_names)
  if (anyDuplicated(lhs_names[valid_lhs])) {
    dups <- unique(lhs_names[valid_lhs][duplicated(lhs_names[valid_lhs])])
    stop("Duplicate model response (LHS) names detected: ",
         paste(dups, collapse = ", "),
         ". Please rename the response variables and/or name the objects list uniquely.")
  }
  obj_names <- names(objects)
  if (is.null(obj_names)) obj_names <- rep("", length(objects))
  obj_names[is.na(obj_names)] <- ""
  empty_idx <- which(!nzchar(obj_names) & valid_lhs)
  if (length(empty_idx)) obj_names[empty_idx] <- lhs_names[empty_idx]
  names(objects) <- obj_names
  resp_names <- names(objects)
  if (any(!nzchar(resp_names)))
    stop("Could not infer response names for all models; please name the objects ",
         "list or ensure each model has a formula with a left-hand side.")
  if (anyDuplicated(resp_names))
    stop("Response identifiers must be unique. Duplicates: ",
         paste(unique(resp_names[duplicated(resp_names)]), collapse = ", "))

  key_to_resp <- .svem_response_key_map(resp_names, lhs_names, valid_lhs)

  # ---- require bootstrap members for every model ----
  has_members <- vapply(objects, function(o) {
    !is.null(o$coef_matrix) && is.matrix(o$coef_matrix) && nrow(o$coef_matrix) >= 2L
  }, logical(1))
  if (!all(has_members)) {
    stop("Thompson sampling requires bootstrap member coefficients ",
         "(coef_matrix with nBoot >= 2) for every model. Missing for: ",
         paste(resp_names[!has_members], collapse = ", "))
  }

  # ---- goals: unnamed positional or named by object name / LHS ----
  goals <- .svem_align_goals(goals, resp_names, key_to_resp)
  goal_info <- .svem_validate_goals(goals, resp_names)
  goals <- goal_info$goals
  goal_type <- goal_info$goal
  wts <- goal_info$normalized_weight

  # ---- candidate table ----
  if (is.null(candidates)) {
    n <- .svem_integer_scalar(n, "n", min = 1L)
    sampled <- svem_random_table_multi(
      objects         = objects,
      n               = n,
      mixture_groups  = mixture_groups,
      debias          = FALSE,
      numeric_sampler = numeric_sampler
    )
    cand_df <- as.data.frame(sampled$data)
  } else {
    if (!is.data.frame(candidates))
      stop("'candidates' must be NULL or a data.frame.")
    if (!nrow(candidates))
      stop("'candidates' has zero rows.")
    cand_df <- as.data.frame(candidates)
  }
  n_cand <- nrow(cand_df)
  if (batch_size > n_cand)
    stop("batch_size (", batch_size, ") exceeds the number of candidates (",
         n_cand, ").")

  # ---- member-level predictions per response ----
  # eta_list[[r]]: n_cand x nBoot_r matrix of member predictions on the
  # response scale (probability scale for binomial responses), via the
  # shared internal helper .svem_member_predictions().
  eta_list  <- stats::setNames(vector("list", length(resp_names)), resp_names)
  mean_pred <- stats::setNames(vector("list", length(resp_names)), resp_names)
  bad_any   <- rep(FALSE, n_cand)

  for (r in resp_names) {
    mp <- .svem_member_predictions(objects[[r]], cand_df, type = "response")
    eta_list[[r]]  <- mp$pred_mat
    mean_pred[[r]] <- rowMeans(mp$pred_mat)
    bad_any <- bad_any | mp$bad_rows
  }
  if (all(bad_any))
    stop("All candidate rows contain unseen or missing levels for at least ",
         "one model; no eligible candidates remain.")

  # ---- Derringer-Suich desirability functions (as in svem_score_random) ----
  .q      <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  .clip01 <- function(z) pmin(pmax(z, 0), 1)
  .ds_max <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    .clip01(ifelse(y <= L, 0, ifelse(y >= U, 1, ((y - L) / (U - L))^s)))
  }
  .ds_min <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    .clip01(ifelse(y <= L, 1, ifelse(y >= U, 0, ((U - y) / (U - L))^s)))
  }
  .ds_target <- function(y, T0, L, U, sL = 1, sR = 1) {
    if (!is.finite(T0) || !is.finite(L) || !is.finite(U) || !(L < T0 && T0 < U))
      return(rep(0, length(y)))
    left  <- ifelse(y <  L, 0, ifelse(y <= T0, ((y - L) / (T0 - L))^sL, NA_real_))
    right <- ifelse(y >  U, 0, ifelse(y >= T0, ((U - y) / (U - T0))^sR, NA_real_))
    .clip01(ifelse(is.na(left), right, ifelse(is.na(right), left, 1)))
  }

  # ---- resolve DS parameters per response (user anchors or robust defaults
  #      from ensemble-mean predictions over the candidate table) ----
  ds_params <- stats::setNames(vector("list", length(resp_names)), resp_names)
  for (r in resp_names) {
    gi <- goals[[r]]
    g  <- goal_type[[r]]
    y  <- mean_pred[[r]]

    L_def <- .q(y, .q_lo); U_def <- .q(y, .q_hi)
    if (!is.finite(L_def)) L_def <- min(y, na.rm = TRUE)
    if (!is.finite(U_def)) U_def <- max(y, na.rm = TRUE)
    span_def  <- U_def - L_def
    full_rng  <- range(y, na.rm = TRUE)
    full_span <- if (!all(is.finite(full_rng))) NA_real_ else diff(full_rng)
    span_floor <- .ds_span_abs_floor
    if (is.finite(full_span) && full_span > 0)
      span_floor <- max(span_floor, .ds_span_rel_floor * full_span)
    span_use <- if (!is.finite(span_def) || span_def <= 0) span_floor else max(span_def, span_floor)
    mid_def <- 0.5 * (L_def + U_def)
    L_def <- mid_def - 0.5 * span_use
    U_def <- mid_def + 0.5 * span_use

    if (g %in% c("max", "min")) {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      # mirror svem_score_random(): degenerate explicit anchors would silently
      # neutralize this response (constant desirability)
      if ((!is.null(gi$lower_acceptable) || !is.null(gi$upper_acceptable)) &&
          (!is.finite(L) || !is.finite(U) || U <= L)) {
        stop("Desirability anchors for response '", r, "' are degenerate ",
             "(need finite lower_acceptable < upper_acceptable; got L = ", L,
             ", U = ", U, ").", call. = FALSE)
      }
      ds_params[[r]] <- list(type = g, L = L, U = U, s = s)
    } else {
      T0 <- as.numeric(gi$target)
      tol       <- if (!is.null(gi$tol)) as.numeric(gi$tol) else NA_real_
      tol_left  <- if (!is.null(gi$tol_left))  as.numeric(gi$tol_left)  else tol
      tol_right <- if (!is.null(gi$tol_right)) as.numeric(gi$tol_right) else tol
      if (!is.finite(tol_left)  || tol_left  <= 0) tol_left  <- .tol_frac * (U_def - L_def)
      if (!is.finite(tol_right) || tol_right <= 0) tol_right <- .tol_frac * (U_def - L_def)
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else (T0 - tol_left)
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else (T0 + tol_right)
      sL <- if (!is.null(gi$shape_left))  as.numeric(gi$shape_left)  else 1
      sR <- if (!is.null(gi$shape_right)) as.numeric(gi$shape_right) else 1
      if (!(is.finite(L) && is.finite(U) && L < T0 && T0 < U)) {
        stop("Resolved target desirability anchors for response '", r,
             "' must satisfy finite L < target < U.", call. = FALSE)
      }
      ds_params[[r]] <- list(type = "target", T0 = T0, L = L, U = U, sL = sL, sR = sR)
    }
  }

  .des <- function(y, p) {
    switch(p$type,
           max    = .ds_max(y, p$L, p$U, p$s),
           min    = .ds_min(y, p$L, p$U, p$s),
           target = .ds_target(y, p$T0, p$L, p$U, p$sL, p$sR))
  }

  # ---- Thompson draws: one member per response per slot, argmax per slot ----
  members <- matrix(NA_integer_, nrow = batch_size, ncol = length(resp_names),
                    dimnames = list(NULL, resp_names))
  picked   <- integer(0)
  draw_val <- matrix(NA_real_, nrow = batch_size, ncol = length(resp_names),
                     dimnames = list(NULL, resp_names))
  draw_des <- matrix(NA_real_, nrow = batch_size, ncol = length(resp_names),
                     dimnames = list(NULL, resp_names))
  ts_score <- numeric(batch_size)

  for (b in seq_len(batch_size)) {
    if (combine == "geom") {
      logscore <- rep(0, n_cand)
    } else {
      score <- rep(0, n_cand)
    }
    d_by_resp <- vector("list", length(resp_names))
    names(d_by_resp) <- resp_names
    for (r in resp_names) {
      m <- sample.int(nrow(objects[[r]]$coef_matrix), 1L)
      members[b, r] <- m
      yv <- eta_list[[r]][, m]
      dv <- .des(yv, ds_params[[r]])
      d_by_resp[[r]] <- list(y = yv, d = dv)
      if (combine == "geom") {
        logscore <- logscore + wts[[r]] * log((1 - geom_floor) * dv + geom_floor)
      } else {
        score <- score + wts[[r]] * dv
      }
    }
    s <- if (combine == "geom") exp(logscore) else score
    s[bad_any] <- -Inf
    s[picked]  <- -Inf
    s[!is.finite(s)] <- -Inf  # NA desirabilities are never eligible
    idx <- which.max(s)
    if (!length(idx) || !is.finite(s[idx]))
      stop("No eligible candidate remains for batch slot ", b, ".")
    picked <- c(picked, idx)
    ts_score[b] <- s[idx]
    for (r in resp_names) {
      draw_val[b, r] <- d_by_resp[[r]]$y[idx]
      draw_des[b, r] <- d_by_resp[[r]]$d[idx]
    }
  }

  # ---- assemble proposals ----
  proposals <- cand_df[picked, , drop = FALSE]
  rownames(proposals) <- NULL
  proposals$slot     <- seq_len(batch_size)
  proposals$ts_score <- ts_score
  for (r in resp_names) {
    proposals[[paste0(r, "_draw")]] <- draw_val[, r]
    proposals[[paste0(r, "_des")]]  <- draw_des[, r]
    proposals[[paste0(r, "_pred")]] <- mean_pred[[r]][picked]
  }

  if (isTRUE(verbose)) {
    cat("svem_thompson_batch: proposed ", batch_size, " run(s) from ",
        n_cand, " candidate(s) for responses: ",
        paste(resp_names, collapse = ", "), "\n", sep = "")
    print(utils::head(proposals[, c("slot", "ts_score",
                                    paste0(resp_names, "_draw"))],
                      batch_size))
  }

  structure(
    list(
      proposals       = proposals,
      members         = members,
      candidate_index = picked,
      candidates      = cand_df,
      weights         = wts,
      ds_params       = ds_params,
      call            = mc
    ),
    class = "svem_thompson_batch"
  )
}

#' Print method for svem_thompson_batch objects
#'
#' @param x An object of class \code{"svem_thompson_batch"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
#' @method print svem_thompson_batch
print.svem_thompson_batch <- function(x, ...) {
  cat("Thompson-sampling batch proposal (svem_thompson_batch)\n")
  cat("  responses:", paste(names(x$weights), collapse = ", "), "\n")
  cat("  batch size:", nrow(x$proposals),
      " | candidates:", nrow(x$candidates), "\n")
  print(x$proposals)
  invisible(x)
}
