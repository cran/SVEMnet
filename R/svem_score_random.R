#' Random-search scoring for SVEM models
#'
#' @description
#' Draw random points from the SVEM sampling schema, compute multi-response
#' desirability scores and (optionally) whole-model-test (WMT) reweighted
#' scores, and attach a scalar uncertainty measure based on percentile CI
#' widths. This function does *not* choose candidates; see
#' \code{\link{svem_select_from_score_table}} for selection and clustering.
#'
#' When \code{specs} is supplied, the function also attempts to append
#' mean-level "in spec" probabilities and related joint indicators using the
#' SVEM bootstrap ensemble via \code{svem_append_design_space_cols}.
#' These quantities reflect uncertainty on the *process mean* at each sampled
#' setting under the fitted SVEM models, not unit-level predictive
#' probabilities. If any error occurs in this spec-limit augmentation, it is
#' caught; a message may be issued when \code{verbose = TRUE}, and the
#' affected table(s) are returned without the spec-related columns.
#'
#' @details
#'
#' \subsection{Typical workflow}{
#' A common pattern is:
#' \enumerate{
#'   \item Fit one or more \code{SVEMnet()} models for the responses of interest.
#'   \item Call \code{svem_score_random()} to:
#'     \itemize{
#'       \item draw candidate settings in factor space,
#'       \item compute Derringer–Suich (DS) desirabilities and a combined
#'             multi-response score, and
#'       \item attach a scalar uncertainty measure derived from percentile CI
#'             widths.
#'     }
#'   \item Optionally provide \code{specs} to append mean-level "in spec"
#'         probabilities and joint indicators based on the SVEM bootstrap
#'         ensemble (process-mean assurance).
#'   \item Use \code{\link{svem_select_from_score_table}} to:
#'     \itemize{
#'       \item select one "best" row (e.g., maximizing \code{score} or
#'             \code{wmt_score}), and
#'       \item pick a small, diverse set of medoid candidates for optimality or
#'             exploration (e.g. high \code{uncertainty_measure}).
#'     }
#'   \item Run selected candidates, append the new data, refit the SVEM models,
#'         and repeat as needed.
#' }
#' }
#'
#' \subsection{Multi-response desirability scoring}{
#' Each response is mapped to a Derringer–Suich desirability
#' \eqn{d_r \in [0,1]} according to its goal:
#' \itemize{
#'   \item \code{goal = "max"}: larger values are better;
#'   \item \code{goal = "min"}: smaller values are better;
#'   \item \code{goal = "target"}: values near a target are best.
#' }
#' Per-response anchors (acceptable lower/upper limits or target-band
#' tolerances) can be supplied in \code{goals}; when not provided, robust
#' defaults are inferred from the sampled responses using the q0.02–q0.98
#' span.
#'
#' Per-response desirabilities are combined into a single scalar \code{score}
#' using either:
#' \itemize{
#'   \item a weighted arithmetic mean (\code{combine = "mean"}), or
#'   \item a weighted geometric mean (\code{combine = "geom"}), with a small
#'         floor applied inside the log to avoid \code{log(0)}.
#' }
#' User-provided weights in \code{goals[[resp]]$weight} are normalized to sum
#' to one and always define \code{weights_original} and the user-weighted
#' \code{score}.
#' }
#'
#' \subsection{Whole-model reweighting (WMT)}{
#' When a WMT object from \code{\link{svem_wmt_multi}} is supplied via the
#' \code{wmt} argument, each response receives a multiplier derived from its
#' whole-model p-value. Final WMT weights are proportional to the product of
#' the user weight and the multiplier, then renormalized to sum to one:
#' \deqn{w_r^{(\mathrm{final})} \propto w_r^{(\mathrm{user})} \times m_r,}
#' where \eqn{m_r} comes from \code{wmt$multipliers}. The user weights always
#' define \code{score}; the WMT-adjusted weights define \code{wmt_score} and
#' the uncertainty weighting when WMT is enabled.
#'
#' \strong{Binomial responses.} If any responses are fitted with
#' \code{family = "binomial"}, supplying a non-\code{NULL} \code{wmt} object
#' is not allowed and the function stops with a clear error. Predictions and
#' CI bounds for binomial responses are interpreted on the probability
#' (response) scale and clamped to \code{[0, 1]} before desirability and
#' uncertainty calculations.
#' }
#'
#' \subsection{Uncertainty measure}{
#' The \code{uncertainty_measure} is a weighted sum of robustly normalized
#' percentile CI widths across responses. For each response, we compute the
#' bootstrap percentile CI width
#' \eqn{\mathrm{CIwidth}_r(x) = u_r(x) - \ell_r(x)} and then map it to the
#' unit interval using an affine rescaling based on the empirical q0.02 and
#' q0.98 quantiles of the CI widths for that response (computed from the
#' table being scored):
#' \deqn{
#'   \tilde W_r(x) =
#'   \frac{
#'     \min\{\max(\mathrm{CIwidth}_r(x), q_{0.02}(r)), q_{0.98}(r)\}
#'     - q_{0.02}(r)
#'   }{
#'     q_{0.98}(r) - q_{0.02}(r)
#'   }.
#' }
#' The scalar \code{uncertainty_measure} is then
#' \deqn{
#'   \text{uncertainty}(x) = \sum_r w_r \, \tilde W_r(x),
#' }
#' where \eqn{w_r} are the user-normalized response weights derived from
#' \code{goals[[resp]]$weight}, regardless of whether a WMT object is
#' supplied. Larger values of \code{uncertainty_measure} indicate settings
#' where the ensemble CI is relatively wide compared to the response's typical
#' scale and are natural targets for exploration.
#' }
#'
#' \subsection{Spec-limit mean-level probabilities}{
#' If \code{specs} is provided, \code{svem_score_random()} attempts to pass
#' the scored table and models to
#' \code{svem_append_design_space_cols} to compute, for each response
#' with an active spec:
#' \itemize{
#'   \item \code{<resp>_p_in_spec_mean}: estimated probability (under the
#'         SVEM bootstrap ensemble) that the process mean at a setting lies
#'         within the specified interval;
#'   \item \code{<resp>_in_spec_point}: 0/1 indicator that the point
#'         prediction lies within the same interval.
#' }
#' and joint quantities:
#' \itemize{
#'   \item \code{p_joint_mean}: product of per-response mean-level
#'         probabilities over responses with active specs;
#'   \item \code{joint_in_spec_point}: 0/1 indicator that all point
#'         predictions are in spec across responses with active specs.
#' }
#' Names in \code{specs} may refer either to \code{names(objects)} or to the
#' model response names; they are automatically aligned to the fitted models.
#'
#' These probabilities are defined on the *conditional means* at each sampled
#' setting, not on individual units or lots, and are best interpreted as
#' ensemble-based assurance measures under the SVEM + FRW pipeline. If the
#' augmentation step fails for any reason (for example, missing predictor
#' columns or incompatible models), the error is caught; a message may be
#' issued when \code{verbose = TRUE}, and \code{score_table} and/or
#' \code{original_data_scored} are returned without the spec-related columns.
#' }
#'
#' @param objects List of \code{svem_model} objects (from
#'   \code{\link{SVEMnet}}). When unnamed, \code{svem_score_random()} attempts
#'   to infer response names from the left-hand sides of the model formulas.
#'   Names (when present) are treated as response identifiers and should
#'   typically match the model response names. All models must share a common
#'   sampling schema (predictor set, factor levels, numeric ranges)
#'   compatible with \code{\link{svem_random_table_multi}}.
#' @param goals List of per-response goal specifications. Either:
#'   \itemize{
#'     \item a named list, where names may be either \code{names(objects)} or
#'           the left-hand-side response names from the fitted models; or
#'     \item an unnamed list with the same length as \code{objects}, in which
#'           case entries are matched to models by position.
#'   }
#'   Each \code{goals[[response]]} must be a list with at least:
#'   \itemize{
#'     \item \code{goal}: one of \code{"max"}, \code{"min"}, \code{"target"};
#'     \item \code{weight}: nonnegative numeric weight.
#'   }
#'   For \code{goal = "target"}, also provide \code{target}. Optional
#'   Derringer–Suich controls:
#'   \itemize{
#'     \item For \code{"max"} or \code{"min"}: \code{lower_acceptable},
#'           \code{upper_acceptable}, \code{shape}.
#'     \item For \code{"target"}: \code{tol} (symmetric), or
#'           \code{tol_left} / \code{tol_right}, and
#'           \code{shape_left} / \code{shape_right}.
#'   }
#'   When anchors/tolerances are not supplied, robust defaults are inferred
#'   from the sampled table using the q0.02–q0.98 span.
#' @param data Optional data frame. When supplied (regardless of whether
#'   \code{wmt} is used), it is scored and returned as
#'   \code{original_data_scored}, with predictions (in \code{<resp>_pred}
#'   columns), per-response desirabilities, \code{score} (and
#'   \code{wmt_score} if applicable), and \code{uncertainty_measure}
#'   appended. When \code{specs} is supplied and the spec-limit augmentation
#'   succeeds, the same mean-level spec columns as in \code{score_table}
#'   (per-response \code{<resp>_p_in_spec_mean}, \code{<resp>_in_spec_point},
#'   and joint \code{p_joint_mean}, \code{joint_in_spec_point}) are appended
#'   as well.
#' @param n Number of random samples to draw in the predictor space. This is
#'   the number of rows in the sampled table used for scoring.
#' @param mixture_groups Optional mixture and simplex constraints passed to
#'   \code{\link{svem_random_table_multi}}. Each group typically specifies
#'   mixture variable names, bounds, and a total.
#' @param level Confidence level for percentile intervals used in the CI width
#'   and uncertainty calculations. Default \code{0.95}.
#' @param combine How to combine per-response desirabilities into a scalar
#'   score. One of:
#'   \itemize{
#'     \item \code{"geom"}: weighted geometric mean (default);
#'     \item \code{"mean"}: weighted arithmetic mean.
#'   }
#' @param numeric_sampler Character string controlling how numeric predictors
#'   are sampled inside \code{\link{svem_random_table_multi}}. One of:
#'   \itemize{
#'     \item \code{"random"}: Latin hypercube sampling when
#'           \pkg{lhs} is available, otherwise independent uniforms;
#'     \item \code{"uniform"}: independent uniforms over stored numeric ranges.
#'   }
#' @param wmt Optional object returned by \code{\link{svem_wmt_multi}}.
#'   When non-\code{NULL}, its \code{multipliers} (and \code{p_values}, if
#'   present) are aligned to \code{names(objects)} and used to define WMT
#'   weights, \code{wmt_score}, and the weighting of the uncertainty measure.
#'   When \code{NULL}, only user weights are used and no WMT reweighting is
#'   applied.
#' @param verbose Logical; if \code{TRUE}, print a compact summary of the run
#'   (and any WMT diagnostics from upstream) to the console.
#' @param specs Optional named list of specification objects, one per response
#'   in \code{objects} for which you want to define a mean-level spec
#'   constraint. Each entry should be either \code{NULL} (no specs for that
#'   response) or a list with components:
#'   \itemize{
#'     \item \code{lower}: numeric lower limit (may be \code{-Inf}, \code{NA},
#'           or \code{NULL} for a one-sided upper spec);
#'     \item \code{upper}: numeric upper limit (may be \code{Inf}, \code{NA},
#'           or \code{NULL} for a one-sided lower spec).
#'   }
#'   Names of \code{specs}, when provided, should be a subset of
#'   \code{names(objects)} or of the model response names (left-hand sides).
#'   The specification structure matches that used by
#'   \code{svem_append_design_space_cols}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{score_table}}{Data frame with predictors, predicted responses
#'     (columns \code{<resp>_pred} for each \code{resp} in \code{names(objects)}),
#'     per-response desirabilities, \code{score}, optional \code{wmt_score},
#'     and \code{uncertainty_measure}. For each response \code{r} in
#'     \code{names(objects)}, additional columns \code{r_lwr}, \code{r_upr}
#'     (percentile CI bounds at level \code{level}) and \code{r_ciw_w}
#'     (weighted, normalized CI width contribution to
#'     \code{uncertainty_measure}) are appended. When \code{specs} is supplied
#'     and the spec-limit augmentation succeeds, additional columns
#'     \code{<resp>_p_in_spec_mean}, \code{<resp>_in_spec_point},
#'     \code{p_joint_mean}, and \code{joint_in_spec_point} are appended.}
#'   \item{\code{original_data_scored}}{If \code{data} is supplied, that data
#'     augmented with prediction columns \code{<resp>_pred}, per-response
#'     desirabilities, \code{score}, optional \code{wmt_score}, and
#'     \code{uncertainty_measure}; otherwise \code{NULL}. When \code{specs} is
#'     supplied and the spec-limit augmentation succeeds, the same mean-level
#'     spec columns as in \code{score_table} are appended to
#'     \code{original_data_scored} as well.}
#'   \item{\code{weights_original}}{User-normalized response weights.}
#'   \item{\code{weights_final}}{Final weights after WMT, if \code{wmt} is
#'     supplied; otherwise equal to \code{weights_original}.}
#'   \item{\code{wmt_p_values}}{Named vector of per-response whole-model
#'     p-values when \code{wmt} is supplied and contains \code{p_values};
#'     otherwise \code{NULL}.}
#'   \item{\code{wmt_multipliers}}{Named vector of per-response WMT multipliers
#'     when \code{wmt} is supplied; otherwise \code{NULL}.}
#' }
#'
#' @seealso
#' \code{\link{SVEMnet}},
#' \code{\link{svem_random_table_multi}},
#' \code{\link{svem_select_from_score_table}},
#' \code{svem_append_design_space_cols()},
#' \code{\link{svem_wmt_multi}}
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------
#' ## Multi-response SVEM scoring with Derringer–Suich desirabilities
#' ## ------------------------------------------------------------------------
#'
#' data(lipid_screen)
#'
#' # Build a deterministic expansion once and reuse for all responses
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_2way  = TRUE,
#'   include_pc_3way  = FALSE
#' )
#'
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' set.seed(1)
#' fit_pot <- SVEMnet(form_pot, lipid_screen)
#' fit_siz <- SVEMnet(form_siz, lipid_screen)
#' fit_pdi <- SVEMnet(form_pdi, lipid_screen)
#'
#' # Collect SVEM models in a named list by response
#' objs <- list(Potency = fit_pot, Size = fit_siz, PDI = fit_pdi)
#'
#' # Targets and user weights for Derringer–Suich desirabilities
#' goals <- list(
#'   Potency = list(goal = "max", weight = 0.6),
#'   Size    = list(goal = "min", weight = 0.3),
#'   PDI     = list(goal = "min", weight = 0.1)
#' )
#'
#' # Optional mixture constraints (composition columns sum to 1)
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' # Basic random-search scoring without WMT or design-space specs
#' set.seed(3)
#' scored_basic <- svem_score_random(
#'   objects         = objs,
#'   goals           = goals,
#'   n               = 10000,          # number of random candidates
#'   mixture_groups  = mix,
#'   combine         = "geom",
#'   numeric_sampler = "random",
#'   verbose         = FALSE
#' )
#'
#' # Scored candidate table: predictors, <resp>_pred, <resp>_des, score, uncertainty
#' names(scored_basic$score_table)
#' head(scored_basic$score_table)
#'
#' # Scored original data (if 'data' is supplied)
#' # scored_basic$original_data_scored contains predictions + desirabilities
#'
#' ## ------------------------------------------------------------------------
#' ## With whole-model tests (WMT) and process-mean specifications
#' ## ------------------------------------------------------------------------
#'
#' set.seed(123)
#' wmt_out <- svem_wmt_multi(
#'   formulas       = list(Potency = form_pot,
#'                         Size    = form_siz,
#'                         PDI     = form_pdi),
#'   data           = lipid_screen,
#'   mixture_groups = mix,
#'   wmt_control    = list(seed = 123),
#'   plot           = FALSE,
#'   verbose        = FALSE
#' )
#'
#' # Simple process-mean specs for a joint design space:
#' #   Potency >= 78, Size <= 100, PDI <= 0.25
#' specs_ds <- list(
#'   Potency = list(lower = 78),
#'   Size    = list(upper = 100),
#'   PDI     = list(upper = 0.25)
#' )
#'
#' set.seed(4)
#' scored_full <- svem_score_random(
#'   objects         = objs,
#'   goals           = goals,
#'   data            = lipid_screen,  # score the original runs as well
#'   n               = 25000,
#'   mixture_groups  = mix,
#'   level           = 0.95,
#'   combine         = "geom",
#'   numeric_sampler = "random",
#'   wmt             = wmt_out,       # optional: WMT reweighting
#'   specs           = specs_ds,      # optional: design-space columns
#'   verbose         = TRUE
#' )
#'
#' # The scored table now includes:
#' #  * score, wmt_score, uncertainty_measure
#' #  * per-response CIs: <resp>_lwr, <resp>_upr
#' #  * design-space columns, e.g. Potency_p_in_spec_mean, p_joint_mean
#' names(scored_full$score_table)
#'
#' ## ------------------------------------------------------------------------
#' ## Positional (unnamed) goals matched to objects by position
#' ## ------------------------------------------------------------------------
#'
#' data(lipid_screen)
#'
#' # Build a deterministic expansion once and reuse for all responses
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_2way  = TRUE,
#'   include_pc_3way  = FALSE
#' )
#'
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' set.seed(1)
#' fit_pot <- SVEMnet(form_pot, lipid_screen)
#' fit_siz <- SVEMnet(form_siz, lipid_screen)
#' fit_pdi <- SVEMnet(form_pdi, lipid_screen)
#'
#' # Collect SVEM models in a list.
#' # Here goals will be matched by position: Potency, Size, PDI.
#' objs <- list(fit_pot, fit_siz, fit_pdi)
#'
#' # Positional goals (unnamed list): must have same length as 'objects'
#' goals_positional <- list(
#'   list(goal = "max", weight = 0.6),  # for Potency (objs[[1]])
#'   list(goal = "min", weight = 0.3),  # for Size    (objs[[2]])
#'   list(goal = "min", weight = 0.1)   # for PDI     (objs[[3]])
#' )
#'
#' set.seed(5)
#' scored_pos <- svem_score_random(
#'   objects         = objs,
#'   goals           = goals_positional,
#'   n               = 5000,
#'   numeric_sampler = "random",
#'   verbose         = FALSE
#' )
#'
#' names(scored_pos$score_table)
#'
#' }

#' @export

svem_score_random <- function(objects,
                              goals,
                              data = NULL,
                              n = 50000,
                              mixture_groups = NULL,
                              level = 0.95,
                              combine = c("geom", "mean"),
                              numeric_sampler = c("random", "uniform"),
                              wmt = NULL,
                              verbose = TRUE,
                              specs   = NULL) {

  # ---- constants / defaults ----
  debias_flag         <- FALSE
  geom_floor          <- 1e-6
  .q_lo               <- 0.02
  .q_hi               <- 0.98
  .tol_frac           <- 0.10
  .exp_left           <- 1.0
  .exp_right          <- 1.0
  .ds_span_abs_floor  <- 1e-6
  .ds_span_rel_floor  <- 0.05


  combine         <- match.arg(combine)
  numeric_sampler <- match.arg(numeric_sampler)

  # ---- basic validation ----
  if (!is.list(objects) || !length(objects))
    stop("objects must be a nonempty list of svem_model objects.")
  if (!all(vapply(objects, inherits, logical(1), what = "svem_model")))
    stop("All elements of 'objects' must be svem_model objects.")

  ## --- infer response names from model formulas ---
  lhs_names <- vapply(objects, function(o) {
    if (!is.null(o$formula) && inherits(o$formula, "formula")) {
      as.character(o$formula[[2L]])
    } else {
      NA_character_
    }
  }, character(1))

  obj_names <- names(objects)
  if (is.null(obj_names)) obj_names <- rep("", length(objects))

  ## fill missing/empty object names from LHS if available
  empty_idx <- which(!nzchar(obj_names) & nzchar(lhs_names))
  if (length(empty_idx)) {
    obj_names[empty_idx] <- lhs_names[empty_idx]
  }
  names(objects) <- obj_names

  ## warn if nonempty list names disagree with LHS
  mismatch <- nzchar(obj_names) & nzchar(lhs_names) & obj_names != lhs_names
  if (any(mismatch)) {
    warning(
      "names(objects) differ from model response names for: ",
      paste0(
        obj_names[mismatch], " (lhs: ", lhs_names[mismatch], ")",
        collapse = ", "
      ),
      ". You can refer to these responses in goals/specs by either the ",
      "list name or the response name."
    )
  }

  resp_names <- names(objects)
  if (any(!nzchar(resp_names))) {
    stop(
      "Could not infer response names for all models; please name the ",
      "objects list or ensure each model has a formula with a left-hand side."
    )
  }

  ## --- goals: allow unnamed (positional) or named by object name or LHS ---
  if (!is.list(goals) || !length(goals))
    stop("goals must be a list of per-response goal specifications.")

  if (is.null(names(goals)) || any(!nzchar(names(goals)))) {
    ## unnamed goals -> match by position
    if (length(goals) != length(objects)) {
      stop(
        "When 'goals' is unnamed, it must have the same length as 'objects'.\n",
        "Either name goals by response, or provide one goal per model."
      )
    }
    names(goals) <- resp_names
  } else {
    ## named goals -> map either object-name or LHS name -> object-name
    gnames <- names(goals)

    ## key -> canonical response name mapping
    key_to_resp <- setNames(resp_names, resp_names)
    valid_lhs   <- nzchar(lhs_names)
    key_to_resp[lhs_names[valid_lhs]] <- resp_names[valid_lhs]

    remapped <- gnames
    for (i in seq_along(gnames)) {
      key <- gnames[i]
      if (!key %in% resp_names && key %in% names(key_to_resp)) {
        remapped[i] <- key_to_resp[[key]]
      }
    }
    names(goals) <- remapped
  }

  miss_goals <- setdiff(resp_names, names(goals))
  if (length(miss_goals)) {
    stop(
      "Missing goals for responses: ",
      paste(miss_goals, collapse = ", "),
      ". You can name goals by the objects list names or by the model ",
      "response names."
    )
  }


  # detect families
  resp_family <- vapply(objects, function(o) {
    fam <- tryCatch(o$family, error = function(e) NA_character_)
    if (!is.character(fam) || !length(fam)) NA_character_ else fam[1]
  }, character(1))
  is_binomial_resp <- tolower(resp_family) %in% "binomial"
  names(is_binomial_resp) <- resp_names

  # WMT not allowed with binomial responses
  if (!is.null(wmt) && any(is_binomial_resp, na.rm = TRUE)) {
    offending <- paste(resp_names[which(is_binomial_resp)], collapse = ", ")
    stop("Supplying `wmt` is not supported when any responses are binomial.\n",
         "Binomial responses detected: ", offending)
  }

  # ---- goals & weights ----
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
    if (!g %in% c("max", "min", "target"))
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

  sw <- sum(goal_df$weight)
  weights_user <- if (sw > 0) goal_df$weight / sw else rep(1 / nrow(goal_df), nrow(goal_df))
  names(weights_user) <- goal_df$response

  # ---- optional WMT: extract from svem_wmt_multi result ----
  wmt_p_vals <- setNames(rep(NA_real_, length(resp_names)), resp_names)
  wmt_mult   <- setNames(rep(1.0,      length(resp_names)), resp_names)

  if (!is.null(wmt)) {
    if (!is.list(wmt) || is.null(wmt$multipliers)) {
      stop("`wmt` must be an object returned by svem_wmt_multi() ",
           "or at least a list with a named numeric 'multipliers' component.")
    }
    mult_vec <- wmt$multipliers
    if (is.null(names(mult_vec))) {
      stop("`wmt$multipliers` must be a named numeric vector.")
    }

    common <- intersect(names(mult_vec), resp_names)
    if (!length(common)) {
      stop("No overlap between names(objects) and names(wmt$multipliers).")
    }
    wmt_mult[common] <- as.numeric(mult_vec[common])

    if (!is.null(wmt$p_values)) {
      p_vec <- wmt$p_values
      if (!is.null(names(p_vec))) {
        common_p <- intersect(names(p_vec), resp_names)
        wmt_p_vals[common_p] <- as.numeric(p_vec[common_p])
      }
    }

    if (isTRUE(verbose)) {
      cat("Using precomputed WMT multipliers for responses:\n")
      print(wmt_mult)
    }
  }

  reweight_flag <- !is.null(wmt)

  # final weights
  weights_final_raw <- weights_user * wmt_mult
  sw2 <- sum(weights_final_raw)
  weights_final <- if (sw2 > 0 && is.finite(sw2)) weights_final_raw / sw2 else weights_user

  wts_user  <- stats::setNames(weights_user, goal_df$response)
  wts_final <- stats::setNames(weights_final, goal_df$response)
  wts_rank  <- if (reweight_flag) wts_final else wts_user
  wts_uncert <- wts_user  # always use user weights for uncertainty


  # ---- sample & predict ----
  sampled_raw <- svem_random_table_multi(
    objects         = objects,
    n               = n,
    mixture_groups  = mixture_groups,
    debias          = debias_flag,
    numeric_sampler = numeric_sampler
  )
  if (!is.list(sampled_raw) || is.null(sampled_raw$data) || is.null(sampled_raw$pred)) {
    stop("svem_random_table_multi must return list(data, pred).")
  }

  data_df <- as.data.frame(sampled_raw$data)
  pred_df <- as.data.frame(sampled_raw$pred, check.names = FALSE, optional = FALSE)

  # map response labels -> prediction columns in sampled_raw$pred
  # and enforce <resp>_pred naming in the score_table
  resp_table_cols <- setNames(character(length(resp_names)), resp_names)

  for (resp in resp_names) {
    obj <- objects[[resp]]

    lhs <- NA_character_
    if (!is.null(obj$formula) && inherits(obj$formula, "formula")) {
      lhs <- tryCatch(
        as.character(obj$formula[[2L]]),
        error = function(e) NA_character_
      )
    }

    candidates <- character(0L)
    if (!is.na(lhs) && nzchar(lhs)) {
      candidates <- c(candidates, paste0(lhs, "_pred"))
    }
    candidates <- unique(c(candidates, paste0(resp, "_pred")))

    found <- candidates[candidates %in% colnames(pred_df)]

    if (!length(found)) {
      stop(
        "Could not match response '", resp,
        "' to prediction columns from svem_random_table_multi(). ",
        "Expected to find one of: ",
        paste(candidates, collapse = ", "),
        " in names(sampled_raw$pred)."
      )
    }

    resp_table_cols[resp] <- found[1L]
  }

  # build prediction matrix with standardized <resp>_pred names
  pred_mat <- as.data.frame(
    lapply(resp_names, function(r) {
      as.numeric(pred_df[[resp_table_cols[[r]]]])
    }),
    check.names = FALSE,
    optional    = FALSE
  )
  colnames(pred_mat) <- paste0(resp_names, "_pred")

  score_table <- cbind(data_df, pred_mat)

  resp_cols       <- resp_names
  resp_pred_cols  <- setNames(paste0(resp_cols, "_pred"), resp_cols)
  predictor_cols  <- colnames(data_df)

  # ---- Derringer–Suich desirabilities ----
  .q       <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  .clip01  <- function(z) pmin(pmax(z, 0), 1)

  .ds_max <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    z <- ifelse(y <= L, 0,
                ifelse(y >= U, 1, ((y - L) / (U - L))^s))
    .clip01(z)
  }
  .ds_min <- function(y, L, U, s = 1) {
    if (!is.finite(L) || !is.finite(U) || U <= L) return(rep(0.5, length(y)))
    z <- ifelse(y <= L, 1,
                ifelse(y >= U, 0, ((U - y) / (U - L))^s))
    .clip01(z)
  }
  .ds_target <- function(y, T0, L, U, sL = 1, sR = 1) {
    if (!is.finite(T0) || !is.finite(L) || !is.finite(U) || !(L < T0 && T0 < U))
      return(rep(0, length(y)))
    left  <- ifelse(y <  L, 0,
                    ifelse(y <= T0, ((y - L) / (T0 - L))^sL, NA_real_))
    right <- ifelse(y >  U, 0,
                    ifelse(y >= T0, ((U - y) / (U - T0))^sR, NA_real_))
    z <- ifelse(is.na(left), right,
                ifelse(is.na(right), left, 1))
    .clip01(z)
  }

  ds_params    <- setNames(vector("list", length(resp_cols)), resp_cols)
  contrib_cols <- list()

  for (r in resp_cols) {
    gi <- goals[[r]]
    g  <- goal_df$goal[goal_df$response == r]

    # predictions used for DS are always <resp>_pred
    y  <- score_table[[resp_pred_cols[[r]]]]

    # ---- validate Derringer–Suich shape parameters (if supplied) ----
    if (!is.null(gi$shape)) {
      if (!is.numeric(gi$shape) || length(gi$shape) != 1L || gi$shape <= 0) {
        stop("For response '", r, "', DS 'shape' must be a single positive number.")
      }
    }
    if (!is.null(gi$shape_left)) {
      if (!is.numeric(gi$shape_left) || length(gi$shape_left) != 1L || gi$shape_left <= 0) {
        stop("For response '", r, "', DS 'shape_left' must be a single positive number.")
      }
    }
    if (!is.null(gi$shape_right)) {
      if (!is.numeric(gi$shape_right) || length(gi$shape_right) != 1L || gi$shape_right <= 0) {
        stop("For response '", r, "', DS 'shape_right' must be a single positive number.")
      }
    }

    if (isTRUE(is_binomial_resp[r])) {
      y <- pmin(pmax(as.numeric(y), 0), 1)
    }

    L_def <- .q(y, .q_lo)
    U_def <- .q(y, .q_hi)
    if (!is.finite(L_def)) L_def <- min(y, na.rm = TRUE)
    if (!is.finite(U_def)) U_def <- max(y, na.rm = TRUE)

    span_def <- U_def - L_def
    full_rng <- range(y, na.rm = TRUE)
    if (!all(is.finite(full_rng))) {
      full_span <- NA_real_
    } else {
      full_span <- full_rng[2L] - full_rng[1L]
    }

    span_floor <- .ds_span_abs_floor
    if (is.finite(full_span) && full_span > 0) {
      span_floor <- max(span_floor, .ds_span_rel_floor * full_span)
    }

    if (!is.finite(span_def) || span_def <= 0) {
      span_use <- span_floor
    } else {
      span_use <- max(span_def, span_floor)
    }

    mid_def <- 0.5 * (L_def + U_def)
    L_def <- mid_def - 0.5 * span_use
    U_def <- mid_def + 0.5 * span_use

    if (g == "max") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      z <- .ds_max(y, L, U, s)
      ds_params[[r]] <- list(type = "max", L = L, U = U, s = s)

    } else if (g == "min") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      z <- .ds_min(y, L, U, s)
      ds_params[[r]] <- list(type = "min", L = L, U = U, s = s)

    } else {  # target
      T0 <- as.numeric(goal_df$target[goal_df$response == r])
      tol       <- if (!is.null(gi$tol)) as.numeric(gi$tol) else NA_real_
      tol_left  <- if (!is.null(gi$tol_left))  as.numeric(gi$tol_left)  else tol
      tol_right <- if (!is.null(gi$tol_right)) as.numeric(gi$tol_right) else tol
      if (!is.finite(tol_left) || tol_left <= 0)  tol_left  <- .tol_frac * (U_def - L_def)
      if (!is.finite(tol_right)|| tol_right <= 0) tol_right <- .tol_frac * (U_def - L_def)

      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else (T0 - tol_left)
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else (T0 + tol_right)

      sL <- if (!is.null(gi$shape_left))  as.numeric(gi$shape_left)  else .exp_left
      sR <- if (!is.null(gi$shape_right)) as.numeric(gi$shape_right) else .exp_right
      if (!(is.finite(L) && is.finite(U) && L < T0 && T0 < U)) {
        band <- .tol_frac * (U_def - L_def)
        L <- T0 - band
        U <- T0 + band
      }
      z <- .ds_target(y, T0, L, U, sL, sR)
      ds_params[[r]] <- list(type = "target", T0 = T0, L = L, U = U, sL = sL, sR = sR)
    }
    contrib_cols[[paste0(r, "_des")]] <- z
  }

  for (nm in names(contrib_cols)) {
    score_table[[nm]] <- contrib_cols[[nm]]
  }

  # ---- combine desirabilities into score + wmt_score ----
  if (combine == "mean") {
    score_user <- Reduce(`+`, lapply(resp_cols, function(r) {
      as.numeric(wts_user[r]) * contrib_cols[[paste0(r, "_des")]]
    }))
    score_rank <- score_user
    wmt_score <- NULL
    if (reweight_flag) {
      wmt_score <- Reduce(`+`, lapply(resp_cols, function(r) {
        as.numeric(wts_final[r]) * contrib_cols[[paste0(r, "_des")]]
      }))
      score_rank <- wmt_score
    }
    score_table$score <- score_user
    if (!is.null(wmt_score)) score_table$wmt_score <- wmt_score

  } else {  # geometric
    log_terms_user <- list()
    log_terms_rank <- list()
    for (r in resp_cols) {
      z    <- contrib_cols[[paste0(r, "_des")]]
      zadj <- (1 - geom_floor) * z + geom_floor
      log_terms_user[[r]] <- as.numeric(wts_user[r]) * log(zadj)
      log_terms_rank[[r]] <- as.numeric(wts_rank[r]) * log(zadj)
    }
    score_user <- exp(Reduce(`+`, log_terms_user))
    score_rank <- exp(Reduce(`+`, log_terms_rank))
    score_table$score <- score_user
    if (reweight_flag) score_table$wmt_score <- score_rank
  }

  # ---- uncertainty_measure from CI widths + store per-response CIs ----
  .normalize01_robust <- function(x) {
    a <- .q(x, .q_lo); b <- .q(x, .q_hi)
    if (!is.finite(a) || !is.finite(b) || b <= a) return(rep(0.5, length(x)))
    (pmin(pmax(x, a), b) - a) / (b - a)
  }

  ciw_w_cols  <- list()
  ci_lwr_cols <- list()
  ci_upr_cols <- list()

  for (r in resp_cols) {
    obj <- objects[[r]]
    pr <- try(
      predict(obj,
              newdata  = score_table[, predictor_cols, drop = FALSE],
              debias   = debias_flag,
              interval = TRUE,
              level    = level),
      silent = TRUE
    )

    lwr <- upr <- rep(NA_real_, nrow(score_table))
    if (!inherits(pr, "try-error") && is.list(pr) &&
        !is.null(pr$lwr) && !is.null(pr$upr)) {
      lwr <- as.numeric(pr$lwr)
      upr <- as.numeric(pr$upr)
    }

    if (isTRUE(is_binomial_resp[r])) {
      lwr <- pmin(pmax(lwr, 0), 1)
      upr <- pmin(pmax(upr, 0), 1)
    }

    # store CI columns
    ci_lwr_cols[[paste0(r, "_lwr")]] <- lwr
    ci_upr_cols[[paste0(r, "_upr")]] <- upr

    width <- upr - lwr
    zc    <- .normalize01_robust(width)
    ciw_w_cols[[paste0(r, "_ciw_w")]] <- as.numeric(wts_uncert[r]) * zc
  }

  if (length(ciw_w_cols)) {
    # attach CIs first, then weighted CI width contributors, then scalar measure
    for (nm in names(ci_lwr_cols)) score_table[[nm]] <- ci_lwr_cols[[nm]]
    for (nm in names(ci_upr_cols)) score_table[[nm]] <- ci_upr_cols[[nm]]
    for (nm in names(ciw_w_cols))  score_table[[nm]] <- ciw_w_cols[[nm]]
    score_table$uncertainty_measure <- Reduce(`+`, ciw_w_cols)
  } else {
    score_table$uncertainty_measure <- NA_real_
  }

  # ---- score original data (optional) ----
  original_data_scored <- NULL
  if (!is.null(data) && is.data.frame(data)) {
    pred_df <- data.frame(row_id = seq_len(nrow(data)))
    ci_lwr  <- list()
    ci_upr  <- list()

    for (r in resp_names) {
      pr <- try(
        predict(objects[[r]], newdata = data, debias = debias_flag,
                interval = TRUE, level = level),
        silent = TRUE
      )
      fit_col <- if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$fit)) {
        as.numeric(pr$fit)
      } else {
        as.numeric(try(
          predict(objects[[r]], newdata = data, debias = debias_flag),
          silent = TRUE
        ))
      }
      if (isTRUE(is_binomial_resp[r])) {
        fit_col <- pmin(pmax(fit_col, 0), 1)
      }
      pred_df[[r]] <- fit_col

      if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
        lwr <- as.numeric(pr$lwr)
        upr <- as.numeric(pr$upr)
        if (isTRUE(is_binomial_resp[r])) {
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

    des_cols <- list()
    for (r in resp_names) {
      y <- pred_df[[r]]
      p <- ds_params[[r]]
      if (p$type == "max") {
        des_cols[[paste0(r, "_des")]] <- .ds_max(y, p$L, p$U, p$s)
      } else if (p$type == "min") {
        des_cols[[paste0(r, "_des")]] <- .ds_min(y, p$L, p$U, p$s)
      } else {
        des_cols[[paste0(r, "_des")]] <- .ds_target(y, p$T0, p$L, p$U, p$sL, p$sR)
      }
    }

    # scores
    if (combine == "mean") {
      score_user <- Reduce(`+`, lapply(resp_names, function(r) {
        as.numeric(wts_user[r]) * des_cols[[paste0(r, "_des")]]
      }))
      score_rank     <- score_user
      wmt_score_vec  <- NULL
      if (reweight_flag) {
        wmt_score_vec <- Reduce(`+`, lapply(resp_names, function(r) {
          as.numeric(wts_final[r]) * des_cols[[paste0(r, "_des")]]
        }))
        score_rank <- wmt_score_vec
      }
    } else {
      logsum_user <- Reduce(`+`, lapply(resp_names, function(r) {
        z  <- des_cols[[paste0(r, "_des")]]
        zf <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts_user[r]) * log(zf)
      }))
      score_user <- exp(logsum_user)

      logsum_rank <- Reduce(`+`, lapply(resp_names, function(r) {
        z  <- des_cols[[paste0(r, "_des")]]
        zf <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts_rank[r]) * log(zf)
      }))
      score_rank    <- exp(logsum_rank)
      wmt_score_vec <- if (reweight_flag) score_rank else NULL
    }

    # uncertainty for original data + per-response CI width contributions
    unc_w      <- list()
    ciw_w_orig <- list()

    for (r in resp_names) {
      width <- ci_upr[[r]] - ci_lwr[[r]]

      a <- .q(width, .q_lo)
      b <- .q(width, .q_hi)
      zc <- if (!is.finite(a) || !is.finite(b) || b <= a) {
        rep(0.5, length(width))
      } else {
        (pmin(pmax(width, a), b) - a) / (b - a)
      }

      unc_w[[r]]      <- as.numeric(wts_uncert[r]) * zc
      ciw_w_orig[[r]] <- as.numeric(wts_uncert[r]) * zc


    }
    unc_vec <- Reduce(`+`, unc_w)

    pred_mat <- pred_df[resp_names]
    colnames(pred_mat) <- paste0(resp_names, "_pred")

    original_data_scored <- cbind(
      data,
      pred_mat,
      as.data.frame(des_cols, stringsAsFactors = FALSE)
    )

    # append CI bounds and CI width contributions for each response
    for (r in resp_names) {
      original_data_scored[[paste0(r, "_lwr")]]   <- ci_lwr[[r]]
      original_data_scored[[paste0(r, "_upr")]]   <- ci_upr[[r]]
      original_data_scored[[paste0(r, "_ciw_w")]] <- ciw_w_orig[[r]]
    }

    original_data_scored$score <- score_user
    if (reweight_flag && !is.null(wmt_score_vec)) {
      original_data_scored$wmt_score <- score_rank
    }
    original_data_scored$uncertainty_measure <- unc_vec
    attr(original_data_scored, "svem_predictor_cols") <- colnames(data)
    attr(original_data_scored, "svem_resp_cols")      <- resp_cols
  }


  # ---- optional spec-limit augmentation (safe) ----
  if (!is.null(specs)) {

    ## allow specs to be named by objects list names OR by formula LHS
    if (!is.null(names(specs))) {
      s_names <- names(specs)

      # mapping: key -> canonical response name
      #   * resp_names (names(objects))
      #   * lhs_names from model formulas, where available
      key_to_resp <- setNames(resp_names, resp_names)

      valid_lhs <- nzchar(lhs_names)
      if (any(valid_lhs)) {
        key_to_resp[lhs_names[valid_lhs]] <- resp_names[valid_lhs]
      }

      remapped <- s_names
      for (i in seq_along(s_names)) {
        key <- s_names[i]
        if (!key %in% resp_names && key %in% names(key_to_resp)) {
          remapped[i] <- key_to_resp[[key]]
        }
      }
      names(specs) <- remapped
    }

    # augment score_table
    score_table <- tryCatch(
      {
        svem_append_design_space_cols(
          score_table = score_table,
          objects     = objects,
          specs       = specs
        )
      },
      error = function(e) {
        if (isTRUE(verbose)) {
          message(
            "svem_score_random(): spec-limit augmentation failed for score_table; ",
            "returning score_table without spec columns. Error: ",
            conditionMessage(e)
          )
        }
        score_table
      }
    )

    # augment original_data_scored, if present
    if (!is.null(original_data_scored) && is.data.frame(original_data_scored)) {
      original_data_scored <- tryCatch(
        {
          svem_append_design_space_cols(
            score_table = original_data_scored,
            objects     = objects,
            specs       = specs
          )
        },
        error = function(e) {
          if (isTRUE(verbose)) {
            message(
              "svem_score_random(): spec-limit augmentation failed for original_data_scored; ",
              "returning original_data_scored without spec columns. Error: ",
              conditionMessage(e)
            )
          }
          original_data_scored
        }
      )
    }
  }


  # store basic attributes for downstream helpers
  attr(score_table, "svem_predictor_cols") <- predictor_cols
  attr(score_table, "svem_resp_cols")      <- resp_cols

  if (!is.null(original_data_scored) && is.data.frame(original_data_scored)) {
    if (!is.null(data) && is.data.frame(data)) {
      attr(original_data_scored, "svem_predictor_cols") <- colnames(data)
    }
    attr(original_data_scored, "svem_resp_cols") <- resp_cols
  }

  list(
    score_table          = score_table,
    original_data_scored = original_data_scored,
    weights_original     = wts_user,
    weights_final        = wts_final,
    wmt_p_values         = if (reweight_flag) wmt_p_vals else NULL,
    wmt_multipliers      = if (reweight_flag) wmt_mult   else NULL
  )
}
