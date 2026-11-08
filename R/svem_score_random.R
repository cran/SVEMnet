.svem_generated_score_names <- function(responses) {
  unique(c(unlist(lapply(responses, function(r) paste0(r, c(
    "_pred", "_des", "_lwr", "_upr", "_ciw_w", "_p_in_spec_mean",
    "_in_spec_point"
  ))), use.names = FALSE), "score", "wmt_score", "uncertainty_measure",
  "p_joint_mean", "joint_in_spec_point"))
}

.svem_check_output_collisions <- function(protected, generated, context) {
  clash <- intersect(protected, generated)
  if (length(clash)) {
    stop(context, ": generated output column(s) would overwrite input variable(s): ",
         paste(clash, collapse = ", "),
         ". Rename these inputs or the response identifiers before proceeding.",
         call. = FALSE)
  }
  invisible(NULL)
}

.svem_align_goals <- function(goals, resp_names, key_to_resp) {
  if (!is.list(goals) || !length(goals)) {
    stop("`goals` must be a nonempty list of per-response goal specifications.",
         call. = FALSE)
  }

  goal_names <- names(goals)
  if (is.null(goal_names)) {
    if (length(goals) != length(resp_names)) {
      stop("When `goals` is unnamed, it must have the same length as `objects`.",
           call. = FALSE)
    }
    names(goals) <- resp_names
    return(goals)
  }

  blank <- is.na(goal_names) | !nzchar(goal_names)
  if (all(blank)) {
    if (length(goals) != length(resp_names)) {
      stop("When `goals` is unnamed, it must have the same length as `objects`.",
           call. = FALSE)
    }
    names(goals) <- resp_names
    return(goals)
  }
  if (any(blank)) {
    stop("`goals` must be either fully named or fully unnamed; blank/NA names cannot be mixed with response names.",
         call. = FALSE)
  }

  unknown <- setdiff(goal_names, names(key_to_resp))
  if (length(unknown)) {
    stop("Unknown goal response name(s): ", paste(unknown, collapse = ", "),
         ". Use an `objects` list name or a fitted-model response name.",
         call. = FALSE)
  }

  remapped <- unname(key_to_resp[goal_names])
  if (anyDuplicated(remapped)) {
    duplicate_responses <- unique(remapped[duplicated(remapped)])
    stop("Goal names/aliases map to duplicate responses: ",
         paste(duplicate_responses, collapse = ", "),
         ". Supply exactly one goal specification per response.",
         call. = FALSE)
  }

  missing_responses <- setdiff(resp_names, remapped)
  if (length(missing_responses)) {
    stop("Missing goals for responses: ", paste(missing_responses, collapse = ", "),
         ". You can name goals by the objects list names or by the model response names.",
         call. = FALSE)
  }

  goals <- goals[match(resp_names, remapped)]
  names(goals) <- resp_names
  goals
}

.svem_response_key_map <- function(resp_names, lhs_names, valid_lhs) {
  if (length(resp_names) != length(lhs_names) ||
      length(valid_lhs) != length(lhs_names)) {
    stop("Internal response-name alignment failure.", call. = FALSE)
  }

  # A fitted-model LHS must not also be another model's canonical list name:
  # in that case a goal key would have two plausible meanings.
  conflict <- valid_lhs & lhs_names %in% resp_names & lhs_names != resp_names
  if (any(conflict)) {
    pairs <- paste0(
      "'", lhs_names[conflict], "' (LHS of '", resp_names[conflict], "')"
    )
    stop(
      "Ambiguous response alias(es): ", paste(pairs, collapse = ", "),
      " also match another object-list name. Rename the conflicting list entry.",
      call. = FALSE
    )
  }

  key_to_resp <- stats::setNames(resp_names, resp_names)
  lhs_alias <- valid_lhs & !(lhs_names %in% resp_names)
  if (any(lhs_alias)) key_to_resp[lhs_names[lhs_alias]] <- resp_names[lhs_alias]
  key_to_resp
}

.svem_validate_goals <- function(goals, resp_names) {
  goal_type <- stats::setNames(character(length(resp_names)), resp_names)
  weights   <- stats::setNames(numeric(length(resp_names)), resp_names)
  targets   <- stats::setNames(rep(NA_real_, length(resp_names)), resp_names)

  for (r in resp_names) {
    gi <- goals[[r]]
    if (!is.list(gi) || is.null(gi$goal) || is.null(gi$weight)) {
      stop("Each goals[[response]] must have goal and weight. Offender: ", r,
           call. = FALSE)
    }

    if (!is.character(gi$goal) || length(gi$goal) != 1L ||
        is.na(gi$goal) || !nzchar(gi$goal)) {
      stop("goal for ", r, " must be a single character value: 'max', 'min', or 'target'.",
           call. = FALSE)
    }
    g <- tolower(gi$goal)
    if (!(g %in% c("max", "min", "target"))) {
      stop("goal for ", r, " must be 'max', 'min', or 'target'.", call. = FALSE)
    }

    gi$weight <- .svem_numeric_scalar(
      gi$weight, paste0("goals[['", r, "']]$weight"), lower = 0
    )

    for (field in c("lower_acceptable", "upper_acceptable")) {
      if (!is.null(gi[[field]])) {
        gi[[field]] <- .svem_numeric_scalar(
          gi[[field]], paste0("goals[['", r, "']]$", field)
        )
      }
    }
    for (field in c("tol", "tol_left", "tol_right", "shape",
                    "shape_left", "shape_right")) {
      if (!is.null(gi[[field]])) {
        gi[[field]] <- .svem_numeric_scalar(
          gi[[field]], paste0("goals[['", r, "']]$", field),
          lower = 0, lower_open = TRUE
        )
      }
    }

    if (!is.null(gi$lower_acceptable) && !is.null(gi$upper_acceptable) &&
        gi$lower_acceptable >= gi$upper_acceptable) {
      stop("Desirability anchors for response '", r,
           "' are degenerate; they must satisfy lower_acceptable < upper_acceptable.",
           call. = FALSE)
    }

    if (!is.null(gi$target)) {
      gi$target <- .svem_numeric_scalar(
        gi$target, paste0("goals[['", r, "']]$target")
      )
    }

    if (g == "target") {
      if (is.null(gi$target)) {
        stop("target must be provided for ", r, " when goal = 'target'.",
             call. = FALSE)
      }
      if (!is.null(gi$lower_acceptable) && gi$lower_acceptable >= gi$target) {
        stop("Target desirability for response '", r,
             "' requires lower_acceptable < target.", call. = FALSE)
      }
      if (!is.null(gi$upper_acceptable) && gi$upper_acceptable <= gi$target) {
        stop("Target desirability for response '", r,
             "' requires target < upper_acceptable.", call. = FALSE)
      }
      targets[[r]] <- gi$target
    }

    gi$goal <- g
    goals[[r]] <- gi
    goal_type[[r]] <- g
    weights[[r]] <- gi$weight
  }

  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("At least one response weight must be strictly positive, and response weights must have a finite sum.",
         call. = FALSE)
  }

  list(
    goals = goals,
    goal = goal_type,
    weight = weights,
    target = targets,
    normalized_weight = weights / weight_sum
  )
}

#' Random-search scoring for SVEM models
#'
#' @description
#' Draw random points from the SVEM sampling schema, compute multi-response
#' desirability scores and (optionally) whole-model-test (WMT) reweighted
#' scores, and attach a scalar uncertainty measure based on percentile CI
#' widths. Companion fits or an explicit preparation budget also return a
#' worklist of recipes to test. For other selections from the scored table,
#' see \code{\link{svem_select_from_score_table}}.
#'
#' Predictions used inside this scorer are always generated with debiasing
#' disabled (i.e., \code{debias = FALSE}) regardless of whether the underlying
#' SVEM fits support calibration.
#'
#' When \code{specs} is supplied, the function also attempts to append
#' mean-level "in spec" probabilities and related joint indicators using the
#' SVEM bootstrap ensemble via \code{svem_append_design_space_cols}.
#' These quantities reflect uncertainty on the *process mean* at each sampled
#' setting under the fitted SVEM models, not unit-level predictive
#' probabilities. If any error occurs in this spec-limit augmentation, it is
#' caught with a warning regardless of \code{verbose}, and the
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
#'     The robust bounded summaries \code{score}, per-response desirabilities,
#'     per-response weighted CI-width summaries, and
#'     \code{uncertainty_measure} may contain ties at 0 or 1 because of the
#'     percentile clipping described below. If the intended single "best" row is
#'     the largest raw predicted response, select on the corresponding
#'     \code{<response>_pred} column. If the intended single exploration row is
#'     the widest interval, create and select on a raw width column such as
#'     \code{<response>_ciw = <response>_upr - <response>_lwr}.
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
#' span. As a result, the largest few percent of predictions may share
#' desirability 1 (and the smallest few percent may share desirability 0) when
#' default anchors are used. This is intentional for robust desirability
#' scoring. Use \code{<response>_pred} directly when selecting the single
#' largest or smallest predicted response.
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
#' define \code{score}; the WMT-adjusted weights define \code{wmt_score}.
#' The uncertainty measure is always weighted using the user weights, even
#' when WMT is supplied.
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
#' q0.98 quantiles of the CI widths for that response. The quantile anchors
#' are computed once from the sampled \code{score_table} and reused when
#' scoring \code{data}, so \code{uncertainty_measure} in
#' \code{original_data_scored} is on the same scale as in \code{score_table}:
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
#' \code{goals[[resp]]$weight}. Larger values of \code{uncertainty_measure}
#' indicate settings where the ensemble CI is relatively wide compared to the
#' response's typical scale and are natural targets for exploration. Because
#' the normalization clips at empirical q0.02 and q0.98, the widest few percent
#' of intervals may share \code{uncertainty_measure = 1}. Use a raw interval
#' width column such as
#' \code{<response>_ciw = <response>_upr - <response>_lwr} when selecting the
#' single widest interval.
#' If a model has fewer than two bootstrap members or its interval prediction
#' fails, the CI columns are \code{NA}. Its uncertainty contribution is the
#' neutral value \code{0.5} times its user weight, and a warning identifies the
#' affected response. These fallback values do not measure fitted uncertainty.
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
#' columns or incompatible models), the error is caught with a warning
#' regardless of \code{verbose}, and \code{score_table} and/or
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
#'   as well. Recognized columns from an earlier scoring round are replaced
#'   with a warning. A generated output name that collides with a model
#'   predictor or an observed response variable is rejected before sampling.
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
#'   weights, \code{wmt_score}. When \code{NULL}, only
#'   user weights are used and no WMT reweighting is applied.
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
#'     \code{p_joint_mean}, and \code{joint_in_spec_point} are appended.
#'     A supplied candidate pool or a worklist adds \code{recipe_id}.
#'     A worklist also includes \code{<resp>_<engine>_pred} for each
#'     requested companion, \code{<engine>_score} for its joint goal and,
#'     when requested, \code{<resp>_gp_latent_sd} for GP exploration.}
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
#' With companions, \code{budget} or \code{gp_endpoint}, the result also includes:
#' \describe{
#'   \item{\code{worklist}}{Scored recipes in preparation order, with \code{slot},
#'     \code{role}, \code{nominated_by} (all endorsing engines) and
#'     \code{nearest_distance} to a previously prepared or earlier selected recipe.
#'     SVEM spread exploration uses \code{maximum_spread} for a unique maximum and
#'     \code{best_predicted_tie} when predicted desirability breaks a spread tie.}
#'   \item{\code{nominations}}{One row per engine with its unrestricted and best
#'     unprepared recipe IDs and tie counts. Status is \code{selected},
#'     \code{dropped_by_budget}, \code{dropped_by_model_count},
#'     \code{not_prioritized}, \code{not_supplied}, \code{no_finite_score} or
#'     \code{no_unprepared_recipe}.
#'     The \code{slot} column identifies the selected preparation and is missing
#'     for omitted optima. Engines endorsing the same recipe share a slot.}
#'   \item{\code{engine_status}}{Fit or prediction status and a diagnostic
#'     message for each engine and response. A companion needs finite predictions
#'     for every positively weighted response on the same recipe to nominate it.}
#'   \item{\code{exploration}}{The best predicted and widest-spread recipes among
#'     the maximum SVEM spread ties, with tie counts and desirability ranges.}
#'   \item{\code{tie_endpoint_ranges}}{Ranges of endpoint predictions, raw
#'     interval widths and weighted spread contributions among those ties.}
#'   \item{\code{selections}}{Worklist rows formatted for \code{svem_export_candidates_csv()}.}
#'   \item{\code{shortfall}}{Without explicit counts, the number of budgeted
#'     preparations left unfilled. With \code{n_optima} or \code{n_explore},
#'     the total model and exploration shortfall. Only an explicitly requested
#'     model count can contribute a model shortfall. Unused physical capacity
#'     is reported separately in \code{allocation_counts$unused_capacity}.}
#'   \item{\code{historical_unmatched}}{Number of history rows with missing or
#'     off-increment settings that could not be matched.}
#'   \item{\code{history_source}}{Whether previous preparations came from supplied
#'     data, all fitted response histories, some histories, or no history.}
#' }
#' With explicit counts, the result also includes:
#' \describe{
#'   \item{\code{allocation_counts}}{Requested, targeted and selected model
#'     and exploration counts, distinct required repeats, role shortfalls,
#'     unfilled model capacity and unused physical capacity.}
#'   \item{\code{merged_nominations}}{Endorsements sharing a physical
#'     preparation with an earlier role.}
#' }
#'
#' @seealso
#' \code{\link{SVEMnet}},
#' \code{\link{svem_random_table_multi}},
#' \code{\link{svem_select_from_score_table}},
#' \code{svem_append_design_space_cols()},
#' \code{\link{svem_wmt_multi}}
#' @template ref-svem
#' @examples
#' ## Small runnable two-response example
#' set.seed(11)
#' n_toy <- 36
#' toy <- data.frame(x1 = runif(n_toy), x2 = runif(n_toy))
#' toy$yield <- 1 + 2 * toy$x1 - toy$x2 + rnorm(n_toy, sd = 0.15)
#' toy$impurity <- 1.5 - toy$x1 + 0.5 * toy$x2 + rnorm(n_toy, sd = 0.10)
#'
#' fit_yield <- SVEMnet(yield ~ x1 + x2, toy, nBoot = 5,
#'                      glmnet_alpha = 1, relaxed = FALSE)
#' fit_impurity <- SVEMnet(impurity ~ x1 + x2, toy, nBoot = 5,
#'                         glmnet_alpha = 1, relaxed = FALSE)
#'
#' toy_scored <- svem_score_random(
#'   objects = list(yield = fit_yield, impurity = fit_impurity),
#'   goals = list(
#'     yield = list(goal = "max", weight = 0.7),
#'     impurity = list(goal = "min", weight = 0.3)
#'   ),
#'   n = 100,
#'   numeric_sampler = "uniform",
#'   verbose = FALSE
#' )
#' head(toy_scored$score_table)
#'
#' \dontrun{
#' ## Production-scale lipid workflow. This is not run automatically because
#' ## it fits three wide ensembles and optionally performs permutation WMTs.
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
#'   wmt_control    = list(seed = 123, nCore = 1L),
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
#'
#' @export
#' @param candidates Optional feasible recipe data frame to score in place of
#'   random sampling. Include all predictors. A unique `recipe_id` is optional.
#'   Repeated physical recipes collapse to the smallest supplied ID (numeric
#'   order or bytewise order for character IDs). Without IDs the scorer assigns
#'   ranks of sorted unique recipe keys, independent of the input row order.
#'   For a worklist, hold each block at its planned operating level. Physical
#'   recipe keys exclude those nuisance blocks.
#' @param budget Optional positive physical preparation count for a worklist.
#'   Under `allocation = "legacy"` it counts distinct new preparations. The
#'   order is elastic-net SVEM, GauPro, Extra Trees, forward and spline optima
#'   when supplied, then SVEM spread, GP latent SD and high-spread medoids.
#'   Under `allocation = "portfolio"` it counts all physical preparations,
#'   including incumbent and benchmark repeats. Shared optima use one place.
#'   With unspecified exploration count, exploration can fill remaining slots.
#'   Supply `data` to exclude previously prepared recipes across all blocks,
#'   or use the recorded preparation settings in companion fits. A worklist's
#'   `history_source` reports the coverage available for this exclusion.
#'   For the legacy allocation, budget required controls separately.
#' @param gp_endpoint Optional response name for GP exploration by latent SD.
#'   A single Gaussian response is inferred automatically when GauPro was fitted.
#' @param increments Optional named numeric vector of positive preparation
#'   increments for an already discretized supplied pool, measured from zero.
#'   For example, use 0.0001 for component fractions in basis points.
#' @param forward Optional named list of Gaussian forward fits for the same
#'   responses. They use the same desirability anchors and scientific weights.
#' @param spline Optional named list of [svem_spline()] fits, one per response.
#'   Scoring uses the same physical factors, block declarations, mixture
#'   constraints, scientific weights and desirability anchors as the primary
#'   fits. Supplying a spline adds its fitted means, joint score and nomination.
#'   A spline treats known block levels as categories and is incompatible with
#'   a primary fit that declares a continuous companion block.
#'   This scorer uses previously fitted objects. Use [svem_plan_runs()] to fit
#'   models and prepare a worklist together.
#' @param n_optima Optional number of distinct new model nominations in
#'   portfolio allocation. Duplicate model optima share one preparation and
#'   the next distinct prioritized model can fill the requested count. With
#'   `n_optima = NULL`, available models fill places left after required
#'   repeats and reserved exploration. An unfilled automatic model place is
#'   recorded as capacity.
#' @param n_explore Optional exact number of uncertainty-medoid preparations.
#'   If both counts are supplied, an omitted budget is their sum plus the
#'   distinct required reference repeats. Supply a total `budget` when only
#'   one count is specified. Reference repeats consume budget first. With `n_explore = NULL`, remaining
#'   eligible slots are filled by exploration. Explicit counts can leave
#'   capacity unused. Inspect `allocation_counts` and `merged_nominations`
#'   for requested and obtained roles, shortages and shared physical recipes.
#' @param allocation `"legacy"` preserves the existing new-recipe allocation.
#'   `"portfolio"` uses the declared priority, then uncertainty medoids. The
#'   default formulation workflow in [svem_plan_runs()] selects this policy.
#' @param incumbent For portfolio allocation, `TRUE` repeats the best observed
#'   preparation in the reference block context; `FALSE` omits it. A positive
#'   integer specifies a row of `data` explicitly. The repeat consumes one slot.
#'   Each historical row must represent an independent preparation; summarize
#'   technical readings first. Selection uses observed outcomes, not fitted
#'   historical means. Missing positively weighted outcomes are ineligible.
#' @param benchmark Optional fixed historical benchmark, carried in addition to
#'   the selected incumbent. Supply one row index of \code{data}, or a named
#'   one-element list identifying a unique row, such as
#'   \code{list(RunID = "original-control")}. The identifier form is convenient
#'   across refits. The recipe must satisfy the current preparation constraints;
#'   it need not have complete response observations. A requested benchmark is
#'   mandatory and consumes a slot. If it is the same physical recipe as the
#'   incumbent, one preparation carries both roles. Too small a budget errors.
#' @param incumbent_summary How repeated independent preparations within a
#'   reference context are ranked. The default \code{"best"} retains the
#'   highest individual observation. A poor repeat cannot demote that stored
#'   maximum. Use \code{"mean"} to score each physical recipe's mean observed
#'   responses over complete preparations; technical readings must be summarized
#'   beforehand. An explicit incumbent row always uses that row alone.
#' @param reference Named list with one historical level of every blocking
#'   factor for automatic incumbent selection. By default uses the planned
#'   levels in `candidates`; raw observations from different blocks are not
#'   pooled. An explicit incumbent row instead supplies its own source context.
#' @param priority Portfolio nomination order. The usual order is incumbent,
#'   elastic net, GauPro, Extra Trees and forward. The default order appends
#'   spline after forward only when a spline fit list is supplied to this scorer.
#'   A custom order can omit or reorder models. Shared optima consume one place.
#' @param exploration_fraction Fraction of remaining eligible recipes offered
#'   to uncertainty medoid selection in portfolio allocation, default 0.1.
#'   The shortlist expands only enough to fill available slots when necessary;
#'   `allocation_details` records that expansion. Use 0.25 for a wider sensitivity.
#' @details When companions or a budget are requested, `score_table` also
#'   contains one set of endpoint means and one comparable desirability per
#'   available engine. The returned `worklist` has one row per distinct
#'   nomination after excluding known preparations. Its
#'   `selections` can be passed to `svem_export_candidates_csv()`.
#'   `nominations` reports each engine's unrestricted and best unprepared
#'   candidate. `exploration` contrasts best predicted desirability among tied
#'   maximum SVEM spread scores with an unclipped standardized-width tie-break.
#'   All engines use the original scientific weights and the same SVEM-resolved
#'   desirability anchors. A sole max/min endpoint ranks by the raw fitted mean.
#'   For several endpoints, equally desirable candidates remain tied even when
#'   their predictions differ outside the specified desirability limits; a
#'   stable recipe key resolves those ties. The scorer does not invent extra
#'   utility beyond the declared goals. Mean-interval, WMT and specification
#'   diagnostics remain specific to SVEM.
#'
#'   Portfolio results also return \code{uncertainty_anchors} and a
#'   \code{uncertainty_eligibility} table with recipe IDs, informative-width
#'   flags and unclipped standardized widths. Constant or unavailable interval
#'   widths cannot support exploration. The worklist records its prediction
#'   context and the incumbent's observed source context and rows. The
#'   incumbent's \code{by_context} table reports observed winners separately,
#'   without pooling or adjusting outcomes across blocks. The incumbent object
#'   (not its \code{by_context} table) also records
#'   \code{recipe_n_preparations} (all attempted preparations of the selected
#'   physical recipe in its source context), \code{recipe_n_complete}, and
#'   \code{recipe_mean_score}, the desirability of its complete-preparation
#'   mean responses. These diagnostics do not change the default best-row
#'   selection. \code{n_preparations} counts the source rows used by the selected
#'   statistic, so it remains one for a best individual observation. The worklist
#'   exposes the recipe counts and mean score with the prefix
#'   \code{incumbent_recipe_}.
#'
#'   The companion worklist is a custom experimental confirmation workflow.
#'   GauPro nominates an optimum by the fitted response mean or joint
#'   desirability. A separate optional slot uses the largest latent standard
#'   deviation for a chosen endpoint. A full Bayesian optimization workflow
#'   specifies an acquisition rule and updates the choice as data arrive.
#'   Expected improvement, for example, combines improvement and uncertainty
#'   in one acquisition value (Jones, Schonlau and Welch, 1998,
#'   \doi{10.1023/A:1008306431147}; Frazier, 2018,
#'   \url{https://arxiv.org/abs/1807.02811}). Here the scientist allocates
#'   preparations, measures the responses and requests the refit.
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
                              specs   = NULL,
                              candidates = NULL,
                              budget = NULL,
                              gp_endpoint = NULL,
                              increments = NULL,
                              forward = NULL,
                              allocation = c("legacy", "portfolio"),
                              incumbent = FALSE,
                              reference = NULL,
                              priority = c("incumbent", "svemnet", "gaupro", "extra_trees", "forward"),
                              exploration_fraction = .1,
                              incumbent_summary = c("best", "mean"),
                              benchmark = NULL,
                              spline = NULL, n_optima = NULL, n_explore = NULL) {

  # ---- constants / defaults ----
  debias_flag         <- FALSE  # always off here by design
  geom_floor          <- 1e-6
  .q_lo               <- 0.02
  .q_hi               <- 0.98
  .tol_frac           <- 0.10
  .exp_left           <- 1.0
  .exp_right          <- 1.0
  .ds_span_abs_floor  <- 1e-6
  .ds_span_rel_floor  <- 0.05

  use_default_priority <- missing(priority)
  allocation <- match.arg(allocation)
  # A scorer cannot fit: only supplied spline fits extend its default priority.
  priority <- .svem_resolve_plan_priority(priority, use_default_priority && is.list(spline),
                                          spline, n_optima, n_explore)
  if ((!is.null(n_optima) || !is.null(n_explore)) && allocation != "portfolio") {
    stop("Explicit allocation counts require allocation = 'portfolio'.", call. = FALSE)
  }
  .svem_validate_allocation_counts(budget, n_optima, n_explore, priority)
  incumbent_summary <- match.arg(incumbent_summary)
  carry <- !identical(incumbent, FALSE)
  benchmark_row <- .svem_resolve_benchmark(benchmark, data)
  if (!is.null(benchmark_row) && allocation != "portfolio") {
    stop("A fixed benchmark requires allocation = 'portfolio'.")
  }
  if (carry && allocation != "portfolio") stop("An incumbent repeat requires allocation = 'portfolio'.")
  exploration_fraction <- .svem_numeric_scalar(exploration_fraction, "exploration_fraction",
                                               lower = 0, upper = 1, lower_open = TRUE)
  combine         <- match.arg(combine)
  numeric_sampler <- match.arg(numeric_sampler)

  # ---- basic validation ----
  if (inherits(objects, "svem_model")) objects <- list(objects)

  if (!is.list(objects) || !length(objects))
    stop("objects must be a nonempty list of svem_model objects.")
  if (!all(vapply(objects, inherits, logical(1), what = "svem_model")))
    stop("All elements of 'objects' must be svem_model objects.")

  n <- .svem_integer_scalar(n, "n", min = 1L)
  level <- .svem_numeric_scalar(level, "level", lower = 0, upper = 1,
                                lower_open = TRUE, upper_open = TRUE)
  verbose <- .svem_logical_scalar(verbose, "verbose")

  if (!is.null(data) && !is.data.frame(data)) {
    stop("'data' must be NULL or a data.frame.")
  }

  # ---- infer response names from model formulas ----
  lhs_names <- vapply(objects, function(o) {
    if (!is.null(o$formula) && inherits(o$formula, "formula")) {
      .svem_response_label(o$formula)
    } else {
      NA_character_
    }
  }, character(1))

  # ---- HARD CHECK: duplicate LHS names create ambiguous key mapping ----
  valid_lhs <- !is.na(lhs_names) & nzchar(lhs_names)
  if (anyDuplicated(lhs_names[valid_lhs])) {
    dups <- unique(lhs_names[valid_lhs][duplicated(lhs_names[valid_lhs])])
    stop(
      "Duplicate model response (LHS) names detected: ",
      paste(dups, collapse = ", "),
      ". Please rename the response variables and/or name the objects list uniquely."
    )
  }

  obj_names <- names(objects)
  if (is.null(obj_names)) obj_names <- rep("", length(objects))
  obj_names[is.na(obj_names)] <- ""

  empty_idx <- which(!nzchar(obj_names) & valid_lhs)
  if (length(empty_idx)) obj_names[empty_idx] <- lhs_names[empty_idx]
  names(objects) <- obj_names

  engines <- intersect(c("gaupro", "extra_trees"),
                       unique(unlist(lapply(objects, function(o) names(o$companions)))))
  if (!is.null(forward)) {
    if (!is.list(forward) || is.null(names(forward)) || anyDuplicated(names(forward)) ||
        !setequal(names(forward), names(objects)) ||
        !all(vapply(forward, inherits, logical(1), what = "svem_model"))) {
      stop("forward must be a named list of fitted models for exactly the same responses.")
    }
    forward <- forward[names(objects)]
    if (!all(vapply(forward, inherits, logical(1), what = "svem_forward"))) {
      stop("forward must contain Gaussian svem_forward() fits, not other SVEM engines.")
    }
    matching <- vapply(names(objects), function(r) {
      identical(.svem_response_label(forward[[r]]$formula),
                .svem_response_label(objects[[r]]$formula))
    }, logical(1))
    if (!all(matching)) stop("Forward response labels must match their corresponding primary model responses.")
    engines <- c(engines, "forward")
  }
  if (!is.null(spline)) {
    spline <- .svem_spline_model_list(spline, names(objects), labels = unname(lhs_names))
    .svem_spline_check_schema(spline, objects)
    engines <- c(engines, "spline")
  }
  make_worklist <- length(engines) > 0L || !is.null(budget) || !is.null(gp_endpoint) || allocation == "portfolio"
  if (!is.null(budget)) budget <- .svem_integer_scalar(budget, "budget", min = 1L)

  resp_names <- names(objects)
  if (any(!nzchar(resp_names))) {
    stop(
      "Could not infer response names for all models; please name the objects list ",
      "or ensure each model has a formula with a left-hand side."
    )
  }

  if (anyDuplicated(resp_names)) {
    dups <- unique(resp_names[duplicated(resp_names)])
    stop(
      "Response identifiers must be unique. Duplicates: ",
      paste(dups, collapse = ", "),
      ". Please rename the 'objects' list entries to unique names."
    )
  }

  # mapping: unambiguous object/LHS key -> canonical response name
  key_to_resp <- .svem_response_key_map(resp_names, lhs_names, valid_lhs)
  if (!is.null(gp_endpoint)) {
    if (!is.character(gp_endpoint) || length(gp_endpoint) != 1L ||
        !gp_endpoint %in% names(key_to_resp)) stop("Unknown gp_endpoint response.")
    gp_endpoint <- unname(key_to_resp[[gp_endpoint]])
  } else if ("gaupro" %in% engines && length(resp_names) == 1L) {
    gp_endpoint <- resp_names
  }

  # Generated names must never replace design factors or observed responses.
  # Check both canonical response names and formula labels before drawing RNG.
  generated_names <- .svem_generated_score_names(unique(c(
    resp_names, lhs_names[valid_lhs]
  )))
  predictor_names <- unique(unlist(lapply(objects, function(o) {
    o$sampling_schema$predictor_vars
  }), use.names = FALSE))
  if (make_worklist || !is.null(candidates)) {
    .svem_check_output_collisions(c(predictor_names, resp_names, generated_names),
                                  .svem_worklist_names(resp_names, engines),
                                  "svem_score_random() worklist")
  }
  .svem_check_output_collisions(predictor_names, generated_names,
                                "svem_score_random()")
  if (!is.null(data)) {
    response_vars <- unique(unlist(lapply(objects, function(o) {
      .svem_response_vars(o$formula)
    }), use.names = FALSE))
    .svem_check_output_collisions(
      intersect(names(data), c(predictor_names, response_vars)),
      generated_names, "svem_score_random() supplied data"
    )
    stale <- intersect(names(data), generated_names)
    if (length(stale)) {
      warning("svem_score_random(): replacing generated column(s) from a previous scoring round: ",
              paste(stale, collapse = ", "), ".", call. = FALSE)
      data <- data[, !names(data) %in% stale, drop = FALSE]
    }
  }

  if (!is.null(specs) && length(specs)) {
    if (!is.list(specs) || is.null(names(specs)) ||
        anyNA(names(specs)) || any(!nzchar(names(specs)))) {
      stop("`specs` must be a named list when provided.", call. = FALSE)
    }
    mapped <- names(specs)
    known <- mapped %in% names(key_to_resp)
    mapped[known] <- unname(key_to_resp[mapped[known]])
    if (anyDuplicated(mapped)) {
      stop("Spec names/aliases map to duplicate responses: ",
           paste(unique(mapped[duplicated(mapped)]), collapse = ", "),
           ". Supply at most one specification per response.", call. = FALSE)
    }
    names(specs) <- mapped
  }

  mismatch <- nzchar(obj_names) & valid_lhs & (obj_names != lhs_names)
  if (any(mismatch, na.rm = TRUE)) {
    warning(
      "names(objects) differ from model response names for: ",
      paste0(obj_names[mismatch], " (lhs: ", lhs_names[mismatch], ")", collapse = ", "),
      ". You can refer to these responses in goals/specs by either the list name or the response name."
    )
  }

  # ---- goals: allow unnamed (positional) or named by object name or LHS ----
  goals <- .svem_align_goals(goals, resp_names, key_to_resp)

  # ---- detect families ----
  resp_family <- vapply(objects, function(o) {
    fam <- tryCatch(o$family, error = function(e) NA_character_)
    if (!is.character(fam) || !length(fam)) NA_character_ else fam[1]
  }, character(1))
  is_binomial_resp <- tolower(resp_family) %in% "binomial"
  names(is_binomial_resp) <- resp_names

  if (!is.null(wmt) && any(is_binomial_resp, na.rm = TRUE)) {
    offending <- paste(resp_names[which(is_binomial_resp)], collapse = ", ")
    stop("Supplying `wmt` is not supported when any responses are binomial.\n",
         "Binomial responses detected: ", offending)
  }

  # ---- goals & weights ----
  goal_info <- .svem_validate_goals(goals, resp_names)
  goals <- goal_info$goals
  goal_df <- data.frame(
    response = resp_names,
    goal = unname(goal_info$goal),
    weight = unname(goal_info$weight),
    target = unname(goal_info$target),
    stringsAsFactors = FALSE
  )
  weights_user <- goal_info$normalized_weight

  # ---- optional WMT: extract from svem_wmt_multi result ----
  wmt_p_vals <- setNames(rep(NA_real_, length(resp_names)), resp_names)
  wmt_mult   <- setNames(rep(1.0,      length(resp_names)), resp_names)

  if (!is.null(wmt)) {
    if (!is.list(wmt) || is.null(wmt$multipliers)) {
      stop("`wmt` must be an object returned by svem_wmt_multi() ",
           "or at least a list with a named numeric 'multipliers' component.")
    }
    mult_vec <- wmt$multipliers
    if (is.null(names(mult_vec))) stop("`wmt$multipliers` must be a named numeric vector.")

    # allow WMT names keyed by LHS
    mult_names2 <- vapply(names(mult_vec), function(k) {
      if (k %in% names(key_to_resp)) key_to_resp[[k]] else k
    }, character(1))
    if (anyDuplicated(mult_names2)) {
      dups <- unique(mult_names2[duplicated(mult_names2)])
      stop("wmt$multipliers names map to duplicate responses: ", paste(dups, collapse = ", "))
    }
    names(mult_vec) <- mult_names2

    common <- intersect(names(mult_vec), resp_names)
    if (!length(common)) {
      stop("No overlap between response names and names(wmt$multipliers) after name alignment.")
    }
    # Surface partial mismatches: a misspelled multiplier name would otherwise
    # silently leave that response at full weight (1.0) and drop the stray entry
    no_mult <- setdiff(resp_names, names(mult_vec))
    if (length(no_mult)) {
      warning("No wmt multiplier found for response(s): ",
              paste(no_mult, collapse = ", "),
              "; using multiplier 1.0 (full weight) for them.", call. = FALSE)
    }
    stray_mult <- setdiff(names(mult_vec), resp_names)
    if (length(stray_mult)) {
      warning("Ignoring wmt multiplier name(s) that match no response: ",
              paste(stray_mult, collapse = ", "), call. = FALSE)
    }
    wmt_mult[common] <- as.numeric(mult_vec[common])

    if (any(!is.finite(wmt_mult)) || any(wmt_mult < 0)) {
      stop("wmt multipliers must be finite and nonnegative after alignment.")
    }

    if (!is.null(wmt$p_values)) {
      p_vec <- wmt$p_values
      if (!is.null(names(p_vec))) {
        p_names2 <- vapply(names(p_vec), function(k) {
          if (k %in% names(key_to_resp)) key_to_resp[[k]] else k
        }, character(1))
        if (anyDuplicated(p_names2)) {
          dups <- unique(p_names2[duplicated(p_names2)])
          stop("wmt$p_values names map to duplicate responses: ", paste(dups, collapse = ", "))
        }
        names(p_vec) <- p_names2
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

  weights_final_raw <- weights_user * wmt_mult
  sw2 <- sum(weights_final_raw)
  weights_final <- if (sw2 > 0 && is.finite(sw2)) weights_final_raw / sw2 else weights_user

  wts_user   <- stats::setNames(weights_user,  goal_df$response)
  wts_final  <- stats::setNames(weights_final, goal_df$response)
  wts_rank   <- if (reweight_flag) wts_final else wts_user
  wts_uncert <- wts_user  # ALWAYS user weights for uncertainty

  # ---- warn once if any models cannot produce bootstrap intervals ----
  has_coef_mat <- vapply(objects, function(o) {
    !is.null(o$coef_matrix) && is.matrix(o$coef_matrix) && nrow(o$coef_matrix) >= 2L
  }, logical(1))
  missing_ci_models <- names(objects)[!has_coef_mat]
  if (length(missing_ci_models)) {
    warning(
      "At least two bootstrap members are needed for intervals. CI columns will be NA for: ",
      paste(missing_ci_models, collapse = ", "),
      "; each missing response contributes a constant 0.5 times its user weight to uncertainty_measure.",
      call. = FALSE
    )
  }

  # ---- sample & predict ----
  # NOTE: svem_random_table_multi() auto-respects discrete numeric supports
  # stored in sampling_schema; no discrete argument is passed here.
  prepared <- NULL
  if (carry) {
    prepared <- .svem_observed_pool(objects, data, candidates,
      .svem_pool_groups(objects, mixture_groups), increments, incumbent, reference)
    candidates <- prepared$candidates
  }
  benchmark_prepared <- NULL
  if (!is.null(benchmark_row)) {
    benchmark_prepared <- tryCatch(.svem_observed_pool(objects, data, candidates,
      .svem_pool_groups(objects, mixture_groups), increments, benchmark_row, NULL),
      error = function(e) stop("Benchmark cannot be prepared: ", conditionMessage(e), call. = FALSE))
    candidates <- benchmark_prepared$candidates
  }
  sampled_raw <- svem_random_table_multi(
    objects         = objects,
    n               = n,
    mixture_groups  = mixture_groups,
    debias          = debias_flag,
    numeric_sampler = numeric_sampler,
    candidates      = candidates
  )

  if (!is.list(sampled_raw) || is.null(sampled_raw$data) || is.null(sampled_raw$pred)) {
    stop("svem_random_table_multi must return list(data, pred).")
  }

  data_df <- as.data.frame(sampled_raw$data)
  pred_df <- as.data.frame(sampled_raw$pred, check.names = FALSE, optional = FALSE)
  if (length(engines)) {
    deployment <- .svem_worklist_numeric_blocks(objects, data_df, !is.null(candidates))
    if (!identical(deployment, data_df)) {
      data_df <- deployment
      pred_df <- .svem_pool_predictions(objects, data_df, debias_flag)$pred
    }
  }
  indexed <- make_worklist || !is.null(candidates) || !is.null(increments)
  if (indexed) {
    pool_index <- .svem_pool_index(data_df, candidates, increments,
                                   if (make_worklist) objects[[1L]]$sampling_schema$blocking)
    data_df <- data_df[pool_index$keep, , drop = FALSE]
    pred_df <- pred_df[pool_index$keep, , drop = FALSE]
  }

  # map response labels -> prediction columns in sampled_raw$pred
  resp_table_cols <- setNames(character(length(resp_names)), resp_names)

  for (resp in resp_names) {
    obj <- objects[[resp]]

    lhs <- NA_character_
    if (!is.null(obj$formula) && inherits(obj$formula, "formula")) {
      lhs <- .svem_response_label(obj$formula)
    }

    matching_cols <- unique(c(if (!is.na(lhs) && nzchar(lhs)) paste0(lhs, "_pred"),
                              paste0(resp, "_pred")))

    found <- matching_cols[matching_cols %in% colnames(pred_df)]
    if (!length(found)) {
      stop(
        "Could not match response '", resp,
        "' to prediction columns from svem_random_table_multi(). Expected one of: ",
        paste(matching_cols, collapse = ", "),
        " in names(sampled_raw$pred)."
      )
    }
    resp_table_cols[resp] <- found[1L]
  }

  pred_mat <- as.data.frame(
    lapply(resp_names, function(r) as.numeric(pred_df[[resp_table_cols[[r]]]])),
    check.names = FALSE,
    optional    = FALSE
  )
  colnames(pred_mat) <- paste0(resp_names, "_pred")

  score_table <- cbind(data_df, pred_mat)
  if (indexed) score_table$recipe_id <- pool_index$ids

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
    y  <- score_table[[resp_pred_cols[[r]]]]

    if (isTRUE(is_binomial_resp[r])) y <- pmin(pmax(as.numeric(y), 0), 1)

    L_def <- .q(y, .q_lo)
    U_def <- .q(y, .q_hi)
    if (!is.finite(L_def)) L_def <- min(y, na.rm = TRUE)
    if (!is.finite(U_def)) U_def <- max(y, na.rm = TRUE)

    span_def <- U_def - L_def
    full_rng <- range(y, na.rm = TRUE)
    full_span <- if (!all(is.finite(full_rng))) NA_real_ else full_rng[2L] - full_rng[1L]

    span_floor <- .ds_span_abs_floor
    if (is.finite(full_span) && full_span > 0) span_floor <- max(span_floor, .ds_span_rel_floor * full_span)
    span_use <- if (!is.finite(span_def) || span_def <= 0) span_floor else max(span_def, span_floor)

    mid_def <- 0.5 * (L_def + U_def)
    L_def <- mid_def - 0.5 * span_use
    U_def <- mid_def + 0.5 * span_use

    # Explicitly supplied anchors that collapse the desirability band would
    # silently neutralize the response (constant desirability for every row);
    # surface the input error instead.
    .check_user_anchors <- function(L, U) {
      if ((!is.null(gi$lower_acceptable) || !is.null(gi$upper_acceptable)) &&
          (!is.finite(L) || !is.finite(U) || U <= L)) {
        stop("Desirability anchors for response '", r, "' are degenerate ",
             "(need finite lower_acceptable < upper_acceptable; got L = ", L,
             ", U = ", U, ").", call. = FALSE)
      }
    }

    if (g == "max") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      .check_user_anchors(L, U)
      z <- .ds_max(y, L, U, s)
      ds_params[[r]] <- list(type = "max", L = L, U = U, s = s)

    } else if (g == "min") {
      L <- if (!is.null(gi$lower_acceptable)) as.numeric(gi$lower_acceptable) else L_def
      U <- if (!is.null(gi$upper_acceptable)) as.numeric(gi$upper_acceptable) else U_def
      s <- if (!is.null(gi$shape))            as.numeric(gi$shape)            else 1
      .check_user_anchors(L, U)
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
        stop("Resolved target desirability anchors for response '", r,
             "' must satisfy finite L < target < U.", call. = FALSE)
      }
      z <- .ds_target(y, T0, L, U, sL, sR)
      ds_params[[r]] <- list(type = "target", T0 = T0, L = L, U = U, sL = sL, sR = sR)
    }

    contrib_cols[[paste0(r, "_des")]] <- z
  }

  for (nm in names(contrib_cols)) score_table[[nm]] <- contrib_cols[[nm]]

  # ---- combine desirabilities into score + wmt_score ----
  if (combine == "mean") {
    score_user <- Reduce(`+`, lapply(resp_cols, function(r) as.numeric(wts_user[r]) * contrib_cols[[paste0(r, "_des")]]))
    score_table$score <- score_user
    if (reweight_flag) {
      score_table$wmt_score <- Reduce(`+`, lapply(resp_cols, function(r) as.numeric(wts_final[r]) * contrib_cols[[paste0(r, "_des")]]))
    }
  } else {
    score_table$score <- exp(Reduce(`+`, lapply(resp_cols, function(r) {
      z <- contrib_cols[[paste0(r, "_des")]]
      zadj <- (1 - geom_floor) * z + geom_floor
      as.numeric(wts_user[r]) * log(zadj)
    })))
    if (reweight_flag) {
      score_table$wmt_score <- exp(Reduce(`+`, lapply(resp_cols, function(r) {
        z <- contrib_cols[[paste0(r, "_des")]]
        zadj <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts_rank[r]) * log(zadj)
      })))
    }
  }

  # ---- uncertainty_measure from CI widths + store per-response CIs ----
  # Anchors (q0.02/q0.98 of the CI widths) are computed once per response from
  # the sampled score_table and reused when scoring 'data', so that
  # uncertainty_measure in original_data_scored is on the same scale as in
  # score_table (mirroring how ds_params are reused for desirabilities).
  .ciw_anchors <- function(x) {
    a <- as.numeric(stats::quantile(x, probs = .q_lo, na.rm = TRUE, names = FALSE))
    b <- as.numeric(stats::quantile(x, probs = .q_hi, na.rm = TRUE, names = FALSE))
    if (!is.finite(a) || !is.finite(b) || b <= a) return(NULL)
    c(a, b)
  }

  .normalize01_robust <- function(x, anchors = NULL) {
    if (is.null(anchors)) anchors <- .ciw_anchors(x)

    # If we can't compute a spread (all NA, constant widths, etc.), use neutral 0.5
    if (is.null(anchors)) {
      return(rep(0.5, length(x)))
    }
    a <- anchors[1L]; b <- anchors[2L]

    out <- (pmin(pmax(x, a), b) - a) / (b - a)

    # Any NA/Inf -> neutral (this keeps uncertainty usable downstream)
    out[!is.finite(out)] <- 0.5
    out
  }

  ciw_anchors <- setNames(vector("list", length(resp_cols)), resp_cols)

  ciw_w_cols  <- list()
  ci_lwr_cols <- list()
  ci_upr_cols <- list()

  for (r in resp_cols) {
    obj <- objects[[r]]
    pr <- if (!has_coef_mat[[r]]) NULL else tryCatch(
      predict(obj,
              newdata  = score_table[, predictor_cols, drop = FALSE],
              debias   = debias_flag,
              interval = TRUE,
              level    = level),
      error = function(e) {
        warning("svem_score_random(): interval prediction failed for response '", r,
                "' in score_table: ", conditionMessage(e),
                "; CI columns are NA and its uncertainty contribution is neutral.",
                call. = FALSE)
        NULL
      }
    )

    lwr <- upr <- rep(NA_real_, nrow(score_table))
    if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$lwr) && !is.null(pr$upr)) {
      lwr <- as.numeric(pr$lwr)
      upr <- as.numeric(pr$upr)
    }

    if (isTRUE(is_binomial_resp[r])) {
      lwr <- pmin(pmax(lwr, 0), 1)
      upr <- pmin(pmax(upr, 0), 1)
    }

    ci_lwr_cols[[paste0(r, "_lwr")]] <- lwr
    ci_upr_cols[[paste0(r, "_upr")]] <- upr

    width <- upr - lwr
    ciw_anchors[[r]] <- .ciw_anchors(width)
    zc    <- .normalize01_robust(width, anchors = ciw_anchors[[r]])
    ciw_w_cols[[paste0(r, "_ciw_w")]] <- as.numeric(wts_uncert[r]) * zc
  }

  for (nm in names(ci_lwr_cols)) score_table[[nm]] <- ci_lwr_cols[[nm]]
  for (nm in names(ci_upr_cols)) score_table[[nm]] <- ci_upr_cols[[nm]]
  for (nm in names(ciw_w_cols))  score_table[[nm]] <- ciw_w_cols[[nm]]

  ciw_mat <- as.data.frame(ciw_w_cols, check.names = FALSE)
  if (ncol(ciw_mat)) {
    score_table$uncertainty_measure <- rowSums(ciw_mat, na.rm = TRUE)
  } else {
    score_table$uncertainty_measure <- rep(0.5, nrow(score_table))  # or NA_real_
  }


  # ---- score original data (optional) ----
  original_data_scored <- NULL
  if (!is.null(data) && is.data.frame(data)) {

    pred_df0 <- data.frame(row_id = seq_len(nrow(data)))
    ci_lwr  <- list()
    ci_upr  <- list()

    for (r in resp_names) {
      pr <- if (!has_coef_mat[[r]]) NULL else tryCatch(
        predict(objects[[r]], newdata = data, debias = debias_flag,
                interval = TRUE, level = level),
        error = function(e) {
          warning("svem_score_random(): interval prediction failed for response '", r,
                  "' in supplied data: ", conditionMessage(e), ".", call. = FALSE)
          NULL
        }
      )

      fit_col <- if (!inherits(pr, "try-error") && is.list(pr) && !is.null(pr$fit)) {
        as.numeric(pr$fit)
      } else {
        tmp <- tryCatch(predict(objects[[r]], newdata = data, debias = debias_flag),
                        error = function(e) {
                          warning("svem_score_random(): point prediction failed for response '", r,
                                  "' in supplied data: ", conditionMessage(e),
                                  "; returning NA predictions.", call. = FALSE)
                          NULL
                        })
        if (is.null(tmp)) rep(NA_real_, nrow(data)) else as.numeric(tmp)
      }

      if (isTRUE(is_binomial_resp[r])) fit_col <- pmin(pmax(fit_col, 0), 1)
      pred_df0[[r]] <- fit_col

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
      y <- pred_df0[[r]]
      p <- ds_params[[r]]
      if (p$type == "max") {
        des_cols[[paste0(r, "_des")]] <- .ds_max(y, p$L, p$U, p$s)
      } else if (p$type == "min") {
        des_cols[[paste0(r, "_des")]] <- .ds_min(y, p$L, p$U, p$s)
      } else {
        des_cols[[paste0(r, "_des")]] <- .ds_target(y, p$T0, p$L, p$U, p$sL, p$sR)
      }
    }

    if (combine == "mean") {
      score_user <- Reduce(`+`, lapply(resp_names, function(r) as.numeric(wts_user[r]) * des_cols[[paste0(r, "_des")]]))
      score_rank <- score_user
      if (reweight_flag) {
        score_rank <- Reduce(`+`, lapply(resp_names, function(r) as.numeric(wts_final[r]) * des_cols[[paste0(r, "_des")]]))
      }
    } else {
      score_user <- exp(Reduce(`+`, lapply(resp_names, function(r) {
        z  <- des_cols[[paste0(r, "_des")]]
        zf <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts_user[r]) * log(zf)
      })))
      score_rank <- exp(Reduce(`+`, lapply(resp_names, function(r) {
        z  <- des_cols[[paste0(r, "_des")]]
        zf <- (1 - geom_floor) * z + geom_floor
        as.numeric(wts_rank[r]) * log(zf)
      })))
    }

    ciw_w_orig <- list()
    for (r in resp_names) {
      width <- ci_upr[[r]] - ci_lwr[[r]]
      # reuse the anchors from the sampled score_table for comparability;
      # if the score_table anchors were degenerate (all-NA/constant widths),
      # the score_table values were neutral 0.5, so match that here rather
      # than recomputing local anchors
      zc <- if (is.null(ciw_anchors[[r]])) {
        rep(0.5, length(width))
      } else {
        .normalize01_robust(width, anchors = ciw_anchors[[r]])
      }
      ciw_w_orig[[paste0(r, "_ciw_w")]] <- as.numeric(wts_uncert[r]) * zc
    }

    ciw_mat0 <- as.data.frame(ciw_w_orig, check.names = FALSE)
    unc_vec <- if (ncol(ciw_mat0)) rowSums(ciw_mat0, na.rm = TRUE) else rep(0.5, nrow(data))


    pred_mat <- pred_df0[resp_names]
    colnames(pred_mat) <- paste0(resp_names, "_pred")

    original_data_scored <- cbind(
      data,
      pred_mat,
      as.data.frame(des_cols, stringsAsFactors = FALSE)
    )

    for (r in resp_names) {
      original_data_scored[[paste0(r, "_lwr")]]   <- ci_lwr[[r]]
      original_data_scored[[paste0(r, "_upr")]]   <- ci_upr[[r]]
      original_data_scored[[paste0(r, "_ciw_w")]] <- ciw_mat0[[paste0(r, "_ciw_w")]]
    }

    original_data_scored$score <- score_user
    if (reweight_flag) original_data_scored$wmt_score <- score_rank
    original_data_scored$uncertainty_measure <- unc_vec
  }

  # ---- optional spec-limit augmentation (safe) ----
  if (!is.null(specs)) {

    if (!is.null(names(specs))) {
      s_names <- names(specs)
      remapped <- s_names
      for (i in seq_along(s_names)) {
        key <- s_names[i]
        if (!key %in% resp_names && key %in% names(key_to_resp)) {
          remapped[i] <- key_to_resp[[key]]
        }
      }
      names(specs) <- remapped
    }

    score_table <- tryCatch(
      svem_append_design_space_cols(score_table = score_table, objects = objects, specs = specs),
      error = function(e) {
        # a warning (not a verbose-only message): the user asked for spec
        # columns and would otherwise discover their absence only downstream
        warning(
          "svem_score_random(): spec-limit augmentation failed for score_table; returning without spec columns. Error: ",
          conditionMessage(e), call. = FALSE
        )
        score_table
      }
    )

    if (!is.null(original_data_scored) && is.data.frame(original_data_scored)) {
      original_data_scored <- tryCatch(
        svem_append_design_space_cols(score_table = original_data_scored, objects = objects, specs = specs),
        error = function(e) {
          warning(
            "svem_score_random(): spec-limit augmentation failed for original_data_scored; returning without spec columns. Error: ",
            conditionMessage(e), call. = FALSE
          )
          original_data_scored
        }
      )
    }
  }

  # store attributes for downstream helpers
  attr(score_table, "svem_predictor_cols") <- predictor_cols
  attr(score_table, "svem_resp_cols")      <- resp_cols
  attr(score_table, "svem_response_vars") <- unique(unlist(lapply(
    objects, function(o) .svem_response_vars(o$formula)
  ), use.names = FALSE))
  if (!is.null(original_data_scored) && is.data.frame(original_data_scored)) {
    # restrict to the model's design factors: colnames(data) would include the
    # measured responses and any ID/notes columns, silently contaminating
    # Gower-distance diversity clustering in svem_select_from_score_table()
    attr(original_data_scored, "svem_predictor_cols") <-
      intersect(predictor_cols, colnames(original_data_scored))
    attr(original_data_scored, "svem_resp_cols")      <- resp_cols
    attr(original_data_scored, "svem_response_vars") <-
      attr(score_table, "svem_response_vars")
  }

  out <- list(
    score_table          = score_table,
    original_data_scored = original_data_scored,
    weights_original     = wts_user,
    weights_final        = wts_final,
    wmt_p_values         = if (reweight_flag) wmt_p_vals else NULL,
    wmt_multipliers      = if (reweight_flag) wmt_mult   else NULL
  )
  if (!make_worklist) return(out)
  des <- function(y, p) {
    switch(p$type,
           max = .ds_max(y, p$L, p$U, p$s),
           min = .ds_min(y, p$L, p$U, p$s),
           target = .ds_target(y, p$T0, p$L, p$U, p$sL, p$sR))
  }
  cross <- .svem_engine_scores(out$score_table, data_df, objects, setdiff(engines, c("forward", "spline")),
                               ds_params, des, wts_user, combine, geom_floor, gp_endpoint)
  out$score_table <- cross$table
  out$engine_status <- cross$status
  if (!is.null(forward)) {
    fw <- .svem_forward_scores(out$score_table, data_df, forward, ds_params, des,
                                wts_user, combine, geom_floor)
    out$score_table <- fw$table
    out$engine_status <- rbind(out$engine_status, fw$status)
  }
  if (!is.null(spline)) {
    sm <- .svem_spline_scores(out$score_table, data_df, spline, ds_params, des,
                              wts_user, combine, geom_floor)
    out$score_table <- sm$table
    out$engine_status <- rbind(out$engine_status, sm$status)
  }
  incumbent_info <- NULL
  if (carry) {
    incumbent_info <- .svem_observed_choice(objects, data, prepared, ds_params, des,
      wts_user, combine, geom_floor, goals, pool_index$keys, incumbent_summary)
    incumbent_info$recipe_id <- out$score_table$recipe_id[incumbent_info$index]
    incumbent_info$by_context$recipe_id <- as.character(out$score_table$recipe_id[
      match(incumbent_info$by_context$recipe_key, pool_index$keys)])
    if (isTRUE(verbose) && incumbent_info$higher_raw_in_other_context) {
      message("A higher raw observed value exists in another block context. ",
              "Raw values across blocks are not automatically comparable; ",
              "inspect $incumbent$by_context or choose an explicit incumbent row.")
    }
  }
  benchmark_info <- NULL
  if (!is.null(benchmark_row)) {
    benchmark_info <- list(index = match(benchmark_prepared$keys[1L], pool_index$keys),
      source_row = benchmark_row, source_data = data[benchmark_row, , drop = FALSE],
      source_context = .svem_context_label(data[benchmark_row, , drop = FALSE],
        objects[[1L]]$sampling_schema$blocking), selector = benchmark)
    if (is.na(benchmark_info$index)) stop("Fixed benchmark was not retained in the candidate pool.")
    benchmark_info$recipe_id <- out$score_table$recipe_id[benchmark_info$index]
  }
  history <- data
  history_source <- "supplied_data"
  if (is.null(history)) {
    tables <- lapply(objects, function(obj) attr(obj$companions, "prepared_recipes", exact = TRUE))
    available <- vapply(tables, is.data.frame, logical(1))
    history_source <- if (all(available)) "fitted_preparations" else if (any(available)) {
      "partial_fitted_preparations"
    } else "no_preparation_history"
    if (any(available)) history <- unique(do.call(rbind, tables[available]))
  }
  out$history_source <- history_source
  work <- .svem_worklist(out$score_table, history, pool_index$keys, increments,
                         objects[[1L]]$sampling_schema, engines, goals, wts_user,
                         wts_uncert, ciw_anchors, gp_endpoint, budget,
                         allocation = allocation, priority = priority,
                         exploration_fraction = exploration_fraction,
                         incumbent = if (carry) incumbent_info$index else NULL,
                         benchmark = if (!is.null(benchmark_info)) benchmark_info$index else NULL,
                         n_optima = n_optima, n_explore = n_explore)
  if (history_source == "partial_fitted_preparations") {
    verified <- work$nominations$best_prepared
    work$nominations$best_prepared[!is.na(verified) & !verified] <- NA
  }
  if (allocation == "portfolio") {
    if (carry) work$nominations$best_ties[work$nominations$engine == "incumbent"] <- incumbent_info$ties
    out$incumbent <- incumbent_info
    out$benchmark <- benchmark_info
    out$desirability_parameters <- ds_params
    out$uncertainty_anchors <- ciw_anchors
    context <- .svem_context_label(data_df, objects[[1L]]$sampling_schema$blocking)[1L]
    work$worklist$prediction_context <- rep(context, nrow(work$worklist))
    work$worklist$reference_context <- rep(NA_character_, nrow(work$worklist))
    work$worklist$incumbent_source_rows <- rep(NA_character_, nrow(work$worklist))
    work$worklist$incumbent_n_preparations <- rep(NA_integer_, nrow(work$worklist))
    work$worklist$incumbent_recipe_n_preparations <- rep(NA_integer_, nrow(work$worklist))
    work$worklist$incumbent_recipe_n_complete <- rep(NA_integer_, nrow(work$worklist))
    work$worklist$incumbent_recipe_mean_score <- rep(NA_real_, nrow(work$worklist))
    work$worklist$incumbent_summary <- rep(NA_character_, nrow(work$worklist))
    work$worklist$benchmark_source_row <- rep(NA_integer_, nrow(work$worklist))
    work$worklist$benchmark_source_context <- rep(NA_character_, nrow(work$worklist))
    if (carry) {
      i <- which(work$worklist$role %in% c("incumbent_repeat", "incumbent_benchmark_repeat"))
      work$worklist$reference_context[i] <- incumbent_info$source_context
      work$worklist$incumbent_source_rows[i] <- paste(incumbent_info$source_rows, collapse = ",")
      work$worklist$incumbent_n_preparations[i] <- incumbent_info$n_preparations
      work$worklist$incumbent_recipe_n_preparations[i] <- incumbent_info$recipe_n_preparations
      work$worklist$incumbent_recipe_n_complete[i] <- incumbent_info$recipe_n_complete
      work$worklist$incumbent_recipe_mean_score[i] <- incumbent_info$recipe_mean_score
      work$worklist$incumbent_summary[i] <- incumbent_info$summary
    }
    if (!is.null(benchmark_info)) {
      i <- which(work$worklist$recipe_id == benchmark_info$recipe_id)
      work$worklist$benchmark_source_row[i] <- benchmark_row
      work$worklist$benchmark_source_context[i] <- benchmark_info$source_context
    }
  }
  c(out, work)
}
