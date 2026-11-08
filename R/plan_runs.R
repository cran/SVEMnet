#' Plan the default formulation follow-up workflow
#'
#' Fit complementary surfaces, carry forward the best observed preparation,
#' and allocate a fixed preparation budget. This is the package's default
#' formulation workflow, composed from [svem_fit_responses()], [svem_forward()]
#' and [svem_score_random()]. Existing individual fitting/scoring calls retain
#' their behavior. Files are written only by an explicit call to
#' [svem_export_candidates_csv()].
#'
#' @param spec Shared [bigexp_terms()] specification.
#' @param data Measured data, one row per independent preparation. Summarize
#'   technical readings within preparation before calling. Use response columns
#'   on the scale on which the scientific goals are defined.
#' @param goals Named response goals as in [svem_score_random()]. Their names
#'   are response columns in `data`. Fixed scientific desirability anchors are
#'   recommended so decisions across rounds have the same interpretation.
#' @param candidates Explicit feasible recipe pool, including the planned
#'   operator/day values and any `recipe_id` labels. Feasible historical recipes
#'   may be appended for the incumbent repeat, without rounding their settings.
#' @param budget Total physical preparation budget, including required repeats.
#'   With both `n_optima` and `n_explore` supplied, omit `budget` to use
#'   their sum plus the distinct required reference repeats. With only one
#'   count, supply a budget to allocate the remaining capacity. With neither
#'   count, the default budget is five preparations.
#' @param fits Optional ordinary named list of already fitted primary SVEM
#'   models with any requested companions. Otherwise fitted sequentially.
#' @param forward Optional named list of forward fits. Otherwise fitted for
#'   each response, using the Gaussian validation-likelihood default.
#' @param fit_control Named arguments passed to [svem_fit_responses()]. Defaults
#'   request `companions = c("gaupro", "extra_trees")` with
#'   `companion_blocks = "interacting"`. Optional engines need the suggested
#'   packages; an unavailable fit remains visible in the scoring diagnostics.
#'   Supply `companions = FALSE` for a core-only workflow. All primary SVEM
#'   algorithm defaults, including relaxation and penalized blocks, are retained.
#' @param forward_control Named additional arguments passed to each
#'   [svem_forward()] fit. Explicit choices such as `objective = "wAIC"` are
#'   preserved. `nBoot` controls its members independently of `fit_control`.
#' @param score_control Named additional arguments to [svem_score_random()],
#'   such as `combine`, `increments` and `mixture_groups`. WMT and specification
#'   probabilities are not computed by default.
#' @param incumbent `TRUE` carries the best observed preparation forward;
#'   `FALSE` omits the repeat. An integer chooses a historical row explicitly.
#' @param benchmark Optional fixed benchmark to carry forward in addition to
#'   the incumbent: one historical row number, or a named list with one unique
#'   preparation identifier, for example `list(RunID = "02SEP22-2")`.
#'   The identifier form stays meaningful when data rows are reordered or added.
#'   A requested benchmark is mandatory and consumes a preparation slot. If it
#'   is the same physical recipe as the incumbent, they share a single repeat.
#' @param spline Optional named response list from [svem_spline()], `TRUE` to fit
#'   the optional smoother, or `FALSE` to disable it. Default `NULL` preserves the
#'   ordinary five-role workflow. For a newly fitted plan with explicit counts
#'   and default priority, spline follows forward. The planner reports cases
#'   outside the declared surface domain or its resource budget. Explicit
#'   `spline = TRUE`, a spline priority, or `n_optima = 5` enables the general
#'   surface for numeric, nominal and mixture/process factors. Unmeasured block
#'   levels have no spline prediction. A pool fixing a numeric coordinate that
#'   varied in the data can make the spline unavailable. For pools above 65,536
#'   distinct numeric sites, supply saved spline fits built on a smaller
#'   response-free reference grid. Smaller count-only
#'   calls retain the former automatic numeric-spline eligibility.
#'   When the planner fits the spline,
#'   `candidates` supplies its response-free reference grid. This grid determines
#'   the basis and can limit its available rank. The declared physical bounds
#'   separately determine recipe eligibility. Choose the reference grid before
#'   observing outcomes and keep it fixed across refits. To score a different
#'   pool with the same basis, pass saved fits from
#'   `svem_fit_responses(..., method = "spline", domain_grid = reference_grid)`.
#'   Declare mixture groups in the shared [bigexp_terms()] specification when
#'   asking the planner to fit the spline. Numeric block values denote categorical
#'   levels. A declared continuous companion block makes the spline unavailable.
#' @param spline_control Named spline fitting controls (`nBoot`, `lambda`,
#'   `store_member_weights`) and the advanced basis controls documented in
#'   [spline_terms()], including `rank`, `null_weight`, `level_modes`,
#'   `level_penalty_ratio`, `max_columns`, `interaction_order`,
#'   `marginal_modes` and `marginal_weight`. A
#'   nonempty control list requires spline fitting to be enabled.
#' @param n_optima Optional count of distinct new model-nominated preparations,
#'   excluding incumbent and benchmark repeats. If five models nominate five
#'   distinct new recipes, they and one incumbent make six preparations.
#'   Duplicate nominations share one preparation and the next model is considered.
#'   With `n_optima = NULL`, available models fill capacity left after repeats
#'   and reserved exploration. An unfilled automatic model place is capacity
#'   rather than a requested-model shortfall.
#' @param n_explore Optional exact exploration-medoid count, including zero.
#'   Default `NULL` fills remaining capacity. This count is distinct from
#'   `exploration_fraction`, which controls the uncertainty shortlist. Requested
#'   counts and required repeats must fit the total `budget`; unused capacity
#'   and genuine model/exploration shortages remain inspectable.
#' @param incumbent_summary `"best"` (default) retains the best observed
#'   individual preparation. `"mean"` first averages repeated independent
#'   preparations of the same recipe within the selected reference context.
#'   Technical readings must already be summarized within preparation.
#' @param reference Historical block context for automatic incumbent selection;
#'   see [svem_score_random()]. Defaults to the planned block levels. Raw
#'   observations from different days/operators are not compared automatically.
#' @param priority Ordered roles, default incumbent, elastic net, GauPro,
#'   Extra Trees, then forward. Omit a role or reorder explicitly to customize.
#'   Listing `"spline"` requests that fit and conflicts with `spline = FALSE`.
#'   With `spline = TRUE`, include `"spline"` in any custom priority.
#' @param exploration_fraction Top fraction of eligible SVEM uncertainty scores
#'   supplied to medoid selection. Default 0.1; 0.25 is a wider sensitivity.
#'
#' @details Five slots normally provide an incumbent repeat and up to four
#'   model optima. Each optimum is the best eligible *unprepared* recipe from
#'   that model. Shared nominees consume one slot and retain all endorsements.
#'   A requested fixed benchmark is reserved before optional model nominations;
#'   the required incumbent and benchmark must both fit within the budget.
#'   Definite budget shortfalls are rejected before fitting. When automatic
#'   incumbent selection needs fitted desirability anchors, the exact required
#'   number is checked after scoring.
#'   Their shared recipe has role `incumbent_benchmark_repeat`; a distinct fixed
#'   benchmark has role `benchmark_repeat`. Remaining slots are diverse medoids
#'   of high ensemble spread. That spread
#'   is a useful exploration diagnostic, not automatically calibrated uncertainty.
#'   If a model also endorses the incumbent, its unrestricted nomination records
#'   that fact. No maximum-GP-uncertainty slot is reserved in this allocation.
#'
#'   After required repeats, the budget truncates the optional model priority.
#'   Budget 1 with the incumbent enabled and no distinct fixed benchmark
#'   repeats the incumbent and produces zero new recipes.
#'   This nested order is a practical default, not a universally optimal ranking.
#'   A scientist needing one new recipe can set `incumbent = FALSE`.
#'   If nominations coincide, unavailable engines leave gaps, or the budget is
#'   larger, exploration fills the remaining slots when informative spread and
#'   eligible recipes are available. A shortfall is reported without inventing
#'   settings. Preparation IDs should be assigned by the laboratory on execution;
#'   `recipe_id` identifies settings within the pool. Preserve supplied stable
#'   IDs and historical recipes in subsequent pools to carry these labels across
#'   rounds; otherwise match physical settings, not automatically generated IDs.
#'
#'   The incumbent is the best *observed* preparation, not necessarily the best
#'   true recipe. Its selection uses measured outcomes with the same goal maps
#'   and weights as candidate scoring. Missing required outcomes are ineligible.
#'   Ties use a stable recipe key. For repeated independent preparations this
#'   intentionally retains the best observed row, even after a poorer repeat.
#'   Set `incumbent_summary = "mean"` to average independent preparations of
#'   each recipe within the reference context before scoring, or explicitly
#'   choose a historical row. Saturated multi-response goals remain exact
#'   desirability ties, resolved by a stable recipe key; no extra utility beyond
#'   the declared goal limits is assumed. Inspect `best_ties` and `new_ties` in
#'   the nominations: shared nominations resolved by the same tie key are not
#'   evidence of independent model agreement.
#'
#'   Refit by calling this function on the augmented measured data with the next
#'   explicit pool. The planner checks the measured predictor levels against
#'   `spec` before fitting. If a new measured day, operator or ingredient level
#'   is absent from the stored specification, rebuild it with
#'   `bigexp_terms(..., data = augmented_data)`, keeping the intended formula,
#'   expansion settings and constraints, before calling the planner again.
#'   This prevents an old specification from mapping new measured levels to
#'   missing values and dropping those preparations. Missing response values
#'   remain subject to the individual endpoint fits' usual missing-data rules.
#'   To rescore different scientific weights or allocations without refitting,
#'   reuse `fits` and `forward` and pass `spline = plan$spline` if the earlier
#'   plan fitted a spline. With both primary and forward fits supplied, changing
#'   counts alone does not request a fresh spline. Request it explicitly with
#'   `spline = TRUE` or a custom priority containing `"spline"`. Reused fits
#'   draw no random numbers. Supplied models must reproduce the complete
#'   training response/design rows, allowing row permutations and each endpoint's
#'   missing-response pattern. The predictors, nominal factors and block roles
#'   must agree with `spec`. A declared discrete support must match the fitted
#'   support. A declared mixture geometry must match any fitted geometry and is
#'   required at fit time for companions. For fits without companions and
#'   without fitted mixture geometry, the declared geometry constrains scoring.
#'   Changed measurements require refitting. Older objects without verifiable
#'   stored training rows cannot be reused here.
#'   Check plausible weight alternatives when they
#'   change scientific priorities; the model portfolio does not remove that
#'   sensitivity. The planner does not automatically add a round or batch factor.
#'   Only explicitly declared blocking factors enter the model. New measured
#'   levels require a rebuilt specification and refit; primary SVEM represents
#'   them as additive penalized effects without recipe-by-block interactions.
#'   Companion blocks may interact with recipe factors. A repeated fixed
#'   benchmark can connect contexts, and exploration extends recipe coverage,
#'   but neither guarantees separation of block and recipe effects. Compare new
#'   recipes with repeats prepared in the same batch; under an additive shift,
#'   their contrast removes that shift without relying on a fitted block term.
#'   Declare blocks from physical conditions before seeing outcomes. A fixed
#'   benchmark chosen in advance avoids selecting that reference for an unusually
#'   high observed response, but its repeat measurements are still noisy. Deliberate
#'   challenge preparations to revisit poor responses remain a manual choice.
#'
#'   GauPro contributes fitted means or joint desirability within this custom
#'   confirmation workflow. This is not a full Bayesian optimizer. For Gaussian
#'   process optimization and acquisition policies see Jones, Schonlau and
#'   Welch (1998), \doi{10.1023/A:1008306431147}, and Frazier (2018),
#'   \url{https://arxiv.org/abs/1807.02811}.
#'   The living vignette extends the workflow and SVEMnet papers at
#'   \url{https://arxiv.org/abs/2212.11264} and
#'   \url{https://arxiv.org/abs/2511.20968}. Those author-maintained versions
#'   may be updated; published journal versions remain static.
#'
#' @return A `svem_run_plan` with ordinary `fits` and `forward` model lists,
#'   full `scores`, `worklist`, `nominations`, observed `incumbent` and optional
#'   `benchmark` provenance, `shortfall` and `policy`. A fitted or supplied spline
#'   adds `spline`. Explicit counts add `allocation_counts` and `merged_nominations`.
#'   Printing shows the achieved budget and roles.
#'   Plotting returns a score/spread view or endpoint prediction panels.
#' @examples
#' d <- data.frame(x = seq(-1, 1, length.out = 16))
#' d$y <- 2 + d$x - d$x^2 + 0.05 * sin(seq_len(nrow(d)))
#' spec <- bigexp_terms(y ~ x, d, factorial_order = 1,
#'                     polynomial_order = 2, report = FALSE)
#' goals <- list(y = list(goal = "max", weight = 1, lower_acceptable = 0, upper_acceptable = 3))
#' pool <- data.frame(x = seq(-1, 1, length.out = 41))
#' set.seed(7)
#' plan <- svem_plan_runs(spec, d, goals, pool, budget = 5,
#'   fit_control = list(nBoot = 3, companions = FALSE),
#'   forward_control = list(nBoot = 3), score_control = list(verbose = FALSE))
#' plan
#' plan$worklist[, c("recipe_id", "slot", "role", "nominated_by", "new_recipe")]
#' # Explicit file writing: svem_export_candidates_csv(plan, file = path)
#' @export
svem_plan_runs <- function(spec, data, goals, candidates, budget = 5L,
                           fits = NULL, forward = NULL,
                           fit_control = list(), forward_control = list(),
                           score_control = list(), incumbent = TRUE, reference = NULL,
                           priority = c("incumbent", "svemnet", "gaupro", "extra_trees", "forward"),
                           exploration_fraction = .1,
                           incumbent_summary = c("best", "mean"), benchmark = NULL,
                           spline = NULL, spline_control = list(),
                           n_optima = NULL, n_explore = NULL) {
  priority_missing <- missing(priority)
  if (!inherits(spec, "bigexp_spec")) stop("spec must come from bigexp_terms().")
  if (!is.data.frame(data)) stop("data must be a measured preparation table.")
  if (!is.list(goals) || !length(goals) || is.null(names(goals)) ||
      anyNA(names(goals)) || any(!nzchar(names(goals))) || anyDuplicated(names(goals)) ||
      !all(names(goals) %in% names(data))) {
    stop("goals must be named by unique response columns in data.")
  }
  reserved_columns <- c("new_recipe", "prediction_context", "reference_context",
    "incumbent_source_rows", "incumbent_n_preparations", "incumbent_summary",
    "incumbent_recipe_n_preparations", "incumbent_recipe_n_complete",
    "incumbent_recipe_mean_score", "benchmark_source_row", "benchmark_source_context")
  collisions <- intersect(reserved_columns, c(spec$vars, names(goals)))
  if (length(collisions)) stop("Reserved portfolio worklist name: ",
                               paste(collisions, collapse = ", "), call. = FALSE)
  prepared_data <- tryCatch(bigexp_prepare(spec, data, unseen = "error")$data, error = function(e) {
    stop("Measured data are incompatible with `spec`: ", conditionMessage(e),
         ". If new factor or blocking levels have been measured, rebuild the ",
         "specification with bigexp_terms(..., data = augmented_data) before ",
         "planning again.", call. = FALSE)
  })
  # Keep the caller's valid factor, character and logical representations.
  # Numeric labels for a declared nominal factor must reach the fit as nominal
  # too, so the saved sampling classes agree with the expanded design.
  nominal_labels <- spec$vars[spec$is_cat[spec$vars] &
    vapply(data[spec$vars], is.numeric, logical(1))]
  data[nominal_labels] <- prepared_data[nominal_labels]
  incumbent_summary <- match.arg(incumbent_summary)
  if (!is.null(spline) && !identical(spline, TRUE) && !identical(spline, FALSE) && !is.list(spline)) {
    stop("spline must be NULL, TRUE, FALSE, or a named list of spline fits.", call. = FALSE)
  }
  explicit_counts <- !is.null(n_optima) || !is.null(n_explore)
  if (missing(budget) && explicit_counts) budget <- NULL
  budget <- .svem_integer_scalar(budget, "budget", min = 1L, allow_null = explicit_counts)
  reuse_only <- !is.null(fits) && !is.null(forward) && is.null(spline)
  priority <- .svem_resolve_plan_priority(priority, priority_missing && !reuse_only,
                                          spline, n_optima, n_explore)
  allowed <- c("incumbent", "svemnet", "gaupro", "extra_trees", "forward", "spline")
  if (!is.character(priority) || !length(priority) || anyNA(priority) ||
      anyDuplicated(priority) || any(!priority %in% allowed)) {
    stop("priority must contain distinct supported role names.")
  }
  if (!priority_missing && identical(spline, FALSE) && "spline" %in% priority) {
    stop("spline = FALSE conflicts with an explicit spline role in priority. Remove that role or enable the spline.", call. = FALSE)
  }
  if (isTRUE(spline) && !"spline" %in% priority) {
    stop("spline = TRUE requires 'spline' in priority; add that role or omit the spline request.", call. = FALSE)
  }
  if (reuse_only && priority_missing && !is.null(n_optima)) {
    requested <- .svem_integer_scalar(n_optima, "n_optima", min = 0L)
    if (requested == length(setdiff(priority, "incumbent")) + 1L) {
      stop("n_optima exceeds the number of prioritized model roles during fit reuse. ",
           "A fifth model candidate requires a supported spline. Pass spline = previous_plan$spline ",
           "to reuse it, or spline = TRUE to fit one. Otherwise request at most four model ",
           "candidates and use n_explore for more preparations.", call. = FALSE)
    }
  }
  .svem_validate_allocation_counts(budget, n_optima, n_explore, priority)
  exploration_fraction <- .svem_numeric_scalar(exploration_fraction,
    "exploration_fraction", lower = 0, upper = 1, lower_open = TRUE)
  reserved <- c("formula", "spec", "data", "response", "responses")
  fit_control <- .svem_plan_control(fit_control, c(reserved, "method", "domain_grid"), "fit_control")
  forward_control <- .svem_plan_control(forward_control, reserved, "forward_control")
  spline_control <- .svem_plan_control(spline_control,
    c(reserved, "domain_grid", "blocking", "mixture_groups", "bounds"), "spline_control")
  .svem_spline_check_control(spline_control, .svem_spline_dimension(spec),
    original = .svem_spline_original_domain(spec, candidates, spline_control))
  score_control <- .svem_plan_control(score_control,
    c("objects", "data", "goals", "candidates", "budget", "forward", "spline", "allocation",
      "incumbent", "reference", "priority", "exploration_fraction", "incumbent_summary", "benchmark",
      "n_optima", "n_explore"), "score_control")
  if (!is.null(fit_control$family) && !identical(fit_control$family, "gaussian")) {
    stop("The default formulation workflow requires Gaussian response fits.")
  }
  .svem_plan_preflight(spec, data, candidates, incumbent, reference, score_control,
    benchmark = benchmark, budget = budget, priority = priority,
    n_optima = n_optima, n_explore = n_explore)
  if (!is.null(fits)) {
    if (length(fit_control)) stop("fit_control cannot modify supplied fits; refit or omit fit_control.")
    fits <- .svem_plan_model_list(fits, names(goals), "fits")
    .svem_plan_check_training(fits, data, "fits", spec)
  }
  if (!is.null(forward)) {
    if (length(forward_control)) stop("forward_control cannot modify supplied fits; refit or omit forward_control.")
    forward <- .svem_plan_model_list(forward, names(goals), "forward")
    .svem_plan_check_training(forward, data, "forward", spec)
  }
  supplied_spline <- is.list(spline)
  if (supplied_spline) {
    if (length(spline_control)) stop("spline_control cannot modify supplied spline fits; refit or omit spline_control.")
    spline <- .svem_spline_model_list(spline, names(goals), data, check_training = TRUE)
    for (r in names(spline)) {
      ss <- spline[[r]]$spec
      declared_groups <- if (!is.null(fits)) fits[[r]]$sampling_schema$mixture_groups else
        if (!is.null(fit_control$mixture_groups)) fit_control$mixture_groups else spec$settings$mixture_groups
      if (!isTRUE(all.equal(.svem_mixture_groups(ss$mixture_groups),
                            .svem_mixture_groups(declared_groups), check.attributes = FALSE))) {
        stop("Spline and primary fits must use the same mixture constraints.", call. = FALSE)
      }
      if (!is.null(ss) && (!setequal(c(ss$predictors, ss$blocking), spec$vars) ||
          !setequal(ss$blocking, spec$settings$blocking))) {
        stop("Supplied spline fits must use the full declared factor and blocking specification.", call. = FALSE)
      }
    }
  }
  want_spline <- !identical(spline, FALSE) && "spline" %in% priority &&
    (isTRUE(spline) || supplied_spline || (explicit_counts && !reuse_only) || !priority_missing)
  if (length(spline_control) && !want_spline) {
    stop("spline_control has no effect because no spline will be fitted.", call. = FALSE)
  }
  continuous_blocks <- fit_control$companion_continuous_blocks
  if (!is.null(fits)) for (f in fits) {
    continuous_blocks <- union(continuous_blocks, attr(f$companions, "input_map")$continuous_blocks)
  }
  spline_basis <- NULL
  declared_groups <- if (!is.null(fits)) fits[[1L]]$sampling_schema$mixture_groups else
    if (!is.null(fit_control$mixture_groups)) fit_control$mixture_groups else spec$settings$mixture_groups
  groups_in_spec <- isTRUE(all.equal(.svem_mixture_groups(declared_groups),
    .svem_mixture_groups(spec$settings$mixture_groups), check.attributes = FALSE))
  if (want_spline && !supplied_spline) {
    spline_basis <- if (length(continuous_blocks)) {
      simpleError("The spline supports additive block IDs, not declared continuous companion block covariates.")
    } else if (!groups_in_spec) {
      simpleError("Declare mixture_groups in the shared bigexp_terms() specification for automatic spline nomination, or supply matching fitted splines.")
    } else .svem_spline_basis(spec, candidates, names(goals)[1L], spline_control,
      general = isTRUE(spline) || (!priority_missing && "spline" %in% priority) ||
        (!is.null(n_optima) && n_optima >= 5L))
    if (inherits(spline_basis, "error")) {
      explanation <- paste0("Spline candidate unavailable: ",
        sub("[[:punct:][:space:]]+$", "", conditionMessage(spline_basis)),
        ". Other models remain available. Set spline = FALSE to omit this optional nominee.")
      if (isTRUE(spline) || !priority_missing) warning(explanation, call. = FALSE)
      else message(explanation)
    }
  } else if (supplied_spline && length(continuous_blocks)) {
    stop("A spline with block IDs cannot follow a declared continuous companion block structure.", call. = FALSE)
  }
  if (is.null(fits)) {
    controls <- utils::modifyList(list(companions = c("gaupro", "extra_trees"),
                                       companion_blocks = "interacting"), fit_control)
    fits <- do.call(svem_fit_responses, c(list(spec = spec, data = data,
                                               responses = names(goals)), controls))
  }
  fits <- .svem_plan_model_list(fits, names(goals), "fits")
  if (is.null(forward)) {
    forward <- stats::setNames(lapply(names(goals), function(r) {
      do.call(svem_forward, c(list(formula = spec, data = data, response = r), forward_control))
    }), names(goals))
  }
  if (want_spline && !supplied_spline) {
    spline <- .svem_fit_spline_responses(spec, data, names(goals), candidates, spline_control,
      allow_unavailable = TRUE, basis = spline_basis,
      warn_unavailable = !inherits(spline_basis, "error"))
  } else if (!want_spline && !supplied_spline) spline <- NULL
  if (!is.null(spline)) .svem_spline_check_schema(spline, fits)
  score_args <- score_control
  if (is.null(score_args$mixture_groups)) score_args$mixture_groups <- spec$settings$mixture_groups
  scores <- do.call(svem_score_random, c(list(objects = fits, goals = goals,
    data = data, candidates = candidates, budget = budget, forward = forward, spline = spline,
    allocation = "portfolio", incumbent = incumbent, reference = reference,
    incumbent_summary = incumbent_summary, benchmark = benchmark,
    priority = priority, exploration_fraction = exploration_fraction,
    n_optima = n_optima, n_explore = n_explore), score_args))
  if (is.null(budget)) budget <- scores$allocation_counts$budget
  result <- list(fits = fits, forward = forward, scores = scores,
    worklist = scores$worklist, nominations = .svem_plan_nominations(scores),
    incumbent = scores$incumbent, benchmark = scores$benchmark, shortfall = scores$shortfall,
    policy = list(name = "default formulation workflow", budget = budget,
      specification = spec$settings, goals = goals,
      version = "1", priority = priority, incumbent = incumbent, reference = reference,
      incumbent_summary = incumbent_summary, benchmark = benchmark,
      exploration_fraction = exploration_fraction, fit_control = fit_control,
      forward_control = forward_control, score_control = score_control,
      package_version = as.character(utils::packageVersion("SVEMnet"))))
  if (!is.null(spline)) {
    result$spline <- spline
    result$policy$spline_control <- spline_control
  }
  if (!is.null(n_optima) || !is.null(n_explore)) {
    result$allocation_counts <- scores$allocation_counts
    result$merged_nominations <- scores$merged_nominations
    result$policy$n_optima <- n_optima
    result$policy$n_explore <- n_explore
  }
  structure(result, class = "svem_run_plan")
}

#' @rdname svem_plan_runs
#' @param x A `svem_run_plan`.
#' @param ... Additional arguments (reserved).
#' @export
print.svem_run_plan <- function(x, ...) {
  cat("SVEMnet default formulation workflow\n")
  if (is.null(x$allocation_counts)) {
    cat(nrow(x$worklist), "of", x$policy$budget, "preparations;",
        sum(x$worklist$new_recipe), "new recipes;", x$shortfall, "unfilled slots\n")
  } else {
    cat(nrow(x$worklist), "of", x$policy$budget, "planned preparations;",
        x$allocation_counts$unused_capacity, "unused places;", x$shortfall,
        "unmet model/exploration requests\n")
  }
  cat(sprintf("Reference repeats: %d | New model candidates: %d | Exploration runs: %d\n",
              sum(!x$worklist$new_recipe), sum(x$worklist$role == "optimum"),
              sum(x$worklist$role == "high_spread_medoid")))
  print(x$worklist[, c("slot", "recipe_id", "role", "nominated_by"), drop = FALSE], row.names = FALSE)
  if (isTRUE(x$incumbent$higher_raw_in_other_context)) {
    cat("A higher raw observed score occurs in another block context. Raw block values\n",
        "are not automatically comparable; inspect incumbent$by_context before changing the reference.\n", sep = "")
  }
  cat("Nomination status (including unavailable and budget-dropped roles):\n")
  columns <- intersect(c("engine", "status", "fit_status", "best_ties", "new_ties",
                         "new_recipe_id", "slot"), names(x$nominations))
  print(x$nominations[, columns, drop = FALSE], row.names = FALSE)
  notes <- x$scores$engine_status
  if (is.data.frame(notes) && all(c("engine", "message") %in% names(notes))) {
    notes <- notes[notes$engine == "spline" & nzchar(notes$message), , drop = FALSE]
  } else notes <- data.frame()
  if (nrow(notes)) {
    cat("Spline notes:\n")
    cat(paste0("  ", notes$response, ": ", notes$message), sep = "\n")
    cat("\n")
  }
  if (!is.null(x$allocation_counts)) {
    cat("Requested and attained preparation counts:\n")
    print(x$allocation_counts)
    if (nrow(x$merged_nominations)) {
      cat("Shared model nominees were merged; remaining models were considered for distinct recipes.\n")
      print(x$merged_nominations, row.names = FALSE)
    }
  }
  cat("Priority:", paste(x$policy$priority, collapse = " > "), "\n")
  cat("High-spread exploration fraction:", x$policy$exploration_fraction, "\n")
  if (!is.null(x$policy$incumbent_summary)) cat("Observed incumbent rule:", x$policy$incumbent_summary, "\n")
  if (!is.null(x$policy$goals)) {
    cat("Scientific goals and weights (full settings retained in policy):\n")
    print(x$policy$goals)
  }
  invisible(x)
}

#' @rdname svem_plan_runs
#' @param type `"tradeoff"` plots SVEM desirability and ensemble spread for the
#'   pool and selected runs; `"responses"` shows predicted endpoint means at
#'   the selected preparations. These are fitted predictions, not observations.
#' @export
plot.svem_run_plan <- function(x, type = c("tradeoff", "responses"), ...) {
  type <- match.arg(type)
  if (type == "responses" && !nrow(x$worklist)) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No eligible preparations: inspect shortfall and nominations") +
      ggplot2::theme_void())
  }
  if (type == "tradeoff") {
    pool <- x$scores$score_table
    work <- x$worklist
    ggplot2::ggplot(pool, ggplot2::aes(x = score, y = uncertainty_measure)) +
      ggplot2::geom_point(alpha = .2, colour = "grey50") +
      ggplot2::geom_point(data = work, ggplot2::aes(colour = role), size = 3) +
      ggplot2::geom_text(data = work, ggplot2::aes(label = slot), vjust = -1,
                        show.legend = FALSE) +
      ggplot2::labs(x = "SVEM fitted desirability", y = "SVEM ensemble spread score",
                    colour = "Preparation role") + ggplot2::theme_minimal()
  } else {
    work <- x$worklist
    rows <- lapply(names(x$fits), function(r) {
      cols <- c(SVEMnet = paste0(r, "_pred"), GauPro = paste0(r, "_gaupro_pred"),
                `Extra Trees` = paste0(r, "_extra_trees_pred"), Forward = paste0(r, "_forward_pred"), Spline = paste0(r, "_spline_pred"))
      do.call(rbind, lapply(names(cols)[cols %in% names(work)], function(e) {
        data.frame(slot = work$slot, response = rep(r, nrow(work)), engine = rep(e, nrow(work)),
                   prediction = work[[cols[[e]]]])
      }))
    })
    tab <- do.call(rbind, rows)
    ggplot2::ggplot(tab, ggplot2::aes(x = slot, y = prediction, colour = engine)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = .4), na.rm = TRUE) +
      ggplot2::facet_wrap(~response, scales = "free_y") +
      ggplot2::labs(x = "Preparation slot", y = "Fitted endpoint mean", colour = "Model") +
      ggplot2::theme_minimal()
  }
}

utils::globalVariables(c("score", "uncertainty_measure", "role", "slot", "prediction", "engine"))
