#' Lipid formulation screening data
#'
#' An example dataset for modeling Potency, Size, and PDI as functions of
#' formulation and process settings. Percent composition columns are stored as
#' proportions in \code{[0, 1]} (for example, 4.19 percent is 0.0419). This
#' table is intended for demonstration of SVEMnet multi-response modeling,
#' desirability-based random-search optimization, and probabilistic
#' design-space construction.
#'
#' @name lipid_screen
#' @docType data
#' @usage data(lipid_screen)
#' @keywords datasets
#'
#' @format A data frame with one row per experimental run and the following columns:
#' \describe{
#'   \item{RunID}{character. Optional identifier for each run.}
#'   \item{PEG}{numeric. Proportion (0-1).}
#'   \item{Helper}{numeric. Proportion (0-1).}
#'   \item{Ionizable}{numeric. Proportion (0-1).}
#'   \item{Cholesterol}{numeric. Proportion (0-1).}
#'   \item{Ionizable_Lipid_Type}{factor. Categorical identity of the ionizable lipid.}
#'   \item{N_P_ratio}{numeric. Molar or mass \eqn{N:P} ratio (unitless).}
#'   \item{flow_rate}{numeric. Process flow rate (arbitrary units).}
#'   \item{Operator}{factor. Categorical blocking factor.}
#'   \item{Potency}{numeric. Response (for example, normalized activity).}
#'   \item{Size}{numeric. Response (for example, particle size in nm).}
#'   \item{PDI}{numeric. Response (polydispersity index).}
#'   \item{Notes}{character. Optional free-text notes.}
#' }
#'
#' @details
#' The four composition columns \code{PEG}, \code{Helper}, \code{Ionizable},
#' and \code{Cholesterol} are stored as proportions in \code{[0,1]}, and in many
#' rows they sum (approximately) to 1, making them natural candidates for
#' mixture constraints in optimization and design-space examples.
#'
#' This dataset accompanies examples showing:
#' \itemize{
#'   \item fitting three SVEM models (Potency, Size, PDI) on a shared expanded
#'         factor space via \code{\link{bigexp_terms}} and \code{\link{bigexp_formula}},
#'   \item random design generation using SVEM random-table helpers (for use
#'         with multi-response optimization),
#'   \item multi-response scoring and candidate selection with
#'         \code{\link{svem_score_random}} (Derringer–Suich desirabilities,
#'         weights, uncertainty) and \code{\link{svem_select_from_score_table}}
#'         (optimal and high-uncertainty medoid candidates),
#'   \item returning both high-score optimal candidates and
#'         high-uncertainty exploration candidates from the same scored table
#'         by changing the target column (for example \code{score} vs
#'         \code{uncertainty_measure} or \code{wmt_score}),
#'   \item optional whole-model reweighting (WMT) of response weights via
#'         \code{\link{svem_wmt_multi}} (for p-values and multipliers) together
#'         with \code{\link{svem_score_random}} (via its \code{wmt} argument),
#'   \item constructing a probabilistic design space in one step by passing
#'         process-mean specifications via the \code{specs} argument of
#'         \code{\link{svem_score_random}} (internally using
#'         \code{svem_append_design_space_cols} when needed).
#' }
#'
#' @source Simulated screening table supplied by the package author.
#' @template ref-svem
#'
#' @examples
#' \donttest{
#' # 1) Load the bundled dataset
#' data(lipid_screen)
#' str(lipid_screen)
#'
#' #2) Build a deterministic expansion using bigexp_terms()
#' #  Provide main effects only on the right-hand side; expansion width
#' #  is controlled via arguments. Here Operator is treated as a blocking
#' #  factor: additive only, no interactions or polynomial terms.
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,   # up to 3-way interactions
#'   polynomial_order = 3,   # include up to cubic terms I(X^2), I(X^3)
#'   include_pc_2way  = TRUE,  # partial-cubic two-way terms Z:I(X^2)
#'   include_pc_3way  = FALSE,  # no partial-cubic three-way terms I(X^2):Z:W
#'    blocking         = "Operator"
#' )
#'
#' # 3) Reuse the same locked expansion for other responses
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' # 4) Fit SVEM models with the shared factor space and expansion
#' set.seed(1)
#' fit_pot <- SVEMnet(form_pot, lipid_screen)
#' fit_siz <- SVEMnet(form_siz, lipid_screen)
#' fit_pdi <- SVEMnet(form_pdi, lipid_screen)
#'
#' # 5) Collect models in a named list by response
#' objs <- list(Potency = fit_pot, Size = fit_siz, PDI = fit_pdi)
#'
#' # 6) Define multi-response goals and weights (DS desirabilities under the hood)
#' #    Maximize Potency (0.6), minimize Size (0.3), minimize PDI (0.1)
#' goals <- list(
#'   Potency = list(goal = "max", weight = 0.6),
#'   Size    = list(goal = "min", weight = 0.3),
#'   PDI     = list(goal = "min", weight = 0.1)
#' )
#'
#' # Mixture constraints: components sum to total, with bounds
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' # 7) Optional: whole-model tests and WMT multipliers via svem_wmt_multi()
#' #    This wrapper runs svem_significance_test_parallel() for each response,
#' #    plots the distance distributions, and prints p-values and multipliers.
#'
#' set.seed(123)
#' wmt_out <- svem_wmt_multi(
#'   formulas       = list(Potency = form_pot,
#'                         Size    = form_siz,
#'                         PDI     = form_pdi),
#'   data           = lipid_screen,
#'   mixture_groups = mix,
#'   wmt_control    = list(seed = 123),
#'   plot           = TRUE
#' )
#'
#' # Inspect WMT p-values and multipliers (also printed by the wrapper)
#' wmt_out$p_values
#' wmt_out$multipliers
#'
#' # 8) Optional: define process-mean specifications for a joint design space.
#' #    Potency at least 78, Size no more than 100, PDI less than 0.25.
#' #    Here we only specify the bounded side; the unbounded side defaults to
#' #    lower = -Inf or upper = Inf inside svem_score_random().
#' specs_ds <- list(
#'   Potency = list(lower = 78),
#'   Size    = list(upper = 100),
#'   PDI     = list(upper = 0.25)
#' )
#'
#' # 9) Random-search scoring in one step via svem_score_random()
#' #    This draws a random candidate table, computes DS desirabilities,
#' #    a combined multi-response score, WMT-adjusted wmt_score (if `wmt`
#' #    is supplied), CI-based uncertainty, and (when `specs` is supplied)
#' #    appends mean-level design-space columns.
#' #
#' #    The `wmt` and `specs` arguments are optional:
#' #      - Omit `wmt` for no whole-model reweighting.
#' #      - Omit `specs` if you do not need design-space probabilities.
#'
#' set.seed(3)
#' scored <- svem_score_random(
#'   objects         = objs,
#'   goals           = goals,
#'   data            = lipid_screen,  # scored and returned as original_data_scored
#'   n               = 25000,
#'   mixture_groups  = mix,
#'   level           = 0.95,
#'   combine         = "geom",
#'   numeric_sampler = "random",
#'   wmt             = wmt_out,   # optional: NULL for no WMT
#'   specs           = specs_ds,  # optional: NULL for no design-space columns
#'   verbose         = TRUE
#' )
#'
#' # 10) Select optimal and exploration sets from the same scored table
#'
#' # Optimal medoid candidates (maximizing DS score)
#' opt_sel <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "score",      # score column is maximized
#'   direction   = "max",
#'   k           = 5,
#'   top_type    = "frac",
#'   top         = 0.1,
#'   label       = "round1_score_optimal"
#' )
#'
#' # Optimal medoid candidates (maximizing WMT-adjusted wmt_score)
#' opt_sel_wmt <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "wmt_score",  # wmt_score column is maximized
#'   direction   = "max",
#'   k           = 5,
#'   top_type    = "frac",
#'   top         = 0.1,
#'   label       = "round1_wmt_optimal"
#' )
#'
#' # Exploration medoid candidates (highest uncertainty_measure)
#' explore_sel <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "uncertainty_measure",  # uncertainty_measure column is maximized
#'   direction   = "max",
#'   k           = 5,
#'   top_type    = "frac",
#'   top         = 0.1,
#'   label       = "round1_explore"
#' )
#'
#' # In-spec medoid candidates (highest joint mean-level assurance)
#' inspec_sel <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "p_joint_mean",         # p_joint_mean column is maximized
#'   direction   = "max",
#'   k           = 5,
#'   top_type    = "frac",
#'   top         = 0.10,
#'   label       = "round1_inspec"
#' )
#'
#' # Single best by score, including per-response CIs
#' opt_sel$best
#'
#' # Single best by WMT-adjusted score, including per-response CIs
#' opt_sel_wmt$best
#'
#' # Diverse high-score candidates (medoids)
#' head(opt_sel_wmt$candidates)
#'
#' # Highest-uncertainty setting and its medoid candidates
#' explore_sel$best
#' head(explore_sel$candidates)
#'
#' # Highest probability mean-in-spec setting and its medoid candidates
#' inspec_sel$best
#' head(inspec_sel$candidates)
#'
#' # 11) Scored original data (predictions, desirabilities, score, wmt_score, uncertainty)
#' head(scored$original_data_scored)
#'
#' # 12) Example: combine new candidate settings with the best existing run
#' #     and (optionally) export a CSV for the next experimental round.
#'
#' # Best existing run from the original scored data (no new medoids; k = 0)
#' best_existing <- svem_select_from_score_table(
#'   score_table = scored$original_data_scored,
#'   target      = "score",
#'   direction   = "max",
#'   k           = 0,               # k <= 0 => only the best row, no medoids
#'   top_type    = "frac",
#'   top         = 1.0,
#'   label       = "round1_existing_best"
#' )
#'
#' # 13) Prepare candidate export tables for the next experimental round.
#' #     svem_export_candidates_csv() accepts individual selection objects.
#' #     Calls are commented out so examples/tests do not write files.
#'
#' # Export a design-style candidate table (factors + responses + predictions)
#' # out.df <- svem_export_candidates_csv(
#' #   opt_sel,
#' #   opt_sel_wmt,
#' #   explore_sel,
#' #   inspec_sel,
#' #   best_existing,
#' #   file      = "lipid_screen_round1_candidates.csv",
#' #   overwrite = TRUE
#' # )
#' # head(out.df)
#'
#' # Export all columns including desirabilities, CI widths, and design-space columns
#' # out.df2 <- svem_export_candidates_csv(
#' #   opt_sel,
#' #   opt_sel_wmt,
#' #   explore_sel,
#' #   inspec_sel,
#' #   best_existing,
#' #   file      = "lipid_screen_round1_candidates_all.csv",
#' #   overwrite = TRUE
#' # )
#' # head(out.df2)
#' }
"lipid_screen"
