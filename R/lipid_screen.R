#' Lipid formulation screening data
#'
#' An example dataset for modeling Potency, Size, and PDI as functions of
#' formulation and process settings. Percent composition columns are stored as
#' proportions in \code{[0, 1]} (e.g., 4.19\% is 0.0419). This table is intended
#' for demonstration of SVEMnet multi-response modeling and desirability-based
#' random-search optimization.
#'
#' @name lipid_screen
#' @docType data
#' @usage data(lipid_screen)
#' @keywords datasets
#'
#' @format A data frame with N rows and the following columns:
#' \describe{
#'   \item{RunID}{character. Optional identifier.}
#'   \item{PEG}{numeric. Proportion (0–1).}
#'   \item{Helper}{numeric. Proportion (0–1).}
#'   \item{Ionizable}{numeric. Proportion (0–1).}
#'   \item{Cholesterol}{numeric. Proportion (0–1).}
#'   \item{Ionizable_Lipid_Type}{factor.}
#'   \item{N_P_ratio}{numeric.}
#'   \item{flow_rate}{numeric.}
#'   \item{Potency}{numeric. Response.}
#'   \item{Size}{numeric. Response (e.g., nm).}
#'   \item{PDI}{numeric. Response (polydispersity index).}
#'   \item{Notes}{character. Optional free-text notes.}
#' }
#'
#' @details
#' This dataset accompanies examples showing:
#' \itemize{
#'   \item fitting three SVEM models (Potency, Size, PDI) on a shared expanded
#'         factor space via \code{bigexp_terms()} and \code{bigexp_formula()},
#'   \item random design generation using SVEM random-table helpers (for use
#'         with multi-response optimization),
#'   \item multi-response optimization with \code{svem_optimize_random()} using
#'         Derringer–Suich desirabilities and weighted combining
#'         (\code{combine = "geom"} or \code{combine = "mean"}),
#'   \item returning both high-score \emph{optimal} candidates and
#'         high-uncertainty \emph{exploration} candidates,
#'   \item optional whole-model reweighting (WMT) of response weights via
#'         p-values, and exporting per-row scores to the original data with
#'         \code{original_data_scored}.
#' }
#'
#' @source Simulated screening table supplied by the package author.
#' @template ref-svem
#' @examples
#' \donttest{
#' # 1) Load the bundled dataset
#' data(lipid_screen)
#' str(lipid_screen)
#'
#' # 2) Build a deterministic expansion using bigexp_terms()
#' #    Provide main effects only on the RHS; expansion width is controlled via arguments.
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,   # up to 3-way interactions
#'   polynomial_order = 3,   # include up to cubic terms I(X^2), I(X^3),
#'   include_pc_2way  = TRUE, # include partial cubic terms
#'   include_pc_3way  = FALSE
#' )
#'
#' # 3) Reuse the same locked expansion for other responses
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' # 4) Fit SVEM models with the shared factor space/expansion
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
#' # Mixture constraints: components sum to 1, with bounds
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' # 7) Run random-search optimization (DS + optimal & exploration candidates)
#' set.seed(2)
#' opt_out <- svem_optimize_random(
#'   objects                  = objs,
#'   goals                    = goals,
#'   n                        = 25000,
#'   mixture_groups           = mix,
#'   level                    = 0.95,
#'   k_candidates             = 5,
#'   top_frac                 = 0.02,
#'   k_exploration_candidates = 5,
#'   exploration_top_frac     = 0.05,
#'   numeric_sampler          = "random",
#'   verbose                  = TRUE
#' )
#'
#' # Inspect optimal solution and candidates (scores and uncertainty included)
#' opt_out$best
#' opt_out$best_pred
#' opt_out$best_ci
#' head(opt_out$candidates)
#'
#' # Inspect exploration target and candidates
#' opt_out$exploration_best
#' head(opt_out$exploration_candidates)
#'
#' # 8) Repeat with WMT reweighting using the original data (requires 'data')
#' #    Choose either "neglog10" (aggressive) or "one_minus_p" (conservative).
#' set.seed(3)
#' opt_wmt <- svem_optimize_random(
#'   objects                  = objs,
#'   goals                    = goals,
#'   data                     = lipid_screen,  # used for WMT and original_data_scored
#'   n                        = 25000,
#'   mixture_groups           = mix,
#'   level                    = 0.95,
#'   k_candidates             = 5,
#'   top_frac                 = 0.02,
#'   k_exploration_candidates = 5,
#'   exploration_top_frac     = 0.05,
#'   numeric_sampler          = "random",
#'   reweight_by_wmt          = TRUE,
#'   wmt_transform            = "neglog10",
#'   verbose                  = TRUE
#' )
#'
#' # Compare weights and look at candidates under WMT
#' opt_wmt$weights_original
#' opt_wmt$weights_final
#' opt_wmt$wmt_p_values
#' head(opt_wmt$candidates)
#' head(opt_wmt$exploration_candidates)
#'
#' # Scored original data (predictions, desirabilities, score, uncertainty)
#' head(opt_wmt$original_data_scored)
#' }
"lipid_screen"
