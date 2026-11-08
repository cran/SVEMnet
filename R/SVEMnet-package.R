#' SVEMnet: Self-Validated Ensemble Models with Relaxed Lasso and Elastic-Net Regression
#'
#' The \code{SVEMnet} package implements Self-Validated Ensemble Models (SVEM)
#' using Elastic Net (including lasso and ridge) regression via \code{glmnet}.
#' SVEM averages predictions from multiple models fitted to fractionally
#' weighted bootstraps of the data, tuned with anti-correlated validation
#' weights. The package supports multi-response optimization with
#' uncertainty-aware candidate generation for iterative formulation and
#' process development.
#'
#' The primary package reference is Karl (2026),
#' \doi{10.1016/j.chemolab.2026.105660}, published in
#' \emph{Chemometrics and Intelligent Laboratory Systems}. It describes the
#' elastic-net and relaxed elastic-net ensemble workflow, shared factor
#' expansions, and multi-response optimization that form the package's core.
#' The analyses in that paper used SVEMnet version 3.2.0; see the package
#' NEWS for subsequent changes.
#' Use \code{citation("SVEMnet")} for the package paper and related references.
#' The LNP workflow paper by Karl et al. (2023) is available in JoVE
#' (\doi{10.3791/65200}) and on arXiv (\url{https://arxiv.org/abs/2212.11264}).
#' The SVEMnet package paper is also on arXiv
#' (\url{https://arxiv.org/abs/2511.20968}). The journal articles are fixed
#' versions of record. The authors update the arXiv manuscripts occasionally.
#' Each arXiv link opens the latest posted revision.
#'
#' A Gaussian response fit can also request GauPro and ranger Extra Trees
#' companions. The worked example builds a shared feasible recipe table,
#' compares the models' nominees, and updates their fits after simulated
#' preparations. Run \code{vignette("lipid-screen-follow-up", package = "SVEMnet")}.
#' Astra assisted with development of the version 4.0.0 companion fitting and
#' prediction code in \code{SVEMnet()} and \code{predict()}, the candidate
#' worklist in \code{svem_score_random()}, and \code{lipid_screen_simulate()}.
#' Astra also assisted with implementation of the default workflow and the
#' optional Gaussian spline. The cited literature provides the statistical
#' foundations of these methods.
#'
#' @section Default formulation workflow:
#' Start with \code{vignette("lipid-screen-follow-up", package = "SVEMnet")}.
#' The same interface handles process-only factors, mixtures and mixture-process
#' factors, including ingredient categories, through a shared specification and
#' an explicit feasible recipe pool. The default portfolio requires numeric
#' Gaussian responses on the declared analysis scale. Binomial models remain
#' available through the individual fitting and scoring functions.
#' \enumerate{
#'   \item Declare factors with \code{\link{bigexp_terms}}, scientific response
#'         goals, mixture constraints, a feasible pool and the planned block.
#'         The existing cubic default includes pure powers, two-factor partial
#'         cubics and three-way interactions. It is a rich candidate basis for
#'         model reduction, not a claim that every term is independently estimable.
#'   \item Call \code{\link{svem_plan_runs}} with the total preparation budget.
#'         The default carries the best observed reference forward, then proposes
#'         new optima from elastic-net SVEM, GauPro, Extra Trees and forward SVEM.
#'         Shared optima use one preparation. Remaining slots are uncertainty medoids.
#'   \item Inspect fits, scores, nominations, reference provenance and any
#'         shortfall. Plot the plan and explicitly export its worklist with
#'         \code{\link{svem_export_candidates_csv}}.
#'   \item Measure every selected preparation, append the outcomes, and call the
#'         same planner for the next round. Use the same scientific response scale
#'         and anchors when comparing rounds.
#' }
#' Whole-model tests and specification-probability diagnostics are optional;
#' neither is required by this streamlined workflow. Users can instead compose
#' the individual functions, including \code{\link{svem_fit_responses}} and
#' \code{\link{svem_score_random}}, or customize the named controls.
#'
#' @section Planning the initial experiment:
#' Karl's \emph{A Sphericity Filter for Space-Filling Designs in Nondeterministic
#' Experiments} (\url{https://arxiv.org/abs/2610.36012}) is a useful starting
#' point for selecting an initial space-filling experiment with protection for
#' low-order effects. Its working model for design protection is separate from
#' the richer candidate model used for response analysis. SVEMnet currently
#' consumes measured data and feasible recipes; it does not implement that
#' initial-design filter. Declare the physical region, prediction target and
#' allowed categorical combinations before generating an initial experiment.
#'
#' @section Core modeling and summaries:
#' \describe{
#'   \item{\code{\link{SVEMnet}}}{Fit an SVEMnet model using Elastic Net regression
#'         (including relaxed elastic net) on fractionally weighted bootstraps.}
#'   \item{\code{\link{predict.svem_model}}}{Predict method for SVEM models
#'         (ensemble-mean aggregation by default, optional debiasing, and
#'         percentile intervals across bootstrap member predictions).}
#'   \item{\code{\link{coef.svem_model}}}{Averaged (optionally debiased)
#'         coefficients from an SVEM model.}
#'   \item{\code{\link{svem_nonzero}}}{Bootstrap nonzero percentages for each
#'         coefficient, with an optional quick plot.}
#'   \item{\code{\link{plot.svem_model}}}{Quick actual-versus-predicted plot for a
#'         fitted model (with ROC, precision-recall, and calibration displays
#'         for binomial fits).}
#' }
#'
#' @section Deterministic wide expansions (bigexp helpers):
#' The \code{bigexp_*} helpers build and reuse a locked polynomial/interaction
#' expansion across multiple responses and datasets:
#' \describe{
#'   \item{\code{\link{bigexp_terms}}}{Build a deterministic expanded RHS
#'         (polynomials, interactions, optional partial-cubic terms) with
#'         locked factor levels and numeric ranges.}
#'   \item{\code{\link{bigexp_prepare}}}{Coerce new data to match a stored
#'         \code{bigexp_spec}, including factor levels and numeric types.}
#'   \item{\code{\link{bigexp_formula}}}{Reuse a locked expansion for another
#'         response to ensure an identical factor space across models.}
#'   \item{\code{\link{with_bigexp_contrasts}}}{Temporarily restore the
#'         contrast options used when a \code{bigexp_spec} was built.}
#'   \item{\code{\link{bigexp_train}}}{Convenience wrapper that builds a
#'         \code{bigexp_spec} and prepares training data in one call.}
#' }
#'
#' @section Random tables, optimization, and candidate generation:
#' \describe{
#'   \item{\code{\link{svem_random_table_multi}}}{Generate one shared random
#'         predictor table (with optional mixture constraints) from cached
#'         factor-space information and obtain predictions from multiple SVEM
#'         models at those points. Supports both Gaussian and binomial models;
#'         binomial predictions are returned on the probability scale. This is
#'         the lower-level sampler used by \code{\link{svem_score_random}}.}
#'   \item{\code{\link{svem_score_random}}}{Random-search scoring for multiple
#'         responses with Derringer–Suich desirabilities, user weights,
#'         optional whole-model-test (WMT) reweighting, percentile CI-based
#'         uncertainty, and (optionally) scoring of existing experimental data.
#'         Returns a scored random-search table and, when \code{data} is
#'         supplied, an augmented copy of the original data with
#'         \code{<resp>_pred}, desirabilities, scores, and an
#'         \code{uncertainty_measure}.}
#'   \item{\code{\link{svem_select_from_score_table}}}{Given a scored table
#'         (typically \code{svem_score_random()$score_table}), select one
#'         "best" row under a chosen objective and a small, diverse set of
#'         medoid candidates via PAM clustering on predictors.}
#'   \item{\code{\link{svem_export_candidates_csv}}}{Concatenate one or more
#'         selection objects from \code{\link{svem_select_from_score_table}}
#'         and export candidate tables (with metadata, predictions, and
#'         response columns) to CSV or return them in-memory for
#'         inspection.}
#' }
#'
#' @section Whole-model testing and plotting:
#' \describe{
#'   \item{\code{\link{svem_significance_test_parallel}}}{Whole-model
#'         significance test with serial execution by default and an optional
#'         private PSOCK cluster for explicit multi-core runs. Supports
#'         mixture-constrained sampling and reuse of a locked
#'         \code{bigexp_spec}; designed for continuous (Gaussian) responses.}
#'   \item{\code{\link{svem_wmt_multi}}}{Helper to run
#'         \code{svem_significance_test_parallel} across multiple responses and
#'         construct whole-model p-values and reweighting multipliers for use
#'         in \code{\link{svem_score_random}}.}
#'   \item{\code{\link[=plot.svem_significance_test]{plot.svem_significance_test}}}{Plot helper
#'         for visualizing one or more significance-test outputs (boxplots of
#'         observed vs permutation distances).}
#' }
#'
#' @section Auxiliary utilities and data:
#' \describe{
#'   \item{\code{\link{glmnet_with_cv}}}{Convenience wrapper around repeated
#'         \code{cv.glmnet()} selection for robust lambda (and optional alpha)
#'         choice.}
#'   \item{\code{\link{lipid_screen}}}{Example dataset for multi-response
#'         modeling, whole-model testing, and mixture-constrained optimization
#'         demonstrations.}
#' }
#'
#' @section Extensions beyond the package paper:
#' GauPro and ranger Extra Trees can nominate additional recipes from the same
#' physical candidate pool. \code{\link{svem_score_random}} can allocate a
#' finite preparation budget to these nominees and to exploration recipes.
#' The package also includes newer methods for further study:
#' \code{\link{svem_forward}} and its single-model benchmark
#' \code{\link{forward_aicc}}, the bootstrap-member batch proposer
#' \code{\link{svem_thompson_batch}}, and the FRW variance diagnostic
#' \code{\link{svem_ij_variance}}. The optional Gaussian thin-plate smoother,
#' \code{\link{spline_terms}} and \code{\link{svem_spline}}, uses the suggested
#' \pkg{mgcv} package. It can nominate an additional recipe on a supported
#' numeric domain after the four default model roles. The opt-in
#' \code{complexity = "edf"} selector in \code{\link{SVEMnet}} is also
#' experimental; the published package workflow uses
#' \code{complexity = "support"}. These extensions are not established by
#' the results in the package paper.
#' Opt-in polynomial centering in \code{\link{bigexp_terms}} is a later
#' design-matrix construction option and was not evaluated in that paper.
#'
#' Bootstrap member-percentile intervals describe variation among member
#' predictions. They do not include new-observation noise and do not by
#' themselves provide calibrated confidence or prediction intervals.
#'
#' @section Families:
#' SVEMnet currently supports:
#' \itemize{
#'   \item Gaussian responses (\code{family = "gaussian"}) with identity link
#'         and optional debiasing / bootstrap member-percentile intervals.
#'   \item Binomial responses (\code{family = "binomial"}) with logit link.
#'         The response must be 0/1 numeric or a two-level factor (first level
#'         treated as 0). Use \code{predict(..., type = "response")} for event
#'         probabilities or \code{type = "class"} for 0/1 labels
#'         (threshold = 0.5 by default).
#' }
#' Some higher-level utilities place additional constraints:
#' \itemize{
#'   \item \code{\link{svem_significance_test_parallel}} is designed and
#'         interpreted for continuous (Gaussian) responses.
#'   \item \code{\link{svem_score_random}} supports mixed Gaussian + binomial
#'         response sets, treating binomial predictions and CIs on the
#'         probability scale, but WMT-based goal reweighting (via
#'         \code{\link{svem_wmt_multi}} and the \code{wmt} argument) is only
#'         allowed when all responses are Gaussian.
#' }
#'
#' @examples
#' if (requireNamespace("GauPro", quietly = TRUE) &&
#'     requireNamespace("ranger", quietly = TRUE)) {
#'   set.seed(41)
#'   design <- data.frame(x = seq(-1, 1, length.out = 16))
#'   design$y <- 2 - (design$x - 0.3)^2 + rnorm(nrow(design), sd = 0.05)
#'   # A small bootstrap count keeps this example quick.
#'   fit <- SVEMnet(y ~ x, design, nBoot = 4L, companions = TRUE)
#'   predict(fit, data.frame(x = c(-0.5, 0, 0.5)), engine = "all")
#' }
#' @keywords package
#' @template ref-svem
#'
#' @docType package
#' @name SVEMnet-package
"_PACKAGE"
