#' SVEMnet: Self-Validated Ensemble Models with Relaxed Lasso and Elastic-Net Regression
#'
#' The \code{SVEMnet} package implements Self-Validated Ensemble Models (SVEM) using Elastic Net (including lasso and ridge) regression via \code{glmnet}. SVEM averages predictions from multiple models fitted to fractionally weighted bootstraps of the data, tuned with anti-correlated validation weights.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{SVEMnet}}}{Fit an SVEMnet model using Elastic Net regression.}
#'   \item{\code{\link{predict.svem_model}}}{Predict method for SVEM models (optionally debiased, with intervals when available).}
#'   \item{\code{\link{coef.svem_model}}}{Bootstrap nonzero percentages and summary of coefficients.}
#'   \item{\code{\link{plot.svem_model}}}{Quick actual-versus-predicted plot for a fitted model.}
#'
#'   \item{\code{\link{bigexp_terms}}}{Build a deterministic expanded RHS (polynomials, interactions) with locked levels/ranges.}
#'   \item{\code{\link{bigexp_formula}}}{Reuse a locked expansion for another response to ensure identical factor space.}
#'
#'   \item{\code{\link{svem_random_table_multi}}}{Generate one shared random predictor table (with optional mixture constraints) and get predictions from multiple SVEM models at those points.}
#'   \item{\code{\link{svem_optimize_random}}}{Random-search optimizer for multiple responses with goals, weights, optional CIs, and diverse PAM-medoids candidates.}
#'
#'   \item{\code{\link{svem_significance_test}}}{Whole-model significance test with optional mixture-group sampling.}
#'   \item{\code{\link{svem_significance_test_parallel}}}{Parallel whole-model significance test (foreach + doParallel).}
#'
#'   \item{\code{\link{plot_svem_significance_tests}}}{Plot helper for visualizing multiple significance-test outputs.}
#'   \item{\code{\link{glmnet_with_cv}}}{Convenience wrapper around repeated \code{cv.glmnet()} selection.}
#'
#'   \item{\code{\link{lipid_screen}}}{Example dataset for multi-response modeling and mixture-constrained optimization demos.}
#' }

#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @keywords package
#' @template ref-svem
#'
#' @docType package
#' @name SVEMnet-package
"_PACKAGE"
