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
#' A typical workflow is to fit models with \code{SVEMnet()}, optionally run a
#' whole-model test (for example, for response reweighting), and then call
#' \code{svem_optimize_random()} to propose optimal and exploration candidates
#' for the next experimental round.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{SVEMnet}}}{Fit an SVEMnet model using Elastic Net regression.}
#'   \item{\code{\link{predict.svem_model}}}{Predict method for SVEM models (optionally debiased, with intervals when available).}
#'   \item{\code{\link{coef.svem_model}}}{Averaged (optionally debiased) coefficients from an SVEM model.}
#'   \item{\code{\link{svem_nonzero}}}{Bootstrap nonzero percentages for each coefficient, with an optional quick plot.}
#'   \item{\code{\link{plot.svem_model}}}{Quick actual-versus-predicted plot for a fitted model.}
#'
#'   \item{\code{\link{bigexp_terms}}}{Build a deterministic expanded RHS (polynomials, interactions) with locked levels/ranges.}
#'   \item{\code{\link{bigexp_formula}}}{Reuse a locked expansion for another response to ensure identical factor space.}
#'
#'   \item{\code{\link{svem_random_table_multi}}}{Generate one shared random predictor table (with optional mixture constraints) and get predictions from multiple SVEM models at those points.}
#'   \item{\code{\link{svem_optimize_random}}}{Random-search optimizer for multiple responses with goals, weights, optional CIs, and diverse PAM-medoids candidates.}
#'
#'   \item{\code{\link{svem_significance_test_parallel}}}{Parallel whole-model significance test (foreach + doParallel) with optional mixture-group sampling.}
#'   \item{\code{\link[=plot.svem_significance_test]{plot.svem_significance_test}}}{Plot helper for visualizing multiple significance-test outputs.}
#'
#'   \item{\code{\link{glmnet_with_cv}}}{Convenience wrapper around repeated \code{cv.glmnet()} selection.}
#'
#'   \item{\code{\link{lipid_screen}}}{Example dataset for multi-response modeling and mixture-constrained optimization demos.}
#' }
#'
#' @section Families:
#' Supports Gaussian and binomial responses (logit link). For \code{family = "binomial"},
#' the response must be 0/1 numeric or a two-level factor (first level treated as 0).
#' Use \code{predict(..., type = "response")} for event probabilities or
#' \code{type = "class"} for 0/1 labels (threshold = 0.5 by default).
#'
#' @keywords package
#' @template ref-svem
#'
#' @docType package
#' @name SVEMnet-package
"_PACKAGE"
