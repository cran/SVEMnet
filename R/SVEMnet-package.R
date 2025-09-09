#' SVEMnet: Self-Validated Ensemble Models with Relaxed Lasso and Elastic-Net Regression
#'
#' The \code{SVEMnet} package implements Self-Validated Ensemble Models (SVEM) using Elastic Net (including lasso and ridge) regression via \code{glmnet}. SVEM averages predictions from multiple models fitted to fractionally weighted bootstraps of the data, tuned with anti-correlated validation weights.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{SVEMnet}}}{Fit an SVEMnet model using Elastic Net regression.}
#'   \item{\code{\link{predict.svem_model}}}{Predict method for SVEM models.}
#'   \item{\code{\link{svem_significance_test}}}{Whole-model significance test.}
#'   \item{\code{\link{svem_significance_test_parallel}}}{Parallel whole-model test.}
#'   \item{\code{\link{plot.svem_model}}}{Plot actual vs predicted for a model.}
#'   \item{\code{\link{coef.svem_model}}}{Coefficient nonzero percentages (bootstrap).}
#'   \item{\code{\link{plot_svem_significance_tests}}}{Plotting helper for multiple tests.}
#'   \item{\code{\link{glmnet_with_cv}}}{Wrapper around \code{cv.glmnet()} with repeats.}
#' }
#'
#' @section Acknowledgments:
#' Development of this package was assisted by GPT o1-preview, which helped in constructing the structure of some of the code and the roxygen documentation. The code for the significance test is taken from the supplementary material of Karl (2024) (it was handwritten by that author).
#'
#' @keywords package
#' @references
#' Gotwalt, C., & Ramsey, P. (2018). Model Validation Strategies for Designed Experiments Using Bootstrapping Techniques With Applications to Biopharmaceuticals. \emph{JMP Discovery Conference}. \url{https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849873/redirect_from_archived_page/true}
#'
#' Karl, A. T. (2024). A randomized permutation whole-model test heuristic for Self-Validated Ensemble Models (SVEM). \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{249}, 105122. \doi{10.1016/j.chemolab.2024.105122}
#'
#' Karl, A., Wisnowski, J., & Rushing, H. (2022). JMP Pro 17 Remedies for Practical Struggles with Mixture Experiments. JMP Discovery Conference. \doi{10.13140/RG.2.2.34598.40003/1}
#'
#' Lemkus, T., Gotwalt, C., Ramsey, P., & Weese, M. L. (2021). Self-Validated Ensemble Models for Design of Experiments. \emph{Chemometrics and Intelligent Laboratory Systems}, 219, 104439. \doi{10.1016/j.chemolab.2021.104439}
#'
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020). Applications of the Fractional-Random-Weight Bootstrap. \emph{The American Statistician}, 74(4), 345–358. \doi{10.1080/00031305.2020.1731599}
#'
#' Ramsey, P., Gaudard, M., & Levin, W. (2021). Accelerating Innovation with Space Filling Mixture Designs, Neural Networks and SVEM. \emph{JMP Discovery Conference}. \url{ https://community.jmp.com/t5/Abstracts/Accelerating-Innovation-with-Space-Filling-Mixture-Designs/ev-p/756841}
#'
#' Ramsey, P., & Gotwalt, C. (2018). Model Validation Strategies for Designed Experiments Using Bootstrapping Techniques With Applications to Biopharmaceuticals. \emph{JMP Discovery Conference - Europe}. \url{https://community.jmp.com/t5/Abstracts/Model-Validation-Strategies-for-Designed-Experiments-Using/ev-p/849647/redirect_from_archived_page/true}
#'
#' Ramsey, P., Levin, W., Lemkus, T., & Gotwalt, C. (2021). SVEM: A Paradigm Shift in Design and Analysis of Experiments. \emph{JMP Discovery Conference - Europe}. \url{https://community.jmp.com/t5/Abstracts/SVEM-A-Paradigm-Shift-in-Design-and-Analysis-of-Experiments-2021/ev-p/756634}
#'
#' Ramsey, P., & McNeill, P. (2023). CMC, SVEM, Neural Networks, DOE, and Complexity: It's All About Prediction. \emph{JMP Discovery Conference}.
#'
#' @docType package
#' @name SVEMnet-package
"_PACKAGE"
