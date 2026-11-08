#' Fit an SVEMnet model (Self-Validated Ensemble Elastic Net)
#'
#' Fit a Self-Validated Ensemble Model (SVEM) with elastic net or relaxed
#' elastic net base learners using \code{glmnet}. Fractional random-weight
#' (FRW) train/validation weights are drawn on each bootstrap replicate,
#' a validation-weighted information criterion (wAIC, wBIC, or wSSE) is
#' minimized to select the penalty, and predictions are ensembled across
#' replicates. Gaussian and binomial responses are supported.
#'
#' This is the core ensemble fitter described in Karl (2026),
#' \doi{10.1016/j.chemolab.2026.105660}. The published package workflow uses
#' the default \code{complexity = "support"}; the optional
#' \code{complexity = "edf"} selector is an experimental extension.
#'
#' You can pass either:
#' \itemize{
#'   \item a standard model formula, e.g.,
#'     \code{y ~ X1 + X2 + X3 + I(X1^2) + (X1 + X2 + X3)^2}, or
#'   \item a \code{bigexp_spec} created by \code{bigexp_terms()}, in which
#'     case SVEMnet will build the design matrix deterministically (locked
#'     types, levels, and contrasts) and, if requested, swap the response
#'     to fit multiple independent responses over the same expansion.
#' }
#'
#' @param formula A formula specifying the model to be fitted, or a
#'   \code{bigexp_spec} created by \code{bigexp_terms()}.
#' @param data A data frame containing the variables in the model.
#' @param nBoot Integer. Number of bootstrap replicates (default \code{200}).
#'   Each replicate draws FRW weights and fits a glmnet path.
#' @param glmnet_alpha Numeric vector of elastic net mixing parameters
#'   \code{alpha} in \code{[0, 1]}. \code{alpha = 1} is lasso, \code{alpha = 0}
#'   is ridge. Defaults to \code{c(0.5, 1)}. When \code{relaxed = TRUE},
#'   \code{alpha = 0} is automatically dropped (ridge + relaxed is not used).
#' @param weight_scheme Character. Weighting scheme for train/validation
#'   copies. One of:
#'   \itemize{
#'     \item \code{"SVEM"} (default): Self-Validated Ensemble Model weights.
#'       For each replicate and row, a shared uniform draw \code{U_i ~ Unif(0, 1)}
#'       is converted to anti-correlated FRW weights
#'       \code{w_train_i = -log(U_i)} and \code{w_valid_i = -log(1 - U_i)}.
#'       Each weight vector is then rescaled to have mean 1 (sum \code{n}).
#'     \item \code{"FRW_plain"}: Fractional random-weight regression without a
#'       separate validation copy. A single FRW vector
#'       \code{w_i = -log(U_i)} is used for both training and validation and
#'       rescaled to have mean 1 (sum \code{n}). This reproduces the FRW
#'       bootstrap regression of Xu et al. (2020) and related work.
#'     \item \code{"Identity"}: Uses unit weights for both training and
#'       validation (no resampling). In combination with \code{nBoot = 1} this
#'       wraps a single \code{glmnet} fit and selects the penalty by the chosen
#'       information criterion, while still using SVEMnet's expansion and
#'       diagnostics.
#'   }
#' @param objective Character. One of \code{"auto"}, \code{"wAIC"},
#'   \code{"wBIC"}, or \code{"wSSE"}.
#'   See Details for the automatic choices, the information criteria and the
#'   loss-only criterion for Gaussian and binomial responses.
#' @param relaxed Logical or character. Default \code{"auto"}.
#'   If \code{TRUE}, use glmnet's relaxed elastic-net path and select both
#'   the penalty \code{lambda} and the relaxed refit parameter \code{gamma}
#'   on each bootstrap. If \code{FALSE}, fit the standard glmnet path without
#'   the relaxed step. If \code{"auto"} (default), SVEMnet uses
#'   \code{relaxed = TRUE} for \code{family = "gaussian"} and
#'   \code{relaxed = FALSE} for \code{family = "binomial"}.
#' @param response Optional character. When \code{formula} is a \code{bigexp_spec},
#'   this names the response column to use on the left-hand side. Defaults to
#'   the response stored in the spec.
#' @param unseen How to treat factor levels not seen in the original
#'   \code{bigexp_spec} when \code{formula} is a \code{bigexp_spec}. One of
#'   \code{"warn_na"} (default; convert unseen levels to \code{NA} with a
#'   warning) or \code{"error"} (stop with an error).
#' @param family Character. One of \code{"gaussian"} (default) or
#'   \code{"binomial"}. For Gaussian models SVEMnet uses the identity link; for
#'   binomial it uses the canonical logit link. The binomial response must be
#'   numeric 0/1, logical, or a factor with exactly two levels (the second
#'   level is treated as 1).
#' @param complexity Character. Complexity measure used in the wAIC/wBIC
#'   penalty term. One of:
#'   \itemize{
#'     \item \code{"support"} (default): the number of nonzero coefficients
#'       (including the intercept), as described in the package references.
#'       It is a model-size proxy. For a fixed lasso tuning parameter it
#'       estimates degrees of freedom under the usual full-rank assumptions;
#'       collinearity and adaptive tuning require additional care.
#'     \item \code{"edf"} (experimental, provided for further study): a
#'       gamma-blended effective degrees of freedom. For an elastic-net path
#'       point with active set \eqn{A}, penalty \eqn{\lambda}, and mixing
#'       parameter \eqn{\alpha < 1}, the penalized fit's effective df is the
#'       ridge-type trace
#'       \eqn{\sum_j e_j / (e_j + n(\lambda/s_{y,w})(1-\alpha))}, where \eqn{e_j} are
#'       the eigenvalues of the weighted, standardized Gram matrix of the
#'       active columns, and \eqn{s_{y,w}} is the training-weighted response
#'       standard deviation. This conditional trace holds the active set and
#'       tuning choices fixed; it does not include their adaptation cost.
#'       Nonuniform penalty factors enter the corresponding diagonal ridge
#'       matrix. Finite coefficient bounds are not supported for EDF selection;
#'       support-selected bounded fits report EDF as unavailable.
#'       The relaxed refit's df is the rank of the active
#'       design. The blended complexity for relaxed mixing parameter
#'       \eqn{\gamma} is
#'       \eqn{1 + \gamma \, df_{enet} + (1-\gamma) \, rank(A)}.
#'       For \eqn{\alpha = 1} (lasso) this reduces to one plus the rank of
#'       the centered active design (Tibshirani and Taylor, 2012), equal to
#'       support when those columns are linearly independent. The numerical
#'       rank uses Gram eigenvalues exceeding \eqn{10^{-10}} times the largest
#'       eigenvalue; very ill-conditioned columns can therefore count as
#'       dependent. Support also overstates the conditional trace for
#'       \eqn{\alpha < 1}. Implemented for \code{family = "gaussian"}
#'       only; binomial fits fall back to \code{"support"} with a warning.
#'   }
#'   The admissibility guardrail (the nonzero non-intercept count must be
#'   strictly less than \eqn{n_{eff,adm}}) and
#'   the reported model-size diagnostics always use the support count,
#'   regardless of this setting. \code{"edf"} changes only the complexity
#'   value inserted into the wAIC/wBIC penalty and is intended for simulation
#'   study of the alpha/gamma comparison; the published SVEM results use
#'   \code{"support"}.
#'   This separate conservative constraint can exclude elastic-net models
#'   with many nonzeros even when their conditional EDF is small.
#' @param store_member_weights Logical; default \code{FALSE}. If \code{TRUE},
#'   retain the exact per-member uniforms and the mean-one training and
#'   validation weights in \code{member_weights}. This can materially increase
#'   object size. The stored weights support reproducibility and experimental
#'   ensemble-mean diagnostics; they do not make the member-percentile limits
#'   calibrated confidence or prediction intervals.
#' @param ... Additional arguments passed to \code{glmnet()}, such as
#'   \code{penalty.factor}, \code{lower.limits}, \code{upper.limits}, or
#'   \code{standardize.response}. An \code{exclude} callback receives the
#'   original design columns and the current fit's response and weights;
#'   an internal compatibility column is never passed to it. Numerically
#'   constant columns remain excluded in addition to user exclusions. When
#'   applying this mask or handling a one-column design, the package checks
#'   that callback returns index the original columns (or form an aligned
#'   logical vector); other callbacks retain \code{glmnet}'s validation.
#'   Any user-supplied \code{weights} are ignored
#'   (SVEM supplies its own fractional random weights), any user-supplied
#'   \code{offset} is ignored (SVEMnet does not support offsets), and any
#'   user-supplied \code{standardize} is ignored; SVEMnet always calls
#'   \code{glmnet} with \code{standardize = TRUE}. Glmnet algorithm-control
#'   values (for example \code{maxit} or \code{thresh}) supplied directly or
#'   through \code{control = list(...)} are routed in a version-compatible way
#'   when possible. The arguments \code{nlambda}, \code{lambda}, and
#'   \code{lambda.min.ratio} are reserved and managed by SVEMnet
#'   (\code{nlambda} is fixed at 500); user-supplied values are ignored with
#'   a warning. Argument names that the installed \code{glmnet} does not
#'   recognize are also ignored with a warning (misspelling protection).
#'   Relaxed fits additionally accept `maxp` and `path`.
#'
#' @details The Gaussian SVEMnet elastic-net search requests 500 lambda values with
#'   glmnet's deviance-based early path stopping disabled (`fdev = 0`,
#'   `devmax = 1`). Other solver limits can still shorten the requested path.
#'   Disabling these rules preserves solver convergence checks, coefficient limits and
#'   SVEM's validation-size support guard. Explicit algorithm controls are
#'   honored. These call-level defaults override session-level stopping settings.
#'   Use `control = list(fdev = 1e-5, devmax = 0.999)` to recover the
#'   former stopping policy. Modern glmnet uses per-call controls; older
#'   supported versions use a scoped override restored even on error.
#'   This default applies to Gaussian SVEMnet members, including fits inside
#'   the whole-model test. Binomial members and the separate `glmnet_with_cv()`
#'   comparison wrapper retain their former stopping defaults. To request the
#'   extended binomial path, supply `control = list(fdev = 0, devmax = 1)`. Retained
#'   solver warnings, their member count, the glmnet version and the effective
#'   call-level controls are available in `diagnostics`. Changing the path can
#'   change fits at any noise level and whole-model p-values or multipliers.
#'   Results are not uniformly better on new data: a small stress check saw
#'   some small-sample and wSSE fits worsen.
#'
#'   GauPro companions normalize the response internally and return means and
#'   latent standard errors in the original response units. An unsupported
#'   single-factor input or a constant response is reported before claiming a
#'   usable GP fit. Engine messages are retained with companion warnings.
#'
#' @param companions Optional Gaussian-response companion fits. \code{FALSE}
#'   fits SVEM alone. \code{TRUE} requests both \pkg{GauPro} and \pkg{ranger}
#'   Extra Trees. A character vector can select \code{"gaupro"} or
#'   \code{"extra_trees"} individually. Companions use the original factors
#'   on the SVEM analysis rows. Requested fits and diagnostics are stored in
#'   \code{$companions}; an unavailable engine leaves the primary fit usable.
#'   The companion component retains the unique physical predictor rows used
#'   for that endpoint to identify previously prepared recipes at nomination.
#'   The Gaussian process uses GauPro's Matern 5/2 kernel on continuous
#'   coordinates and latent-factor kernels on nominal inputs. It estimates
#'   the nugget with one seeded optimization start. The fit records its seed,
#'   package version, kernel settings and attained GauPro deviance on its
#'   internally normalized response scale. That deviance is not comparable
#'   with earlier unnormalized fits. Printing the companion component reports
#'   how many fit notes were retained in its `warnings` field.
#' @param companion_blocks With \code{"respect"} (default), a varying block
#'   declared in a BigExp specification makes the companion fits unsupported.
#'   \code{"interacting"} explicitly permits each companion to learn a
#'   different recipe surface per block. Numeric day, plate or batch IDs
#'   are nominal labels by default. Supply observed IDs when scoring them.
#'   An unreplicated numeric block triggers a diagnostic at fitting time.
#' @param companion_continuous_blocks Names of measured numeric blocking
#'   covariates, such as ambient temperature, to treat as continuous inputs in
#'   an interacting companion. Such a covariate can be scored between observed
#'   values within its fitted range. Other blocks retain their nominal encoding.
#'   The SVEM formula and its block penalties are unaffected.
#' @param mixture_groups Optional list of fixed-sum mixture groups in
#'   \code{\link{svem_random_table_multi}} format. A BigExp fit inherits its
#'   groups from \code{\link{bigexp_terms}}. Companions use one fewer
#'   varying coordinate per group. They drop the last varying component in
#'   the declared order and check all components at prediction. The GP and forest
#'   can give different fits in a different choice of free coordinates.
#' @return An object of class \code{"svem_model"} (and \code{"svem_binomial"}
#'   when \code{family = "binomial"}) with components:
#' \itemize{
#'   \item \code{parms}: Vector of ensemble-averaged coefficients, including
#'     the intercept.
#'   \item \code{parms_debiased}: Vector of coefficients after optional
#'     debiasing (see Details; Gaussian only).
#'   \item \code{debias_fit}: If debiasing was performed, the calibration
#'     model \code{lm(y ~ y_pred)}; otherwise \code{NULL}.
#'   \item \code{coef_matrix}: Matrix of per-bootstrap coefficients
#'     (rows = bootstraps, columns = intercept and predictors).
#'   \item \code{nBoot}: Number of bootstrap replicates actually used.
#'   \item \code{glmnet_alpha}: Vector of alpha values considered.
#'   \item \code{best_alphas}: Per-bootstrap alpha selected by the criterion.
#'   \item \code{best_lambdas}: Per-bootstrap lambda selected by the criterion.
#'   \item \code{best_relax_gammas}: Per-bootstrap relaxed gamma selected when
#'     \code{relaxed = TRUE}; \code{NA} otherwise.
#'   \item \code{complexity}: The complexity measure actually used
#'     (\code{"support"} or \code{"edf"}).
#'   \item \code{best_edfs}: Per-bootstrap blended effective df of the
#'     selected Gaussian path point. This diagnostic is computed even when
#'     \code{complexity = "support"}; it is \code{NA} for binomial fits.
#'   \item \code{member_diagnostics}: Per-member Gaussian validation SSE,
#'     selected objective value, support size, effective df, raw and bounded
#'     Kish validation size, support-based residual df and scale, fallback
#'     status, and admissibility warnings. Gaussian residual-scale fields are
#'     \code{NA} for binomial fits.
#'   \item \code{member_weights}: If \code{store_member_weights = TRUE}, a
#'     list containing matrices \code{uniforms}, \code{training}, and
#'     \code{validation}, with rows aligned to \code{coef_matrix}; otherwise
#'     \code{NULL}.
#'   \item \code{weight_scheme}: The weighting scheme that was used.
#'   \item \code{relaxed}: Logical flag indicating whether relaxed paths were
#'     used.
#'   \item \code{relaxed_input}: The user-supplied value for \code{relaxed}
#'     (one of \code{TRUE}, \code{FALSE}, or \code{"auto"}). The resolved flag
#'     actually used is reported in \code{relaxed}.
#'   \item \code{dropped_alpha0_for_relaxed}: Logical; \code{TRUE} if
#'     \code{alpha = 0} was dropped because \code{relaxed = TRUE}.
#'   \item \code{objective_input}: The objective requested by the user.
#'   \item \code{objective_used}: The objective actually used after applying
#'     the "auto" rule (for example \code{"wAIC"} or \code{"wBIC"}).
#'   \item \code{objective}: Same as \code{objective_used} (for convenience).
#'   \item \code{auto_used}: Logical; \code{TRUE} if \code{objective = "auto"}.
#'   \item \code{auto_decision}: The objective selected by the auto rule
#'     (wAIC or wBIC) when \code{auto_used = TRUE}.
#'   \item \code{diagnostics}: List with summary information, including:
#'     \itemize{
#'       \item \code{k_summary}: Median and IQR of selected model size
#'         (number of nonzero coefficients including intercept).
#'       \item \code{fallback_rate}: Proportion of bootstraps that fell back
#'         to an intercept-only fit.
#'       \item \code{fallback_reasons}: Named counts summarizing why any
#'         non-structural bootstrap fits required that fallback.
#'        \item \code{n_eff_summary}: Summary of raw Kish effective validation
#'         sizes \eqn{n_eff = (\sum_i w_i^{valid})^2 / \sum_i (w_i^{valid})^2}
#'         across bootstraps (before truncation to form \eqn{n_eff_adm}).
#'       \item \code{alpha_freq}: Relative frequency of selected alpha values
#'         (if any).
#'       \item \code{relax_gamma_freq}: Relative frequency of selected relaxed
#'         gamma values (if \code{relaxed = TRUE} and any were selected).
#'     }
#'   \item \code{actual_y}: Numeric response vector used in fitting (0/1 for
#'     binomial).
#'   \item \code{training_X}: Numeric model matrix without the intercept
#'     column used for training.
#'   \item \code{y_pred}: Fitted values from the ensemble on the training
#'     data. For Gaussian this is on the response scale; for binomial it is
#'     on the probability scale.
#'   \item \code{y_pred_debiased}: Debiased fitted values on the training
#'     data (Gaussian only); \code{NULL} otherwise.
#'   \item \code{nobs}: Number of observations used in fitting.
#'   \item \code{nparm}: Number of parameters in the full expansion
#'     (intercept plus predictors).
#'   \item \code{formula}: The formula used for fitting (possibly derived
#'     from a \code{bigexp_spec}).
#'   \item \code{terms}: \code{terms} object used for building the design
#'     matrix, with a compact evaluation environment containing any required
#'     formula transformations.
#'   \item \code{xlevels}: Factor levels recorded at training time.
#'   \item \code{contrasts}: Contrasts used for building the design matrix.
#'   \item \code{schema}: Compact description for safe prediction, including
#'     \code{feature_names}, \code{terms_str}, \code{xlevels},
#'     \code{contrasts}, \code{contrasts_options}, and a simple hash.
#'   \item \code{sampling_schema}: Schema used to generate random candidate
#'     tables, including predictor names, variable classes, numeric ranges,
#'     and factor levels.
#'   \item \code{used_bigexp_spec}: Logical flag indicating whether a
#'     \code{bigexp_spec} was used.
#'   \item \code{family}: The fitted family ("gaussian" or "binomial").
#' }
#'
#' @details
#' SVEMnet implements Self-Validated Ensemble Models using elastic net and
#' relaxed elastic net base learners from glmnet. Each bootstrap replicate
#' draws fractional random weights, builds a train and validation copy, fits
#' a path over lambda (and optionally over alpha and relaxed gamma), and
#' selects a path point by minimizing a validation-weighted criterion. Final
#' predictions are obtained by averaging replicate predictions on the chosen
#' scale.
#'
#' By default, \code{relaxed = "auto"} resolves to \code{TRUE} for Gaussian
#' fits and \code{FALSE} for binomial fits.
#'
#' The function is typically used in small-n design-of-experiments (DOE)
#' workflows where classical train/validation splits and cross-validation
#' can be unstable. A common pattern is: (1) build a deterministic expansion
#' with \code{bigexp_terms()}, (2) fit SVEM models via \code{SVEMnet()},
#' (3) perform whole-model significance testing, and (4) call
#' \code{svem_score_random()} for constrained multi-response optimization.
#'
#' Weighting schemes:
#' \itemize{
#'   \item With \code{weight_scheme = "SVEM"}, SVEMnet uses a pair of
#'     anti-correlated FRW vectors for train and validation. All rows appear
#'     in every replicate, but train and validation contributions are
#'     separated through the shared uniform draws.
#'   \item With \code{weight_scheme = "FRW_plain"}, a single FRW vector is
#'     used for both train and validation, which reproduces FRW regression
#'     without a self-validation split. This is mainly provided for method
#'     comparison and teaching.
#'   \item With \code{weight_scheme = "Identity"}, both train and validation
#'     weights are 1. Setting \code{nBoot = 1} in this mode yields a single
#'     glmnet fit whose penalty is chosen by the selected information
#'     criterion, without any bootstrap variation.
#' }
#'
#' Selection criteria (Gaussian):
#' For \code{family = "gaussian"}, the validation loss is a weighted sum of
#' squared errors on the validation copy. Let
#' \eqn{\mathrm{SSE}_w = \sum_i w_i^{valid} r_i^2} denote the weighted SSE.
#' The criteria are:
#' \itemize{
#'   \item \code{"wSSE"}: loss-only selector that minimizes the weighted SSE
#'     \eqn{\mathrm{SSE}_w},
#'   \item \code{"wAIC"}: Gaussian AIC analog
#'     \eqn{C(\lambda) = n \log\{\mathrm{SSE}_w(\lambda)/n\} + 2 k},
#'   \item \code{"wBIC"}: Gaussian BIC analog
#'     \eqn{C(\lambda) = n \log\{\mathrm{SSE}_w(\lambda)/n\} +
#'       \log(n_eff_adm)\, k}.
#' }
#' The FRW validation weights are rescaled to have mean one, so that their sum
#' is always \eqn{\sum_i w_i^{valid} = n}. The AIC/BIC analogs therefore use
#' \eqn{n \log(\mathrm{SSE}_w/n)} as the Gaussian loss term, while the
#' \code{"wSSE"} selector uses \eqn{\mathrm{SSE}_w} directly.
#'
#' The effective validation size is computed from the FRW weights using Kish's
#' effective sample size \eqn{n_eff = (\sum_i w_i^{valid})^2 /
#'   \sum_i (w_i^{valid})^2} and then truncated to lie between 2 and \code{n}
#' to form \eqn{n_eff_adm}. The AIC-style selector uses a \eqn{2 k} penalty;
#' the BIC-style selector uses a \eqn{\log(n_eff_adm) k} penalty, so that the
#' loss term is scaled by total validation weight while the complexity penalty
#' uses a weight-concentration adjustment. For \code{"wAIC"} and
#' \code{"wBIC"}, the number of non-intercept coefficients must be strictly
#' less than \eqn{n_eff_adm}; equality is inadmissible.
#'
#' Because the FRW validation weights are random rather than fixed design
#' weights, these information-criterion scores are used heuristically for
#' relative model comparison within each FRW replicate, rather than as exact
#' AIC/BIC values.
#' The training and validation copies reuse the same responses; they are not
#' independent held-out samples. Kish size measures weight concentration for
#' a weighted mean, not regression rank, leverage or selection complexity.
#' The name \code{wAIC} does not refer to Bayesian WAIC. Also,
#' \code{wBIC} has a smaller per-coefficient penalty than \code{wAIC} when
#' \eqn{n_eff_adm < \exp(2)}. Normalizing weights to total \code{n} defines
#' the loss-to-penalty scale; substituting Kish size in the loss would be a
#' different selection rule.
#'
#' For diagnostics, SVEMnet reports the raw Kish effective sizes across
#' bootstraps (see \code{diagnostics$n_eff_summary}), while \eqn{n_eff_adm}
#' is used internally in the penalty and model-size guardrail. Near-interpolating
#' path points are screened out via a simple model size guardrail before
#' minimization. When \code{objective = "auto"}, Gaussian fits use
#' \code{"wAIC"} and binomial fits use \code{"wBIC"}.
#' This guard is a heuristic cap, not a guarantee of positive residual
#' degrees of freedom. With \code{weight_scheme = "Identity"}, support may
#' reach \code{n}; shrunken fits can still have positive SSE, while exact
#' interpolation makes the Gaussian profile likelihood singular.
#' Nonpositive support-based residual denominators are flagged in
#' \code{member_diagnostics}.
#'
#' This structure (pseudo-likelihood using total weight and BIC penalty using a
#' Kish-type effective sample size) parallels survey-weighted information
#' criteria as in Lumley and Scott (2015) and Kish (1965).
#'
#' Selection criteria (binomial):
#' For \code{family = "binomial"}, the validation loss is the weighted
#' negative log-likelihood on the FRW validation copy (equivalently,
#' proportional to the binomial deviance up to a constant factor). Let
#' \eqn{\mathrm{NLL}} denote the weighted negative log-likelihood. The same
#' labels are used:
#' \itemize{
#'   \item \code{"wSSE"}: loss-only selector based on \eqn{\mathrm{NLL}}
#'     (the name is retained for backward compatibility),
#'   \item \code{"wAIC"}: deviance-style criterion
#'     \eqn{C(\lambda) = 2\,\mathrm{NLL}(\lambda) + 2 k},
#'   \item \code{"wBIC"}: deviance-style criterion
#'     \eqn{C(\lambda) = 2\,\mathrm{NLL}(\lambda) + \log(n_eff_adm)\, k}.
#' }
#' The effective validation size \eqn{n_eff_adm} and the model size guardrail
#' are handled as in the Gaussian case: for \code{"wAIC"} and \code{"wBIC"}
#' we compute a Kish effective size from the FRW validation weights, truncate
#' it to lie between 2 and \code{n}, and require the number of nonzero
#' coefficients (excluding the intercept) to be less than this effective size
#' when evaluating the criterion.
#'
#' Auto rule:
#' When \code{objective = "auto"}, SVEMnet selects the criterion by family:
#' \itemize{
#'   \item \code{family = "gaussian"} -> \code{"wAIC"}
#'   \item \code{family = "binomial"} -> \code{"wBIC"}
#' }
#'
#' Relaxed elastic net:
#' When \code{relaxed = TRUE}, SVEMnet calls glmnet with
#' \code{relax = TRUE} and traverses a small grid of relaxed refit values
#' (gamma). For each alpha and gamma, SVEMnet evaluates all lambda path
#' points on the validation copy and records the combination that minimizes
#' the selected criterion. Model size is always defined as the number of
#' nonzero coefficients including the intercept, so standard and relaxed
#' paths are scored on the same scale.
#'
#' Gaussian debiasing:
#' For Gaussian models, SVEMnet optionally performs a simple linear
#' calibration of ensemble predictions on the training data. When there is
#' sufficient variation in the fitted values and \code{nBoot} is at least
#' 10, the function fits \code{lm(y ~ y_pred)} and uses the coefficients to
#' construct debiased coefficients and debiased fitted values. Binomial
#' fits do not use debiasing; predictions are ensembled on the probability
#' or link scale directly.
#'
#' Implementation notes:
#' \itemize{
#'   \item Predictors are always standardized internally via
#'     \code{glmnet(..., standardize = TRUE)}.
#'   \item The terms object is stored with a compact evaluation environment so
#'     that prediction supports standard and local transformation functions
#'     without accidentally capturing the full calling environment.
#'   \item A compact schema (feature names, factor levels, contrasts, and
#'     a simple hash) is stored to allow \code{predict()} and companion
#'     functions to rebuild model matrices deterministically, even when the
#'     original data frame is not available.
#'   \item A separate sampling schema stores raw predictor ranges and
#'     factor levels for use in random candidate generation for optimization.
#' }
#'
#' @template ref-svem
#'
#' @references
#' Zou, H., Hastie, T., & Tibshirani, R. (2007). On the "degrees of freedom"
#' of the lasso. \emph{The Annals of Statistics}, 35(5), 2173–2192.
#' \doi{10.1214/009053607000000127}
#'
#' @importFrom stats runif lm predict coef var model.frame model.matrix model.response delete.response IQR median
#' @importFrom glmnet glmnet
#'
#' @examples
#' set.seed(42)
#'
#' n  <- 30
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X3 <- rnorm(n)
#' eps <- rnorm(n, sd = 0.5)
#' y <- 1 + 2 * X1 - 1.5 * X2 + 0.5 * X3 + 1.2 * (X1 * X2) +
#'   0.8 * (X1^2) + eps
#' dat <- data.frame(y, X1, X2, X3)
#'
#' # Minimal hand-written expansion
#' mod_relax <- SVEMnet(
#'   y ~ (X1 + X2 + X3)^2 + I(X1^2) + I(X2^2),
#'   data          = dat,
#'   glmnet_alpha  = 1,
#'   nBoot         = 10,
#'   objective     = "auto",
#'   weight_scheme = "SVEM",
#'   relaxed       = FALSE
#' )
#'
#' pred_in_raw <- predict(mod_relax, dat, debias = FALSE)
#' pred_in_db  <- predict(mod_relax, dat, debias = TRUE)
#'
#' # ---------------------------------------------------------------------------
#' # Big expansion (full factorial + polynomial surface + partial-cubic crosses)
#' # Build once, reuse for one or more responses
#' # ---------------------------------------------------------------------------
#' spec <- bigexp_terms(
#'   y ~ X1 + X2 + X3,
#'   data             = dat,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_3way  = FALSE
#' )
#'
#' # Fit using the spec (auto-prepares data)
#' fit_y <- SVEMnet(
#'   spec, dat,
#'   glmnet_alpha  = 1,
#'   nBoot         = 8,
#'   objective     = "auto",
#'   weight_scheme = "SVEM"
#' )
#'
#' # A second, independent response over the same expansion
#' set.seed(99)
#' dat$y2 <- 0.5 + 1.4 * X1 - 0.6 * X2 + 0.2 * X3 + rnorm(n, 0, 0.4)
#' fit_y2 <- SVEMnet(
#'   spec, dat, response = "y2",
#'   glmnet_alpha  = 1,
#'   nBoot         = 8,
#'   objective     = "auto",
#'   weight_scheme = "SVEM"
#' )
#'
#' svem_nonzero(fit_y2)
#'
#' p1 <- predict(fit_y,  dat)
#' p2 <- predict(fit_y2, dat, debias = TRUE)
#'
#' # Show that a new batch expands identically under the same spec
#' newdat <- data.frame(
#'   y  = y,
#'   X1 = X1 + rnorm(n, 0, 0.05),
#'   X2 = X2 + rnorm(n, 0, 0.05),
#'   X3 = X3 + rnorm(n, 0, 0.05)
#' )
#' prep_new <- bigexp_prepare(spec, newdat)
#' stopifnot(identical(
#'   colnames(model.matrix(spec$formula, bigexp_prepare(spec, dat)$data)),
#'   colnames(model.matrix(spec$formula, prep_new$data))
#' ))
#' preds_new <- predict(fit_y, prep_new$data)
#'
#' ## Binomial example
#' set.seed(2)
#' n  <- 80
#' X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
#' eta <- -0.3 + 1.1 * X1 - 0.8 * X2 + 0.5 * X1 * X3
#' p   <- plogis(eta)
#' yb  <- rbinom(n, 1, p)
#' db  <- data.frame(yb = yb, X1 = X1, X2 = X2, X3 = X3)
#'
#' fit_b <- SVEMnet(
#'   yb ~ (X1 + X2 + X3)^2, db,
#'   nBoot        = 10,
#'   glmnet_alpha = 1,
#'   family       = "binomial"
#' )
#'
#' ## Probabilities, link, and classes
#' p_resp <- predict(fit_b, db, type = "response")
#' p_link <- predict(fit_b, db, type = "link")
#' y_hat  <- predict(fit_b, db, type = "class")  # 0/1 labels
#'
#' ## Mean aggregation with uncertainty on probability scale
#' out_b <- predict(
#'   fit_b, db,
#'   type     = "response",
#'   se.fit   = TRUE,
#'   interval = TRUE,
#'   level    = 0.9
#' )
#' str(out_b)
#'
#' ## Example with blocking (requires SVEMnet to store sampling_schema$blocking)
#' set.seed(2)
#' df_block <- data.frame(
#'   y1         = rnorm(40),
#'   y2         = rnorm(40),
#'   X1         = runif(40),
#'   X2         = runif(40),
#'   Operator   = factor(sample(paste0("Op", 1:3), 40, TRUE)),
#'   AmbientTmp = rnorm(40, mean = 22, sd = 2)
#' )
#'
#' spec_block <- bigexp_terms(
#'   y1 ~ X1 + X2,
#'   data             = df_block,
#'   factorial_order  = 2,
#'   polynomial_order = 2,
#'   blocking         = c("Operator", "AmbientTmp")
#' )
#'
#' fit_b1 <- SVEMnet(spec_block, df_block, response = "y1", nBoot = 8,
#'                   glmnet_alpha = 1, relaxed = FALSE)
#' fit_b2 <- SVEMnet(spec_block, df_block, response = "y2", nBoot = 8,
#'                   glmnet_alpha = 1, relaxed = FALSE)
#'
#' tab_block <- svem_random_table_multi(list(fit_b1, fit_b2), n = 200)
#' @export
SVEMnet <- function(formula, data,
                    nBoot = 200,
                    glmnet_alpha = c( 0.5, 1),
                    weight_scheme = c("SVEM", "FRW_plain", "Identity"),
                    objective = c("auto", "wAIC", "wBIC", "wSSE"),
                    relaxed = "auto",
                    response = NULL,
                    unseen = c("warn_na","error"),
                    family = c("gaussian", "binomial"),
                    complexity = c("support", "edf"),
                    ...,
                    store_member_weights = FALSE,
                    companions = FALSE,
                    companion_blocks = c("respect", "interacting"),
                    mixture_groups = NULL,
                    companion_continuous_blocks = character()) {

  engines <- .svem_companion_names(companions)
  companion_blocks <- match.arg(companion_blocks)
  if (is.null(companion_continuous_blocks)) companion_continuous_blocks <- character()
  if (!is.character(companion_continuous_blocks) || anyNA(companion_continuous_blocks) ||
      any(!nzchar(companion_continuous_blocks)) || anyDuplicated(companion_continuous_blocks)) {
    stop("companion_continuous_blocks must be a character vector of unique block names.",
         call. = FALSE)
  }
  mixture_groups <- .svem_mixture_groups(mixture_groups)

  ## ----------------- family handling -----------------
  if (inherits(family, "family")) {
    stop(
      "SVEMnet currently accepts family only as a character string: ",
      "family = 'gaussian' or family = 'binomial'."
    )
  }

  family   <- match.arg(family)
  fam_name <- family

  if (!fam_name %in% c("gaussian", "binomial")) {
    stop("SVEMnet currently supports family = 'gaussian' or 'binomial' only.")
  }


  ## ---- argument checks ----
  objective     <- match.arg(objective)
  weight_scheme <- match.arg(weight_scheme)
  unseen        <- match.arg(unseen)
  complexity    <- match.arg(complexity)
  store_member_weights <- .svem_logical_scalar(
    store_member_weights, "store_member_weights"
  )

  if (identical(complexity, "edf") && fam_name == "binomial") {
    warning("complexity = 'edf' is implemented for family = 'gaussian' only; ",
            "using complexity = 'support'.")
    complexity <- "support"
  }
  use_edf <- identical(complexity, "edf")

  ## ---- relaxed handling (supports logical or "auto") ----
  relaxed_input <- relaxed
  relax_flag <- NA

  if (is.character(relaxed)) {
    relaxed <- match.arg(relaxed, c("auto"))
    # family has already been parsed into fam_name above
    relax_flag <- if (relaxed == "auto") (fam_name == "gaussian") else NA
  } else if (is.logical(relaxed) && length(relaxed) == 1L && !is.na(relaxed)) {
    relax_flag <- relaxed
  } else {
    stop("'relaxed' must be TRUE, FALSE, or the character string \"auto\" (default).")
  }


  nBoot <- .svem_integer_scalar(nBoot, "nBoot", min = 1L)
  nBoot_input <- nBoot  # keep what user asked for (fix #2)

  if (!is.numeric(glmnet_alpha) || any(!is.finite(glmnet_alpha))) {
    stop("'glmnet_alpha' must be numeric and finite.")
  }
  if (any(glmnet_alpha < 0 | glmnet_alpha > 1)) {
    stop("'glmnet_alpha' values must be in [0, 1].")
  }
  glmnet_alpha <- unique(as.numeric(glmnet_alpha))
  if (!length(glmnet_alpha)) glmnet_alpha <- 1



  ## ---- model frame / matrix construction ----
  data <- as.data.frame(data)
  contrasts_opts_used <- NULL

  # Detect bigexp_spec either passed directly or attached to the formula
  using_spec     <- inherits(formula, "bigexp_spec")
  spec           <- NULL
  spec_from_attr <- FALSE

  if (using_spec) {
    # Case 1: SVEMnet(spec, data, ...)
    spec <- formula
  } else {
    # Case 2: SVEMnet(form_with_attr, data, ...)
    spec_attr <- attr(formula, "bigexp_spec", exact = TRUE)
    if (!is.null(spec_attr) && inherits(spec_attr, "bigexp_spec")) {
      using_spec     <- TRUE
      spec           <- spec_attr
      spec_from_attr <- TRUE
    }
  }

  if (!is.null(spec$settings$mixture_groups)) {
    inherited <- spec$settings$mixture_groups
    if (!is.null(mixture_groups) && !identical(mixture_groups, inherited))
      stop("mixture_groups conflicts with the groups declared in bigexp_terms().", call. = FALSE)
    mixture_groups <- inherited
  }
  if (using_spec && length(mixture_groups) &&
      (length(spec$settings$mixture_vars) || isTRUE(spec$settings$center_polynomials)) &&
      !setequal(spec$settings$mixture_vars, unlist(lapply(mixture_groups, `[[`, "vars"), use.names = FALSE))) {
    stop("Declare consistent mixture_groups and mixture membership in bigexp_terms().",
         call. = FALSE)
  }
  contrasts_used_full <- NULL  # will capture pre-intercept-drop contrasts (fix #1)
  schema_data <- data

  if (using_spec) {
    # Decide which formula to use:
    # - If 'response' is supplied, always swap LHS to that name.
    # - Else, if the spec came via attribute, respect the formula's own LHS.
    # - Else (direct spec), use spec$formula as stored.
    if (!is.null(response)) {
      f_use <- stats::as.formula(paste(response, "~", spec$rhs), env = baseenv())
    } else if (spec_from_attr) {
      f_use <- formula
    } else {
      f_use <- spec$formula
    }

    prep <- bigexp_prepare(spec, data, unseen = unseen)
    schema_data <- prep$data

    spec_settings <- spec$settings
    if (is.null(spec_settings)) spec_settings <- list()

    # Blocking variables from the bigexp_spec (may be NULL / empty)
    blocking <- spec_settings$blocking
    if (is.null(blocking)) blocking <- character(0L)

    # Contrast options + per-factor contrasts from the spec
    if (!is.null(spec_settings$contrasts_options)) {
      spec_contrasts_opts <- spec_settings$contrasts_options
    } else {
      spec_contrasts_opts <- getOption("contrasts")
    }

    if (!is.null(spec_settings$contrasts)) {
      spec_contrasts <- spec_settings$contrasts
    } else {
      spec_contrasts <- NULL
    }

    contrasts_opts_used <- spec_contrasts_opts
    old_opts <- options(contrasts = spec_contrasts_opts)
    on.exit(options(old_opts), add = TRUE)

    mf <- stats::model.frame(f_use, prep$data, na.action = stats::na.omit)
    if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")
    X  <- stats::model.matrix(f_use, mf, contrasts.arg = spec_contrasts)

    # capture contrasts BEFORE any subsetting that might drop attributes (fix #1)
    contrasts_used_full <- attr(X, "contrasts")

  } else {
    # Plain formula path, no bigexp_spec
    f_use <- formula
    mf    <- stats::model.frame(f_use, data, na.action = stats::na.omit)
    mf <- .svem_drop_unused_frame_levels(mf)
    if (nrow(mf) < 2L) stop("Not enough complete cases after NA removal.")
    contrasts_opts_used <- getOption("contrasts")
    X       <- stats::model.matrix(f_use, mf)
    blocking <- character(0L)   # no blocking when not using bigexp_spec

    # capture contrasts BEFORE any subsetting that might drop attributes (fix #1)
    contrasts_used_full <- attr(X, "contrasts")
  }

  y <- stats::model.response(mf)

  terms_full <- attr(mf, "terms")
  .svem_reject_formula_offset(terms_full, "SVEMnet")
  model_frame_schema <- .svem_model_frame_schema(mf)
  if (!identical(attr(terms_full, "intercept"), 1L) ||
      (using_spec && identical(spec$settings$intercept, FALSE))) {
    stop(
      "SVEMnet fits follow the published intercept convention; formulas using `~ 0`, ",
      "`- 1`, or `bigexp_terms(intercept = FALSE)` are not supported.",
      call. = FALSE
    )
  }

  if (is.matrix(y) || is.array(y) || length(y) != nrow(mf)) {
    stop("SVEMnet requires a single univariate response; matrix/multivariate responses are not supported.",
         call. = FALSE)
  }

  ## ---- response handling by family ----
  .coerce_binomial_01 <- function(y) {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) {
        stop("Binomial SVEMnet requires a factor with exactly 2 levels or 0/1 numeric y.")
      }
      return(as.integer(y == levels(y)[2L]))
    } else if (is.logical(y)) {
      return(as.integer(y))
    } else if (is.numeric(y)) {
      uy <- sort(unique(y))
      if (!all(uy %in% c(0,1))) {
        stop("Binomial SVEMnet requires numeric y coded as 0/1.")
      }
      return(as.integer(y))
    } else {
      stop("Unsupported y type for binomial; use 0/1 numeric, logical, or a 2-level factor.")
    }
  }

  if (fam_name == "gaussian") {
    if (!is.numeric(y)) stop("Gaussian SVEMnet requires numeric y.")
    y_vec <- as.numeric(y)
  } else if (fam_name == "binomial") {
    y_vec <- .coerce_binomial_01(y)
  }

  ## drop intercept column (glmnet adds its own); match only the exact
  ## model.matrix name so a real predictor named "Intercept" is not deleted
  int_idx <- which(colnames(X) == "(Intercept)")
  if (length(int_idx)) X <- X[, -int_idx, drop = FALSE]

  # re-attach contrasts explicitly in case subsetting dropped them (fix #1)
  if (!is.null(contrasts_used_full)) {
    attr(X, "contrasts") <- contrasts_used_full
  }

  if (any(!is.finite(y_vec)) || any(!is.finite(X))) {
    stop("Non-finite values in response/predictors after NA handling.")
  }
  storage.mode(X) <- "double"

  n <- nrow(X); p <- ncol(X); nobs <- n; nparm <- p + 1L
  numerically_constant <- .svem_numerically_constant_columns(X)
  structural_intercept_only <- length(numerically_constant) == p ||
    (fam_name == "gaussian" && length(unique(y_vec)) < 2L) ||
    (fam_name == "binomial" && length(unique(y_vec)) < 2L)

  ## Capture a compact, self-contained formula environment for prediction.
  compact_formula <- .svem_compact_formula_terms(
    f_use, terms_full, data_names = names(schema_data)
  )
  f_store <- compact_formula$formula
  terms_clean <- .svem_freeze_factor_predvars(compact_formula$terms, schema_data)
  schema_data_used <- schema_data
  omitted_rows <- as.integer(attr(mf, "na.action"))
  if (length(omitted_rows) && nrow(schema_data_used) >= max(omitted_rows)) {
    schema_data_used <- schema_data_used[-omitted_rows, , drop = FALSE]
  }

  predictor_vars <- .svem_raw_predictor_vars(terms_full, names(schema_data))
  if (length(mixture_groups)) {
    mixture_vars <- unlist(lapply(mixture_groups, `[[`, "vars"), use.names = FALSE)
    if (length(setdiff(mixture_vars, predictor_vars)) ||
        length(intersect(mixture_vars, c(blocking, spec$settings$discrete_numeric))) ||
        !all(vapply(schema_data_used[mixture_vars], is.numeric, logical(1)))) {
      stop("Mixture components must be continuous, non-blocking model predictors.", call. = FALSE)
    }
    .svem_check_mixture_totals(schema_data_used, mixture_groups)
  }
  if (!using_spec) {
    schema_data_used <- .svem_drop_unused_raw_levels(schema_data_used, predictor_vars)
  }

  ## ---- guard: response must not also be a predictor ----
  ## .svem_response_vars() handles transformed responses such as log(y),
  ## for which as.character(f_use[[2L]]) is length 2 and would crash the
  ## '&&' below on R >= 4.3.
  resp_vars <- .svem_response_vars(f_use)
  resp_in_rhs <- intersect(resp_vars, predictor_vars)
  if (length(resp_in_rhs)) {
    stop(
      "SVEMnet does not allow a predictor with the same name as the response ('",
      paste(resp_in_rhs, collapse = "', '"),
      "'). Please rename your variables or adjust the formula."
    )
  }

  var_classes <- setNames(vapply(predictor_vars, function(v) {
    if (v %in% names(data)) class(data[[v]])[1] else NA_character_
  }, character(1)), predictor_vars)

  xlevels <- list()
  factor_levels <- list()
  for (v in predictor_vars) {
    if (v %in% names(schema_data_used)) {
      if (is.factor(schema_data_used[[v]])) {
        xlevels[[v]] <- levels(schema_data_used[[v]])
        lev <- levels(schema_data_used[[v]])
        factor_levels[[v]] <- lev[lev %in% as.character(schema_data_used[[v]])]
      } else if (is.character(schema_data_used[[v]])) {
        lev <- sort(unique(as.character(schema_data_used[[v]])))
        xlevels[[v]] <- lev
        factor_levels[[v]] <- lev
      } else if (is.logical(schema_data_used[[v]])) {
        factor_levels[[v]] <- c("FALSE", "TRUE")[c(FALSE, TRUE) %in% schema_data_used[[v]]]
      }
    }
  }

  # Use the contrasts captured before subsetting (fix #1)
  contrasts_used <- contrasts_used_full

  # Helper for categorical mode (used for blocking factors)
  .mode_level <- function(x) {
    if (is.null(x)) return(NA_character_)
    x_chr <- as.character(x)
    x_chr <- x_chr[!is.na(x_chr)]
    if (!length(x_chr)) return(NA_character_)
    lev <- if (is.factor(x)) levels(x) else sort(unique(x_chr))
    tab <- table(factor(x_chr, levels = lev))
    if (!length(tab)) return(NA_character_)
    maxfreq <- max(tab)
    modes <- names(tab)[tab == maxfreq]
    modes[1L]  # break ties arbitrarily but deterministically
  }

  # Normalize blocking to be only predictors actually present
  blocking <- intersect(blocking, predictor_vars)

  # Modes for blocking categorical variables, if any
  block_cat_modes <- NULL
  if (length(blocking)) {
    # "Categorical" here means we stored factor_levels
    block_cat <- intersect(blocking, names(factor_levels))
    if (length(block_cat)) {
      block_cat_modes <- vector("list", length(block_cat))
      names(block_cat_modes) <- block_cat
      for (v in block_cat) {
        if (v %in% colnames(mf)) {
          mode_v <- .mode_level(mf[[v]])
          if (!is.na(mode_v)) {
            block_cat_modes[[v]] <- mode_v
          }
        }
      }
      # Drop if everything came back NA / empty
      if (!length(block_cat_modes) ||
          all(vapply(block_cat_modes, is.null, logical(1)))) {
        block_cat_modes <- NULL
      }
    }
  }

  # Numeric ranges for sampling:
  # If a bigexp_spec was used and it stored num_range, prefer the locked ranges
  # (this does not affect fitting; it only affects downstream random sampling).
  if (using_spec &&
      !is.null(spec$num_range) &&
      is.matrix(spec$num_range) &&
      ncol(spec$num_range) > 0L &&
      !is.null(colnames(spec$num_range))) {

    # Determine numeric vars using the spec typing (if available), else fall back to mf
    if (!is.null(spec$is_cat) && length(spec$is_cat)) {
      num_vars <- predictor_vars[predictor_vars %in% names(spec$is_cat) & !spec$is_cat[predictor_vars]]
    } else {
      num_vars <- predictor_vars[vapply(predictor_vars, function(v) {
        v %in% names(schema_data_used) && is.numeric(schema_data_used[[v]])
      }, logical(1))]
    }

    cols <- intersect(num_vars, colnames(spec$num_range))
    if (length(cols)) {
      num_ranges <- spec$num_range[, cols, drop = FALSE]
      storage.mode(num_ranges) <- "double"
    } else {
      num_ranges <- matrix(numeric(0), nrow = 2, ncol = 0,
                           dimnames = list(c("min","max"), NULL))
    }

  } else {
    num_vars <- predictor_vars[vapply(predictor_vars, function(v) {
      v %in% names(schema_data_used) && is.numeric(schema_data_used[[v]])
    }, logical(1))]
    if (length(num_vars)) {
      rng_mat <- vapply(num_vars, function(v) {
        r <- range(schema_data_used[[v]], na.rm = TRUE)
        if (!all(is.finite(r)) || r[1] == r[2]) {
          r <- c(min(schema_data_used[[v]], na.rm = TRUE), max(schema_data_used[[v]], na.rm = TRUE))
        }
        r
      }, numeric(2))
      rownames(rng_mat) <- c("min","max")
      num_ranges <- as.matrix(rng_mat)
    } else {
      num_ranges <- matrix(numeric(0), nrow = 2, ncol = 0,
                           dimnames = list(c("min","max"), NULL))
    }
  }

  # Locked categorical support applies even to factor-only specifications.
  if (using_spec && !is.null(spec$levels) && length(spec$levels)) {
    fl <- spec$levels[intersect(names(spec$levels), predictor_vars)]
    fl <- fl[!vapply(fl, is.null, logical(1))]
    if (length(fl)) factor_levels <- fl
  }

  # Discrete numeric levels recorded in a bigexp_spec (for downstream sampling only)
  discrete_numeric_vars   <- character(0L)
  discrete_numeric_levels <- list()
  if (using_spec && !is.null(spec$settings)) {
    dn <- spec$settings$discrete_numeric
    dl <- spec$settings$discrete_levels

    if (is.null(dn)) dn <- character(0L)
    if (!is.character(dn)) dn <- character(0L)
    dn <- intersect(unique(dn), predictor_vars)

    if (!is.null(dl) && is.list(dl) && length(dn)) {
      # keep only named entries matching dn
      dl <- dl[intersect(names(dl), dn)]
      # coerce each to numeric unique sorted finite
      for (nm in names(dl)) {
        vv <- as.numeric(dl[[nm]])
        vv <- sort(unique(vv[is.finite(vv)]))
        if (length(vv)) discrete_numeric_levels[[nm]] <- vv
      }
    }

    if (length(dn)) discrete_numeric_vars <- dn
  }

  direct_support <- .svem_direct_factor_support(terms_clean, schema_data_used)
  for (v in setdiff(names(direct_support), names(discrete_numeric_levels))) {
    discrete_numeric_levels[[v]] <- direct_support[[v]]
  }
  discrete_numeric_vars <- union(discrete_numeric_vars, names(direct_support))



  ## ---- objective selection (auto) ----
  auto_used     <- identical(objective, "auto")
  auto_decision <- NA_character_

  objective_used <- if (auto_used) {
    if (fam_name == "gaussian") "wAIC" else "wBIC"
  } else {
    objective
  }
  auto_decision <- if (auto_used) objective_used else NA_character_


  ## ---- drop ridge when relaxed ----
  dropped_alpha0_for_relaxed <- FALSE
  if (isTRUE(relax_flag) && any(glmnet_alpha == 0)) {
    warning("Dropping alpha = 0 (ridge) for relaxed fits; ridge + relaxed is not supported here.")
    glmnet_alpha <- glmnet_alpha[glmnet_alpha != 0]
    if (!length(glmnet_alpha)) glmnet_alpha <- 1
    dropped_alpha0_for_relaxed <- TRUE
  }

  ## ---- containers ----
  coef_matrix  <- matrix(NA_real_, nrow = nBoot, ncol = p + 1L)
  colnames(coef_matrix) <- c("(Intercept)", colnames(X))
  best_alphas       <- rep(NA_real_, nBoot)
  best_lambdas      <- rep(NA_real_, nBoot)
  best_relax_gammas <- rep(NA_real_, nBoot)
  k_sel_vec         <- rep(NA_integer_, nBoot)
  edf_sel_vec       <- rep(NA_real_, nBoot)
  fallbacks         <- integer(nBoot)
  fallback_reason   <- rep(NA_character_, nBoot)
  n_eff_keep        <- numeric(nBoot)
  objective_value   <- rep(NA_real_, nBoot)
  validation_sse    <- rep(NA_real_, nBoot)
  solver_warnings   <- rep(NA_character_, nBoot)
  frw_uniforms      <- if (store_member_weights) matrix(NA_real_, nBoot, n) else NULL
  training_weights  <- if (store_member_weights) matrix(NA_real_, nBoot, n) else NULL
  validation_weights <- if (store_member_weights) matrix(NA_real_, nBoot, n) else NULL

  setup <- .svem_glmnet_setup(X, list(...), fam_name, relax_flag,
                             numerically_constant = numerically_constant)
  glmnet_design <- setup$design
  X_glm <- glmnet_design$x
  glmnet_base_args <- setup$args

  ## ---- bootstrap loop ----
  for (i in seq_len(nBoot)) {
    eps <- .Machine$double.eps
    U <- rep(NA_real_, n)
    if (weight_scheme == "SVEM") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- -log1p(-U)
    } else if (weight_scheme == "FRW_plain") {
      U <- pmin(pmax(stats::runif(n), eps), 1 - eps)
      w_train <- -log(U); w_valid <- w_train
    } else {
      w_train <- rep(1, n); w_valid <- rep(1, n)
    }
    w_train <- w_train * (n / sum(w_train))
    w_valid <- w_valid * (n / sum(w_valid))
    if (store_member_weights) {
      frw_uniforms[i, ]       <- U
      training_weights[i, ]   <- w_train
      validation_weights[i, ] <- w_valid
    }

    sumw  <- sum(w_valid); sumw2 <- sum(w_valid^2)
    n_eff_raw <- (sumw^2) / (sumw2 + eps)
    n_eff_adm <- max(2, min(n, n_eff_raw))
    n_eff_keep[i] <- n_eff_raw

    if (structural_intercept_only) {
      mu_w <- sum(w_train * y_vec) / sum(w_train)
      intercept <- if (fam_name == "gaussian") {
        mu_w
      } else {
        stats::qlogis(pmin(pmax(mu_w, 1e-6), 1 - 1e-6))
      }
      coef_matrix[i, ] <- c(intercept, rep(0, p))
      best_alphas[i] <- NA_real_
      best_lambdas[i] <- NA_real_
      best_relax_gammas[i] <- if (isTRUE(relax_flag)) NA_real_ else 1
      k_sel_vec[i] <- 1L
      edf_sel_vec[i] <- if (fam_name == "gaussian") 1 else NA_real_
      if (fam_name == "gaussian") {
        sse_i <- sum(w_valid * (y_vec - intercept)^2)
        validation_sse[i] <- sse_i
        objective_value[i] <- switch(
          objective_used,
          "wSSE" = sse_i,
          "wAIC" = sumw * log(max(sse_i, .Machine$double.xmin) / sumw) + 2,
          "wBIC" = sumw * log(max(sse_i, .Machine$double.xmin) / sumw) + log(n_eff_adm)
        )
      }
      next
    }

    if (fam_name == "gaussian") {
      edf_context <- .svem_edf_context(X, y_vec, w_train,
                                     glmnet_base_args, ncol(X_glm))
      if (use_edf && !is.null(edf_context$reason)) stop(edf_context$reason)
    }

    best_val_score <- Inf
    best_alpha     <- NA_real_
    best_lambda    <- NA_real_
    best_beta_hat  <- rep(NA_real_, p + 1L)
    best_k         <- NA_integer_
    best_kpen      <- NA_real_
    best_edf       <- NA_real_
    best_val_core  <- NA_real_
    best_relax_g   <- if (isTRUE(relax_flag)) NA_real_ else 1
    attempt_reasons <- character(0L)

    relax_gamma_grid <- if (isTRUE(relax_flag)) c( 0.2, 0.6,  1) else 1

    for (alpha in glmnet_alpha) {
      bundle <- .svem_glmnet_path(X_glm, y_vec, w_train, alpha,
                                  glmnet_base_args, relax_gamma_grid,
                                  glmnet_design$sentinel)
      attempt_reasons <- c(attempt_reasons, bundle$errors)
      notes <- unique(c(if (!is.na(solver_warnings[i])) solver_warnings[i], bundle$warnings))
      if (length(notes)) solver_warnings[i] <- paste(notes, collapse = " | ")
      fit <- bundle$fit
      for (path in bundle$paths) {
        relax_g <- path$gamma
        coef_path <- path$coefficients
        if (nrow(coef_path) != p + 1L) next
        L <- ncol(coef_path); if (L == 0L) next

        pred_link <- tryCatch({
          XB <- X %*% coef_path[-1L, , drop = FALSE]
          sweep(XB, 2, coef_path[1L, ], FUN = "+")
        }, error = function(e) NULL)
        if (is.null(pred_link) || nrow(pred_link) != n) next

        ## family-specific validation loss
        if (fam_name == "gaussian") {
          res <- pred_link - y_vec
          val_core <- as.vector(crossprod(w_valid, res^2))  # weighted SSE
          val_core[!is.finite(val_core)] <- Inf
          adj_val_core <- pmax(val_core, .Machine$double.xmin)
        } else {  # binomial: negative log-likelihood
          mu <- 1 / (1 + exp(-pred_link))
          mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
          yy <- matrix(y_vec, nrow = n, ncol = L)
          loglik <- colSums(w_valid * (yy * log(mu) + (1 - yy) * log(1 - mu)))
          val_core <- -loglik   # NLL
          val_core[!is.finite(val_core)] <- Inf
          adj_val_core <- pmax(val_core, .Machine$double.eps)
        }

        k_raw <- as.integer(.svem_solver_support_path(coef_path))
        if (length(k_raw) != length(adj_val_core)) next

        k_slope <- pmax(0L, k_raw - 1L)
        k_eff   <- 1L + k_slope
        logN_pen <- log(n_eff_adm)

        ## The "support" penalty uses k_eff. Conditional "edf" uses one
        ## plus the numerical rank of the active weighted, centered design
        ## at alpha = 1, equal to support only for independent active columns.
        ## At alpha < 1 it uses the gamma-blended conditional trace.
        k_pen <- k_eff
        edf_path <- if (fam_name == "gaussian") as.numeric(k_eff) else
          rep(NA_real_, L)
        if (use_edf && length(fit$lambda) == L) {
          edf_path <- .svem_conditional_edf(edf_context, coef_path,
                                           fit$lambda, alpha, relax_g)
          k_pen <- edf_path
        }

        if (fam_name == "gaussian") {
          n_like <- sumw
          mse_w  <- adj_val_core / n_like
          metric <- switch(objective_used,
                           "wSSE" = val_core,
                           "wAIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- n_like * log(mse_w[mask]) + 2 * k_pen[mask]
                             out
                           },
                           "wBIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- n_like * log(mse_w[mask]) +
                               logN_pen * k_pen[mask]
                             out
                           },
                           stop("Unknown objective: ", objective_used))
        } else {
          ## binomial: objective based on NLL
          metric <- switch(objective_used,
                           "wSSE" = adj_val_core,
                           "wAIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- 2 * adj_val_core[mask] + 2 * k_eff[mask]
                             out
                           },
                           "wBIC" = {
                             out <- rep(Inf, L)
                             mask <- (k_slope < n_eff_adm)
                             out[mask] <- 2 * adj_val_core[mask] +
                               logN_pen * k_eff[mask]
                             out
                           },
                           stop("Unknown objective: ", objective_used))
        }

        metric[!is.finite(metric)] <- Inf
        if (!any(is.finite(metric))) next

        idx_min    <- which.min(metric)
        lambda_opt <- fit$lambda[idx_min]
        val_score  <- metric[idx_min]
        selected_edf <- if (fam_name == "gaussian") {
          as.numeric(edf_path[idx_min])
        } else NA_real_
        if (fam_name == "gaussian" && !use_edf && length(fit$lambda) == L) {
          selected_edf <- .svem_conditional_edf(
            edf_context, coef_path[, idx_min, drop = FALSE],
            fit$lambda[idx_min], alpha, relax_g)
        }

        if (is.finite(val_score) && val_score < best_val_score) {
          best_val_score <- val_score
          best_alpha     <- alpha
          best_lambda    <- lambda_opt
          best_beta_hat  <- coef_path[, idx_min]
          best_k         <- k_raw[idx_min]
          best_kpen      <- as.numeric(k_pen[idx_min])
          best_edf       <- selected_edf
          best_val_core  <- as.numeric(val_core[idx_min])
          best_relax_g   <- relax_g
        }
      }
    }

    if (anyNA(best_beta_hat) || !all(is.finite(best_beta_hat))) {
      fallbacks[i]  <- 1L
      attempt_reasons <- unique(attempt_reasons[nzchar(attempt_reasons)])
      fallback_reason[i] <- if (length(attempt_reasons)) {
        paste(attempt_reasons, collapse = " | ")
      } else {
        "no admissible finite glmnet path point"
      }
      mu_w          <- sum(w_train * y_vec) / sum(w_train)
      if (fam_name == "gaussian") {
        best_beta_hat <- c(mu_w, rep(0, p))
      } else {
        mu_w <- pmin(pmax(mu_w, 1e-6), 1 - 1e-6)
        best_beta_hat <- c(qlogis(mu_w), rep(0, p))
      }
      best_alpha    <- NA_real_
      best_lambda   <- NA_real_
      best_k        <- 1L
      best_kpen     <- 1
      best_edf      <- if (fam_name == "gaussian") 1 else NA_real_
      if (fam_name == "gaussian") {
        best_val_core <- sum(w_valid * (y_vec - best_beta_hat[1L])^2)
      } else {
        best_val_core <- NA_real_
      }
      best_relax_g  <- if (isTRUE(relax_flag)) NA_real_ else 1
    }

    coef_matrix[i, ]     <- best_beta_hat
    best_alphas[i]       <- best_alpha
    best_lambdas[i]      <- best_lambda
    best_relax_gammas[i] <- best_relax_g
    k_sel_vec[i]         <- best_k
    edf_sel_vec[i]       <- if (fam_name == "gaussian") best_edf else NA_real_
    objective_value[i]   <- if (fallbacks[i] == 1L) NA_real_ else best_val_score
    validation_sse[i]    <- if (fam_name == "gaussian") best_val_core else NA_real_
  }

  fallback_summary <- if (any(fallbacks == 1L)) {
    sort(table(fallback_reason[fallbacks == 1L]), decreasing = TRUE)
  } else {
    integer(0L)
  }
  if (!structural_intercept_only && all(fallbacks == 1L)) {
    stop(
      "All bootstrap fits fell back to intercept-only models. Reasons: ",
      paste(names(fallback_summary), as.integer(fallback_summary),
            sep = " (n=", collapse = "); "),
      ").",
      call. = FALSE
    )
  }
  if (any(fallbacks == 1L) && !all(fallbacks == 1L)) {
    warning(
      sum(fallbacks == 1L), " of ", nBoot,
      " bootstrap fits fell back to intercept-only models. See `diagnostics$fallback_reasons`.",
      call. = FALSE
    )
  }

  ## ---- finalize ----
  valid_rows <- rowSums(!is.finite(coef_matrix)) == 0
  if (!any(valid_rows)) stop("All bootstrap iterations failed to produce valid coefficients.")
  coef_matrix        <- coef_matrix[valid_rows, , drop = FALSE]
  best_alphas       <- best_alphas[valid_rows]
  best_lambdas      <- best_lambdas[valid_rows]
  best_relax_gammas <- best_relax_gammas[valid_rows]
  k_sel_vec         <- k_sel_vec[valid_rows]
  edf_sel_vec       <- edf_sel_vec[valid_rows]
  fallbacks         <- fallbacks[valid_rows]
  fallback_reason   <- fallback_reason[valid_rows]
  n_eff_keep        <- n_eff_keep[valid_rows]
  objective_value   <- objective_value[valid_rows]
  validation_sse    <- validation_sse[valid_rows]
  solver_warnings   <- solver_warnings[valid_rows]
  if (store_member_weights) {
    frw_uniforms       <- frw_uniforms[valid_rows, , drop = FALSE]
    training_weights   <- training_weights[valid_rows, , drop = FALSE]
    validation_weights <- validation_weights[valid_rows, , drop = FALSE]
  }

  nBoot_used <- nrow(coef_matrix)  # fix #2

  avg_coefficients <- colMeans(coef_matrix)
  lin_pred <- as.vector(X %*% avg_coefficients[-1L] + avg_coefficients[1L])

  if (fam_name == "gaussian") {
    y_pred <- lin_pred
  } else {
    y_pred <- tryCatch({
      cm <- coef_matrix
      if (!is.matrix(cm) || nrow(cm) < 1L || ncol(cm) != (ncol(X) + 1L)) {
        stop("coef_matrix shape mismatch")
      }

      intercepts <- cm[, 1L]
      betas      <- cm[, -1L, drop = FALSE]

      eta_mat <- X %*% t(betas) + matrix(intercepts, nrow = nrow(X), ncol = length(intercepts), byrow = TRUE)

      # prevent exp overflow/underflow
      eta_mat <- pmin(pmax(eta_mat, -30), 30)

      p_mat <- 1 / (1 + exp(-eta_mat))
      p_mat <- pmin(pmax(p_mat, 1e-8), 1 - 1e-8)

      rowMeans(p_mat)
    }, error = function(e) {
      p <- 1 / (1 + exp(-lin_pred))
      pmin(pmax(p, 1e-8), 1 - 1e-8)
    })
  }


  ## debias: only for gaussian (leave binomial as-is for now)
  debias_fit <- NULL
  y_pred_debiased <- NULL
  parms_debiased  <- avg_coefficients
  calibration <- list(fit = NULL, reason = NULL)

  if (fam_name == "gaussian" && nBoot_used >= 10 && stats::var(y_pred) > 0) {  # fix #2
    calibration <- .svem_fit_calibration(y_vec, y_pred, X, avg_coefficients)
    debias_fit <- calibration$fit
    y_pred_debiased <- if (!is.null(debias_fit)) stats::predict(debias_fit) else NULL
    parms_debiased <- avg_coefficients
    ab <- try(.svem_calibration_coefficients(debias_fit), silent = TRUE)
    if (!inherits(ab, "try-error") &&
        length(ab) >= 2 &&
        is.finite(ab[1]) && is.finite(ab[2])) {
      a <- unname(ab[1]); b <- unname(ab[2])
      int_name <- "(Intercept)"
      if (!is.null(names(parms_debiased)) && int_name %in% names(parms_debiased)) {
        parms_debiased[int_name] <- a + b * parms_debiased[int_name]
        if (length(parms_debiased) > 1L) {
          slope_names <- setdiff(names(parms_debiased), int_name)
          parms_debiased[slope_names] <- b * parms_debiased[slope_names]
        }
      } else {
        parms_debiased[1] <- a + b * parms_debiased[1]
        if (length(parms_debiased) > 1L) parms_debiased[-1] <- b * parms_debiased[-1]
      }
    }
  }

  diagnostics <- list(
    numerically_constant_columns = colnames(X)[numerically_constant],
    calibration_skip_reason = calibration$reason,
    solver_warnings = solver_warnings,
    solver_warning_count = sum(!is.na(solver_warnings)),
    glmnet_version = as.character(utils::packageVersion("glmnet")),
    glmnet_control = if (!is.null(glmnet_base_args$control)) glmnet_base_args$control else
      c(glmnet_base_args[intersect(names(glmnet_base_args), .glmnet_direct_control_args())],
        attr(glmnet_base_args, "svem_legacy_control", exact = TRUE)),
    k_summary = c(k_median = stats::median(k_sel_vec),
                  k_iqr    = stats::IQR(k_sel_vec)),
    fallback_rate = mean(fallbacks),
    fallback_reasons = fallback_summary,
    n_eff_summary = summary(n_eff_keep),
    edf_summary = if (fam_name == "gaussian") summary(edf_sel_vec) else NULL,
    relax_gamma_freq =
      if (isTRUE(relax_flag) && any(is.finite(best_relax_gammas))) {
        prop.table(table(round(best_relax_gammas, 2)))
      } else NULL
  )
  tf <- table(best_alphas[is.finite(best_alphas)])
  diagnostics$alpha_freq <- if (length(tf)) {
    out <- as.numeric(tf) / sum(tf)
    names(out) <- names(tf)
    out
  } else numeric()

  n_eff_adm_keep <- pmax(2, pmin(n, n_eff_keep))
  residual_df_support <- n - k_sel_vec
  residual_variance_support <- if (fam_name == "gaussian") {
    ifelse(residual_df_support > 0 & is.finite(validation_sse),
           validation_sse / residual_df_support, NA_real_)
  } else rep(NA_real_, nBoot_used)
  residual_scale_support <- sqrt(residual_variance_support)
  member_warning <- rep(NA_character_, nBoot_used)
  for (ii in seq_len(nBoot_used)) {
    wi <- character(0L)
    if (fallbacks[ii] == 1L) {
      wi <- c(wi, paste0("fallback: ", fallback_reason[ii]))
    }
    if (fam_name == "gaussian" && residual_df_support[ii] <= 0) {
      wi <- c(wi, "nonpositive n - support residual denominator")
    }
    if (objective_used %in% c("wAIC", "wBIC") &&
        (n_eff_adm_keep[ii] - (k_sel_vec[ii] - 1L)) <= 1) {
      wi <- c(wi, "selected support is within one coefficient of the Kish guardrail")
    }
    if (length(wi)) member_warning[ii] <- paste(wi, collapse = "; ")
  }
  member_diagnostics <- data.frame(
    member = seq_len(nBoot_used),
    objective = rep(objective_used, nBoot_used),
    validation_sse = validation_sse,
    objective_value = objective_value,
    support_size = k_sel_vec,
    effective_df = edf_sel_vec,
    kish_validation_n_eff = n_eff_keep,
    kish_validation_n_eff_adm = n_eff_adm_keep,
    residual_df_support = residual_df_support,
    residual_variance_support = residual_variance_support,
    residual_scale_support = residual_scale_support,
    admissible_for_residual_scale =
      fam_name == "gaussian" & residual_df_support > 0 &
      is.finite(residual_scale_support) & fallbacks == 0L,
    fallback = fallbacks == 1L,
    admissibility_warning = member_warning,
    stringsAsFactors = FALSE
  )
  member_weights <- if (store_member_weights) {
    list(
      uniforms = frw_uniforms,
      training = training_weights,
      validation = validation_weights,
      normalization = "each weight row rescaled to mean one"
    )
  } else NULL

  feature_names <- colnames(X)
  terms_str <- tryCatch(
    paste(deparse(stats::delete.response(terms_clean)), collapse = " "),
    error = function(e) NA_character_
  )
  safe_hash <- function(s) {
    if (!is.character(s) || !length(s) || is.na(s)) return(NA_character_)
    bytes <- charToRaw(paste0(s, collapse = ""))
    sprintf("h%08x_%d", sum(as.integer(bytes)), length(bytes))
  }
  schema <- list(
    model_frame       = model_frame_schema,
    feature_names     = feature_names,
    terms_str         = terms_str,
    xlevels           = xlevels,
    contrasts         = contrasts_used,
    contrasts_options = contrasts_opts_used,
    terms_hash        = safe_hash(terms_str)
  )

  sampling_schema <- list(
    predictor_vars  = predictor_vars,
    var_classes     = var_classes,
    num_ranges      = num_ranges,
    factor_levels   = factor_levels,
    blocking        = blocking,        # character(0) when no blocking
    block_cat_modes = block_cat_modes, # NULL when none / not available
    discrete_numeric = discrete_numeric_vars,
    discrete_levels  = discrete_numeric_levels
  )
  if (length(mixture_groups)) sampling_schema$mixture_groups <- mixture_groups


  result <- list(
    parms             = avg_coefficients,
    parms_debiased    = parms_debiased,
    debias_fit        = debias_fit,
    coef_matrix       = coef_matrix,
    nBoot             = nBoot_used,     # fix #2 (actual used)
    nBoot_input       = nBoot_input,    # keep what user asked for (added, non-breaking)
    glmnet_alpha      = glmnet_alpha,
    best_alphas       = best_alphas,
    best_lambdas      = best_lambdas,
    best_relax_gammas = if (isTRUE(relax_flag)) best_relax_gammas else
      rep(NA_real_, length(best_alphas)),
    weight_scheme     = weight_scheme,
    complexity        = complexity,
    best_edfs         = edf_sel_vec,
    member_diagnostics = member_diagnostics,
    member_weights    = member_weights,
    store_member_weights = store_member_weights,
    relaxed           = isTRUE(relax_flag),
    relaxed_input     = relaxed_input,
    dropped_alpha0_for_relaxed = dropped_alpha0_for_relaxed,
    objective_input   = objective,
    objective_used    = objective_used,
    objective         = objective_used,
    auto_used         = auto_used,
    auto_decision     = if (auto_used) auto_decision else NA_character_,
    diagnostics       = diagnostics,
    actual_y          = y_vec,
    training_X        = X,
    y_pred            = y_pred,
    y_pred_debiased   = y_pred_debiased,
    nobs              = nobs,
    nparm             = nparm,
    formula           = f_store,
    terms             = terms_clean,
    xlevels           = xlevels,
    contrasts         = contrasts_used,
    schema            = schema,
    sampling_schema   = sampling_schema,
    used_bigexp_spec  = using_spec,
    family            = fam_name
  )

  cls <- c("svem_model", "SVEMnet")
  if (fam_name == "binomial") {
    cls <- c("svem_binomial", cls)
  }
  class(result) <- cls
  if (length(engines)) {
    result$companions <- .svem_fit_companions(result, schema_data_used, engines,
                                             companion_blocks, unseen, using_spec,
                                             companion_continuous_blocks)
  }
  result
}
