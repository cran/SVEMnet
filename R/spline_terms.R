#' Specify a Gaussian smooth surface on a declared recipe domain
#'
#' Construct the surface basis once, before fitting any responses. Each distinct
#' numeric site in the reference grid receives equal weight in the shared
#' numeric basis centering and scaling. Nominal coverage does not reweight it.
#' The grid need not contain a response or block column.
#'
#' @param formula A two-sided formula with one response name and raw study
#'   predictor names joined by `+`, for example `y ~ lipid + ratio + type`.
#'   With `y ~ .`, include only intended predictors, an optional response and
#'   declared blocks in `domain_grid`. Numerical coordinates enter a joint
#'   surface; factor, character and logical predictors are nominal. Ordered
#'   factors are also treated as nominal; use explicit numeric values to model
#'   an ordered physical scale.
#' @param domain_grid Feasible reference recipes defined independently of
#'   response values. Exact duplicate recipes count once. Set this grid before
#'   observing outcomes and retain the saved specification for later fits.
#'   The numerical surface uses distinct numerical sites; adding another
#'   nominal level does not change its numeric basis at a fixed surface rank. A smaller reference grid
#'   can be used to build a basis for a much larger prediction pool. The
#'   constructor caps distinct numeric reference sites at 65,536; use a smaller
#'   response-free reference grid when the scoring pool is larger.
#' @param blocking Optional nuisance-factor names. Known measured levels
#'   contribute ordinary penalized additive shifts, without interactions with
#'   study factors. Numeric block IDs are matched by their exact values.
#'   Unreplicated numeric IDs warn. Continuous blocking trends are unsupported.
#' @param mixture_groups Fixed-sum mixture groups declared as in [bigexp_terms()].
#'   Include every component. Multiple disjoint groups, legal zeros and fixed
#'   components are supported. Each varying group is mapped to its independent
#'   closure coordinates. The reference grid must span these coordinates.
#' @param bounds Optional two-row numeric matrix of lower and upper legal
#'   limits, with columns named for numeric study predictors. Defaults to the
#'   reference ranges. A coordinate fixed on the reference grid remains fixed
#'   in this saved model. Different legal limits do not change coordinate
#'   scaling or the curvature penalty. Prediction beyond reference coverage is
#'   extrapolation. Eligibility uses the package's absolute physical tolerance
#'   of `1e-8`, scaled by `max(1, abs(total))` for each mixture total.
#'   Additional preparation constraints belong in the candidate pool.
#' @param rank Advanced shared numeric-surface rank including its intercept. Default
#'   `NULL` retains rank 26 for the original one-to-three-coordinate surface.
#'   General surfaces use a quadratic-sized rank target, bounded by distinct
#'   reference sites and `max_columns`, while reserving nominal main effects.
#'   Marginal curvature columns are additional and share the overall column
#'   budget. This target does not imply an exact quadratic basis. With too few sites
#'   or a limited budget, the numeric surface can reduce to linear directions.
#' @param null_weight Positive penalty on linear directions and nominal main
#'   effects relative to the smallest positive curvature penalty; default 0.05.
#' @param level_modes Advanced maximum number of smooth curvature modes for
#'   each nominal level's deviation from the common surface; default 3.
#' @param level_penalty_ratio Positive multiplier for nominal-level deviation
#'   penalties, including level-specific numeric slopes; default 1.
#' @param max_columns Advanced budget for non-block model columns; default 200.
#'   Main effects are retained. Complete interaction families and curvature
#'   ranks are reduced when necessary and reported in `basis_notes`. Increase
#'   this explicit resource budget when more columns are needed. Memory and
#'   fitting cost increase with the column count. There is no fixed
#'   three-coordinate or twelve-coordinate restriction.
#' @param interaction_order Advanced maximum order, 1 to 3, of nominal cell
#'   effects; default 3. Third-order structure can also include two nominal
#'   factors times a numeric linear direction when the budget permits.
#' @param marginal_modes Advanced maximum curvature modes per physical numeric
#'   factor in general surfaces with at least four independent coordinates;
#'   default 4. These one-factor smooths complement the joint surface. Zero
#'   disables them. Whole ranks are reduced uniformly by `max_columns`.
#' @param marginal_weight Positive multiplier for marginal curvature penalties;
#'   default 0.25 relative to the first shared joint-curvature penalty.
#' @details The original numeric surface is retained exactly for one to three
#'   independent coordinates with at most one mixture and no fixed component.
#'   General specifications store version 2. Their shared numeric surface is
#'   thin-plate for one to three coordinates and a low-rank Duchon spline in
#'   higher dimensions, with a linear penalty null space. Higher-dimensional
#'   surfaces also include modest one-factor curvature to avoid spending every
#'   curvature direction on joint geometry. Its positive penalty allows overlap
#'   with joint modes. This is a practical basis choice, not a universal optimum.
#'   Nominal factors add
#'   penalized level shifts, low-rank deviations from the common surface and
#'   cell effects. These deviations shrink toward the common surface rather
#'   than fitting an independent high-dimensional surface at every level.
#'   Fixed coordinates use explicit point bounds or a scale-relative
#'   floating-point tolerance of eight machine epsilons. This removes rounding
#'   noise without using an absolute width threshold on a narrow varying factor.
#'
#'   With no varying numeric coordinate, the general model reduces to
#'   penalized categorical effects, or an intercept if every factor is fixed;
#'   no numeric spline is claimed in that limit. Such limits need no `mgcv`.
#'   Interaction effects specific to unmeasured combinations of otherwise
#'   observed levels are not identified; predictions there use shared and
#'   lower-order shrinkage assumptions. New or unmeasured nominal
#'   study levels have no prediction. Study factors can interact; declared
#'   blocks remain additive. The selected member means and their variation
#'   are not calibrated prediction intervals. Factor support does not imply
#'   that a sparse design identifies all interactions or narrow response peaks.
#'
#'   The optional `mgcv` package supplies the numeric spline construction.
#'   The stored basis, penalties and charts are reused in prediction and after
#'   serialization. General specifications retain their basis controls,
#'   legality tolerance and construction-package version. Numerical prediction
#'   at new sites requires `mgcv`.
#'   Non-default general basis controls can select version 2 even on a
#'   low-dimensional numeric problem; saved version 1 fits remain unchanged.
#'   The implementation was developed with assistance from Astra.
#' @return A serializable `svem_spline_spec` containing physical predictors,
#'   bounds, reference grid, coordinate transforms, basis and diagonal penalty.
#'   General specifications additionally record nominal levels, included
#'   interactions, the actual column count and any automatic budget reductions.
#'   Pass the specification to [svem_spline()].
#' @references Wood, S. N. (2003). Thin plate regression splines.
#'   Journal of the Royal Statistical Society B, 65, 95--114.
#'   \doi{10.1111/1467-9868.00374}.
#'   Duchon, J. (1977). Splines minimizing rotation-invariant semi-norms in
#'   Sobolev spaces. In Construction Theory of Functions of Several Variables,
#'   85--100. See also the `mgcv` documentation for
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.ds.smooth.spec.html}.
#' @examples
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   recipes <- expand.grid(u = seq(0, 1, length.out = 7),
#'                          v = seq(0, 1, length.out = 7))
#'   spec <- spline_terms(y ~ u + v, recipes)
#'   measured <- recipes[seq(1, 49, by = 2), ]
#'   measured$y <- with(measured, u * (1 - v))
#'   fit <- svem_spline(spec, measured, nBoot = 3)
#'   predict(fit, recipes[1:2, ])
#' }
#' @export
spline_terms <- function(formula, domain_grid, blocking = NULL,
                         mixture_groups = NULL, bounds = NULL,
                         rank = NULL, null_weight = 0.05,
                         level_modes = 3L, level_penalty_ratio = 1,
                         max_columns = 200L, interaction_order = 3L,
                         marginal_modes = 4L, marginal_weight = 0.25) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) || !is.data.frame(domain_grid) ||
      !nrow(domain_grid) || anyDuplicated(names(domain_grid))) {
    stop("Supply a two-sided formula with a named response and a nonempty domain grid.", call. = FALSE)
  }
  response <- as.character(formula[[2L]])
  if (!is.null(blocking) && (!is.character(blocking) || anyNA(blocking) ||
       any(!nzchar(blocking)) || anyDuplicated(blocking) || response %in% blocking)) {
    stop("blocking must name distinct nuisance factors, separate from the response.", call. = FALSE)
  }
  blocks <- as.character(blocking)
  raw_names <- function(expr) {
    if (is.symbol(expr) && !identical(expr, quote(.))) return(as.character(expr))
    if (is.call(expr) && identical(expr[[1L]], quote(`+`)) && length(expr) == 3L) {
      return(c(raw_names(expr[[2L]]), raw_names(expr[[3L]])))
    }
    stop("Use raw predictor names joined by +, or a single dot, for the spline surface.", call. = FALSE)
  }
  rhs <- formula[[3L]]
  vars <- if (identical(rhs, quote(.))) setdiff(names(domain_grid), response) else raw_names(rhs)
  if (response %in% vars || anyDuplicated(vars)) stop("Repeated predictors or response leakage in spline formula.", call. = FALSE)
  factors <- setdiff(vars, blocks)
  if (!length(factors) || !all(factors %in% names(domain_grid))) {
    stop("The domain grid must contain each declared surface predictor.", call. = FALSE)
  }
  max_columns <- .svem_integer_scalar(max_columns, "max_columns", min = 1L)
  level_modes <- .svem_integer_scalar(level_modes, "level_modes", min = 0L)
  interaction_order <- .svem_integer_scalar(interaction_order, "interaction_order", min = 1L, max = 3L)
  level_penalty_ratio <- .svem_numeric_scalar(level_penalty_ratio, "level_penalty_ratio", lower = 0, lower_open = TRUE)
  marginal_modes <- .svem_integer_scalar(marginal_modes, "marginal_modes", min = 0L)
  marginal_weight <- .svem_numeric_scalar(marginal_weight, "marginal_weight", lower = 0, lower_open = TRUE)
  groups <- .svem_mixture_groups(mixture_groups)
  if (!.svem_spline_v1_eligible(domain_grid[factors], groups) ||
      .svem_spline_nondefault(list(max_columns = max_columns, level_modes = level_modes,
        level_penalty_ratio = level_penalty_ratio, interaction_order = interaction_order,
        marginal_modes = marginal_modes, marginal_weight = marginal_weight))) {
    return(.svem_spline_general_terms(response, factors, blocks, domain_grid, groups,
      bounds = bounds, rank = rank, null_weight = null_weight,
      level_modes = level_modes, level_penalty_ratio = level_penalty_ratio,
      max_columns = max_columns, interaction_order = interaction_order,
      marginal_modes = marginal_modes, marginal_weight = marginal_weight))
  }
  if (is.null(rank)) rank <- 26L
  grid <- unique(domain_grid[factors])
  if (!all(vapply(grid, .svem_spline_numeric, logical(1))) || any(!is.finite(as.matrix(grid)))) {
    stop("Spline surface predictors must be finite, ordinary numeric physical coordinates.", call. = FALSE)
  }
  .svem_spline_check_sites(nrow(grid))
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the optional 'mgcv' package to use spline_terms().", call. = FALSE)
  }
  if (length(groups) > 1L) stop("The joint spline supports at most one fixed-sum mixture.", call. = FALSE)
  group <- if (length(groups)) groups[[1L]] else NULL
  parts <- group$vars
  if (length(parts) && (length(parts) < 2L || !all(parts %in% factors) || group$total <= 0 ||
                        any(group$lower < 0))) {
    stop("Specify all nonnegative mixture components in the surface and a positive mixture total.", call. = FALSE)
  }
  d <- length(factors) - as.integer(length(parts) > 0L)
  rank <- .svem_integer_scalar(rank, "rank", min = d + 2L)
  null_weight <- .svem_numeric_scalar(null_weight, "null_weight", lower = 0, lower_open = TRUE)
  if (!d %in% 1:3 || nrow(grid) < d + 2L) {
    stop("The spline needs one to three independent coordinates and at least dimension + 2 distinct reference recipes.", call. = FALSE)
  }
  chart <- apply(as.matrix(grid), 2L, range)
  if (is.null(bounds)) bounds <- chart
  if (!is.matrix(bounds) || !is.numeric(bounds) || nrow(bounds) != 2L ||
      anyDuplicated(colnames(bounds)) || !setequal(colnames(bounds), factors)) {
    stop("bounds must have two numeric rows and one named column per surface predictor.", call. = FALSE)
  }
  bounds <- bounds[, factors, drop = FALSE]
  if (any(!is.finite(bounds)) || any(bounds[2L, ] <= bounds[1L, ])) {
    stop("Each predictor needs distinct finite legal limits.", call. = FALSE)
  }
  widths <- chart[2L, ] - chart[1L, ]
  if (any(!is.finite(widths)) || any(widths <= 0)) {
    stop("Each spline coordinate must vary over a finite range on the reference grid.", call. = FALSE)
  }
  spec <- list(version = 1L, response = response, predictors = factors, blocking = blocks,
               process = setdiff(factors, parts), mixture_groups = groups,
               lower = stats::setNames(as.numeric(chart[1L, ]), factors),
               legal_lower = stats::setNames(as.numeric(bounds[1L, ]), factors),
               legal_upper = stats::setNames(as.numeric(bounds[2L, ]), factors),
               widths = stats::setNames(as.numeric(widths), factors), dimension = d,
               domain_grid = grid, requested_rank = rank, null_weight = null_weight)
  if (length(parts)) {
    a <- as.matrix(grid[parts])
    spec$mix_center <- colMeans(a)
    normal <- .svem_spline_mixture_normal(spec$widths[parts])
    if (any(!is.finite(normal)) || !any(normal != 0)) stop("Mixture widths are numerically unsupported.", call. = FALSE)
    v <- normal
    v[length(v)] <- v[length(v)] + 1
    h <- diag(length(v)) - 2 * tcrossprod(v) / sum(v^2)
    tangent <- h[, seq_len(length(v)-1L), drop = FALSE]
    centered <- sweep(sweep(a, 2L, spec$mix_center, "-"), 2L, spec$widths[parts], "/")
    axes <- eigen(crossprod(centered %*% tangent), symmetric = TRUE)
    spec$tangent <- tangent %*% axes$vectors
    spec$mix_radius <- max(sqrt(rowSums((centered %*% spec$tangent)^2)))
    if (!is.finite(spec$mix_radius) || spec$mix_radius <= 0) stop("Degenerate mixture reference domain.", call. = FALSE)
  }
  z <- .svem_spline_coordinates(spec, grid)
  if (qr(scale(as.matrix(z), center = TRUE, scale = FALSE))$rank < d) {
    stop("The reference grid does not span the declared independent coordinates.", call. = FALSE)
  }
  k <- min(rank, nrow(grid))
  # The smooth specification captures these names for lookup in z.
  .u <- .v <- .w <- NULL
  xt <- list(max.knots = 2000, seed = 17091)
  sm_term <- switch(d, mgcv::s(.u, bs = "tp", m = 2, k = k, xt = xt),
      mgcv::s(.u, .v, bs = "tp", m = 2, k = k, xt = xt),
      mgcv::s(.u, .v, .w, bs = "tp", m = 2, k = k, xt = xt))
  sm <- mgcv::smoothCon(sm_term, data = z, knots = NULL, absorb.cons = TRUE,
                        scale.penalty = FALSE)[[1L]]
  spec$center <- colMeans(sm$X)
  x <- sweep(sm$X, 2L, spec$center, "-")
  gram <- eigen(crossprod(x) / nrow(x), symmetric = TRUE)
  if (min(gram$values) <= max(gram$values) * max(dim(x)) * .Machine$double.eps) {
    stop("The reference grid does not identify the requested thin-plate basis.", call. = FALSE)
  }
  whiten <- sweep(gram$vectors, 2L, sqrt(gram$values), "/")
  rough <- crossprod(whiten, sm$S[[1L]] %*% whiten)
  eig <- eigen((rough + t(rough))/2, symmetric = TRUE)
  tolerance <- max(abs(eig$values)) * sqrt(.Machine$double.eps)
  linear <- abs(eig$values) <= tolerance
  if (any(eig$values < -tolerance) || sum(linear) != d || all(linear)) {
    stop("The roughness penalty does not identify the expected linear and curvature modes.", call. = FALSE)
  }
  spec$transform <- whiten %*% eig$vectors
  spec$penalty <- eig$values / min(eig$values[!linear])
  spec$penalty[linear] <- null_weight
  # Prediction evaluates the stored knots and thin-plate basis definition.
  # The centered grid design is already retained in reference_X.
  sm$X <- NULL
  spec$smooth <- sm
  spec$reference_X <- x %*% spec$transform
  spec$rank <- ncol(spec$reference_X) + 1L
  class(spec) <- "svem_spline_spec"
  spec
}

.svem_spline_numeric <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x))

.svem_spline_mixture_normal <- function(widths) {
  norm <- sqrt(sum(widths^2))
  if (!is.finite(norm) || norm < sqrt(.Machine$double.xmin)) {
    widths <- widths / max(widths)
    norm <- sqrt(sum(widths^2))
  }
  widths / norm
}

.svem_spline_legal <- function(spec, data) {
  if (identical(spec$version, 2L)) return(.svem_spline_general_legal(spec, data))
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(spec$predictors %in% names(data)) ||
      !all(vapply(data[spec$predictors], .svem_spline_numeric, logical(1)))) {
    stop("Supply all named numeric physical coordinates.", call. = FALSE)
  }
  x <- as.matrix(data[spec$predictors])
  tol <- 1e-8
  # Earlier saved specifications used the chart bounds as their legal limits.
  lower <- if (is.null(spec$legal_lower)) spec$lower else spec$legal_lower
  upper <- if (is.null(spec$legal_upper)) spec$lower + spec$widths else spec$legal_upper
  ok <- rowSums(!is.finite(x) | sweep(x, 2L, lower, "-") < -tol |
                sweep(x, 2L, upper, "-") > tol) == 0L
  if (length(spec$mixture_groups)) {
    g <- spec$mixture_groups[[1L]]
    a <- x[, g$vars, drop = FALSE]
    ok <- ok & abs(rowSums(a) - g$total) <= tol * max(1, abs(g$total)) &
      rowSums(sweep(a, 2L, g$lower, "-") < -tol |
              sweep(a, 2L, g$upper, "-") > tol) == 0L
  }
  ok[is.na(ok)] <- FALSE
  unname(ok)
}

.svem_spline_coordinates <- function(spec, data) {
  if (!all(.svem_spline_legal(spec, data))) {
    stop("Recipe violates the declared spline bounds or mixture total, or contains nonfinite settings.", call. = FALSE)
  }
  x <- as.matrix(data[spec$predictors])
  u <- sweep(sweep(x, 2L, spec$lower, "-"), 2L, spec$widths, "/")
  if (any(!is.finite(u))) stop("Physical coordinates exceed the numerical chart range.", call. = FALSE)
  z <- 2 * u[, spec$process, drop = FALSE] - 1
  if (length(spec$mixture_groups)) {
    g <- spec$mixture_groups[[1L]]
    a <- x[, g$vars, drop = FALSE]
    centered <- sweep(sweep(a, 2L, spec$mix_center, "-"), 2L, spec$widths[g$vars], "/")
    z <- cbind(z, centered %*% spec$tangent / spec$mix_radius)
  }
  colnames(z) <- c(".u", ".v", ".w")[seq_len(spec$dimension)]
  as.data.frame(z)
}

.svem_spline_surface <- function(spec, data) {
  if (identical(spec$version, 2L)) return(.svem_spline_general_design(spec, data))
  z <- .svem_spline_coordinates(spec, data)
  if (!nrow(data)) return(spec$reference_X[integer(), , drop = FALSE])
  if (identical(unname(as.matrix(data[spec$predictors])), unname(as.matrix(spec$domain_grid)))) {
    return(spec$reference_X)
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the optional 'mgcv' package to evaluate this saved spline at new recipes.", call. = FALSE)
  }
  mat <- mgcv::PredictMat(spec$smooth, z)
  sweep(mat, 2L, spec$center, "-") %*% spec$transform
}

#' @export
print.svem_spline_spec <- function(x, ...) {
  if (identical(x$version, 2L)) return(.svem_spline_general_print_spec(x, ...))
  cat("Joint thin-plate spline specification for ", x$response, "\n",
      "  ", nrow(x$domain_grid), " reference recipes, ", x$dimension,
      " independent coordinates, surface rank ", x$rank, "\n", sep = "")
  invisible(x)
}

