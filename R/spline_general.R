# Saved general Gaussian surface. Developed with assistance from Astra.
# Uses mgcv thin-plate/Duchon bases, with penalized nominal deviations and
# ordinary additive penalized blocks. Existing version-1 arithmetic is separate.
.svem_spline_general_constants <- function() {
  list(max_columns   = 200L,    # cap on model columns before block contrasts
       max_sites     = 65536L,  # cap on distinct reference coordinate sites
       tp_knots      = 2000,    # matches the version 1 construction
       ds_knots      = 1000,    # bounds the kernel memory at higher dimension
       knot_seed     = 17091,   # matches the version 1 construction
       tol           = 1e-8)
}

.svem_spline_check_sites <- function(n) {
  limit <- .svem_spline_general_constants()$max_sites
  if (n > limit) stop("The reference grid has more than ", limit,
    " distinct coordinate sites. Choose a smaller response-free reference grid.", call. = FALSE)
}

# Quadratic-complete budget: intercept, d linear and d(d+1)/2 curvature modes.
# This budget is a target. A truncated kernel basis does not exactly span the
# quadratics, so the probe measures the shortfall.
.svem_spline_general_default_rank <- function(d, max_columns = 200L) {
  if (!d) return(1L)
  max(d + 2L, min(max_columns + 1L, max(26L, (d + 1L) * (d + 2L) / 2L)))
}

.svem_spline_general_is_nominal <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

# Radix sorting keeps the saved column layout independent of the locale.
.svem_spline_general_levels <- function(x) {
  if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])), method = "radix")
}

.svem_spline_general_preserving_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  force(expr)
}

# The same closure chart as version 1, restricted to the components that vary.
.svem_spline_general_mixture_chart <- function(sites, vars, widths) {
  a <- sites[, vars, drop = FALSE]
  center <- colMeans(a)
  normal <- .svem_spline_mixture_normal(widths[vars])
  v <- normal
  v[length(v)] <- v[length(v)] + 1
  h <- diag(length(v)) - 2 * tcrossprod(v) / sum(v^2)
  tangent <- h[, seq_len(length(v) - 1L), drop = FALSE]
  centered <- sweep(sweep(a, 2L, center, "-"), 2L, widths[vars], "/")
  axes <- eigen(crossprod(centered %*% tangent), symmetric = TRUE)
  tangent <- tangent %*% axes$vectors
  radius <- max(sqrt(rowSums((centered %*% tangent)^2)))
  if (!is.finite(radius) || radius <= 0) stop("Degenerate mixture reference domain.", call. = FALSE)
  list(vars = vars, center = center, tangent = tangent, radius = radius)
}

# d <= 3 reuses the version 1 construction. Higher dimension uses a Duchon spline
# with a linear null space and the kernel r^2 log(r) (m = 2, s = d/2 - 1).
.svem_spline_general_smooth_term <- function(coordinates, k, d) {
  constants <- .svem_spline_general_constants()
  args <- if (d <= 3L) {
    list(bs = "tp", m = 2, k = k,
         xt = list(max.knots = constants$tp_knots, seed = constants$knot_seed))
  } else {
    list(bs = "ds", m = c(2, d / 2 - 1), k = k,
         xt = list(max.knots = constants$ds_knots, seed = constants$knot_seed))
  }
  eval(as.call(c(list(quote(mgcv::s)), lapply(coordinates, as.name), args)))
}

.svem_spline_general_terms <- function(response, predictors, blocks, domain_grid, groups,
                                       bounds = NULL, rank = NULL, null_weight = 0.05,
                                       level_modes = 3L, level_penalty_ratio = 1, max_columns = 200L,
                                       interaction_order = 3L, marginal_modes = 4L, marginal_weight = 0.25) {
  constants <- .svem_spline_general_constants()
  max_columns <- .svem_integer_scalar(max_columns, "max_columns", min = 1L)
  level_modes <- .svem_integer_scalar(level_modes, "level_modes", min = 0L)
  interaction_order <- .svem_integer_scalar(interaction_order, "interaction_order", min = 1L, max = 3L)
  level_penalty_ratio <- .svem_numeric_scalar(level_penalty_ratio, "level_penalty_ratio", lower = 0, lower_open = TRUE)
  marginal_modes <- .svem_integer_scalar(marginal_modes, "marginal_modes", min = 0L)
  marginal_weight <- .svem_numeric_scalar(marginal_weight, "marginal_weight", lower = 0, lower_open = TRUE)
  grid <- unique(domain_grid[predictors])
  is_num <- vapply(grid, .svem_spline_numeric, logical(1))
  is_nom <- vapply(grid, .svem_spline_general_is_nominal, logical(1))
  if (!all(is_num | is_nom)) {
    stop("Spline predictors must be ordinary numeric coordinates or factor, character or logical levels.", call. = FALSE)
  }
  nominal <- predictors[is_nom]
  continuous <- predictors[is_num]
  if (anyNA(grid[nominal])) stop("Reference nominal levels must be nonmissing.", call. = FALSE)
  sites_df <- unique(grid[continuous])
  if (!length(continuous)) sites_df <- grid[1L, continuous, drop = FALSE]
  sites <- as.matrix(sites_df)
  if (any(!is.finite(sites))) stop("Spline numeric coordinates must be finite.", call. = FALSE)
  chart <- if (length(continuous)) apply(sites, 2L, range) else matrix(numeric(), 2L, 0L)
  colnames(chart) <- continuous
  widths <- chart[2L, ] - chart[1L, ]
  # Use relative floating-point precision, not an absolute physical width cutoff.
  if (any(!is.finite(widths))) stop("Reference coordinate ranges exceed numerical precision.", call. = FALSE)
  fixed <- continuous[.svem_spline_general_fixed(chart)]
  if (is.null(bounds)) bounds <- chart
  if (!is.matrix(bounds) || !is.numeric(bounds) || nrow(bounds) != 2L ||
      anyDuplicated(colnames(bounds)) || !setequal(colnames(bounds), continuous)) {
    stop("bounds must have two numeric rows and one named column per numeric predictor.", call. = FALSE)
  }
  bounds <- bounds[, continuous, drop = FALSE]
  # Explicit point bounds identify a fixed setting even when a computed
  # balance differs by rounding noise. Other tiny varying factors stay varying.
  fixed <- union(fixed, continuous[bounds[2L, ] == bounds[1L, ]])
  if (any(!is.finite(bounds)) || any(bounds[2L, ] < bounds[1L, ])) {
    stop("Each numeric predictor needs finite ordered legal limits.", call. = FALSE)
  }
  # Rounding noise in a fixed coordinate cannot create additional sites or
  # alter the common numeric basis weighting.
  if (length(fixed)) {
    sites_df[fixed] <- as.list(chart[1L, fixed] / 2 + chart[2L, fixed] / 2)
    sites_df <- unique(sites_df)
    sites <- as.matrix(sites_df)
  }
  .svem_spline_check_sites(nrow(sites))
  parts <- unlist(lapply(groups, `[[`, "vars"), use.names = FALSE)
  charts <- list()
  for (g in groups) {
    if (length(g$vars) < 2L || !all(g$vars %in% continuous) || g$total <= 0 || any(g$lower < 0)) {
      stop("Specify all nonnegative mixture components in the surface and a positive mixture total.",
           call. = FALSE)
    }
    free <- setdiff(g$vars, fixed)
    # One varying component cannot satisfy a fixed total. Report that as closure.
    if (length(free) == 1L) {
      stop("Mixture components ", paste(g$vars, collapse = ", "),
           " must sum to their total; only ", free, " varies on the reference grid.", call. = FALSE)
    }
    if (length(free)) charts[[length(charts) + 1L]] <- .svem_spline_general_mixture_chart(sites, free, widths)
  }
  process <- setdiff(continuous, c(parts, fixed))
  d <- length(process) + sum(vapply(charts, function(m) length(m$vars) - 1L, integer(1)))
  d <- as.integer(d)
  if (d && nrow(sites) < d + 1L) {
    stop("The reference grid does not span the declared independent coordinates.", call. = FALSE)
  }
  if (is.null(rank)) {
    nominal_main <- sum(vapply(nominal, function(v) {
      n <- length(.svem_spline_general_levels(grid[[v]])); if (n > 1L) n else 0L
    }, integer(1)))
    if (d + nominal_main > max_columns)
      stop("Linear directions and nominal main effects exceed max_columns; increase that resource budget.", call. = FALSE)
    rank <- min(.svem_spline_general_default_rank(d, max_columns), max_columns - nominal_main + 1L)
  }
  rank <- .svem_integer_scalar(rank, "rank", min = if (d) d + 1L else 1L)
  if (d && rank - 1L > max_columns) {
    stop("The numeric surface exceeds max_columns; increase that explicit resource budget or lower rank.", call. = FALSE)
  }
  null_weight <- .svem_numeric_scalar(null_weight, "null_weight", lower = 0, lower_open = TRUE)
  spec <- list(version = 2L, response = response, predictors = predictors, blocking = blocks,
               continuous = continuous, nominal = nominal, process = process, fixed = fixed,
               fixed_value = stats::setNames(chart[1L, fixed] / 2 + chart[2L, fixed] / 2, fixed),
               mixture_groups = groups, charts = charts,
               lower = stats::setNames(as.numeric(chart[1L, ]), continuous),
               legal_lower = stats::setNames(as.numeric(bounds[1L, ]), continuous),
               legal_upper = stats::setNames(as.numeric(bounds[2L, ]), continuous),
               widths = stats::setNames(as.numeric(widths), continuous),
               levels = stats::setNames(lapply(nominal, function(v) .svem_spline_general_levels(grid[[v]])), nominal),
               dimension = d, coordinates = if (d) paste0(".z", seq_len(d)) else character(), domain_grid = grid,
               requested_rank = rank, null_weight = null_weight,
               controls = list(rank = rank, null_weight = null_weight, level_modes = level_modes,
                 level_penalty_ratio = level_penalty_ratio, max_columns = max_columns,
                 interaction_order = interaction_order, marginal_modes = marginal_modes,
                 marginal_weight = marginal_weight),
               legal_tol = constants$tol, mgcv_version = NA_character_,
               max_columns = max_columns, interaction_order = interaction_order,
               basis_notes = character(),
               smoother = if (!d) "categorical shrinkage" else if (d <= 3L) "thin-plate" else "Duchon")
  if (!all(.svem_spline_general_legal(spec, grid))) {
    stop("The reference grid violates the declared spline bounds or mixture total.", call. = FALSE)
  }
  if (!d) {
    spec["smooth"] <- list(NULL)
    spec$surface_penalty <- numeric()
    spec$rank <- 1L
    return(.svem_spline_general_nominal_layout(spec, level_modes, level_penalty_ratio))
  }
  z <- .svem_spline_general_chart_coordinates(spec, sites_df)
  if (qr(scale(as.matrix(z), center = TRUE, scale = FALSE))$rank < d) {
    stop("The reference grid does not span the declared independent coordinates. ",
         "Declare any fixed-sum closure in mixture_groups.", call. = FALSE)
  }
  # One distinct coordinate site counts once, however many levels it admits.
  # Adding a nominal level therefore cannot change the shared surface basis.
  k <- min(rank, nrow(sites))
  if (k <= d + 1L) {
    spec$center <- colMeans(z)
    x <- sweep(as.matrix(z), 2L, spec$center, "-")
    gram <- eigen(crossprod(x) / nrow(x), symmetric = TRUE)
    spec$transform <- sweep(gram$vectors, 2L, sqrt(gram$values), "/")
    spec["smooth"] <- list(NULL)
    spec$surface_penalty <- rep(null_weight, d)
    spec$rank <- d + 1L
    spec$smoother <- "linear limit"
    spec$basis_notes <- "Reference sites or the basis budget support linear directions only."
    return(.svem_spline_general_nominal_layout(spec, level_modes, level_penalty_ratio))
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the optional 'mgcv' package to construct a numeric spline.", call. = FALSE)
  }
  spec$mgcv_version <- as.character(utils::packageVersion("mgcv"))
  sm <- .svem_spline_general_preserving_rng(mgcv::smoothCon(.svem_spline_general_smooth_term(spec$coordinates, k, d), data = z,
                                           knots = NULL, absorb.cons = TRUE,
                                           scale.penalty = FALSE)[[1L]])
  spec$center <- colMeans(sm$X)
  x <- sweep(sm$X, 2L, spec$center, "-")
  gram <- eigen(crossprod(x) / nrow(x), symmetric = TRUE)
  if (min(gram$values) <= max(gram$values) * max(dim(x)) * .Machine$double.eps) {
    stop("The reference grid does not identify the requested surface basis.", call. = FALSE)
  }
  whiten <- sweep(gram$vectors, 2L, sqrt(gram$values), "/")
  rough <- crossprod(whiten, sm$S[[1L]] %*% whiten)
  eig <- eigen((rough + t(rough)) / 2, symmetric = TRUE)
  tolerance <- max(abs(eig$values)) * sqrt(.Machine$double.eps)
  linear <- abs(eig$values) <= tolerance
  if (any(eig$values < -tolerance) || sum(linear) != d || all(linear)) {
    stop("The roughness penalty does not identify the expected linear and curvature modes.", call. = FALSE)
  }
  spec$transform <- whiten %*% eig$vectors
  surface_penalty <- eig$values / min(eig$values[!linear])
  surface_penalty[linear] <- null_weight
  spec$surface_penalty <- surface_penalty
  sm$X <- NULL
  spec[["smooth"]] <- sm
  spec$rank <- ncol(x) + 1L
  spec <- .svem_spline_general_marginals(spec, marginal_modes, marginal_weight)
  .svem_spline_general_nominal_layout(spec, level_modes, level_penalty_ratio)
}

# Eigenvalues decrease, so the linear modes come last and the smoothest
# curvature modes immediately precede them.
.svem_spline_general_nominal_layout <- function(spec, level_modes, level_penalty_ratio) {
  r <- length(spec$surface_penalty)
  d <- spec$dimension
  active <- spec$nominal[lengths(spec$levels[spec$nominal]) > 1L]
  if (!d && !length(active)) spec$smoother <- "intercept limit"
  levels_count <- sum(lengths(spec$levels[active]))
  budget <- spec$max_columns
  if (r + levels_count > budget) {
    stop("Shared surface and nominal main effects exceed max_columns; increase that explicit resource budget or lower rank.", call. = FALSE)
  }
  # Reserve all main effects. Reduce whole interaction families, preserving
  # permutation symmetry across equally supported category labels and factors.
  linear_modes <- if (d) tail(seq_len(r), d) else integer()
  if (r + levels_count * (1L + length(linear_modes)) > budget) {
    linear_modes <- integer()
    spec$basis_notes <- c(spec$basis_notes, "Level-specific numeric slopes omitted by max_columns.")
  }
  used <- r + levels_count * (1L + length(linear_modes))
  cells <- list()
  for (ord in seq_len(spec$interaction_order)[-1L]) {
    if (length(active) < ord) next
    proposed <- utils::combn(active, ord, simplify = FALSE)
    need <- sum(vapply(proposed, function(p) prod(lengths(spec$levels[p])), numeric(1)))
    if (used + need <= budget) {cells <- c(cells, proposed); used <- used + need}
    else spec$basis_notes <- c(spec$basis_notes, paste0(ord, "-factor nominal cells omitted by max_columns."))
  }
  available_curvature <- if (!is.null(spec$nominal_curvature)) spec$nominal_curvature else head(seq_len(r), r - d)
  q <- min(level_modes, length(available_curvature), if (levels_count) floor((budget - used) / levels_count) else level_modes)
  curvature <- if (q > 0L) tail(available_curvature, q) else integer()
  if (q < min(level_modes, length(available_curvature))) spec$basis_notes <- c(spec$basis_notes, "Level curvature rank reduced by max_columns.")
  modes <- c(curvature, linear_modes)
  used <- used + levels_count * length(curvature)
  # A pair of nominal factors times a numeric linear direction is a third
  # order study interaction. Allocate the complete family when it fits.
  pair_cells <- cells[lengths(cells) == 2L]
  pair_modes <- integer()
  need <- d * sum(vapply(pair_cells, function(p) prod(lengths(spec$levels[p])), numeric(1)))
  if (spec$interaction_order >= 3L && length(pair_cells) && d) {
    if (used + need <= budget) pair_modes <- tail(seq_len(r), d)
    else spec$basis_notes <- c(spec$basis_notes, "Nominal-pair by numeric slopes omitted by max_columns.")
  }
  penalty <- spec$surface_penalty
  columns <- if (r) paste0("surface", seq_len(r)) else character()
  families <- if (!is.null(spec$numeric_families)) spec$numeric_families else c(rep("joint_curvature", r-d), rep("linear", d))
  for (v in active) {
    lev <- spec$levels[[v]]
    penalty <- c(penalty, rep(spec$null_weight, length(lev)))
    columns <- c(columns, paste0("nominal_", v, "_", lev))
    families <- c(families, rep("nominal_main", length(lev)))
    for (m in modes) {
      penalty <- c(penalty, rep(level_penalty_ratio * max(1, spec$surface_penalty[m]), length(lev)))
      columns <- c(columns, paste0("nominal_", v, "_", lev, "_mode", m))
      families <- c(families, rep("level_deviation", length(lev)))
    }
  }
  for (p in cells) {
    size <- prod(lengths(spec$levels[p]))
    label <- paste(p, collapse = ":")
    penalty <- c(penalty, rep(1, size))
    columns <- c(columns, paste0("cell_", label, "_", seq_len(size)))
    families <- c(families, rep("nominal_cell", size))
    if (length(p) == 2L) for (m in pair_modes) {
      penalty <- c(penalty, rep(level_penalty_ratio, size))
      columns <- c(columns, paste0("cell_", label, "_", seq_len(size), "_mode", m))
      families <- c(families, rep("nominal_numeric", size))
    }
  }
  spec$active_nominal <- active
  spec$cells <- cells
  spec$pairs <- pair_cells
  spec$pair_modes <- pair_modes
  spec$level_modes <- modes
  spec$level_penalty_ratio <- level_penalty_ratio
  spec$penalty <- penalty
  spec$columns <- columns
  spec$column_families <- families
  stopifnot(length(penalty) <= budget, all(penalty > 0))
  class(spec) <- "svem_spline_spec"
  spec
}

.svem_spline_general_legal <- function(spec, data) {
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(spec$predictors %in% names(data)) ||
      !all(vapply(data[spec$continuous], .svem_spline_numeric, logical(1))) ||
      !all(vapply(data[spec$nominal], .svem_spline_general_is_nominal, logical(1)))) {
    stop("Supply all named physical coordinates and nominal levels.", call. = FALSE)
  }
  tol <- if (is.null(spec$legal_tol)) .svem_spline_general_constants()$tol else spec$legal_tol
  x <- as.matrix(data[spec$continuous])
  ok <- rowSums(!is.finite(x) | sweep(x, 2L, spec$legal_lower, "-") < -tol |
                  sweep(x, 2L, spec$legal_upper, "-") > tol) == 0L
  # The saved model has no direction for a coordinate that was fixed on the grid.
  if (length(spec$fixed)) {
    off <- abs(sweep(x[, spec$fixed, drop = FALSE], 2L, spec$fixed_value, "-")) > tol
    ok <- ok & rowSums(off) == 0L
  }
  for (g in spec$mixture_groups) {
    a <- x[, g$vars, drop = FALSE]
    ok <- ok & abs(rowSums(a) - g$total) <= tol * max(1, abs(g$total)) &
      rowSums(sweep(a, 2L, g$lower, "-") < -tol | sweep(a, 2L, g$upper, "-") > tol) == 0L
  }
  for (v in spec$nominal) {
    ok <- ok & !is.na(match(as.character(data[[v]]), spec$levels[[v]]))
  }
  ok[is.na(ok)] <- FALSE
  unname(ok)
}

.svem_spline_general_coordinates <- function(spec, data) {
  if (!all(.svem_spline_general_legal(spec, data))) {
    stop("Recipe violates the declared spline bounds, mixture total, fixed reference coordinates or ",
         "nominal levels, or contains nonfinite settings.", call. = FALSE)
  }
  .svem_spline_general_chart_coordinates(spec, data)
}

.svem_spline_general_chart_coordinates <- function(spec, data) {
  x <- as.matrix(data[spec$continuous])
  # Fixed columns divide by a zero width here, and the subset below excludes them.
  if (!spec$dimension) return(as.data.frame(matrix(numeric(), nrow(data), 0L)))
  safe_widths <- spec$widths
  safe_widths[safe_widths == 0] <- 1
  u <- sweep(sweep(x, 2L, spec$lower, "-"), 2L, safe_widths, "/")
  z <- 2 * u[, spec$process, drop = FALSE] - 1
  for (m in spec$charts) {
    centered <- sweep(sweep(x[, m$vars, drop = FALSE], 2L, m$center, "-"), 2L, spec$widths[m$vars], "/")
    z <- cbind(z, centered %*% m$tangent / m$radius)
  }
  if (any(!is.finite(z))) stop("Physical coordinates exceed the numerical chart range.", call. = FALSE)
  colnames(z) <- spec$coordinates
  as.data.frame(z)
}

# Version 2 has no cached grid design. All prediction paths evaluate the same
# saved basis, so no batch can switch between a cached and a computed result.
.svem_spline_general_design <- function(spec, data) {
  empty <- matrix(numeric(), 0L, length(spec$penalty), dimnames = list(NULL, spec$columns))
  if (!nrow(data)) return(empty)
  z <- .svem_spline_general_coordinates(spec, data)
  if (!is.null(spec[["smooth"]]) && !requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the optional 'mgcv' package to evaluate this saved spline.", call. = FALSE)
  }
  phi <- if (!spec$dimension) matrix(numeric(), nrow(data), 0L) else .svem_spline_general_numeric(spec, data, z)
  n <- nrow(phi)
  parts <- list(phi)
  index <- list()
  for (v in spec$active_nominal) {
    i <- match(as.character(data[[v]]), spec$levels[[v]])
    index[[v]] <- i
    h <- matrix(0, n, length(spec$levels[[v]]))
    h[cbind(seq_len(n), i)] <- 1
    parts[[length(parts) + 1L]] <- h
    for (m in spec$level_modes) parts[[length(parts) + 1L]] <- h * phi[, m]
  }
  for (p in spec$cells) {
    sizes <- lengths(spec$levels[p])
    cell <- rep(0L, n)
    for (j in seq_along(p)) cell <- cell * sizes[j] + index[[p[j]]] - 1L
    h <- matrix(0, n, prod(sizes))
    h[cbind(seq_len(n), cell + 1L)] <- 1
    parts[[length(parts) + 1L]] <- h
    if (length(p) == 2L) for (m in spec$pair_modes) parts[[length(parts) + 1L]] <- h * phi[, m]
  }
  x <- do.call(cbind, parts)
  dimnames(x) <- list(NULL, spec$columns)
  x
}

# Measured nominal levels are per response, as measured block levels are.
.svem_spline_general_covered <- function(spec, measured_levels, data) {
  ok <- .svem_spline_general_legal(spec, data)
  for (v in spec$nominal) ok <- ok & as.character(data[[v]]) %in% measured_levels[[v]]
  ok
}

# Deterministic capacity report at equal weights. It consumes no random draws.
.svem_spline_general_capacity <- function(a, inverse_penalty, z, lambda, cap) {
  unit <- rep(1, length(z))
  member <- .svem_spline_member(a, inverse_penalty, z, unit, unit, lambda)
  data.frame(lambda = lambda, slope_edf = member$slope_edf,
             median_validation_kish_cap = cap,
             admissible_at_median_cap = member$slope_edf < cap)
}

.svem_spline_general_print_spec <- function(x, ...) {
  cat("General-factor ", x$smoother, " spline specification for ", x$response, "\n",
      "  ", x$dimension, " independent coordinates, surface rank ", x$rank, ", ",
      length(x$penalty), " model columns\n",
      "  Fixed coordinates: ", if (length(x$fixed)) paste(x$fixed, collapse = ", ") else "none", "\n",
      "  Nominal factors: ", if (length(x$nominal)) paste(x$nominal, collapse = ", ") else "none", "\n",
      if (length(x$active_nominal)) paste0("  Modelled: level shifts, ", length(x$level_modes),
        " level deviation modes, nominal cells up to order ", x$interaction_order, ".\n") else "",
      if (length(x$basis_notes)) paste0("  Basis notes: ", paste(x$basis_notes, collapse = " "), "\n") else "",
      sep = "")
  invisible(x)
}


.svem_spline_v1_eligible <- function(grid, groups) {
  grid <- unique(grid)
  if (!all(vapply(grid, .svem_spline_numeric, logical(1))) || length(groups) > 1L) return(FALSE)
  x <- as.matrix(grid)
  if (any(!is.finite(x))) return(TRUE)
  if (any(.svem_spline_general_fixed(apply(x, 2L, range)))) return(FALSE)
  d <- ncol(x) - length(groups)
  d %in% 1:3 && nrow(grid) >= d + 2L
}

.svem_spline_general_fixed <- function(chart) {
  width <- chart[2L, ] - chart[1L, ]
  scale <- pmax(abs(chart[1L, ]), abs(chart[2L, ]))
  width <= 8 * .Machine$double.eps * scale
}


# Physical-factor curvature complements the joint surface in dimensions >= 4.
# Its projection removes shared linear directions, not joint curvature; the
# positive diagonal penalties make the resulting overlapping basis identifiable.
.svem_spline_general_marginals <- function(spec, modes, weight) {
  spec$marginal_modes <- modes
  spec$marginal_weight <- weight
  if (spec$dimension < 4L || !modes) return(spec)
  r <- length(spec$surface_penalty); d <- spec$dimension
  ref <- spec$domain_grid[!duplicated(spec$domain_grid[setdiff(spec$continuous, spec$fixed)]), , drop = FALSE]
  z <- .svem_spline_general_chart_coordinates(spec, ref)
  phi <- sweep(mgcv::PredictMat(spec[["smooth"]], z), 2L, spec$center, "-") %*% spec$transform
  lin <- tail(seq_len(r), d); linear <- phi[, lin, drop = FALSE]
  vars <- setdiff(spec$continuous, spec$fixed)
  nominal_main <- sum(lengths(spec$levels)[lengths(spec$levels) > 1L])
  q <- min(modes, floor((spec$max_columns - r - nominal_main) / length(vars)))
  if (q < modes) spec$basis_notes <- c(spec$basis_notes, "Marginal curvature rank reduced by max_columns.")
  if (q < 1L) return(spec)
  .axis <- NULL # Captured for lookup in axis data by the smooth specification.
  marginal <- list(); penalty <- numeric()
  for (v in vars) {
    a <- data.frame(.axis = 2 * (ref[[v]] - spec$lower[[v]]) / spec$widths[[v]] - 1)
    k <- min(q + 2L, length(unique(a$.axis)))
    if (k < 3L) next
    sm <- .svem_spline_general_preserving_rng(mgcv::smoothCon(
      mgcv::s(.axis, bs = "tp", m = 2, k = k, xt = list(max.knots = 1000, seed = 17091)),
      data = a, knots = NULL, absorb.cons = TRUE, scale.penalty = FALSE)[[1L]])
    center <- colMeans(sm$X); x <- sweep(sm$X, 2L, center, "-")
    gram <- eigen(crossprod(x) / nrow(x), symmetric = TRUE)
    if (min(gram$values) <= max(gram$values) * max(dim(x)) * .Machine$double.eps)
      stop("The reference grid does not identify marginal curvature for ", v, ".", call. = FALSE)
    whiten <- sweep(gram$vectors, 2L, sqrt(gram$values), "/")
    rough <- crossprod(whiten, sm$S[[1L]] %*% whiten)
    eig <- eigen((rough + t(rough)) / 2, symmetric = TRUE)
    positive <- which(eig$values > max(abs(eig$values)) * sqrt(.Machine$double.eps))
    take <- tail(positive, min(q, length(positive)))
    if (!length(take)) next
    transform <- whiten %*% eig$vectors[, take, drop = FALSE]
    projection <- qr.solve(linear, x %*% transform)
    sm$X <- NULL
    marginal[[v]] <- list(var = v, lower = spec$lower[[v]], width = spec$widths[[v]],
      smooth = sm, center = center, transform = transform, projection = projection)
    penalty <- c(penalty, weight * eig$values[take] / min(eig$values[positive]))
  }
  if (!length(penalty)) return(spec)
  spec$marginal_specs <- marginal
  spec$shared_linear <- lin
  spec$nominal_curvature <- setdiff(seq_len(r), lin)
  spec$numeric_order <- c(spec$nominal_curvature, r + seq_along(penalty), lin)
  spec$surface_penalty <- c(spec$surface_penalty, penalty)[spec$numeric_order]
  spec$numeric_families <- c(rep("joint_curvature", r-d), rep("marginal_curvature", length(penalty)), rep("linear", d))
  spec$rank <- length(spec$surface_penalty) + 1L
  spec
}

.svem_spline_general_numeric <- function(spec, data, z) {
  mat <- if (is.null(spec[["smooth"]])) as.matrix(z) else mgcv::PredictMat(spec[["smooth"]], z)
  phi <- sweep(mat, 2L, spec$center, "-") %*% spec$transform
  if (is.null(spec$marginal_specs)) return(phi)
  linear <- phi[, spec$shared_linear, drop = FALSE]
  parts <- list(phi)
  for (a in spec$marginal_specs) {
    u <- data.frame(.axis = 2 * (data[[a$var]] - a$lower) / a$width - 1)
    curved <- sweep(mgcv::PredictMat(a$smooth, u), 2L, a$center, "-") %*% a$transform
    parts[[length(parts) + 1L]] <- curved - linear %*% a$projection
  }
  do.call(cbind, parts)[, spec$numeric_order, drop = FALSE]
}

# Shared dispatch defaults; signatures remain literal for readable help.
.svem_spline_general_defaults <- function() {
  list(max_columns = 200, level_modes = 3, level_penalty_ratio = 1,
       interaction_order = 3, marginal_modes = 4, marginal_weight = .25)
}
.svem_spline_nondefault <- function(control) {
  defaults <- .svem_spline_general_defaults()
  any(vapply(intersect(names(control), names(defaults)), function(v)
    !identical(as.numeric(control[[v]]), defaults[[v]]), logical(1)))
}

# Disclosure only: these notes never change scores or selection.
.svem_spline_general_notes <- function(fit) {
  spec <- fit$spec
  notes <- spec$basis_notes
  fraction <- mean(fit$member_diagnostics$intercept_only)
  if (fraction >= .5 && fraction < 1) notes <- c(notes, sprintf(
    "%.0f%% of members selected the intercept-only fit; the directional nomination comes from the remaining members.", 100 * fraction))
  vars <- setdiff(spec$continuous, spec$fixed)
  if (length(vars) && nrow(fit$training_data)) {
    x <- as.matrix(fit$training_data[vars])
    tol <- if (is.null(spec$legal_tol)) .svem_spline_general_constants()$tol else spec$legal_tol
    outside <- rowSums(sweep(x, 2L, spec$lower[vars], "-") < -tol |
      sweep(x, 2L, spec$lower[vars] + spec$widths[vars], "-") > tol) > 0L
    if (any(outside)) notes <- c(notes, sprintf(
      "%d of %d measured preparations lie outside the reference grid range.", sum(outside), length(outside)))
  }
  paste(notes, collapse = " ")
}
