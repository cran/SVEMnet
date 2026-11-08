#' Fit a Gaussian self-validated smooth surface
#'
#' Each member fits a training-weighted penalized mean and residual variance
#' and scores these on its validation-weighted copy. For training and
#' validation mean squared residuals `t` and `v`, the selection score is
#' `log(t) + v/t`. A training fit that understates residual variation can
#' receive a poor validation score. Predictions average the selected means.
#'
#' @param spec A saved [spline_terms()] specification.
#' @param data Measured rows containing the physical predictors, optional block
#'   IDs and a numeric response. Rows missing one of these fields are omitted.
#'   Study factors can be numeric or nominal as declared in the saved specification.
#' @param response A single response column name. Defaults to the specification
#'   response. Reuse the same physical specification for multiple endpoints.
#' @param nBoot Number of fractional random-weight members (default 200).
#' @param store_member_weights Store the exact member uniforms and mean-one
#'   training and validation weights. Use `set.seed()` for reproducibility.
#' @param lambda Advanced smoothing path starting with `Inf` for the
#'   intercept-only model, followed by strictly decreasing positive values.
#' @details One uniform draw per preparation gives antithetic exponential
#'   training and validation weights. Minimize weighted training mean squared
#'   error plus lambda times the quadratic penalty. Each member considers an
#'   intercept-only fit, by default followed by nine smoothing strengths decreasing from
#'   10 to 0.0001 on a log scale. The selected nonintercept effective degrees
#'   of freedom must be strictly below the validation Kish size bounded between
#'   2 and n. Both linear directions and measured additive block contrasts
#'   contribute to this bound. Each block contrast has penalty one.
#'
#'   A constant response returns that constant exactly. General fits also treat
#'   variation within eight scale-relative machine epsilons as rounding-only,
#'   returning a constant up to that rounding. Otherwise the score
#'   uses responses centered at their midrange and divided by their largest
#'   absolute centered value, with training variance floored at squared
#'   machine epsilon in those units. For positive variance above the floor
#'   this conditioning changes all candidate scores by one additive constant.
#'   The first minimum on the ordered admissible path wins exact ties.
#'
#'   General fits record per-response measured nominal levels and an equal-weight
#'   path-capacity diagnostic beside the actual member Kish caps. Their
#'   member diagnostics also report effective degrees of freedom by basis
#'   family. These add to slope EDF but are a descriptive allocation among
#'   overlapping penalized bases, not a scientific variance decomposition. A constant
#'   general surface has no directional optimum. When the joint fit is flat,
#'   portfolio scoring omits the spline nominee and reports `no_signal` for
#'   fully covered endpoints. Endpoints with missing predictions retain their
#'   coverage status. Ordinary predictions remain available on covered recipes.
#'   A constant endpoint can still contribute
#'   to a joint goal when other weighted endpoints vary.
#'
#' @return A `svem_spline` fit. Its `member_diagnostics` record smoothing
#'   strengths, conditional effective degrees of freedom, Kish caps, selected
#'   residual variances and Gaussian scores in conditioned response units.
#'   `path_scores` contains all scores with excluded points set to infinity.
#'   `coef_matrix` contains conditioned-response coefficients in the saved
#'   basis, one row per member. Ordinary predictions use original response
#'   units. Optional member spread describes variation among fitted members.
#' @seealso [spline_terms()]
#' @export
svem_spline <- function(spec, data, response = NULL, nBoot = 200L,
                        store_member_weights = FALSE,
                        lambda = c(Inf, 10^seq(1, -4, length.out = 9L))) {
  if (!inherits(spec, "svem_spline_spec") || !isTRUE(spec$version %in% c(1L, 2L))) {
    stop("Pass a supported spline_terms() specification.", call. = FALSE)
  }
  nBoot <- .svem_integer_scalar(nBoot, "nBoot", min = 1L)
  store_member_weights <- .svem_logical_scalar(store_member_weights, "store_member_weights")
  if (is.null(response)) response <- spec$response
  if (!is.character(response) || length(response) != 1L || is.na(response) ||
      !nzchar(response) || response %in% c(spec$predictors, spec$blocking)) {
    stop("response must name one numeric measured endpoint separate from the predictors.", call. = FALSE)
  }
  measured <- .svem_spline_measured(spec, data, response)
  y <- measured$y
  raw <- measured$x
  n <- length(y)
  if (n < 2L) stop("At least two complete measured preparations are required.", call. = FALSE)
  unreplicated <- spec$blocking[vapply(raw[spec$blocking], function(x) {
    .svem_spline_numeric(x) && !anyDuplicated(x)
  }, logical(1))]
  if (length(unreplicated)) warning("Treating unreplicated numeric block(s) as nominal IDs: ",
      paste(unreplicated, collapse = ", "), ". Continuous spline blocking trends are unsupported.", call. = FALSE)
  blocks <- lapply(raw[spec$blocking], function(x) {
    lev <- if (is.factor(x)) levels(droplevels(x)) else unique(x)
    size <- length(lev)
    contrasts <- if (size > 1L) qr.Q(qr(stats::contr.helmert(size))) * sqrt(size) else matrix(numeric(), 1L, 0L)
    list(levels = lev, contrasts = contrasts)
  })
  x <- .svem_spline_design(spec, raw, blocks)
  penalty <- c(spec$penalty, rep(1, ncol(x) - length(spec$penalty)))
  constant <- all(y == y[1L]) || (identical(spec$version, 2L) &&
    max(y) - min(y) <= 8 * .Machine$double.eps * max(abs(y)))
  center <- if (constant) y[1L] else min(y)/2 + max(y)/2
  centered <- y - center
  scale <- if (constant) 1 else max(abs(centered))
  if (!is.finite(scale) || scale <= 0) stop("Response scale is numerically unsupported.", call. = FALSE)
  z <- centered / scale
  inverse_penalty <- 1 / sqrt(penalty)
  a <- sweep(x, 2L, inverse_penalty, "*")
  lambda <- .svem_spline_check_lambda(lambda)
  groups <- if (identical(spec$version, 2L)) c(spec$column_families, rep("block", ncol(x)-length(spec$penalty))) else NULL
  coefficients <- matrix(0, nBoot, ncol(x) + 1L)
  colnames(coefficients) <- c("(Intercept)", colnames(x))
  scores <- matrix(NA_real_, nBoot, length(lambda))
  diagnostics <- vector("list", nBoot)
  if (store_member_weights) uniforms <- training <- validation <- matrix(0, nBoot, n)
  eps <- .Machine$double.eps
  for (i in seq_len(nBoot)) {
    u <- pmin(pmax(stats::runif(n), eps), 1 - eps)
    w <- -log(u)
    v <- -log1p(-u)
    w <- w * (n / sum(w))
    v <- v * (n / sum(v))
    if (store_member_weights) {
      uniforms[i, ] <- u
      training[i, ] <- w
      validation[i, ] <- v
    }
    member <- .svem_spline_member(a, inverse_penalty, z, w, v, lambda, groups = groups)
    index <- if (constant) 1L else which.min(member$score)
    coefficients[i, ] <- member$coefficients[, index]
    scores[i, ] <- member$score
    diagnostics[[i]] <- c(member = i, grid_index = index, lambda = lambda[index],
      effective_df = 1 + member$slope_edf[index], slope_edf = member$slope_edf[index],
      kish_validation_n_eff = member$kish, kish_validation_cap = member$cap,
      excluded_points = sum(member$slope_edf >= member$cap),
      intercept_only = as.numeric(index == 1L), training_mse = member$t[index],
      validation_mse = member$v[index], gaussian_likelihood = member$score[index])
    # Component EDF is descriptive: overlapping penalized basis families
    # do not constitute a unique scientific variance decomposition.
    if (length(groups)) diagnostics[[i]] <- c(diagnostics[[i]],
      stats::setNames(member$family_edf[, index], paste0("edf_", rownames(member$family_edf))))
  }
  details <- as.data.frame(do.call(rbind, diagnostics))
  details$intercept_only <- as.logical(details$intercept_only)
  rhs <- Reduce(function(a, b) call("+", a, b), lapply(c(spec$predictors, spec$blocking), as.name))
  formula <- stats::as.formula(call("~", as.name(response), rhs), env = baseenv())
  fit <- structure(list(formula = formula, family = "gaussian",
      response = response, spec = spec, actual_y = y, training_data = raw,
      training_rows = measured$rows, nBoot = nBoot, objective_used = "gaussian_likelihood",
      block_levels = lapply(blocks, `[[`, "levels"), blocks = blocks,
      coef_matrix = coefficients, mean_coefficients = colMeans(coefficients),
      response_center = center, response_scale = scale, lambda_grid = lambda,
      penalty = penalty, path_scores = scores, member_diagnostics = details,
      member_weights = if (store_member_weights) list(training = training,
        validation = validation, uniforms = uniforms) else NULL), class = "svem_spline")
  if (identical(spec$version, 2L)) {
    fit$nominal_levels <- lapply(raw[spec$nominal], function(x) unique(as.character(x)))
    fit$member_diagnostics$admissible_points <- length(lambda) - 1L - details$excluded_points
    fit$capacity <- .svem_spline_general_capacity(a, inverse_penalty, z, lambda,
      stats::median(details$kish_validation_cap))
  }
  fit
}

.svem_spline_measured <- function(spec, data, response) {
  vars <- c(spec$predictors, spec$blocking)
  fields <- c(vars, response)
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(fields %in% names(data)) ||
      !.svem_spline_numeric(data[[response]]) ||
      !all(vapply(data[setdiff(spec$predictors, spec$nominal)], .svem_spline_numeric, logical(1)))) {
    stop("Measured data must contain the numeric response and all declared predictor and block columns.", call. = FALSE)
  }
  for (v in spec$nominal) if (!.svem_spline_general_is_nominal(data[[v]])) {
    stop("Supply declared nominal predictors as factor, character or logical levels.", call. = FALSE)
  }
  for (b in spec$blocking) {
    x <- data[[b]]
    if (!(is.factor(x) || is.character(x) || .svem_spline_numeric(x))) {
      stop("A spline block must contain factor, character or numeric level identifiers.", call. = FALSE)
    }
  }
  rows <- which(stats::complete.cases(data[fields]))
  raw <- data[rows, vars, drop = FALSE]
  y <- as.numeric(data[[response]][rows])
  if (any(!is.finite(y)) || any(vapply(raw[spec$blocking], function(x) {
    is.numeric(x) && any(!is.finite(x))
  }, logical(1)))) stop("Complete measured rows contain a nonfinite response or block.", call. = FALSE)
  list(x = raw, y = y, rows = rows)
}

.svem_spline_design <- function(spec, data, blocks) {
  x <- .svem_spline_surface(spec, data)
  if (!identical(spec$version, 2L)) colnames(x) <- paste0("surface", seq_len(ncol(x)))
  for (b in names(blocks)) {
    if (!b %in% names(data)) stop("Missing declared spline block: ", b, call. = FALSE)
    value <- data[[b]]
    support <- blocks[[b]]$levels
    if (is.numeric(support)) {
      if (!.svem_spline_numeric(value)) stop("Expected numeric spline block IDs in ", b, ".", call. = FALSE)
      index <- match(value, support)
    } else {
      index <- match(as.character(value), support)
    }
    if (anyNA(index)) stop("Unknown or missing measured spline block level in ", b, ".", call. = FALSE)
    extra <- blocks[[b]]$contrasts[index, , drop = FALSE]
    if (ncol(extra)) {
      colnames(extra) <- paste0("block", match(b, names(blocks)), "_", seq_len(ncol(extra)))
      x <- cbind(x, extra)
    }
  }
  x
}

.svem_spline_check_lambda <- function(lambda) {
  lambda <- unname(lambda)
  if (!is.numeric(lambda) || !is.null(dim(lambda)) || length(lambda) < 2L ||
      !identical(as.numeric(lambda[1L]), Inf) || any(!is.finite(lambda[-1L])) ||
      any(lambda[-1L] <= 0) || any(diff(lambda) >= 0)) {
    stop("lambda must start at Inf and decrease strictly through positive finite values.", call. = FALSE)
  }
  lambda
}

.svem_spline_member <- function(a, inverse_penalty, y, w, v, lambda, groups = NULL) {
  wt <- w / sum(w)
  vv <- v / sum(v)
  if (!ncol(a)) {
    mean <- sum(wt * y)
    t <- sum(wt * (y - mean)^2)
    valid <- sum(vv * (y - mean)^2)
    kish <- sum(v)^2 / (sum(v^2) + .Machine$double.eps)
    out <- list(coefficients = matrix(mean, 1L, length(lambda)),
      score = rep(log(max(t, .Machine$double.eps^2)) + valid / max(t, .Machine$double.eps^2), length(lambda)),
      slope_edf = rep(0, length(lambda)), kish = kish, cap = max(2, min(length(y), kish)),
      t = rep(t, length(lambda)), v = rep(valid, length(lambda)))
    if (!is.null(groups)) out$family_edf <- matrix(numeric(), 0L, length(lambda))
    return(out)
  }
  xbar <- drop(crossprod(a, wt))
  ybar <- sum(wt * y)
  centered <- sweep(a, 2L, xbar, "-")
  gram <- eigen(crossprod(centered * sqrt(wt)), symmetric = TRUE)
  if (min(gram$values) < -max(1, max(gram$values)) * ncol(a) * 1e-12) {
    stop("The training smoother has a numerically indefinite Gram matrix.", call. = FALSE)
  }
  eigenvalues <- pmax(gram$values, 0)
  inverse <- 1 / outer(eigenvalues, lambda, "+")
  rotated_y <- drop(crossprod(gram$vectors, crossprod(centered, wt * (y - ybar))))
  theta <- gram$vectors %*% (rotated_y * inverse)
  intercept <- ybar - drop(crossprod(xbar, theta))
  residual <- (y - ybar) - centered %*% theta
  t <- drop(crossprod(wt, residual^2))
  valid <- drop(crossprod(vv, residual^2))
  safe <- pmax(t, .Machine$double.eps^2)
  score <- log(safe) + valid / safe
  slope_edf <- colSums(eigenvalues * inverse)
  kish <- sum(v)^2 / (sum(v^2) + .Machine$double.eps)
  cap <- max(2, min(length(y), kish))
  score[slope_edf >= cap] <- Inf
  if (anyNA(score) || !is.finite(score[1L])) stop("Invalid Gaussian smoothing score.", call. = FALSE)
  out <- list(coefficients = rbind(intercept, theta * inverse_penalty), score = score,
       slope_edf = slope_edf, kish = kish, cap = cap, t = t, v = valid)
  if (!is.null(groups)) {
    contribution <- gram$vectors^2 %*% (eigenvalues * inverse)
    out$family_edf <- rowsum(contribution, groups, reorder = FALSE)
  }
  out
}

#' Predict a self-validated Gaussian spline
#'
#' @param object A fitted [svem_spline()] object.
#' @param newdata Physical recipes with the named surface variables and all
#'   declared blocking factors at known measured levels. Nominal study factors
#'   must also be at levels measured for this response. Keep `mgcv` installed
#'   to predict at new sites after reloading a serialized fit.
#' @param type Use `"response"` for the ensemble mean in original units and
#'   `"members"` for the matrix of member predictions (rows are recipes).
#' @param se.fit Return the mean and the sample standard deviation across
#'   fitted members. This spread is a diagnostic, not a calibrated standard
#'   error or prediction interval.
#' @param interval Return empirical percentile bounds across fitted members.
#'   These describe member variation and have no predictive coverage guarantee.
#' @param level Central fraction of member predictions for percentile bounds.
#' @param debias Must be `FALSE`. Predictions use the original fitted means.
#' @param ... Reserved for extensions.
#' @return A numeric vector, a recipe-by-member matrix, or a list with
#'   `fit`, optionally `se.fit`, and optionally `lwr` and `upr`.
#' @export
predict.svem_spline <- function(object, newdata, type = c("response", "members"),
                                se.fit = FALSE, interval = FALSE, level = 0.95,
                                debias = FALSE, ...) {
  if (length(list(...))) stop("Unused spline prediction arguments.", call. = FALSE)
  type <- match.arg(type)
  se.fit <- .svem_logical_scalar(se.fit, "se.fit")
  interval <- .svem_logical_scalar(interval, "interval")
  level <- .svem_numeric_scalar(level, "level", lower = 0, upper = 1,
                                lower_open = TRUE, upper_open = TRUE)
  if (.svem_logical_scalar(debias, "debias")) stop("The spline uses its fitted means; debias must be FALSE.", call. = FALSE)
  if ((se.fit || interval) && (type == "members" || object$nBoot < 2L)) {
    stop("Spread and percentile bounds require response predictions and at least two fitted spline members.", call. = FALSE)
  }
  if (identical(object$spec$version, 2L) &&
      !all(.svem_spline_general_covered(object$spec, object$nominal_levels, newdata))) {
    stop("Unknown or unmeasured nominal level, or a recipe outside the saved spline domain.", call. = FALSE)
  }
  surface <- .svem_spline_design(object$spec, newdata, object$blocks)
  if (!nrow(newdata)) {
    if (type == "members") return(matrix(numeric(), 0L, object$nBoot))
    if (se.fit || interval) {
      answer <- list(fit = numeric())
      if (se.fit) answer$se.fit <- numeric()
      if (interval) answer[c("lwr", "upr")] <- list(numeric(), numeric())
      return(answer)
    }
    return(numeric())
  }
  x <- cbind("(Intercept)" = 1, surface)
  if (type == "response" && !se.fit && !interval) {
    return(as.numeric(object$response_center + object$response_scale * (x %*% object$mean_coefficients)))
  }
  m <- object$response_center + object$response_scale * tcrossprod(x, object$coef_matrix)
  if (type == "members") return(m)
  answer <- list(fit = rowMeans(m))
  if (se.fit) answer$se.fit <- apply(m, 1L, stats::sd)
  if (interval) {
    alpha <- (1 - level)/2
    q <- apply(m, 1L, stats::quantile, probs = c(alpha, 1 - alpha), names = FALSE)
    answer$lwr <- q[1L, ]
    answer$upr <- q[2L, ]
  }
  answer
}

#' @export
print.svem_spline <- function(x, ...) {
  d <- x$member_diagnostics
  cat("Gaussian self-validated ", if (identical(x$spec$version, 2L)) x$spec$smoother else "thin-plate",
      if (identical(x$spec$version, 2L)) " surface for " else " spline for ", x$response, "\n",
      "  ", nrow(x$training_data), " measured preparations, ", x$nBoot, " members\n",
      "  Surface rank ", x$spec$rank, ", median effective degrees of freedom ",
      format(stats::median(d$effective_df), digits = 4), "\n",
      "  ", sum(d$intercept_only), " intercept-only members, ",
      sum(d$excluded_points), " path points excluded by the validation Kish cap\n", sep = "")
  if (identical(x$spec$version, 2L)) {
    notes <- .svem_spline_general_notes(x)
    if (nzchar(notes)) cat("  ", notes, "\n", sep = "")
  }
  invisible(x)
}

.svem_spline_check_training <- function(fit, data) {
  rows <- .svem_spline_measured(fit$spec, data, fit$response)
  if (nrow(rows$x) != nrow(fit$training_data)) stop("The complete spline training rows differ.", call. = FALSE)
  keyed <- function(x, y) {
    for (b in c(fit$spec$blocking, fit$spec$nominal)) if (!is.numeric(x[[b]])) x[[b]] <- as.character(x[[b]])
    x <- cbind(as.data.frame(y), x)
    ord <- do.call(order, c(unname(as.list(x)), list(method = "radix")))
    lapply(x, function(v) if (is.numeric(v)) as.numeric(v[ord]) else as.character(v[ord]))
  }
  if (!identical(keyed(rows$x, rows$y), keyed(fit$training_data, fit$actual_y))) {
    stop("The measured outcomes or physical settings differ from the saved spline fit.", call. = FALSE)
  }
  invisible(TRUE)
}
