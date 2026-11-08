# Optional spline integration. These helpers do not modify primary SVEM fits.
.svem_spline_spec_from_bigexp <- function(spec, domain_grid, response, rank = NULL,
    null_weight = 0.05, level_modes = 3L, level_penalty_ratio = 1,
    max_columns = 200L, interaction_order = 3L, marginal_modes = 4L, marginal_weight = 0.25, general = TRUE) {
  if (!inherits(spec, "bigexp_spec")) stop("spec must come from bigexp_terms().", call. = FALSE)
  blocks <- spec$settings$blocking
  surface <- setdiff(spec$vars, blocks)
  if (!length(surface) || !is.data.frame(domain_grid) || !all(surface %in% names(domain_grid))) {
    stop("Supply the physical science factors and a declared feasible domain grid.", call. = FALSE)
  }
  groups <- spec$settings$mixture_groups
  if (!general && !.svem_spline_original_domain(spec, domain_grid, list())) {
    stop("The optional general spline is fitted when explicitly requested with spline = TRUE or n_optima = 5; this factor space is outside the original automatic numeric-spline domain.", call. = FALSE)
  }
  nominal <- surface[spec$is_cat[surface]]
  numeric_surface <- setdiff(surface, nominal)
  grid <- domain_grid
  for (v in nominal) grid[[v]] <- factor(as.character(grid[[v]]), levels = spec$levels[[v]])
  for (v in intersect(numeric_surface, names(spec$settings$discrete_levels))) {
    if (any(!grid[[v]] %in% spec$settings$discrete_levels[[v]])) {
      stop("Spline reference grid violates the declared discrete numeric support of ", v, ".", call. = FALSE)
    }
  }
  rhs <- Reduce(function(a, b) call("+", a, b), lapply(c(surface, blocks), as.name))
  form <- stats::as.formula(call("~", as.name(response), rhs), env = baseenv())
  bounds <- spec$num_range[, numeric_surface, drop = FALSE]
  for (g in .svem_mixture_groups(groups)) {
    bounds[1L, g$vars] <- g$lower
    bounds[2L, g$vars] <- g$upper
  }
  spline_terms(form, domain_grid = grid, blocking = blocks, mixture_groups = groups,
    bounds = bounds, rank = rank, null_weight = null_weight, level_modes = level_modes,
    level_penalty_ratio = level_penalty_ratio, max_columns = max_columns,
    interaction_order = interaction_order, marginal_modes = marginal_modes, marginal_weight = marginal_weight)
}

.svem_spline_nominal_labels <- function(data, nominal) {
  for (v in intersect(nominal, names(data))) {
    if (.svem_spline_numeric(data[[v]])) data[[v]] <- as.character(data[[v]])
  }
  data
}

.svem_spline_model_list <- function(objects, responses, data = NULL, check_training = FALSE,
                                    labels = responses) {
  if (!is.list(objects) || is.null(names(objects)) || anyDuplicated(names(objects)) ||
      !setequal(names(objects), responses)) {
    stop("spline must be a named list of fits for exactly the named response goals.", call. = FALSE)
  }
  objects <- objects[responses]
  for (i in seq_along(responses)) {
    r <- responses[i]
    obj <- objects[[r]]
    if (!inherits(obj, "svem_spline") && !inherits(obj, "svem_spline_unavailable")) {
      stop("spline must contain svem_spline() fits.", call. = FALSE)
    }
    if (!identical(.svem_response_label(obj$formula), unname(labels[i]))) {
      stop("Spline response labels must match their named goal columns; do not relabel fitted responses.", call. = FALSE)
    }
    if (check_training && inherits(obj, "svem_spline_unavailable")) {
      stop("Cannot reuse a previously unavailable spline for '", r,
           "'. Refit on the current measured data or set spline = TRUE to retry. ",
           "For an unsupported study, set spline = FALSE and reduce n_optima if necessary.", call. = FALSE)
    }
    if (check_training && inherits(obj, "svem_spline")) {
      measured <- .svem_spline_nominal_labels(data, obj$spec$nominal)
      tryCatch(.svem_spline_check_training(obj, measured), error = function(e) {
        stop("Cannot reuse supplied spline model for '", r, "': ", conditionMessage(e),
             " Refit on current measured data.", call. = FALSE)
      })
    }
  }
  objects
}

.svem_spline_check_schema <- function(fits, objects) {
  for (r in names(objects)) {
    fit <- fits[[r]]
    if (!inherits(fit, "svem_spline")) next
    primary <- objects[[r]]
    if (length(attr(primary$companions, "input_map")$continuous_blocks)) {
      stop("A spline with block IDs cannot follow declared continuous companion block covariates.", call. = FALSE)
    }
    if (!is.null(primary$family) && !identical(primary$family, "gaussian")) {
      stop("A Gaussian spline nominee requires Gaussian primary response fits.", call. = FALSE)
    }
    schema <- primary$sampling_schema
    if (!all(setdiff(fit$spec$predictors, fit$spec$nominal) %in% colnames(schema$num_ranges)) ||
        !all(fit$spec$nominal %in% names(schema$factor_levels))) {
      stop("Spline surface predictors must match the numeric and nominal physical factors in the primary model schema.", call. = FALSE)
    }
    if (!setequal(c(fit$spec$predictors, fit$spec$blocking), schema$predictor_vars) ||
        !setequal(fit$spec$blocking, schema$blocking)) {
      stop("Spline and primary fits must use the same raw science and blocking factors; a conditional slice cannot replace a full-domain model.", call. = FALSE)
    }
    if (!isTRUE(all.equal(.svem_mixture_groups(fit$spec$mixture_groups),
                          .svem_mixture_groups(schema$mixture_groups), check.attributes = FALSE))) {
      stop("Spline and primary fits must use the same mixture constraints.", call. = FALSE)
    }
  }
  invisible(NULL)
}

.svem_spline_basis_controls <- c("rank", "null_weight", "level_modes", "level_penalty_ratio", "max_columns", "interaction_order", "marginal_modes", "marginal_weight")
.svem_spline_fit_controls <- c("nBoot", "store_member_weights", "lambda")

.svem_spline_check_control <- function(control, dimension = NULL, original = FALSE) {
  allowed <- c(.svem_spline_basis_controls, .svem_spline_fit_controls)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || anyDuplicated(names(control)) || !all(names(control) %in% allowed)))) {
    stop("Spline controls must contain distinct named spline fitting or basis controls: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  if ("nBoot" %in% names(control)) .svem_integer_scalar(control$nBoot, "nBoot", min = 1L)
  if ("store_member_weights" %in% names(control))
    .svem_logical_scalar(control$store_member_weights, "store_member_weights")
  if ("lambda" %in% names(control)) .svem_spline_check_lambda(control$lambda)
  if ("rank" %in% names(control)) {
    minimum <- if (length(dimension) == 1L && is.finite(dimension) && dimension >= 1L) dimension + (if (isTRUE(original)) 2L else 1L) else 1L
    .svem_integer_scalar(control$rank, "rank", min = minimum)
  }
  if ("null_weight" %in% names(control))
    .svem_numeric_scalar(control$null_weight, "null_weight", lower = 0, lower_open = TRUE)
  if ("level_modes" %in% names(control)) .svem_integer_scalar(control$level_modes, "level_modes", min = 0L)
  if ("level_penalty_ratio" %in% names(control)) .svem_numeric_scalar(control$level_penalty_ratio, "level_penalty_ratio", lower = 0, lower_open = TRUE)
  if ("max_columns" %in% names(control)) .svem_integer_scalar(control$max_columns, "max_columns", min = 1L)
  if ("interaction_order" %in% names(control)) .svem_integer_scalar(control$interaction_order, "interaction_order", min = 1L, max = 3L)
  if ("marginal_modes" %in% names(control)) .svem_integer_scalar(control$marginal_modes, "marginal_modes", min = 0L)
  if ("marginal_weight" %in% names(control)) .svem_numeric_scalar(control$marginal_weight, "marginal_weight", lower = 0, lower_open = TRUE)
  invisible(control)
}

.svem_spline_dimension <- function(spec) {
  surface <- setdiff(spec$vars, spec$settings$blocking)
  numeric <- surface[!spec$is_cat[surface]]
  varying <- numeric[spec$num_range[2L, numeric] > spec$num_range[1L, numeric]]
  groups <- .svem_mixture_groups(spec$settings$mixture_groups)
  as.integer(length(varying) - sum(vapply(groups, function(g) any(g$vars %in% varying), logical(1))))
}

.svem_spline_basis <- function(spec, domain_grid, response, control = list(), general = TRUE) {
  geom <- intersect(names(control), .svem_spline_basis_controls)
  tryCatch(do.call(.svem_spline_spec_from_bigexp,
    c(list(spec = spec, domain_grid = domain_grid, response = response, general = general), control[geom])),
    error = identity)
}

.svem_fit_spline_responses <- function(spec, data, responses, domain_grid, control = list(),
                                      allow_unavailable = FALSE, basis = NULL,
                                      warn_unavailable = TRUE) {
  .svem_spline_check_control(control, .svem_spline_dimension(spec),
    original = .svem_spline_original_domain(spec, domain_grid, control))
  geom <- intersect(names(control), .svem_spline_basis_controls)
  if (is.null(basis)) basis <- .svem_spline_basis(spec, domain_grid, responses[1L], control)
  if (!inherits(basis, "error")) data <- .svem_spline_nominal_labels(data, basis$nominal)
  fit_one <- function(r) {
    if (inherits(basis, "error")) stop(basis)
    do.call(svem_spline, c(list(spec = basis, data = data, response = r),
                          control[setdiff(names(control), geom)]))
  }
  result <- lapply(responses, function(r) {
    if (!allow_unavailable) return(fit_one(r))
    tryCatch(fit_one(r), error = function(e) {
      if (warn_unavailable) warning("Spline candidate unavailable for '", r, "': ", conditionMessage(e),
              ". See the spline engine status and allocation counts.", call. = FALSE)
      structure(list(formula = stats::as.formula(call("~", as.name(r), 1), env = baseenv()),
                     message = conditionMessage(e), status = "fit_failed"),
                class = "svem_spline_unavailable")
    })
  })
  stats::setNames(result, responses)
}

.svem_spline_scores <- function(tab, x, fits, parameters, des, weights, combine, floor) {
  contributions <- list()
  statuses <- list()
  for (r in names(weights)) {
    obj <- fits[[r]]
    pr <- if (inherits(obj, "svem_spline_unavailable")) {
      simpleError(obj$message)
    } else tryCatch({
      legal <- if (identical(obj$spec$version, 2L))
        .svem_spline_general_covered(obj$spec, obj$nominal_levels, x) else .svem_spline_legal(obj$spec, x)
      for (b in names(obj$blocks)) {
        if (!b %in% names(x)) stop("Missing declared spline block: ", b, call. = FALSE)
        support <- obj$blocks[[b]]$levels
        value <- x[[b]]
        if (is.numeric(support) && !.svem_spline_numeric(value)) {
          stop("Expected numeric spline block IDs in ", b, ".", call. = FALSE)
        }
        index <- match(if (is.numeric(support)) value else as.character(value), support)
        legal <- legal & !is.na(index)
      }
      out <- rep(NA_real_, nrow(x))
      if (any(legal)) out[legal] <- as.numeric(stats::predict(obj, newdata = x[legal, , drop = FALSE]))
      out
    }, error = identity)
    state <- if (inherits(pr, "error")) "prediction_failed" else "ok"
    msg <- if (inherits(pr, "error")) conditionMessage(pr) else ""
    if (inherits(pr, "error")) pr <- rep(NA_real_, nrow(x))
    if (length(pr) != nrow(x)) stop("Spline returned an invalid number of predictions.", call. = FALSE)
    if (inherits(obj, "svem_spline_unavailable")) state <- "fit_failed"
    if (state == "ok") {
      uncovered <- sum(!is.finite(pr))
      if (uncovered == length(pr)) state <- "no_finite_predictions"
      else if (uncovered > 0L) state <- "partial_coverage"
      if (uncovered > 0L) {
        msg <- paste0(uncovered, if (identical(obj$spec$version, 2L))
          " recipe(s) have no prediction within the saved domain, measured nominal levels and measured block levels." else
          " recipe(s) have no prediction within the spline's domain and measured block levels.")
      }
    }
    if (state %in% c("ok", "partial_coverage") && identical(obj$spec$version, 2L))
      msg <- trimws(paste(msg, .svem_spline_general_notes(obj)))
    tab[[paste0(r, "_spline_pred")]] <- pr
    statuses[[r]] <- data.frame(engine = "spline", response = r, status = state, message = msg)
    if (weights[[r]] > 0) {
      d <- des(pr, parameters[[r]])
      d[!is.finite(pr)] <- NA_real_
      if (combine == "geom") d <- log(floor + (1 - floor) * d)
      contributions[[r]] <- weights[[r]] * d
    }
  }
  score <- Reduce("+", contributions)
  if (combine == "geom") score <- exp(score)
  # Flat general surfaces provide no direction for a new optimum. A constant
  # endpoint still contributes to a joint goal when another endpoint varies.
  positive <- names(weights)[weights > 0]
  flat <- vapply(fits[positive], function(f) inherits(f, "svem_spline") &&
    all(f$mean_coefficients[-1L] == 0), logical(1))
  no_signal <- all(flat) && any(vapply(fits[positive], function(f) identical(f$spec$version, 2L), logical(1)))
  if (no_signal) {
    score[] <- NA_real_
    for (r in positive) {
      if (statuses[[r]]$status == "ok") statuses[[r]]$status <- "no_signal"
      statuses[[r]]$message <- trimws(paste(statuses[[r]]$message,
        "The fitted spline surface is constant; no directional optimum is available."))
    }
  }
  tab$spline_score <- score
  unavailable <- vapply(fits[names(weights)[weights > 0]], inherits,
                        logical(1), what = "svem_spline_unavailable")
  if (!any(is.finite(score)) && !any(unavailable) && !no_signal) {
    warning("Spline fits cannot score the joint goal; inspect engine_status.", call. = FALSE)
  }
  list(table = tab, status = do.call(rbind, statuses))
}

.svem_resolve_plan_priority <- function(priority, use_default, spline, n_optima, n_explore) {
  explicit_counts <- !is.null(n_optima) || !is.null(n_explore)
  if (use_default && !identical(spline, FALSE) && (explicit_counts || isTRUE(spline) || is.list(spline))) {
    priority <- c(priority, "spline")
  }
  priority
}


# Match the public constructor's original/general dispatch when validating rank
# before any optional or primary fit. A lower-rank general linear limit remains
# admissible; the original numeric path retains its former rank requirement.
.svem_spline_original_domain <- function(spec, grid, control) {
  surface <- setdiff(spec$vars, spec$settings$blocking)
  if (!is.data.frame(grid) || !all(surface %in% names(grid)) || any(spec$is_cat[surface])) return(FALSE)
  if (.svem_spline_nondefault(control)) return(FALSE)
  .svem_spline_v1_eligible(grid[surface], .svem_mixture_groups(spec$settings$mixture_groups))
}
