# Planner-only validation: no fitting, prediction, or random draws.
.svem_plan_preflight <- function(spec, data, candidates, incumbent, reference,
                                 score_control, benchmark = NULL, budget = NULL,
                                 priority = NULL, n_optima = NULL, n_explore = NULL) {
  if (!is.data.frame(candidates) || !nrow(candidates) ||
      anyDuplicated(names(candidates))) {
    stop("candidates must be a nonempty data frame with unique column names.", call. = FALSE)
  }
  pool <- tryCatch(bigexp_prepare(spec, candidates, unseen = "error")$data,
    error = function(e) stop("Candidate pool is incompatible with `spec`: ",
      conditionMessage(e), ". For an unmeasured future day/operator, score at a known ",
      "reference level and record the actual execution context separately; rebuild ",
      "the specification and refit after new levels have been measured.", call. = FALSE))
  vars <- spec$vars
  if (anyNA(pool[vars]) || any(vapply(pool[vars], function(z) {
    is.numeric(z) && any(!is.finite(z))
  }, logical(1)))) stop("Candidate predictor settings must be finite and nonmissing.", call. = FALSE)
  blocks <- spec$settings$blocking
  for (b in blocks) {
    if (length(unique(pool[[b]])) != 1L) {
      stop("Supply one planned level of each blocking factor per worklist.", call. = FALSE)
    }
    if (!any(.svem_matches_block(data[[b]], pool[[b]][1L]))) {
      stop("Candidate blocking level in '", b, "' has no measured preparations. ",
           "For a future day/operator, score at a known reference level and record ",
           "the actual execution context separately; rebuild the specification and ",
           "refit after the new level has been measured.", call. = FALSE)
    }
  }
  .svem_pool_index(pool[vars], candidates, score_control$increments, blocks)
  declared_groups <- .svem_mixture_groups(spec$settings$mixture_groups)
  groups <- .svem_mixture_groups(score_control$mixture_groups)
  if (!is.null(declared_groups) && !is.null(score_control$mixture_groups) &&
      !identical(groups, declared_groups)) {
    stop("score_control mixture_groups must agree with the declared specification.", call. = FALSE)
  }
  if (is.null(groups)) groups <- declared_groups
  .svem_check_candidate_mixtures(pool, groups, 1e-8)
  benchmark_row <- .svem_resolve_benchmark(benchmark, data)
  check_budget <- function(rows) {
    carry <- !identical(incumbent, FALSE) && "incumbent" %in% priority
    if (!is.null(n_optima) || !is.null(n_explore)) {
      required_min <- as.integer(carry) + as.integer(!is.null(benchmark_row))
      if (carry && !is.null(benchmark_row) && length(rows)) {
        design <- setdiff(vars, blocks)
        refkeys <- .svem_recipe_keys(data[c(benchmark_row, rows), design, drop = FALSE],
                                      design, score_control$increments)
        if (!anyNA(refkeys) && refkeys[1L] %in% refkeys[-1L]) required_min <- 1L
      }
      .svem_count_targets(budget, n_optima, n_explore, seq_len(required_min), priority)
    }
    if (carry && !is.null(benchmark_row) && !is.null(budget) && budget < 2L) {
      design <- setdiff(vars, blocks)
      keys <- .svem_recipe_keys(data[c(benchmark_row, rows), design, drop = FALSE],
                                design, score_control$increments)
      if (!anyNA(keys) && !keys[1L] %in% keys[-1L]) {
        stop("budget is too small for the distinct required incumbent and fixed benchmark; allow at least two preparations.", call. = FALSE)
      }
    }
  }
  if (identical(incumbent, FALSE)) {
    check_budget(integer())
    return(invisible(NULL))
  }
  explicit <- is.numeric(incumbent) && length(incumbent) == 1L
  if (explicit) {
    row <- .svem_integer_scalar(incumbent, "incumbent", min = 1L, max = nrow(data))
    check_budget(row)
    if (!is.null(reference)) stop("Use either a historical row index or reference block values, not both.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!isTRUE(incumbent)) stop("incumbent must be TRUE, FALSE, or a historical row index.", call. = FALSE)
  if (is.null(reference)) reference <- as.list(pool[1L, blocks, drop = FALSE])
  if (!is.list(reference) || (length(reference) && (is.null(names(reference)) ||
      anyNA(names(reference)) || anyDuplicated(names(reference)))) ||
      !setequal(names(reference), blocks)) {
    stop("reference must name exactly the blocking factors, each with one historical value.", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(data))
  for (b in names(reference)) {
    value <- reference[[b]]
    if (length(value) != 1L || is.na(value)) stop("Each reference block value must be a single nonmissing value.", call. = FALSE)
    keep <- keep & .svem_matches_block(data[[b]], value)
  }
  if (!any(keep)) {
    stop("No historical preparations in the reference block context. Supply an incumbent row or set incumbent = FALSE.", call. = FALSE)
  }
  check_budget(which(keep))
  invisible(NULL)
}

# Shared with the scorer: resolution is pure and keeps the caller's stable
# identifier selector separate from its current row position.
.svem_resolve_benchmark <- function(benchmark, data) {
  if (is.null(benchmark)) return(NULL)
  if (!is.data.frame(data) || !nrow(data)) {
    stop("A fixed benchmark requires nonempty measured data.", call. = FALSE)
  }
  if (is.numeric(benchmark) && length(benchmark) == 1L) {
    return(.svem_integer_scalar(benchmark, "benchmark", min = 1L, max = nrow(data)))
  }
  if (!is.list(benchmark) || length(benchmark) != 1L || is.null(names(benchmark)) ||
      is.na(names(benchmark)) || !nzchar(names(benchmark)) ||
      !names(benchmark) %in% names(data)) {
    stop("benchmark must be one historical row number or a named list with one identifier column and value.", call. = FALSE)
  }
  value <- benchmark[[1L]]
  if (!is.atomic(value) || length(value) != 1L || !is.null(dim(value)) || is.na(value)) {
    stop("The benchmark identifier must be one nonmissing value.", call. = FALSE)
  }
  values <- data[[names(benchmark)]]
  rows <- which(!is.na(values) & as.character(values) == as.character(value))
  if (length(rows) != 1L) {
    stop("The benchmark identifier must match exactly one measured preparation; found ",
         length(rows), ". Use a unique preparation identifier.", call. = FALSE)
  }
  as.integer(rows)
}

.svem_plan_model_list <- function(objects, responses, label) {
  if (!is.list(objects) || is.null(names(objects)) || anyNA(names(objects)) ||
      anyDuplicated(names(objects)) || !setequal(names(objects), responses) ||
      !all(vapply(objects, inherits, logical(1), what = "svem_model"))) {
    stop(label, " must name one fitted model for each goal response.", call. = FALSE)
  }
  objects <- objects[responses]
  if (identical(label, "fits") &&
      any(vapply(objects, inherits, logical(1), what = c("svem_forward", "svem_forward_ic")))) {
    stop("fits is reserved for primary SVEM models. Use the forward argument for svem_forward() fits. ",
         "Single forward information-criterion fits are unsupported in this planner.", call. = FALSE)
  }
  if (identical(label, "fits") && any(vapply(objects, inherits, logical(1), what = "svem_cv"))) {
    stop("fits is reserved for primary self-validated SVEM models. A cross-validated fit cannot supply the primary SVEM nominee.", call. = FALSE)
  }
  if (identical(label, "forward") &&
      !all(vapply(objects, inherits, logical(1), what = "svem_forward"))) {
    stop("forward must contain Gaussian svem_forward() fits, not other SVEM engines.", call. = FALSE)
  }
  labels <- vapply(objects, function(o) .svem_response_label(o$formula), character(1))
  if (!identical(unname(labels), responses)) {
    stop("Supplied model response labels must match the named goal columns; do not relabel or swap fitted responses.", call. = FALSE)
  }
  if (any(vapply(objects, function(o) !is.null(o$family) && o$family != "gaussian", logical(1)))) {
    stop("The default formulation workflow requires Gaussian response fits.", call. = FALSE)
  }
  objects
}

.svem_plan_training_rows <- function(object, data) {
  if (!is.matrix(object$training_X) || !is.numeric(object$actual_y) ||
      is.null(object$formula) || is.null(object$terms) ||
      length(object$actual_y) != nrow(object$training_X) ||
      (ncol(object$training_X) && is.null(colnames(object$training_X)))) {
    stop("Stored training rows are unavailable or incomplete.")
  }
  data <- .svem_validate_newdata_classes(object, data)
  response <- eval(object$formula[[2L]], data, environment(object$formula))
  raw_vars <- intersect(all.vars(stats::delete.response(object$terms)), names(data))
  potentially_complete <- !is.na(response) & stats::complete.cases(data[raw_vars])
  # Never turn newly measured complete rows into omissions merely because an
  # old supplied model lacks their factor level. Missing endpoint rows remain
  # eligible for the endpoint-specific omission convention.
  for (v in intersect(names(object$xlevels), names(data))) {
    unknown <- !is.na(data[[v]]) & !as.character(data[[v]]) %in% object$xlevels[[v]]
    if (any(unknown & potentially_complete)) {
      stop("New measured factor levels are absent from the stored training levels in '", v, "'.")
    }
  }
  # A direct factor(x) in a plain formula can map new numeric settings to NA.
  # BigExp discrete supports restrict planned preparations, not measured rows.
  support <- if (!isTRUE(object$used_bigexp_spec)) .svem_discrete_map(object$sampling_schema)
  for (v in intersect(names(support), names(data))) {
    unknown <- !is.na(data[[v]]) & !data[[v]] %in% support[[v]]
    if (any(unknown & potentially_complete)) {
      stop("New measured discrete predictor settings are absent from the stored training support in '", v, "'.")
    }
  }
  data <- .svem_harmonize_raw_levels(object, data)
  tt <- object$terms
  if (identical(environment(tt), baseenv())) environment(tt) <- asNamespace("stats")
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass)
  mf <- .svem_harmonize_model_frame(object, mf)
  mf <- mf[stats::complete.cases(mf), , drop = FALSE]
  ctr <- object$contrasts
  if (is.null(ctr)) {
    saved <- object$schema$contrasts_options
    if (length(object$xlevels) && is.null(saved)) stop("Stored contrast settings are unavailable.")
    old <- options("contrasts")
    on.exit(options(old), add = TRUE)
    if (!is.null(saved)) options(contrasts = saved)
  } else {
    if (!is.list(ctr)) ctr <- as.list(ctr)
    ctr <- lapply(ctr, function(z) if (is.character(z)) get(z, mode = "function") else z)
  }
  X <- stats::model.matrix(stats::delete.response(tt), data = mf, contrasts.arg = ctr)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  y <- stats::model.response(mf)
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || any(!is.finite(X))) {
    stop("Measured data do not reproduce finite training responses and predictors.")
  }
  old_X <- object$training_X
  if (!setequal(colnames(X), colnames(old_X)) || nrow(X) != nrow(old_X)) {
    stop("The complete training rows or design columns differ.")
  }
  X <- X[, colnames(old_X), drop = FALSE]
  old_y <- as.numeric(object$actual_y)
  if (any(!is.finite(old_X)) || any(!is.finite(old_y))) stop("Stored training values are not finite.")
  # Sorting paired response/design rows preserves multiplicity and allows a
  # harmless data-row permutation. Response values must agree exactly; only
  # model-matrix reconstruction receives a machine-roundoff allowance.
  rows <- function(y, X) {
    a <- cbind(y, X)
    a[do.call(order, c(unname(as.data.frame(a)), list(method = "radix"))), , drop = FALSE]
  }
  now <- rows(as.numeric(y), X)
  prior <- rows(old_y, old_X)
  if (!identical(as.numeric(now[, 1L]), as.numeric(prior[, 1L])) ||
      any(abs(now[, -1L, drop = FALSE] - prior[, -1L, drop = FALSE]) >
          64 * .Machine$double.eps * pmax(1, abs(now[, -1L, drop = FALSE]),
                                         abs(prior[, -1L, drop = FALSE])))) {
    stop("The measured outcomes or design settings differ from the stored training rows.")
  }
  invisible(NULL)
}

.svem_plan_check_training <- function(objects, data, label, spec) {
  nominal <- spec$vars[spec$is_cat[spec$vars]]
  groups <- .svem_mixture_groups(spec$settings$mixture_groups)
  discrete <- spec$settings$discrete_levels
  for (response in names(objects)) {
    object <- objects[[response]]
    tryCatch({
      schema <- object$sampling_schema
      if (!setequal(schema$predictor_vars, spec$vars) ||
          !setequal(schema$blocking, spec$settings$blocking) ||
          !setequal(names(schema$factor_levels), nominal)) {
        stop("The fitted predictors, nominal factors or blocking roles differ from the declared specification.")
      }
      fitted_groups <- .svem_mixture_groups(schema$mixture_groups)
      if (!is.null(groups) && (!is.null(fitted_groups) || length(object$companions)) &&
          !identical(fitted_groups, groups)) {
        stop("The fitted mixture constraints differ from the declared specification.")
      }
      if (length(discrete)) {
        fitted_discrete <- .svem_discrete_map(schema)
        for (v in names(discrete)) {
          if (!setequal(fitted_discrete[[v]], discrete[[v]])) {
            stop("The fitted discrete support for '", v, "' differs from the declared specification.")
          }
        }
      }
      .svem_plan_training_rows(object, data)
    }, error = function(e) {
      stop("Cannot reuse supplied ", label, " model for '", response, "': ",
           conditionMessage(e), " Refit with the current specification and measured data or omit `", label,
           "` so the planner can fit it. Older objects without verifiable training ",
           "rows cannot be reused by this planner.", call. = FALSE)
    })
  }
  invisible(NULL)
}

.svem_plan_nominations <- function(scores) {
  nominations <- scores$nominations
  expected <- c("incumbent", "svemnet", "gaupro", "extra_trees", "forward")
  if ("spline" %in% scores$engine_status$engine || "spline" %in% nominations$engine) expected <- c(expected, "spline")
  if (!is.null(scores$benchmark)) expected <- c("benchmark", expected)
  for (engine in setdiff(expected, nominations$engine)) {
    row <- nominations[NA_integer_, , drop = FALSE]
    row$engine <- engine
    row$status <- "not_requested"
    nominations <- rbind(nominations, row)
  }
  nominations$fit_status <- vapply(nominations$engine, function(engine) {
    if (engine %in% c("incumbent", "benchmark")) return("observed")
    states <- scores$engine_status[scores$engine_status$engine == engine, , drop = FALSE]
    if (is.null(states) || !nrow(states)) return("not_requested")
    if (all(states$status == "ok")) return("ok")
    paste(paste0(states$response, ": ", states$status), collapse = "; ")
  }, character(1))
  rownames(nominations) <- NULL
  nominations[match(expected, nominations$engine), , drop = FALSE]
}
