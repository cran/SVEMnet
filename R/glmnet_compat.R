.glmnet_direct_control_args <- function() {
  c("thresh", "maxit", "dfmax", "pmax", "trace.it")
}

# Variable names appearing in the response (LHS) of a two-sided formula.
# For a call LHS such as log(y), as.character(formula[[2L]]) returns a
# length-2 vector ("log", "y"), which breaks scalar contexts ('&&', 'if',
# '[[<-') on R >= 4.2; all.vars() is robust to any LHS expression.
.svem_response_vars <- function(formula) {
  tryCatch(all.vars(formula[[2L]]), error = function(e) character(0L))
}

# Single-string label for the response expression (e.g. "y" or "log(y)"),
# for display and column naming.
.svem_response_label <- function(formula) {
  tryCatch(paste(deparse(formula[[2L]]), collapse = " "),
           error = function(e) NA_character_)
}

# Internal selection uses the solver's exact active set. glmnet returns exact
# zeros at stored path knots; a raw-coefficient tolerance changes support when
# predictor or response units change.
.svem_solver_active <- function(beta) {
  slope <- beta[-1L]
  which(is.finite(slope) & slope != 0)
}

.svem_solver_support <- function(beta, count_intercept = TRUE) {
  if (!is.numeric(beta) || !length(beta)) return(NA_integer_)
  as.integer(length(.svem_solver_active(beta)) + as.integer(count_intercept))
}

.svem_solver_support_path <- function(coefficients, count_intercept = TRUE) {
  slopes <- coefficients[-1L, , drop = FALSE]
  as.integer(count_intercept) + colSums(is.finite(slopes) & slopes != 0)
}

.svem_reject_formula_offset <- function(terms, fun_name) {
  if (length(attr(terms, "offset"))) {
    stop(fun_name, " does not support formula offset() terms.", call. = FALSE)
  }
  invisible(NULL)
}

.svem_raw_predictor_vars <- function(terms, data_names) {
  intersect(all.vars(stats::delete.response(terms)), data_names)
}

# Plain formulas estimate only levels present among complete training rows.
# Frozen bigexp specifications deliberately bypass this preparation.
.svem_drop_unused_frame_levels <- function(mf) {
  vars <- names(mf)
  if (isTRUE(attr(attr(mf, "terms"), "response") == 1L)) vars <- vars[-1L]
  removed <- character()
  for (v in vars) {
    if (!is.factor(mf[[v]])) next
    value <- droplevels(mf[[v]])
    absent <- setdiff(levels(mf[[v]]), levels(value))
    if (length(absent)) {
      if (!is.null(attr(mf[[v]], "contrasts"))) {
        stop("Cannot remove unobserved levels from predictor '", v,
             "' with explicit contrasts. Drop unused levels before defining ",
             "the contrasts, or use a locked bigexp specification.", call. = FALSE)
      }
      removed <- c(removed, paste0(v, " (", paste(absent, collapse = ", "), ")"))
      mf[[v]] <- value
    }
  }
  if (length(removed)) {
    warning("Removed unobserved categorical level(s) from formula predictors: ",
            paste(removed, collapse = "; "), ".", call. = FALSE)
  }
  mf
}

.svem_drop_unused_raw_levels <- function(data, predictor_vars) {
  for (v in intersect(predictor_vars, names(data))) {
    if (is.factor(data[[v]])) data[[v]] <- droplevels(data[[v]])
  }
  data
}

# Model-frame levels also cover expressions such as factor(x) and paste(a, b).
# Keep these distinct from raw-variable levels and sampling supports.
.svem_model_frame_schema <- function(mf) {
  tt <- attr(mf, "terms")
  vars <- names(mf)
  if (isTRUE(attr(tt, "response") == 1L)) vars <- vars[-1L]
  vars <- vars[vapply(mf[vars], function(x) is.factor(x) || is.character(x), logical(1))]
  levels <- lapply(mf[vars], function(x) {
    if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
  })
  list(levels = levels,
       ordered = vars[vapply(mf[vars], is.ordered, logical(1))])
}

.svem_harmonize_raw_levels <- function(object, newdata) {
  xlev <- object$xlevels
  if (is.null(xlev) || !is.list(xlev)) return(newdata)
  classes <- .svem_training_var_classes(object)
  for (v in intersect(names(xlev), names(newdata))) {
    ordered <- !isTRUE(object$used_bigexp_spec) &&
      (is.ordered(newdata[[v]]) || identical(unname(classes[v]), "ordered"))
    newdata[[v]] <- factor(as.character(newdata[[v]]), levels = xlev[[v]],
                           ordered = ordered,
                           exclude = if (anyNA(xlev[[v]])) NULL else NA)
  }
  newdata
}

.svem_harmonize_model_frame <- function(object, mf) {
  saved <- object$schema$model_frame
  # Old serialized objects have no evaluated-frame metadata.
  if (is.null(saved) || is.null(saved$levels)) return(mf)
  for (v in intersect(names(saved$levels), names(mf))) {
    mf[[v]] <- factor(as.character(mf[[v]]), levels = saved$levels[[v]],
                       ordered = v %in% saved$ordered,
                       exclude = if (anyNA(saved$levels[[v]])) NULL else NA)
  }
  mf
}

# Only direct categorical conversions establish discrete numeric support.
# A general transform such as cut(x, ...) does not make x itself discrete.
.svem_base_factor_call <- function(expr, env) {
  if (!is.call(expr)) return(NULL)
  head <- expr[[1L]]
  name <- if (is.symbol(head)) as.character(head) else ""
  qualified <- is.call(head) && length(head) == 3L &&
    as.character(head[[1L]]) %in% c("::", ":::") &&
    identical(as.character(head[[2L]]), "base")
  if (qualified) name <- as.character(head[[3L]])
  if (!name %in% c("factor", "as.factor")) return(NULL)
  base_fun <- get(name, baseenv())
  if (!qualified && !identical(get0(name, envir = env, mode = "function",
                                    inherits = TRUE), base_fun)) return(NULL)
  list(name = name, call = match.call(base_fun, expr))
}

.svem_freeze_factor_predvars <- function(terms, data) {
  env <- attr(terms, ".Environment")
  walk <- function(expr) {
    if (!is.call(expr)) return(expr)
    conversion <- .svem_base_factor_call(expr, env)
    if (!is.null(conversion)) {
      matched <- conversion$call
      if (is.symbol(matched$x) && as.character(matched$x) %in% names(data)) {
        if (identical(conversion$name, "factor") &&
            !"levels" %in% names(matched)) {
          raw_call <- matched
          raw_call[[1L]] <- quote(base::factor)
          raw_call$labels <- NULL
          raw_levels <- levels(eval(raw_call, data, env))
          expr <- matched
          expr[[1L]] <- quote(base::factor)
          expr$levels <- raw_levels
        } else if (identical(conversion$name, "as.factor")) {
          value <- eval(matched, data, env)
          expr <- as.call(list(quote(base::factor), x = matched$x,
                              levels = levels(value), ordered = is.ordered(value),
                              exclude = if (anyNA(levels(value))) NULL else NA))
        }
      }
    }
    # relevel(g, ref) also needs an omitted reference level to exist while
    # the call is evaluated; the saved final frame still rejects that level.
    head <- expr[[1L]]
    is_relevel <- identical(head, as.name("relevel")) &&
      identical(get0("relevel", envir = env, mode = "function", inherits = TRUE),
                stats::relevel)
    qualified_relevel <- is.call(head) && length(head) == 3L &&
      as.character(head[[1L]]) %in% c("::", ":::") &&
      identical(as.character(head[[2L]]), "stats") &&
      identical(as.character(head[[3L]]), "relevel")
    if (is_relevel || qualified_relevel) {
      matched <- match.call(stats::relevel, expr)
      if (is.symbol(matched$x) && as.character(matched$x) %in% names(data)) {
        raw <- data[[as.character(matched$x)]]
        if (identical(class(raw), "factor")) {
          expr <- matched
          expr$x <- as.call(list(quote(base::factor), x = matched$x,
                                 levels = levels(raw),
                                 exclude = if (anyNA(levels(raw))) NULL else NA))
        }
      }
    }
    for (j in seq_along(expr)[-1L]) expr[j] <- list(walk(expr[[j]]))
    expr
  }
  pv <- attr(terms, "predvars")
  if (!is.null(pv)) attr(terms, "predvars") <- walk(pv)
  terms
}

.svem_direct_factor_support <- function(terms, data) {
  found <- list()
  env <- attr(terms, ".Environment")
  walk <- function(expr) {
    if (!is.call(expr)) return(invisible(NULL))
    conversion <- .svem_base_factor_call(expr, env)
    if (!is.null(conversion) && is.symbol(conversion$call$x)) {
      v <- as.character(conversion$call$x)
      if (v %in% names(data) && is.numeric(data[[v]])) {
        converted <- eval(conversion$call, data, env)
        support <- sort(unique(data[[v]][is.finite(data[[v]]) & !is.na(converted)]))
        found[[v]] <<- if (is.null(found[[v]])) support else intersect(found[[v]], support)
      }
    }
    for (j in seq_along(expr)[-1L]) walk(expr[[j]])
    invisible(NULL)
  }
  variables <- as.list(attr(stats::delete.response(terms), "variables"))[-1L]
  for (expr in variables) walk(expr)
  found
}

.svem_numerically_constant_columns <- function(X) {
  if (!ncol(X)) return(integer())
  center <- colMeans(X)
  variance <- colMeans(sweep(X, 2L, center, "-")^2)
  neps <- nrow(X) * .Machine$double.eps
  bound <- neps * variance + (neps * abs(center))^2
  unname(which(is.finite(variance) & variance <= bound))
}

.svem_union_exclude <- function(exclude, fixed, p) {
  if (!length(fixed)) return(exclude)
  normalize <- function(value) {
    if (is.logical(value)) {
      if (length(value) != p || anyNA(value)) {
        stop("Logical `exclude` must have one nonmissing value per design column.", call. = FALSE)
      }
      which(value)
    } else value
  }
  if (is.function(exclude)) {
    user_exclude <- exclude
    return(function(...) unique(c(normalize(user_exclude(...)), fixed)))
  }
  unique(c(normalize(exclude), fixed))
}

.svem_fit_calibration <- function(y, y_pred, X, beta) {
  center <- mean(y_pred)
  spread <- max(abs(y_pred - center))
  m_eps <- length(beta) * .Machine$double.eps
  gamma <- m_eps / (1 - m_eps)
  magnitude <- abs(beta[1L]) + as.vector(abs(X) %*% abs(beta[-1L]))
  if (!is.finite(spread) || spread <= 2 * gamma * max(magnitude)) {
    return(list(fit = NULL, reason = "prediction spread is at its rounding-error bound"))
  }
  form <- y ~ y_pred
  environment(form) <- baseenv()
  fit <- stats::lm(form, data = data.frame(y = y, y_pred = y_pred))
  if (fit$rank == 2L && all(is.finite(stats::coef(fit)))) {
    return(list(fit = fit, reason = NULL))
  }
  form <- stats::as.formula(substitute(y ~ I((y_pred - CENTER) / SCALE),
                                       list(CENTER = center, SCALE = spread)),
                            env = baseenv())
  fit <- stats::lm(form, data = data.frame(y = y, y_pred = y_pred))
  if (fit$rank != 2L || any(!is.finite(stats::coef(fit)))) {
    return(list(fit = NULL, reason = "calibration has no finite rank-two fit"))
  }
  attr(fit, "svem_calibration_transform") <- c(center = center, scale = spread)
  list(fit = fit, reason = NULL)
}

.svem_calibration_coefficients <- function(fit) {
  ab <- stats::coef(fit)
  transform <- attr(fit, "svem_calibration_transform")
  if (!is.null(transform)) {
    ab <- c(ab[1L] - ab[2L] * transform[["center"]] / transform[["scale"]],
            ab[2L] / transform[["scale"]])
  }
  ab
}

.svem_apply_calibration <- function(fit, values) {
  ab <- stats::coef(fit)
  if (length(ab) < 2L || any(!is.finite(ab[1:2]))) return(values)
  transform <- attr(fit, "svem_calibration_transform")
  if (is.null(transform)) return(ab[1L] + ab[2L] * values)
  ab[1L] + ab[2L] * ((values - transform[["center"]]) / transform[["scale"]])
}

# Build a small, self-contained evaluation environment for stored formulas and
# terms objects.  The parent supplies the standard formula helpers (poly,
# contrasts, etc.); only objects genuinely supplied by the formula environment
# are copied.  Local functions are cloned with similarly compact environments
# so serializing a fitted model does not also serialize the caller's data frame.
.svem_compact_formula_env <- function(formula, data_names = character(0L)) {
  source_env <- environment(formula)
  if (is.null(source_env)) source_env <- parent.frame()
  target_env <- new.env(parent = asNamespace("stats"))

  assigned_names <- function(expr) {
    out <- character(0L)
    walk <- function(e) {
      if (!is.call(e)) return(invisible(NULL))
      head <- as.character(e[[1L]])[1L]
      if (head %in% c("<-", "<<-", "=", "->", "->>")) {
        lhs <- if (head %in% c("->", "->>")) e[[3L]] else e[[2L]]
        if (is.symbol(lhs)) out <<- c(out, as.character(lhs))
      } else if (identical(head, "for") && length(e) >= 2L && is.symbol(e[[2L]])) {
        out <<- c(out, as.character(e[[2L]]))
      }
      for (j in seq_along(e)[-1L]) walk(e[[j]])
      invisible(NULL)
    }
    walk(expr)
    unique(out)
  }

  copy_binding <- function(name, from, into, depth = 0L) {
    if (depth > 12L || exists(name, envir = into, inherits = FALSE) ||
        !exists(name, envir = from, inherits = TRUE)) {
      return(invisible(NULL))
    }
    value <- get(name, envir = from, inherits = TRUE)
    inherited <- tryCatch(get(name, envir = parent.env(into), inherits = TRUE),
                          error = function(e) NULL)
    if (!is.null(inherited) && identical(value, inherited)) {
      return(invisible(NULL))
    }

    if (is.function(value) && !is.primitive(value)) {
      old_env <- environment(value)
      if (!is.null(old_env) && !isNamespace(old_env) &&
          !identical(old_env, baseenv()) && !identical(old_env, emptyenv())) {
        fn_env <- new.env(parent = into)
        cloned <- value
        environment(cloned) <- fn_env
        assign(name, cloned, envir = into)

        locals <- unique(c(names(formals(value)), assigned_names(body(value))))
        formal_expr <- as.call(c(quote(list), formals(value)))
        deps <- unique(c(
          all.names(body(value), functions = FALSE, unique = TRUE),
          all.names(body(value), functions = TRUE, unique = TRUE),
          all.names(formal_expr, functions = FALSE, unique = TRUE),
          all.names(formal_expr, functions = TRUE, unique = TRUE)
        ))
        deps <- setdiff(deps, c(locals, "...", "", "function", "{") )
        for (dep in deps) copy_binding(dep, old_env, fn_env, depth + 1L)
        return(invisible(NULL))
      }
    }
    assign(name, value, envir = into)
    invisible(NULL)
  }

  expr <- formula
  symbols <- unique(c(
    all.names(expr, functions = FALSE, unique = TRUE),
    all.names(expr, functions = TRUE, unique = TRUE)
  ))
  symbols <- setdiff(symbols, c(data_names, "", ".", "~", "+", "-", "*", "/",
                                "^", ":", "::", ":::", "(", "[", "[[",
                                "$", "@", "=", "<-", "|"))
  for (name in symbols) copy_binding(name, source_env, target_env)
  target_env
}

.svem_compact_formula_terms <- function(formula, terms, data_names = character(0L)) {
  env <- .svem_compact_formula_env(formula, data_names = data_names)
  environment(formula) <- env
  environment(terms) <- env
  list(formula = formula, terms = terms)
}

# glmnet 4.x rejects a one-column x matrix.  Add an excluded, identically-zero
# column only for the call into glmnet; user-visible designs and coefficients
# never contain it.  penalty.factor is normalized so the real feature retains
# the same effective penalty it has in a genuine one-feature fit.
.svem_glmnet_callback_exclude <- function(callback, p, has_sentinel) {
  force(callback); force(p); force(has_sentinel)
  function(x, y, weights, ...) {
    callback_x <- if (has_sentinel) x[, seq_len(p), drop = FALSE] else x
    value <- callback(x = callback_x, y = y, weights = weights, ...)
    if (is.null(value)) return(integer())
    message <- paste0("`exclude` callback must return indices in 1..", p,
                      " or a logical vector of length ", p, ".")
    if (is.logical(value)) {
      if (length(value) != p || anyNA(value)) stop(message, call. = FALSE)
      return(which(value))
    }
    if (!is.numeric(value) || is.complex(value) || any(!is.finite(value)) ||
        any(value != floor(value)) || any(value < 1 | value > p)) {
      stop(message, call. = FALSE)
    }
    as.integer(value)
  }
}

.svem_glmnet_callback_design <- function(design, callback, p) {
  pf <- design$args$penalty.factor
  if (is.null(pf)) pf <- rep(1, ncol(design$x))
  if (!is.numeric(pf) || is.complex(pf) || length(pf) != ncol(design$x) || anyNA(pf)) {
    stop("For an `exclude` callback, `penalty.factor` must be numeric, have ",
         "one value per design column, and contain no missing values.", call. = FALSE)
  }
  has_sentinel <- !is.null(design$sentinel)
  fixed <- c(design$numerically_constant_columns,
             if (has_sentinel) ncol(design$x) else integer())
  # glmnet checks callback exclusions before adding infinite-penalty exclusions,
  # then resets every excluded penalty to one before normalization. This keeps
  # our fixed exclusions out of its callback-only minimum-column check.
  pf[fixed] <- Inf
  design$args$penalty.factor <- pf
  design$exclude <- .svem_glmnet_callback_exclude(callback, p, has_sentinel)
  design
}

.svem_prepare_glmnet_design <- function(x, args = list(), exclude = NULL,
                                        numerically_constant = NULL) {
  if (is.null(numerically_constant)) {
    numerically_constant <- .svem_numerically_constant_columns(x)
  }
  callback <- NULL
  if (is.function(exclude) && (length(numerically_constant) || ncol(x) == 1L)) {
    callback <- exclude
    exclude <- NULL
  } else {
    exclude <- .svem_union_exclude(exclude, numerically_constant, ncol(x))
  }
  if (!is.matrix(x) || ncol(x) != 1L) {
    design <- list(x = x, args = args, exclude = exclude, sentinel = NULL,
                   numerically_constant_columns = numerically_constant)
    if (!is.null(callback)) design <- .svem_glmnet_callback_design(design, callback, ncol(x))
    return(design)
  }

  sentinel <- ".svem_internal_zero_sentinel"
  while (sentinel %in% colnames(x)) sentinel <- paste0(sentinel, "_")
  x_fit <- cbind(x, stats::setNames(list(rep(0, nrow(x))), sentinel)[[1L]])
  colnames(x_fit)[ncol(x_fit)] <- sentinel
  storage.mode(x_fit) <- "double"

  if (!is.null(args$penalty.factor)) {
    pf <- args$penalty.factor
    if (!is.numeric(pf) || is.complex(pf) || length(pf) != 1L || !is.finite(pf) || pf <= 0) {
      stop("For a one-predictor fit, `penalty.factor` must be one positive finite number.")
    }
    args$penalty.factor <- c(1, 1)
  }

  sentinel_index <- ncol(x_fit)
  if (is.null(exclude)) {
    exclude <- sentinel_index
  } else if (is.function(exclude)) {
    user_exclude <- exclude
    exclude <- function(...) unique(c(user_exclude(...), sentinel_index))
  } else if (is.logical(exclude)) {
    if (length(exclude) != 1L || is.na(exclude)) {
      stop("For a one-predictor fit, logical `exclude` must have length one.")
    }
    exclude <- unique(c(if (exclude) 1L else integer(0L), sentinel_index))
  } else if (is.numeric(exclude)) {
    if (any(!is.finite(exclude)) || any(exclude != floor(exclude)) ||
        any(exclude < 1L) || any(exclude > 1L)) {
      stop(
        "For a one-predictor fit, numeric `exclude` must contain only the integer index 1.",
        call. = FALSE
      )
    }
    exclude <- unique(c(as.integer(exclude), sentinel_index))
  } else {
    stop("`exclude` must be NULL, numeric, logical, or a function.")
  }

  design <- list(x = x_fit, args = args, exclude = exclude, sentinel = sentinel,
                 numerically_constant_columns = numerically_constant)
  if (!is.null(callback)) design <- .svem_glmnet_callback_design(design, callback, ncol(x))
  design
}

.svem_strip_glmnet_sentinel <- function(coef_mat, sentinel = NULL) {
  if (is.null(sentinel) || is.null(rownames(coef_mat)) ||
      !sentinel %in% rownames(coef_mat)) return(coef_mat)
  coef_mat[setdiff(rownames(coef_mat), sentinel), , drop = FALSE]
}

.svem_training_var_classes <- function(object) {
  out <- tryCatch(object$sampling_schema$var_classes, error = function(e) NULL)
  if (is.null(out)) out <- tryCatch(object$schema$var_classes, error = function(e) NULL)
  if (is.null(out)) character(0L) else out
}

# Enforce the raw-variable type contract before model.frame() performs implicit
# coercions.  Integer and double predictors are intentionally interchangeable;
# character input is accepted for a predictor trained as a factor because it
# can be mapped safely to the saved training levels.
.svem_validate_newdata_classes <- function(object, newdata) {
  classes <- .svem_training_var_classes(object)
  if (!length(classes)) return(newdata)

  required <- names(classes)[!is.na(classes) & nzchar(classes)]
  missing <- setdiff(required, names(newdata))
  if (length(missing)) {
    stop(
      "`newdata` is missing required predictor(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  for (name in names(classes)) {
    expected <- classes[[name]]
    if (is.na(expected) || !nzchar(expected) || !name %in% names(newdata)) next
    value <- newdata[[name]]
    ok <- switch(
      expected,
      numeric = is.numeric(value),
      integer = is.numeric(value),
      factor = (is.factor(value) && !is.ordered(value)) || is.character(value),
      ordered = is.ordered(value) || is.character(value),
      character = is.character(value) || is.factor(value),
      logical = is.logical(value),
      inherits(value, expected)
    )
    # bigexp_prepare() intentionally turns categorical inputs into plain
    # factors. Accept both the original raw class and that prepared class.
    categorical <- name %in% names(object$sampling_schema$factor_levels)
    if (isTRUE(object$used_bigexp_spec) && categorical &&
        expected %in% c("factor", "ordered", "character", "logical")) {
      ok <- ok || (is.factor(value) && !is.ordered(value)) || is.character(value)
    }
    if (!isTRUE(ok)) {
      stop(
        "Predictor '", name, "' was fitted as ", expected,
        " but `newdata` supplies class ", paste(class(value), collapse = "/"),
        ". Integer and double are interchangeable; categorical predictors ",
        "may be supplied as factor or character.",
        call. = FALSE
      )
    }
  }
  newdata
}

# Names accepted inside control = list(...): the legacy direct-control set plus
# whatever the installed glmnet's glmnet.control() understands.
.svem_known_control_keys <- function() {
  ctl <- tryCatch(names(formals(glmnet::glmnet.control)),
                  error = function(e) character(0))
  unique(c(.glmnet_direct_control_args(), ctl))
}

# Warn about (and return) argument names that none of the target glmnet
# functions understand. glmnet() itself silently swallows unknown arguments
# through its '...', so misspelled arguments would otherwise vanish without
# any feedback.
.svem_unknown_glmnet_args <- function(dots,
                                      funs = list(glmnet::glmnet),
                                      extra_allowed = character(0),
                                      context = "glmnet",
                                      warn = TRUE) {
  nms <- names(dots)
  if (is.null(nms) || !length(nms)) return(character(0))
  nms <- nms[nzchar(nms)]
  if (!length(nms)) return(character(0))

  allowed <- unique(c(
    unlist(lapply(funs, function(f) names(formals(f)))),
    "control",
    .svem_known_control_keys(),
    extra_allowed
  ))
  unknown <- setdiff(nms, allowed)
  if (length(unknown) && isTRUE(warn)) {
    warning(
      "Ignoring unrecognized argument(s) passed to ", context, ": ",
      paste(unknown, collapse = ", "),
      ". Check for misspellings; glmnet silently ignores unknown arguments.",
      call. = FALSE
    )
  }
  unknown
}

.glmnet_supports_control <- function(fun) {
  "control" %in% names(formals(fun))
}

.glmnet_check_control_list <- function(control) {
  if (is.null(control)) return(control)
  if (!is.list(control)) {
    stop("'control' must be a named list of glmnet algorithm-control values.")
  }
  nms <- names(control)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("'control' must be a named list of glmnet algorithm-control values.")
  }
  control
}

.glmnet_prepare_call_args <- function(args, fun, default_control = list(),
                                      direct_control = .glmnet_direct_control_args(),
                                      old_control = .glmnet_direct_control_args(),
                                      scoped_legacy = FALSE) {
  if (!length(args)) args <- list()

  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- rep("", length(args))

  user_control <- NULL
  if ("control" %in% arg_names) {
    user_control <- .glmnet_check_control_list(args[["control"]])
    args[["control"]] <- NULL
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
  }

  direct_names <- intersect(direct_control, arg_names)
  direct_values <- args[direct_names]
  args <- args[setdiff(seq_along(args), match(direct_names, arg_names))]
  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- rep("", length(args))

  # glmnet.control()-only keys (epsnr, mxitnr, fdev, ...) passed directly:
  # glmnet()'s '...' would swallow them silently, so either route them via the
  # control list (new glmnet) or warn that they cannot take effect (old glmnet).
  ctl_only_names  <- intersect(setdiff(.svem_known_control_keys(), direct_control), arg_names)
  ctl_only_values <- args[ctl_only_names]
  if (length(ctl_only_names)) {
    args <- args[setdiff(seq_along(args), match(ctl_only_names, arg_names))]
  }

  if (.glmnet_supports_control(fun)) {
    control <- default_control
    if (!is.null(user_control)) {
      control <- utils::modifyList(control, user_control, keep.null = TRUE)
    }
    if (length(ctl_only_values)) {
      control <- utils::modifyList(control, ctl_only_values, keep.null = TRUE)
    }
    if (length(direct_values)) {
      control <- utils::modifyList(control, direct_values, keep.null = TRUE)
    }
    if (length(control)) args$control <- control
    return(args)
  }

  if (isTRUE(scoped_legacy)) {
    control <- default_control
    for (overrides in list(user_control, ctl_only_values, direct_values))
      if (length(overrides)) control <- utils::modifyList(control, overrides, keep.null = TRUE)
    global_keys <- setdiff(names(formals(glmnet::glmnet.control)), c("factory", old_control))
    unsupported <- setdiff(names(control), c(old_control, global_keys))
    if (length(unsupported)) stop("Unsupported glmnet control entr",
      if (length(unsupported) == 1L) "y: " else "ies: ",
      paste(unsupported, collapse = ", "), ".", call. = FALSE)
    for (nm in intersect(names(control), old_control)) args[[nm]] <- control[[nm]]
    global <- control[intersect(names(control), global_keys)]
    if (length(global)) attr(args, "svem_legacy_control") <- global
    return(args)
  }

  if (length(ctl_only_values)) {
    warning(
      "The installed glmnet cannot route these algorithm-control argument(s) ",
      "per call; they are ignored: ", paste(ctl_only_names, collapse = ", "),
      ". Set them globally with glmnet::glmnet.control() if needed.",
      call. = FALSE
    )
  }

  translated <- list()
  if (!is.null(user_control)) {
    extra <- setdiff(names(user_control), old_control)
    if (length(extra)) {
      stop(
        "The installed glmnet does not support control = list(...). ",
        "Only these control entries can be translated for older glmnet: ",
        paste(old_control, collapse = ", "), ". Unsupported entr",
        if (length(extra) == 1L) "y: " else "ies: ",
        paste(extra, collapse = ", "), "."
      )
    }
    translated <- user_control[intersect(names(user_control), old_control)]
  }

  old_values <- default_control
  if (length(translated)) {
    old_values <- utils::modifyList(old_values, translated, keep.null = TRUE)
  }
  if (length(direct_values)) {
    old_values <- utils::modifyList(old_values, direct_values, keep.null = TRUE)
  }

  if (length(old_values)) {
    for (nm in names(old_values)) {
      if (!nm %in% names(args)) args[[nm]] <- old_values[[nm]]
    }
  }

  args
}
