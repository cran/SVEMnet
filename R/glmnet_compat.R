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
.svem_prepare_glmnet_design <- function(x, args = list(), exclude = NULL) {
  if (!is.matrix(x) || ncol(x) != 1L) {
    return(list(x = x, args = args, exclude = exclude, sentinel = NULL))
  }

  sentinel <- ".svem_internal_zero_sentinel"
  while (sentinel %in% colnames(x)) sentinel <- paste0(sentinel, "_")
  x_fit <- cbind(x, stats::setNames(list(rep(0, nrow(x))), sentinel)[[1L]])
  colnames(x_fit)[ncol(x_fit)] <- sentinel
  storage.mode(x_fit) <- "double"

  if (!is.null(args$penalty.factor)) {
    pf <- args$penalty.factor
    if (!is.numeric(pf) || length(pf) != 1L || !is.finite(pf) || pf <= 0) {
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

  list(x = x_fit, args = args, exclude = exclude, sentinel = sentinel)
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
                                      old_control = .glmnet_direct_control_args()) {
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
