# Solver invocation and coefficient extraction for SVEMnet.
# No selection criterion is evaluated here; solver diagnostics are retained.
.svem_glmnet_path <- function(X, y, weights, alpha, args, gamma,
                               sentinel = NULL) {
  messages <- errors <- character()
  fit <- tryCatch(withCallingHandlers(
    do.call(glmnet::glmnet, c(list(x = X, y = y, alpha = alpha,
                                   weights = weights), args)),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { errors <<- conditionMessage(e); NULL })
  paths <- list()
  if (!is.null(fit) && length(fit$lambda)) for (g in gamma) {
    cf <- tryCatch({
      if (isTRUE(args$relax)) as.matrix(stats::coef(fit, s = fit$lambda, gamma = g))
      else as.matrix(stats::coef(fit, s = fit$lambda))
    }, error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
    if (!is.null(cf)) paths[[length(paths) + 1L]] <- list(
      coefficients = .svem_strip_glmnet_sentinel(cf, sentinel),
      lambda = fit$lambda, gamma = g)
  }
  if (!is.null(fit) && !length(fit$lambda)) errors <- c(errors, "glmnet returned an empty lambda path")
  list(fit = fit, paths = paths, warnings = unique(messages), errors = unique(errors))
}

.svem_glmnet_setup <- function(X, dots, fam_name, relax_flag) {
  ## ---- capture and sanitize user '...' ----
  if (length(dots)) {
    if ("weights" %in% names(dots)) {
      warning("Ignoring user 'weights'; SVEM uses bootstrap weights.")
      dots$weights <- NULL
    }
    if ("standardize" %in% names(dots)) {
      warning("Ignoring user 'standardize'; SVEMnet always standardizes.")
      dots$standardize <- NULL
    }
    if ("offset" %in% names(dots)) {
      warning("Ignoring user 'offset'; SVEMnet does not support offsets.")
      dots$offset <- NULL
    }
    protected <- c("x","y","alpha","intercept","relax","standardize",
                   "nlambda","lambda.min.ratio","lambda","family")
    dropped_reserved <- intersect(names(dots), protected)
    if (length(dropped_reserved)) {
      warning("Ignoring reserved glmnet argument(s) managed by SVEMnet: ",
              paste(dropped_reserved, collapse = ", "), ".")
    }
    dots <- dots[setdiff(names(dots), protected)]

    # glmnet() silently swallows unknown arguments; catch misspellings here
    unknown <- .svem_unknown_glmnet_args(dots, funs = list(glmnet::glmnet),
                                         context = "glmnet via SVEMnet")
    if (length(unknown)) dots <- dots[setdiff(names(dots), unknown)]
  }

  exclude_user <- dots$exclude
  dots$exclude <- NULL
  glmnet_design <- .svem_prepare_glmnet_design(
    X, args = dots, exclude = exclude_user
  )
  X_glm <- glmnet_design$x
  dots <- glmnet_design$args

  glmnet_base_args <- .glmnet_prepare_call_args(
    c(list(intercept = TRUE,
           standardize = TRUE,
           nlambda = 500,
           relax = isTRUE(relax_flag),
           family = fam_name,
           exclude = glmnet_design$exclude),
      dots),
    glmnet::glmnet,
    default_control = list(maxit = 1e6)
  )

  list(design = glmnet_design, args = glmnet_base_args, dots = dots)
}
