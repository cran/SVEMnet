# Solver invocation and coefficient extraction for SVEMnet.
# No selection criterion is evaluated here; solver diagnostics are retained.
.svem_glmnet_path <- function(X, y, weights, alpha, args, gamma,
                               sentinel = NULL) {
  messages <- errors <- character()
  fit <- tryCatch(withCallingHandlers(
    .svem_glmnet_call(X, y, weights, alpha, args),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { errors <<- conditionMessage(e); NULL })
  paths <- list()
  base <- if (!is.null(fit) && length(fit$lambda))
    tryCatch(.svem_glmnet_coefficient_matrix(fit),
      error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL }) else NULL
  if (!is.null(base)) for (g in gamma) {
    cf <- tryCatch(.svem_glmnet_blend(fit, base, if (isTRUE(args$relax)) g else 1),
      error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
    if (!is.null(cf)) paths[[length(paths) + 1L]] <- list(
      coefficients = .svem_strip_glmnet_sentinel(cf, sentinel),
      lambda = fit$lambda, gamma = g)
  }
  if (!is.null(fit) && !length(fit$lambda)) errors <- c(errors, "glmnet returned an empty lambda path")
  list(fit = fit, paths = paths, warnings = unique(messages), errors = unique(errors))
}

.svem_glmnet_setup <- function(X, dots, fam_name, relax_flag,
                               numerically_constant = NULL) {
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
                                         extra_allowed = if (isTRUE(relax_flag)) c("maxp", "path") else character(),
                                         context = "glmnet via SVEMnet")
    if (length(unknown)) dots <- dots[setdiff(names(dots), unknown)]
  }

  exclude_user <- dots$exclude
  dots$exclude <- NULL
  glmnet_design <- .svem_prepare_glmnet_design(
    X, args = dots, exclude = exclude_user,
    numerically_constant = numerically_constant
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
    default_control = c(list(maxit = 1e6),
      if (identical(fam_name, "gaussian")) list(fdev = 0, devmax = 1)),
    scoped_legacy = TRUE
  )

  list(design = glmnet_design, args = glmnet_base_args, dots = dots)
}

# glmnet >= 5 has per-call controls. Older supported releases require a scoped
# session override; restore every setting even when the solver raises an error.
.svem_glmnet_call <- function(X, y, weights, alpha, args) {
  controls <- attr(args, "svem_legacy_control", exact = TRUE)
  if (length(controls)) {
    old <- glmnet::glmnet.control()
    on.exit(do.call(glmnet::glmnet.control, old), add = TRUE)
    do.call(glmnet::glmnet.control, controls)
  }
  do.call(glmnet::glmnet, c(list(x = X, y = y, alpha = alpha, weights = weights), args))
}

# At the stored lambda values coef() only retrieves these matrices. Avoid its
# interpolation and sparse reconstruction overhead, including singleton drops.
.svem_glmnet_coefficient_matrix <- function(fit) {
  beta <- fit$beta
  if (is.null(dim(beta))) stop("Inconsistent glmnet coefficient dimensions.", call. = FALSE)
  beta <- as.matrix(beta)
  if (ncol(beta) != length(fit$lambda) || length(fit$a0) != ncol(beta))
    stop("Inconsistent glmnet coefficient path.", call. = FALSE)
  out <- rbind(as.numeric(fit$a0), beta)
  rownames(out) <- c("(Intercept)", rownames(beta))
  colnames(out) <- if (is.null(names(fit$lambda)))
    paste0("s=", format(fit$lambda)) else names(fit$lambda)
  out
}

.svem_glmnet_blend <- function(fit, base, gamma) {
  if (gamma == 1 || !inherits(fit, "relaxed")) return(base)
  relaxed <- fit$relaxed
  if (is.null(relaxed) || !length(relaxed$lambda))
    stop("glmnet returned no relaxed refit path.", call. = FALSE)
  # Preserve glmnet's documented gamma floor and its maxp endpoint mapping.
  # Do not replace its relaxed coordinate-descent estimator by exact OLS.
  index <- match(fit$lambda, relaxed$lambda, nomatch = 0L)
  if (!any(index > 0L)) stop("glmnet returned no matching relaxed refit.", call. = FALSE)
  index[index == 0L] <- max(index)
  beta <- relaxed$beta
  if (is.null(dim(beta))) beta <- matrix(beta, nrow = nrow(base) - 1L)
  beta <- as.matrix(beta)
  if (nrow(beta) != nrow(base) - 1L || ncol(beta) != length(relaxed$lambda) ||
      length(relaxed$a0) != ncol(beta))
    stop("Inconsistent glmnet relaxed coefficient path.", call. = FALSE)
  g <- max(gamma, 1e-5)
  out <- rbind(g * base[1L, ] + (1 - g) * relaxed$a0[index],
               g * base[-1L, , drop = FALSE] +
                 (1 - g) * beta[, index, drop = FALSE])
  dimnames(out) <- dimnames(base)
  out
}
