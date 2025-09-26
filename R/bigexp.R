# bigexp.R
# Wide polynomial/interaction expansion, deterministic across datasets

# ---- internal utility ---------------------------------------------------------

#' Build RHS for factorial/response-surface/partial-cubic expansion
#' @keywords internal
#' @noRd
.bigexp_build_rhs <- function(vars, cont_vars,
                              factorial_order = 3L,
                              include_pure_cubic = TRUE,
                              include_pc_3way   = FALSE,
                              intercept = TRUE) {
  stopifnot(length(vars) >= 1L, factorial_order %in% c(2L, 3L))

  main_block <- if (length(vars) == 1L) vars else
    paste0("(", paste(vars, collapse = " + "), ")^", factorial_order)

  rhs_parts <- c(main_block)

  if (length(cont_vars)) {
    # Response surface (squares)
    rhs_parts <- c(rhs_parts, paste(paste0("I(", cont_vars, "^2)"), collapse = " + "))

    # Optional pure cubic I(X^3)
    if (isTRUE(include_pure_cubic)) {
      rhs_parts <- c(rhs_parts, paste(paste0("I(", cont_vars, "^3)"), collapse = " + "))
    }

    # Partial-cubic 2-way: I(X^2):Z
    pc2 <- character()
    for (xi in cont_vars) {
      others <- setdiff(vars, xi)
      if (length(others)) pc2 <- c(pc2, paste0("I(", xi, "^2):", others))
    }
    if (length(pc2)) rhs_parts <- c(rhs_parts, paste(pc2, collapse = " + "))

    # Optional partial-cubic 3-way: I(X^2):Z:W
    if (isTRUE(include_pc_3way)) {
      pc3 <- character()
      for (xi in cont_vars) {
        others <- setdiff(vars, xi)
        if (length(others) >= 2L) {
          for (pair in utils::combn(others, 2L, simplify = FALSE)) {
            pc3 <- c(pc3, paste0("I(", xi, "^2):", pair[[1]], ":", pair[[2]]))
          }
        }
      }
      if (length(pc3)) rhs_parts <- c(rhs_parts, paste(pc3, collapse = " + "))
    }
  }

  rhs <- paste(rhs_parts, collapse = " + ")
  if (!intercept) rhs <- paste0(rhs, " - 1")
  rhs
}

# ---- 1) spec builder ----------------------------------------------------------

#' Create a deterministic expansion spec for wide polynomial/interaction models
#'
#' Builds a specification that locks variable types and factor levels from `data`,
#' then encodes a large expansion:
#' - Full factorial up to 2- or 3-way among the listed main effects
#' - Response surface (squares of continuous predictors)
#' - Partial cubic crosses (I(X^2):Z, optionally 3-way I(X^2):Z:W)
#'
#' Provide a main-effects formula (e.g., `y ~ X1 + X2 + G` or `y ~ .`).
#' The function constructs the complex RHS and records factor levels so that
#' future datasets expand identically.
#'
#' @param formula Main-effects formula on the RHS (no `:`/`^`/`I()` there). `y ~ .` allowed.
#' @param data Data frame used to decide types and lock factor levels.
#' @param factorial_order Integer, 2 or 3 (default 3).
#' @param discrete_threshold Numeric predictors with <= this many unique finite
#'   values are treated as categorical (default 2).
#' @param include_pure_cubic Logical; include `I(X^3)` for continuous predictors (default FALSE).
#' @param include_pc_3way Logical; include partial-cubic 3-ways `I(X^2):Z:W` (default FALSE).
#' @param intercept Include intercept (default TRUE).
#'
#' @return An object of class `"bigexp_spec"` with components:
#' \itemize{
#' \item `formula` — expanded formula (`y ~ ...`) using the training response.
#' \item `rhs` — the right-hand-side expansion string, reusable for any response.
#' \item `vars` — predictor names (in the order discovered from `formula`/`data`).
#' \item `is_cat` — named logical: categorical vs continuous.
#' \item `levels` — list of locked factor levels (level order preserved).
#' \item `num_range` — 2 x p matrix of ranges for continuous vars (informational).
#' \item `settings` — list of expansion settings, including saved contrasts.
#' }
#'
#' @examples
#' # spec <- bigexp_terms(y ~ X1 + X2 + G, data = df, factorial_order = 3)
#' @export
#' @importFrom stats model.frame na.pass as.formula model.matrix
#' @importFrom utils combn
bigexp_terms <- function(formula, data,
                         factorial_order = 3L,
                         discrete_threshold = 2L,
                         include_pure_cubic = FALSE,
                         include_pc_3way   = FALSE,
                         intercept = TRUE) {
  stopifnot(is.data.frame(data), factorial_order %in% c(2L, 3L))

  ftxt <- paste(deparse(formula), collapse = "")
  if (grepl("[:^]|I\\(", ftxt)) {
    stop("Provide main effects only on the RHS (or use y ~ .). The helper will generate interactions and powers.")
  }

  mf   <- stats::model.frame(formula, data, na.action = stats::na.pass)
  tt   <- attr(mf, "terms")
  resp <- as.character(attr(tt, "variables"))[2L]
  vars <- attr(tt, "term.labels")
  if (length(vars) == 0L && grepl("~\\s*\\.", ftxt)) {
    vars <- setdiff(names(data), resp)
  }
  if (!length(vars)) stop("No predictors found on the RHS of formula.")

  is_cat <- setNames(logical(length(vars)), vars)
  levels_list <- vector("list", length(vars)); names(levels_list) <- vars
  num_range   <- matrix(NA_real_, nrow = 2, ncol = 0,
                        dimnames = list(c("min","max"), character()))
  dat0 <- as.data.frame(data)

  for (v in vars) {
    x <- dat0[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      is_cat[v] <- TRUE
      fx <- factor(x)
      levels_list[[v]] <- levels(fx)
      dat0[[v]] <- fx
    } else if (is.numeric(x)) {
      ux <- unique(x[is.finite(x)])
      if (length(ux) <= discrete_threshold) {
        is_cat[v] <- TRUE
        fx <- factor(x)
        levels_list[[v]] <- levels(fx)
        dat0[[v]] <- fx
      } else {
        is_cat[v] <- FALSE
        r <- range(x, finite = TRUE, na.rm = TRUE)
        num_range <- cbind(num_range, setNames(matrix(r, nrow = 2), v))
      }
    } else {
      is_cat[v] <- TRUE
      fx <- factor(x)
      levels_list[[v]] <- levels(fx)
      dat0[[v]] <- fx
    }
  }

  cont_vars <- vars[!is_cat]
  rhs <- .bigexp_build_rhs(vars, cont_vars,
                           factorial_order, include_pure_cubic,
                           include_pc_3way, intercept)
  form_expanded <- stats::as.formula(paste(resp, "~", rhs))

  # Record factor-contrast mapping as used today (and current global options)
  mm_tmp <- stats::model.matrix(stats::as.formula(paste("~", rhs)), dat0)
  contrasts_used <- attr(mm_tmp, "contrasts")
  contrasts_opts <- getOption("contrasts")

  structure(list(
    formula   = form_expanded,
    rhs       = rhs,
    vars      = vars,
    is_cat    = is_cat,
    levels    = levels_list,
    num_range = num_range,
    settings  = list(
      factorial_order     = factorial_order,
      discrete_threshold  = discrete_threshold,
      include_pure_cubic  = include_pure_cubic,
      include_pc_3way     = include_pc_3way,
      intercept           = intercept,
      contrasts           = contrasts_used,
      contrasts_options   = contrasts_opts
    )
  ), class = "bigexp_spec")
}

# ---- 2) data preparer ---------------------------------------------------------

#' Prepare data to match a bigexp_spec (stable expansion across datasets)
#'
#' Coerces new data to the types/levels locked in `spec` so that
#' `model.matrix(spec$formula, data)` yields the same columns and order
#' every time, even if some factor levels are absent in the new batch.
#'
#' @param spec Object returned by [bigexp_terms()].
#' @param data New data frame (train/test/future batches).
#' @param unseen How to handle unseen factor levels in `data`:
#'   `"warn_na"` (default; convert to `NA` with a warning) or `"error"` (stop).
#'
#' @return `list(formula = spec$formula, data = coerced_data)`
#'
#' @examples
#' # spec <- bigexp_terms(y ~ X1 + X2 + G, train)
#' # new_in <- bigexp_prepare(spec, newdata)
#' # fit <- some_model(new_in$formula, data = new_in$data, ...)
#' @export
bigexp_prepare <- function(spec, data, unseen = c("warn_na","error")) {
  stopifnot(inherits(spec, "bigexp_spec"), is.data.frame(data))
  unseen <- match.arg(unseen)

  dat2 <- as.data.frame(data)
  vars <- spec$vars

  unseen_hits <- list()
  for (v in vars) {
    if (isTRUE(spec$is_cat[[v]])) {
      lv <- spec$levels[[v]]
      if (is.null(lv)) lv <- sort(unique(as.character(dat2[[v]])))
      vals <- as.character(dat2[[v]])
      bad  <- setdiff(unique(vals[!is.na(vals)]), lv)
      if (length(bad)) {
        if (identical(unseen, "error")) stop("Unseen level(s) in ", v, ": ", paste(bad, collapse = ", "))
        unseen_hits[[v]] <- bad
      }
      dat2[[v]] <- factor(vals, levels = lv)
    } else {
      if (!is.numeric(dat2[[v]])) suppressWarnings(dat2[[v]] <- as.numeric(dat2[[v]]))
    }
  }

  if (length(unseen_hits)) {
    msg <- paste0(
      "Unseen levels mapped to NA (rows with NA may be dropped by na.omit):\n",
      paste(vapply(names(unseen_hits), function(nm)
        paste0("  * ", nm, ": ", paste(unseen_hits[[nm]], collapse = ", ")), ""), collapse = "\n")
    )
    warning(msg, call. = FALSE)
  }

  list(formula = spec$formula, data = dat2)
}

# ---- 3) helpers for multi-response & consistency ------------------------------

#' Construct a formula for a new response using a bigexp_spec
#'
#' Useful when you want to fit multiple responses with the same expansion.
#'
#' @param spec A "bigexp_spec".
#' @param response Character scalar: the new response (column name in your data).
#' @return A formula like `response ~ <locked big expansion>`.
#' @examples
#' # f2 <- bigexp_formula(spec, "y2")
#' @export
bigexp_formula <- function(spec, response) {
  stopifnot(inherits(spec, "bigexp_spec"))
  if (missing(response)) return(spec$formula)
  stopifnot(is.character(response), length(response) == 1L, nzchar(response))
  stats::as.formula(paste(response, "~", spec$rhs))
}

#' Build a model matrix using the spec's stored contrasts (if present)
#' @param spec A "bigexp_spec".
#' @param data A data frame to prepare and encode.
#' @return The design matrix returned by `model.matrix()`.
#' @examples
#' # MM <- bigexp_model_matrix(spec, newdata)
#' @export
bigexp_model_matrix <- function(spec, data) {
  inp <- bigexp_prepare(spec, data)

  # Ensure contrasts are functions if stored as character names
  cts <- spec$settings$contrasts
  if (!is.null(cts) && length(cts)) {
    cts <- lapply(cts, function(ci) {
      if (is.character(ci) && length(ci) == 1L && nzchar(ci)) {
        get(ci, mode = "function")
      } else {
        ci
      }
    })
  } else {
    cts <- NULL
  }

  stats::model.matrix(
    spec$formula,
    inp$data,
    contrasts.arg = cts
  )
}

#' Evaluate an expression with the spec's recorded contrast options
#' @param spec A "bigexp_spec".
#' @param code Code to evaluate with temporarily restored `options(contrasts=...)`.
#' @examples
#' # with_bigexp_contrasts(spec, { fit <- SVEMnet(f, data = d, ...) })
#' @export
with_bigexp_contrasts <- function(spec, code) {
  opts <- spec$settings$contrasts_options
  if (is.null(opts)) return(eval.parent(substitute(code)))
  old <- options(contrasts = opts); on.exit(options(old), add = TRUE)
  eval.parent(substitute(code))
}




# ---- 4) train-time convenience wrapper (optional) -----------------------------

#' Build a spec and prepare training data in one call
#'
#' Convenience wrapper around [bigexp_terms()] and [bigexp_prepare()]. It creates
#' the deterministic expansion spec from the training data and immediately
#' prepares that same training data to match the locked types/levels.
#' Prefer the two-step API in documentation, but this is handy for quick scripts.
#'
#' @param formula Main-effects formula (e.g., `y ~ X1 + X2 + G` or `y ~ .`).
#' @param data    Training data frame (used to lock types/levels and prepare).
#' @param ...     Additional arguments forwarded to [bigexp_terms()]
#'                (e.g., `factorial_order`, `discrete_threshold`,
#'                `include_pure_cubic`, `include_pc_3way`, `intercept`).
#'
#' @return A list with components:
#' \itemize{
#'   \item `spec` — the "bigexp_spec" object.
#'   \item `formula` — the expanded formula (`spec$formula`).
#'   \item `data` — the coerced training data ready for modeling.
#' }
#'
#' @examples
#' # tr <- bigexp_train(y ~ X1 + X2 + X3, dat, factorial_order = 3)
#' # fit <- SVEMnet(tr$formula, data = tr$data, ...)
#' @export
bigexp_train <- function(formula, data, ...) {
  stopifnot(is.data.frame(data))
  spec <- bigexp_terms(formula, data, ...)
  prep <- bigexp_prepare(spec, data)
  structure(list(spec = spec, formula = spec$formula, data = prep$data),
            class = "bigexp_train")
}

# ---- S3 print method ----------------------------------------------------------

#' @export
print.bigexp_spec <- function(x, ...) {
  cat("bigexp_spec: ",
      "(factorial_order = ", x$settings$factorial_order,
      ", partial_cubic_3way = ", x$settings$include_pc_3way,
      ", pure_cubic = ", x$settings$include_pure_cubic, ")\n", sep = "")
  cat("  Predictors (", length(x$vars), "): ", paste(x$vars, collapse = ", "), "\n", sep = "")
  if (any(!x$is_cat)) {
    cat("  Continuous: ", paste(names(x$is_cat)[!x$is_cat], collapse = ", "), "\n", sep = "")
  }
  if (any(x$is_cat)) {
    cat("  Categorical: ", paste(names(x$is_cat)[x$is_cat], collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$settings$contrasts_options)) {
    cat("  Stored contrasts options: ",
        paste(x$settings$contrasts_options, collapse = ", "), "\n", sep = "")
  }
  cat("  Formula:\n  ", deparse(x$formula), "\n", sep = "")
  invisible(x)
}
