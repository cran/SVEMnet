# bigexp.R
# Wide polynomial and interaction expansion, deterministic across datasets

# ---- internal utility ---------------------------------------------------------

#' Build RHS for factorial and polynomial expansion
#'
#' Internal helper that constructs a right-hand-side (RHS) expression string
#' containing:
#' \itemize{
#'   \item Full factorial interactions (using \code{:}) among the listed main
#'         effects up to \code{factorial_order}.
#'   \item Polynomial powers \code{I(X^k)} for continuous predictors up to
#'         \code{polynomial_order}.
#'   \item Optional partial-cubic interactions of the form \code{Z:I(X^2)} and
#'         \code{I(X^2):Z:W}.
#' }
#'
#' This function is not exported. It is used by \code{bigexp_terms()} to build
#' the reusable RHS string that can be applied to multiple responses.
#'
#' @param vars Character vector of predictor names (main effects).
#' @param cont_vars Character subset of \code{vars} that are treated as continuous.
#' @param factorial_order Integer >= 1. Maximum order of factorial interactions
#'   among the main effects (for example, 1 means main effects only, 2 adds
#'   two-way interactions, 3 adds three-way interactions, and so on).
#' @param polynomial_order Integer >= 1. Maximum polynomial degree for continuous
#'   predictors. A value of 1 means only linear terms; 2 adds squares
#'   \code{I(X^2)}; 3 adds cubes \code{I(X^3)}; in general, all powers
#'   \code{I(X^k)} for \code{k} from 2 up to \code{polynomial_order} are added.
#' @param include_pc_2way Logical. If \code{TRUE} (default) and
#'   \code{polynomial_order >= 2}, include partial-cubic two-way terms of the
#'   form \code{Z:I(X^2)}, where \code{X} is a continuous predictor and \code{Z}
#'   is any other predictor in \code{vars}.
#' @param include_pc_3way Logical. If \code{TRUE} and \code{polynomial_order >= 2},
#'   include partial-cubic three-way terms of the form \code{I(X^2):Z:W}, where
#'   \code{X} is a continuous predictor and \code{Z}, \code{W} are other
#'   predictors in \code{vars}.
#' @param intercept Logical. If \code{FALSE}, the final RHS explicitly drops the
#'   intercept using \code{" - 1"}.
#'
#' @keywords internal
#' @noRd
.bigexp_build_rhs <- function(vars, cont_vars,
                              factorial_order  = 3L,
                              polynomial_order = 3L,
                              include_pc_2way  = TRUE,
                              include_pc_3way  = FALSE,
                              intercept        = TRUE) {
  stopifnot(length(vars) >= 1L)

  if (!is.numeric(factorial_order) || length(factorial_order) != 1L ||
      !is.finite(factorial_order) || factorial_order < 1) {
    stop("factorial_order must be a single finite integer >= 1.")
  }
  if (!is.numeric(polynomial_order) || length(polynomial_order) != 1L ||
      !is.finite(polynomial_order) || polynomial_order < 1) {
    stop("polynomial_order must be a single finite integer >= 1.")
  }

  factorial_order  <- as.integer(factorial_order)
  polynomial_order <- as.integer(polynomial_order)

  # Explicit factorial block: main effects plus all : interactions up to order k
  build_factorial_terms <- function(vs, k) {
    k <- min(k, length(vs))
    terms <- vs
    if (k >= 2L) {
      for (deg in 2L:k) {
        cmb <- utils::combn(vs, deg, simplify = FALSE)
        terms <- c(terms, vapply(cmb, function(x) paste(x, collapse = ":"), ""))
      }
    }
    terms
  }

  fac_terms <- build_factorial_terms(vars, factorial_order)
  rhs_parts <- c(paste(fac_terms, collapse = " + "))

  # Polynomial powers and partial-cubic crosses
  if (length(cont_vars) > 0L && polynomial_order >= 2L) {
    # I(X^2), I(X^3), ..., I(X^p)
    poly_terms <- character()
    for (deg in 2L:polynomial_order) {
      poly_terms <- c(poly_terms, paste0("I(", cont_vars, "^", deg, ")"))
    }
    rhs_parts <- c(rhs_parts, paste(poly_terms, collapse = " + "))

    # Optional Z:I(X^2) (two-way partial-cubic)
    if (isTRUE(include_pc_2way)) {
      pc2 <- character()
      for (xi in cont_vars) {
        others <- setdiff(vars, xi)
        if (length(others)) {
          # form: X3:I(X2^2)
          pc2 <- c(pc2, paste0(others, ":I(", xi, "^2)"))
        }
      }
      if (length(pc2)) rhs_parts <- c(rhs_parts, paste(pc2, collapse = " + "))
    }

    # Optional I(X^2):Z:W (three-way partial-cubic)
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

#' Create a deterministic expansion spec for wide polynomial and interaction models
#'
#' \code{bigexp_terms()} builds a specification object that:
#' \itemize{
#'   \item decides which predictors are treated as continuous or categorical,
#'   \item locks factor levels from the supplied data,
#'   \item records the contrast settings used when the model matrix is first built, and
#'   \item constructs a reusable right-hand side (RHS) string for a large expansion
#'         that can be shared across multiple responses and datasets.
#' }
#'
#' The expansion can include:
#' \itemize{
#'   \item full factorial interactions among the listed main effects, up to a
#'         chosen order;
#'   \item polynomial terms \code{I(X^k)} for continuous predictors up to a
#'         chosen degree; and
#'   \item optional partial-cubic interactions of the form \code{Z:I(X^2)} and
#'         \code{I(X^2):Z:W}.
#' }
#'
#' Predictor types are inferred from \code{data}:
#' \itemize{
#'   \item factors, characters, and logicals are treated as categorical;
#'   \item numeric predictors with at most \code{discrete_threshold} distinct
#'         finite values are treated as categorical; and
#'   \item all other numeric predictors are treated as continuous, and their
#'         observed ranges are stored.
#' }
#'
#' Once built, a \code{"bigexp_spec"} can be reused to create consistent
#' expansions for new datasets via \code{\link{bigexp_prepare}}, and
#' \code{\link{bigexp_formula}}. The RHS
#' and contrast settings are locked, so the same spec applied to different data
#' produces design matrices with the same columns in the same order (up to
#' missing levels for specific batches).
#'
#' @section Typical workflow:
#' In a typical multi-response workflow you:
#' \enumerate{
#'   \item Call \code{bigexp_terms()} once on your training data to build and
#'         lock the expansion (types, levels, contrasts, RHS).
#'   \item Fit models using \code{spec$formula} and the original data
#'         (for example, \code{SVEMnet(spec$formula, data, ...)} or
#'         \code{lm(spec$formula, data)}).
#'   \item For new batches, call \code{\link{bigexp_prepare}} with the same \code{spec} so that
#'         design matrices have exactly the same columns and coding.
#'   \item For additional responses on the same factor space, use
#'         \code{\link{bigexp_formula}} to swap the left-hand side while
#'         reusing the locked expansion.
#' }
#'
#' @param formula Main-effects formula of the form \code{y ~ X1 + X2 + G} or
#'   \code{y ~ .}. The right-hand side should contain main effects only; do not
#'   include \code{:} (interactions), \code{^} (factorial shortcuts),
#'   \code{I()} powers, or inline polynomial expansions. The helper will
#'   generate interactions and polynomial terms automatically.
#' @param data Data frame used to decide types and lock factor levels.
#' @param factorial_order Integer >= 1. Maximum order of factorial interactions
#'   among the main effects. For example, 1 gives main effects only, 2 gives
#'   up to two-way interactions, 3 gives up to three-way interactions, and so on.
#' @param discrete_threshold Numeric predictors with at most this many unique
#'   finite values are treated as categorical. Default is \code{2}.
#' @param polynomial_order Integer >= 1. Maximum polynomial degree for continuous
#'   predictors. A value of 1 means only linear terms; 2 adds squares \code{I(X^2)};
#'   3 adds cubes \code{I(X^3)}; in general, all powers \code{I(X^k)} for
#'   \code{k} from 2 up to \code{polynomial_order} are added.
#' @param include_pc_2way Logical. If \code{TRUE} (default) and
#'   \code{polynomial_order >= 2}, include partial-cubic two-way terms of the
#'   form \code{Z:I(X^2)}.
#' @param include_pc_3way Logical. If \code{TRUE} and \code{polynomial_order >= 2},
#'   include partial-cubic three-way terms \code{I(X^2):Z:W}.
#' @param intercept Logical. If \code{TRUE} (default), include an intercept in the
#'   expansion; if \code{FALSE}, the generated RHS drops the intercept.
#'
#' @return An object of class \code{"bigexp_spec"} with components:
#' \itemize{
#'   \item \code{formula}: expanded formula of the form \code{y ~ <big expansion>},
#'         using the response from the input formula.
#'   \item \code{rhs}: right-hand side expansion string (reusable for any response).
#'   \item \code{vars}: character vector of predictor names in the order discovered
#'         from the formula and data.
#'   \item \code{is_cat}: named logical vector indicating which predictors are
#'         treated as categorical (\code{TRUE}) versus continuous (\code{FALSE}).
#'   \item \code{levels}: list of locked factor levels for categorical predictors.
#'   \item \code{num_range}: 2 x p numeric matrix of ranges for continuous variables
#'         (rows \code{c("min","max")}).
#'   \item \code{settings}: list of expansion settings, including
#'         \code{factorial_order}, \code{polynomial_order},
#'         \code{discrete_threshold}, \code{include_pc_2way},
#'         \code{include_pc_3way}, \code{intercept}, and stored contrast
#'         information.
#' }
#'
#' @seealso \code{\link{bigexp_prepare}}, \code{\link{bigexp_formula}},
#'   \code{\link{bigexp_train}}
#'
#' @examples
#' ## Example 1: small design with one factor
#' set.seed(1)
#' df <- data.frame(
#'   y  = rnorm(20),
#'   X1 = rnorm(20),
#'   X2 = rnorm(20),
#'   G  = factor(sample(c("A", "B"), 20, replace = TRUE))
#' )
#'
#' ## Two-way interactions and up to cubic terms in X1 and X2
#' spec <- bigexp_terms(
#'   y ~ X1 + X2 + G,
#'   data             = df,
#'   factorial_order  = 2,
#'   polynomial_order = 3
#' )
#'
#' print(spec)
#'
#'
#' ## Example 2: pure main effects (no interactions, no polynomial terms)
#' spec_main <- bigexp_terms(
#'   y ~ X1 + X2 + G,
#'   data             = df,
#'   factorial_order  = 1,  # main effects only
#'   polynomial_order = 1   # no I(X^2) or higher
#' )
#'
#'
#' @export
#' @importFrom stats model.frame na.pass as.formula model.matrix
#' @importFrom utils combn
bigexp_terms <- function(formula, data,
                         factorial_order    = 3L,
                         discrete_threshold = 2L,
                         polynomial_order   = 3L,
                         include_pc_2way    = TRUE,
                         include_pc_3way    = FALSE,
                         intercept          = TRUE) {
  stopifnot(is.data.frame(data))

  if (!is.numeric(factorial_order) || length(factorial_order) != 1L ||
      !is.finite(factorial_order) || factorial_order < 1) {
    stop("factorial_order must be a single finite integer >= 1.")
  }
  if (!is.numeric(polynomial_order) || length(polynomial_order) != 1L ||
      !is.finite(polynomial_order) || polynomial_order < 1) {
    stop("polynomial_order must be a single finite integer >= 1.")
  }
  if (!is.numeric(discrete_threshold) || length(discrete_threshold) != 1L ||
      !is.finite(discrete_threshold) || discrete_threshold < 0) {
    stop("discrete_threshold must be a single non-negative number.")
  }

  factorial_order    <- as.integer(factorial_order)
  polynomial_order   <- as.integer(polynomial_order)
  discrete_threshold <- as.integer(discrete_threshold)

  ftxt <- paste(deparse(formula), collapse = "")
  if (grepl("[:^]|I\\(", ftxt)) {
    stop(
      "Provide main effects only on the RHS (or use y ~ .). ",
      "The helper will generate interactions and powers."
    )
  }

  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  tt <- attr(mf, "terms")

  if (attr(tt, "response") != 1L) {
    stop("Formula must be of the form y ~ X1 + X2 + ... (two-sided).")
  }

  # Require a simple symbol on the LHS (no transformations)
  lhs_call <- attr(tt, "variables")[[2L]]
  if (!is.symbol(lhs_call)) {
    stop(
      "Left-hand side must be a single variable name (e.g., y ~ ...). ",
      "Transformations like log(y) are not supported here."
    )
  }
  resp <- as.character(lhs_call)

  vars <- attr(tt, "term.labels")
  if (length(vars) == 0L && grepl("~\\s*\\.", ftxt)) {
    # Use the terms object: all RHS vars that actually entered the terms
    vars <- all.vars(stats::delete.response(tt))
  }
  if (!length(vars)) {
    stop("No predictors found on the right hand side of formula.")
  }

  is_cat      <- setNames(logical(length(vars)), vars)
  levels_list <- vector("list", length(vars)); names(levels_list) <- vars
  num_range   <- matrix(
    NA_real_, nrow = 2, ncol = 0,
    dimnames = list(c("min", "max"), character())
  )
  dat0 <- as.data.frame(data)

  # Decide which predictors are categorical vs continuous, and lock levels/ranges
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
        if (!all(is.finite(r))) {
          stop(
            "Continuous predictor '", v,
            "' has no finite values in the data used to build the expansion. ",
            "Either drop this predictor or clean/impute the data before calling bigexp_terms()."
          )
        }
        col_mat <- matrix(r, nrow = 2)
        colnames(col_mat) <- v
        num_range <- cbind(num_range, col_mat)
      }
    } else {
      is_cat[v] <- TRUE
      fx <- factor(x)
      levels_list[[v]] <- levels(fx)
      dat0[[v]] <- fx
    }
  }


  cont_vars <- vars[!is_cat]
  rhs <- .bigexp_build_rhs(
    vars             = vars,
    cont_vars        = cont_vars,
    factorial_order  = factorial_order,
    polynomial_order = polynomial_order,
    include_pc_2way  = include_pc_2way,
    include_pc_3way  = include_pc_3way,
    intercept        = intercept
  )
  form_expanded <- stats::as.formula(paste(resp, "~", rhs))

  # Record factor-contrast mapping as used today (and current global options)
  mm_tmp         <- stats::model.matrix(stats::as.formula(paste("~", rhs)), dat0)
  contrasts_used <- attr(mm_tmp, "contrasts")
  contrasts_opts <- getOption("contrasts")

  structure(
    list(
      formula   = form_expanded,
      rhs       = rhs,
      vars      = vars,
      is_cat    = is_cat,
      levels    = levels_list,
      num_range = num_range,
      settings  = list(
        factorial_order    = factorial_order,
        discrete_threshold = discrete_threshold,
        polynomial_order   = polynomial_order,
        include_pc_2way    = include_pc_2way,
        include_pc_3way    = include_pc_3way,
        intercept          = intercept,
        contrasts          = contrasts_used,
        contrasts_options  = contrasts_opts
      )
    ),
    class = "bigexp_spec"
  )
}

# ---- 2) data preparer ---------------------------------------------------------

#' Prepare data to match a \code{bigexp_spec}
#'
#' \code{bigexp_prepare()} coerces a new data frame so that it matches a
#' previously built \code{\link{bigexp_terms}} spec. It:
#' \itemize{
#'   \item applies the locked factor levels for categorical predictors,
#'   \item enforces that continuous variables remain numeric (and errors
#'         if they are not), and
#'   \item optionally warns about or errors on unseen factor levels.
#' }
#'
#' Columns that are not listed in \code{spec$vars} (for example, the response
#' or extra metadata columns) are left unchanged.
#'
#' The goal is that \code{model.matrix(spec$formula, data)} will produce the same set of columns in
#' the same order across all datasets prepared with the same spec, even if some
#' levels are missing in a particular batch.
#'
#' @param spec Object returned by \code{\link{bigexp_terms}}.
#' @param data New data frame (for example, training, test, or future batches).
#' @param unseen How to handle unseen factor levels in \code{data}:
#'   \code{"warn_na"} (default) maps unseen levels to \code{NA} and issues a
#'   warning, or \code{"error"} stops with an error if any unseen levels are
#'   encountered.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{formula}: the expanded formula stored in the spec
#'         (same as \code{spec$formula}).
#'   \item \code{data}: a copy of the input data with predictor columns coerced
#'         to match the spec (types and levels), suitable for
#'         \code{model.frame()} / \code{model.matrix()}.
#' }
#'
#' @seealso \code{\link{bigexp_terms}}
#' @examples
#' set.seed(1)
#' train <- data.frame(
#'   y  = rnorm(10),
#'   X1 = rnorm(10),
#'   X2 = rnorm(10),
#'   G  = factor(sample(c("A", "B"), 10, replace = TRUE))
#' )
#'
#' spec <- bigexp_terms(
#'   y ~ X1 + X2 + G,
#'   data             = train,
#'   factorial_order  = 2,
#'   polynomial_order = 2
#' )
#'
#' newdata <- data.frame(
#'   y  = rnorm(5),
#'   X1 = rnorm(5),
#'   X2 = rnorm(5),
#'   G  = factor(sample(c("A", "B"), 5, replace = TRUE))
#' )
#'
#' prep <- bigexp_prepare(spec, newdata)
#' str(prep$data)
#'
#' @export
bigexp_prepare <- function(spec, data, unseen = c("warn_na", "error")) {
  stopifnot(inherits(spec, "bigexp_spec"), is.data.frame(data))
  unseen <- match.arg(unseen)

  dat2 <- as.data.frame(data)
  vars <- spec$vars

  # Ensure all required predictors are present
  missing <- setdiff(vars, names(dat2))
  if (length(missing)) {
    stop(
      "Missing required predictor(s) in data: ",
      paste(missing, collapse = ", ")
    )
  }

  unseen_hits <- list()
  for (v in vars) {
    if (isTRUE(spec$is_cat[[v]])) {
      lv <- spec$levels[[v]]
      if (is.null(lv)) {
        lv <- sort(unique(as.character(dat2[[v]])))
      }
      vals <- as.character(dat2[[v]])
      bad  <- setdiff(unique(vals[!is.na(vals)]), lv)
      if (length(bad)) {
        if (identical(unseen, "error")) {
          stop("Unseen level(s) in ", v, ": ", paste(bad, collapse = ", "))
        }
        unseen_hits[[v]] <- bad
      }
      dat2[[v]] <- factor(vals, levels = lv)
    } else {
      # Continuous variable expected; enforce numeric type in new data
      if (!is.numeric(dat2[[v]])) {
        stop(
          "Variable '", v,
          "' was treated as numeric when the spec was built, ",
          "but is not numeric in new data."
        )
      }
    }
  }

  if (length(unseen_hits)) {
    msg <- paste0(
      "Unseen levels mapped to NA; downstream model.frame/model.matrix ",
      "may drop or fail on rows containing NA depending on na.action.\n",
      paste(
        vapply(
          names(unseen_hits),
          function(nm) {
            paste0("  * ", nm, ": ", paste(unseen_hits[[nm]], collapse = ", "))
          },
          character(1)
        ),
        collapse = "\n"
      )
    )
    warning(msg, call. = FALSE)
  }

  list(formula = spec$formula, data = dat2)
}

# ---- 3) helpers for multi-response and consistency ----------------------------

#' Construct a formula for a new response using a bigexp_spec
#'
#' bigexp_formula() lets you reuse an existing expansion spec for multiple
#' responses. It keeps the right hand side locked but changes the response
#' variable on the left hand side.
#'
#' This is useful when you want to fit separate models for several responses
#' on the same factor space while guaranteeing that they all use exactly the
#' same design columns and coding.
#'
#' @param spec A "bigexp_spec" object created by bigexp_terms().
#' @param response Character scalar giving the name of the new response column
#'   in your data. If omitted, the original formula is returned unchanged.
#'
#' @return A formula of the form \code{response ~ rhs}, where the right-hand side
#'   is taken from the locked expansion stored in \code{spec}.
#'
#' @examples
#' set.seed(1)
#' df2 <- data.frame(
#'   y1 = rnorm(10),
#'   y2 = rnorm(10),
#'   X1 = rnorm(10),
#'   X2 = rnorm(10)
#' )
#'
#' spec2 <- bigexp_terms(
#'   y1 ~ X1 + X2,
#'   data             = df2,
#'   factorial_order  = 2,
#'   polynomial_order = 2
#' )
#'
#' f2 <- bigexp_formula(spec2, "y2")
#' f2
#'
#' @export
bigexp_formula <- function(spec, response) {
  stopifnot(inherits(spec, "bigexp_spec"))
  if (missing(response)) {
    return(spec$formula)
  }
  stopifnot(is.character(response), length(response) == 1L, nzchar(response))
  stats::as.formula(paste(response, "~", spec$rhs))
}


#' Evaluate code with the spec's recorded contrast options
#'
#' with_bigexp_contrasts() temporarily restores the contrasts options that
#' were active when the spec was built, runs a block of code, and then
#' restores the original options. This is useful when a modeling function
#' uses the global \code{options("contrasts")} to decide how to encode factors
#' (for example, \code{lm()}, \code{glm()}, or other modeling functions that
#' call \code{model.matrix()} internally).
#'
#' @param spec A "bigexp_spec" object with stored contrasts_options in settings.
#' @param code Code to evaluate with temporarily restored options.
#'
#' @examples
#' set.seed(1)
#' df4 <- data.frame(
#'   y  = rnorm(10),
#'   X1 = rnorm(10),
#'   G  = factor(sample(c("A", "B"), 10, replace = TRUE))
#' )
#'
#' spec4 <- bigexp_terms(
#'   y ~ X1 + G,
#'   data             = df4,
#'   factorial_order  = 2,
#'   polynomial_order = 2
#' )
#'
#' with_bigexp_contrasts(spec4, {
#'   mm4 <- model.matrix(spec4$formula, df4)
#'   head(mm4)
#' })
#'
#' @export
with_bigexp_contrasts <- function(spec, code) {
  opts <- spec$settings$contrasts_options
  if (is.null(opts)) {
    return(eval.parent(substitute(code)))
  }
  old <- options(contrasts = opts)
  on.exit(options(old), add = TRUE)
  eval.parent(substitute(code))
}

# ---- 4) train-time convenience wrapper ----------------------------------------

#' Build a spec and prepare training data in one call
#'
#' bigexp_train() is a convenience wrapper around \code{\link{bigexp_terms}} and
#' \code{\link{bigexp_prepare}}. It:
#' \itemize{
#'   \item builds a deterministic expansion spec from the training data; and
#'   \item immediately prepares that same data to match the locked types and levels.
#' }
#'
#' This is handy when you want a single object that contains both the spec
#' and the training data in a form that is ready to pass into a modeling
#' function. For more control, you can call \code{bigexp_terms()} and
#' \code{bigexp_prepare()} explicitly instead.
#'
#' @param formula Main-effects formula such as \code{y ~ X1 + X2 + G} or \code{y ~ .}.
#'   Only main effects should appear on the right hand side.
#' @param data Training data frame used to lock types and levels.
#' @param ... Additional arguments forwarded to \code{bigexp_terms()}, such as
#'   \code{factorial_order}, \code{discrete_threshold}, \code{polynomial_order},
#'   \code{include_pc_2way}, \code{include_pc_3way}, and \code{intercept}.
#'
#' @return An object of class \code{"bigexp_train"} which is a list with components:
#' \itemize{
#'   \item \code{spec}: the \code{"bigexp_spec"} object returned by
#'         \code{bigexp_terms()}.
#'   \item \code{formula}: the expanded formula \code{spec$formula}.
#'   \item \code{data}: the prepared training data (predictors coerced to match
#'         \code{spec}), suitable for passing directly to modeling functions
#'         such as \code{lm()}, \code{glm()}, or \code{SVEMnet()}.
#' }
#'
#' @examples
#' set.seed(1)
#' df5 <- data.frame(
#'   y  = rnorm(20),
#'   X1 = rnorm(20),
#'   X2 = rnorm(20)
#' )
#'
#' tr <- bigexp_train(
#'   y ~ X1 + X2,
#'   data             = df5,
#'   factorial_order  = 2,
#'   polynomial_order = 3
#' )
#'
#' ## Prepared training data and expanded formula:
#' str(tr$data)
#' tr$formula
#'
#' ## Example: fit a model using the expanded formula
#' fit_lm <- lm(tr$formula, data = tr$data)
#' summary(fit_lm)
#'
#' @export
bigexp_train <- function(formula, data, ...) {
  stopifnot(is.data.frame(data))
  spec <- bigexp_terms(formula, data, ...)
  prep <- bigexp_prepare(spec, data)
  structure(
    list(
      spec    = spec,
      formula = spec$formula,
      data    = prep$data
    ),
    class = "bigexp_train"
  )
}

# ---- S3 print method ----------------------------------------------------------

#' Print method for bigexp_spec objects
#'
#' This print method shows a compact summary of the expansion settings and the
#' predictors that are treated as continuous or categorical.
#'
#' @param x A "bigexp_spec" object.
#' @param ... Unused.
#'
#' @export
print.bigexp_spec <- function(x, ...) {
  cat(
    "bigexp_spec: ",
    "(factorial_order = ", x$settings$factorial_order,
    ", polynomial_order = ", x$settings$polynomial_order,
    ", partial_cubic_2way = ", x$settings$include_pc_2way,
    ", partial_cubic_3way = ", x$settings$include_pc_3way,
    ")\n",
    sep = ""
  )
  cat(
    "  Predictors (", length(x$vars), "): ",
    paste(x$vars, collapse = ", "), "\n",
    sep = ""
  )
  if (any(!x$is_cat)) {
    cat(
      "  Continuous: ",
      paste(names(x$is_cat)[!x$is_cat], collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (any(x$is_cat)) {
    cat(
      "  Categorical: ",
      paste(names(x$is_cat)[x$is_cat], collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$settings$contrasts_options)) {
    co <- x$settings$contrasts_options
    cat(
      "  Stored contrasts options: ",
      paste(co, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  cat("  Formula:\n  ", deparse(x$formula), "\n", sep = "")
  invisible(x)
}
