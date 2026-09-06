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
#' @param polynomial_centers Named numeric vector of fixed centers for powered
#'   predictors. Unlisted predictors retain their raw powers.
#'
#' @keywords internal
#' @noRd
.bigexp_build_rhs <- function(vars, cont_vars,
                              factorial_order  = 3L,
                              polynomial_order = 3L,
                              include_pc_2way  = TRUE,
                              include_pc_3way  = FALSE,
                              intercept        = TRUE,
                              polynomial_centers = numeric(0L)) {
  stopifnot(length(vars) >= 1L)

  factorial_order  <- .svem_integer_scalar(factorial_order, "factorial_order", min = 1L)
  polynomial_order <- .svem_integer_scalar(polynomial_order, "polynomial_order", min = 1L)

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
    power_term <- function(x, degree) {
      if (x %in% names(polynomial_centers)) {
        # Embed a round-trip numeric literal: no training-data environment or
        # prediction-batch mean is needed to evaluate the resulting formula.
        center <- format(polynomial_centers[[x]], digits = 17L,
                         scientific = TRUE, trim = TRUE, decimal.mark = ".")
        paste0("I((", x, " - (", center, "))^", degree, ")")
      } else {
        paste0("I(", x, "^", degree, ")")
      }
    }
    # I(X^2), I(X^3), ..., I(X^p)
    poly_terms <- character()
    for (deg in 2L:polynomial_order) {
      poly_terms <- c(poly_terms, vapply(cont_vars, power_term, "", degree = deg,
                                        USE.NAMES = FALSE))
    }
    rhs_parts <- c(rhs_parts, paste(poly_terms, collapse = " + "))

    # Optional Z:I(X^2) (two-way partial-cubic); these are cross terms, so they
    # require the corresponding factorial interaction order to be enabled
    if (isTRUE(include_pc_2way) && factorial_order >= 2L) {
      pc2 <- character()
      for (xi in cont_vars) {
        others <- setdiff(vars, xi)
        if (length(others)) {
          # form: X3:I(X2^2)
          pc2 <- c(pc2, paste0(others, ":", power_term(xi, 2L)))
        }
      }
      if (length(pc2)) rhs_parts <- c(rhs_parts, paste(pc2, collapse = " + "))
    }

    # Optional I(X^2):Z:W (three-way partial-cubic)
    if (isTRUE(include_pc_3way) && factorial_order >= 3L) {
      pc3 <- character()
      for (xi in cont_vars) {
        others <- setdiff(vars, xi)
        if (length(others) >= 2L) {
          for (pair in utils::combn(others, 2L, simplify = FALSE)) {
            pc3 <- c(pc3, paste0(power_term(xi, 2L), ":", pair[[1]], ":", pair[[2]]))
          }
        }
      }
      if (length(pc3)) rhs_parts <- c(rhs_parts, paste(pc3, collapse = " + "))
    }
  }

  rhs <- paste(rhs_parts, collapse = " + ")
  if (!isTRUE(intercept)) rhs <- paste0(rhs, " - 1")
  rhs
}

# ---- 1) spec builder ----------------------------------------------------------
#' Create a deterministic expansion spec for wide polynomial and interaction models
#'
#' \code{bigexp_terms()} builds a specification object that:
#' \itemize{
#'   \item decides which predictors are treated as continuous or categorical,
#'   \item optionally treats selected variables as blocking factors that enter
#'         the model only additively and never in interactions or polynomials,
#'   \item locks factor levels from the supplied data,
#'   \item records the contrast settings used when the model matrix is first built, and
#'   \item constructs a reusable right-hand side (RHS) expression string for a large
#'         expansion that can be shared across multiple responses and datasets.
#' }
#'
#' The expansion for non-blocking predictors can include:
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
#'   \item all other numeric predictors are treated as continuous, and their
#'         observed ranges are stored.
#' }
#'
#' Variables listed in \code{blocking} are included in the spec and are typed
#' using the same rules as other predictors (for example, a numeric blocking
#' variable with many distinct values is treated as continuous). However,
#' blocking variables enter the model only as additive main effects, without
#' interactions or polynomial terms, regardless of \code{factorial_order} or
#' \code{polynomial_order}.
#'
#' Once built, a \code{"bigexp_spec"} can be reused to create consistent
#' expansions for new datasets via \code{\link{bigexp_prepare}} and
#' \code{\link{bigexp_formula}}. The RHS and contrast settings are locked, so
#' the same spec applied to different data produces design matrices with the
#' same columns in the same order (up to missing levels for specific batches).
#'
#' @section Centered polynomial powers:
#' With \code{center_polynomials = TRUE}, a power \code{I(X^k)} is replaced
#' by \code{I((X - c)^k)}, where \code{c} is the arithmetic mean of the
#' finite values of \code{X} in the data supplied to \code{bigexp_terms()}.
#' The same replacement is used wherever that power occurs in a partial-cubic
#' term. Main effects, ordinary factorial interactions, and the non-powered
#' factors in partial-cubic terms remain on their original scales.
#'
#' Each center is computed separately from its predictor, without excluding
#' rows because another predictor or response is missing. This keeps the
#' expansion independent of the response used to build the spec. The constants
#' are recorded in \code{spec$settings$polynomial_centers} and embedded in the
#' formula; \code{bigexp_prepare()}, \code{bigexp_formula()}, and prediction
#' reuse them without recentering new data. Predictor ranges, discrete sampling
#' levels, and supplied data values are unchanged. Only design-matrix
#' construction changes; fitting algorithms are unchanged, although selecting
#' or penalizing terms in the new basis can change the fitted model.
#'
#' Predictors named in \code{mixture_vars} retain raw powers. Mixture membership
#' is explicit and is not inferred from names, ranges, or observed sums.
#'
#' @section Typical workflow:
#' In a typical multi-response workflow you:
#' \enumerate{
#'   \item Call \code{bigexp_terms()} once on your training data to build and
#'         lock the expansion (types, levels, contrasts, RHS).
#'   \item Fit models using \code{spec$formula} and the original data
#'         (for example, \code{SVEMnet(spec$formula, data, ...)} or
#'         \code{lm(spec$formula, data)}).
#'   \item For new batches, call \code{\link{bigexp_prepare}} with the same
#'         \code{spec} so that design matrices have exactly the same columns
#'         and coding.
#'   \item For additional responses on the same factor space, use
#'         \code{\link{bigexp_formula}} to swap the left-hand side while
#'         reusing the locked expansion.
#' }
#'
#' @param formula Main-effects formula of the form \code{y ~ X1 + X2 + G} or
#'   \code{y ~ .}. The right-hand side should contain main effects only; do not
#'   include \code{:} (interactions), \code{^} (factorial shortcuts),
#'   \code{I()} powers, transformed terms such as \code{log(X1)}, or inline
#'   polynomial expansions. The helper will generate interactions and
#'   polynomial terms automatically. Dot exclusions such as
#'   \code{y ~ . - X1} are honored. Note that with \code{y ~ .} every
#'   non-response column of \code{data} becomes a predictor, including any
#'   other response columns present in the data frame; list the main effects
#'   explicitly (or use exclusions) when the data holds several responses.
#' @param data Data frame used to decide types and lock factor levels.
#' @param factorial_order Integer >= 1. Maximum order of factorial interactions
#'   among the non-blocking main effects. For example, 1 gives main effects
#'   only, 2 gives up to two-way interactions, 3 gives up to three-way
#'   interactions, and so on.
#' @param polynomial_order Integer >= 1. Maximum polynomial degree for continuous
#'   non-blocking predictors. A value of 1 means only linear terms; 2 adds
#'   squares \code{I(X^2)}; 3 adds cubes \code{I(X^3)}; in general, all powers
#'   \code{I(X^k)} for \code{k} from 2 up to \code{polynomial_order} are added.
#' @param include_pc_2way Logical. If \code{TRUE} (default),
#'   \code{polynomial_order >= 2}, and \code{factorial_order >= 2}, include
#'   partial-cubic two-way terms of the form \code{Z:I(X^2)} where \code{X} is
#'   continuous and \code{Z} is another non-blocking predictor. (These are
#'   cross terms, so \code{factorial_order = 1} — "main effects only" —
#'   suppresses them.)
#' @param include_pc_3way Logical. If \code{TRUE},
#'   \code{polynomial_order >= 2}, and \code{factorial_order >= 3}, include
#'   partial-cubic three-way terms \code{I(X^2):Z:W} among non-blocking
#'   predictors.
#' @param intercept Logical. If \code{TRUE} (default), include an intercept in the
#'   expansion; if \code{FALSE}, the generated RHS drops the intercept. The
#'   no-intercept form is available for standalone design-matrix workflows;
#'   \code{SVEMnet()} and \code{glmnet_with_cv()} require the intercept convention
#'   used by the published SVEM methods and reject such a specification.
#' @param blocking Optional character vector of column names in \code{data} to
#'   treat as blocking factors. These variables are included in the spec and
#'   typed like other predictors (categorical vs continuous), but they enter the
#'   model only as additive main effects and never appear in interactions,
#'   polynomials, or partial-cubic terms.
#'   \strong{Important:} when using \code{y ~ .}, blocking variables are
#'   automatically excluded from the "non-blocking" predictor set so they do not
#'   trigger a conflict error.
#'   When using an explicit RHS (for example \code{y ~ X1 + X2}), blocking
#'   variables must not also be explicitly listed on the right-hand side.
#' @param discrete_numeric Optional specification of "discrete numeric" predictors
#'   for downstream sampling (for example in \code{svem_random_table_multi()}).
#'   These predictors are still treated as numeric for modeling and expansion
#'   (that is, they remain continuous in the design matrix and may participate
#'   in polynomial and interaction terms). This option only records a finite set
#'   of preferred numeric levels to be used when randomly generating recipes.
#'   Supply either:
#'   \itemize{
#'     \item a character vector of predictor names, in which case the allowed
#'           levels are inferred as the sorted unique finite values observed in
#'           \code{data}; or
#'     \item a named list mapping predictor names to numeric vectors of allowed
#'           levels. If an entry is \code{NULL} or length zero, levels are
#'           inferred from \code{data} for that predictor.
#'   }
#' @param audit How to handle suspicious typing / high-cardinality issues when
#'   building the spec. One of \code{"warn"} (default), \code{"error"}, or
#'   \code{"none"}. Audits cover numeric-like character/factor columns (including
#'   percent strings like \code{"25\%"}), and very high-cardinality categorical
#'   predictors that are likely IDs or mis-typed numerics.
#' @param audit_numeric_rate Numeric in (0,1). If at least this fraction of
#'   non-missing values parse as numeric (after stripping commas and an optional
#'   trailing \code{\%}), the column is flagged as numeric-like when stored as
#'   character/factor.
#' @param audit_unique_ratio Numeric in (0, 1). For categorical predictors, warn/error
#'   if \code{unique(non-missing) / n_nonmissing >= audit_unique_ratio}.
#' @param audit_min_n Integer >= 1. Minimum number of non-missing values required
#'   before audits are applied.
#' @param report Logical. If \code{TRUE} (default), print a compact summary of the
#'   inferred predictor types and settings (via \code{print.bigexp_spec}) when
#'   \code{bigexp_terms()} returns.
#' @param center_polynomials Logical; default \code{FALSE}. If \code{TRUE},
#'   subtract each eligible predictor's fixed training-data mean before
#'   raising it to a power. Applies to numeric, non-blocking predictors except
#'   those in \code{mixture_vars}. See \emph{Centered polynomial powers}.
#' @param mixture_vars Optional character vector naming numeric, non-blocking
#'   predictors whose polynomial powers must remain uncentered. This only
#'   excludes them from polynomial centering; it does not impose mixture sum
#'   constraints or replace \code{mixture_groups} in downstream sampling.
#'
#' @return An object of class \code{"bigexp_spec"} with components:
#' \itemize{
#'   \item \code{formula}: expanded formula of the form \code{y ~ <big expansion>},
#'         using the response from the input formula.
#'   \item \code{rhs}: right-hand side expansion string (reusable for any response).
#'   \item \code{vars}: character vector of predictor names (including blocking
#'         variables) in the order discovered from the formula and data.
#'   \item \code{is_cat}: named logical vector indicating which predictors are
#'         treated as categorical (\code{TRUE}) versus continuous (\code{FALSE}).
#'   \item \code{levels}: list of locked factor levels for categorical predictors.
#'   \item \code{num_range}: 2 x p numeric matrix of ranges for continuous variables
#'         (rows \code{c("min","max")}).
#'   \item \code{settings}: list of expansion settings, including
#'         \code{factorial_order}, \code{polynomial_order},
#'         \code{include_pc_2way},
#'         \code{include_pc_3way}, \code{intercept}, \code{blocking}, and
#'         stored contrast information. When polynomial centering or mixture
#'         exclusions are requested, also records \code{center_polynomials},
#'         the named numeric vector \code{polynomial_centers}, and
#'         \code{mixture_vars}.
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
#' ## Optional centering of powers using fixed training-data means
#' spec_centered <- bigexp_terms(
#'   y ~ X1 + X2 + G, data = df,
#'   factorial_order = 2, polynomial_order = 2,
#'   center_polynomials = TRUE
#' )
#' spec_centered$settings$polynomial_centers
#' ## For mixture-process data, also supply mixture_vars = c("A", "B", "C")
#' ## to keep the named mixture components' powers on their raw scales.
#'
#' ## Example 2: pure main effects (no interactions, no polynomial terms)
#' spec_main <- bigexp_terms(
#'   y ~ X1 + X2 + G,
#'   data             = df,
#'   factorial_order  = 1,  # main effects only
#'   polynomial_order = 1   # no I(X^2) or higher
#' )
#'
#' ## Example 3: blocking factors (categorical and continuous)
#' set.seed(2)
#' df_block <- data.frame(
#'   y           = rnorm(30),
#'   X1          = rnorm(30),
#'   X2          = rnorm(30),
#'   G           = factor(sample(c("A", "B"), 30, replace = TRUE)),
#'   Operator    = factor(sample(paste0("Op", 1:3), 30, replace = TRUE)),
#'   AmbientTemp = rnorm(30, mean = 22, sd = 2)  # continuous blocking covariate
#' )
#'
#' ## Here Operator (categorical) and AmbientTemp (continuous) are blocking factors:
#' ## they enter additively, but do not appear in interactions or polynomials.
#' spec_block <- bigexp_terms(
#'   y ~ X1 + X2 + G,
#'   data             = df_block,
#'   factorial_order  = 2,
#'   polynomial_order = 3,
#'   blocking         = c("Operator", "AmbientTemp")
#' )
#'
#' print(spec_block)
#' spec_block$rhs
#'
#' ## Example 4: discrete numeric predictors (finite numeric support)
#' ## A common case is a numeric process setting that only takes a small set
#' ## of allowed values (e.g., 0.5, 1, 2, 4). Use `discrete_numeric` in
#' ## bigexp_terms() so downstream sampling respects those levels automatically.
#' set.seed(3)
#' D_allowed <- c(0.5, 1, 2, 4)
#' df_disc <- data.frame(
#'   y  = rnorm(40),
#'   D  = sample(D_allowed, 40, replace = TRUE),   # numeric with discrete support
#'   X1 = rnorm(40),
#'   G  = factor(sample(c("A", "B"), 40, replace = TRUE))
#' )
#'
#' # Record that D should be treated as "discrete numeric" for downstream sampling.
#' # Levels are inferred automatically from the training data.
#' spec_disc <- bigexp_terms(
#'   y ~ D + X1 + G,
#'   data             = df_disc,
#'   factorial_order  = 2,
#'   polynomial_order = 2,
#'   discrete_numeric = "D"
#' )
#'
#'
#' # Fit. The discrete support is retained in fit$sampling_schema.
#' fit_disc <- SVEMnet(spec_disc, df_disc, nBoot = 8,
#'                     glmnet_alpha = 1, relaxed = FALSE)
#'
#' # Score random candidates; sampled D values stay in D_allowed
#' scored <- svem_score_random(
#'   objects         = list(y = fit_disc),
#'   goals           = list(y = list(goal = "max", weight = 1)),
#'   n               = 200,
#'   numeric_sampler = "random",
#'   verbose         = FALSE
#' )
#'
#' table(scored$score_table$D)
#' stopifnot(all(scored$score_table$D %in% D_allowed))
#' @export
#' @importFrom stats model.frame na.pass as.formula model.matrix
#' @importFrom utils combn
bigexp_terms <- function(formula, data,
                         factorial_order    = 3L,
                         polynomial_order   = 3L,
                         include_pc_2way    = TRUE,
                         include_pc_3way    = FALSE,
                         intercept          = TRUE,
                         blocking           = NULL,
                         discrete_numeric   = NULL,
                         audit              = c("warn", "error", "none"),
                         audit_numeric_rate = 0.90,
                         audit_unique_ratio = 0.80,
                         audit_min_n        = 12L,
                         report             = TRUE,
                         center_polynomials = FALSE,
                         mixture_vars       = NULL) {
  stopifnot(is.data.frame(data))

  factorial_order  <- .svem_integer_scalar(factorial_order, "factorial_order", min = 1L)
  polynomial_order <- .svem_integer_scalar(polynomial_order, "polynomial_order", min = 1L)
  audit_numeric_rate <- .svem_numeric_scalar(
    audit_numeric_rate, "audit_numeric_rate", lower = 0, upper = 1,
    lower_open = TRUE
  )
  audit_unique_ratio <- .svem_numeric_scalar(
    audit_unique_ratio, "audit_unique_ratio", lower = 0, upper = 1,
    lower_open = TRUE
  )
  audit_min_n <- .svem_integer_scalar(audit_min_n, "audit_min_n", min = 1L)
  report <- .svem_logical_scalar(report, "report")

  ## Validate logical flags up front so that the stored settings, the printed
  ## spec, and the generated RHS cannot disagree about what was requested.
  intercept <- .svem_logical_scalar(intercept, "intercept")
  include_pc_2way <- .svem_logical_scalar(include_pc_2way, "include_pc_2way")
  include_pc_3way <- .svem_logical_scalar(include_pc_3way, "include_pc_3way")
  center_polynomials <- .svem_logical_scalar(center_polynomials, "center_polynomials")
  if (is.null(mixture_vars)) {
    mixture_vars <- character(0L)
  } else if (!is.character(mixture_vars) || anyNA(mixture_vars) ||
             any(!nzchar(mixture_vars))) {
    stop("mixture_vars must be a character vector of nonempty predictor names, or NULL.",
         call. = FALSE)
  } else {
    mixture_vars <- unique(mixture_vars)
  }

  ## Validate blocking
  if (is.null(blocking)) {
    blocking <- character(0L)
  } else {
    if (!is.character(blocking)) {
      stop("blocking must be a character vector of column names, or NULL.")
    }
    blocking <- unique(blocking)
    missing_blocking <- setdiff(blocking, names(data))
    if (length(missing_blocking)) {
      stop(
        "Blocking variable(s) not found in 'data': ",
        paste(missing_blocking, collapse = ", ")
      )
    }
  }

  ftxt <- paste(deparse(formula), collapse = "")
  if (grepl("[:^]|I\\(", ftxt)) {
    stop(
      "Provide main effects only on the RHS (or use y ~ .). ",
      "The helper will generate interactions and powers."
    )
  }
  dot_rhs <- "." %in% all.vars(formula)

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

  # Disallow response as a blocking variable
  if (length(blocking) && resp %in% blocking) {
    stop(
      "Blocking variables must not include the response '", resp, "'. ",
      "Remove it from the 'blocking' argument."
    )
  }

  vars <- attr(tt, "term.labels")

  ## Robust handling for y ~ . (and any case where '.' leaks into term labels)
  if (length(vars) && any(vars == ".")) {
    vars <- vars[vars != "."]
  }
  if (dot_rhs) {
    # '.' was expanded by terms() against 'data', so term.labels already honor
    # exclusions such as 'y ~ . - X1'. Blocking variables are excluded here
    # because they are added separately via the 'blocking' argument.
    if (length(vars)) {
      vars <- setdiff(vars, blocking)
    } else {
      # Fallback when '.' did not expand: use training-data columns,
      # excluding the response and blocking variables.
      vars <- setdiff(names(data), c(resp, blocking))
    }
  }

  if (!length(vars) && !length(blocking)) {
    stop("No predictors found on the right hand side of formula, and no blocking variables supplied.")
  }

  # Disallow variables listed both in the RHS and in 'blocking' for explicit RHS.
  # For y ~ ., blocking variables are excluded above.
  if (length(blocking) && !dot_rhs) {
    conflict_blocking <- intersect(blocking, vars)
    if (length(conflict_blocking)) {
      stop(
        "The following variables are listed both in the model RHS and in 'blocking': ",
        paste(conflict_blocking, collapse = ", "),
        ". Blocking variables are included automatically via the 'blocking' argument; ",
        "remove them from the formula RHS."
      )
    }
  }

  ## All predictors in the spec: union of RHS vars and blocking cols
  vars_all <- if (length(blocking)) unique(c(vars, blocking)) else vars

  ## Every predictor must be an actual column of 'data'. Transformed RHS terms
  ## such as log(X1) or poly(X1, 2) produce term labels that are not columns;
  ## catch them here with a clear message instead of failing later in
  ## bigexp_prepare() with a confusing "missing predictor" error.
  missing_cols <- setdiff(vars_all, names(data))
  if (length(missing_cols)) {
    # Non-syntactic column names (e.g. "Flow Rate (mL/min)") arrive from
    # terms() wrapped in backticks and would otherwise be misdiagnosed as
    # transformed terms; give the real reason instead.
    debacktick <- gsub("^`|`$", "", missing_cols)
    nonsyntactic <- missing_cols[debacktick %in% names(data)]
    if (length(nonsyntactic)) {
      stop(
        "Column name(s) not syntactically valid in R: ",
        paste(gsub("^`|`$", "", nonsyntactic), collapse = ", "),
        ". bigexp_terms() builds formula text from column names, so they must be ",
        "syntactic; rename the columns first (e.g. names(data) <- make.names(names(data)))."
      )
    }
    stop(
      "Predictor(s) not found as columns in 'data': ",
      paste(missing_cols, collapse = ", "),
      ". bigexp_terms() supports only untransformed main effects on the RHS; ",
      "apply any transformations to the data before calling bigexp_terms()."
    )
  }

  # ---- discrete-numeric parsing (record-only; modeling stays numeric) ----
  discrete_numeric_vars   <- character(0L)
  discrete_numeric_levels <- list()

  if (!is.null(discrete_numeric)) {
    if (is.character(discrete_numeric)) {
      discrete_numeric_vars <- unique(discrete_numeric)
      if (any(!nzchar(discrete_numeric_vars))) {
        stop("discrete_numeric must not contain empty names.")
      }
      bad <- setdiff(discrete_numeric_vars, vars_all)
      if (length(bad)) {
        stop(
          "discrete_numeric variable(s) not found among predictors: ",
          paste(bad, collapse = ", ")
        )
      }
      # levels inferred from data below (after dat0 is created)
    } else if (is.list(discrete_numeric)) {
      if (is.null(names(discrete_numeric)) || any(!nzchar(names(discrete_numeric)))) {
        stop(
          "When discrete_numeric is a list, it must be a *named* list mapping ",
          "predictor names to numeric vectors of allowed levels."
        )
      }
      discrete_numeric_vars <- unique(names(discrete_numeric))
      bad <- setdiff(discrete_numeric_vars, vars_all)
      if (length(bad)) {
        stop(
          "discrete_numeric variable(s) not found among predictors: ",
          paste(bad, collapse = ", ")
        )
      }
      discrete_numeric_levels <- discrete_numeric
    } else {
      stop(
        "discrete_numeric must be NULL, a character vector of predictor names, ",
        "or a named list mapping predictor names to numeric vectors of allowed levels."
      )
    }
  }

  is_cat      <- setNames(logical(length(vars_all)), vars_all)
  levels_list <- vector("list", length(vars_all)); names(levels_list) <- vars_all
  num_range   <- matrix(
    NA_real_, nrow = 2, ncol = 0,
    dimnames = list(c("min", "max"), character())
  )
  dat0  <- as.data.frame(data)
  audit <- match.arg(audit)

  .emit <- function(...) {
    if (identical(audit, "warn")) warning(..., call. = FALSE)
    if (identical(audit, "error")) stop(...)
    invisible(NULL) # "none" => do nothing
  }

  # ---- preflight audit for common "CSV typing" pitfalls ----
  for (v in vars_all) {
    x <- dat0[[v]]

    if (is.character(x) || is.factor(x)) {
      # requires helper defined elsewhere in file:
      # .bigexp_numeric_like_info()
      info <- .bigexp_numeric_like_info(x)

      if (is.finite(info$rate) && info$rate >= audit_numeric_rate) {
        vals <- as.character(x)
        nonmiss <- !is.na(vals) & nzchar(trimws(vals))
        n <- sum(nonmiss)
        if (n >= audit_min_n) {
          u <- length(unique(vals[nonmiss]))
          prop_unique <- u / n

          if (prop_unique >= audit_unique_ratio) {
            .emit(
              "Predictor '", v, "' is ",
              if (info$kind == "percent") "percent-like" else "numeric-like",
              " but stored as ", class(x)[1L], " with very high cardinality (",
              u, " unique / ", n, " non-missing; ", sprintf("%.0f%%", 100 * prop_unique), ").\n",
              "This will be treated as *categorical* and can explode the design matrix.\n",
              "Example values: ", paste(info$example, collapse = ", "), "\n",
              "Optional Fix (if current behavior not intended): convert it to numeric before bigexp_terms() (e.g., strip '%' and divide by 100), ",
              "or exclude/block it if it's an ID."
            )
          } else {
            .emit(
              "Predictor '", v, "' looks ",
              if (info$kind == "percent") "percent-like" else "numeric-like",
              " but is stored as ", class(x)[1L], ". bigexp_terms() will treat it as *categorical*.\n",
              "Example values: ", paste(info$example, collapse = ", "), "\n",
              "Optional Fix (if current behavior not intended): convert it to numeric before bigexp_terms()"
            )
          }
        }
      }
    }
  }

  blocking_vars <- intersect(vars_all, blocking)

  ## Decide which predictors are categorical vs continuous, and lock levels/ranges
  for (v in vars_all) {
    x <- dat0[[v]]

    if (is.factor(x) || is.character(x) || is.logical(x)) {
      is_cat[v] <- TRUE
      fx <- factor(x)
      levels_list[[v]] <- levels(fx)
      dat0[[v]] <- fx

    } else if (is.numeric(x)) {
      # ALWAYS treat numeric as continuous; no auto-discretization
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

    } else {
      is_cat[v] <- TRUE
      fx <- factor(x)
      levels_list[[v]] <- levels(fx)
      dat0[[v]] <- fx
    }
  }

  ## Categorical predictors need at least 2 observed levels; otherwise
  ## model.matrix() fails later with an unhelpful base-R contrasts error.
  for (v in vars_all) {
    if (isTRUE(is_cat[[v]]) && length(levels_list[[v]]) < 2L) {
      stop(
        "Categorical predictor '", v, "' has fewer than 2 levels in the ",
        "training data (it is constant or all-missing). Remove it from the ",
        "formula/blocking, or supply data in which it varies."
      )
    }
  }

  # ---- high-cardinality categorical guardrail (IDs / accidental strings) ----
  for (v in vars_all) {
    if (isTRUE(is_cat[[v]])) {
      vals <- as.character(dat0[[v]])
      nonmiss <- !is.na(vals) & nzchar(trimws(vals))
      n <- sum(nonmiss)
      if (n >= audit_min_n) {
        u <- length(unique(vals[nonmiss]))
        prop_unique <- u / n
        if (prop_unique >= audit_unique_ratio) {
          .emit(
            "Categorical predictor '", v, "' has very high cardinality (",
            u, " unique / ", n, " non-missing; ", sprintf("%.0f%%", 100 * prop_unique), ").\n",
            "This often indicates an ID column or a numeric column imported as text, and can greatly expand the model matrix.\n",
            "Consider excluding it, making it a blocking factor, or converting to numeric if appropriate."
          )
        }
      }
    }
  }

  ## Finalize discrete numeric levels (either provided or inferred from data)
  if (length(discrete_numeric_vars)) {
    for (v in discrete_numeric_vars) {
      x <- dat0[[v]]

      if (!is.numeric(x)) {
        stop(
          "discrete_numeric variable '", v, "' must be numeric in the training data."
        )
      }

      lv_user <- NULL
      if (length(discrete_numeric_levels) && !is.null(discrete_numeric_levels[[v]])) {
        lv_user <- discrete_numeric_levels[[v]]
      }

      if (is.null(lv_user) || !length(lv_user)) {
        lv <- sort(unique(as.numeric(x[is.finite(x)])))
      } else {
        lv <- sort(unique(as.numeric(lv_user)))
      }

      lv <- lv[is.finite(lv)]
      if (!length(lv)) {
        stop(
          "Could not determine any finite discrete levels for '", v, "'. ",
          "Provide finite numeric values in data or supply them via discrete_numeric."
        )
      }

      discrete_numeric_levels[[v]] <- lv
    }

    # ensure list is named and limited to discrete_numeric_vars
    discrete_numeric_levels <- discrete_numeric_levels[discrete_numeric_vars]
  } else {
    discrete_numeric_levels <- list()
  }

  ## Non-blocking predictors are eligible for factorial / polynomial terms
  main_vars <- setdiff(vars_all, blocking_vars)
  cont_vars <- main_vars[!is_cat[main_vars]]

  invalid_mixture <- setdiff(mixture_vars, cont_vars)
  if (length(invalid_mixture)) {
    stop("mixture_vars must name numeric, non-blocking predictors in the expansion: ",
         paste(invalid_mixture, collapse = ", "), call. = FALSE)
  }
  polynomial_centers <- numeric(0L)
  if (center_polynomials && polynomial_order >= 2L) {
    eligible <- setdiff(cont_vars, mixture_vars)
    polynomial_centers <- vapply(eligible, function(v) {
      x <- dat0[[v]]
      mean(x[is.finite(x)])
    }, numeric(1L))
    if (any(!is.finite(polynomial_centers))) {
      stop("Could not compute finite polynomial centers for all eligible predictors.",
           call. = FALSE)
    }
  }

  rhs_terms <- character()

  if (length(main_vars) > 0L) {
    rhs_main <- .bigexp_build_rhs(
      vars             = main_vars,
      cont_vars        = cont_vars,
      factorial_order  = factorial_order,
      polynomial_order = polynomial_order,
      include_pc_2way  = include_pc_2way,
      include_pc_3way  = include_pc_3way,
      intercept        = intercept,
      polynomial_centers = polynomial_centers
    )
    rhs_terms <- c(rhs_terms, rhs_main)
  }

  if (length(blocking_vars) > 0L) {
    rhs_terms <- c(rhs_terms, blocking_vars)
  }

  rhs <- paste(rhs_terms, collapse = " + ")

  ## If there are no main_vars, handle intercept explicitly
  if (!length(main_vars)) {
    if (!isTRUE(intercept)) {
      if (nzchar(rhs)) {
        rhs <- paste(rhs, "- 1")
      } else {
        rhs <- "- 1"
      }
    }
  }

  # env = baseenv(): the default would capture this function's evaluation
  # frame (including the training data and the temporary model matrix), so
  # every saved spec/model would embed a hidden copy of the data.
  form_expanded <- stats::as.formula(paste(resp, "~", rhs), env = baseenv())

  # Record factor-contrast mapping as used today (and current global options)
  mm_tmp         <- stats::model.matrix(stats::as.formula(paste("~", rhs), env = baseenv()), dat0)
  contrasts_used <- attr(mm_tmp, "contrasts")
  contrasts_opts <- getOption("contrasts")

  spec <- structure(
    list(
      formula   = form_expanded,
      rhs       = rhs,
      vars      = vars_all,
      is_cat    = is_cat,
      levels    = levels_list,
      num_range = num_range,
      settings  = list(
        factorial_order    = factorial_order,
        polynomial_order   = polynomial_order,
        include_pc_2way    = include_pc_2way,
        include_pc_3way    = include_pc_3way,
        intercept          = intercept,
        blocking           = blocking_vars,
        discrete_numeric   = discrete_numeric_vars,
        discrete_levels    = discrete_numeric_levels,
        contrasts          = contrasts_used,
        contrasts_options  = contrasts_opts
      )
    ),
    class = "bigexp_spec"
  )

  # Leave the default spec unchanged, including its metadata and formula text.
  if (center_polynomials || length(mixture_vars)) {
    spec$settings$center_polynomials <- center_polynomials
    spec$settings$polynomial_centers <- polynomial_centers
    spec$settings$mixture_vars <- mixture_vars
  }

  # Attach the spec to its own stored formula so that the documented workflow
  # SVEMnet(spec$formula, data, ...) carries the full spec (blocking,
  # discrete-numeric supports, locked levels/contrasts) exactly like
  # SVEMnet(bigexp_formula(spec, "y"), data, ...). The attached copy's own
  # $formula is the bare pre-attachment formula, so this does not recurse.
  # The bigexp_formula class only adds a print method that hides the
  # attached spec (default formula printing would dump it in full).
  attr(spec$formula, "bigexp_spec") <- spec
  class(spec$formula) <- c("bigexp_formula", "formula")

  if (isTRUE(report)) {
    print(spec)
  }

  spec
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
      # Continuous variable expected; enforce numeric type in new data.
      # NOTE: "discrete_numeric" predictors (when recorded in spec$settings)
      # remain numeric here as well; intermediate values are allowed.
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
  f <- stats::as.formula(paste(response, "~", spec$rhs), env = baseenv())
  attr(f, "bigexp_spec") <- spec
  class(f) <- c("bigexp_formula", "formula")
  f
}

#' Print a formula that carries a \code{bigexp_spec} attribute
#'
#' Prints the formula itself and a one-line note, instead of the default
#' formula print, which would dump the entire attached spec.
#'
#' @param x A formula with an attached \code{bigexp_spec}
#'   (from \code{\link{bigexp_terms}} or \code{\link{bigexp_formula}}).
#' @param ... Passed to the default formula print method.
#' @return \code{x}, invisibly.
#' @export
print.bigexp_formula <- function(x, ...) {
  y <- x
  attr(y, "bigexp_spec") <- NULL
  class(y) <- "formula"
  print(y, ...)
  cat("<locked bigexp_spec attached>\n")
  invisible(x)
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
#'   \code{factorial_order},  \code{polynomial_order},
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
#' predictors that are treated as continuous or categorical. It also reports any
#' variables that were designated as blocking factors and therefore enter the
#' model only additively (no interactions, no polynomials).
#'
#' @param x A "bigexp_spec" object.
#' @param ... Unused.
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
#' print(spec4)
#'
#' ## Example with a blocking factor:
#' set.seed(2)
#' df_block2 <- data.frame(
#'   y           = rnorm(12),
#'   X1          = rnorm(12),
#'   G           = factor(sample(c("A", "B"), 12, replace = TRUE)),
#'   Operator    = factor(sample(letters[1:3], 12, replace = TRUE)),
#'   AmbientTemp = rnorm(12, mean = 22, sd = 1.5)
#' )
#'
#' spec_block2 <- bigexp_terms(
#'   y ~ X1 + G,
#'   data             = df_block2,
#'   factorial_order  = 2,
#'   polynomial_order = 3,
#'   blocking         = c("Operator", "AmbientTemp")
#' )
#'
#' print(spec_block2)
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
  if (!is.null(x$settings$blocking) && length(x$settings$blocking)) {
    cat(
      "  Blocking (additive only): ",
      paste(x$settings$blocking, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$settings$discrete_numeric) && length(x$settings$discrete_numeric)) {
    cat(
      "  Discrete numeric (sampling levels recorded): ",
      paste(x$settings$discrete_numeric, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (isTRUE(x$settings$center_polynomials)) {
    centers <- x$settings$polynomial_centers
    cat("  Centered polynomial powers: ",
        if (length(centers)) paste(names(centers), collapse = ", ") else "none",
        "\n", sep = "")
    if (length(x$settings$mixture_vars)) {
      cat("  Mixture predictors (raw powers): ",
          paste(x$settings$mixture_vars, collapse = ", "), "\n", sep = "")
    }
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

#' @keywords internal
#' @noRd
.bigexp_numeric_like_info <- function(x) {
  vals <- as.character(x)
  vals <- trimws(vals)

  nonmiss <- !is.na(vals) & nzchar(vals)
  if (!any(nonmiss)) {
    return(list(rate = 0, kind = NA_character_, example = character(0)))
  }

  v <- vals[nonmiss]

  # detect percent-style strings
  has_pct <- grepl("%", v, fixed = TRUE)
  pct_rate <- mean(has_pct)

  # clean numeric-ish text
  v_clean <- gsub("%", "", v, fixed = TRUE)
  v_clean <- gsub(",", "", v_clean, fixed = TRUE)
  v_clean <- trimws(v_clean)

  num <- suppressWarnings(as.numeric(v_clean))

  # if mostly percent strings, interpret as percent scale (0-1)
  kind <- "number"
  if (pct_rate >= 0.80) {
    num <- num / 100
    kind <- "percent"
  }

  rate <- mean(is.finite(num))

  # examples for messaging
  ex <- unique(v)
  ex <- ex[seq_len(min(3L, length(ex)))]

  list(rate = rate, kind = kind, example = ex)
}

