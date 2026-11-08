#' Fit several responses over one shared factor specification
#'
#' Fit one independent SVEM model per response using the same specification
#' from \code{\link{bigexp_terms}}. The returned ordinary named list works
#' directly with \code{\link{svem_random_table_multi}} and
#' \code{\link{svem_score_random}} when using the default SVEM method. Optional
#' spline lists are passed through the scorer's separate `spline` argument.
#'
#' Fits run sequentially in the requested response order. With the default
#' method, this helper neither sets nor restores the random seed: it uses the
#' same random draws as equivalent sequential calls to
#' \code{SVEMnet(spec, data, response = ...)}. Fitting choices in \code{...}
#' are passed to the selected fitting method.
#' For a single response already named in the specification,
#' \code{SVEMnet(spec, data)} is sufficient.
#'
#' @param spec A \code{bigexp_spec} from \code{\link{bigexp_terms}}, or a saved
#'   [spline_terms()] specification when `method = "spline"`.
#' @param data A data frame containing the shared predictors and responses.
#' @param responses Character vector of unique response column names, in the
#'   order to fit. Names on this vector are ignored; the returned list is
#'   named by the response column names themselves.
#' @param ... Arguments for the selected fitting method. The default passes
#'   them unchanged to \code{\link{SVEMnet}}, for example \code{nBoot},
#'   \code{objective}, or \code{companions}. The spline method accepts
#'   \code{nBoot}, \code{lambda} and \code{store_member_weights}; with a
#'   \code{bigexp_spec} it also accepts the advanced basis controls documented
#'   in [spline_terms()]. Use \code{responses}, rather than a separate
#'   \code{response} argument.
#'
#' @param method `"svemnet"` (default) keeps the original sequential fitting
#'   path. `"spline"` fits optional Gaussian smooth surfaces, including nominal
#'   study factors and multiple mixture groups, with the same
#'   raw-factor specification. Spline fits can be passed as the separate
#'   `spline` list to [svem_plan_runs()] or [svem_score_random()].
#' @param domain_grid Explicit immutable response-free feasible recipe grid
#'   required for `method = "spline"` with a bigexp specification. A saved
#'   [spline_terms()] specification already contains its domain.
#' @return An ordinary named list in response order: `svem_model` entries for
#'   the default method, or `svem_spline` entries for the optional spline method.
#'   No joint-response fit or additional list wrapper class is introduced.
#'
#' @examples
#' dat <- data.frame(x = seq(-1, 1, length.out = 20), z = rep(c(0, 1), 10))
#' dat$y1 <- 2 + dat$x + 0.1 * sin(seq_len(nrow(dat)))
#' dat$y2 <- 3 - dat$z + 0.1 * cos(seq_len(nrow(dat)))
#' spec <- bigexp_terms(y1 ~ x + z, dat, factorial_order = 1,
#'                     polynomial_order = 1, report = FALSE)
#' set.seed(12)
#' fits <- svem_fit_responses(spec, dat, c("y1", "y2"), nBoot = 3,
#'                            glmnet_alpha = 1, relaxed = FALSE)
#' predict(fits$y1, dat[1:2, ])
#' # One response needs no list helper:
#' one_fit <- SVEMnet(spec, dat, nBoot = 3, glmnet_alpha = 1, relaxed = FALSE)
#'
#' @export
svem_fit_responses <- function(spec, data, responses, ...,
                               method = c("svemnet", "spline"), domain_grid = NULL) {
  method <- match.arg(method)
  if ("response" %in% names(sys.call()) ||
      "response" %in% names(match.call(expand.dots = FALSE)$...)) {
    stop("Use `responses`, not a separate `response` argument.", call. = FALSE)
  }
  if (!inherits(spec, "bigexp_spec") && !(method == "spline" && inherits(spec, "svem_spline_spec"))) {
    stop("`spec` must be a bigexp_spec from bigexp_terms().", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame containing the responses and predictors.",
         call. = FALSE)
  }
  if (!is.character(responses) || !length(responses) || anyNA(responses) ||
      any(!nzchar(responses))) {
    stop("`responses` must be a nonempty character vector of response column names.",
         call. = FALSE)
  }
  if (anyDuplicated(responses)) {
    stop("`responses` must contain unique names; duplicate responses are not allowed.",
         call. = FALSE)
  }
  unknown <- setdiff(responses, names(data))
  if (length(unknown)) {
    stop("Unknown response column(s) in `data`: ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  responses <- unname(responses)
  if (method == "spline") {
    if (inherits(spec, "svem_spline_spec")) {
      if (!is.null(domain_grid)) stop("A saved spline specification already contains its domain_grid.", call. = FALSE)
      return(stats::setNames(lapply(responses, function(r) svem_spline(spec, data, response = r, ...)), responses))
    }
    if (is.null(domain_grid)) stop("Supply an explicit response-free domain_grid for spline fits.", call. = FALSE)
    return(.svem_fit_spline_responses(spec, data, responses, domain_grid, list(...)))
  }
  if (!is.null(domain_grid)) stop("domain_grid is used only with method = 'spline'.", call. = FALSE)
  fits <- lapply(responses, function(response) {
    SVEMnet(spec, data, response = response, ...)
  })
  stats::setNames(fits, responses)
}
