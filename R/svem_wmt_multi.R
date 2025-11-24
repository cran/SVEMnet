#' Whole-model tests for multiple SVEM responses (WMT wrapper)
#'
#' @description
#' Convenience wrapper around \code{\link{svem_significance_test_parallel}}
#' for running whole-model tests (WMT) on multiple responses that share the
#' same dataset and mixture constraints. This helper:
#' \itemize{
#'   \item takes a formula or a list of formulas and a single \code{data} frame,
#'   \item calls \code{svem_significance_test_parallel()} for each response,
#'   \item extracts per-response p-values and converts them to WMT
#'         multipliers via a chosen transform, and
#'   \item optionally plots the WMT objects together using
#'         \code{\link{plot.svem_significance_test}} and prints a compact
#'         summary of p-values and multipliers.
#' }
#'
#' The resulting \code{multipliers} vector is designed to be passed directly
#' to downstream scoring functions (for example, as an optional WMT argument
#' to \code{svem_score_random()}), with response names matched by
#' \code{names()}.
#'
#' @param formulas A single formula or a (preferably named) list of formulas,
#'   one per response. If unnamed, response names are inferred from the
#'   left-hand side of each formula; non-unique names are made unique.
#' @param data Data frame containing the predictors and responses referenced
#'   in \code{formulas}.
#' @param mixture_groups Optional mixture and simplex constraints passed to
#'   \code{\link{svem_significance_test_parallel}}.
#' @param wmt_transform Character; transformation used to convert WMT
#'   p-values into multipliers. One of:
#'   \itemize{
#'     \item \code{"neglog10"}: \eqn{f(p) = [-\log_{10}(p)]^{\text{strength}}},
#'     \item \code{"one_minus_p"}: \eqn{f(p) = (1 - p)^{\text{strength}}}.
#'   }
#'   Currently, \code{strength = 1} is used internally.
#' @param wmt_control Optional list of extra arguments passed directly to
#'   \code{\link{svem_significance_test_parallel}}. By default this is
#'   \code{list(seed = 123)} so that WMT calls are reproducible; you may
#'   override or extend this (e.g. \code{list(seed = 999, nPerm = 300)}).
#'   Any entries not recognized by \code{svem_significance_test_parallel}
#'   are ignored by that function.
#' @param plot Logical; if \code{TRUE} (default), attempt to plot all
#'   successfully computed WMT objects together via
#'   \code{\link{plot.svem_significance_test}}.
#' @param verbose Logical; if \code{TRUE} (default), print progress and a
#'   compact summary of p-values and multipliers.
#'
#' @return
#' A list of class \code{"svem_wmt_multi"} with components:
#' \describe{
#'   \item{\code{wmt_objects}}{Named list of WMT objects (one per response),
#'         as returned by \code{svem_significance_test_parallel()}. Entries
#'         are \code{NULL} where a WMT call failed.}
#'   \item{\code{p_values}}{Named numeric vector of per-response p-values
#'         (bounded away from 0/1), or \code{NA} when unavailable.}
#'   \item{\code{multipliers}}{Named numeric vector of per-response WMT
#'         multipliers derived from the p-values using \code{wmt_transform}.}
#'   \item{\code{wmt_transform}}{The transformation used.}
#'   \item{\code{wmt_control}}{The list of arguments passed through to
#'         \code{svem_significance_test_parallel()}.}
#' }
#'
#' @seealso
#' \code{\link{svem_significance_test_parallel}},
#' \code{\link{plot.svem_significance_test}}
#' @importFrom stats sd ave
#' @examples
#' \donttest{
#' data(lipid_screen)
#'
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_2way  = TRUE,
#'   include_pc_3way  = FALSE
#' )
#'
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' set.seed(123)
#' wmt_out <- svem_wmt_multi(
#'   formulas       = list(Potency = form_pot,
#'                         Size    = form_siz,
#'                         PDI     = form_pdi),
#'   data           = lipid_screen,
#'   mixture_groups = mix,
#'   wmt_transform  = "neglog10",
#'   wmt_control    = list(seed = 123),
#'   plot           = TRUE
#' )
#'
#' wmt_out$p_values
#' wmt_out$multipliers
#'
#' ## later: pass wmt_out$multipliers into svem_score_random()
#' }
#'
#' @export
svem_wmt_multi <- function(formulas,
                           data,
                           mixture_groups = NULL,
                           wmt_transform = c("neglog10", "one_minus_p"),
                           wmt_control = list(seed = 123),
                           plot = TRUE,
                           verbose = TRUE) {

  # ---- basic checks ----
  if (inherits(formulas, "formula")) {
    formulas <- list(formulas)
  }
  if (!is.list(formulas) || !length(formulas)) {
    stop("`formulas` must be a formula or a non-empty list of formulas.")
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  wmt_transform <- match.arg(wmt_transform)

  # ---- infer / clean response names ----
  f_names <- names(formulas)
  if (is.null(f_names) || any(!nzchar(f_names))) {
    f_names <- vapply(formulas, function(fml) {
      if (!inherits(fml, "formula")) return(NA_character_)
      lhs <- tryCatch(as.character(fml[[2L]]),
                      error = function(e) NA_character_)
      if (!is.na(lhs) && nzchar(lhs)) lhs else "response"
    }, character(1))

    # make unique if needed
    if (anyDuplicated(f_names)) {
      idx <- ave(seq_along(f_names), f_names, FUN = seq_along)
      dup <- duplicated(f_names) | (idx > 1L)
      f_names[dup] <- paste0(f_names[dup], "_", idx[dup])
    }
    names(formulas) <- f_names
  }

  # transform constants (aligned with svem_score_random)
  wmt_strength <- 1
  wmt_floor    <- 1e-3

  res_list <- vector("list", length(formulas))
  names(res_list) <- names(formulas)

  p_vals <- setNames(rep(NA_real_, length(formulas)), names(formulas))
  mults  <- setNames(rep(1.0,      length(formulas)), names(formulas))

  # ---- loop over responses ----
  for (nm in names(formulas)) {
    fml <- formulas[[nm]]
    if (!inherits(fml, "formula")) {
      stop("All entries in `formulas` must be formulas. Offender: ", nm)
    }

    if (isTRUE(verbose)) {
      cat("Running whole-model test (WMT) for", nm, "...\n")
    }

    # Build argument list:
    # - always supply formula, data, mixture_groups
    # - allow wmt_control to override (e.g., seed, nPerm, nPoint, spec, ...)
    args <- c(
      list(
        formula        = fml,
        data           = data,
        mixture_groups = mixture_groups
      ),
      wmt_control
    )

    wmt_obj <- tryCatch(
      do.call(svem_significance_test_parallel, args),
      error = function(e) {
        if (isTRUE(verbose)) {
          cat(" -", nm, ": WMT error ->", conditionMessage(e), "\n")
        }
        NULL
      }
    )
    res_list[[nm]] <- wmt_obj

    p <- NA_real_
    if (!is.null(wmt_obj) && is.finite(wmt_obj$p_value)) {
      p <- as.numeric(wmt_obj$p_value)
      # bound away from 0 and 1 to stabilize transforms
      p <- min(max(p, 1e-16), 1 - 1e-12)
    }
    p_vals[nm] <- p

    mult_raw <- if (is.na(p)) {
      1.0
    } else {
      switch(
        wmt_transform,
        neglog10    = (-log10(p))^wmt_strength,
        one_minus_p = (1 - p)^wmt_strength
      )
    }
    mult <- max(mult_raw, wmt_floor)
    mults[nm] <- mult

    if (isTRUE(verbose)) {
      cat(sprintf(" - %s: p = %s, multiplier = %.4g\n",
                  nm,
                  if (is.na(p)) "NA" else format(p, digits = 4),
                  mult))
    }
  }

  # ---- optional combined plot using plot.svem_significance_test ----
  if (isTRUE(plot)) {
    non_null <- Filter(Negate(is.null), res_list)
    if (length(non_null) >= 1L) {
      lbls <- names(non_null)
      # plot(wmt1, wmt2, ..., labels = lbls)
      try(
        {
          args_plot <- c(non_null[1L],
                         non_null[-1L],
                         list(labels = lbls))
          print(do.call(plot, args_plot))
        },
        silent = TRUE
      )
    }
  }

  if (isTRUE(verbose)) {
    cat("\nWMT p-values:\n")
    print(p_vals)
    cat("\nWMT multipliers:\n")
    print(mults)
  }

  out <- list(
    wmt_objects   = res_list,
    p_values      = p_vals,
    multipliers   = mults,
    wmt_transform = wmt_transform,
    wmt_control   = wmt_control
  )
  class(out) <- c("svem_wmt_multi", class(out))
  out
}
