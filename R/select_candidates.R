#' Select diverse candidate rows via PAM medoids
#'
#' @description
#' Internal utility to select a small set of diverse candidate rows from a
#' scored table. The function:
#' \enumerate{
#'   \item Orders rows according to a score column \code{by}, with larger or
#'         smaller values treated as more desirable depending on
#'         \code{direction}.
#'   \item Retains the top fraction \code{top_frac} of rows.
#'   \item Computes a dissimilarity matrix on \code{predictor_cols} using
#'         \code{cluster::daisy()} with the requested metric.
#'   \item Runs \code{cluster::pam()} on the dissimilarities and returns the
#'         row indices of the selected medoids, in the original table.
#' }
#'
#' This helper is used inside optimization routines to pick a small, diverse
#' subset of high scoring or high uncertainty points while respecting the
#' original sampled candidate set.
#'
#' @param table A \code{data.frame} containing candidate rows, including the
#'   score column referenced by \code{by} and the predictor columns listed in
#'   \code{predictor_cols}.
#' @param by Character scalar. Name of the column in \code{table} used for
#'   ranking.
#' @param top_frac Numeric scalar in the open interval \code{(0, 1]} giving
#'   the fraction of top ranked rows to retain before clustering.
#' @param k Integer scalar giving the desired number of medoids (candidates)
#'   to select. If \code{k <= 0} or \code{nrow(table) == 0}, the function
#'   returns \code{integer(0)}. Internally, \code{k} is truncated to the
#'   number of available rows in the top fraction.
#' @param predictor_cols Character vector of column names in \code{table} that
#'   define the space on which diversity is measured. These columns are passed
#'   to \code{cluster::daisy()}.
#' @param direction Character scalar, either \code{"max"} or \code{"min"},
#'   indicating whether larger or smaller values of \code{by} are treated as
#'   more desirable when ranking rows and defining the top fraction.
#' @param metric Character scalar giving the dissimilarity metric passed to
#'   \code{cluster::daisy()}. Defaults to \code{"gower"}.
#'
#' @return
#' An integer vector of row indices (referring to \code{table}) corresponding
#' to the selected PAM medoids drawn from the top \code{top_frac} fraction of
#' rows by \code{by}. If no candidates can be selected (for example because
#' \code{k <= 0} or \code{nrow(table) == 0}), returns \code{integer(0)}.
#'
#' @details
#' Rows are ranked according to \code{by}, with larger or smaller values
#' treated as better depending on \code{direction} ("max" or "min").
#' Rows whose \code{by} value is \code{NA}, \code{NaN}, \code{Inf}, or
#' \code{-Inf} are excluded before ranking and cannot enter the top fraction.
#' The function errors when no finite ranking values remain.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' tab <- data.frame(
#'   score = rnorm(n),
#'   x1    = runif(n),
#'   x2    = runif(n)
#' )
#'
#' # Select 3 medoids from the top 20% by score (maximize score)
#' idx_max <- svem_select_candidates(
#'   table          = tab,
#'   by             = "score",
#'   top_frac       = 0.2,
#'   k              = 3,
#'   predictor_cols = c("x1", "x2"),
#'   direction      = "max"
#' )
#'
#' # Select 3 medoids from the bottom 20% by score (minimize score)
#' idx_min <- svem_select_candidates(
#'   table          = tab,
#'   by             = "score",
#'   top_frac       = 0.2,
#'   k              = 3,
#'   predictor_cols = c("x1", "x2"),
#'   direction      = "min"
#' )
#'
#' @importFrom cluster daisy pam
#' @keywords internal
#' @noRd
svem_select_candidates <- function(table,
                                   by,
                                   top_frac,
                                   k,
                                   predictor_cols,
                                   direction = c("max", "min"),
                                   metric = "gower") {
  # ---- basic checks ----
  if (!is.data.frame(table)) {
    stop("`table` must be a data.frame.")
  }
  n <- nrow(table)
  if (n == 0L) {
    return(integer(0L))
  }

  direction <- match.arg(direction)

  if (!is.character(by) || length(by) != 1L || !nzchar(by)) {
    stop("`by` must be a nonempty character scalar naming a column in `table`.")
  }
  if (!(by %in% colnames(table))) {
    stop("Column `", by, "` not found in `table`.")
  }
  if (!is.numeric(table[[by]])) {
    stop("Column `", by, "` must be numeric for ranking.")
  }

  top_frac <- .svem_numeric_scalar(top_frac, "top_frac", lower = 0, upper = 1,
                                   lower_open = TRUE)

  k <- .svem_integer_scalar(k, "k", min = -.Machine$integer.max)
  if (k <= 0L) {
    return(integer(0L))
  }

  if (!is.character(predictor_cols) || !length(predictor_cols)) {
    stop("`predictor_cols` must be a nonempty character vector of column names.")
  }
  if (anyNA(predictor_cols) || any(!nzchar(predictor_cols))) {
    stop("`predictor_cols` cannot contain NA or empty column names.")
  }
  missing_pred <- setdiff(predictor_cols, colnames(table))
  if (length(missing_pred)) {
    stop("The following `predictor_cols` are not present in `table`: ",
         paste(missing_pred, collapse = ", "))
  }

  # Rank only finite values.  Inf/-Inf are sentinel values in several scoring
  # workflows and must not become the best row or enter the PAM candidate set.
  vals <- table[[by]]
  finite_idx <- which(is.finite(vals))
  if (!length(finite_idx)) {
    stop("Column `", by, "` has no finite values; cannot rank candidates.")
  }
  n_rank <- length(finite_idx)

  # Determine how many top rows to keep; the small epsilon prevents
  # floating-point noise (e.g. 0.2 * 50 = 10.000000000000002) from
  # inflating the count by one
  m_top <- max(1L, min(n_rank, ceiling(top_frac * n_rank - 1e-9)))

  # Deterministic tie-breaker: retain original row order among equal scores.
  rank_key <- if (direction == "max") -vals[finite_idx] else vals[finite_idx]
  ord <- order(rank_key, finite_idx)
  top_idx <- finite_idx[ord[seq_len(m_top)]]

  top_X <- table[top_idx, predictor_cols, drop = FALSE]

  # Coerce problematic types for daisy() (character -> factor; integer64 -> numeric)
  for (nm in names(top_X)) {
    if (is.character(top_X[[nm]])) {
      top_X[[nm]] <- factor(top_X[[nm]])
    } else if (inherits(top_X[[nm]], "integer64")) {
      top_X[[nm]] <- as.numeric(top_X[[nm]])
    }
  }

  # If k is larger than available rows, truncate
  k <- min(k, nrow(top_X))
  if (k <= 0L) {
    return(integer(0L))
  }

  # cluster::pam() requires k <= n - 1 (and daisy() needs n >= 2). When the
  # requested count covers the whole top set there is nothing to cluster:
  # return all top rows directly instead of crashing.
  if (k >= nrow(top_X)) {
    return(top_idx)
  }

  # Critical: reset rownames so pam id.med maps to positions 1..m_top reliably
  rownames(top_X) <- as.character(seq_len(nrow(top_X)))

  # Dissimilarities and PAM medoids (cluster is a hard dependency)
  # daisy() warns when a numeric column has only 2 unique values; still interval-scaled.

  d <- withCallingHandlers(
    cluster::daisy(top_X, metric = metric),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("binary variable\\(s\\).*treated as interval scaled", msg, ignore.case = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )


  pam_fit <- cluster::pam(d, k = k, diss = TRUE)

  # pam_fit$id.med may be labels; map back robustly
  med_id <- pam_fit$id.med
  med_pos <- suppressWarnings(as.integer(med_id))

  if (anyNA(med_pos)) {
    # fall back to matching by labels if coercion fails
    labs <- attr(d, "Labels")
    med_pos <- match(as.character(med_id), labs)
  }

  med_pos <- med_pos[!is.na(med_pos)]
  if (!length(med_pos)) {
    return(integer(0L))
  }

  top_idx[med_pos]
}
