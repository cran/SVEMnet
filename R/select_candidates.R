#' Select diverse candidate rows via PAM medoids
#'
#' @description
#' Internal utility to select a small set of diverse candidate rows from a
#' scored table. The function:
#' \enumerate{
#'   \item Orders rows in decreasing order of a score column \code{by}.
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
#'   ranking. Rows are sorted in decreasing order of this column; the highest
#'   values are treated as most desirable.
#' @param top_frac Numeric scalar in the open interval \code{(0, 1]} giving
#'   the fraction of top ranked rows to retain before clustering.
#' @param k Integer scalar giving the desired number of medoids (candidates)
#'   to select. If \code{k <= 0} or \code{nrow(table) == 0}, the function
#'   returns \code{integer(0)}. Internally, \code{k} is truncated to the
#'   number of available rows in the top fraction.
#' @param predictor_cols Character vector of column names in \code{table} that
#'   define the space on which diversity is measured. These columns are passed
#'   to \code{cluster::daisy()}.
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
#' The function assumes that larger values of \code{by} are better and always
#' ranks in decreasing order. Any rows with missing values in the \code{by}
#' column are placed at the end by \code{order()} and may or may not enter the
#' top fraction, depending on \code{top_frac} and the number of nonmissing
#' rows.
#' @importFrom cluster daisy pam
#' @keywords internal
#' @noRd
svem_select_candidates <- function(table,
                                   by,
                                   top_frac,
                                   k,
                                   predictor_cols,
                                   metric = "gower") {
  # Basic checks
  if (!is.data.frame(table)) {
    stop("`table` must be a data.frame.")
  }
  n <- nrow(table)
  if (n == 0L || k <= 0L) {
    return(integer(0L))
  }

  if (!is.character(by) || length(by) != 1L || !nzchar(by)) {
    stop("`by` must be a nonempty character scalar naming a column in `table`.")
  }
  if (!(by %in% colnames(table))) {
    stop("Column `", by, "` not found in `table`.")
  }
  if (!is.numeric(table[[by]])) {
    stop("Column `", by, "` must be numeric for ranking.")
  }

  if (!is.numeric(top_frac) || length(top_frac) != 1L ||
      !is.finite(top_frac) || top_frac <= 0 || top_frac > 1) {
    stop("`top_frac` must be a single finite number in (0, 1].")
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) {
    stop("`k` must be a single finite numeric value.")
  }
  k <- as.integer(k)

  if (!is.character(predictor_cols) || !length(predictor_cols)) {
    stop("`predictor_cols` must be a nonempty character vector of column names.")
  }
  missing_pred <- setdiff(predictor_cols, colnames(table))
  if (length(missing_pred)) {
    stop("The following `predictor_cols` are not present in `table`: ",
         paste(missing_pred, collapse = ", "))
  }

  # Determine how many top rows to keep
  m_top <- max(1L, min(n, ceiling(top_frac * n)))

  # Order by score, decreasing; NA go to the end
  ord     <- order(table[[by]], decreasing = TRUE, na.last = TRUE)
  top_idx <- ord[seq_len(m_top)]

  top_X <- table[top_idx, predictor_cols, drop = FALSE]

  # If k is larger than available rows, truncate
  k <- min(k, m_top)
  if (k <= 0L) {
    return(integer(0L))
  }

  # Dissimilarities and PAM medoids (cluster is a hard dependency)
  d <- cluster::daisy(top_X, metric = metric)
  pam_fit <- cluster::pam(d, k = k, diss = TRUE)

  med_top_pos    <- pam_fit$id.med
  med_global_idx <- top_idx[med_top_pos]

  med_global_idx
}
