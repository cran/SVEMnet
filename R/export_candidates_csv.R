#' Export SVEM candidate sets to CSV
#'
#' @description
#' Given one or more selection objects returned by
#' \code{\link{svem_select_from_score_table}}, concatenate their
#' \code{$best} rows and \code{$candidates} and export a CSV suitable for
#' planning new experimental runs.
#'
#' Each row is tagged with:
#' \itemize{
#'   \item \code{candidate_type}: \code{"best"} or \code{"medoid"}.
#'   \item \code{selection_label}: derived from the \code{label} argument
#'         used in \code{svem_select_from_score_table()} when available.
#' }
#'
#' The function does not modify any response or prediction columns (for
#' example, \code{Potency}, \code{Potency_pred}); it simply harmonizes
#' columns across inputs (adding \code{NA}-filled columns where necessary),
#' concatenates rows, and reorders a few metadata columns for readability.
#'
#' Any columns named \code{candidate_type}, \code{selection_label}, or
#' \code{Notes_from_SVEMnet} that are present in the final data frame are moved to the
#' leftmost positions in that order.
#'
#' @param ... One or more objects returned by
#'   \code{\link{svem_select_from_score_table}}. You may also pass a single
#'   list of such objects.
#' @param file Character scalar; path to the CSV file to be written.
#'   Required only when \code{write_file = TRUE}.
#' @param overwrite Logical; if \code{FALSE} (default) and \code{file}
#'   already exists, an error is thrown. If \code{TRUE}, any existing file
#'   at \code{file} is overwritten. Only used when \code{write_file = TRUE}.
#' @param write_file Logical; if \code{TRUE} (default), write the combined
#'   table to \code{file} as CSV and print the full path. If \code{FALSE},
#'   no file is written and \code{file} may be \code{NULL}; the concatenated
#'   \code{data.frame} is still returned (invisibly).
#'
#' @return Invisibly, the \code{data.frame} that was written to CSV (or
#'   would be written, when \code{write_file = FALSE}).
#'
#' @examples
#' # Fabricated scored candidates keep this example independent of model fitting.
#' scored <- data.frame(
#'   x1 = c(0.1, 0.4, 0.7, 0.9),
#'   x2 = factor(c("A", "B", "A", "B")),
#'   score = c(0.55, 0.82, 0.91, 0.73),
#'   uncertainty_measure = c(0.15, 0.31, 0.20, 0.44)
#' )
#'
#' optimal <- svem_select_from_score_table(
#'   scored, target = "score", direction = "max", k = 0,
#'   label = "round1_optimal"
#' )
#' exploratory <- svem_select_from_score_table(
#'   scored, target = "uncertainty_measure", direction = "max", k = 0,
#'   label = "round1_exploratory"
#' )
#'
#' # Export to disk when desired:
#' # svem_export_candidates_csv(
#' #   optimal, exploratory, file = "round1_candidates.csv",
#' #   overwrite  = FALSE,
#' #   write_file = TRUE
#' # )
#'
#' # Inspect the combined table without writing a file.
#' cand_tbl <- svem_export_candidates_csv(optimal, exploratory, write_file = FALSE)
#' head(cand_tbl)

#'
#' @export
svem_export_candidates_csv <- function(...,
                                       file       = NULL,
                                       overwrite  = FALSE,
                                       write_file = TRUE) {

  # Validate file only if we are actually writing
  if (isTRUE(write_file)) {
    if (is.null(file) || !is.character(file) ||
        length(file) != 1L || !nzchar(file)) {
      stop("`file` must be a nonempty character scalar giving the CSV path ",
           "when write_file = TRUE.")
    }

    if (file.exists(file) && !isTRUE(overwrite)) {
      stop("File '", file, "' already exists. Set `overwrite = TRUE` to overwrite.")
    }
  }

  sels_raw <- list(...)

  if (!length(sels_raw)) {
    stop("Provide at least one svem_select_from_score_table() object.")
  }

  # helper: is this a selection object?
  is_sel <- function(x) {
    is.list(x) && all(c("best", "candidates") %in% names(x))
  }

  # Allow a single list of selection objects as well:
  # if we got exactly one argument, it's a list, and *all* of its elements
  # are selection objects, treat that as "a list of selections".
  if (length(sels_raw) == 1L &&
      is.list(sels_raw[[1L]]) &&
      all(vapply(sels_raw[[1L]], is_sel, logical(1L)))) {
    sels <- sels_raw[[1L]]
  } else {
    sels <- sels_raw
  }

  if (!all(vapply(sels, is_sel, logical(1L)))) {
    stop("All inputs must be objects returned by svem_select_from_score_table().")
  }

  rows_list <- list()
  idx       <- 0L

  for (sel in sels) {
    # Prefer the evaluated label stored in the selection object; fall back to
    # deparsing the call for objects created by older SVEMnet versions (note
    # that the deparse fallback yields the expression text, not its value,
    # when a variable was passed as `label`).
    sel_label <- NA_character_
    if (!is.null(sel$label) && is.character(sel$label) && length(sel$label) == 1L) {
      sel_label <- sel$label
    } else if (!is.null(sel$call) && is.language(sel$call)) {
      call_obj <- sel$call
      if (!is.null(call_obj$label)) {
        lab_raw <- call_obj$label
        if (is.character(lab_raw) && length(lab_raw) == 1L) {
          sel_label <- lab_raw
        } else if (is.symbol(lab_raw)) {
          sel_label <- as.character(lab_raw)
        } else {
          sel_label <- paste(deparse(lab_raw), collapse = " ")
        }
      }
    }

    best <- sel$best
    if (!is.null(best) && nrow(best) > 0L) {
      best$candidate_type  <- "best"
      best$selection_label <- sel_label
      idx                  <- idx + 1L
      rows_list[[idx]]     <- best
    }

    cand <- sel$candidates
    if (!is.null(cand) && nrow(cand) > 0L) {
      cand$candidate_type  <- "medoid"
      cand$selection_label <- sel_label
      idx                  <- idx + 1L
      rows_list[[idx]]     <- cand
    }
  }

  if (!length(rows_list)) {
    stop("No rows found in the supplied selection objects (no best/candidates).")
  }

  ## ---- Harmonize columns across all selections before rbind ----
  all_cols <- unique(unlist(lapply(rows_list, names)))

  # helper: prototype column for a given name
  get_proto_col <- function(col_name) {
    for (df in rows_list) {
      if (col_name %in% names(df)) {
        return(df[[col_name]])
      }
    }
    NULL
  }

  rows_list <- lapply(rows_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    if (length(missing)) {
      for (m in missing) {
        proto <- get_proto_col(m)
        if (is.null(proto)) {
          # Fallback: default to logical NA when no prototype is found
          df[[m]] <- NA
        } else if (is.factor(proto)) {
          df[[m]] <- factor(NA, levels = levels(proto))
        } else if (is.integer(proto)) {
          df[[m]] <- NA_integer_
        } else if (is.numeric(proto)) {
          df[[m]] <- NA_real_
        } else if (is.logical(proto)) {
          df[[m]] <- NA
        } else {
          df[[m]] <- NA_character_
        }
      }
    }
    df[, all_cols, drop = FALSE]
  })

  # Bind all rows together (now column-aligned)
  out <- do.call(rbind, rows_list)
  rownames(out) <- NULL

  # Move metadata columns to the left (if present)
  meta_order   <- c("candidate_type", "selection_label", "Notes_from_SVEMnet","Notes")
  meta_present <- meta_order[meta_order %in% names(out)]
  if (length(meta_present)) {
    other_cols <- setdiff(names(out), meta_present)
    out <- out[, c(meta_present, other_cols), drop = FALSE]
  }

  # Only write to disk if requested
  if (isTRUE(write_file)) {
    utils::write.csv(out, file = file, row.names = FALSE)

    # Informative message with the normalized (absolute) path
    full_path <- tryCatch(
      normalizePath(file, winslash = "/", mustWork = FALSE),
      error = function(e) file
    )
    message("svem_export_candidates_csv(): wrote candidates to: ", full_path)
  }

  invisible(out)
}
