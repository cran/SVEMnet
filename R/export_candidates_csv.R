#' Export SVEM candidate sets to CSV
#'
#' @description
#' A plan from \code{\link{svem_plan_runs}}, or a scored result from
#' \code{\link{svem_score_random}} with a budgeted
#' \code{$worklist}, or that worklist data frame itself, exports the selected
#' preparations in \code{slot} order. This path preserves \code{recipe_id},
#' \code{slot}, \code{role}, \code{nominated_by}, and all other columns as
#' supplied. It does not relabel exploration rows as best candidates. A direct
#' worklist must have all four identity columns, unique nonmissing recipe IDs,
#' positive unique integer-valued slots, and nonempty roles and endorsements.
#' A short worklist exports only the available preparations; an empty worklist
#' exports column headers. No rows are invented to fill a budget shortfall.
#'
#' Alternatively, given one or more selection objects returned by
#' \code{\link{svem_select_from_score_table}}, concatenate their
#' \code{$best} rows and \code{$candidates} and export a CSV suitable for
#' planning new experimental runs.
#'
#' On the legacy selection-object path, each row is tagged with:
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
#' On that legacy path, columns named \code{candidate_type}, \code{selection_label}, or
#' \code{Notes_from_SVEMnet} that are present in the final data frame are moved to the
#' leftmost positions in that order.
#'
#' @param ... One or more objects returned by
#'   \code{\link{svem_select_from_score_table}}. You may also pass a single
#'   list of such objects. Alternatively, pass one run plan, budgeted scored result or
#'   its recognizable worklist data frame. These worklist forms cannot be
#'   combined with other arguments in \code{...}. A score result without a
#'   worklist is not exported: request \code{budget} when scoring or select
#'   rows with \code{svem_select_from_score_table()} first.
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
#' @param planned_date Optional planned preparation date: one \code{Date} or
#'   an ISO date string such as \code{"2026-11-02"}. Adds a
#'   \code{planned_date} column. This is scheduling metadata, not the actual
#'   execution date, a predictor, or an automatically added blocking factor.
#' @param preparation_prefix Optional nonempty string identifying this batch,
#'   for example \code{"studyA-20261102-round1"}. Adds \code{preparation_id}
#'   values formed from the prefix and \code{-001}, \code{-002}, and so on,
#'   in exported row order. Use a distinct prefix for every physical batch,
#'   including batches on the same date. IDs are unique within this export;
#'   the function does not maintain a global laboratory ID registry. Re-exporting
#'   the same plan with the same prefix reproduces its labels. Exporting an
#'   unlabeled subset renumbers its rows; to preserve previously assigned IDs,
#'   filter the returned labeled table and omit these labeling options on export.
#'   For a laboratory's
#'   own IDs, add a \code{preparation_id} column to the worklist instead.
#'   Existing columns are preserved; supplying either option when its output
#'   column already exists errors rather than overwriting it. With both options
#'   omitted, export columns and labels are unchanged, including legacy exports.
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
#' # Optional scheduling labels, separate from the formulation's recipe_id:
#' labelled <- svem_export_candidates_csv(optimal, exploratory, write_file = FALSE,
#'   planned_date = as.Date("2026-11-02"),
#'   preparation_prefix = "studyA-20261102-round1")
#' labelled[c("planned_date", "preparation_id", "selection_label")]

#'
#' @export
svem_export_candidates_csv <- function(...,
                                       file       = NULL,
                                       overwrite  = FALSE,
                                       write_file = TRUE,
                                       planned_date = NULL,
                                       preparation_prefix = NULL) {

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

  write_output <- function(out) {
    out <- .svem_export_labels(out, planned_date, preparation_prefix)
    if (isTRUE(write_file)) {
      utils::write.csv(out, file = file, row.names = FALSE)
      full_path <- tryCatch(
        normalizePath(file, winslash = "/", mustWork = FALSE),
        error = function(e) file
      )
      message("svem_export_candidates_csv(): wrote candidates to: ", full_path)
    }
    invisible(out)
  }

  sels_raw <- list(...)
  if (length(sels_raw) == 1L) {
    input <- sels_raw[[1L]]
    if (inherits(input, "svem_run_plan")) {
      if (is.null(input$worklist)) stop("This run plan has no worklist to export.", call. = FALSE)
      return(write_output(.svem_export_worklist(input$worklist)))
    }
    if (is.list(input) && !is.data.frame(input) && "score_table" %in% names(input)) {
      if (is.null(input$worklist)) {
        stop("This scored result has no worklist. Request `budget` in ",
             "svem_score_random(), or select rows with svem_select_from_score_table() ",
             "and export that selection. The entire candidate pool is not exported.",
             call. = FALSE)
      }
      return(write_output(.svem_export_worklist(input$worklist)))
    }
    if (is.data.frame(input)) {
      return(write_output(.svem_export_worklist(input)))
    }
  }

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
    protected <- unique(c(attr(sel$best, "svem_predictor_cols"),
                          attr(sel$best, "svem_response_vars"),
                          attr(sel$candidates, "svem_predictor_cols"),
                          attr(sel$candidates, "svem_response_vars")))
    .svem_check_output_collisions(protected,
                                  c("candidate_type", "selection_label"),
                                  "svem_export_candidates_csv()")
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

  write_output(out)
}

# Recognize a worklist by its complete identity schema, not by prediction or
# score columns that also occur in an unselected candidate pool.
.svem_export_worklist <- function(work) {
  identity <- c("recipe_id", "slot", "role", "nominated_by")
  if (!is.data.frame(work) || anyDuplicated(names(work)) ||
      !all(identity %in% names(work))) {
    stop("A worklist must be a data frame with unique column names including ",
         "recipe_id, slot, role, and nominated_by. Export a budgeted scored ",
         "result or select rows with svem_select_from_score_table() first.",
         call. = FALSE)
  }
  ids <- work$recipe_id
  if (!is.atomic(ids) || !is.null(dim(ids)) || anyNA(ids) ||
      any(!nzchar(as.character(ids))) || anyDuplicated(as.character(ids))) {
    stop("Worklist `recipe_id` values must be unique, nonmissing, and nonempty.",
         call. = FALSE)
  }
  slots <- work$slot
  if (!is.numeric(slots) || !is.null(dim(slots)) || anyNA(slots) ||
      any(!is.finite(slots)) || any(slots < 1 | slots != floor(slots)) ||
      anyDuplicated(slots)) {
    stop("Worklist `slot` values must be unique positive whole numbers.", call. = FALSE)
  }
  for (name in c("role", "nominated_by")) {
    value <- work[[name]]
    if ((!is.character(value) && !is.factor(value)) || !is.null(dim(value)) ||
        anyNA(value) || any(!nzchar(trimws(as.character(value))))) {
      stop("Worklist `", name, "` values must be nonmissing, nonempty labels.",
           call. = FALSE)
    }
  }
  out <- work[order(slots, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Explicit laboratory labels are applied only on export, never to model input.
.svem_export_labels <- function(out, planned_date, preparation_prefix) {
  if (is.null(planned_date) && is.null(preparation_prefix)) return(out)
  if (!is.null(planned_date)) {
    if ("planned_date" %in% names(out)) {
      stop("The export already has planned_date; preserve it or remove that column explicitly.", call. = FALSE)
    }
    if (length(planned_date) != 1L || is.na(planned_date)) {
      stop("planned_date must be one nonmissing Date or ISO date string.", call. = FALSE)
    }
    if (is.character(planned_date)) {
      if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", planned_date)) {
        stop("planned_date must use the unambiguous YYYY-MM-DD format.", call. = FALSE)
      }
      value <- tryCatch(as.Date(planned_date, format = "%Y-%m-%d"), error = function(e) NA)
      if (is.na(value) || !identical(format(value, "%Y-%m-%d"), planned_date)) {
        stop("planned_date must use the unambiguous YYYY-MM-DD format.", call. = FALSE)
      }
      planned_date <- value
    }
    if (!inherits(planned_date, "Date") || !is.finite(as.numeric(planned_date))) {
      stop("planned_date must be one nonmissing Date or ISO date string.", call. = FALSE)
    }
    out$planned_date <- rep(planned_date, nrow(out))
  }
  if (!is.null(preparation_prefix)) {
    if (!is.character(preparation_prefix) || length(preparation_prefix) != 1L ||
        is.na(preparation_prefix) || !nzchar(trimws(preparation_prefix)) ||
        grepl("[\r\n]", preparation_prefix)) {
      stop("preparation_prefix must be one nonempty single-line string.", call. = FALSE)
    }
    if ("preparation_id" %in% names(out)) {
      stop("The export already has preparation_id; preserve it or remove that column explicitly.", call. = FALSE)
    }
    out$preparation_id <- if (nrow(out)) paste0(preparation_prefix, "-",
      sprintf("%03d", seq_len(nrow(out)))) else character()
  }
  out
}
